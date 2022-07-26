MODULE UFEMISM_main_model
  ! The main regional ice-sheet model

#include <petsc/finclude/petscksp.h>

! ===== USEs =====
! ================

  USE mpi
  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: c_backspace
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list
  USE data_types_module,               ONLY: type_model_region, type_grid, type_remapping_mesh_mesh
  USE reference_fields_module,         ONLY: initialise_reference_geometries, map_reference_geometries_to_mesh
  USE mesh_memory_module,              ONLY: deallocate_mesh_all
  USE mesh_creation_module,            ONLY: create_mesh_from_cart_data
  USE mesh_mapping_module,             ONLY: calc_remapping_operators_mesh_mesh, deallocate_remapping_operators_mesh_mesh, &
                                             calc_remapping_operator_mesh2grid, deallocate_remapping_operators_mesh2grid, &
                                             calc_remapping_operator_grid2mesh, deallocate_remapping_operators_grid2mesh
  USE mesh_update_module,              ONLY: determine_mesh_fitness, create_new_mesh
  use mesh_single_module,              only: create_new_mesh_single, create_single_mesh_from_cart_data
  USE netcdf_module,                   ONLY: initialise_debug_fields, create_output_files, associate_debug_fields, &
                                             write_to_output_files, create_debug_file, reallocate_debug_fields
  USE ice_dynamics_module,             ONLY: initialise_ice_model, remap_ice_model, run_ice_model
  use reallocate_mod,                  only: reallocate
  use utilities_module,                only: time_display, inverse_oblique_sg_projection

! ===== Preamble =====
! ====================

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  subroutine run_model( region, t_end)
    ! Run the model region until the next coupling time

    implicit none

    ! In/output variables:
    type(type_model_region), intent(inout) :: region
    real(dp),                intent(in)    :: t_end

    ! Local variables:
    character(len=256)                     :: routine_name
    real(dp)                               :: meshfitness
    integer                                :: it
    real(dp)                               :: t1, t2, dt_ave

    ! Add routine to path
    routine_name = 'run_model('  //  region%name  //  ')'
    call init_routine( routine_name)

    if (par%master) then
      write(*,"(A)") ''
      write(*,"(A,A,A,A,A,F9.3,A,F9.3,A)") &
            ' Running model region ', region%name, ' (', TRIM(region%long_name), &
            ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'
    end if

    ! Set the intermediary pointers in "debug" to this region's debug data fields
    call associate_debug_fields(  region)

    ! Computation time tracking
    region%tcomp_total          = 0._dp
    region%tcomp_ice            = 0._dp
    region%tcomp_thermo         = 0._dp
    region%tcomp_climate        = 0._dp
    region%tcomp_GIA            = 0._dp
    region%tcomp_mesh           = 0._dp

    t1 = MPI_WTIME()
    t2 = 0._dp

    ! Initialise iteration counter
    it = 0
    ! Initialise averaged time step
    dt_ave = 0._dp

    ! ====================================
    ! ===== The main model time loop =====
    ! ====================================

    do while (region%time < t_end)

      ! Update iteration counter
      it = it + 1

      ! Mesh update
      ! ===========

      ! Check if the mesh needs to be updated
      t2 = MPI_WTIME()
      meshfitness = 1._dp
      if (region%time > region%t_last_mesh + C%dt_mesh_min) then
        call determine_mesh_fitness(region%mesh, region%ice, meshfitness)
      end if
      region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2

      ! If required, update the mesh
      if (meshfitness < C%mesh_fitness_threshold) then
        region%t_last_mesh = region%time
        t2 = MPI_WTIME()
        call run_model_update_mesh( region)
        region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2
      end if

      ! Ice dynamics
      ! ============

      ! Calculate ice velocities and the resulting change in ice geometry
      ! NOTE: geometry is not updated yet; this happens at the end of the time loop
      t1 = MPI_WTIME()
      call run_ice_model( region, t_end)
      t2 = MPI_WTIME()
      region%tcomp_ice = region%tcomp_ice + t2 - t1

      ! == Time display
      ! ===============

      if (par%master .AND. C%do_time_display) then
        call time_display(region, t_end, dt_ave, it)
      end if

      ! == Output
      ! =========

      ! Write NetCDF output
      if (region%do_output) then
        ! If the mesh has been updated, create a new NetCDF file
        if (.not. region%output_file_exists) then
          call create_output_files( region)
          call sync
          region%output_file_exists = .true.
        end if
        ! Write to regional NetCDF output files
        call write_to_output_files( region)
      end if

      ! Advance region time
      ! ===================

      region%time = region%time + region%dt
      dt_ave = dt_ave + region%dt
      call sync

    end do ! while (region%time < t_end)

    ! ===========================================
    ! ===== End of the main model time loop =====
    ! ===========================================

    ! Write to NetCDF output one last time at the end of the simulation
    if (region%time == C%end_time_of_run) then
      ! If the mesh has been updated, create a new NetCDF file
      if (.not. region%output_file_exists) then
        call create_output_files( region)
        call sync
        region%output_file_exists = .TRUE.
      end if
      call write_to_output_files( region)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_model

  subroutine run_model_update_mesh( region)
    ! Perform a mesh update: create a new mesh based on the current modelled ice-sheet
    ! geometry, map all the model data from the old to the new mesh. Deallocate the
    ! old mesh, update the output files and the square grid maps.

    implicit none

    ! In- and output variables
    type(type_model_region),             intent(inout) :: region

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'run_model_update_mesh'
    type(type_remapping_mesh_mesh)                     :: map

    ! Add routine to path
    call init_routine( routine_name)

    ! Create a new mesh
    if (C%use_submesh) then
      call create_new_mesh( region)
    else
      call create_new_mesh_single( region)
    end if

    if (par%master) then
      write(*,"(A)") '  Reallocating and remapping after mesh update...'
    end if

    ! Update the mapping operators between the new mesh and the fixed square grids
    call deallocate_remapping_operators_mesh2grid( region%grid_output)
    call deallocate_remapping_operators_mesh2grid( region%grid_GIA   )
    call deallocate_remapping_operators_mesh2grid( region%grid_smooth)

    call deallocate_remapping_operators_grid2mesh( region%grid_output)
    call deallocate_remapping_operators_grid2mesh( region%grid_GIA   )
    call deallocate_remapping_operators_grid2mesh( region%grid_smooth)

    call calc_remapping_operator_mesh2grid( region%mesh_new, region%grid_output)
    call calc_remapping_operator_mesh2grid( region%mesh_new, region%grid_GIA   )
    call calc_remapping_operator_mesh2grid( region%mesh_new, region%grid_smooth)

    call calc_remapping_operator_grid2mesh( region%grid_output, region%mesh_new)
    call calc_remapping_operator_grid2mesh( region%grid_GIA   , region%mesh_new)
    call calc_remapping_operator_grid2mesh( region%grid_smooth, region%mesh_new)

    ! Calculate the mapping arrays
    call calc_remapping_operators_mesh_mesh( region%mesh, region%mesh_new, map)

    ! Reallocate memory for reference geometries
    if (allocated(region%refgeo_init%Hi)) then
      deallocate( region%refgeo_init%Hi )
      deallocate( region%refgeo_init%Hb )
      deallocate( region%refgeo_init%Hs )

      deallocate( region%refgeo_PD%Hi)
      deallocate( region%refgeo_PD%Hb)
      deallocate( region%refgeo_PD%Hs)

      deallocate( region%refgeo_GIAeq%Hi)
      deallocate( region%refgeo_GIAeq%Hb)
      deallocate( region%refgeo_GIAeq%Hs)
    end if

    allocate( region%refgeo_init%Hi (region%mesh_new%vi1:region%mesh_new%vi2), source=0.0_dp)
    allocate( region%refgeo_init%Hb (region%mesh_new%vi1:region%mesh_new%vi2), source=0.0_dp)
    allocate( region%refgeo_init%Hs (region%mesh_new%vi1:region%mesh_new%vi2), source=0.0_dp)

    allocate( region%refgeo_PD%Hi   (region%mesh_new%vi1:region%mesh_new%vi2), source=0.0_dp)
    allocate( region%refgeo_PD%Hb   (region%mesh_new%vi1:region%mesh_new%vi2), source=0.0_dp)
    allocate( region%refgeo_PD%Hs   (region%mesh_new%vi1:region%mesh_new%vi2), source=0.0_dp)

    allocate( region%refgeo_GIAeq%Hi(region%mesh_new%vi1:region%mesh_new%vi2), source=0.0_dp)
    allocate( region%refgeo_GIAeq%Hb(region%mesh_new%vi1:region%mesh_new%vi2), source=0.0_dp)
    allocate( region%refgeo_GIAeq%Hs(region%mesh_new%vi1:region%mesh_new%vi2), source=0.0_dp)

    ! Map reference geometries from the square grids to the mesh
    call map_reference_geometries_to_mesh( region, region%mesh_new)

    ! Remap all the submodels
    call remap_ice_model( region%mesh, region%mesh_new, map, region%ice, region%refgeo_PD, region%time)

    ! Deallocate shared memory for the mapping arrays
    call deallocate_remapping_operators_mesh_mesh( map)

    ! Deallocate the old mesh, bind the region%mesh pointers to the new mesh.
    !call deallocate_mesh_all( region%mesh) TOdo make this automatic, fix it
    region%mesh = region%mesh_new

    ! When the next output is written, new output files must be created.
    region%output_file_exists = .FALSE.

    ! Reallocate the debug fields for the new mesh, create a new debug file
    call reallocate_debug_fields( region)
    call create_debug_file(       region)

    ! Run all model components again after updating the mesh
    region%t_next_SIA     = region%time
    region%t_next_SSA     = region%time
    region%t_next_DIVA    = region%time
    region%t_next_thermo  = region%time
    region%t_next_climate = region%time
    region%t_next_ocean   = region%time
    region%t_next_SMB     = region%time
    region%t_next_BMB     = region%time
    region%t_next_ELRA    = region%time

    region%do_SIA         = .true.
    region%do_SSA         = .true.
    region%do_DIVA        = .true.
    region%do_thermo      = .true.
    region%do_climate     = .true.
    region%do_ocean       = .true.
    region%do_SMB         = .true.
    region%do_BMB         = .true.
    region%do_ELRA        = .true.

    if (par%master) then
      write(*,"(A)") '  Finished reallocating and remapping.'
      write(*,"(A)") '  Running again now...'
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_model_update_mesh

  SUBROUTINE initialise_model( region, name)
    ! Initialise the entire model region - read initial and PD data, create the mesh,
    ! initialise the ice dynamics, climate, ocean, and SMB sub models

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),          INTENT(INOUT)     :: region
    CHARACTER(LEN=3),                 INTENT(IN)        :: name

    ! Local variables:
    CHARACTER(LEN=256)                                  :: routine_name

    ! Add routine to path
    routine_name = 'initialise_model('  //  name  //  ')'
    call init_routine( routine_name)

    ! ===== Basic initialisation =====
    ! ================================

    ! Region name
    region%name      = name
    if (region%name == 'NAM') THEN
      region%long_name = 'North America'
    else if (region%name == 'EAS') THEN
      region%long_name = 'Eurasia'
    else if (region%name == 'GRL') THEN
      region%long_name = 'Greenland'
    else if (region%name == 'ANT') THEN
      region%long_name = 'Antarctica'
    end if

    if (par%master) then
      write(*,"(A)") ''
      write(*,"(5A)") ' Initialising model region ', region%name, &
                      ' (', TRIM(region%long_name), ')...'
    end if

    ! ===== Allocate memory for timers and scalars =====
    ! ==================================================

    call allocate_region_timers_and_scalars( region)

    ! ===== PD and init reference data fields =====
    ! =============================================

    call initialise_reference_geometries( region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq, region%name)

    ! ===== The mesh =====
    ! ====================

    if (C%use_submesh) then
      call create_mesh_from_cart_data( region)
    else
      call create_single_mesh_from_cart_data( region)
    endif

    ! ===== Map reference geometries to the mesh =====
    ! ================================================

    ! Allocate memory for reference data on the mesh
    allocate( region%refgeo_init%Hi (region%mesh%vi1:region%mesh%vi2), source=0.0_dp)
    allocate( region%refgeo_init%Hb (region%mesh%vi1:region%mesh%vi2), source=0.0_dp)
    allocate( region%refgeo_init%Hs (region%mesh%vi1:region%mesh%vi2), source=0.0_dp)

    allocate( region%refgeo_PD%Hi   (region%mesh%vi1:region%mesh%vi2), source=0.0_dp)
    allocate( region%refgeo_PD%Hb   (region%mesh%vi1:region%mesh%vi2), source=0.0_dp)
    allocate( region%refgeo_PD%Hs   (region%mesh%vi1:region%mesh%vi2), source=0.0_dp)

    allocate( region%refgeo_GIAeq%Hi(region%mesh%vi1:region%mesh%vi2), source=0.0_dp)
    allocate( region%refgeo_GIAeq%Hb(region%mesh%vi1:region%mesh%vi2), source=0.0_dp)
    allocate( region%refgeo_GIAeq%Hs(region%mesh%vi1:region%mesh%vi2), source=0.0_dp)

    ! Map data from the square grids to the mesh
    call map_reference_geometries_to_mesh( region, region%mesh)

    ! ===== The different square grids =====
    ! ======================================

    if (par%master) then
      write(*,"(A)") '  Initialising square grids for output, GIA, and data smoothing...'
    end if

    call initialise_model_square_grid( region, region%grid_output, C%dx_grid_output)
    call initialise_model_square_grid( region, region%grid_GIA,    C%dx_grid_GIA   )
    call initialise_model_square_grid( region, region%grid_smooth, C%dx_grid_smooth)

    ! ===== Initialise dummy fields for debugging =====
    ! =================================================

    if (par%master) then
      write(*,"(A)") '  Initialising debug fields...'
    end if

    call initialise_debug_fields( region)

    ! ===== Output files =====
    ! ========================

    call create_output_files(    region)
    call associate_debug_fields( region)
    region%output_file_exists = .TRUE.
    call sync

    ! ===== The ice dynamics model =====
    ! ==================================

    call initialise_ice_model( region%mesh, region%ice, region%refgeo_init)

    if (par%master) then
      write(*,"(3A)") ' Finished initialising model region ', region%name, '.'
    end if

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = HUGE( 1))

  END SUBROUTINE initialise_model

! ===== Auxiliary routines =====
! ==============================

  SUBROUTINE allocate_region_timers_and_scalars( region)
    ! Allocate shared memory for this region's timers (used for the asynchronous coupling between the
    ! ice dynamics and the secondary model components), and for the scalars (integrated ice volume and
    ! area, SMB components, computation times, etc.)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_region_timers_and_scalars'

    ! Add routine to path
    call init_routine( routine_name)

    ! Timers and time steps
    ! =====================
    region%time           = C%start_time_of_run
    region%dt             = C%dt_min
    region%dt_prev        = C%dt_min

    region%t_last_mesh    = C%start_time_of_run
    region%t_next_mesh    = C%start_time_of_run + C%dt_mesh_min
    region%do_mesh        = .FALSE.

    region%t_last_SIA     = C%start_time_of_run
    region%t_next_SIA     = C%start_time_of_run
    region%do_SIA         = .TRUE.

    region%t_last_SSA     = C%start_time_of_run
    region%t_next_SSA     = C%start_time_of_run
    region%do_SSA         = .TRUE.

    region%t_last_DIVA    = C%start_time_of_run
    region%t_next_DIVA    = C%start_time_of_run
    region%do_DIVA        = .TRUE.

    region%t_last_thermo  = C%start_time_of_run
    region%t_next_thermo  = C%start_time_of_run + C%dt_thermo
    region%do_thermo      = .FALSE.

    region%t_last_climate = C%start_time_of_run
    region%t_next_climate = C%start_time_of_run
    region%do_climate     = .TRUE.

    region%t_last_ocean   = C%start_time_of_run
    region%t_next_ocean   = C%start_time_of_run
    region%do_ocean       = .TRUE.

    region%t_last_SMB     = C%start_time_of_run
    region%t_next_SMB     = C%start_time_of_run
    region%do_SMB         = .TRUE.

    region%t_last_BMB     = C%start_time_of_run
    region%t_next_BMB     = C%start_time_of_run
    region%do_BMB         = .TRUE.

    region%t_last_ELRA    = C%start_time_of_run
    region%t_next_ELRA    = C%start_time_of_run
    region%do_ELRA        = .TRUE.

    region%t_last_output  = C%start_time_of_run
    region%t_next_output  = C%start_time_of_run
    region%do_output      = .TRUE.

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = 65)

  END SUBROUTINE allocate_region_timers_and_scalars

  SUBROUTINE initialise_model_square_grid( region, grid, dx)
    ! Initialise a regular square grid enveloping this model region

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    TYPE(type_grid),            INTENT(INOUT)     :: grid
    REAL(dp),                   INTENT(IN)        :: dx

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_model_square_grid'
    REAL(dp)                                      :: xmid, ymid
    INTEGER                                       :: nsx, nsy, i, j, n
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp

    ! Add routine to path
    call init_routine( routine_name)

    nsx = 0
    nsy = 0

    ! Resolution
    grid%dx = dx

    ! Determine the center of the model domain
    xmid = (region%mesh%xmin + region%mesh%xmax) / 2._dp
    ymid = (region%mesh%ymin + region%mesh%ymax) / 2._dp

    nsx = CEILING((region%mesh%xmax - xmid - grid%dx/2._dp) / grid%dx)
    nsy = CEILING((region%mesh%ymax - ymid - grid%dx/2._dp) / grid%dx)

    grid%nx = 2*nsx + 1
    grid%ny = 2*nsy + 1

    allocate( grid%x ( grid%nx ))
    allocate( grid%y ( grid%ny ))

    do i = 1, grid%nx
      grid%x( i) = -nsx*grid%dx + (i-1)*grid%dx
    end do
    do j = 1, grid%ny
      grid%y( j) = -nsy*grid%dx + (j-1)*grid%dx
    end do

    grid%xmin = grid%x(1      )
    grid%xmax = grid%x(grid%nx)
    grid%ymin = grid%y(1      )
    grid%ymax = grid%y(grid%ny)

    ! Tolerance; points lying within this distance of each other are treated as identical
    grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

    ! Set up grid-to-vector translation tables
    grid%n  = grid%nx * grid%ny
    allocate ( grid%ij2n( grid%nx, grid%ny ))
    allocate ( grid%n2ij( grid%n , 2       ))
    n = 0
    do i = 1, grid%nx
        do j = 1, grid%ny
          n = n+1
          grid%ij2n( i,j) = n
          grid%n2ij( n,:) = [i,j]
        end do
    end do

    ! Assign range to each processor
    call partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    call partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

    ! Calculate lat-lon coordinates
    allocate( grid%lat(grid%nx, grid%ny))
    allocate( grid%lon(grid%nx, grid%ny))

    do i = 1, grid%nx
    do j = 1, grid%ny
      call inverse_oblique_sg_projection( grid%x( i), grid%y( j), region%mesh%lambda_M, region%mesh%phi_M, region%mesh%alpha_stereo, grid%lon( i,j), grid%lat( i,j))
    end do
    end do

    ! Calculate mapping arrays between the mesh and the grid
    call calc_remapping_operator_mesh2grid( region%mesh, grid)
    call calc_remapping_operator_grid2mesh( grid, region%mesh)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = 15)

  END SUBROUTINE initialise_model_square_grid

END MODULE UFEMISM_main_model
