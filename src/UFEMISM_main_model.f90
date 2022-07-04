MODULE UFEMISM_main_model

  ! The main regional ice-sheet model

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list

  ! ! Import specific functionality
  USE data_types_module,               ONLY: type_model_region, type_grid, type_remapping_mesh_mesh
  USE reference_fields_module,         ONLY: initialise_reference_geometries, map_reference_geometries_to_mesh
  USE mesh_memory_module,              ONLY: deallocate_mesh_all
  USE mesh_help_functions_module,      ONLY: inverse_oblique_sg_projection
  USE mesh_creation_module,            ONLY: create_mesh_from_cart_data
  USE mesh_mapping_module,             ONLY: calc_remapping_operators_mesh_mesh, deallocate_remapping_operators_mesh_mesh, &
                                             calc_remapping_operator_mesh2grid, deallocate_remapping_operators_mesh2grid, &
                                             calc_remapping_operator_grid2mesh, deallocate_remapping_operators_grid2mesh
  USE mesh_update_module,              ONLY: determine_mesh_fitness, create_new_mesh
  USE netcdf_module,                   ONLY: initialise_debug_fields, create_output_files, associate_debug_fields, &
                                             write_to_output_files, create_debug_file, reallocate_debug_fields
  USE ice_dynamics_module,             ONLY: initialise_ice_model, remap_ice_model, determine_timesteps_and_actions
  use reallocate_mod,                  only: reallocate

  IMPLICIT NONE

CONTAINS

  SUBROUTINE run_model( region, t_end)
    ! Run the model until t_end (usually a 100 years further)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),           INTENT(INOUT)     :: region
    REAL(dp),                          INTENT(IN)        :: t_end

    ! Local variables:
    CHARACTER(LEN=256)                                   :: routine_name
    INTEGER                                              :: it
    REAL(dp)                                             :: meshfitness
    REAL(dp)                                             :: t1, t2

    ! Add routine to path
    routine_name = 'run_model('  //  region%name  //  ')'
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE (0,'(A,A,A,A,A,F9.3,A,F9.3,A)') '  Running model region ', region%name, ' (', TRIM(region%long_name), &
                                                          ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'

    ! Set the intermediary pointers in "debug" to this region's debug data fields
    CALL associate_debug_fields(  region)

    ! Computation time tracking
    region%tcomp_total          = 0._dp
    region%tcomp_ice            = 0._dp
    region%tcomp_thermo         = 0._dp
    region%tcomp_climate        = 0._dp
    region%tcomp_GIA            = 0._dp
    region%tcomp_mesh           = 0._dp

    t1 = MPI_WTIME()
    t2 = 0._dp

  ! ====================================
  ! ===== The main model time loop =====
  ! ====================================

    it = 0
    DO WHILE (region%time < t_end)
      it = it + 1

    ! Mesh update
    ! ===========

      ! Check if the mesh needs to be updated
      t2 = MPI_WTIME()
      meshfitness = 1._dp
      IF (region%time > region%t_last_mesh + C%dt_mesh_min) THEN
        CALL determine_mesh_fitness(region%mesh, region%ice, meshfitness)
      END IF
      region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2

      ! If required, update the mesh
      IF (meshfitness < C%mesh_fitness_threshold) THEN
      ! IF (.FALSE.) THEN
      ! IF (.TRUE.) THEN
        region%t_last_mesh = region%time
        t2 = MPI_WTIME()
        CALL run_model_update_mesh( region)
        region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2
      END IF

    ! Determine time step
    ! ===================

      ! Adjust the time step to prevent overshooting other model components
      CALL determine_timesteps_and_actions( region, t_end)

    ! Output
    ! ======

      ! Write output
      IF (region%do_output) THEN
        ! If the mesh has been updated, create a new NetCDF file
        IF (.NOT. region%output_file_exists) THEN
          CALL create_output_files( region)
          CALL sync
          region%output_file_exists = .TRUE.
        END IF
        CALL write_to_output_files( region)
      END IF

    ! Advance region time
    ! ===================

      region%time = region%time + region%dt
      CALL sync

    END DO ! DO WHILE (region%time < t_end)

  ! ===========================================
  ! ===== End of the main model time loop =====
  ! ===========================================

    ! Write to NetCDF output one last time at the end of the simulation
    IF (region%time == C%end_time_of_run) THEN
      ! If the mesh has been updated, create a new NetCDF file
      IF (.NOT. region%output_file_exists) THEN
        CALL create_output_files( region)
        CALL sync
        region%output_file_exists = .TRUE.
      END IF
      CALL write_to_output_files( region)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model

  ! Update the mesh and everything else that needs it
  SUBROUTINE run_model_update_mesh( region)
    ! Perform a mesh update: create a new mesh based on the current modelled ice-sheet
    ! geometry, map all the model data from the old to the new mesh. Deallocate the
    ! old mesh, update the output files and the square grid maps.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_model_update_mesh'
    TYPE(type_remapping_mesh_mesh)                     :: map

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new mesh
    CALL create_new_mesh( region)

    ! Update the mapping operators between the new mesh and the fixed square grids
    CALL deallocate_remapping_operators_mesh2grid(           region%grid_output)
    CALL deallocate_remapping_operators_mesh2grid(           region%grid_GIA   )
    CALL deallocate_remapping_operators_mesh2grid(           region%grid_smooth)

    CALL deallocate_remapping_operators_grid2mesh(           region%grid_output)
    CALL deallocate_remapping_operators_grid2mesh(           region%grid_GIA   )
    CALL deallocate_remapping_operators_grid2mesh(           region%grid_smooth)

    CALL calc_remapping_operator_mesh2grid( region%mesh_new, region%grid_output)
    CALL calc_remapping_operator_mesh2grid( region%mesh_new, region%grid_GIA   )
    CALL calc_remapping_operator_mesh2grid( region%mesh_new, region%grid_smooth)

    CALL calc_remapping_operator_grid2mesh( region%grid_output, region%mesh_new)
    CALL calc_remapping_operator_grid2mesh( region%grid_GIA   , region%mesh_new)
    CALL calc_remapping_operator_grid2mesh( region%grid_smooth, region%mesh_new)

    ! Calculate the mapping arrays
    CALL calc_remapping_operators_mesh_mesh( region%mesh, region%mesh_new, map)

    ! Reallocate memory for reference geometries
    if (allocated(region%refgeo_init%Hi)) then
      deallocate( region%refgeo_init%Hi )
      deallocate( region%refgeo_init%Hb )
      deallocate( region%refgeo_init%Hs )

      deallocate( region%refgeo_PD%Hi   )
      deallocate( region%refgeo_PD%Hb   )
      deallocate( region%refgeo_PD%Hs   )

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
    CALL map_reference_geometries_to_mesh( region, region%mesh_new)

    ! Remap all the submodels
    CALL remap_ice_model(      region%mesh, region%mesh_new, map, region%ice, region%refgeo_PD, region%time)

    ! Deallocate shared memory for the mapping arrays
    CALL deallocate_remapping_operators_mesh_mesh( map)

    ! Deallocate the old mesh, bind the region%mesh pointers to the new mesh.
    !CALL deallocate_mesh_all( region%mesh) TODO make this automatic, fix it
    region%mesh = region%mesh_new

    ! When the next output is written, new output files must be created.
    region%output_file_exists = .FALSE.

    ! Reallocate the debug fields for the new mesh, create a new debug file
    CALL reallocate_debug_fields( region)
    CALL create_debug_file(       region)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model_update_mesh

  ! Initialise the entire model region - read initial and PD data, create the mesh,
  ! initialise the ice dynamics, climate, ocean, and SMB sub models
  SUBROUTINE initialise_model( region, name)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),          INTENT(INOUT)     :: region
    CHARACTER(LEN=3),                 INTENT(IN)        :: name

    ! Local variables:
    CHARACTER(LEN=256)                                  :: routine_name

    ! Add routine to path
    routine_name = 'initialise_model('  //  name  //  ')'
    CALL init_routine( routine_name)

    ! ===== Basic initialisation =====
    ! ================================

    ! Region name
    region%name      = name
    IF (region%name == 'NAM') THEN
      region%long_name = 'North America'
    ELSE IF (region%name == 'EAS') THEN
      region%long_name = 'Eurasia'
    ELSE IF (region%name == 'GRL') THEN
      region%long_name = 'Greenland'
    ELSE IF (region%name == 'ANT') THEN
      region%long_name = 'Antarctica'
    END IF

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Initialising model region ', region%name, ' (', TRIM(region%long_name), ')...'

    ! ===== Allocate memory for timers and scalars =====
    ! ==================================================

    CALL allocate_region_timers_and_scalars( region)

    ! ===== PD and init reference data fields =====
    ! =============================================

    CALL initialise_reference_geometries( region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq, region%name)

    ! ===== The mesh =====
    ! ====================

    CALL create_mesh_from_cart_data( region)

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
    CALL map_reference_geometries_to_mesh( region, region%mesh)

    ! ===== The different square grids =====
    ! ======================================

    IF (par%master) WRITE(0,*) '  Initialising square grids for output, GIA, and data smoothing...'

    CALL initialise_model_square_grid( region, region%grid_output, C%dx_grid_output)
    CALL initialise_model_square_grid( region, region%grid_GIA,    C%dx_grid_GIA   )
    CALL initialise_model_square_grid( region, region%grid_smooth, C%dx_grid_smooth)

    ! ===== Initialise dummy fields for debugging =====
    ! =================================================

    IF (par%master) WRITE(0,*) '  Initialising debug fields...'

    CALL initialise_debug_fields( region)

    ! ===== Output files =====
    ! ========================

    CALL create_output_files(    region)
    CALL associate_debug_fields( region)
    region%output_file_exists = .TRUE.
    CALL sync

    ! ===== The ice dynamics model =====
    ! ==================================

    CALL initialise_ice_model( region%mesh, region%ice, region%refgeo_init)

    IF (par%master) WRITE (0,*) ' Finished initialising model region ', region%name, '.'

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = HUGE( 1))

  END SUBROUTINE initialise_model
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
    CALL init_routine( routine_name)

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
    CALL finalise_routine( routine_name, n_extra_windows_expected = 65)
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
    CALL init_routine( routine_name)

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

    DO i = 1, grid%nx
      grid%x( i) = -nsx*grid%dx + (i-1)*grid%dx
    END DO
    DO j = 1, grid%ny
      grid%y( j) = -nsy*grid%dx + (j-1)*grid%dx
    END DO

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
    DO i = 1, grid%nx
      IF (MOD(i,2) == 1) THEN
        DO j = 1, grid%ny
          n = n+1
          grid%ij2n( i,j) = n
          grid%n2ij( n,:) = [i,j]
        END DO
      ELSE
        DO j = grid%ny, 1, -1
          n = n+1
          grid%ij2n( i,j) = n
          grid%n2ij( n,:) = [i,j]
        END DO
      END IF
    END DO

    ! Assign range to each processor
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

    ! Calculate lat-lon coordinates
    allocate( grid%lat(grid%nx, grid%ny))
    allocate( grid%lon(grid%nx, grid%ny))

    DO i = 1, grid%nx
    DO j = 1, grid%ny
      CALL inverse_oblique_sg_projection( grid%x( i), grid%y( j), region%mesh%lambda_M, region%mesh%phi_M, region%mesh%alpha_stereo, grid%lon( i,j), grid%lat( i,j))
    END DO
    END DO

    ! Calculate mapping arrays between the mesh and the grid
    CALL calc_remapping_operator_mesh2grid( region%mesh, grid)
    CALL calc_remapping_operator_grid2mesh( grid, region%mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 15)

  END SUBROUTINE initialise_model_square_grid

END MODULE UFEMISM_main_model
