module UFEMISM_main_model
  ! The main regional ice-sheet model

! ===== Uses =====
! ================

  use mpi
  use configuration_module,          only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parallel_module,               only : par, sync, ierr, cerr, partition_list
  use data_types_module,             only : type_model_region, type_grid, type_remapping_mesh_mesh, &
                                            type_climate_matrix_global, type_ocean_matrix_global
  use reference_fields_module,       only : initialise_reference_geometries, map_reference_geometries_to_mesh
  use mesh_memory_module,            only : deallocate_mesh_all
  use mesh_creation_module,          only : create_mesh_from_cart_data
  use mesh_mapping_module,           only : calc_remapping_operators_mesh_mesh, deallocate_remapping_operators_mesh_mesh, &
                                            calc_remapping_operator_mesh2grid, deallocate_remapping_operators_mesh2grid, &
                                            calc_remapping_operator_grid2mesh, deallocate_remapping_operators_grid2mesh
  use mesh_update_module,            only : determine_mesh_fitness, create_new_mesh
  use mesh_single_module,            only : create_new_mesh_single, create_single_mesh_from_cart_data
  use netcdf_module,                 only : initialise_debug_fields, create_output_files, associate_debug_fields, &
                                            write_to_output_files, create_debug_file, reallocate_debug_fields
  use ice_dynamics_module,           only : initialise_ice_model, remap_ice_model, run_ice_model, update_ice_thickness
  use climate_module,                only : initialise_climate_model_regional, remap_climate_model, run_climate_model
  use ocean_module,                  only : initialise_ocean_model_regional, remap_ocean_model, run_ocean_model
  use BMB_module,                    only : initialise_BMB_model, remap_bmb_model, run_BMB_model
  use SMB_module,                    only : initialise_SMB_model, remap_smb_model, run_SMB_model
  use general_ice_model_data_module, only : initialise_mask_noice, initialise_basins
  use utilities_module,              only : time_display, inverse_oblique_sg_projection
  use thermodynamics_module,         only : initialise_thermo_model, run_thermo_model
  use general_sea_level_module,      only : calculate_PD_sealevel_contribution, calculate_icesheet_volume_and_area
  use scalar_data_output_module,     only : initialise_regional_scalar_data, write_regional_scalar_data
  use forcing_module,                only : forcing, update_sealevel_at_model_time

! ===== Preamble =====
! ====================

  implicit none

contains

! ===== Main routines =====
! =========================

  subroutine run_model( region, climate_matrix_global, t_end)
    ! Run the model region until the next coupling time

    implicit none

    ! In/output variables:
    type(type_model_region),          intent(inout) :: region
    type(type_climate_matrix_global), intent(inout) :: climate_matrix_global
    real(dp),                         intent(in)    :: t_end

    ! Local variables:
    character(len=256)                              :: routine_name
    real(dp)                                        :: meshfitness
    integer                                         :: it
    real(dp)                                        :: tstart, tstop, t1, t2, dt_ave

    ! === Initialisation ===
    ! ======================

    ! Set routine's name
    routine_name = 'run_model('  //  region%name  //  ')'

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise coupling-interval timer
    tstart = MPI_WTIME()

    ! Screen message
    if (par%master) then
      write(*,"(A)") ''
      write(*,"(5A,F9.3,A,F9.3,A)") &
            ' Running model region ', region%name, ' (', trim(region%long_name), &
            ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'
    end if
    call sync

    ! Set the intermediary pointers in "debug" to this region's debug data fields
    call associate_debug_fields(  region)

    ! Computation time tracking
    region%tcomp_total   = 0._dp
    region%tcomp_ice     = 0._dp
    region%tcomp_thermo  = 0._dp
    region%tcomp_climate = 0._dp
    region%tcomp_GIA     = 0._dp
    region%tcomp_mesh    = 0._dp

    ! Initialise sub-model timers
    t1 = 0._dp
    t2 = 0._dp

    ! Initialise iteration counter
    it = 0
    ! Initialise averaged time step
    dt_ave = 0._dp

    ! === The main model time loop ===
    ! ================================

    do while (region%time < t_end)

      ! Update iteration counter
      it = it + 1

      ! Sea-level
      ! =========

      if (C%choice_sealevel_model == 'prescribed') then
        ! Update global sea level based on record
        call update_sealevel_at_model_time( region%time)
        ! Update regional sea level based on record
        region%ice%SL_a( region%mesh%vi1:region%mesh%vi2) = forcing%sealevel_obs
      else
        ! Updated during coupling interval, or update not needed
      end if

      ! GIA
      ! ===

      ! Set timer
      t1 = MPI_WTIME()

      ! Pick the GIA model
      select case (C%choice_GIA_model)
        case ('none')
          ! Nothing to be done
        case ('ELRA')
          ! ELRA model
          call crash('ELRA model not implemented yet!')
        case ('SELEN')
          ! SELEN model
          call crash('SELEN model not implemented yet!')
        case default
          ! Unknown case
          call crash('unknown choice_GIA_model "' // &
                      trim(C%choice_GIA_model) // '"!')
      end select

      ! Record computation time
      t2 = MPI_WTIME()
      region%tcomp_GIA = region%tcomp_GIA + t2 - t1

      ! Mesh update
      ! ===========

      ! Set timer
      t1 = MPI_WTIME()

      !Default value
      meshfitness = 1._dp

      ! Check if the mesh needs to be updated
      if (region%time > region%t_last_mesh + C%dt_mesh_min) then
        call determine_mesh_fitness(region%mesh, region%ice, meshfitness)
      end if

      ! If required, update the mesh
      if (meshfitness < C%mesh_fitness_threshold) then
        region%t_last_mesh = region%time
        call run_model_update_mesh( region, climate_matrix_global)
      end if

      ! Record computation time
      t2 = MPI_WTIME()
      region%tcomp_mesh = region%tcomp_mesh + t2 - t1

      ! Ice dynamics
      ! ============

      ! Set timer
      t1 = MPI_WTIME()

      ! Calculate ice velocities and the resulting change in ice geometry
      ! NOTE: geometry is not updated yet; this happens at the end of the time loop
      call run_ice_model( region, t_end)

      ! Record computation time
      t2 = MPI_WTIME()
      region%tcomp_ice = region%tcomp_ice + t2 - t1

      ! Time display
      ! ============

      if (par%master .and. C%do_time_display) then
        call time_display(region, t_end, dt_ave, it)
      end if
      call sync

      ! Climate, ocean, SMB and BMB
      ! ===========================

      ! Set timer
      t1 = MPI_WTIME()

      ! Run the climate model
      if (region%do_climate) then
        call run_climate_model( region, climate_matrix_global, region%time)
      end if

      ! Run the ocean model
      if (region%do_ocean) then
        call run_ocean_model( region%mesh, region%ocean_matrix)
      end if

      ! Run the SMB model
      if (region%do_SMB) then
        call run_SMB_model( region%mesh, region%ice, region%climate_matrix, region%time, region%SMB, region%mask_noice)
      end if

      ! Run the BMB model
      if (region%do_BMB) then
        call run_BMB_model( region%mesh, region%ice, region%ocean_matrix%applied, region%BMB)
      end if

      ! Record computation time
      t2 = MPI_WTIME()
      region%tcomp_climate = region%tcomp_climate + t2 - t1

      ! Thermodynamics
      ! ==============

      ! Set timer
      t1 = MPI_WTIME()

      ! Run the thermodynamics model
      CALL run_thermo_model( region%mesh, region%ice, region%climate_matrix%applied, region%ocean_matrix%applied, region%SMB, region%time, region%do_thermo)

      ! Record computation time
      t2 = MPI_WTIME()
      region%tcomp_thermo = region%tcomp_thermo + t2 - t1

      ! Output
      ! ======

      ! Write NetCDF output
      if (region%do_output) then
        ! Check if the mesh has been updated
        if (.not. region%output_file_exists) then
          ! Create a new NetCDF file
          call create_output_files( region)
        end if
        ! Write to regional NetCDF output files
        call write_to_output_files( region)
      end if

      ! Update ice sheet area and volume
      call calculate_icesheet_volume_and_area( region)
      ! Write to regional scalar output
      call write_regional_scalar_data( region, region%time)

      ! == Update ice geometry
      ! ======================

      call update_ice_thickness( region%mesh, region%ice, region%refgeo_PD)

      ! Advance region time
      ! ===================

      region%time = region%time + region%dt
      dt_ave = dt_ave + region%dt
      call sync

    end do ! while (region%time < t_end)

    ! ===========================================
    ! ===== End of the main model time loop =====
    ! ===========================================

    ! Determine total ice sheet area, volume, volume-above-flotation
    ! and GMSL contribution, used for global scalar output
    call calculate_icesheet_volume_and_area( region)

    ! Write to NetCDF output one last time at the end of the simulation
    if (region%time == C%end_time_of_run) then
      ! If the mesh has been updated, create a new NetCDF file
      if (.not. region%output_file_exists) then
        call create_output_files( region)
      end if
      ! Write to regional 2-/3-D output
      call write_to_output_files( region)
      ! Write to regional scalar output
      call write_regional_scalar_data( region, region%time)
    end if

    ! Record coupling-interval computation time
    tstop = MPI_WTIME()
    region%tcomp_total = tstop - tstart

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_model

  subroutine run_model_update_mesh( region, climate_matrix_global)
    ! Perform a mesh update: create a new mesh based on the current modelled ice-sheet
    ! geometry, map all the model data from the old to the new mesh. Deallocate the
    ! old mesh, update the output files and the square grid maps.

    implicit none

    ! In- and output variables
    type(type_model_region),          intent(inout) :: region
    type(type_climate_matrix_global), intent(inout) :: climate_matrix_global

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'run_model_update_mesh'
    type(type_remapping_mesh_mesh)                :: map

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
    call sync

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

    ! Recalculate the "no ice" mask
    call initialise_mask_noice( region, region%mesh_new)

    ! Redefine the ice basins
    call initialise_basins( region%mesh, region%ice)

    ! Map reference geometries from the square grids to the mesh
    call map_reference_geometries_to_mesh( region, region%mesh_new)

    ! Remap de GIA submodel
    select case (C%choice_GIA_model)
      case ('none')
        ! Nothing to do
      case ('ELRA')
        ! ELRA model
        call crash('choice_GIA_model "' // trim(C%choice_GIA_model) // '" not implemented yet!')
      case default
        ! Unknown case
        call crash('unknown choice_GIA_model "' // trim(C%choice_GIA_model) // '"!')
    end select

    ! Remap all other submodels
    call remap_ice_model( region%mesh, region%mesh_new, map, region%ice, region%refgeo_PD, region%time)
    call remap_climate_model( region%mesh_new, region%climate_matrix, climate_matrix_global)
    CALL remap_ocean_model( region%mesh_new, region%ocean_matrix)
    call remap_SMB_model( region%mesh, region%mesh_new, map, region%SMB)
    call remap_BMB_model( region%mesh, region%mesh_new, map, region%BMB)

    ! Deallocate shared memory for the mapping arrays
    call deallocate_remapping_operators_mesh_mesh( map)

    ! Deallocate the old mesh, bind the region%mesh pointers to the new mesh.
    !call deallocate_mesh_all( region%mesh) TOdo make this automatic, fix it
    region%mesh = region%mesh_new

    ! When the next output is written, new output files must be created.
    region%output_file_exists = .false.

    ! Reallocate the debug fields for the new mesh, create a new debug file
    call reallocate_debug_fields( region)
    call create_debug_file(       region)

    ! Recompute the PD sea level contribution on the new mesh
    call calculate_PD_sealevel_contribution( region)

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
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_model_update_mesh

  subroutine initialise_model( region, name, climate_matrix_global, ocean_matrix_global)
    ! Initialise the entire model region - read initial and PD data, create the mesh,
    ! initialise the ice dynamics, climate, ocean, and SMB sub models

    implicit none

    ! In/output variables:
    type(type_model_region),          intent(inout) :: region
    character(len=3),                 intent(in)    :: name
    type(type_climate_matrix_global), intent(inout) :: climate_matrix_global
    type(type_ocean_matrix_global),   intent(inout) :: ocean_matrix_global

    ! Local variables:
    character(len=256)                              :: routine_name

    ! Add routine to path
    routine_name = 'initialise_model('  //  name  //  ')'
    call init_routine( routine_name)

    ! ===== Basic initialisation =====
    ! ================================

    ! Region name
    region%name = name

    if (region%name == 'NAM') then
      region%long_name = 'North America'
    else if (region%name == 'EAS') then
      region%long_name = 'Eurasia'
    else if (region%name == 'GRL') then
      region%long_name = 'Greenland'
    else if (region%name == 'ANT') then
      region%long_name = 'Antarctica'
    end if

    if (par%master) then
      write(*,"(A)") ''
      write(*,"(5A)") ' Initialising model region ', region%name, &
                      ' (', trim(region%long_name), ')...'
    end if
    call sync

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
    call sync

    call initialise_model_square_grid( region, region%grid_output, C%dx_grid_output)
    call initialise_model_square_grid( region, region%grid_GIA,    C%dx_grid_GIA   )
    call initialise_model_square_grid( region, region%grid_smooth, C%dx_grid_smooth)

    ! ===== Dummy fields for debugging =====
    ! ======================================

    call initialise_debug_fields( region)
    call associate_debug_fields( region)

    ! ===== The "no ice" mask =====
    ! =============================

    call initialise_mask_noice(region, region%mesh)

    ! ===== The ice dynamics model =====
    ! ==================================

    call initialise_ice_model( region%mesh, region%ice, region%refgeo_init, region%refgeo_PD)

    ! ===== Ice basins =====
    ! ======================

    call initialise_basins( region%mesh, region%ice)

    ! ===== The climate model =====
    ! =============================

    call initialise_climate_model_regional( region, climate_matrix_global)
    call run_climate_model( region, climate_matrix_global, C%start_time_of_run)

    ! ===== The ocean model =====
    ! ===========================

    call initialise_ocean_model_regional( region, ocean_matrix_global)

    ! ===== The SMB model =====
    ! =========================

    call initialise_SMB_model( region%mesh, region%ice, region%SMB, region%name)
    call run_SMB_model( region%mesh, region%ice, region%climate_matrix, C%start_time_of_run, region%SMB, region%mask_noice)

    ! ===== The BMB model =====
    ! =========================

    call initialise_BMB_model( region%mesh, region%ice, region%BMB, region%name)

    ! ===== The GIA model =====
    ! =========================

    select case (C%choice_GIA_model)

      case ('none')
        ! Nothing to do

      case ('ELRA')
        ! ELRA model
        call crash('choice_GIA_model "' // trim(C%choice_GIA_model) // '" not implemented yet!')

      case default
        ! Unknown case
        call crash('unknown choice_GIA_model "' // trim(C%choice_GIA_model) // '"!')

    end select

    ! ===== Ice temperature =====
    ! ===========================

    ! Initialise the ice temperature profile
    call initialise_thermo_model( region%mesh, region%ice, region%climate_matrix%applied, region%ocean_matrix%applied, region%SMB, region%name)

    ! ===== Scalar ice data =====
    ! ===========================

    ! Calculate ice sheet metadata (volume, area, GMSL contribution),
    ! for writing to the first time point of the output file
    call calculate_PD_sealevel_contribution( region)
    call calculate_icesheet_volume_and_area( region)

    ! ===== Output files =====
    ! ========================

    ! Create output file for regional scalar data
    call initialise_regional_scalar_data( region)

    ! Write regional scalar data at time t=0
    CALL write_regional_scalar_data( region, C%start_time_of_run)

    ! Create output file for regional 2-/3-D data
    call create_output_files( region)

    ! Write 2-/3-D data at time t=0
    call write_to_output_files( region)

    ! ===== Finalisation =====
    ! ========================

    if (par%master) then
      write(*,"(3A)") ' Finished initialising model region ', region%name, '.'
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = huge( 1))

  end subroutine initialise_model

! ===== Auxiliary routines =====
! ==============================

  subroutine allocate_region_timers_and_scalars( region)
    ! Allocate shared memory for this region's timers (used for the asynchronous coupling between the
    ! ice dynamics and the secondary model components), and for the scalars (integrated ice volume and
    ! area, SMB components, computation times, etc.)

    implicit none

    ! In/output variables:
    type(type_model_region),      intent(inout) :: region

    ! Local variables:
    character(len=256), parameter               :: routine_name = 'allocate_region_timers_and_scalars'

    ! Add routine to path
    call init_routine( routine_name)

    ! Timers and time steps
    ! =====================

    region%time           = C%start_time_of_run
    region%dt             = C%dt_min
    region%dt_prev        = C%dt_min

    region%t_last_mesh    = C%start_time_of_run
    region%t_next_mesh    = C%start_time_of_run + C%dt_mesh_min
    region%do_mesh        = .false.

    region%t_last_SIA     = C%start_time_of_run
    region%t_next_SIA     = C%start_time_of_run
    region%do_SIA         = .true.

    region%t_last_SSA     = C%start_time_of_run
    region%t_next_SSA     = C%start_time_of_run
    region%do_SSA         = .true.

    region%t_last_DIVA    = C%start_time_of_run
    region%t_next_DIVA    = C%start_time_of_run
    region%do_DIVA        = .true.

    region%t_last_thermo  = C%start_time_of_run
    region%t_next_thermo  = C%start_time_of_run + C%dt_thermo
    region%do_thermo      = .false.

    region%t_last_climate = C%start_time_of_run
    region%t_next_climate = C%start_time_of_run
    region%do_climate     = .true.

    region%t_last_ocean   = C%start_time_of_run
    region%t_next_ocean   = C%start_time_of_run
    region%do_ocean       = .true.

    region%t_last_SMB     = C%start_time_of_run
    region%t_next_SMB     = C%start_time_of_run
    region%do_SMB         = .true.

    region%t_last_BMB     = C%start_time_of_run
    region%t_next_BMB     = C%start_time_of_run
    region%do_BMB         = .true.

    region%t_last_ELRA    = C%start_time_of_run
    region%t_next_ELRA    = C%start_time_of_run
    region%do_ELRA        = .true.

    region%t_last_output  = C%start_time_of_run
    region%t_next_output  = C%start_time_of_run
    region%do_output      = .true.

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = 65)

  end subroutine allocate_region_timers_and_scalars

  subroutine initialise_model_square_grid( region, grid, dx)
    ! Initialise a regular square grid enveloping this model region

    implicit none

    ! In/output variables:
    type(type_model_region),       intent(inout) :: region
    type(type_grid),               intent(inout) :: grid
    real(dp),                      intent(in)    :: dx

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'initialise_model_square_grid'
    real(dp)                                     :: xmid, ymid
    integer                                      :: nsx, nsy, i, j, n
    real(dp), parameter                          :: tol = 1E-9_dp

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Basic grid data ===
    ! =======================

    ! Projection parameters for this region
    select case (region%name)

      case ('NAM')
        ! North America
        grid%lambda_M     = C%lambda_M_NAM
        grid%phi_M        = C%phi_M_NAM
        grid%alpha_stereo = C%alpha_stereo_NAM

      case ('EAS')
        ! Eurasia
        grid%lambda_M     = C%lambda_M_EAS
        grid%phi_M        = C%phi_M_EAS
        grid%alpha_stereo = C%alpha_stereo_EAS

      case ('GRL')
        ! Greenland
        grid%lambda_M     = C%lambda_M_GRL
        grid%phi_M        = C%phi_M_GRL
        grid%alpha_stereo = C%alpha_stereo_GRL

      case ('ANT')
        ! Antarctica
        grid%lambda_M     = C%lambda_M_ANT
        grid%phi_M        = C%phi_M_ANT
        grid%alpha_stereo = C%alpha_stereo_ANT

    end select

    nsx = 0
    nsy = 0

    ! Resolution
    grid%dx = dx

    ! Determine the center of the model domain
    xmid = (region%mesh%xmin + region%mesh%xmax) / 2._dp
    ymid = (region%mesh%ymin + region%mesh%ymax) / 2._dp

    ! Number of points at each side of domain center
    nsx = floor((region%mesh%xmax - xmid) / grid%dx)
    nsy = floor((region%mesh%ymax - ymid) / grid%dx)

    ! Determine total number of points per dimension as twice the
    ! number per side, plus the middle point, plus 1 grid cell to
    ! make sure the mesh lies completely inside the grid for any
    ! grid resolution
    grid%nx = (2*nsx + 1) + 1
    grid%ny = (2*nsy + 1) + 1

    ! Determine total number of grid points
    grid%n  = grid%nx * grid%ny

    ! Assign range to each processor
    call partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    call partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

    ! === Grid generation ===
    ! =======================

    ! Allocate x and y coordinates for square grid
    allocate( grid%x ( grid%nx ))
    allocate( grid%y ( grid%ny ))

    ! Compute corners of grid
    grid%xmin = grid%x(1      )
    grid%xmax = grid%x(grid%nx)
    grid%ymin = grid%y(1      )
    grid%ymax = grid%y(grid%ny)

    ! Compute x coordinates
    grid%xmin = xmid - nsx * grid%dx
    grid%xmax = xmid + nsx * grid%dx
    do i = 1, grid%nx
      grid%x( i) = grid%xmin + (i-1)*grid%dx
    end do

    ! Compute y coordinates
    grid%ymin = ymid - nsy * grid%dx
    grid%ymax = ymid + nsy * grid%dx
    do j = 1, grid%ny
      grid%y( j) = grid%ymin + (j-1)*grid%dx
    end do

    ! === Grid-to-vector tables ===
    ! =============================

    ! Allocate table data
    allocate ( grid%ij2n( grid%nx, grid%ny ))
    allocate ( grid%n2ij( grid%n , 2       ))

    ! Tolerance; points lying within this distance of each other are treated as identical
    grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

    ! Set up grid-to-vector translation tables
    n = 0
    do i = 1, grid%nx
    do j = 1, grid%ny
      n = n+1
      grid%ij2n( i,j) = n
      grid%n2ij( n,:) = [i,j]
    end do
    end do

    ! === Mappings between mesh and grid ===
    ! ======================================

    ! Calculate mapping arrays between the mesh and the grid
    call calc_remapping_operator_mesh2grid( region%mesh, grid)
    call calc_remapping_operator_grid2mesh( grid, region%mesh)

    ! === Geographical coordinates ===
    ! ================================

    ! Calculate lat-lon coordinates
    allocate( grid%lat(grid%nx, grid%ny))
    allocate( grid%lon(grid%nx, grid%ny))

    do i = 1, grid%nx
    do j = 1, grid%ny
      call inverse_oblique_sg_projection( grid%x( i), grid%y( j), region%mesh%lambda_M, &
                                          region%mesh%phi_M, region%mesh%alpha_stereo, &
                                          grid%lon( i,j), grid%lat( i,j))
    end do
    end do

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_model_square_grid

end module UFEMISM_main_model
