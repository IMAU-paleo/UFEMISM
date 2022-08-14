MODULE UFEMISM_main_model
  ! The main regional ice-sheet model

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

  USE mpi
  USE configuration_module,                ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                        ONLY: perr
  USE parallel_module,                     ONLY: par, sync, ierr, cerr, partition_list, &
                                                 allocate_shared_int_0D,   allocate_shared_dp_0D, &
                                                 allocate_shared_int_1D,   allocate_shared_dp_1D, &
                                                 allocate_shared_int_2D,   allocate_shared_dp_2D, &
                                                 allocate_shared_int_3D,   allocate_shared_dp_3D, &
                                                 allocate_shared_bool_0D,  allocate_shared_bool_1D, &
                                                 reallocate_shared_int_0D, reallocate_shared_dp_0D, &
                                                 reallocate_shared_int_1D, reallocate_shared_dp_1D, &
                                                 reallocate_shared_int_2D, reallocate_shared_dp_2D, &
                                                 reallocate_shared_int_3D, reallocate_shared_dp_3D, &
                                                 deallocate_shared
  USE netcdf_module,                       ONLY: debug, write_to_debug_file, create_output_files, write_to_output_files, &
                                                 initialise_debug_fields, associate_debug_fields, reallocate_debug_fields, &
                                                 create_debug_file, write_PETSc_matrix_to_NetCDF, create_regional_scalar_output_file
  USE data_types_module,                   ONLY: type_model_region, type_mesh, type_grid, type_remapping_mesh_mesh, &
                                                 type_climate_matrix_global, type_ocean_matrix_global
  USE reference_fields_module,             ONLY: initialise_reference_geometries, map_reference_geometries_to_mesh
  USE mesh_memory_module,                  ONLY: deallocate_mesh_all
  USE mesh_help_functions_module,          ONLY: inverse_oblique_sg_projection
  USE mesh_creation_module,                ONLY: create_mesh_from_cart_data
  USE mesh_mapping_module,                 ONLY: calc_remapping_operators_mesh_mesh, deallocate_remapping_operators_mesh_mesh, &
                                                 calc_remapping_operator_mesh2grid, deallocate_remapping_operators_mesh2grid, &
                                                 calc_remapping_operator_grid2mesh, deallocate_remapping_operators_grid2mesh
  USE mesh_update_module,                  ONLY: determine_mesh_fitness, create_new_mesh
  USE general_ice_model_data_module,       ONLY: initialise_mask_noice, initialise_basins
  USE forcing_module,                      ONLY: forcing, update_sealevel_record_at_model_time
  USE ice_dynamics_module,                 ONLY: initialise_ice_model,                    remap_ice_model,      run_ice_model,      update_ice_thickness
  USE thermodynamics_module,               ONLY: initialise_ice_temperature,                                    run_thermo_model,   calc_ice_rheology
  USE climate_module,                      ONLY: initialise_climate_model_regional,       remap_climate_model,  run_climate_model
  USE ocean_module,                        ONLY: initialise_ocean_model_regional,         remap_ocean_model,    run_ocean_model
  USE SMB_module,                          ONLY: initialise_SMB_model,                    remap_SMB_model,      run_SMB_model,      SMB_IMAUITM_inversion
  USE BMB_module,                          ONLY: initialise_BMB_model,                    remap_BMB_model,      run_BMB_model
  USE isotopes_module,                     ONLY: initialise_isotopes_model,               remap_isotopes_model, run_isotopes_model, calculate_reference_isotopes
  USE bedrock_ELRA_module,                 ONLY: initialise_ELRA_model,                   remap_ELRA_model,     run_ELRA_model

  USE tests_and_checks_module,             ONLY: run_all_matrix_tests
  USE basal_conditions_and_sliding_module, ONLY: basal_sliding_inversion
  USE restart_module,                      ONLY: read_mesh_from_restart_file, read_init_data_from_restart_file, remap_restart_data
  USE general_sea_level_module,            ONLY: calculate_PD_sealevel_contribution, calculate_icesheet_volume_and_area
  USE ice_velocity_module,                 ONLY: solve_DIVA
  USE utilities_module,                    ONLY: time_display
  USE scalar_data_output_module,           ONLY: write_regional_scalar_data

# if (defined(DO_SELEN))
  USE SELEN_main_module,                   ONLY: apply_SELEN_bed_geoid_deformation_rates, remap_SELEN_model
# endif

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  ! Run the model region until the next coupling time
  SUBROUTINE run_model( region, climate_matrix_global, t_end)
    ! Run the model until t_end (usually a 100 years further)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),          INTENT(INOUT)  :: region
    TYPE(type_climate_matrix_global), INTENT(INOUT)  :: climate_matrix_global
    REAL(dp),                         INTENT(IN)     :: t_end

    ! Local variables:
    CHARACTER(LEN=256)                               :: routine_name
    REAL(dp)                                         :: meshfitness
    INTEGER                                          :: it
    REAL(dp)                                         :: tstart, tstop, t1, t2, dt_ave
    CHARACTER(LEN=9)                                 :: r_time, r_step, r_adv, r_ave

    ! Add routine to path
    routine_name = 'run_model('  //  region%name  //  ')'
    CALL init_routine( routine_name)

    IF (par%master) THEN

      WRITE(0,*) ''
      WRITE (0,'(5A,F9.3,A,F9.3,A)') '  Running model region ', region%name, ' (', TRIM(region%long_name), &
                                     ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'
    END IF

    ! Set the intermediary pointers in "debug" to this region's debug data fields
    CALL associate_debug_fields(  region)

    ! Computation time tracking
    region%tcomp_total          = 0._dp
    region%tcomp_ice            = 0._dp
    region%tcomp_thermo         = 0._dp
    region%tcomp_climate        = 0._dp
    region%tcomp_GIA            = 0._dp
    region%tcomp_mesh           = 0._dp

    ! Initialise coupling-interval timer
    tstart = MPI_WTIME()

    ! Initialise sub-model timers
    t1 = MPI_WTIME()
    t2 = 0._dp

    ! Initialise iteration counter
    it = 0
    ! Initialise averaged time step
    dt_ave = 0._dp

    ! ====================================
    ! ===== The main model time loop =====
    ! ====================================

    DO WHILE (region%time < t_end)

      ! Update iteration counter
      it = it + 1

      ! == Sea-level
      ! ============

      IF (C%choice_sealevel_model == 'prescribed') THEN
        ! Update global sea level based on record
        CALL update_sealevel_record_at_model_time( region%time)
        ! Update regional sea level based on record
        region%ice%SL_a( region%mesh%vi1:region%mesh%vi2) = forcing%sealevel_obs
      ELSE
        ! Updated during coupling interval or update not needed
      END IF

      ! == GIA
      ! ======

      t1 = MPI_WTIME()
      IF     (C%choice_GIA_model == 'none') THEN
        ! Nothing to be done
      ELSEIF (C%choice_GIA_model == 'ELRA') THEN
        CALL run_ELRA_model( region)
      ELSEIF (C%choice_GIA_model == 'SELEN') THEN
#       if (defined(DO_SELEN))
        CALL apply_SELEN_bed_geoid_deformation_rates( region)
#       endif
      ELSE
        CALL crash('unknown choice_GIA_model "' // TRIM(C%choice_GIA_model) // '"!')
      END IF
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_GIA = region%tcomp_GIA + t2 - t1

      ! == Mesh update
      ! ==============

      ! Check if the mesh needs to be updated
      IF (par%master) THEN
        t2 = MPI_WTIME()
      END IF

      meshfitness = 1._dp

      IF (region%time > region%t_last_mesh + C%dt_mesh_min) THEN
        CALL determine_mesh_fitness(region%mesh, region%ice, meshfitness)
      END IF

      IF (par%master) THEN
        region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2
      END IF

      ! If needed (or commanded), update the mesh
      IF ( meshfitness < C%mesh_fitness_threshold .OR. &
          (region%do_mesh .AND. region%time >= C%start_time_of_run + C%do_force_mesh_update_after) ) THEN

        region%t_last_mesh = region%time

        IF (par%master) THEN
          t2 = MPI_WTIME()
        END IF

        CALL run_model_update_mesh( region, climate_matrix_global)

        IF (par%master) THEN
          region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2
        END IF

        ! Make sure a forced mesh update is not triggered again next time step
        region%do_mesh = .FALSE.

      END IF

      ! == Ice dynamics
      ! ===============

      ! Calculate ice velocities and the resulting change in ice geometry
      ! NOTE: geometry is not updated yet; this happens at the end of the time loop
      t1 = MPI_WTIME()
      CALL run_ice_model( region, t_end)
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_ice = region%tcomp_ice + t2 - t1

      ! == Time display
      ! ===============

      if (par%master .AND. C%do_time_display) then
        call time_display(region, t_end, dt_ave, it)
      end if

      ! == Climate, ocean, SMB and BMB
      ! ==============================

      t1 = MPI_WTIME()

      ! Run the climate model
      IF (region%do_climate) THEN
        CALL run_climate_model( region, climate_matrix_global, region%time)
      END IF

      ! Run the ocean model
      IF (region%do_ocean) THEN
        CALL run_ocean_model( region%mesh, region%grid_smooth, region%ice, region%ocean_matrix, region%climate_matrix, region%name, region%time)
      END IF

      ! Run the SMB model
      IF (region%do_SMB) THEN
        CALL run_SMB_model( region%mesh, region%ice, region%climate_matrix, region%time, region%SMB, region%mask_noice)
      END IF

      ! Run the BMB model
      IF (region%do_BMB) THEN
        CALL run_BMB_model( region%mesh, region%ice, region%ocean_matrix%applied, region%BMB, region%name, region%time, region%refgeo_PD)
      END IF

      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_climate = region%tcomp_climate + t2 - t1

      ! == Thermodynamics
      ! =================

      t1 = MPI_WTIME()
      CALL run_thermo_model( region%mesh, region%ice, region%climate_matrix%applied, region%ocean_matrix%applied, region%SMB, region%time, do_solve_heat_equation = region%do_thermo)
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_thermo = region%tcomp_thermo + t2 - t1

      ! == Isotopes
      ! ===========

      CALL run_isotopes_model( region)

      ! == Basal sliding inversion
      ! ==========================

      IF (C%do_basal_sliding_inversion .AND. region%do_basal) THEN
        ! Adjust bed roughness
        CALL basal_sliding_inversion( region%mesh, region%grid_smooth, region%ice, region%refgeo_PD, region%time)
      END IF

      ! == SBM IMAU-ITM inversion
      ! =========================

      IF (C%do_SMB_IMAUITM_inversion .AND. region%do_SMB_inv) THEN
        CALL SMB_IMAUITM_inversion( region%mesh, region%ice, region%climate_matrix%applied, region%SMB, region%refgeo_PD)
      END IF

      ! == Output
      ! =========

      ! Write NetCDF output
      IF (region%do_output) THEN
        ! If the mesh has been updated, create a new NetCDF file
        IF (.NOT. region%output_file_exists) THEN
          CALL create_output_files( region)
          CALL sync
          region%output_file_exists = .TRUE.
        END IF
        ! ! Update ice sheet area and volume
        ! CALL calculate_icesheet_volume_and_area( region)
        ! ! Write to regional scalar output
        ! CALL write_regional_scalar_data( region, region%time)
        ! Write to regional 2-/3-D output
        CALL write_to_output_files( region)
      END IF

      ! Update ice sheet area and volume
      CALL calculate_icesheet_volume_and_area( region)
      ! Write to regional scalar output
      CALL write_regional_scalar_data( region, region%time)

      ! == Update ice geometry
      ! ======================

      CALL update_ice_thickness( region%mesh, region%ice, region%mask_noice, region%refgeo_PD, region%refgeo_GIAeq)
      CALL sync

      ! == Advance region time
      ! ======================

      IF (par%master) region%time = region%time + region%dt
      IF (par%master) dt_ave = dt_ave + region%dt
      CALL sync

    END DO

    ! ===========================================
    ! ===== End of the main model time loop =====
    ! ===========================================

    ! Determine total ice sheet area, volume, volume-above-flotation
    ! and GMSL contribution, used for global scalar output
    CALL calculate_icesheet_volume_and_area( region)

    ! Write to NetCDF output one last time at the end of the simulation
    IF (region%time == C%end_time_of_run) THEN
      ! Update output time of record
      region%t_last_output = C%end_time_of_run
      ! If the mesh has been updated, create a new NetCDF file
      IF (.NOT. region%output_file_exists) THEN
        CALL create_output_files( region)
        CALL sync
        region%output_file_exists = .TRUE.
      END IF
      ! Write to regional 2-/3-D output
      CALL write_to_output_files( region)
      ! Write to regional scalar output
      CALL write_regional_scalar_data( region, region%time)
    END IF

    tstop = MPI_WTIME()
    region%tcomp_total = tstop - tstart

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model

  ! Update the mesh and everything else that needs it
  SUBROUTINE run_model_update_mesh( region, climate_matrix_global)
    ! Perform a mesh update: create a new mesh based on the current modelled ice-sheet
    ! geometry, map all the model data from the old to the new mesh. Deallocate the
    ! old mesh, update the output files and the square grid maps.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_climate_matrix_global),    INTENT(INOUT) :: climate_matrix_global

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_model_update_mesh'
    TYPE(type_remapping_mesh_mesh)                     :: map

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new mesh
    CALL create_new_mesh( region)

    IF (par%master) WRITE(0,*) '  Reallocating and remapping after mesh update...'

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

    ! Recalculate the "no ice" mask (also needed for the climate model)
    CALL reallocate_shared_int_1D( region%mesh_new%nV, region%mask_noice, region%wmask_noice)
    CALL initialise_mask_noice( region, region%mesh_new)

    ! Redefine the ice basins
    CALL reallocate_shared_int_1D( region%mesh_new%nV, region%ice%basin_ID, region%ice%wbasin_ID)
    CALL initialise_basins( region%mesh_new, region%ice%basin_ID, region%ice%nbasins, region%name)

    ! Reallocate memory for reference geometries
    CALL deallocate_shared( region%refgeo_init%wHi )
    CALL deallocate_shared( region%refgeo_init%wHb )
    CALL deallocate_shared( region%refgeo_init%wHs )

    CALL deallocate_shared( region%refgeo_PD%wHi   )
    CALL deallocate_shared( region%refgeo_PD%wHb   )
    CALL deallocate_shared( region%refgeo_PD%wHs   )

    CALL deallocate_shared( region%refgeo_GIAeq%wHi)
    CALL deallocate_shared( region%refgeo_GIAeq%wHb)
    CALL deallocate_shared( region%refgeo_GIAeq%wHs)

    CALL allocate_shared_dp_1D( region%mesh_new%nV, region%refgeo_init%Hi , region%refgeo_init%wHi )
    CALL allocate_shared_dp_1D( region%mesh_new%nV, region%refgeo_init%Hb , region%refgeo_init%wHb )
    CALL allocate_shared_dp_1D( region%mesh_new%nV, region%refgeo_init%Hs , region%refgeo_init%wHs )

    CALL allocate_shared_dp_1D( region%mesh_new%nV, region%refgeo_PD%Hi   , region%refgeo_PD%wHi   )
    CALL allocate_shared_dp_1D( region%mesh_new%nV, region%refgeo_PD%Hb   , region%refgeo_PD%wHb   )
    CALL allocate_shared_dp_1D( region%mesh_new%nV, region%refgeo_PD%Hs   , region%refgeo_PD%wHs   )

    CALL allocate_shared_dp_1D( region%mesh_new%nV, region%refgeo_GIAeq%Hi, region%refgeo_GIAeq%wHi)
    CALL allocate_shared_dp_1D( region%mesh_new%nV, region%refgeo_GIAeq%Hb, region%refgeo_GIAeq%wHb)
    CALL allocate_shared_dp_1D( region%mesh_new%nV, region%refgeo_GIAeq%Hs, region%refgeo_GIAeq%wHs)

    ! Map reference geometries from the square grids to the mesh
    CALL map_reference_geometries_to_mesh( region, region%mesh_new)

    ! Remap the GIA submodel
    IF (C%choice_GIA_model == 'none') THEN
      ! Do nothing
    ELSEIF (C%choice_GIA_model == 'ELRA') THEN
      CALL remap_ELRA_model(   region%mesh, region%mesh_new, map, region%ice, region%refgeo_PD, region%grid_GIA)
    ELSEIF (C%choice_GIA_model == 'SELEN') THEN
#     if (defined(DO_SELEN))
      CALL remap_SELEN_model(  region%mesh_new, region%SELEN)
#     endif
    ELSE
      CALL crash('unknown choice_GIA_model "' // TRIM( C%choice_GIA_model) // '"!')
    END IF

    ! Remap all other submodels
    CALL remap_ice_model(      region%mesh, region%mesh_new, map, region%ice, region%refgeo_PD, region%time)
    CALL remap_climate_model(  region%mesh, region%mesh_new, map, region%climate_matrix, climate_matrix_global, region%refgeo_PD, region%grid_smooth, region%mask_noice, region%name, region%time)
    CALL remap_ocean_model(    region%mesh, region%mesh_new, map, region%ocean_matrix)
    CALL remap_SMB_model(      region%mesh, region%mesh_new, map, region%SMB)
    CALL remap_BMB_model(      region%mesh, region%mesh_new, map, region%BMB)
    CALL remap_isotopes_model( region%mesh, region%mesh_new, map, region)

    ! Deallocate shared memory for the mapping arrays
    CALL deallocate_remapping_operators_mesh_mesh( map)

    ! Deallocate the old mesh, bind the region%mesh pointers to the new mesh.
    CALL deallocate_mesh_all( region%mesh)
    region%mesh = region%mesh_new

    ! Run the sub-models once to fill them in
    CALL run_climate_model( region, climate_matrix_global, region%time)
    CALL run_ocean_model( region%mesh, region%grid_smooth, region%ice, region%ocean_matrix, region%climate_matrix, region%name, region%time)
    CALL run_SMB_model( region%mesh, region%ice, region%climate_matrix, region%time, region%SMB, region%mask_noice)
    CALL run_BMB_model( region%mesh, region%ice, region%ocean_matrix%applied, region%BMB, region%name, region%time, region%refgeo_PD)

    ! Remap key restart data
    IF (C%is_restart) THEN
      CALL remap_restart_data( region)
    END IF

    ! When the next output is written, new output files must be created.
    region%output_file_exists = .FALSE.

    ! Reallocate the debug fields for the new mesh, create a new debug file
    CALL reallocate_debug_fields( region)
    CALL create_debug_file(       region)

    ! Recalculate the reference precipitation isotope content (must be done after region%mesh has been cycled)
    CALL calculate_reference_isotopes( region)

    ! Recompute the PD sea level contribution on the new mesh
    CALL calculate_PD_sealevel_contribution( region)

    region%t_next_SIA     = region%time
    region%t_next_SSA     = region%time
    region%t_next_DIVA    = region%time
    region%t_next_thermo  = region%time
    region%t_next_climate = region%time
    region%t_next_ocean   = region%time
    region%t_next_SMB     = region%time
    region%t_next_BMB     = region%time
    region%t_next_ELRA    = region%time
    region%t_next_basal   = region%time
    region%t_next_SMB_inv = region%time

    region%do_SIA         = .TRUE.
    region%do_SSA         = .TRUE.
    region%do_DIVA        = .TRUE.
    region%do_thermo      = .TRUE.
    region%do_climate     = .TRUE.
    region%do_ocean       = .TRUE.
    region%do_SMB         = .TRUE.
    region%do_BMB         = .TRUE.
    region%do_ELRA        = .TRUE.
    region%do_basal       = .TRUE.
    region%do_SMB_inv     = .TRUE.

    IF (par%master) WRITE(0,*) '  Finished reallocating and remapping.'
    IF (par%master) WRITE(0,*) '  Running again now...'

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model_update_mesh

  ! Initialise the entire model region
  SUBROUTINE initialise_model( region, name, climate_matrix_global, ocean_matrix_global)
    ! Initialise the entire model region - read initial and PD data, create the mesh,
    ! and initialise the ice dynamics, climate, ocean, and SMB sub models

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),          INTENT(INOUT)  :: region
    CHARACTER(LEN=3),                 INTENT(IN)     :: name
    TYPE(type_climate_matrix_global), INTENT(INOUT)  :: climate_matrix_global
    TYPE(type_ocean_matrix_global),   INTENT(INOUT)  :: ocean_matrix_global

    ! Local variables:
    CHARACTER(LEN=256)                               :: routine_name

    ! Add routine to path
    routine_name = 'initialise_model('  //  name  //  ')'
    CALL init_routine( routine_name)

    ! ===== Region name =====
    ! =======================

    ! Region name
    region%name        = name
    IF      (region%name == 'NAM') THEN
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

    ! ===== Timers and scalars =====
    ! ==============================

    CALL allocate_region_timers_and_scalars( region)

    ! ===== Initial geometry =====
    ! ============================

    IF (C%is_restart) THEN

      ! Read mesh from a restart file
      CALL read_mesh_from_restart_file( region%mesh, region%restart, region%name, region%time)
      ! Read data (on the mesh) from a restart file
      CALL read_init_data_from_restart_file( region%restart, region%name)
      ! Initialise topographic data fields (on a square grid), mapping the initial topo from the mesh onto the grid
      CALL initialise_reference_geometries( region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq, region%name, region%mesh, region%restart)

    ELSE

      ! Initialise topographic data fields (on a square grid)
      CALL initialise_reference_geometries( region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq, region%name)
      ! Create a new mesh from the reference initial geometry
      CALL create_mesh_from_cart_data( region)

    END IF

    ! ===== Reference geometries: grid to mesh =====
    ! ==============================================

    IF (par%master) WRITE(0,*) '  Mapping reference geometries onto the initial model mesh...'

    ! Allocate memory for reference data on the mesh
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_init%Hi , region%refgeo_init%wHi )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_init%Hb , region%refgeo_init%wHb )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_init%Hs , region%refgeo_init%wHs )

    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_PD%Hi,    region%refgeo_PD%wHi   )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_PD%Hb,    region%refgeo_PD%wHb   )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_PD%Hs,    region%refgeo_PD%wHs   )

    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_GIAeq%Hi, region%refgeo_GIAeq%wHi)
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_GIAeq%Hb, region%refgeo_GIAeq%wHb)
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_GIAeq%Hs, region%refgeo_GIAeq%wHs)

    ! Map data from the square grids to the mesh
    CALL map_reference_geometries_to_mesh( region, region%mesh)

    IF (par%master) WRITE(0,*) '  Finished mapping reference geometries.'

    ! ===== Square grids =====
    ! ========================

    IF (par%master) WRITE(0,*) '  Initialising square grids for output, GIA, and data smoothing...'

    CALL initialise_model_square_grid( region, region%grid_output, C%dx_grid_output)
    CALL initialise_model_square_grid( region, region%grid_GIA,    C%dx_grid_GIA   )
    CALL initialise_model_square_grid( region, region%grid_smooth, C%dx_grid_smooth)

    ! ===== Debug fields =====
    ! ========================

    IF (par%master) WRITE(0,*) '  Initialising debug fields...'

    CALL initialise_debug_fields( region)
    CALL associate_debug_fields( region)

    ! ===== The "no ice" mask =====
    ! =============================

    CALL allocate_shared_int_1D( region%mesh%nV, region%mask_noice, region%wmask_noice)
    CALL initialise_mask_noice( region, region%mesh)

    ! ===== The ice dynamics model =====
    ! ==================================

    CALL initialise_ice_model( region%mesh, region%ice, region%refgeo_init, region%refgeo_PD, region%restart)

    ! ===== Ice basins =====
    ! ======================

    ! Allocate shared memory
    CALL allocate_shared_int_1D( region%mesh%nV, region%ice%basin_ID, region%ice%wbasin_ID)
    CALL allocate_shared_int_0D(                 region%ice%nbasins,  region%ice%wnbasins )

    ! Define basins
    CALL initialise_basins( region%mesh, region%ice%basin_ID, region%ice%nbasins, region%name)

    ! ===== The climate model =====
    ! =============================

    CALL initialise_climate_model_regional( region, climate_matrix_global)

    ! ===== The ocean model =====
    ! ===========================

    CALL initialise_ocean_model_regional( region, ocean_matrix_global)

    ! ===== The SMB model =====
    ! =========================

    CALL initialise_SMB_model( region%mesh, region%ice, region%SMB, region%name, region%restart)

    ! ===== The BMB model =====
    ! =========================

    CALL initialise_BMB_model( region%mesh, region%ice, region%BMB, region%name, region%restart)

    ! ===== The GIA model =====
    ! =========================

    IF     (C%choice_GIA_model == 'none') THEN
      ! Nothing to be done
    ELSEIF (C%choice_GIA_model == 'ELRA') THEN
      CALL initialise_ELRA_model( region%mesh, region%grid_GIA, region%ice, region%refgeo_PD)
    ELSEIF (C%choice_GIA_model == 'SELEN') THEN
      ! Nothing to be done; SELEN is initialised globally by UFEMISM_program.
    ELSE
      CALL crash('unknown choice_GIA_model "' // TRIM(C%choice_GIA_model) // '"!')
    END IF

    ! ===== The isotopes model =====
    ! ==============================

    CALL initialise_isotopes_model( region)

    ! ===== Geothermal heat flux =====
    ! ================================

    ! This (init regional GHF) is currently done in the ice dynamics module (initialise_ice_model).
    ! Might be good to move it to the forcing module later, a la IMAU-ICE.

    ! ===== Ice temperature profile =====
    ! ===================================

    ! Run the climate and SMB models once, to get the correct surface temperature+SMB fields for the ice temperature initialisation
    CALL run_climate_model( region, climate_matrix_global, C%start_time_of_run)
    CALL run_SMB_model( region%mesh, region%ice, region%climate_matrix, C%start_time_of_run, region%SMB, region%mask_noice)

    ! Initialise the ice temperature profile
    CALL initialise_ice_temperature( region%mesh, region%ice, region%climate_matrix%applied, region%ocean_matrix%applied, region%SMB, region%name, region%restart)

    ! Initialise the rheology
    CALL calc_ice_rheology( region%mesh, region%ice, C%start_time_of_run)

    ! ===== Exception: Initial velocities for choice_ice_dynamics == "none" =====
    ! ===========================================================================

    ! If we're running with choice_ice_dynamics == "none", calculate a velocity field
    ! once during initialisation (so that the thermodynamics are solved correctly)
    IF (C%choice_ice_dynamics == 'none') THEN
      C%choice_ice_dynamics = 'DIVA'
      CALL solve_DIVA( region%mesh, region%ice)
      C%choice_ice_dynamics = 'none'
    END IF

    ! == Model wind-up
    ! ================

    CALL run_model_windup( region)

    ! ===== Scalar ice data =====
    ! ===========================

    ! Calculate ice sheet metadata (volume, area, GMSL contribution),
    ! for writing to the first time point of the output file
    CALL calculate_PD_sealevel_contribution( region)
    CALL calculate_icesheet_volume_and_area( region)

    ! ===== Regional output =====
    ! ===========================

    ! Create output file for regional scalar data
    CALL create_regional_scalar_output_file( region)

    ! Create output file for regional 2-/3-D data
    CALL create_output_files( region)
    region%output_file_exists = .TRUE.

    ! Write scalar data at time t=0
    CALL write_regional_scalar_data( region, C%start_time_of_run)

    ! Write 2-/3-D data at time t=0
    CALL write_to_output_files( region)

    ! === Finalisation ===
    ! ====================

    IF (par%master) WRITE (0,*) ' Finished initialising model region ', region%name, '.'

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = HUGE( 1))

  END SUBROUTINE initialise_model

! ===== Auxiliary routines =====
! ==============================

  ! Allocate some basic variables
  SUBROUTINE allocate_region_timers_and_scalars( region)
    ! Allocate shared memory for this region's timers (used for the asynchronous coupling between the
    ! ice dynamics and the secondary model components), and for the scalars (integrated ice volume and
    ! area, SMB components, computation times, etc.)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),       INTENT(INOUT)  :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'allocate_region_timers_and_scalars'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Timers and time steps
    ! =====================

    CALL allocate_shared_dp_0D(   region%time,             region%wtime            )
    CALL allocate_shared_dp_0D(   region%dt,               region%wdt              )
    CALL allocate_shared_dp_0D(   region%dt_prev,          region%wdt_prev         )
    CALL allocate_shared_dp_0D(   region%dt_crit_SIA,      region%wdt_crit_SIA     )
    CALL allocate_shared_dp_0D(   region%dt_crit_SSA,      region%wdt_crit_SSA     )
    CALL allocate_shared_dp_0D(   region%dt_crit_ice,      region%wdt_crit_ice     )
    CALL allocate_shared_dp_0D(   region%dt_crit_ice_prev, region%wdt_crit_ice_prev)

    CALL allocate_shared_dp_0D(   region%t_last_mesh,      region%wt_last_mesh     )
    CALL allocate_shared_dp_0D(   region%t_next_mesh,      region%wt_next_mesh     )
    CALL allocate_shared_bool_0D( region%do_mesh,          region%wdo_mesh         )

    CALL allocate_shared_dp_0D(   region%t_last_SIA,       region%wt_last_SIA      )
    CALL allocate_shared_dp_0D(   region%t_next_SIA,       region%wt_next_SIA      )
    CALL allocate_shared_bool_0D( region%do_SIA,           region%wdo_SIA          )

    CALL allocate_shared_dp_0D(   region%t_last_SSA,       region%wt_last_SSA      )
    CALL allocate_shared_dp_0D(   region%t_next_SSA,       region%wt_next_SSA      )
    CALL allocate_shared_bool_0D( region%do_SSA,           region%wdo_SSA          )

    CALL allocate_shared_dp_0D(   region%t_last_DIVA,      region%wt_last_DIVA     )
    CALL allocate_shared_dp_0D(   region%t_next_DIVA,      region%wt_next_DIVA     )
    CALL allocate_shared_bool_0D( region%do_DIVA,          region%wdo_DIVA         )

    CALL allocate_shared_dp_0D(   region%t_last_thermo,    region%wt_last_thermo   )
    CALL allocate_shared_dp_0D(   region%t_next_thermo,    region%wt_next_thermo   )
    CALL allocate_shared_bool_0D( region%do_thermo,        region%wdo_thermo       )

    CALL allocate_shared_dp_0D(   region%t_last_climate,   region%wt_last_climate  )
    CALL allocate_shared_dp_0D(   region%t_next_climate,   region%wt_next_climate  )
    CALL allocate_shared_bool_0D( region%do_climate,       region%wdo_climate      )

    CALL allocate_shared_dp_0D(   region%t_last_ocean,     region%wt_last_ocean    )
    CALL allocate_shared_dp_0D(   region%t_next_ocean,     region%wt_next_ocean    )
    CALL allocate_shared_bool_0D( region%do_ocean,         region%wdo_ocean        )

    CALL allocate_shared_dp_0D(   region%t_last_SMB,       region%wt_last_SMB      )
    CALL allocate_shared_dp_0D(   region%t_next_SMB,       region%wt_next_SMB      )
    CALL allocate_shared_bool_0D( region%do_SMB,           region%wdo_SMB          )

    CALL allocate_shared_dp_0D(   region%t_last_BMB,       region%wt_last_BMB      )
    CALL allocate_shared_dp_0D(   region%t_next_BMB,       region%wt_next_BMB      )
    CALL allocate_shared_bool_0D( region%do_BMB,           region%wdo_BMB          )

    CALL allocate_shared_dp_0D(   region%t_last_ELRA,      region%wt_last_ELRA     )
    CALL allocate_shared_dp_0D(   region%t_next_ELRA,      region%wt_next_ELRA     )
    CALL allocate_shared_bool_0D( region%do_ELRA,          region%wdo_ELRA         )

    CALL allocate_shared_dp_0D(   region%t_last_basal,     region%wt_last_basal    )
    CALL allocate_shared_dp_0D(   region%t_next_basal,     region%wt_next_basal    )
    CALL allocate_shared_bool_0D( region%do_basal,         region%wdo_basal        )

    CALL allocate_shared_dp_0D(   region%t_last_SMB_inv,   region%wt_last_SMB_inv  )
    CALL allocate_shared_dp_0D(   region%t_next_SMB_inv,   region%wt_next_SMB_inv  )
    CALL allocate_shared_bool_0D( region%do_SMB_inv,       region%wdo_SMB_inv      )

    CALL allocate_shared_dp_0D(   region%t_last_output,    region%wt_last_output   )
    CALL allocate_shared_dp_0D(   region%t_next_output,    region%wt_next_output   )
    CALL allocate_shared_bool_0D( region%do_output,        region%wdo_output       )

    IF (par%master) THEN

      region%time             = C%start_time_of_run
      region%dt               = C%dt_min
      region%dt_prev          = C%dt_min
      region%dt_crit_SIA      = C%dt_min
      region%dt_crit_SSA      = C%dt_min
      region%dt_crit_ice      = C%dt_min
      region%dt_crit_ice_prev = C%dt_min

      region%t_last_mesh      = C%start_time_of_run
      region%t_next_mesh      = C%start_time_of_run + C%dt_mesh_min
      IF (C%do_force_mesh_update) THEN
        region%do_mesh          = .TRUE.
      END IF

      region%t_last_SIA       = C%start_time_of_run
      region%t_next_SIA       = C%start_time_of_run
      region%do_SIA           = .TRUE.

      region%t_last_SSA       = C%start_time_of_run
      region%t_next_SSA       = C%start_time_of_run
      region%do_SSA           = .TRUE.

      region%t_last_DIVA      = C%start_time_of_run
      region%t_next_DIVA      = C%start_time_of_run
      region%do_DIVA          = .TRUE.

      region%t_last_thermo    = C%start_time_of_run
      region%t_next_thermo    = C%start_time_of_run + C%dt_thermo
      region%do_thermo        = .FALSE.

      region%t_last_climate   = C%start_time_of_run
      region%t_next_climate   = C%start_time_of_run
      region%do_climate       = .TRUE.

      region%t_last_ocean     = C%start_time_of_run
      region%t_next_ocean     = C%start_time_of_run
      region%do_ocean         = .TRUE.

      region%t_last_SMB       = C%start_time_of_run
      region%t_next_SMB       = C%start_time_of_run
      region%do_SMB           = .TRUE.

      region%t_last_BMB       = C%start_time_of_run
      region%t_next_BMB       = C%start_time_of_run
      region%do_BMB           = .TRUE.

      region%t_last_ELRA      = C%start_time_of_run
      region%t_next_ELRA      = C%start_time_of_run
      IF (C%choice_GIA_model == 'ELRA') THEN
        region%do_ELRA        = .TRUE.
      ELSE
        region%do_ELRA        = .FALSE.
      END IF

      region%t_last_basal     = C%start_time_of_run
      region%t_next_basal     = C%start_time_of_run + C%dt_basal
      region%do_basal         = .FALSE.

      region%t_last_SMB_inv   = C%start_time_of_run
      region%t_next_SMB_inv   = C%start_time_of_run + C%dt_SMB_inv
      region%do_SMB_inv       = .FALSE.

      region%t_last_output    = C%start_time_of_run
      region%t_next_output    = C%start_time_of_run
      region%do_output        = .TRUE.
    END IF

    ! ===== Scalars =====
    ! ===================

    ! Ice-sheet volume and area
    CALL allocate_shared_dp_0D( region%ice_area                     , region%wice_area                     )
    CALL allocate_shared_dp_0D( region%ice_volume                   , region%wice_volume                   )
    CALL allocate_shared_dp_0D( region%ice_volume_PD                , region%wice_volume_PD                )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation   , region%wice_volume_above_flotation   )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation_PD, region%wice_volume_above_flotation_PD)

    ! Regionally integrated SMB components
    CALL allocate_shared_dp_0D( region%int_T2m                      , region%wint_T2m                      )
    CALL allocate_shared_dp_0D( region%int_snowfall                 , region%wint_snowfall                 )
    CALL allocate_shared_dp_0D( region%int_rainfall                 , region%wint_rainfall                 )
    CALL allocate_shared_dp_0D( region%int_melt                     , region%wint_melt                     )
    CALL allocate_shared_dp_0D( region%int_refreezing               , region%wint_refreezing               )
    CALL allocate_shared_dp_0D( region%int_runoff                   , region%wint_runoff                   )
    CALL allocate_shared_dp_0D( region%int_SMB                      , region%wint_SMB                      )
    CALL allocate_shared_dp_0D( region%int_BMB                      , region%wint_BMB                      )
    CALL allocate_shared_dp_0D( region%int_MB                       , region%wint_MB                       )

    ! Englacial isotope content
    CALL allocate_shared_dp_0D( region%GMSL_contribution            , region%wGMSL_contribution            )
    CALL allocate_shared_dp_0D( region%mean_isotope_content         , region%wmean_isotope_content         )
    CALL allocate_shared_dp_0D( region%mean_isotope_content_PD      , region%wmean_isotope_content_PD      )
    CALL allocate_shared_dp_0D( region%d18O_contribution            , region%wd18O_contribution            )
    CALL allocate_shared_dp_0D( region%d18O_contribution_PD         , region%wd18O_contribution_PD         )

    ! Computation times
    CALL allocate_shared_dp_0D( region%tcomp_total                  , region%wtcomp_total                  )
    CALL allocate_shared_dp_0D( region%tcomp_ice                    , region%wtcomp_ice                    )
    CALL allocate_shared_dp_0D( region%tcomp_thermo                 , region%wtcomp_thermo                 )
    CALL allocate_shared_dp_0D( region%tcomp_climate                , region%wtcomp_climate                )
    CALL allocate_shared_dp_0D( region%tcomp_GIA                    , region%wtcomp_GIA                    )
    CALL allocate_shared_dp_0D( region%tcomp_mesh                   , region%wtcomp_mesh                   )

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 71)

  END SUBROUTINE allocate_region_timers_and_scalars

  ! Initialise a regular grid
  SUBROUTINE initialise_model_square_grid( region, grid, dx)
    ! Initialise a regular square grid enveloping this model
    ! region based on the mesh data and settings.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region), INTENT(INOUT)  :: region
    TYPE(type_grid),         INTENT(INOUT)  :: grid
    REAL(dp),                INTENT(IN)     :: dx

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER           :: routine_name = 'initialise_model_square_grid'
    REAL(dp)                                :: xmid, ymid
    INTEGER                                 :: nsx, nsy, i, j, n
    REAL(dp), PARAMETER                     :: tol = 1E-9_dp

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate basic grid data
    CALL allocate_shared_int_0D( grid%nx,           grid%wnx          )
    CALL allocate_shared_int_0D( grid%ny,           grid%wny          )
    CALL allocate_shared_int_0D( grid%n,            grid%wn           )
    CALL allocate_shared_dp_0D(  grid%dx,           grid%wdx          )
    CALL allocate_shared_dp_0D(  grid%lambda_M,     grid%wlambda_M    )
    CALL allocate_shared_dp_0D(  grid%phi_M,        grid%wphi_M       )
    CALL allocate_shared_dp_0D(  grid%alpha_stereo, grid%walpha_stereo)
    CALL allocate_shared_dp_0D(  grid%xmin,         grid%wxmin        )
    CALL allocate_shared_dp_0D(  grid%xmax,         grid%wxmax        )
    CALL allocate_shared_dp_0D(  grid%ymin,         grid%wymin        )
    CALL allocate_shared_dp_0D(  grid%ymax,         grid%wymax        )
    CALL allocate_shared_dp_0D(  grid%tol_dist,     grid%wtol_dist    )


    ! === Basic grid data ===
    ! =======================

    ! Let master do this
    IF (par%master) THEN

      ! Resolution
      grid%dx = dx

      ! Projection parameters for this region
      IF     (region%name == 'NAM') THEN
        grid%lambda_M     = C%lambda_M_NAM
        grid%phi_M        = C%phi_M_NAM
        grid%alpha_stereo = C%alpha_stereo_NAM
      ELSEIF (region%name == 'EAS') THEN
        grid%lambda_M     = C%lambda_M_EAS
        grid%phi_M        = C%phi_M_EAS
        grid%alpha_stereo = C%alpha_stereo_EAS
      ELSEIF (region%name == 'GRL') THEN
        grid%lambda_M     = C%lambda_M_GRL
        grid%phi_M        = C%phi_M_GRL
        grid%alpha_stereo = C%alpha_stereo_GRL
      ELSEIF (region%name == 'ANT') THEN
        grid%lambda_M     = C%lambda_M_ANT
        grid%phi_M        = C%phi_M_ANT
        grid%alpha_stereo = C%alpha_stereo_ANT
      END IF

      ! Determine the center of the model domain
      xmid = (region%mesh%xmin + region%mesh%xmax) / 2._dp
      ymid = (region%mesh%ymin + region%mesh%ymax) / 2._dp

      ! Number of points at each side of domain center
      nsx = FLOOR((region%mesh%xmax - xmid) / grid%dx)
      nsy = FLOOR((region%mesh%ymax - ymid) / grid%dx)

      ! Determine total number of points per dimension as twice the
      ! number per side, plus the middle point, plus 1 grid cell to
      ! make sure the mesh lies completely inside the grid for any
      ! grid resolution
      grid%nx = (2*nsx + 1) + 1
      grid%ny = (2*nsy + 1) + 1

      ! Determine total number of grid points
      grid%n  = grid%nx * grid%ny

    END IF
    CALL sync

    ! Assign range to each processor
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

    ! === Grid generation ===
    ! =======================

    ! Allocate x and y coordinates for square grid
    CALL allocate_shared_dp_1D( grid%nx, grid%x, grid%wx)
    CALL allocate_shared_dp_1D( grid%ny, grid%y, grid%wy)

    ! Let master do this
    IF (par%master) THEN

      ! Compute corners of grid
      grid%xmin = grid%x(1      )
      grid%xmax = grid%x(grid%nx)
      grid%ymin = grid%y(1      )
      grid%ymax = grid%y(grid%ny)

      ! Compute x coordinates
      grid%xmin = xmid - nsx * grid%dx
      grid%xmax = xmid + nsx * grid%dx
      DO i = 1, grid%nx
        grid%x( i) = grid%xmin + (i-1)*grid%dx
      END DO

      ! Compute y coordinates
      grid%ymin = ymid - nsy * grid%dx
      grid%ymax = ymid + nsy * grid%dx
      DO j = 1, grid%ny
        grid%y( j) = grid%ymin + (j-1)*grid%dx
      END DO

    END IF
    CALL sync

    ! === Grid-to-vector tables ===
    ! =============================

    ! Allocate table data
    CALL allocate_shared_int_2D( grid%nx, grid%ny, grid%ij2n, grid%wij2n)
    CALL allocate_shared_int_2D( grid%n , 2      , grid%n2ij, grid%wn2ij)

    ! Tolerance; points lying within this distance are treated as identical
    grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

    ! Set up grid-to-vector translation tables
    IF (par%master) THEN
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
    END IF
    CALL sync

    ! === Mappings between mesh and grid ===
    ! ======================================

    ! Calculate mapping arrays between the mesh and the grid
    CALL calc_remapping_operator_mesh2grid( region%mesh, grid)
    CALL calc_remapping_operator_grid2mesh( grid, region%mesh)

    ! === Geographical coordinates ===
    ! ================================

    ! Calculate lat-lon coordinates
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, grid%lat, grid%wlat)
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, grid%lon, grid%wlon)

    ! Compute the lat/lon coordinate for each grid
    ! point using the mesh projection parameters
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL inverse_oblique_sg_projection( grid%x( i), grid%y( j), region%mesh%lambda_M, region%mesh%phi_M, region%mesh%alpha_stereo, grid%lon( i,j), grid%lat( i,j))
    END DO
    END DO
    CALL sync

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 18)

  END SUBROUTINE initialise_model_square_grid

  SUBROUTINE run_model_windup( region)
    ! After a restart, go back in time a set amount of years and run only
    ! the ice dynamics model (without ice thickness updates) until the
    ! start of the run. This is useful to avoid shocks in the initial
    ! evolution of the ice due to an empty, fresh-new velocity solver.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),       INTENT(INOUT)  :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'run_model_windup'
    INTEGER                                       :: it
    REAL(dp)                                      :: dt_ave, t_end

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%windup_total_years <= 0._dp) THEN
      ! Finalise routine path
      CALL finalise_routine( routine_name)
      ! Exit rutine
      RETURN
    END IF

    ! === Wind up the model ===
    ! =========================

    IF (par%master) THEN
      WRITE (0,*) '  Winding up model region ', region%name, ' for a nice restart...'
    END IF

    ! Save current time as end-time for wind-up
    t_end = region%time

    ! Bring the timer C%windup_total_years years back in time
    region%time = region%time - C%windup_total_years

    ! Let the model know we want to run velocities from this point on
    region%t_last_SIA       = region%time
    region%t_next_SIA       = region%time
    region%do_SIA           = .TRUE.

    region%t_last_SSA       = region%time
    region%t_next_SSA       = region%time
    region%do_SSA           = .TRUE.

    region%t_last_DIVA      = region%time
    region%t_next_DIVA      = region%time
    region%do_DIVA          = .TRUE.

    region%t_last_thermo    = region%time
    region%t_next_thermo    = region%time
    region%do_thermo        = .TRUE.

    ! Asume mid-value dt's from previous run
    region%dt               = C%dt_max * .5_dp
    region%dt_prev          = C%dt_max * .5_dp
    region%dt_crit_SIA      = C%dt_max * .5_dp
    region%dt_crit_SSA      = C%dt_max * .5_dp
    region%dt_crit_ice      = C%dt_max * .5_dp
    region%dt_crit_ice_prev = C%dt_max * .5_dp

    ! Initialise iteration counter
    it = 0
    ! Initialise averaged time step
    dt_ave = 0._dp

    ! Run the ice model until coming back to present
    DO WHILE (region%time < t_end)

      ! Update iteration counter
      it = it + 1

      ! Calculate ice velocities
      CALL run_ice_model( region, t_end)

      ! Calculate ice temperatures
      CALL run_thermo_model( region%mesh, region%ice, region%climate_matrix%applied, &
                             region%ocean_matrix%applied, region%SMB, region%time, &
                             do_solve_heat_equation = region%do_thermo)

      ! Display progress
      if (par%master .AND. C%do_time_display) then
        call time_display(region, t_end, dt_ave, it)
      end if

      ! Advance region time and repeat
      IF (par%master) THEN
        region%time = region%time + region%dt
        dt_ave = dt_ave + region%dt
      END IF
      CALL sync

    END DO

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model_windup

END MODULE UFEMISM_main_model