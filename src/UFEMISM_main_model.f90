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
                                                 create_debug_file, write_PETSc_matrix_to_NetCDF
  USE data_types_module,                   ONLY: type_model_region, type_mesh, type_grid, type_remapping_mesh_mesh, &
                                                 type_climate_matrix_global, type_ocean_matrix_global
  USE reference_fields_module,             ONLY: initialise_reference_geometries, map_reference_geometries_to_mesh, remap_restart_init_topo
  USE mesh_memory_module,                  ONLY: deallocate_mesh_all
  USE mesh_help_functions_module,          ONLY: inverse_oblique_sg_projection
  USE mesh_creation_module,                ONLY: create_mesh_from_cart_data
  USE mesh_mapping_module,                 ONLY: calc_remapping_operators_mesh_mesh, deallocate_remapping_operators_mesh_mesh, &
                                                 calc_remapping_operator_mesh2grid, deallocate_remapping_operators_mesh2grid, &
                                                 calc_remapping_operator_grid2mesh, deallocate_remapping_operators_grid2mesh
  USE mesh_update_module,                  ONLY: determine_mesh_fitness, create_new_mesh
  USE general_ice_model_data_module,       ONLY: initialise_mask_noice, initialise_basins
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
  USE restart_module,                      ONLY: read_mesh_from_restart_file, read_init_data_from_restart_file
  USE general_sea_level_module,            ONLY: calculate_PD_sealevel_contribution
  USE ice_velocity_module,                 ONLY: solve_DIVA
  USE text_output_module,                  ONLY: create_regional_text_output, write_regional_text_output

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
    TYPE(type_model_region),           INTENT(INOUT)     :: region
    TYPE(type_climate_matrix_global),  INTENT(INOUT)     :: climate_matrix_global
    REAL(dp),                          INTENT(IN)        :: t_end

    ! Local variables:
    CHARACTER(LEN=256)                                   :: routine_name
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

    DO WHILE (region%time < t_end)

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
      IF (par%master) t2 = MPI_WTIME()
      meshfitness = 1._dp
      IF (region%time > region%t_last_mesh + C%dt_mesh_min) THEN
        CALL determine_mesh_fitness(region%mesh, region%ice, meshfitness)
      END IF
      IF (par%master) region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2

      ! If required, update the mesh
      IF (meshfitness < C%mesh_fitness_threshold) THEN
      ! IF (.FALSE.) THEN
      ! IF (.TRUE.) THEN
        region%t_last_mesh = region%time
        IF (par%master) t2 = MPI_WTIME()
        CALL run_model_update_mesh( region, climate_matrix_global)
        IF (par%master) region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2
      END IF

      ! == Ice dynamics
      ! ===============

      ! Calculate ice velocities and the resulting change in ice geometry
      ! NOTE: geometry is not updated yet; this happens at the end of the time loop
      t1 = MPI_WTIME()
      CALL run_ice_model( region, t_end)
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_ice = region%tcomp_ice + t2 - t1

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

      IF (C%do_basal_sliding_inversion) THEN
        IF (region%do_basal) THEN
          CALL basal_sliding_inversion( region%mesh, region%grid_smooth, region%ice, region%refgeo_PD)
        END IF
      END IF

      ! == SBM IMAU-ITM inversion
      ! =========================

      IF (C%do_SMB_IMAUITM_inversion) THEN
        IF (region%do_SMB_inv) THEN
          CALL SMB_IMAUITM_inversion( region%mesh, region%ice, region%climate_matrix%applied, region%SMB, region%refgeo_PD)
        END IF
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
        ! Write to regional NetCDF output files
        CALL write_to_output_files( region)
      END IF

      ! Write scalar text output
      IF (region%do_output .OR. region%do_SMB .OR. region%do_BMB) THEN
        ! Determine total ice sheet area, volume, and volume-above-flotation,
        CALL calculate_icesheet_volume_and_area( region)
        ! Write to regional text output file
        CALL write_regional_text_output( region)
      END IF

      ! == Update ice geometry
      ! ======================

      CALL update_ice_thickness( region%mesh, region%ice, region%mask_noice, region%refgeo_PD, region%refgeo_GIAeq)
      CALL sync

      ! == Advance region time
      ! ======================

      IF (par%master) region%time = region%time + region%dt
      CALL sync

    END DO

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
      CALL write_regional_text_output( region)
    END IF

    ! Determine total ice sheet area, volume, volume-above-flotation and GMSL contribution,
    ! used for writing to global text output and in the inverse routine
    CALL calculate_icesheet_volume_and_area( region)

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
    TYPE(type_climate_matrix_global),    INTENT(IN)    :: climate_matrix_global

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

    IF (C%is_restart) THEN
      ! Map restart initial geometry from old to new mesh (if needed)
      CALL remap_restart_init_topo( region, region%mesh, region%mesh_new, map)
    END IF

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
    CALL remap_climate_model(  region%mesh, region%mesh_new, map, region%climate_matrix, climate_matrix_global, region%refgeo_PD, region%grid_smooth, region%mask_noice, region%name)
    CALL remap_ocean_model(    region%mesh, region%mesh_new, map, region%ocean_matrix)
    CALL remap_SMB_model(      region%mesh, region%mesh_new, map, region%SMB)
    CALL remap_BMB_model(      region%mesh, region%mesh_new, map, region%BMB)
    CALL remap_isotopes_model( region%mesh, region%mesh_new, map, region)

    ! Deallocate shared memory for the mapping arrays
    CALL deallocate_remapping_operators_mesh_mesh( map)

    ! Deallocate the old mesh, bind the region%mesh pointers to the new mesh.
    CALL deallocate_mesh_all( region%mesh)
    region%mesh = region%mesh_new

    ! When the next output is written, new output files must be created.
    region%output_file_exists = .FALSE.

    ! Reallocate the debug fields for the new mesh, create a new debug file
    CALL reallocate_debug_fields( region)
    CALL create_debug_file(       region)

    ! Recalculate the reference precipitation isotope content (must be done after region%mesh has been cycled)
    CALL calculate_reference_isotopes( region)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model_update_mesh

  ! Initialise the entire model region
  SUBROUTINE initialise_model( region, name, climate_matrix_global, ocean_matrix_global)
    ! Initialise the entire model region - read initial and PD data, create the mesh,
    ! and initialise the ice dynamics, climate, ocean, and SMB sub models

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),          INTENT(INOUT)     :: region
    CHARACTER(LEN=3),                 INTENT(IN)        :: name
    TYPE(type_climate_matrix_global), INTENT(INOUT)     :: climate_matrix_global
    TYPE(type_ocean_matrix_global),   INTENT(INOUT)     :: ocean_matrix_global

    ! Local variables:
    CHARACTER(LEN=256)                                  :: routine_name

    ! Add routine to path
    routine_name = 'initialise_model('  //  name  //  ')'
    CALL init_routine( routine_name)

    ! ===== Basic initialisation =====
    ! ================================

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

    ! ===== Allocate memory for timers and scalars =====
    ! ==================================================

    CALL allocate_region_timers_and_scalars( region)

    ! ===== PD, GIAeq, and init reference data fields =====
    ! =====================================================

    CALL initialise_reference_geometries( region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq, region%name)

    ! ===== The mesh =====
    ! ====================

    IF (C%is_restart) THEN
      ! Read mesh and data (on the mesh) from a restart file
      CALL read_mesh_from_restart_file( region)
      CALL read_init_data_from_restart_file( region)
    ELSE
      ! Create a new mesh from the reference initial geometry
      CALL create_mesh_from_cart_data( region)
    END IF

    ! ===== Map reference geometries to the mesh =====
    ! ================================================

    IF (par%master) WRITE(0,*) '  Mapping reference geometries onto the initial model mesh...'

    ! Allocate memory for reference data on the mesh
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_init%Hi , region%refgeo_init%wHi )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_init%Hb , region%refgeo_init%wHb )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_init%Hs , region%refgeo_init%wHs )

    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_PD%Hi   , region%refgeo_PD%wHi   )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_PD%Hb   , region%refgeo_PD%wHb   )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_PD%Hs   , region%refgeo_PD%wHs   )

    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_GIAeq%Hi, region%refgeo_GIAeq%wHi)
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_GIAeq%Hb, region%refgeo_GIAeq%wHb)
    CALL allocate_shared_dp_1D( region%mesh%nV, region%refgeo_GIAeq%Hs, region%refgeo_GIAeq%wHs)

    ! Map data from the square grids to the mesh
    CALL map_reference_geometries_to_mesh( region, region%mesh)

    IF (par%master) WRITE(0,*) '  Finished mapping reference geometries.'

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

    ! ===== The "no ice" mask =====
    ! =============================

    CALL allocate_shared_int_1D( region%mesh%nV, region%mask_noice, region%wmask_noice)
    CALL initialise_mask_noice( region, region%mesh)

    ! ===== The ice dynamics model =====
    ! ==================================

    CALL initialise_ice_model( region%mesh, region%ice, region%refgeo_init, region%refgeo_PD, region%restart)

    ! ===== Define ice basins =====
    ! =============================

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

    CALL initialise_BMB_model( region%mesh, region%ice, region%BMB, region%name)

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

    ! ===== Initialise the ice temperature profile =====
    ! ==================================================

    ! Run the climate and SMB models once, to get the correct surface temperature+SMB fields for the ice temperature initialisation
    CALL run_climate_model( region, climate_matrix_global, C%start_time_of_run)
    CALL run_SMB_model( region%mesh, region%ice, region%climate_matrix, C%start_time_of_run, region%SMB, region%mask_noice)

    ! Initialise the ice temperature profile
    CALL initialise_ice_temperature( region%mesh, region%ice, region%climate_matrix%applied, region%ocean_matrix%applied, region%SMB, region%name, region%restart)

    ! Initialise the rheology
    CALL calc_ice_rheology( region%mesh, region%ice, C%start_time_of_run)

    ! ===== Scalar output (regionally integrated ice volume, SMB components, etc.) =====
    ! ==================================================================================

    ! Calculate ice sheet metadata (volume, area, GMSL contribution) for writing to the first line of the output file
    CALL calculate_PD_sealevel_contribution( region)
    CALL calculate_icesheet_volume_and_area( region)

    ! ===== Regional ASCII output (Scalar data, POIs) =====
    ! =====================================================

    CALL create_regional_text_output( region)

    ! ===== Exception: Initial velocities for choice_ice_dynamics == "none" =====
    ! ===========================================================================

    ! If we're running with choice_ice_dynamics == "none", calculate a velocity field
    ! once during initialisation (so that the thermodynamics are solved correctly)
    IF (C%choice_ice_dynamics == 'none') THEN
      C%choice_ice_dynamics = 'DIVA'
      CALL solve_DIVA( region%mesh, region%ice)
      C%choice_ice_dynamics = 'none'
    END IF

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
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_region_timers_and_scalars'

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
      IF (C%choice_GIA_model == 'ELRA') THEN
        region%do_ELRA      = .TRUE.
      ELSE
        region%do_ELRA      = .FALSE.
      END IF

      region%t_last_basal   = C%start_time_of_run
      region%t_next_basal   = C%start_time_of_run + C%dt_basal
      region%do_basal       = .FALSE.

      region%t_last_SMB_inv = C%start_time_of_run
      region%t_next_SMB_inv = C%start_time_of_run + C%dt_SMB_inv
      region%do_SMB_inv     = .FALSE.

      region%t_last_output  = C%start_time_of_run
      region%t_next_output  = C%start_time_of_run
      region%do_output      = .TRUE.
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
    ! Initialise a regular square grid enveloping this model region

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    TYPE(type_grid),            INTENT(INOUT)     :: grid
    REAL(dp),                   INTENT(IN)        :: dx

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_model_square_grid'
    REAL(dp)                                      :: xmid, ymid
    INTEGER                                       :: nsx, nsy, i, j, n
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    nsx = 0
    nsy = 0

    CALL allocate_shared_int_0D( grid%nx,           grid%wnx          )
    CALL allocate_shared_int_0D( grid%ny,           grid%wny          )
    CALL allocate_shared_dp_0D(  grid%dx,           grid%wdx          )
    CALL allocate_shared_dp_0D(  grid%lambda_M,     grid%wlambda_M    )
    CALL allocate_shared_dp_0D(  grid%phi_M,        grid%wphi_M       )
    CALL allocate_shared_dp_0D(  grid%alpha_stereo, grid%walpha_stereo)
    CALL allocate_shared_dp_0D(  grid%xmin,         grid%wxmin        )
    CALL allocate_shared_dp_0D(  grid%xmax,         grid%wxmax        )
    CALL allocate_shared_dp_0D(  grid%ymin,         grid%wymin        )
    CALL allocate_shared_dp_0D(  grid%ymax,         grid%wymax        )

    IF (par%master) THEN

      ! Resolution
      grid%dx = dx

      ! Projection parameters for this region (determined from the config)
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

      nsx = CEILING((region%mesh%xmax - xmid - grid%dx/2._dp) / grid%dx)
      nsy = CEILING((region%mesh%ymax - ymid - grid%dx/2._dp) / grid%dx)

      grid%nx = 2*nsx + 1
      grid%ny = 2*nsy + 1

    END IF ! IF (par%master) THEN
    CALL sync

    CALL allocate_shared_dp_1D( grid%nx, grid%x, grid%wx)
    CALL allocate_shared_dp_1D( grid%ny, grid%y, grid%wy)

    ! Fill in x and y
    IF (par%master) THEN
    
      ! x
      grid%xmin = xmid - nsx * grid%dx
      grid%xmax = xmid + nsx * grid%dx
      DO i = 1, grid%nx
        grid%x( i) = grid%xmin + (i-1)*grid%dx
      END DO
    
      ! y
      grid%ymin = ymid - nsy * grid%dx
      grid%ymax = ymid + nsy * grid%dx
      DO j = 1, grid%ny
        grid%y( j) = grid%ymin + (j-1)*grid%dx
      END DO
      
    END IF ! IF (par%master) THEN
    CALL sync

    ! Tolerance; points lying within this distance of each other are treated as identical
    CALL allocate_shared_dp_0D( grid%tol_dist, grid%wtol_dist)
    IF (par%master) grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp
    CALL sync

    ! Set up grid-to-vector translation tables
    CALL allocate_shared_int_0D(                   grid%n           , grid%wn           )
    IF (par%master) grid%n  = grid%nx * grid%ny
    CALL sync
    CALL allocate_shared_int_2D( grid%nx, grid%ny, grid%ij2n        , grid%wij2n        )
    CALL allocate_shared_int_2D( grid%n , 2      , grid%n2ij        , grid%wn2ij        )
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

    ! Assign range to each processor
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

    ! Calculate lat-lon coordinates
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, grid%lat, grid%wlat)
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, grid%lon, grid%wlon)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL inverse_oblique_sg_projection( grid%x( i), grid%y( j), region%mesh%lambda_M, region%mesh%phi_M, region%mesh%alpha_stereo, grid%lon( i,j), grid%lat( i,j))
    END DO
    END DO
    CALL sync

    ! Calculate mapping arrays between the mesh and the grid
    CALL calc_remapping_operator_mesh2grid( region%mesh, grid)
    CALL calc_remapping_operator_grid2mesh( grid, region%mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 18)

  END SUBROUTINE initialise_model_square_grid

  ! Calculate regional ice volume and area
  SUBROUTINE calculate_icesheet_volume_and_area( region)
    ! Calculate this region's ice sheet's volume and area

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calculate_icesheet_volume_and_area'
    INTEGER                                       :: vi
    REAL(dp)                                      :: ice_area, ice_volume, thickness_above_flotation, ice_volume_above_flotation

    ! Add routine to path
    CALL init_routine( routine_name)

    ice_area                   = 0._dp
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp

    ! Calculate ice area and volume for processor domain
    DO vi = region%mesh%vi1, region%mesh%vi2

      IF (region%ice%mask_ice_a( vi) == 1) THEN
        ice_volume = ice_volume + (region%ice%Hi_a( vi) * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area))
        ice_area   = ice_area   + region%mesh%A( vi) * 1.0E-06_dp ! [km^3]

        ! Thickness above flotation
        thickness_above_flotation = MAX(0._dp, region%ice%Hi_a( vi) - MAX(0._dp, (region%ice%SL_a( vi) - region%ice%Hb_a( vi)) * (seawater_density / ice_density)))

        ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area)
      END IF

    END DO
    CALL sync

    CALL MPI_REDUCE( ice_area,                   region%ice_area,                   1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume,                 region%ice_volume,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! Calculate GMSL contribution
    IF (par%master) region%GMSL_contribution = -1._dp * (region%ice_volume_above_flotation - region%ice_volume_above_flotation_PD)
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_icesheet_volume_and_area

END MODULE UFEMISM_main_model