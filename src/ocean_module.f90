MODULE ocean_module

  ! Contains all the routines for calculating the ocean forcing.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             reallocate_shared_dp_2D, deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_mesh, type_grid, type_model_region, type_ice_model, &
                                             type_ocean_matrix_global, type_ocean_snapshot_global, &
                                             type_ocean_matrix_regional, type_ocean_snapshot_regional, &
                                             type_highres_ocean_data, type_climate_matrix_regional, &
                                             type_remapping_mesh_mesh
  USE netcdf_module,                   ONLY: inquire_PD_obs_global_ocean_file, read_PD_obs_global_ocean_file, &
                                             inquire_GCM_global_ocean_file, read_GCM_global_ocean_file, &
                                             inquire_hires_geometry_file, read_hires_geometry_file, &
                                             inquire_extrapolated_ocean_file, read_extrapolated_ocean_file, &
                                             create_extrapolated_ocean_file
  USE forcing_module,                  ONLY: forcing, update_CO2_at_model_time
  USE utilities_module,                ONLY: remap_cons_2nd_order_1D, inverse_oblique_sg_projection, &
                                             map_glob_to_grid_3D, surface_elevation, &
                                             extrapolate_Gaussian_floodfill
  USE mesh_mapping_module,             ONLY: calc_remapping_operator_grid2mesh, map_grid2mesh_3D, &
                                             calc_remapping_operator_mesh2grid, map_mesh2grid_2D, &
                                             deallocate_remapping_operators_mesh2grid, &
                                             deallocate_remapping_operators_grid2mesh, &
                                             smooth_Gaussian_2D, remap_field_dp_3D

  IMPLICIT NONE

CONTAINS

! == The main routines that should be called from the main ice model/program
! ==========================================================================

  SUBROUTINE run_ocean_model( mesh, grid, ice, ocean_matrix, climate_matrix, region_name, time)
    ! Run the regional ocean model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TyPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_matrix_regional),    INTENT(INOUT) :: ocean_matrix
    TYPE(type_climate_matrix_regional),  INTENT(IN)    :: climate_matrix
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ocean_model == 'none') THEN
      ! No ocean data is used at all
    ELSEIF (C%choice_ocean_model == 'idealised') THEN
      ! Assign some idealised temperature/salinity ocean profiles

      CALL run_ocean_model_idealised( mesh, ice, ocean_matrix%applied, region_name, time)

      IF (par%master) WRITE(0,*) 'This subroutine (run_ocean_model) is not ready for this choice of ocean model...'
      IF (par%master) WRITE(0,*) 'Feel free to fill it before running the model :)'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ELSEIF (C%choice_ocean_model == 'uniform_warm_cold') THEN
      ! Uniform warm/cold ocean model (the old ANICE way)

      CALL run_ocean_model_uniform_warm_cold( mesh, ocean_matrix, time)

      IF (par%master) WRITE(0,*) 'This subroutine (run_ocean_model) is not ready for this choice of ocean model...'
      IF (par%master) WRITE(0,*) 'Feel free to fill it before running the model :)'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ELSEIF (C%choice_ocean_model == 'PD_obs') THEN
      ! Keep the ocean fixed to present-day observed conditions

      CALL run_ocean_model_PD_obs( mesh, ocean_matrix)

      IF (par%master) WRITE(0,*) 'This subroutine (run_ocean_model) is not ready for this choice of ocean model...'
      IF (par%master) WRITE(0,*) 'Feel free to fill it before running the model :)'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ELSEIF (C%choice_ocean_model == 'matrix_warm_cold') THEN
      ! Run the warm/cold ocean matrix

      CALL run_ocean_model_matrix_warm_cold( mesh, grid, ocean_matrix, climate_matrix, region_name, time)

    ELSE
      IF (par%master) WRITE(0,*) 'run_ocean_model - ERROR: unknown choice_ocean_model "', TRIM(C%choice_ocean_model), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model

  SUBROUTINE initialise_ocean_model_regional( region, ocean_matrix_global)
    ! Initialise the regional ocean model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_ocean_matrix_global),      INTENT(IN)    :: ocean_matrix_global

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ocean_model_regional'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '  Initialising regional ocean model "', TRIM(C%choice_ocean_model), '"...'

    IF     (C%choice_ocean_model == 'none') THEN
      ! No ocean data is used at all
    ELSEIF (C%choice_ocean_model == 'idealised') THEN
      ! Only initialise the "applied" regional snapshot

      CALL allocate_ocean_snapshot_regional( region%mesh, region%ocean_matrix%applied, name = 'applied')

    ELSEIF (C%choice_ocean_model == 'uniform_warm_cold') THEN
      ! Only uses the T_ocean_mean variables

      CALL allocate_shared_dp_0D( region%ocean_matrix%GCM_cold%T_ocean_mean, region%ocean_matrix%GCM_cold%wT_ocean_mean)
      CALL allocate_shared_dp_0D( region%ocean_matrix%GCM_PI%T_ocean_mean,   region%ocean_matrix%GCM_PI%wT_ocean_mean  )
      CALL allocate_shared_dp_0D( region%ocean_matrix%GCM_warm%T_ocean_mean, region%ocean_matrix%GCM_warm%wT_ocean_mean)

    ELSEIF (C%choice_ocean_model == 'PD_obs') THEN
      ! Allocate both the "PD_obs" and "applied" snapshots, and initialise the present-day observations

      CALL initialise_ocean_model_PD_obs_regional( region, ocean_matrix_global)

    ELSEIF (C%choice_ocean_model == 'matrix_warm_cold') THEN
    !   ! Allocate all the snapshots used in the warm/cold ocean matrix

      CALL initialise_ocean_matrix_regional( region, ocean_matrix_global)

    ELSE
      IF (par%master) WRITE(0,*) 'initialise_ocean_model_regional - ERROR: unknown choice_ocean_model "', TRIM(C%choice_ocean_model), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=37)

  END SUBROUTINE initialise_ocean_model_regional

  SUBROUTINE initialise_ocean_model_global( ocean_matrix)
    ! Initialise the global ocean model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ocean_matrix_global),      INTENT(INOUT) :: ocean_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ocean_model_global'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) ' Initialising global ocean model "', TRIM(C%choice_ocean_model), '"...'

    IF     (C%choice_ocean_model == 'none') THEN
      ! No global ocean data is used at all

    ELSEIF (C%choice_ocean_model == 'idealised') THEN
      ! No global ocean data is used at all

    ELSEIF (C%choice_ocean_model == 'uniform_warm_cold') THEN
      ! No global ocean data is used at all

    ELSEIF (C%choice_ocean_model == 'PD_obs') THEN
      ! Only initialise the "PD_obs" global snapshot

      CALL initialise_ocean_PD_obs_global( ocean_matrix%PD_obs, name = 'WOA18')

    ELSEIF (C%choice_ocean_model == 'matrix_warm_cold') THEN
      ! Allocate all the global snapshots used in the warm/cold ocean matrix

      CALL initialise_ocean_matrix_global( ocean_matrix)

    ELSE
      IF (par%master) WRITE(0,*) 'initialise_ocean_model_global - ERROR: unknown choice_ocean_model "', TRIM(C%choice_ocean_model), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=39)

  END SUBROUTINE initialise_ocean_model_global

  ! == Idealised ocean configurations
  ! =================================

  SUBROUTINE run_ocean_model_idealised( mesh, ice, ocean, region_name, time)
    ! Run the regional ocean model
    !
    ! Assign some idealised temperature/salinity ocean profiles

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TyPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_ocean_model == 'idealised') THEN
      IF (par%master) WRITE(0,*) 'run_ocean_model_idealised - ERROR: choice_ocean_model should be set to "idealised"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Choose an idealised ocean profile
    IF     (C%choice_idealised_ocean == 'MISMIP+_COLD') THEN
      CALL run_ocean_model_idealised_MISMIPplus_COLD( mesh, ocean)
    ELSEIF (C%choice_idealised_ocean == 'MISMIP+_WARM') THEN
      CALL run_ocean_model_idealised_MISMIPplus_WARM( mesh, ocean)
    ELSEIF (C%choice_idealised_ocean == 'MISOMIP1') THEN
      CALL run_ocean_model_idealised_MISOMIP1( mesh, ocean, time)
    ELSEIF (C%choice_idealised_ocean == 'Reese2018_ANT') THEN
      CALL run_ocean_model_idealised_Reese2018_ANT( mesh, ice, ocean, region_name)
    ELSE
      IF (par%master) WRITE(0,*) 'run_ocean_model_idealised - ERROR: unknown choice_idealised_ocean "', TRIM(C%choice_idealised_ocean), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_idealised
  SUBROUTINE run_ocean_model_idealised_MISMIPplus_COLD( mesh, ocean)
    ! Run the regional ocean model
    !
    ! Set the ocean temperature and salinity to the ISOMIP+ "COLD" profile (Asay-Davis et al., 2016, Table 5)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ocean_snapshot_regional),  INTENT(INOUT) :: ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_idealised_MISMIPplus_COLD'
    INTEGER                                            :: vi,k
    REAL(dp)                                           :: w
    REAL(dp), PARAMETER                                :: Tzero     = -1.9_dp    ! Sea surface temperature [degC] (originally T0, but that name is already taken...)
    REAL(dp), PARAMETER                                :: Tbot      = -1.9_dp    ! Sea floor   temperature [degC]
    REAL(dp), PARAMETER                                :: Szero     = 33.8_dp    ! Sea surface salinity    [PSU]
    REAL(dp), PARAMETER                                :: Sbot      = 34.55_dp   ! Sea floor   salinity    [PSU]
    REAL(dp), PARAMETER                                :: depth_max = 720._dp    ! Maximum depth for the profile (constant values below that)

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Fill in the temperature and salinity profiles
    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, C%nz_ocean

      ! Interpolation weight
      w = MIN(1._dp, MAX(0._dp, C%z_ocean( k) / depth_max ))

      ! Temperature
      ocean%T_ocean(          vi,k) = w * Tbot + (1._dp - w) * Tzero
      ocean%T_ocean_ext(      vi,k) = w * Tbot + (1._dp - w) * Tzero
      ocean%T_ocean_corr_ext( vi,k) = w * Tbot + (1._dp - w) * Tzero

      ! Salinity
      ocean%S_ocean(          vi,k) = w * Sbot + (1._dp - w) * Szero
      ocean%S_ocean_ext(      vi,k) = w * Sbot + (1._dp - w) * Szero
      ocean%S_ocean_corr_ext( vi,k) = w * Sbot + (1._dp - w) * Szero

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_idealised_MISMIPplus_COLD
  SUBROUTINE run_ocean_model_idealised_MISMIPplus_WARM( mesh, ocean)
    ! Run the regional ocean model
    !
    ! Set the ocean temperature and salinity to the ISOMIP+ "WARM" profile (Asay-Davis et al., 2016, Table 6)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ocean_snapshot_regional),  INTENT(INOUT) :: ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_idealised_MISMIPplus_WARM'
    INTEGER                                            :: vi,k
    REAL(dp)                                           :: w
    REAL(dp), PARAMETER                                :: Tzero     = -1.9_dp    ! Sea surface temperature [degC] (originally T0, but that name is already taken...)
    REAL(dp), PARAMETER                                :: Tbot      =  1.0_dp    ! Sea floor   temperature [degC]
    REAL(dp), PARAMETER                                :: Szero     = 33.8_dp    ! Sea surface salinity    [PSU]
    REAL(dp), PARAMETER                                :: Sbot      = 34.7_dp    ! Sea floor   salinity    [PSU]
    REAL(dp), PARAMETER                                :: depth_max = 720._dp    ! Maximum depth for the profile (constant values below that)

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Fill in the temperature and salinity profiles
    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, C%nz_ocean

      ! Interpolation weight
      w = MIN(1._dp, MAX(0._dp, C%z_ocean( k) / depth_max ))

      ! Temperature
      ocean%T_ocean(          vi,k) = w * Tbot + (1._dp - w) * Tzero
      ocean%T_ocean_ext(      vi,k) = w * Tbot + (1._dp - w) * Tzero
      ocean%T_ocean_corr_ext( vi,k) = w * Tbot + (1._dp - w) * Tzero

      ! Salinity
      ocean%S_ocean(          vi,k) = w * Sbot + (1._dp - w) * Szero
      ocean%S_ocean_ext(      vi,k) = w * Sbot + (1._dp - w) * Szero
      ocean%S_ocean_corr_ext( vi,k) = w * Sbot + (1._dp - w) * Szero

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_idealised_MISMIPplus_WARM
  SUBROUTINE run_ocean_model_idealised_MISOMIP1( mesh, ocean, time)
    ! Run the regional ocean model
    !
    ! Set the ocean temperature and salinity according to the ISOMIP+ protocol

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ocean_snapshot_regional),  INTENT(INOUT) :: ocean
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_idealised_MISOMIP1'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%MISOMIP1_scenario == 'IceOcean0') THEN
      ! Cold ocean always

      CALL run_ocean_model_idealised_MISMIPplus_COLD( mesh, ocean)

    ELSEIF (C%MISOMIP1_scenario == 'IceOcean1ra' .OR. &
            C%MISOMIP1_scenario == 'IceOcean2ra') THEN
      ! Cold ocean during spin-up; warm ocean for 100 years, then cold ocean again

      IF (time < 0._dp) THEN
        CALL run_ocean_model_idealised_MISMIPplus_COLD( mesh, ocean)
      ELSEIF (time < 100._dp) THEN
        CALL run_ocean_model_idealised_MISMIPplus_WARM( mesh, ocean)
      ELSE
        CALL run_ocean_model_idealised_MISMIPplus_COLD( mesh, ocean)
      END IF

    ELSEIF (C%MISOMIP1_scenario == 'IceOcean1rr' .OR. &
            C%MISOMIP1_scenario == 'IceOcean2rr') THEN
      ! Cold ocean during spin-up; warm ocean after t = 0

      IF (time < 0._dp) THEN
        CALL run_ocean_model_idealised_MISMIPplus_COLD( mesh, ocean)
      ELSE
        CALL run_ocean_model_idealised_MISMIPplus_WARM( mesh, ocean)
      END IF

    ELSE
      IF (par%master) WRITE(0,*) 'run_ocean_model_idealised_MISOMIP1 - ERROR: unknown MISOMIP1_scenario "', TRIM(C%MISOMIP1_scenario), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_idealised_MISOMIP1
  SUBROUTINE run_ocean_model_idealised_Reese2018_ANT( mesh, ice, ocean, region_name)
    ! Run the regional ocean model
    !
    ! Set the ocean temperature and salinity to basin-dependent values
    ! provided by Reese et al. (2018) so that PICO gives realistic present-day melt rates

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_idealised_Reese2018_ANT'
    INTEGER                                            :: vi
    REAL(dp)                                           :: T,S

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. region_name == 'ANT') THEN
      IF (par%master) WRITE(0,*) 'run_ocean_model_idealised_Reese2018_ANT - ERROR: only applicable to Antarctica!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    IF (.NOT. (C%choice_basin_scheme_ANT == 'file' .AND. C%do_merge_basins_ANT)) THEN
      IF (par%master) THEN
        WRITE(0,*) ''
        WRITE(0,*) ' ===== '
        WRITE(0,*) 'run_ocean_model_idealised_Reese2018_ANT - WARNING: This really only works when using the external'
        WRITE(0,*) '  Antarctic ice basins file "ant_full_drainagesystem_polygons.txt". This can be downloaded from:'
        WRITE(0,*) '  https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems'
        WRITE(0,*) '  ...and you will also need to set do_merge_basins_ANT_config = .TRUE.'
        WRITE(0,*) ' ===== '
        WRITE(0,*) ''
      END IF
    END IF

    ! Fill in the temperature and salinity values
    T = 0._dp
    S = 0._dp
    DO vi = mesh%vi1, mesh%vi2

      IF     (ice%basin_ID( vi) == 1) THEN
        T = -1.76_dp
        S = 34.82_dp
      ELSEIF (ice%basin_ID( vi) == 2) THEN
        T = -1.66_dp
        S = 34.70_dp
      ELSEIF (ice%basin_ID( vi) == 3) THEN
        T = -1.65_dp
        S = 34.48_dp
      ELSEIF (ice%basin_ID( vi) == 4) THEN
        T = -1.58_dp
        S = 34.49_dp
      ELSEIF (ice%basin_ID( vi) == 5) THEN
        T = -1.51_dp
        S = 34.50_dp
      ELSEIF (ice%basin_ID( vi) == 6) THEN
        T = -1.73_dp
        S = 34.70_dp
      ELSEIF (ice%basin_ID( vi) == 7) THEN
        T = -1.68_dp
        S = 34.65_dp
      ELSEIF (ice%basin_ID( vi) == 8) THEN
        T = -0.73_dp
        S = 34.73_dp
      ELSEIF (ice%basin_ID( vi) == 9) THEN
        T = -1.61_dp
        S = 34.75_dp
      ELSEIF (ice%basin_ID( vi) == 10) THEN
        T = -1.30_dp
        S = 34.84_dp
      ELSEIF (ice%basin_ID( vi) == 11) THEN
        T = -1.58_dp
        S = 34.79_dp
      ELSEIF (ice%basin_ID( vi) == 12) THEN
        T = -0.36_dp
        S = 34.58_dp
      ELSEIF (ice%basin_ID( vi) == 13) THEN
        T =  0.80_dp
        S = 34.79_dp
      ELSEIF (ice%basin_ID( vi) == 14) THEN
        T =  1.10_dp
        S = 34.85_dp
      ELSEIF (ice%basin_ID( vi) == 15) THEN
        T =  0.23_dp
        S = 34.7_dp
      ELSEIF (ice%basin_ID( vi) == 16) THEN
        T = -1.23_dp
        S = 34.67_dp
      ELSEIF (ice%basin_ID( vi) == 17) THEN
        T = -1.80_dp
        S = 34.84_dp
      END IF

      ! Temperature
      ocean%T_ocean(          vi,:) = T
      ocean%T_ocean_ext(      vi,:) = T
      ocean%T_ocean_corr_ext( vi,:) = T

      ! Salinity
      ocean%S_ocean(          vi,:) = S
      ocean%S_ocean_ext(      vi,:) = S
      ocean%S_ocean_corr_ext( vi,:) = S

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_idealised_Reese2018_ANT

  ! == Uniform warm/cold ocean model (the old ANICE way)
  ! ====================================================

  SUBROUTINE run_ocean_model_uniform_warm_cold( mesh, ocean_matrix, time)
    ! Run the regional ocean model
    !
    ! Uniform warm/cold ocean model (the old ANICE way)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ocean_matrix_regional),    INTENT(INOUT) :: ocean_matrix
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_uniform_warm_cold'
    REAL(dp) :: dummy_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) 'run_ocean_model_uniform_warm_cold - ERROR: need to port this stuff from the ANICE_legacy BMB routine!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    dummy_dp = mesh%nV
    dummy_dp = ocean_matrix%applied%T_ocean_mean
    dummy_dp = time

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_uniform_warm_cold

  ! == Static present-day observed ocean
  ! ====================================

  SUBROUTINE run_ocean_model_PD_obs( mesh, ocean_matrix)
    ! Run the regional ocean model
    !
    ! Just use the present-day observed ocean data

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ocean_matrix_regional),    INTENT(INOUT) :: ocean_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_PD_obs'
    INTEGER                                            :: vi,k

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi2, mesh%vi2
    DO k = 1, C%nz_ocean
      ocean_matrix%applied%T_ocean(          vi,k) = ocean_matrix%PD_obs%T_ocean(          vi,k)
      ocean_matrix%applied%T_ocean_ext(      vi,k) = ocean_matrix%PD_obs%T_ocean_ext(      vi,k)
      ocean_matrix%applied%T_ocean_corr_ext( vi,k) = ocean_matrix%PD_obs%T_ocean_corr_ext( vi,k)
      ocean_matrix%applied%S_ocean(          vi,k) = ocean_matrix%PD_obs%S_ocean(          vi,k)
      ocean_matrix%applied%S_ocean_ext(      vi,k) = ocean_matrix%PD_obs%S_ocean_ext(      vi,k)
      ocean_matrix%applied%S_ocean_corr_ext( vi,k) = ocean_matrix%PD_obs%S_ocean_corr_ext( vi,k)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_PD_obs

  SUBROUTINE initialise_ocean_model_PD_obs_regional( region, ocean_matrix_global)
    ! Initialise the regional ocean model
    !
    ! Just use the present-day observed ocean data

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_ocean_matrix_global),      INTENT(IN)    :: ocean_matrix_global

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ocean_model_PD_obs_regional'
    INTEGER                                            :: vi,k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate all the snapshots
    CALL allocate_ocean_snapshot_regional( region%mesh, region%ocean_matrix%PD_obs,  name = 'PD_obs' )
    CALL allocate_ocean_snapshot_regional( region%mesh, region%ocean_matrix%applied, name = 'applied')

    ! Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid,
    ! extrapolate mapped ocean data to cover the entire 3D domain, and finally
    ! map to the actual ice model resolution
    CALL get_extrapolated_ocean_data( region, ocean_matrix_global%PD_obs, region%ocean_matrix%PD_obs, C%filename_PD_obs_ocean)

    DO k = 1, C%nz_ocean
    DO vi = region%mesh%vi1, region%mesh%vi2

      ! PD_obs doesn't have a bias-corrected version
      region%ocean_matrix%PD_obs%T_ocean_corr_ext(  vi,k) = region%ocean_matrix%PD_obs%T_ocean_ext( vi,k)
      region%ocean_matrix%PD_obs%S_ocean_corr_ext(  vi,k) = region%ocean_matrix%PD_obs%S_ocean_ext( vi,k)

      ! Initialise applied ocean forcing with present-day observations
      region%ocean_matrix%applied%T_ocean(          vi,k) = region%ocean_matrix%PD_obs%T_ocean(          vi,k)
      region%ocean_matrix%applied%T_ocean_ext(      vi,k) = region%ocean_matrix%PD_obs%T_ocean_ext(      vi,k)
      region%ocean_matrix%applied%T_ocean_corr_ext( vi,k) = region%ocean_matrix%PD_obs%T_ocean_corr_ext( vi,k)
      region%ocean_matrix%applied%S_ocean(          vi,k) = region%ocean_matrix%PD_obs%S_ocean(          vi,k)
      region%ocean_matrix%applied%S_ocean_ext(      vi,k) = region%ocean_matrix%PD_obs%S_ocean_ext(      vi,k)
      region%ocean_matrix%applied%S_ocean_corr_ext( vi,k) = region%ocean_matrix%PD_obs%S_ocean_corr_ext( vi,k)

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ocean_model_PD_obs_regional

! == Ocean matrix with warm and cold snapshots
! ============================================

  SUBROUTINE run_ocean_model_matrix_warm_cold( mesh, grid, ocean_matrix, climate_matrix, region_name, time)
    ! Run the regional ocean model
    !
    ! Run the warm/cold ocean matrix

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ocean_matrix_regional),    INTENT(INOUT) :: ocean_matrix
    TYPE(type_climate_matrix_regional),  INTENT(IN)    :: climate_matrix
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_matrix_warm_cold'
    INTEGER                                            :: vi,k
    REAL(dp)                                           :: a
    REAL(dp)                                           :: CO2
    REAL(dp)                                           :: w_CO2
    REAL(dp), DIMENSION(:    ), POINTER                :: w_ins, w_ins_smooth,  w_ice,  w_tot
    REAL(dp), DIMENSION(:,:  ), POINTER                :: w_tot_final
    INTEGER                                            :: ww_ins, ww_ins_smooth, ww_ice, ww_tot, ww_tot_final
    REAL(dp)                                           :: w_ins_av
    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, w_ins,        ww_ins                  )
    CALL allocate_shared_dp_1D( mesh%nV, w_ins_smooth, ww_ins_smooth           )
    CALL allocate_shared_dp_1D( mesh%nV, w_ice,        ww_ice                  )
    CALL allocate_shared_dp_1D( mesh%nV, w_tot,        ww_tot                  )
    CALL allocate_shared_dp_2D( mesh%nV, C%nz_ocean, w_tot_final,  ww_tot_final)

    ! Find CO2 interpolation weight (use either prescribed or modelled CO2)
    ! =====================================================================

    CALL update_CO2_at_model_time( time)

    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      CO2 = forcing%CO2_obs
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CO2 = forcing%CO2_mod
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CO2 = 0._dp
      WRITE(0,*) '  ERROR - run_ocean_model_matrix_warm_cold must only be called with the correct forcing method, check your code!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSE
      CO2 = 0._dp
      WRITE(0,*) '  ERROR - choice_forcing_method "', C%choice_forcing_method, '" not implemented in run_ocean_model_matrix_warm_cold!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    w_CO2 = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (CO2 - C%matrix_low_CO2_level) / (C%matrix_high_CO2_level - C%matrix_low_CO2_level) ))

    ! Find ice interpolation weight
    ! =============================

    ! Calculate weighting field
    DO vi = mesh%vi1, mesh%vi2
      w_ins( vi) = MAX( -w_cutoff, MIN( 1._dp + w_cutoff,  (   climate_matrix%applied%I_abs(  vi) -     climate_matrix%GCM_cold%I_abs( vi)) / &  ! Berends et al., 2018 - Eq. 3
                                                           (   climate_matrix%GCM_warm%I_abs( vi) -     climate_matrix%GCM_cold%I_abs( vi)) ))
    END DO
    CALL sync
    w_ins_av     = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (SUM(climate_matrix%applied%I_abs )     - SUM(climate_matrix%GCM_cold%I_abs)     ) / &
                                                          (SUM(climate_matrix%GCM_warm%I_abs)     - SUM(climate_matrix%GCM_cold%I_abs)     ) ))

    ! Smooth the weighting field
    w_ins_smooth( mesh%vi1:mesh%vi2) = w_ins( mesh%vi1:mesh%vi2)
    CALL smooth_Gaussian_2D( mesh, grid, w_ins_smooth, 200000._dp)

    ! Combine unsmoothed, smoothed, and regional average weighting fields (Berends et al., 2018, Eq. 4)
    w_ice( mesh%vi1:mesh%vi2) = (1._dp * w_ins_smooth( mesh%vi1:mesh%vi2) + 6._dp * w_ins_av) / 7._dp

    ! Combine weigths CO2 and ice
    ! ===========================

    IF         (region_name == 'NAM') THEN
      w_tot( mesh%vi1:mesh%vi2) = (C%ocean_matrix_CO2vsice_NAM * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_NAM) * w_ice( mesh%vi1:mesh%vi2))
    ELSEIF     (region_name == 'EAS') THEN
      w_tot( mesh%vi1:mesh%vi2) = (C%ocean_matrix_CO2vsice_EAS * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_EAS) * w_ice( mesh%vi1:mesh%vi2))
    ELSEIF     (region_name == 'GRL') THEN
      w_tot( mesh%vi1:mesh%vi2) = (C%ocean_matrix_CO2vsice_GRL * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_GRL) * w_ice( mesh%vi1:mesh%vi2))
    ELSEIF     (region_name == 'ANT') THEN
      w_tot( mesh%vi1:mesh%vi2) = (C%ocean_matrix_CO2vsice_ANT * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_ANT) * w_ice( mesh%vi1:mesh%vi2))
    END IF

    ! Update the history of the weighing fields
    ! =========================================

    ! 1st entry is the current value, 2nd is 1*dt_ocean ago, 3d is 2*dt_ocean ago, etc.
    ocean_matrix%applied%w_tot_history( mesh%vi1:mesh%vi2, 2:ocean_matrix%applied%nw_tot_history) = ocean_matrix%applied%w_tot_history( mesh%vi1:mesh%vi2, 1:ocean_matrix%applied%nw_tot_history-1)
    ocean_matrix%applied%w_tot_history( mesh%vi1:mesh%vi2, 1                                    ) =                      w_tot(         mesh%vi1:mesh%vi2                                         )

    ! Interpolate the GCM ocean snapshots
    ! =============================

    DO k = 1, C%nz_ocean
    DO vi = grid%i1, grid%i2
      a = ( ( C%z_ocean (k) / C%z_ocean (C%nz_ocean) ) * (ocean_matrix%applied%nw_tot_history-1) ) + 1

      w_tot_final (vi,k) = ( SUM(ocean_matrix%applied%w_tot_history(vi,1:FLOOR(a))) + ( (a - FLOOR(a)) * ocean_matrix%applied%w_tot_history(vi,CEILING(a)) ) ) * (1._dp / a)

      ocean_matrix%applied%T_ocean_corr_ext (vi,k) = (           w_tot_final( vi,k)  * ocean_matrix%GCM_warm%T_ocean_corr_ext( vi,k)  ) + &
                                                     (  (1._dp - w_tot_final( vi,k)) * ocean_matrix%GCM_cold%T_ocean_corr_ext( vi,k)  )
      ocean_matrix%applied%S_ocean_corr_ext (vi,k) = (           w_tot_final( vi,k)  * ocean_matrix%GCM_warm%S_ocean_corr_ext( vi,k)  ) + &
                                                     (  (1._dp - w_tot_final( vi,k)) * ocean_matrix%GCM_cold%S_ocean_corr_ext( vi,k)  )
    END DO
    END DO
    CALL sync

    CALL deallocate_shared( ww_ins)
    CALL deallocate_shared( ww_ins_smooth)
    CALL deallocate_shared( ww_ice)
    CALL deallocate_shared( ww_tot)
    CALL deallocate_shared( ww_tot_final)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_matrix_warm_cold

! == Initialise the global ocean matrix
  SUBROUTINE initialise_ocean_matrix_global( ocean_matrix)
    ! Initialise the global warm/cold ocean matrix

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ocean_matrix_global),      INTENT(INOUT) :: ocean_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ocean_matrix_global'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise the present-day observed global ocean (e.g. WOA18)
    CALL initialise_ocean_PD_obs_global(   ocean_matrix%PD_obs,   name = 'WOA18')

    ! Initialise the (GCM-modelled) global ocean snapshots
    CALL initialise_ocean_snapshot_global( ocean_matrix%GCM_PI,   name = 'ref_PI', nc_filename = C%filename_GCM_ocean_snapshot_PI,   CO2 = 280._dp,                 orbit_time = 0._dp                   )
    CALL initialise_ocean_snapshot_global( ocean_matrix%GCM_warm, name = 'warm',   nc_filename = C%filename_GCM_ocean_snapshot_warm, CO2 = C%matrix_high_CO2_level, orbit_time = C%matrix_warm_orbit_time)
    CALL initialise_ocean_snapshot_global( ocean_matrix%GCM_cold, name = 'cold',   nc_filename = C%filename_GCM_ocean_snapshot_cold, CO2 = C%matrix_low_CO2_level,  orbit_time = C%matrix_cold_orbit_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=39)

  END SUBROUTINE initialise_ocean_matrix_global

  SUBROUTINE initialise_ocean_PD_obs_global( PD_obs, name)
    ! Initialise the present-day observed global ocean snapshot

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_ocean_snapshot_global),   INTENT(INOUT) :: PD_obs
    CHARACTER(LEN=*),                   INTENT(IN)    :: name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'initialise_ocean_PD_obs_global'

    ! Add routine to path
    CALL init_routine( routine_name)

    PD_obs%name            = name
    PD_obs%netcdf%filename = C%filename_PD_obs_ocean

    ! Inquire if all required variables are present in the NetCDF file, and read the grid size.
    CALL allocate_shared_int_0D( PD_obs%nlon,         PD_obs%wnlon        )
    CALL allocate_shared_int_0D( PD_obs%nlat,         PD_obs%wnlat        )
    CALL allocate_shared_int_0D( PD_obs%nz_ocean_raw, PD_obs%wnz_ocean_raw)
    IF (par%master) CALL inquire_PD_obs_global_ocean_file( PD_obs)
    CALL sync

    ! Allocate memory
    CALL allocate_shared_dp_1D( PD_obs%nlon,                                   PD_obs%lon,         PD_obs%wlon        )
    CALL allocate_shared_dp_1D(              PD_obs%nlat,                      PD_obs%lat,         PD_obs%wlat        )
    CALL allocate_shared_dp_1D(                           PD_obs%nz_ocean_raw, PD_obs%z_ocean_raw, PD_obs%wz_ocean_raw)

    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, PD_obs%nz_ocean_raw, PD_obs%T_ocean_raw, PD_obs%wT_ocean_raw)
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, PD_obs%nz_ocean_Raw, PD_obs%S_ocean_raw, PD_obs%wS_ocean_raw)

    ! Read data from the NetCDF file
    IF (par%master) WRITE(0,*) '   Reading PD observed ocean data from file ', TRIM(PD_obs%netcdf%filename), '...'
    IF (par%master) CALL read_PD_obs_global_ocean_file( PD_obs)
    CALL sync

    ! Determine process domains
    CALL partition_list( PD_obs%nlon, par%i, par%n, PD_obs%i1, PD_obs%i2)

    ! Map the data to the desired vertical grid
    CALL map_global_ocean_data_to_UFEMISM_vertical_grid( PD_obs)

    ! Deallocate raw ocean data
    CALL deallocate_shared( PD_obs%wnz_ocean_raw)
    CALL deallocate_shared( PD_obs%wz_ocean_raw )
    CALL deallocate_shared( PD_obs%wT_ocean_raw )
    CALL deallocate_shared( PD_obs%wS_ocean_raw )

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=6)

  END SUBROUTINE initialise_ocean_PD_obs_global

  SUBROUTINE initialise_ocean_snapshot_global( snapshot, name, nc_filename, CO2, orbit_time)
    ! Initialise a (GCM-modelled) global ocean snapshot

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ocean_snapshot_global), INTENT(INOUT) :: snapshot
    CHARACTER(LEN=*),                 INTENT(IN)    :: name
    CHARACTER(LEN=*),                 INTENT(IN)    :: nc_filename
    REAL(dp),                         INTENT(IN)    :: CO2
    REAL(dp),                         INTENT(IN)    :: orbit_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                   :: routine_name = 'initialise_ocean_snapshot_global'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Metadata
    snapshot%name            = name
    snapshot%netcdf%filename = nc_filename

    ! General forcing info
    CALL allocate_shared_dp_0D( snapshot%CO2,        snapshot%wCO2       )
    CALL allocate_shared_dp_0D( snapshot%orbit_time, snapshot%worbit_time)
    CALL allocate_shared_dp_0D( snapshot%orbit_ecc,  snapshot%worbit_ecc )
    CALL allocate_shared_dp_0D( snapshot%orbit_obl,  snapshot%worbit_obl )
    CALL allocate_shared_dp_0D( snapshot%orbit_pre,  snapshot%worbit_pre )

    snapshot%CO2        = CO2
    snapshot%orbit_time = orbit_time

    ! Inquire if all required variables are present in the NetCDF file, and read the grid size.
    CALL allocate_shared_int_0D( snapshot%nlon,         snapshot%wnlon        )
    CALL allocate_shared_int_0D( snapshot%nlat,         snapshot%wnlat        )
    CALL allocate_shared_int_0D( snapshot%nz_ocean_raw, snapshot%wnz_ocean_raw)
    IF (par%master) CALL inquire_GCM_global_ocean_file( snapshot)
    CALL sync

    ! Allocate memory
    CALL allocate_shared_dp_1D( snapshot%nlon,                                       snapshot%lon,         snapshot%wlon        )
    CALL allocate_shared_dp_1D(                snapshot%nlat,                        snapshot%lat,         snapshot%wlat        )
    CALL allocate_shared_dp_1D(                               snapshot%nz_ocean_raw, snapshot%z_ocean_raw, snapshot%wz_ocean_raw)
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, snapshot%nz_ocean_raw, snapshot%T_ocean_raw, snapshot%wT_ocean_raw)
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, snapshot%nz_ocean_raw, snapshot%S_ocean_raw, snapshot%wS_ocean_raw)

    ! Read data from the NetCDF file
    IF (par%master) WRITE(0,*) '   Reading GCM ocean snapshot ', TRIM(snapshot%name), ' from file ', TRIM(snapshot%netcdf%filename), '...'
    IF (par%master) CALL read_GCM_global_ocean_file( snapshot)
    CALL sync

    ! Determine process domains
    CALL partition_list( snapshot%nlon, par%i, par%n, snapshot%i1, snapshot%i2)

    ! Map the data to the desired vertical grid
    CALL map_global_ocean_data_to_UFEMISM_vertical_grid( snapshot)

    ! Deallocate raw ocean data
    CALL deallocate_shared( snapshot%wnz_ocean_raw)
    CALL deallocate_shared( snapshot%wz_ocean_raw )
    CALL deallocate_shared( snapshot%wT_ocean_raw )
    CALL deallocate_shared( snapshot%wS_ocean_raw )

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=11)

  END SUBROUTINE initialise_ocean_snapshot_global

  ! == Initialising the regional ocean matrix
  SUBROUTINE initialise_ocean_matrix_regional( region, ocean_matrix_global)
    ! Initialise the regional warm/cold ocean matrix

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_ocean_matrix_global),      INTENT(IN)    :: ocean_matrix_global

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ocean_matrix_regional'
    INTEGER                                            :: vi,k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate all the snapshots
    CALL allocate_ocean_snapshot_regional( region%mesh, region%ocean_matrix%PD_obs,   name = 'PD_obs')
    CALL allocate_ocean_snapshot_regional( region%mesh, region%ocean_matrix%GCM_PI,   name = 'GCM_PI')
    CALL allocate_ocean_snapshot_regional( region%mesh, region%ocean_matrix%GCM_cold, name = 'GCM_cold')
    CALL allocate_ocean_snapshot_regional( region%mesh, region%ocean_matrix%GCM_warm, name = 'GCM_warm')
    CALL allocate_ocean_snapshot_regional( region%mesh, region%ocean_matrix%applied,  name = 'applied')

    ! Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid,
    ! extrapolate mapped ocean data to cover the entire 3D domain, and finally
    ! map to the actual ice model resolution
    CALL get_extrapolated_ocean_data( region, ocean_matrix_global%PD_obs,   region%ocean_matrix%PD_obs,   C%filename_PD_obs_ocean           )
    CALL get_extrapolated_ocean_data( region, ocean_matrix_global%GCM_PI,   region%ocean_matrix%GCM_PI,   C%filename_GCM_ocean_snapshot_PI  )
    CALL get_extrapolated_ocean_data( region, ocean_matrix_global%GCM_cold, region%ocean_matrix%GCM_cold, C%filename_GCM_ocean_snapshot_cold)
    CALL get_extrapolated_ocean_data( region, ocean_matrix_global%GCM_warm, region%ocean_matrix%GCM_warm, C%filename_GCM_ocean_snapshot_warm)

    ! Correct regional ocean data for GCM bias
    CALL correct_GCM_bias_ocean( region%mesh, region%ocean_matrix, region%ocean_matrix%GCM_warm)
    CALL correct_GCM_bias_ocean( region%mesh, region%ocean_matrix, region%ocean_matrix%GCM_cold)

    DO k = 1, C%nz_ocean
    DO vi = region%mesh%vi1, region%mesh%vi2

      ! PI and PD_obs don't have a bias-corrected version
      region%ocean_matrix%GCM_PI%T_ocean_corr_ext(  vi,k) = region%ocean_matrix%GCM_PI%T_ocean_ext(      vi,k)
      region%ocean_matrix%GCM_PI%S_ocean_corr_ext(  vi,k) = region%ocean_matrix%GCM_PI%S_ocean_ext(      vi,k)
      region%ocean_matrix%PD_obs%T_ocean_corr_ext(  vi,k) = region%ocean_matrix%PD_obs%T_ocean_ext(      vi,k)
      region%ocean_matrix%PD_obs%S_ocean_corr_ext(  vi,k) = region%ocean_matrix%PD_obs%S_ocean_ext(      vi,k)

      ! Initialise applied ocean forcing with present-day observations
      region%ocean_matrix%applied%T_ocean(          vi,k) = region%ocean_matrix%PD_obs%T_ocean(          vi,k)
      region%ocean_matrix%applied%T_ocean_ext(      vi,k) = region%ocean_matrix%PD_obs%T_ocean_ext(      vi,k)
      region%ocean_matrix%applied%T_ocean_corr_ext( vi,k) = region%ocean_matrix%PD_obs%T_ocean_corr_ext( vi,k)
      region%ocean_matrix%applied%S_ocean(          vi,k) = region%ocean_matrix%PD_obs%S_ocean(          vi,k)
      region%ocean_matrix%applied%S_ocean_ext(      vi,k) = region%ocean_matrix%PD_obs%S_ocean_ext(      vi,k)
      region%ocean_matrix%applied%S_ocean_corr_ext( vi,k) = region%ocean_matrix%PD_obs%S_ocean_corr_ext( vi,k)

    END DO
    END DO
    CALL sync

    ! Allocate memory for the weighing fields history, and initialise
    CALL allocate_shared_int_0D( region%ocean_matrix%applied%nw_tot_history, region%ocean_matrix%applied%wnw_tot_history)
    region%ocean_matrix%applied%nw_tot_history = CEILING( C%ocean_w_tot_hist_averaging_window / C%dt_ocean) + 1
    CALL allocate_shared_dp_2D(  region%mesh%nV, region%ocean_matrix%applied%nw_tot_history, region%ocean_matrix%applied%w_tot_history, region%ocean_matrix%applied%ww_tot_history)
    region%ocean_matrix%applied%w_tot_history = 0._dp ! Initiate at cold conditions

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=37)

  END SUBROUTINE initialise_ocean_matrix_regional

  SUBROUTINE allocate_ocean_snapshot_regional( mesh, ocean, name)
    ! Allocate shared memory for a regional ocean snapshot

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ocean_snapshot_regional),  INTENT(INOUT) :: ocean
    CHARACTER(LEN=*),                    INTENT(IN)    :: name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_ocean_snapshot_regional'

    ! Add routine to path
    CALL init_routine( routine_name)

    ocean%name = name

    ! Allocate shared memory
    CALL allocate_shared_dp_0D(                      ocean%T_ocean_mean,     ocean%wT_ocean_mean    )
    CALL allocate_shared_dp_2D( mesh%nV, C%nz_ocean, ocean%T_ocean,          ocean%wT_ocean         )
    CALL allocate_shared_dp_2D( mesh%nV, C%nz_ocean, ocean%S_ocean,          ocean%wS_ocean         )
    CALL allocate_shared_dp_2D( mesh%nV, C%nz_ocean, ocean%T_ocean_ext,      ocean%wT_ocean_ext     )
    CALL allocate_shared_dp_2D( mesh%nV, C%nz_ocean, ocean%S_ocean_ext,      ocean%wS_ocean_ext     )
    CALL allocate_shared_dp_2D( mesh%nV, C%nz_ocean, ocean%T_ocean_corr_ext, ocean%wT_ocean_corr_ext)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz_ocean, ocean%S_ocean_corr_ext, ocean%wS_ocean_corr_ext)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=7)

  END SUBROUTINE allocate_ocean_snapshot_regional

  SUBROUTINE correct_GCM_bias_ocean( mesh, ocean_matrix, ocean)
    ! Correct regional ocean data for GCM bias
    ! (must be done on regional grid, since the GCM grid and the World Ocean Atlas grid are generally not the same!)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ocean_matrix_regional),    INTENT(IN)    :: ocean_matrix
    TYPE(type_ocean_snapshot_regional),  INTENT(INOUT) :: ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'correct_GCM_bias_ocean'
    INTEGER                                            :: vi,k

    ! Add routine to path
    CALL init_routine( routine_name)

    DO k = 1, C%nz_ocean
    DO vi = mesh%vi1, mesh%vi2
      ocean%T_ocean_corr_ext( vi,k) = ocean%T_ocean_ext( vi,k) - (ocean_matrix%GCM_PI%T_ocean_ext( vi,k) - ocean_matrix%PD_obs%T_ocean_ext( vi,k))
      ocean%S_ocean_corr_ext( vi,k) = ocean%S_ocean_ext( vi,k) - (ocean_matrix%GCM_PI%S_ocean_ext( vi,k) - ocean_matrix%PD_obs%S_ocean_ext( vi,k))
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE correct_GCM_bias_ocean

! == Regrid 3-D ocean data fields in the vertical direction
! =========================================================

  SUBROUTINE map_global_ocean_data_to_UFEMISM_vertical_grid( ocean_glob)
    ! Map global 3-D ocean temperature and salinity from whatever vertical grid
    ! it has been provided on to the vertical grid used by UFEMISM

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_ocean_snapshot_global), INTENT(INOUT) :: ocean_glob

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                   :: routine_name = 'map_global_ocean_data_to_UFEMISM_vertical_grid'
    INTEGER                                         :: i,j,k
    INTEGER,  DIMENSION(:    ), ALLOCATABLE         :: z_mask_old, z_mask_new
    REAL(dp)                                        :: z_floor
    REAL(dp)                                        :: NaN

    ! Add routine to path
    CALL init_routine( routine_name)

    ! A trick
    NaN = -1._dp
    NaN = SQRT( NaN)

    ! Allocate shared memory for the global 3-D ocean temperature and salinity
    CALL allocate_shared_dp_3D( ocean_glob%nlon, ocean_glob%nlat, C%nz_ocean, ocean_glob%T_ocean, ocean_glob%wT_ocean)
    CALL allocate_shared_dp_3D( ocean_glob%nlon, ocean_glob%nlat, C%nz_ocean, ocean_glob%S_ocean, ocean_glob%wS_ocean)

    ! Use "masked" 2-nd order conservative 1-D remapping
    ALLOCATE( z_mask_old( ocean_glob%nz_ocean_raw))
    ALLOCATE( z_mask_new( C%nz_ocean             ))

    DO i = ocean_glob%i1, ocean_glob%i2
    DO j = 1, ocean_glob%nlat

      ! Determine local depth of the ocean floor, fill in both data masks
      IF (ocean_glob%T_ocean_raw( i,j,ocean_glob%nz_ocean_raw) == ocean_glob%T_ocean_raw( i,j,ocean_glob%nz_ocean_raw)) THEN
        ! Ocean floor lies below the vertical limit of the provided data
        z_mask_old = 1
        z_floor = ocean_glob%z_ocean_raw( ocean_glob%nz_ocean_raw) + (ocean_glob%z_ocean_raw( 2) - ocean_glob%z_ocean_raw( 1))
      ELSEIF (ocean_glob%T_ocean_raw( i,j,1) /= ocean_glob%T_ocean_raw( i,j,1)) THEN
        ! This grid cell isn't ocean at all
        z_mask_old = 0
        z_floor    = 0._dp
      ELSE
        z_mask_old = 1
        k = ocean_glob%nz_ocean_raw
        DO WHILE (ocean_glob%T_ocean_raw( i,j,k) /= ocean_glob%T_ocean_raw( i,j,k))
          z_mask_old( k) = 0
          z_floor = ocean_glob%z_ocean_raw( k)
          k = k - 1
        END DO
      END IF

      z_mask_new = 0
      DO k = 1, C%nz_ocean
        IF (C%z_ocean( k) < z_floor) z_mask_new = 1
      END DO

      ! Regrid vertical column
      CALL remap_cons_2nd_order_1D( ocean_glob%z_ocean_raw, z_mask_old, ocean_glob%T_ocean_raw( i,j,:), C%z_ocean, z_mask_new, ocean_glob%T_ocean( i,j,:))
      CALL remap_cons_2nd_order_1D( ocean_glob%z_ocean_raw, z_mask_old, ocean_glob%S_ocean_raw( i,j,:), C%z_ocean, z_mask_new, ocean_glob%S_ocean( i,j,:))

      ! Fill masked values with NaN
      DO k = 1, C%nz_ocean
        IF (z_mask_new( k) == 0) THEN
          ocean_glob%T_ocean( i,j,k) = NaN
          ocean_glob%S_ocean( i,j,k) = NaN
        END IF
      END DO

    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( z_mask_old)
    DEALLOCATE( z_mask_new)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=2)

  END SUBROUTINE map_global_ocean_data_to_UFEMISM_vertical_grid

  SUBROUTINE initialise_ocean_vertical_grid
    ! Set up the vertical grid used for ocean data

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ocean_vertical_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_ocean_vertical_grid == 'regular') THEN
      CALL initialise_ocean_vertical_grid_regular
    ELSE
      CALL crash('unknown choice_ocean_vertical_grid "' // TRIM( C%choice_ocean_vertical_grid) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ocean_vertical_grid
  SUBROUTINE initialise_ocean_vertical_grid_regular
    ! Set up the vertical grid used for ocean data - regular grid

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ocean_vertical_grid_regular'
    INTEGER                                       :: k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine the number of vertical layers to be used
    C%nz_ocean = 1 + FLOOR( C%ocean_vertical_grid_max_depth / C%ocean_regular_grid_dz)

    ! Allocate memory
    ALLOCATE( C%z_ocean( C%nz_ocean))

    ! Fill in the values
    DO k = 1, C%nz_ocean
      C%z_ocean( k) = REAL(k-1,dp) * C%ocean_regular_grid_dz
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ocean_vertical_grid_regular

! == Extrapolating incomplete ocean data into the complete 3-D model domain
! =========================================================================

  SUBROUTINE get_extrapolated_ocean_data( region, ocean_glob, ocean_reg, filename_ocean_glob)
    ! Check if extrapolated ocean files for the current ice model setting exist. If so,
    ! read those. If not, perform the extrapolation and save the results to a new netCDF file.

    ! When creating a set of extrapolated files, a header file is created that describes
    ! the ice model settings for which those files were created. We check all existing header
    ! files, if any of them match the current settings, we read the extrapolated files listed
    ! there. If none of them match, we create a set of extrapolated files (and the accompanying
    ! header file) from scratch.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_ocean_snapshot_global),    INTENT(IN)    :: ocean_glob
    TYPE(type_ocean_snapshot_regional),  INTENT(INOUT) :: ocean_reg
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename_ocean_glob

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'get_extrapolated_ocean_data'
    LOGICAL                                            :: foundmatch
    TYPE(type_highres_ocean_data)                      :: hires
    INTEGER                                            :: i,j,n
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ===== Initialise the high-resolution extrapolated ocean data =====
    ! ==================================================================

    ! If a valid preprocessed file exists, read data from there. If not, perform
    ! the preprocessing and save the result to a file to save on future work

    ! First, check if any existing header matches the current ice model set-up.
    CALL check_for_matching_ocean_header( region, filename_ocean_glob, foundmatch, ocean_reg%hires_ocean_foldername)

    IF (foundmatch) THEN
      IF (par%master) WRITE(0,*) '   Found valid extrapolated ocean data in folder "', TRIM( ocean_reg%hires_ocean_foldername), '"'
      CALL get_hires_ocean_data_from_file( region%mesh, hires, ocean_reg%hires_ocean_foldername)
    ELSE
      ! No header fitting the current ice model set-up was found. Create a new one describing
      ! the current set-up, and generate extrapolated ocean data files from scratch.
      IF (par%master) WRITE(0,*) '   Creating new extrapolated ocean data in folder "', TRIM( ocean_reg%hires_ocean_foldername), '"'
      CALL map_and_extrapolate_hires_ocean_data( region, ocean_glob, hires)
      CALL write_hires_extrapolated_ocean_data_to_file( hires, filename_ocean_glob, ocean_reg%hires_ocean_foldername)
    END IF ! IF (.NOT. foundmatch) THEN

    ! ===== Map extrapolated data from the high-resolution grid to the actual ice-model mesh =====
    ! ============================================================================================

    IF (par%master) WRITE(0,*) '    Mapping high-resolution extrapolated ocean data unto the ice-model mesh...'

    ! Tolerance; points lying within this distance of each other are treated as identical
    CALL allocate_shared_dp_0D( hires%grid%tol_dist, hires%grid%wtol_dist)
    IF (par%master) hires%grid%tol_dist = ((hires%grid%xmax - hires%grid%xmin) + (hires%grid%ymax - hires%grid%ymin)) * tol / 2._dp

    ! Set up grid-to-vector translation tables
    CALL allocate_shared_int_0D(                               hires%grid%n,    hires%grid%wn)
    IF (par%master) hires%grid%n  = hires%grid%nx * hires%grid%ny
    CALL sync
    CALL allocate_shared_int_2D( hires%grid%nx, hires%grid%ny, hires%grid%ij2n, hires%grid%wij2n)
    CALL allocate_shared_int_2D( hires%grid%n , 2            , hires%grid%n2ij, hires%grid%wn2ij)
    IF (par%master) THEN
      n = 0
      DO i = 1, hires%grid%nx
        IF (MOD(i,2) == 1) THEN
          DO j = 1, hires%grid%ny
            n = n+1
            hires%grid%ij2n( i,j) = n
            hires%grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = hires%grid%ny, 1, -1
            n = n+1
            hires%grid%ij2n( i,j) = n
            hires%grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO
    END IF
    CALL sync

    ! Map high-resolution ocean data to UFEMISM's mesh
    CALL calc_remapping_operator_grid2mesh( hires%grid, region%mesh)
    CALL map_grid2mesh_3D( hires%grid, region%mesh, hires%T_ocean, ocean_reg%T_ocean_ext)
    CALL map_grid2mesh_3D( hires%grid, region%mesh, hires%S_ocean, ocean_reg%S_ocean_ext)
    CALL deallocate_remapping_operators_mesh2grid( hires%grid)

    ! Clean up after yourself
    CALL deallocate_shared( hires%grid%wnx             )
    CALL deallocate_shared( hires%grid%wny             )
    CALL deallocate_shared( hires%grid%wx              )
    CALL deallocate_shared( hires%grid%wy              )
    CALL deallocate_shared( hires%grid%wxmin           )
    CALL deallocate_shared( hires%grid%wxmax           )
    CALL deallocate_shared( hires%grid%wymin           )
    CALL deallocate_shared( hires%grid%wymax           )
    CALL deallocate_shared( hires%grid%wdx             )
    CALL deallocate_shared( hires%grid%wlambda_M       )
    CALL deallocate_shared( hires%grid%wphi_M          )
    CALL deallocate_shared( hires%grid%walpha_stereo   )
    CALL deallocate_shared( hires%grid%wlat            )
    CALL deallocate_shared( hires%grid%wlon            )
    CALL deallocate_shared( hires%grid%wtol_dist       )
    CALL deallocate_shared( hires%grid%wn              )
    CALL deallocate_shared( hires%grid%wij2n           )
    CALL deallocate_shared( hires%grid%wn2ij           )
    CALL deallocate_shared( hires%wT_ocean             )
    CALL deallocate_shared( hires%wS_ocean             )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_extrapolated_ocean_data

  SUBROUTINE get_hires_ocean_data_from_file( mesh, hires, hires_foldername)
    ! Read high-resolution extrapolated ocean data from an external file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    CHARACTER(LEN=256),                  INTENT(IN)    :: hires_foldername

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'get_hires_ocean_data_from_file'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if the NetCDF file has all the required dimensions and variables
    CALL allocate_shared_int_0D( hires%grid%nx,  hires%grid%wnx)
    CALL allocate_shared_int_0D( hires%grid%ny,  hires%grid%wny)
    IF (par%master) THEN
      hires%netcdf%filename = TRIM( hires_foldername)//'/extrapolated_ocean_data.nc'
      CALL inquire_extrapolated_ocean_file( hires)
    END IF
    CALL sync

    ! Allocate shared memory for x,y and the actual data
    CALL allocate_shared_dp_1D( hires%grid%nx,                hires%grid%x,  hires%grid%wx )
    CALL allocate_shared_dp_1D(                hires%grid%ny, hires%grid%y,  hires%grid%wy )
    CALL allocate_shared_dp_3D( hires%grid%nx, hires%grid%ny, C%nz_ocean, hires%T_ocean, hires%wT_ocean)
    CALL allocate_shared_dp_3D( hires%grid%nx, hires%grid%ny, C%nz_ocean, hires%S_ocean, hires%wS_ocean)

    ! Read the data from the NetCDF file
    IF (par%master) THEN
      WRITE(0,*) '    Reading high-resolution extrapolated ocean data from file "', TRIM( hires%netcdf%filename), '"...'
      CALL read_extrapolated_ocean_file( hires)
    END IF
    CALL sync

    ! Allocate shared memory for other grid parameters
    CALL allocate_shared_dp_0D(  hires%grid%dx,           hires%grid%wdx          )
    CALL allocate_shared_dp_0D(  hires%grid%xmin,         hires%grid%wxmin        )
    CALL allocate_shared_dp_0D(  hires%grid%xmax,         hires%grid%wxmax        )
    CALL allocate_shared_dp_0D(  hires%grid%ymin,         hires%grid%wymin        )
    CALL allocate_shared_dp_0D(  hires%grid%ymax,         hires%grid%wymax        )
    CALL allocate_shared_dp_0D(  hires%grid%lambda_M,     hires%grid%wlambda_M    )
    CALL allocate_shared_dp_0D(  hires%grid%phi_M,        hires%grid%wphi_M       )
    CALL allocate_shared_dp_0D(  hires%grid%alpha_stereo, hires%grid%walpha_stereo)

    ! Polar stereographic projection parameters and resolution
    IF (par%master) THEN
      ! Projection parameters are of course identical to those used for this ice model region
      hires%grid%lambda_M     = mesh%lambda_M
      hires%grid%phi_M        = mesh%phi_M
      hires%grid%alpha_stereo = mesh%alpha_stereo
      ! But the resolution is different
      hires%grid%dx           = hires%grid%x( 2) - hires%grid%x( 1)
    END IF
    CALL sync

    ! Assign range to each processor
    CALL partition_list( hires%grid%nx, par%i, par%n, hires%grid%i1, hires%grid%i2)
    CALL partition_list( hires%grid%ny, par%i, par%n, hires%grid%j1, hires%grid%j2)

    ! Lat,lon coordinates
    CALL allocate_shared_dp_2D( hires%grid%nx, hires%grid%ny, hires%grid%lat, hires%grid%wlat)
    CALL allocate_shared_dp_2D( hires%grid%nx, hires%grid%ny, hires%grid%lon, hires%grid%wlon)

    DO j = 1, hires%grid%ny
    DO i = hires%grid%i1, hires%grid%i2
      CALL inverse_oblique_sg_projection( hires%grid%x( i), hires%grid%y( j), hires%grid%lambda_M, hires%grid%phi_M, hires%grid%alpha_stereo, hires%grid%lon( i,j), hires%grid%lat( i,j))
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=16)

  END SUBROUTINE get_hires_ocean_data_from_file

  SUBROUTINE write_hires_extrapolated_ocean_data_to_file( hires, filename_ocean_glob, hires_foldername)
    ! 1. Create a new folder inside the "extrapolated_ocean_files" folder
    ! 2. Create a header file inside this new folder listing the current model settings
    ! 3. Create a NetCDF file inside this new folder
    ! 4. Write the high-resolution extrapolated ocean data to this NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename_ocean_glob
    CHARACTER(LEN=256),                  INTENT(IN)    :: hires_foldername

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_hires_extrapolated_ocean_data_to_file'
    CHARACTER(LEN=256)                                 :: hires_ocean_filename

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new folder where the new extrapolated ocean file will be stored.
    IF (par%master) CALL system('mkdir ' // TRIM(hires_foldername))

    ! Create a header file describing the current ice-model set-up.
    CALL write_ocean_header( hires, filename_ocean_glob, hires_foldername)

    ! Create a NetCDF file and write data to it
    IF (par%master) hires_ocean_filename = TRIM(hires_foldername)//'/extrapolated_ocean_data.nc'
    IF (par%master) WRITE(0,*) '    Writing extrapolated ocean data to file "', TRIM(hires_ocean_filename), '"...'
    CALL create_extrapolated_ocean_file( hires, hires_ocean_filename)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_hires_extrapolated_ocean_data_to_file

  SUBROUTINE check_for_matching_ocean_header( region, filename_ocean_glob, foundmatch, hires_foldername)
    ! Inspect all the folder inside the "extrapolated_ocean_files" folder, read their header files,
    ! and see if any of them match the current model settings. If so, return the name of the folder
    ! where the matching header was found. If not, return the name of the folder where a new header
    ! should be created.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(IN)    :: region
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename_ocean_glob
    LOGICAL,                             INTENT(OUT)   :: foundmatch
    CHARACTER(LEN=256),                  INTENT(OUT)   :: hires_foldername

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_for_matching_ocean_header'
    INTEGER                                            :: cerr, ierr
    CHARACTER(LEN=256)                                 :: header_filename
    LOGICAL                                            :: header_exists
    INTEGER                                            :: folder_i

    CHARACTER(LEN=256)                                 :: original_ocean_filename_read
    CHARACTER(LEN=256)                                 :: choice_ocean_vertical_grid_read
    INTEGER                                            :: nz_ocean_read
    REAL(dp)                                           :: ocean_vertical_grid_max_depth_read
    REAL(dp)                                           :: ocean_extrap_res_read
    REAL(dp)                                           :: ocean_extrap_Gauss_sigma_read
    REAL(dp)                                           :: lambda_M_read
    REAL(dp)                                           :: phi_M_read
    REAL(dp)                                           :: alpha_stereo_read

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      folder_i = 1
      DO WHILE (folder_i < 1000)

        ! Generate a foldername to inspect
        IF     (folder_i < 10)   THEN
          WRITE( hires_foldername,'(A,A,A,A,I1)') TRIM(C%ocean_extrap_dir), '/', region%name, '_00', folder_i
        ELSEIF (folder_i < 100)  THEN
          WRITE( hires_foldername,'(A,A,A,A,I2)') TRIM(C%ocean_extrap_dir), '/', region%name, '_0',  folder_i
        ELSEIF (folder_i < 1000) THEN
          WRITE( hires_foldername,'(A,A,A,A,I3)') TRIM(C%ocean_extrap_dir), '/', region%name, '_',   folder_i
        ELSE
          IF (par%master) WRITE(0,*) ' ERROR: tried a thousand folders!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF

        ! Check if a header in this folder exists. If not, then we've inspected all existing headers
        ! without finding the good one, so we must generate the extrapolated ocean files from scratch.
        header_filename = TRIM( hires_foldername)//'/'//'header.txt'
        INQUIRE( FILE = header_filename, EXIST = header_exists)

        IF (.NOT. header_exists) THEN
          ! No more headers exist to be inspected. Return foundmatch = .FALSE.

          foundmatch = .FALSE.
          EXIT

        ELSE
          ! If the header exists, read it and see if it fits the current ice-model set-up.

          CALL read_ocean_header( &
            header_filename,                    &
            original_ocean_filename_read,       &
            choice_ocean_vertical_grid_read,    &
            nz_ocean_read,                      &
            ocean_vertical_grid_max_depth_read, &
            ocean_extrap_res_read,              &
            ocean_extrap_Gauss_sigma_read,      &
            lambda_M_read,                      &
            phi_M_read,                         &
            alpha_stereo_read)

          IF ( TRIM(original_ocean_filename_read)    == TRIM(filename_ocean_glob)             .AND. &
               TRIM(choice_ocean_vertical_grid_read) == TRIM(choice_ocean_vertical_grid_read) .AND. &
               nz_ocean_read                         == C%nz_ocean                            .AND. &
               ocean_vertical_grid_max_depth_read    == C%ocean_vertical_grid_max_depth       .AND. &
               ocean_extrap_res_read                 == C%ocean_extrap_res                    .AND. &
               ocean_extrap_Gauss_sigma_read         == C%ocean_extrap_Gauss_sigma            .AND. &
               lambda_M_read                         == region%mesh%lambda_M                  .AND. &
               phi_M_read                            == region%mesh%phi_M                     .AND. &
               alpha_stereo_read                     == region%mesh%alpha_stereo) THEN
            ! This header matches the current model set-up!

            foundmatch = .TRUE.
            EXIT

          ELSE
            ! This header doesn't match the current model set-up. Try the next one.
            folder_i = folder_i + 1
          END IF

        END IF

      END DO ! DO WHILE (header_i < 1000)

    END IF ! IF (par%master) THEN
    CALL sync

    CALL MPI_BCAST( hires_foldername,       256,   MPI_CHAR,    0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( foundmatch,               1,   MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( header_filename,          256, MPI_CHAR,    0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( folder_i,                 1,   MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_for_matching_ocean_header
  SUBROUTINE read_ocean_header( &
            header_filename,                    &
            original_ocean_filename,       &
            choice_ocean_vertical_grid,    &
            nz_ocean,                      &
            ocean_vertical_grid_max_depth, &
            ocean_extrap_res,              &
            ocean_extrap_Gauss_sigma,      &
            lambda_M,                      &
            phi_M,                         &
            alpha_stereo                   )
    ! Read a header file listing the model settings that were used to create a high-resolution extrapolated ocean data file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: header_filename
    CHARACTER(LEN=256),                  INTENT(OUT)   :: original_ocean_filename
    CHARACTER(LEN=256),                  INTENT(OUT)   :: choice_ocean_vertical_grid
    INTEGER,                             INTENT(OUT)   :: nz_ocean
    REAL(dp),                            INTENT(OUT)   :: ocean_vertical_grid_max_depth
    REAL(dp),                            INTENT(OUT)   :: ocean_extrap_res
    REAL(dp),                            INTENT(OUT)   :: ocean_extrap_Gauss_sigma
    REAL(dp),                            INTENT(OUT)   :: lambda_M
    REAL(dp),                            INTENT(OUT)   :: phi_M
    REAL(dp),                            INTENT(OUT)   :: alpha_stereo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_ocean_header'
    INTEGER                                            :: ios, cerr, ierr

    ! The NAMELIST that's used to read the external header file.
    NAMELIST /HEADER/original_ocean_filename,             &
                     choice_ocean_vertical_grid,          &
                     nz_ocean,                            &
                     ocean_vertical_grid_max_depth,       &
                     ocean_extrap_res,                    &
                     ocean_extrap_Gauss_sigma,            &
                     lambda_M,                            &
                     phi_M,                               &
                     alpha_stereo

    ! Add routine to path
    CALL init_routine( routine_name)

    OPEN( UNIT = 29, FILE = TRIM(header_filename), STATUS='OLD', ACTION='READ', iostat=ios)
    IF (ios /= 0) THEN
      WRITE(0,*) ' ERROR: could not open ""', TRIM(header_filename), '"'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! In the following statement the entire configuration file is read, using the namelist (NML=HEADER)
    READ(  UNIT = 29, NML = HEADER, IOSTAT = ios)
    CLOSE( UNIT = 29)

    IF (ios /= 0) THEN
      WRITE(0,*) ' ERROR: could not read ""', TRIM(header_filename), '"'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_ocean_header
  SUBROUTINE write_ocean_header( hires, filename_ocean_glob, hires_foldername)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename_ocean_glob
    CHARACTER(LEN=256),                  INTENT(IN)    :: hires_foldername

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_ocean_header'
    INTEGER, DIMENSION(8)                              :: datevec
    CHARACTER(LEN=256)                                 :: header_filename

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Let the Master do the work
    IF (par%master) THEN

      header_filename = TRIM( hires_foldername)//'/header.txt'

      OPEN( UNIT = 1337, FILE = header_filename, STATUS = 'NEW')

      CALL date_and_time( VALUES = datevec)

      WRITE(UNIT = 1337, FMT = '(A)') '&HEADER'
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A,I4,A,I2,A,I2)') '! Icemodel-ocean header file, created on ', datevec(1), '-', datevec(2), '-', datevec(3)
      WRITE(UNIT = 1337, FMT = '(A)') '!'
      WRITE(UNIT = 1337, FMT = '(A)') '! This header describes the icemodel set-up that was used to created this'
      WRITE(UNIT = 1337, FMT = '(A)') '! extrapolated ocean data file. Since creating these is computationally intensive,'
      WRITE(UNIT = 1337, FMT = '(A)') '! reading them from files is preferred. These header files make ice model'
      WRITE(UNIT = 1337, FMT = '(A)') '! a bit more flexible when using different input files for the four model regions.'
      WRITE(UNIT = 1337, FMT = '(A)') ''

      WRITE(UNIT = 1337, FMT = '(A)')       '! The original global ocean file that was extrapolated'
      WRITE(UNIT = 1337, FMT = '(A,A,A)')   'original_ocean_filename       = ''', TRIM(filename_ocean_glob), ''''
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)')       '! The vertical grid the ocean data was projected to'
      WRITE(UNIT = 1337, FMT = '(A,A,A)')   'choice_ocean_vertical_grid    = ''', TRIM(C%choice_ocean_vertical_grid), ''''
      WRITE(UNIT = 1337, FMT = '(A,I5)')    'nz_ocean                      = ', C%nz_ocean
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'ocean_vertical_grid_max_depth = ', C%ocean_vertical_grid_max_depth
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)')       '! Resolution and Gaussian smoothing radius used for the high-resolution extrapolation'
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'ocean_extrap_res              = ', C%ocean_extrap_res
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'ocean_extrap_Gauss_sigma      = ', C%ocean_extrap_Gauss_sigma
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)')       '! Parameters of the high-resolution grid'
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'lambda_M                      = ', hires%grid%lambda_M
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'phi_M                         = ', hires%grid%phi_M
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'alpha_stereo                  = ', hires%grid%alpha_stereo

      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)') '/'

      CLOSE(UNIT = 1337)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_ocean_header

  SUBROUTINE map_and_extrapolate_hires_ocean_data( region, ocean_glob, hires)
    ! Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid,
    ! extrapolate mapped ocean data to cover the entire 3D domain, and finally
    ! map to the actual ice model resolution.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),          INTENT(INOUT) :: region
    TYPE(type_ocean_snapshot_global), INTENT(IN)    :: ocean_glob
    TYPE(type_highres_ocean_data),    INTENT(INOUT) :: hires

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                   :: routine_name = 'map_and_extrapolate_hires_ocean_data'
    INTEGER                                         :: vi,i,j
    REAL(dp), DIMENSION(:    ), POINTER             ::  basin_ID_dp_lores
    REAL(dp), DIMENSION(:,:  ), POINTER             ::  basin_ID_dp_hires,  basin_ID_dp_hires_ext
    INTEGER                                         :: wbasin_ID_dp_lores, wbasin_ID_dp_hires, wbasin_ID_dp_hires_ext
    INTEGER                                         :: ii,jj,n
    LOGICAL                                         :: foundit
    REAL(dp), PARAMETER                             :: tol = 1E-9_dp
    CHARACTER(LEN=256)                              :: choice_basin_scheme
    CHARACTER(LEN=256)                              :: filename_basins

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ===== Read the high-resolution geometry data and generate the hi-res grid based on that=====
    ! ============================================================================================

    ! Determine which file to use for this region
    IF     (region%name == 'NAM') THEN
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_NAM
      choice_basin_scheme       = C%choice_basin_scheme_NAM
      filename_basins           = C%filename_basins_NAM
    ELSEIF (region%name == 'EAS') THEN
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_EAS
      choice_basin_scheme       = C%choice_basin_scheme_EAS
      filename_basins           = C%filename_basins_EAS
    ELSEIF (region%name == 'GRL') THEN
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_GRL
      choice_basin_scheme       = C%choice_basin_scheme_GRL
      filename_basins           = C%filename_basins_GRL
    ELSEIF (region%name == 'ANT') THEN
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_ANT
      choice_basin_scheme       = C%choice_basin_scheme_ANT
      filename_basins           = C%filename_basins_ANT
    END IF

    ! Check if the NetCDF file has all the required dimensions and variables
    CALL allocate_shared_int_0D( hires%grid%nx, hires%grid%wnx)
    CALL allocate_shared_int_0D( hires%grid%ny, hires%grid%wny)
    IF (par%master) THEN
      CALL inquire_hires_geometry_file( hires)
    END IF
    CALL sync

    ! Allocate shared memory for x,y and the actual data
    CALL allocate_shared_dp_1D( hires%grid%nx,                hires%grid%x, hires%grid%wx)
    CALL allocate_shared_dp_1D(                hires%grid%ny, hires%grid%y, hires%grid%wy)
    CALL allocate_shared_dp_2D( hires%grid%nx, hires%grid%ny, hires%Hi,     hires%wHi    )
    CALL allocate_shared_dp_2D( hires%grid%nx, hires%grid%ny, hires%Hb,     hires%wHb    )

    ! Read the data from the NetCDF file
    IF (par%master) THEN
      WRITE(0,*) '    Reading high-resolution geometry for ocean extrapolation from file "', TRIM( hires%netcdf_geo%filename), '"...'
      CALL read_hires_geometry_file( hires)
    END IF
    CALL sync

    ! Allocate shared memory for other grid parameters
    CALL allocate_shared_dp_0D(  hires%grid%dx,           hires%grid%wdx          )
    CALL allocate_shared_dp_0D(  hires%grid%xmin,         hires%grid%wxmin        )
    CALL allocate_shared_dp_0D(  hires%grid%xmax,         hires%grid%wxmax        )
    CALL allocate_shared_dp_0D(  hires%grid%ymin,         hires%grid%wymin        )
    CALL allocate_shared_dp_0D(  hires%grid%ymax,         hires%grid%wymax        )
    CALL allocate_shared_dp_0D(  hires%grid%lambda_M,     hires%grid%wlambda_M    )
    CALL allocate_shared_dp_0D(  hires%grid%phi_M,        hires%grid%wphi_M       )
    CALL allocate_shared_dp_0D(  hires%grid%alpha_stereo, hires%grid%walpha_stereo)

    ! Polar stereographic projection parameters and resolution
    IF (par%master) THEN
      ! Projection parameters are of course identical to those used for this ice model region
      hires%grid%lambda_M     = region%mesh%lambda_M
      hires%grid%phi_M        = region%mesh%phi_M
      hires%grid%alpha_stereo = region%mesh%alpha_stereo
      ! But the resolution is different
      hires%grid%dx           = hires%grid%x( 2) - hires%grid%x( 1)
      hires%grid%xmin         = hires%grid%x( 1            )
      hires%grid%xmax         = hires%grid%x( hires%grid%nx)
      hires%grid%ymin         = hires%grid%y( 1            )
      hires%grid%ymax         = hires%grid%y( hires%grid%ny)
      ! Check if this is the resolution we want
      IF (hires%grid%dx /= C%ocean_extrap_res) THEN
        CALL crash('high-resolution geometry file "' // TRIM( hires%netcdf_geo%filename) // '" has a different resolution from C%ocean_extrap_res = {dp_01}', dp_01 = C%ocean_extrap_res)
      END IF
    END IF
    CALL sync

    ! Tolerance; points lying within this distance of each other are treated as identical
    CALL allocate_shared_dp_0D( hires%grid%tol_dist, hires%grid%wtol_dist)
    IF (par%master) hires%grid%tol_dist = ((hires%grid%xmax - hires%grid%xmin) + (hires%grid%ymax - hires%grid%ymin)) * tol / 2._dp

    ! Set up grid-to-vector translation tables
    CALL allocate_shared_int_0D( hires%grid%n, hires%grid%wn)
    IF (par%master) hires%grid%n = hires%grid%nx * hires%grid%ny
    CALL sync

    CALL allocate_shared_int_2D( hires%grid%nx, hires%grid%ny, hires%grid%ij2n, hires%grid%wij2n)
    CALL allocate_shared_int_2D( hires%grid%n , 2,             hires%grid%n2ij, hires%grid%wn2ij)
    IF (par%master) THEN
      n = 0
      DO i = 1, hires%grid%nx
        IF (MOD(i,2) == 1) THEN
          DO j = 1, hires%grid%ny
            n = n+1
            hires%grid%ij2n( i,j) = n
            hires%grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = hires%grid%ny, 1, -1
            n = n+1
            hires%grid%ij2n( i,j) = n
            hires%grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO
    END IF
    CALL sync

    ! Assign range to each processor
    CALL partition_list( hires%grid%nx, par%i, par%n, hires%grid%i1, hires%grid%i2)
    CALL partition_list( hires%grid%ny, par%i, par%n, hires%grid%j1, hires%grid%j2)

    ! Lat,lon coordinates
    CALL allocate_shared_dp_2D( hires%grid%nx, hires%grid%ny, hires%grid%lat, hires%grid%wlat)
    CALL allocate_shared_dp_2D( hires%grid%nx, hires%grid%ny, hires%grid%lon, hires%grid%wlon)

    DO i = hires%grid%i1, hires%grid%i2
    DO j = 1, hires%grid%ny
      CALL inverse_oblique_sg_projection( hires%grid%x( i), hires%grid%y( j), hires%grid%lambda_M, hires%grid%phi_M, hires%grid%alpha_stereo, hires%grid%lon( i,j), hires%grid%lat( i,j))
    END DO
    END DO
    CALL sync

    ! ===== Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid =====
    ! ================================================================================================

    IF (par%master) WRITE(0,'(A,F4.1,A)') '     Mapping ocean data from the global lat/lon-grid to the ', hires%grid%dx / 1000._dp, ' km regional x/y-grid...'

    ! Allocate shared memory for high-resolution ocean data
    CALL allocate_shared_dp_3D( hires%grid%nx, hires%grid%ny, C%nz_ocean, hires%T_ocean, hires%wT_ocean)
    CALL allocate_shared_dp_3D( hires%grid%nx, hires%grid%ny, C%nz_ocean, hires%S_ocean, hires%wS_ocean)

    ! Map the data from the global lon/lat-grid to the high-resolution regional x/y-grid
    CALL map_glob_to_grid_3D( ocean_glob%nlat, ocean_glob%nlon, ocean_glob%lat, ocean_glob%lon, hires%grid, ocean_glob%T_ocean, hires%T_ocean)
    CALL map_glob_to_grid_3D( ocean_glob%nlat, ocean_glob%nlon, ocean_glob%lat, ocean_glob%lon, hires%grid, ocean_glob%S_ocean, hires%S_ocean)

    ! ===== Perform the extrapolation on the high-resolution grid =====
    ! =================================================================

    IF (par%master) WRITE(0,'(A,F4.1,A)') '     Defining ice basins on the ', hires%grid%dx / 1000._dp, ' km regional x/y-grid...'

    ! Allocate shared memory for ice basins on the high-resolution grid
    CALL allocate_shared_int_2D( hires%grid%nx, hires%grid%ny, hires%basin_ID, hires%wbasin_ID)
    CALL allocate_shared_int_0D(                               hires%nbasins,  hires%wnbasins )

    IF (par%master) hires%nbasins = region%ice%nbasins
    CALL sync

    ! Instead of doing the "proper" basin definition on high resolution (which is insanely slow),
    ! just downscale the basin ID field from the ice model (using some tricks to get accurate values near the boundaries)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( region%mesh%nV,                 basin_ID_dp_lores,     wbasin_ID_dp_lores    )
    CALL allocate_shared_dp_2D( hires%grid%nx,  hires%grid%ny,  basin_ID_dp_hires,     wbasin_ID_dp_hires    )
    CALL allocate_shared_dp_2D( hires%grid%nx,  hires%grid%ny,  basin_ID_dp_hires_ext, wbasin_ID_dp_hires_ext)

    IF     (choice_basin_scheme == 'none') THEN
      ! No basins are defined (i.e. the whole region is one big, big basin)

      DO j = 1, hires%grid%ny
      DO i = hires%grid%i1, hires%grid%i2
        hires%basin_ID( i,j) = 1
      END DO
      END DO
      CALL sync

    ELSEIF (choice_basin_scheme == 'file') THEN

      ! Convert basin ID field to double precision (for remapping)
      DO vi = region%mesh%vi1, region%mesh%vi2
        basin_ID_dp_lores( vi) = REAL( region%ice%basin_ID( vi), dp)
      END DO
      CALL sync

      ! Map double-precision basin ID from ice-model mesh to high-resolution grid
      CALL calc_remapping_operator_mesh2grid( region%mesh, hires%grid)
      CALL map_mesh2grid_2D( region%mesh, hires%grid, basin_ID_dp_lores, basin_ID_dp_hires)
      CALL deallocate_remapping_operators_mesh2grid( hires%grid)

      ! Remove all near-boundary cells
      DO j = 1, hires%grid%ny
      DO i = hires%grid%i1, hires%grid%i2
        IF (MODULO( basin_ID_dp_hires( i,j), 1._dp) > 0.01_dp) THEN
          basin_ID_dp_hires( i,j) = -1._dp
        END IF
      END DO
      END DO
      CALL sync

      ! For those, use extrapolation instead
      basin_ID_dp_hires_ext( hires%grid%i1:hires%grid%i2,:) = basin_ID_dp_hires( hires%grid%i1:hires%grid%i2,:)
      CALL sync

      DO j = 1, hires%grid%ny
      DO i = hires%grid%i1, hires%grid%i2
        IF (basin_ID_dp_hires_ext( i,j) == -1._dp) THEN

            n = 0
            foundit = .FALSE.
            DO WHILE (.NOT. foundit)

              n = n+1

              ! Take the value of the nearest non-boundary cell
              DO jj = MAX(1,j-n), MIN(hires%grid%ny,j+n)
              DO ii = MAX(1,i-n), MIN(hires%grid%nx,i+n)
                IF (basin_ID_dp_hires( ii,jj) > -1._dp) THEN
                  basin_ID_dp_hires_ext( i,j) = basin_ID_dp_hires( ii,jj)
                  foundit = .TRUE.
                  EXIT
                END IF
              END DO
              IF (foundit) EXIT
              END DO

              ! Safety
              IF (n > MAX(hires%grid%nx, hires%grid%ny)) THEN
                WRITE(0,*) 'map_and_extrapolate_ocean_data - ERROR: basin ID downscaling got stuck!'
                CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
              END IF

            END DO ! DO WHILE (.NOT. foundit)

          END IF ! IF (basin_ID_dp_hires_ext( i,j) == -1._dp) THEN
        END DO
        END DO
      CALL sync

      ! Convert hi-resolution basin ID field back to integer precision
      DO j = 1, hires%grid%ny
      DO i = hires%grid%i1, hires%grid%i2
        hires%basin_ID( i,j) = NINT( basin_ID_dp_hires_ext( i,j))
      END DO
      END DO
      CALL sync

    ELSE
      CALL crash('unknown choice_basin_scheme "' // TRIM(choice_basin_scheme) // '"!')
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wbasin_ID_dp_lores)
    CALL deallocate_shared( wbasin_ID_dp_hires)
    CALL deallocate_shared( wbasin_ID_dp_hires_ext)

    ! Perform the extrapolation on the high-resolution grid
    IF (par%master) WRITE(0,'(A,F4.1,A)') '     Performing ocean data extrapolation on the ', hires%grid%dx / 1000._dp, ' km regional x/y-grid...'
    CALL extend_regional_ocean_data_to_cover_domain( hires)

    ! Clean up fields that were needed only for the extrapolation
    CALL deallocate_shared( hires%wHi      )
    CALL deallocate_shared( hires%wHb      )
    CALL deallocate_shared( hires%wbasin_ID)
    CALL deallocate_shared( hires%wnbasins )

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=20)

  END SUBROUTINE map_and_extrapolate_hires_ocean_data
  SUBROUTINE extend_regional_ocean_data_to_cover_domain( hires)
    ! Extend global ocean data over the whole grid, based on the procedure outlined in
    ! Jourdain, N. C., Asay-Davis, X., Hattermann, T., Straneo, F., Seroussi, H., Little, C. M., & Nowicki, S. (2020).
    ! A protocol for calculating basal melt rates in the ISMIP6 Antarctic ice sheet projections. The Cryosphere, 14(9), 3111-3134.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extend_regional_ocean_data_to_cover_domain'
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: NaN, Hs, z_bedrock, z_icebase, z
    INTEGER,  DIMENSION(:,:,:), POINTER                ::  mask_wetdry,  mask_hasdata
    INTEGER                                            :: wmask_wetdry, wmask_hasdata
    INTEGER                                            :: k1,k2,bi
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE            :: mask, mask_filled
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_T, d_S
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  T_ocean_ext,  S_ocean_ext
    INTEGER                                            :: wT_ocean_ext, wS_ocean_ext

    LOGICAL,  PARAMETER                                :: verbose = .FALSE.

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Useful
    NaN = -1._dp
    NaN = SQRT( NaN)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( hires%grid%nx, hires%grid%ny, C%nz_ocean, T_ocean_ext, wT_ocean_ext)
    CALL allocate_shared_dp_3D( hires%grid%nx, hires%grid%ny, C%nz_ocean, S_ocean_ext, wS_ocean_ext)

    ! Initialise the extrapolated product with the provided ocean data
    DO k = 1, C%nz_ocean
    DO j = 1, hires%grid%ny
    DO i = hires%grid%i1, hires%grid%i2
      T_ocean_ext( i,j,k) = hires%T_ocean( i,j,k)
      S_ocean_ext( i,j,k) = hires%S_ocean( i,j,k)
    END DO
    END DO
    END DO
    CALL sync

    ! Define the two masks needed for the four extrapolation steps:
    !
    !  - mask_wetdry:
    !      1 = actually    wet (i.e. open ocean, sub-shelf cavity, above sea floor and beneath ice base)
    !      2 = potentially wet (i.e. grounded marine ice, above sea floor)
    !      3 =             dry (i.e. beneath bedrock                     )
    !
    !  - mask_hasdata:
    !      0 = has no data
    !      1 = has data provided
    !      2 = has data extrapolated

    CALL allocate_shared_int_3D( hires%grid%nx, hires%grid%ny, C%nz_ocean, mask_wetdry,  wmask_wetdry )
    CALL allocate_shared_int_3D( hires%grid%nx, hires%grid%ny, C%nz_ocean, mask_hasdata, wmask_hasdata)

    DO j = 1, hires%grid%ny
    DO i = hires%grid%i1, hires%grid%i2

      Hs = surface_elevation( hires%Hi( i,j), hires%Hb( i,j), 0._dp)
      z_bedrock = hires%Hb( i,j)
      z_icebase = Hs - hires%Hi( i,j)

      DO k = 1, C%nz_ocean

        z = -C%z_ocean( k)

        ! mask_wetdry
        IF (z < z_bedrock) THEN
          ! This 3D hires%grid box is beneath the bedrock surface, so dry
          mask_wetdry( i,j,k) = 3
        ELSE
          ! This 3D hires%grid box is above the bedrock surface so at least potentially wet
          IF (z < z_icebase) THEN
            ! This 3D hires%grid box is above the bedrock surface and below the ice base, so it is actually wet
            mask_wetdry( i,j,k) = 1
          ELSE
            ! This 3D hires%grid box is above the bedrock surface and above the ice base, so it is potentially wet (i.e. inside grounded marine ice)
            mask_wetdry( i,j,k) = 2
          END IF
        END IF

        ! mask_hasdata
        IF (hires%T_ocean( i,j,k) /= hires%T_ocean( i,j,k)) THEN
          ! This 3D hires%grid box has no data (yet)
          mask_hasdata( i,j,k) = 0
        ELSE
          ! Data is already provided for this 3D hires%grid box
          mask_hasdata( i,j,k) = 1
        END IF

      END DO

    END DO
    END DO
    CALL sync

    ! ================================================================
    ! ===== Step 1: horizontal extrapolation into shelf cavities =====
    ! ================================================================

    IF (par%master .AND. verbose) WRITE(0,*) '    extend_regional_ocean_data_to_cover_domain - step 1'

    ! Here, we start with the ocean data as provided (i.e. only for open ocean), and
    ! perform a horizontal extrapolation (so for each vertical layer separately) into
    ! the shelf cavities. Only "actually wet" 3D grid boxes are allowed to be filled,
    ! so the fill is limited by both bedrock sills and grounded ice.

    ! Allocate memory for the mask and data field of a single extrapolation step
    ALLOCATE( mask(        hires%grid%nx, hires%grid%ny))
    ALLOCATE( mask_filled( hires%grid%nx, hires%grid%ny))
    ALLOCATE( d_T(         hires%grid%nx, hires%grid%ny))
    ALLOCATE( d_S(         hires%grid%nx, hires%grid%ny))

    ! Parallelised by partitioning the vertical domain
    CALL partition_list( C%nz_ocean, par%i, par%n, k1, k2)

    DO k = k1, k2

      ! Extrapolate per basin
      DO bi = 1, hires%nbasins

        IF (verbose) WRITE(0,'(A,I2,A,I3,A,I3,A,I3,A,I3)') '        process ', par%i, ': vertical layer ', k, '/', C%nz_ocean, ', basin ', bi, '/', hires%nbasins

        ! Define the mask and initial data fields for this particular flood-fill
        ! (i.e. this vertical layer and this basin)
        mask        = 0
        d_T         = NaN
        d_S         = NaN
        DO j = 1, hires%grid%ny
        DO i = 1, hires%grid%nx
          IF (hires%basin_ID( i,j) == bi) THEN
            IF (mask_hasdata( i,j,k) == 1) THEN
              ! This is where the source data comes from
              mask( i,j) = 2
              d_T(  i,j) = T_ocean_ext( i,j,k)
              d_S(  i,j) = S_ocean_ext( i,j,k)
            ELSEIF (mask_hasdata( i,j,k) == 0 .AND. mask_wetdry( i,j,k) == 1) THEN
              ! This is where we're supposed to fill it in
              mask( i,j) = 1
            END IF
          END IF
        END DO
        END DO

        ! Perform the flood-fill-based Gaussian extrapolation
        CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_T, C%ocean_extrap_Gauss_sigma, mask_filled)
        CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_S, C%ocean_extrap_Gauss_sigma, mask_filled)

        ! Copy extrapolated data to the data structure

        DO j = 1, hires%grid%ny
        DO i = 1, hires%grid%nx
          IF (mask_filled( i,j) == 1) THEN
            T_ocean_ext(  i,j,k) = d_T( i,j)
            S_ocean_ext(  i,j,k) = d_S( i,j)
            mask_hasdata( i,j,k) = 2
          END IF
        END DO
        END DO

      END DO ! DO bi = 1, ice%nbasins

    END DO ! DO k = k1, k2
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( mask       )
    DEALLOCATE( mask_filled)
    DEALLOCATE( d_T        )
    DEALLOCATE( d_S        )

    ! ===========================================================================
    ! ===== Step 2: vertical extrapolation into sill-blocked shelf cavities =====
    ! ===========================================================================

    IF (par%master .AND. verbose) WRITE(0,*) '    extend_regional_ocean_data_to_cover_domain - step 2'

    ! Here, we start with the ocean data that has been horizontally extrapolated into
    ! the shelf cavities, allowing for bedrock topography to block the fill, so that
    ! for example the lower parts of the Filchner-Ronne and Ross cavities have not yet
    ! been filled. We now extrapolate the data vertically from the filled parts to
    ! fill those parts of the cavities. Barring any really weird geometry, the entire
    ! cavities will now be filled.

    DO j = 1, hires%grid%ny
    DO i = hires%grid%i1, hires%grid%i2

      ! Move down through the vertical column
      DO k = 2, C%nz_ocean
        ! If this grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        IF (mask_wetdry( i,j,k) == 1 .AND. mask_hasdata( i,j,k) == 0) THEN
          ! This 3D grid box is wet but has no data
          IF (mask_hasdata( i,j,k-1) == 1 .OR. mask_hasdata( i,j,k-1) == 2) THEN
            ! The one above it has data; copy data
            mask_hasdata( i,j,k) = 2
            T_ocean_ext( i,j,k) = T_ocean_ext( i,j,k-1)
            S_ocean_ext( i,j,k) = S_ocean_ext( i,j,k-1)
          END IF
        END IF
      END DO

    END DO
    END DO
    CALL sync

    ! ===============================================================
    ! ===== Step 3: vertical extrapolation into ice and bedrock =====
    ! ===============================================================

    IF (par%master .AND. verbose) WRITE(0,*) '    extend_regional_ocean_data_to_cover_domain - step 3'

    ! Extrapolate data vertically into 3D grid boxes that are occupied by ice
    ! or bedrock (since they might turn into ocean at some point during a simulation)

    DO j = 1, hires%grid%ny
    DO i = hires%grid%i1, hires%grid%i2

      ! Move down through the vertical column
      DO k = 2, C%nz_ocean
        ! If this grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        IF (mask_hasdata( i,j,k) == 0) THEN
          ! This 3D grid box is wet but has no data
          IF (mask_hasdata( i,j,k-1) == 1 .OR. mask_hasdata( i,j,k-1) == 2) THEN
            ! The one above it has data; copy data
            mask_hasdata( i,j,k) = 2
            T_ocean_ext( i,j,k) = T_ocean_ext( i,j,k-1)
            S_ocean_ext( i,j,k) = S_ocean_ext( i,j,k-1)
          END IF
        END IF
      END DO

      ! Move up through the vertical column
      DO k = C%nz_ocean-1, 1, -1
        ! If this grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        IF (mask_hasdata( i,j,k) == 0) THEN
          ! This 3D grid box is wet but has no data
          IF (mask_hasdata( i,j,k+1) == 1 .OR. mask_hasdata( i,j,k+1) == 2) THEN
            ! The one above it has data; copy data
            mask_hasdata( i,j,k) = 2
            T_ocean_ext(  i,j,k) = T_ocean_ext( i,j,k+1)
            S_ocean_ext(  i,j,k) = S_ocean_ext( i,j,k+1)
          END IF
        END IF
      END DO

    END DO
    END DO
    CALL sync

    ! =================================================================
    ! ===== Step 4: horizontal extrapolation into ice and bedrock =====
    ! =================================================================

    IF (par%master .AND. verbose) WRITE(0,*) '    extend_regional_ocean_data_to_cover_domain - step 4'

    ! In the last step, extrapolate data horizontally into 3D
    ! grid boxes that are occupied by ice or bedrock

    ! Allocate memory for the mask and data field of a single extrapolation step
    ALLOCATE( mask(        hires%grid%nx, hires%grid%ny))
    ALLOCATE( mask_filled( hires%grid%nx, hires%grid%ny))
    ALLOCATE( d_T(         hires%grid%nx, hires%grid%ny))
    ALLOCATE( d_S(         hires%grid%nx, hires%grid%ny))

    ! Parallelised by partitioning the vertical domain
    CALL partition_list( C%nz_ocean, par%i, par%n, k1, k2)

    DO k = k1, k2

      ! Extrapolate per basin
      DO bi = 1, hires%nbasins

        IF (verbose) WRITE(0,'(A,I2,A,I3,A,I3,A,I3,A,I3)') '        process ', par%i, ': vertical layer ', k, '/', C%nz_ocean, ', basin ', bi, '/', hires%nbasins

        ! Define the mask and initial data fields for this particular flood-fill
        ! (i.e. this vertical layer and this basin)
        mask        = 0
        d_T         = NaN
        d_S         = NaN
        DO j = 1, hires%grid%ny
        DO i = 1, hires%grid%nx
          IF (hires%basin_ID( i,j) == bi) THEN
            IF (mask_hasdata( i,j,k) == 1 .OR. mask_hasdata( i,j,k) == 2) THEN
              ! This is where the source data comes from
              mask( i,j) = 2
              d_T(  i,j) = T_ocean_ext( i,j,k)
              d_S(  i,j) = S_ocean_ext( i,j,k)
            ELSEIF (mask_hasdata( i,j,k) == 0) THEN
              ! This is where we're supposed to fill it in
              mask( i,j) = 1
            END IF
          END IF
        END DO
        END DO

        ! Perform the flood-fill-based Gaussian extrapolation

        CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_T, C%ocean_extrap_Gauss_sigma, mask_filled)
        CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_S, C%ocean_extrap_Gauss_sigma, mask_filled)

        ! Copy extrapolated data to the data structure
        DO j = 1, hires%grid%ny
        DO i = 1, hires%grid%nx
          IF (mask_filled( i,j) == 1) THEN
            T_ocean_ext(  i,j,k) = d_T( i,j)
            S_ocean_ext(  i,j,k) = d_S( i,j)
            mask_hasdata( i,j,k) = 2
          END IF
        END DO
        END DO

      END DO ! DO bi = 1, ice%nbasins

      ! One more pass without considering basins (sometimes, a handful of isolated "basin enclaves" can
      ! occur at high resolution, which will not be filled when using the basin-constrained flood-fill)

      ! Define the mask and initial data fields for this particular flood-fill
      ! (i.e. this vertical layer and this basin)
      mask        = 0
      d_T         = NaN
      d_S         = NaN

      DO j = 1, hires%grid%ny
      DO i = 1, hires%grid%nx
        IF (mask_hasdata( i,j,k) == 1 .OR. mask_hasdata( i,j,k) == 2) THEN
          ! This is where the source data comes from
          mask( i,j) = 2
          d_T(  i,j) = T_ocean_ext( i,j,k)
          d_S(  i,j) = S_ocean_ext( i,j,k)
        ELSEIF (mask_hasdata( i,j,k) == 0) THEN
          ! This is where we're supposed to fill it in
          mask( i,j) = 1
        END IF
      END DO
      END DO

      ! Perform the flood-fill-based Gaussian extrapolation
      CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_T, C%ocean_extrap_Gauss_sigma, mask_filled)
      CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_S, C%ocean_extrap_Gauss_sigma, mask_filled)

      ! Copy extrapolated data to the data structure
      DO j = 1, hires%grid%ny
      DO i = 1, hires%grid%nx
        IF (mask_filled( i,j) == 1) THEN
          T_ocean_ext(  i,j,k) = d_T( i,j)
          S_ocean_ext(  i,j,k) = d_S( i,j)
          mask_hasdata( i,j,k) = 2
        END IF
      END DO
      END DO

      ! Check if all pixels have now been filled
      DO j = 1, hires%grid%ny
      DO i = 1, hires%grid%nx
        IF (T_ocean_ext( i,j,k) /= T_ocean_ext( i,j,k)) THEN
          WRITE(0,*) '    extend_regional_ocean_data_to_cover_domain - ERROR: unfilled pixels remains at [i,j] = [', i,',',j,']'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END DO
      END DO

    END DO ! DO k = k1, k2
    CALL sync

    ! Copy data to the hires structure
    DO k = 1, C%nz_ocean
    DO j = 1, hires%grid%ny
    DO i = hires%grid%i1, hires%grid%i2
      hires%T_ocean( i,j,k) = T_ocean_ext( i,j,k)
      hires%S_ocean( i,j,k) = S_ocean_ext( i,j,k)
    END DO
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( mask       )
    DEALLOCATE( mask_filled)
    DEALLOCATE( d_T        )
    DEALLOCATE( d_S        )
    CALL deallocate_shared( wmask_wetdry )
    CALL deallocate_shared( wmask_hasdata)
    CALL deallocate_shared( wT_ocean_ext )
    CALL deallocate_shared( wS_ocean_ext )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extend_regional_ocean_data_to_cover_domain

! == Remap ocean fields after mesh update
! =======================================

  SUBROUTINE remap_ocean_model( mesh_old, mesh_new, map, ocean_matrix)
    ! Reallocate all the data fields (no remapping needed, instead we just run
    ! the climate model immediately after a mesh update)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_old
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ocean_matrix_regional),    INTENT(INOUT) :: ocean_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_ocean_model'
    INTEGER                                            :: vi,k

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '  Reallocating regional ocean model...'

    IF     (C%choice_ocean_model == 'none') THEN
      ! No ocean data is used at all

    ELSEIF (C%choice_ocean_model == 'idealised') THEN
      ! Only reallocate the "applied" regional snapshot

      CALL reallocate_subocean( mesh_new, ocean_matrix%applied)

    ELSEIF (C%choice_ocean_model == 'uniform_warm_cold') THEN
      ! Only uses the 0-D T_ocean_mean variables. Nothing to do here.

    ELSEIF (C%choice_ocean_model == 'PD_obs') THEN

      ! Allocate the PD_obs and applied snapshots
      CALL reallocate_subocean( mesh_new, ocean_matrix%PD_obs )
      CALL reallocate_subocean( mesh_new, ocean_matrix%applied)

      ! Map high-resolution extrapolated ocean data to the new model mesh
      CALL remap_ocean_data( mesh_new, ocean_matrix%PD_obs)

      DO k = 1, C%nz_ocean
      DO vi = mesh_new%vi1, mesh_new%vi2

        ! PD_obs doesn't have a bias-corrected version
        ocean_matrix%PD_obs%T_ocean_corr_ext(  vi,k) = ocean_matrix%PD_obs%T_ocean_ext( vi,k)
        ocean_matrix%PD_obs%S_ocean_corr_ext(  vi,k) = ocean_matrix%PD_obs%S_ocean_ext( vi,k)

        ! Initialise applied ocean forcing with present-day observations
        ocean_matrix%applied%T_ocean(          vi,k) = ocean_matrix%PD_obs%T_ocean(          vi,k)
        ocean_matrix%applied%T_ocean_ext(      vi,k) = ocean_matrix%PD_obs%T_ocean_ext(      vi,k)
        ocean_matrix%applied%T_ocean_corr_ext( vi,k) = ocean_matrix%PD_obs%T_ocean_corr_ext( vi,k)
        ocean_matrix%applied%S_ocean(          vi,k) = ocean_matrix%PD_obs%S_ocean(          vi,k)
        ocean_matrix%applied%S_ocean_ext(      vi,k) = ocean_matrix%PD_obs%S_ocean_ext(      vi,k)
        ocean_matrix%applied%S_ocean_corr_ext( vi,k) = ocean_matrix%PD_obs%S_ocean_corr_ext( vi,k)

      END DO
      END DO
      CALL sync

    ELSEIF (C%choice_ocean_model == 'matrix_warm_cold') THEN

      ! Allocate all snapshots
      CALL reallocate_subocean( mesh_new, ocean_matrix%PD_obs  )
      CALL reallocate_subocean( mesh_new, ocean_matrix%GCM_PI  )
      CALL reallocate_subocean( mesh_new, ocean_matrix%GCM_cold)
      CALL reallocate_subocean( mesh_new, ocean_matrix%GCM_warm)
      CALL reallocate_subocean( mesh_new, ocean_matrix%applied )

      ! Map high-resolution extrapolated ocean data to the new model mesh
      CALL remap_ocean_data( mesh_new, ocean_matrix%PD_obs  )
      CALL remap_ocean_data( mesh_new, ocean_matrix%GCM_PI  )
      CALL remap_ocean_data( mesh_new, ocean_matrix%GCM_cold)
      CALL remap_ocean_data( mesh_new, ocean_matrix%GCM_warm)

      ! Correct regional ocean data for GCM bias
      CALL correct_GCM_bias_ocean( mesh_new, ocean_matrix, ocean_matrix%GCM_warm)
      CALL correct_GCM_bias_ocean( mesh_new, ocean_matrix, ocean_matrix%GCM_cold)

      DO k = 1, C%nz_ocean
      DO vi = mesh_new%vi1, mesh_new%vi2

        ! PI and PD_obs don't have a bias-corrected version
        ocean_matrix%GCM_PI%T_ocean_corr_ext(  vi,k) = ocean_matrix%GCM_PI%T_ocean_ext(      vi,k)
        ocean_matrix%GCM_PI%S_ocean_corr_ext(  vi,k) = ocean_matrix%GCM_PI%S_ocean_ext(      vi,k)
        ocean_matrix%PD_obs%T_ocean_corr_ext(  vi,k) = ocean_matrix%PD_obs%T_ocean_ext(      vi,k)
        ocean_matrix%PD_obs%S_ocean_corr_ext(  vi,k) = ocean_matrix%PD_obs%S_ocean_ext(      vi,k)

        ! Initialise applied ocean forcing with present-day observations
        ocean_matrix%applied%T_ocean(          vi,k) = ocean_matrix%PD_obs%T_ocean(          vi,k)
        ocean_matrix%applied%T_ocean_ext(      vi,k) = ocean_matrix%PD_obs%T_ocean_ext(      vi,k)
        ocean_matrix%applied%T_ocean_corr_ext( vi,k) = ocean_matrix%PD_obs%T_ocean_corr_ext( vi,k)
        ocean_matrix%applied%S_ocean(          vi,k) = ocean_matrix%PD_obs%S_ocean(          vi,k)
        ocean_matrix%applied%S_ocean_ext(      vi,k) = ocean_matrix%PD_obs%S_ocean_ext(      vi,k)
        ocean_matrix%applied%S_ocean_corr_ext( vi,k) = ocean_matrix%PD_obs%S_ocean_corr_ext( vi,k)

      END DO
      END DO
      CALL sync

      ! The weighing fields history must be remapped
      CALL remap_field_dp_3D( mesh_old, mesh_new, map, ocean_matrix%applied%w_tot_history, ocean_matrix%applied%ww_tot_history, 'trilin')

    ELSE
      IF (par%master) WRITE(0,*) 'remap_ocean_model - ERROR: unknown choice_ocean_model "', TRIM(C%choice_ocean_model), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ocean_model

  SUBROUTINE reallocate_subocean( mesh_new, ocean_reg)
    ! Reallocate data fields of a regional ocean snapshot after a mesh update

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ocean_snapshot_regional),  INTENT(INOUT) :: ocean_reg

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'reallocate_subocean'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Reallocate shared memory
    CALL reallocate_shared_dp_2D( mesh_new%nV, C%nz_ocean, ocean_reg%T_ocean,          ocean_reg%wT_ocean         )
    CALL reallocate_shared_dp_2D( mesh_new%nV, C%nz_ocean, ocean_reg%S_ocean,          ocean_reg%wS_ocean         )
    CALL reallocate_shared_dp_2D( mesh_new%nV, C%nz_ocean, ocean_reg%T_ocean_ext,      ocean_reg%wT_ocean_ext     )
    CALL reallocate_shared_dp_2D( mesh_new%nV, C%nz_ocean, ocean_reg%S_ocean_ext,      ocean_reg%wS_ocean_ext     )
    CALL reallocate_shared_dp_2D( mesh_new%nV, C%nz_ocean, ocean_reg%T_ocean_corr_ext, ocean_reg%wT_ocean_corr_ext)
    CALL reallocate_shared_dp_2D( mesh_new%nV, C%nz_ocean, ocean_reg%S_ocean_corr_ext, ocean_reg%wS_ocean_corr_ext)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE reallocate_subocean

  SUBROUTINE remap_ocean_data( mesh_new, ocean_reg)
    ! Check if extrapolated ocean files for the current ice model
    ! setting exist, based on the header found/created during
    ! initialistion. If so, read those. If not, throw an error, as
    ! those files should be there.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_new
    TYPE(type_ocean_snapshot_regional),  INTENT(INOUT) :: ocean_reg

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_ocean_data'
    LOGICAL                                            :: folder_exists
    TYPE(type_highres_ocean_data)                      :: hires
    INTEGER                                            :: i,j,n
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ===== Read the high-resolution extrapolated ocean data =====
    ! ============================================================

    INQUIRE( FILE = ocean_reg%hires_ocean_foldername, EXIST = folder_exists)
    IF (folder_exists) THEN
      IF (par%master) WRITE(0,*) '   Found extrapolated ocean data in folder "', TRIM( ocean_reg%hires_ocean_foldername), '"'
      CALL get_hires_ocean_data_from_file( mesh_new, hires, ocean_reg%hires_ocean_foldername)
    ELSE
      ! No header fitting the current ice model set-up was found.
      ! This should not occur, so throw an error and abort.
      IF (par%master) WRITE(0,*) 'remap_ocean_data - ERROR: no valid extrapolated ocean data found during mesh update!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! ===== Map extrapolated data from the high-resolution grid to the actual ice-model mesh =====
    ! ============================================================================================

    IF (par%master) WRITE(0,*) '    Mapping high-resolution extrapolated ocean data onto the new mesh...'

    ! Tolerance; points lying within this distance of each other are treated as identical
    CALL allocate_shared_dp_0D( hires%grid%tol_dist, hires%grid%wtol_dist)
    IF (par%master) hires%grid%tol_dist = ((hires%grid%xmax - hires%grid%xmin) + (hires%grid%ymax - hires%grid%ymin)) * tol / 2._dp

    ! Set up grid-to-vector translation tables
    CALL allocate_shared_int_0D(                               hires%grid%n,    hires%grid%wn)
    IF (par%master) hires%grid%n  = hires%grid%nx * hires%grid%ny
    CALL sync
    CALL allocate_shared_int_2D( hires%grid%nx, hires%grid%ny, hires%grid%ij2n, hires%grid%wij2n)
    CALL allocate_shared_int_2D( hires%grid%n , 2            , hires%grid%n2ij, hires%grid%wn2ij)
    IF (par%master) THEN
      n = 0
      DO i = 1, hires%grid%nx
        IF (MOD(i,2) == 1) THEN
          DO j = 1, hires%grid%ny
            n = n+1
            hires%grid%ij2n( i,j) = n
            hires%grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = hires%grid%ny, 1, -1
            n = n+1
            hires%grid%ij2n( i,j) = n
            hires%grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO
    END IF
    CALL sync

    ! Map high-resolution ocean data to UFEMISM's mesh
    CALL calc_remapping_operator_grid2mesh( hires%grid, mesh_new)
    CALL map_grid2mesh_3D( hires%grid, mesh_new, hires%T_ocean, ocean_reg%T_ocean_ext)
    CALL map_grid2mesh_3D( hires%grid, mesh_new, hires%S_ocean, ocean_reg%S_ocean_ext)
    CALL deallocate_remapping_operators_mesh2grid( hires%grid)

    ! Clean up after yourself
    CALL deallocate_shared( hires%grid%wnx             )
    CALL deallocate_shared( hires%grid%wny             )
    CALL deallocate_shared( hires%grid%wx              )
    CALL deallocate_shared( hires%grid%wy              )
    CALL deallocate_shared( hires%grid%wxmin           )
    CALL deallocate_shared( hires%grid%wxmax           )
    CALL deallocate_shared( hires%grid%wymin           )
    CALL deallocate_shared( hires%grid%wymax           )
    CALL deallocate_shared( hires%grid%wdx             )
    CALL deallocate_shared( hires%grid%wlambda_M       )
    CALL deallocate_shared( hires%grid%wphi_M          )
    CALL deallocate_shared( hires%grid%walpha_stereo   )
    CALL deallocate_shared( hires%grid%wlat            )
    CALL deallocate_shared( hires%grid%wlon            )
    CALL deallocate_shared( hires%grid%wtol_dist       )
    CALL deallocate_shared( hires%grid%wn              )
    CALL deallocate_shared( hires%grid%wij2n           )
    CALL deallocate_shared( hires%grid%wn2ij           )
    CALL deallocate_shared( hires%wT_ocean             )
    CALL deallocate_shared( hires%wS_ocean             )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ocean_data

END MODULE ocean_module