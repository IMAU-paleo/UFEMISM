MODULE BMB_module
  ! Contains all the routines for calculating the basal mass balance.

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
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
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             interpolate_ocean_depth, is_floating, extrapolate_Gaussian_floodfill_mesh
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_BMB_model, type_remapping_mesh_mesh, &
                                             type_climate_snapshot_regional, type_ocean_snapshot_regional, &
                                             type_remapping_mesh_mesh, type_reference_geometry, type_restart_data, &
                                             type_grid
  USE forcing_module,                  ONLY: forcing, get_insolation_at_time_month_and_lat
  USE mesh_mapping_module,             ONLY: remap_field_dp_2D, remap_field_dp_3D, smooth_Gaussian_2D

  IMPLICIT NONE

CONTAINS

! ===== The main routines that should be called from the main ice model/program =====
! ==========================================================================

  SUBROUTINE run_BMB_model( mesh, ice, ocean, BMB, region_name, time, refgeo)
    ! Run the selected BMB model

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),   INTENT(INOUT) :: ocean
    TYPE(type_BMB_model),                 INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name
    REAL(dp),                             INTENT(IN)    :: time
    TYPE(type_reference_geometry),        INTENT(IN)    :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model'
    INTEGER                                            :: vi
    REAL(dp)                                           :: amplification_factor

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Re-initialise total and grounded BMB
    BMB%BMB(       mesh%vi1:mesh%vi2) = 0._dp
    BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp

    ! Re-initialise BMB over ice shelves
    IF (C%choice_BMB_shelf_model == 'inversion') THEN
      DO vi = mesh%vi1, mesh%vi2
        IF ( ice%mask_sheet_a( vi) == 1 .AND. (.NOT. is_floating( refgeo%Hi( vi), refgeo%Hb( vi), 0._dp)) ) THEN
           BMB%BMB_shelf( vi) = 0._dp
        END IF
      END DO
    ELSE
      BMB%BMB_shelf( mesh%vi1:mesh%vi2) = 0._dp
    END IF
    CALL sync

    ! Run the selected shelf BMB model
    IF     (C%choice_BMB_shelf_model == 'uniform') THEN
      BMB%BMB_shelf( mesh%vi1:mesh%vi2) = C%BMB_shelf_uniform
      CALL sync
    ELSEIF (C%choice_BMB_shelf_model == 'idealised') THEN
      CALL run_BMB_model_idealised(            mesh, ice, BMB, time)
    ELSEIF (C%choice_BMB_shelf_model == 'ANICE_legacy') THEN
      CALL run_BMB_model_ANICE_legacy(         mesh, ice, BMB, region_name, time)
    ELSEIF (C%choice_BMB_shelf_model == 'Favier2019_lin') THEN
      CALL run_BMB_model_Favier2019_linear(    mesh, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'Favier2019_quad') THEN
      CALL run_BMB_model_Favier2019_quadratic( mesh, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'Favier2019_Mplus') THEN
      CALL run_BMB_model_Favier2019_Mplus(     mesh, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'Lazeroms2018_plume') THEN
      CALL run_BMB_model_Lazeroms2018_plume(   mesh, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'PICO') THEN
      CALL run_BMB_model_PICO(                 mesh, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'PICOP') THEN
      CALL run_BMB_model_PICOP(                mesh, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'inversion') THEN
      CALL run_BMB_model_shelf_inversion(      mesh, ice, BMB, refgeo)
    ELSEIF (C%choice_BMB_shelf_model == 'Bernales202X') THEN
      CALL run_BMB_model_Bernales202X( mesh, ice, ocean, BMB, refgeo)
    ELSE
      CALL crash('unknown choice_BMB_shelf_model "' // TRIM(C%choice_BMB_shelf_model) // '"!')
    END IF

    ! Run the selected sheet BMB model
    IF     (C%choice_BMB_sheet_model == 'uniform') THEN
      BMB%BMB_sheet( mesh%vi1:mesh%vi2) = C%BMB_sheet_uniform
    ELSE
      CALL crash('unknown choice_BMB_sheet_model "' // TRIM(C%choice_BMB_sheet_model) // '"!')
    END IF

    ! Extrapolate melt field from the regular (FCMP) mask to the (extended) PMP mask
    IF (C%choice_BMB_subgrid == 'PMP') CALL extrapolate_melt_from_FCMP_to_PMP( mesh, ice, BMB)

    ! If desired, the BMB can be tuned locally by changing the amplification factor
    DO vi = mesh%vi1, mesh%vi2

      amplification_factor = 1._dp

      IF (C%choice_BMB_shelf_amplification == 'uniform') THEN
        amplification_factor = 1._dp
      ELSE IF (C%choice_BMB_shelf_amplification == 'basin') THEN
        IF (region_name == 'ANT') THEN
          amplification_factor = C%basin_BMB_amplification_factor_ANT( ice%basin_ID( vi))
        ELSE IF (region_name == 'GRL') THEN
          amplification_factor = C%basin_BMB_amplification_factor_GRL( ice%basin_ID( vi))
        ELSE
          amplification_factor = 1._dp
        END IF
      ELSE
        CALL crash('unknown choice_BMB_shelf_amplification "' // TRIM(C%choice_BMB_shelf_amplification) // '"!')
      END IF

      BMB%BMB_shelf( vi) = amplification_factor * BMB%BMB_shelf( vi)

    END DO
    CALL sync

    ! Add sheet and shelf melt rates together, applying the selected scheme for sub-grid shelf melt
    ! (see Leguy et al. 2021 for explanations of the three schemes)
    DO vi = mesh%vi1, mesh%vi2

      ! Add sub-sheet melt rates (no sub-grid scaling yet)
      BMB%BMB( vi) = 0._dp
      IF (ice%mask_sheet_a( vi) == 1) THEN
        BMB%BMB( vi) = BMB%BMB_sheet( vi)
      END IF

      ! Add sub-shelf melt rates
      IF (C%choice_BMB_shelf_model == 'inversion') THEN
        ! For the inversion of melt rates or temperatures, add inverted rates
        ! everywhere. BMB_shelf can be non-zero only at points where the model
        ! OR the reference data is shelf or ocean (this helps to account for
        ! 'unwanted' grounding line advance or retreat, and to get an idea of
        ! the melt rates needed at ocean points). Thus, no need for masks here.
        BMB%BMB( vi) = BMB%BMB( vi) + BMB%BMB_shelf( vi)
      ELSE
        ! Different sub-grid schemes for sub-shelf melt. BMB_shelf can be non-zero
        ! only at ice shelf points, so no need for masks here.
        IF     (C%choice_BMB_subgrid == 'FCMP') THEN
          IF (ice%mask_shelf_a( vi) == 1) BMB%BMB( vi) = BMB%BMB( vi) + BMB%BMB_shelf( vi)
        ELSEIF (C%choice_BMB_subgrid == 'PMP') THEN
          BMB%BMB( vi) = BMB%BMB( vi) + (1._dp - ice%f_grndx_a( vi)) * BMB%BMB_shelf( vi)
        ELSEIF (C%choice_BMB_subgrid == 'NMP') THEN
          IF (ice%f_grndx_a( vi) == 0._dp) BMB%BMB( vi) = BMB%BMB( vi) + BMB%BMB_shelf( vi)
        ELSE
          CALL crash('unknown choice_BMB_subgrid "' // TRIM(C%choice_BMB_subgrid) // '"!')
        END IF
      END IF

    END DO
    CALL sync

    ! Limit basal melt
    DO vi = mesh%vi1, mesh%vi2
      BMB%BMB( vi) = MAX( BMB%BMB( vi), C%BMB_min)
      BMB%BMB( vi) = MIN( BMB%BMB( vi), C%BMB_max)
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( BMB%BMB_sheet, 'BMB%BMB_sheet')
    CALL check_for_NaN_dp_1D( BMB%BMB_shelf, 'BMB%BMB_shelf')
    CALL check_for_NaN_dp_1D( BMB%BMB,       'BMB%BMB'      )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model

  SUBROUTINE initialise_BMB_model( mesh, ice, ocean, BMB, region_name, restart)
    ! Allocate memory for the data fields of the SMB model.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_restart_data),             INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '  Initialising BMB model: sheet = "', TRIM(C%choice_BMB_sheet_model), '", shelf = "', TRIM(C%choice_BMB_shelf_model), '"...'

    ! General
    CALL allocate_shared_dp_1D( mesh%nV, BMB%BMB      , BMB%wBMB      )
    CALL allocate_shared_dp_1D( mesh%nV, BMB%BMB_shelf, BMB%wBMB_shelf)
    CALL allocate_shared_dp_1D( mesh%nV, BMB%BMB_sheet, BMB%wBMB_sheet)

    ! Shelf
    IF     (C%choice_BMB_shelf_model == 'uniform' .OR. &
            C%choice_BMB_shelf_model == 'idealised') THEN
      ! Nothing else needs to be done
    ELSEIF (C%choice_BMB_shelf_model == 'ANICE_legacy') THEN
      CALL initialise_BMB_model_ANICE_legacy( mesh, BMB, region_name)
    ELSEIF (C%choice_BMB_shelf_model == 'Favier2019_lin' .OR. &
            C%choice_BMB_shelf_model == 'Favier2019_quad' .OR. &
            C%choice_BMB_shelf_model == 'Favier2019_Mplus') THEN
      CALL initialise_BMB_model_Favier2019( mesh, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'Lazeroms2018_plume') THEN
      CALL initialise_BMB_model_Lazeroms2018_plume( mesh, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'PICO') THEN
      CALL initialise_BMB_model_PICO(  mesh, ice, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'PICOP') THEN
      CALL initialise_BMB_model_PICOP( mesh, ice, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'inversion') THEN
      CALL initialise_BMB_model_inversion( mesh, BMB, restart)
    ELSEIF (C%choice_BMB_shelf_model == 'Bernales202X') THEN
      CALL initialise_BMB_model_Bernales202X( mesh, ice, ocean, BMB)
    ELSE
      CALL crash('unknown choice_BMB_shelf_model "' // TRIM(C%choice_BMB_shelf_model) // '"!')
    END IF

    ! Sheet
    IF     (C%choice_BMB_sheet_model == 'uniform') THEN
      ! Nothing else needs to be done
    ELSE
      CALL crash('unknown choice_BMB_sheet_model "' // TRIM(C%choice_BMB_sheet_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=16)

  END SUBROUTINE initialise_BMB_model

! ===== Idealised BMB schemes =====
! =================================

  SUBROUTINE run_BMB_model_idealised( mesh, ice, BMB, time)
    ! Idealised BMB schemes

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_idealised_BMB_shelf == 'MISMIP+') THEN
      ! The schematic basal melt used in the MISMIPplus experiments

      CALL run_BMB_model_MISMIPplus( mesh, ice, BMB, time)

    ELSE
      CALL crash('unknown choice_idealised_BMB_shelf "' // TRIM(C%choice_idealised_BMB_shelf) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_idealised

  SUBROUTINE run_BMB_model_MISMIPplus( mesh, ice, BMB, time)
    ! The schematic basal melt used in the MISMIPplus experiments

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_MISMIPplus'
    INTEGER                                            :: vi
    REAL(dp)                                           :: zd, cavity_thickness

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%MISMIPplus_scenario == 'ice0') THEN
      ! The reference scenario; no basal melt ever

      BMB%BMB(       mesh%vi1:mesh%vi2) = 0._dp
      BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
      BMB%BMB_shelf( mesh%vi1:mesh%vi2) = 0._dp
      CALL sync

    ELSEIF (C%MISMIPplus_scenario == 'ice1ra') THEN
      ! Increased melt for 100 yr, followed by zero melt for 100 yr

      IF (time < 0._dp) THEN
        ! t = -10,000 to t = 0: spin-up, no melt

        BMB%BMB(       mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_shelf( mesh%vi1:mesh%vi2) = 0._dp
        CALL sync

      ELSEIF (time < 100._dp) THEN
        ! t = 0 to t = 100: melt

        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
        DO vi = mesh%vi1, mesh%vi2

          zd = ice%Hs_a( vi) - ice%Hi_a( vi)
          cavity_thickness = MAX( 0._dp, zd - ice%Hb_a( vi))

          ! Cornford et al. (2020), Eq. 7
          BMB%BMB_shelf( vi) = -0.2_dp * TANH( cavity_thickness / 75._dp) * MAX( -100._dp - zd, 0._dp)

        END DO
        CALL sync

      ELSE ! IF (time < 0._dp) THEN
        ! After t = 100: no melt

        BMB%BMB(       mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_shelf( mesh%vi1:mesh%vi2) = 0._dp
        CALL sync

      END IF ! IF (time < 0._dp) THEN

    ELSEIF (C%MISMIPplus_scenario == 'ice1rr') THEN
      ! Increased melt forever

      IF (time < 0._dp) THEN
        ! t = -10,000 to t = 0: spin-up, no melt

        BMB%BMB(       mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_shelf( mesh%vi1:mesh%vi2) = 0._dp
        CALL sync

      ELSE ! IF (time < 0._dp) THEN
        ! t > 0: melt

        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
        DO vi = mesh%vi1, mesh%vi2

          zd = ice%Hs_a( vi) - ice%Hi_a( vi)
          cavity_thickness = MAX( 0._dp, zd - ice%Hb_a( vi))

          ! Cornford et al. (2020), Eq. 7
          BMB%BMB_shelf( vi) = -0.2_dp * TANH( cavity_thickness / 75._dp) * MAX( -100._dp - zd, 0._dp)

        END DO
        CALL sync

      END IF ! IF (time < 0._dp) THEN

    ELSEIF (C%MISMIPplus_scenario == 'ice2ra') THEN
      ! Increased "calving" for 100 yr, followed by zero "calving" for 100 yr

      IF (time < 0._dp) THEN
        ! t = -10,000 to t = 0: spin-up, no "calving"

        BMB%BMB(       mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_shelf( mesh%vi1:mesh%vi2) = 0._dp
        CALL sync

      ELSEIF (time < 100._dp) THEN
        ! t = 0 to t = 100: "calving"

        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
        DO vi = mesh%vi1, mesh%vi2
          IF (mesh%V( vi,1) > 80000._dp) THEN ! Actually the border is at x = 480 km, but our coordinate system is shifted...
            BMB%BMB_shelf( vi) = -100._dp
          ELSE
            BMB%BMB_shelf( vi) = 0._dp
          END IF
        END DO
        CALL sync

      ELSE ! IF (time < 0._dp) THEN
        ! After t = 100: no "calving"

        BMB%BMB(       mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_shelf( mesh%vi1:mesh%vi2) = 0._dp
        CALL sync

      END IF ! IF (time < 0._dp) THEN

    ELSEIF (C%MISMIPplus_scenario == 'ice2rr') THEN
      ! Increased "calving" forever

      IF (time < 0._dp) THEN
        ! t = -10,000 to t = 0: spin-up, no "calving"

        BMB%BMB(       mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
        BMB%BMB_shelf( mesh%vi1:mesh%vi2) = 0._dp
        CALL sync

      ELSE
        ! t > 0: "calving"

        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
        DO vi = mesh%vi1, mesh%vi2
          IF (mesh%V( vi,1) > 80000._dp) THEN ! Actually the border is at x = 480 km, but our coordinate system is shifted...
            BMB%BMB_shelf( vi) = -100._dp
          ELSE
            BMB%BMB_shelf( vi) = 0._dp
          END IF
        END DO
        CALL sync

      END IF ! IF (time < 0._dp) THEN

    ELSE
      CALL crash('unknown MISMIPplus_scenario "' // TRIM(C%MISMIPplus_scenario) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_MISMIPplus

! ===== The ANICE_legacy sub-shelf melt model =====
! =================================================

  SUBROUTINE run_BMB_model_ANICE_legacy( mesh, ice, BMB, region_name, time)
    ! Calculate sub-shelf melt with the ANICE_legacy model, which is based on the glacial-interglacial
    ! parameterisation by Pollard & DeConto (2012), the distance-to-open-ocean + subtended-angle parameterisation
    ! by Pollard & DeConto (2012), and the linear temperature-melt relation by Martin et al. (2011).

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_ANICE_legacy'
    INTEGER                                            :: vi
    REAL(dp)                                           :: T_ocean_mean, Q_TOA_jun_65N, Q_TOA_jan_80S
    REAL(dp)                                           :: BMB_shelf                             ! Sub-shelf melt rate for non-exposed shelf  [m/year]
    REAL(dp)                                           :: BMB_shelf_exposed                     ! Sub-shelf melt rate for exposed shelf      [m/year]
    REAL(dp)                                           :: BMB_deepocean                         ! Sub-shelf melt rate for deep-ocean areas   [m/year]
    REAL(dp)                                           :: w_ins, w_PD, w_warm, w_cold, w_deep, w_expo
    REAL(dp)                                           :: T_freeze                              ! Freezing temperature at the base of the shelf (Celcius)
    REAL(dp)                                           :: water_depth
    REAL(dp), PARAMETER                                :: cp0     = 3974._dp                    ! specific heat capacity of the ocean mixed layer (J kg-1 K-1)
    REAL(dp), PARAMETER                                :: gamma_T = 1.0E-04_dp                  ! Thermal exchange velocity (m s-1)

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    BMB%sub_angle( mesh%vi1:mesh%vi2) = 360._dp
    BMB%dist_open( mesh%vi1:mesh%vi2) = 0._dp
    w_ins                             = 0._dp
    w_PD                              = 0._dp
    w_warm                            = 0._dp
    w_cold                            = 0._dp
    w_deep                            = 0._dp
    w_expo                            = 0._dp
    BMB_shelf                         = 0._dp
    BMB_shelf_exposed                 = 0._dp
    BMB_deepocean                     = 0._dp
    CALL sync

    ! Find the "subtended angle" and distance-to-open-ocean of all shelf pixels
    DO vi = mesh%vi1, mesh%vi2
      CALL calculate_sub_angle_dist_open( mesh, ice, vi, BMB%sub_angle( vi), BMB%dist_open( vi))
    END DO
    CALL sync

    ! Find the weight from insolation
    IF (region_name == 'NAM' .OR. region_name == 'EAS' .OR. region_name == 'GRL') THEN
      CALL get_insolation_at_time_month_and_lat( time, 6, 65._dp, Q_TOA_jun_65N)
      w_ins = MAX(0._dp, (Q_TOA_jun_65N - 462.29_dp) / 40._dp)
    ELSEIF (region_name == 'ANT') THEN
      CALL get_insolation_at_time_month_and_lat( time, 1, -80._dp, Q_TOA_jan_80S)
      w_ins = MAX(0._dp, (Q_TOA_jan_80S - 532.19_dp) / 40._dp)
    END IF

    ! Initialise to prevent compiler warnings
    T_ocean_mean = 0._dp

    ! Determine mean ocean temperature and basal melt rates for deep ocean and exposed shelves
    IF (C%choice_forcing_method == 'none') THEN

      ! In this case, no CO2/d18O forcing is used; just assume PD weights
      w_warm = 0._dp
      w_cold = 0._dp
      w_PD   = 1._dp

    ELSE IF (C%choice_forcing_method == 'CO2_direct') THEN

      ! Use the prescribed CO2 record as a glacial index
      IF (forcing%CO2_obs > 280._dp) THEN
        ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 400 ppmv
        w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (400._dp - forcing%CO2_obs) / (400._dp - 280._dp) )) - w_ins
        w_warm = 1._dp - w_PD
        w_cold = 0._dp
      ELSE
        ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
        w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (forcing%CO2_obs - 190._dp) / (280._dp - 190._dp) )) + w_ins
        w_cold = 1._dp - w_PD
        w_warm = 0._dp
      END IF

    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN

      ! Use modelled CO2 as a glacial index
      IF (forcing%CO2_mod > 280._dp) THEN
        ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 400 ppmv
        w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (400._dp - forcing%CO2_mod) / (400._dp - 280._dp) )) - w_ins
        w_warm = 1._dp - w_PD
        w_cold = 0._dp
      ELSE
        ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
        w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (forcing%CO2_mod - 190._dp) / (280._dp - 190._dp) )) + w_ins
        w_cold = 1._dp - w_PD
        w_warm = 0._dp
      END IF

    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN

      ! Use modelled global mean annual temperature change as a glacial index
      IF (forcing%dT_glob_inverse > 0._dp) THEN
        ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 405 ppmv
        w_warm = MIN(1.25_dp, MAX(-0.25_dp, forcing%dT_glob_inverse / 12._dp )) - w_ins
        w_PD   = 1._dp - w_warm
        w_cold = 0._dp
      ELSE
        ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
        w_cold = MIN(1.25_dp, MAX(-0.25_dp, -forcing%dT_glob_inverse / 12._dp )) + w_ins
        w_PD   = 1._dp - w_cold
        w_warm = 0._dp
      END IF

    ELSEIF (C%choice_forcing_method == 'climate_direct' .OR. C%choice_forcing_method == 'SMB_direct') THEN

      ! In this case, no CO2/d18O forcing is used; just assume PD weights
      w_warm = 0._dp
      w_cold = 0._dp
      w_PD   = 1._dp

    ELSE ! IF (C%choice_forcing_method == 'CO2_direct') THEN
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF ! IF (C%choice_forcing_method == 'CO2_direct') THEN

    T_ocean_mean      = w_PD * BMB%T_ocean_mean_PD      + w_warm * BMB%T_ocean_mean_warm      + w_cold * BMB%T_ocean_mean_cold
    BMB_deepocean     = w_PD * BMB%BMB_deepocean_PD     + w_warm * BMB%BMB_deepocean_warm     + w_cold * BMB%BMB_deepocean_cold
    BMB_shelf_exposed = w_PD * BMB%BMB_shelf_exposed_PD + w_warm * BMB%BMB_shelf_exposed_warm + w_cold * BMB%BMB_shelf_exposed_cold
    CALL sync

    ! Use the (interpolated, spatially uniform) ocean temperature and the subtended angle + distance-to-open-ocean
    ! to calculate sub-shelf melt rates using the parametrisation from Martin et al., 2011
    DO vi = mesh%vi1, mesh%vi2

      IF (ice%mask_shelf_a( vi) == 1) THEN
        ! Sub-shelf melt

        ! Freezing temperature at the bottom of the ice shelves, scaling with depth below water level
        T_freeze = 0.0939_dp - 0.057_dp * 35._dp - 7.64E-04_dp * ice%Hi_a( vi) * ice_density / seawater_density

        ! Sub-shelf melt rate for non-exposed shelves (Martin, TC, 2011) - melt values, when T_ocean > T_freeze.
        BMB_shelf   = seawater_density * cp0 * sec_per_year * gamma_T * BMB%subshelf_melt_factor * &
                   (T_ocean_mean - T_freeze) / (L_fusion * ice_density)

      ELSE
        BMB_shelf = 0._dp
      END IF

      IF (ice%mask_shelf_a( vi) == 1 .OR. ice%mask_ocean_a( vi) == 1) THEN

        water_depth = ice%SL_a( vi) - ice%Hb_a( vi)
        w_deep = MAX(0._dp, MIN(1._dp, (water_depth - BMB%deep_ocean_threshold_depth) / 200._dp))
        w_expo = MAX(0._dp, MIN(1._dp, (BMB%sub_angle( vi) - 80._dp)/30._dp)) * EXP(-BMB%dist_open( vi) / 100000._dp)

        ! As the subtended angle and distance to open ocean is not straightforward on a triangular mesh,
        ! for now disregard their effects by setting the weight for exposed shelves to 0. Come back to
        ! this later if necessary.
        w_expo = 0._dp

        BMB%BMB_shelf( vi) = (w_deep * BMB_deepocean) + (1._dp - w_deep) * (w_expo * BMB_shelf_exposed + (1._dp - w_expo) * BMB_shelf)

      ELSE
        BMB%BMB_shelf( vi) = 0._dp
      END IF

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( BMB%sub_angle, 'BMB%sub_angle')
    CALL check_for_NaN_dp_1D( BMB%dist_open, 'BMB%dist_open')
    CALL check_for_NaN_dp_1D( BMB%BMB_shelf, 'BMB%BMB_shelf')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_ANICE_legacy

  SUBROUTINE initialise_BMB_model_ANICE_legacy( mesh, BMB, region_name)
    ! Allocate memory for the data fields of the ANICE_legacy shelf BMB model.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_ANICE_legacy'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Variables
    CALL allocate_shared_dp_1D( mesh%nV, BMB%sub_angle, BMB%wsub_angle)
    CALL allocate_shared_dp_1D( mesh%nV, BMB%dist_open, BMB%wdist_open)

    ! Tuning parameters
    CALL allocate_shared_dp_0D( BMB%T_ocean_mean_PD,            BMB%wT_ocean_mean_PD           )
    CALL allocate_shared_dp_0D( BMB%T_ocean_mean_cold,          BMB%wT_ocean_mean_cold         )
    CALL allocate_shared_dp_0D( BMB%T_ocean_mean_warm,          BMB%wT_ocean_mean_warm         )
    CALL allocate_shared_dp_0D( BMB%BMB_deepocean_PD,           BMB%wBMB_deepocean_PD          )
    CALL allocate_shared_dp_0D( BMB%BMB_deepocean_cold,         BMB%wBMB_deepocean_cold        )
    CALL allocate_shared_dp_0D( BMB%BMB_deepocean_warm,         BMB%wBMB_deepocean_warm        )
    CALL allocate_shared_dp_0D( BMB%BMB_shelf_exposed_PD,       BMB%wBMB_shelf_exposed_PD      )
    CALL allocate_shared_dp_0D( BMB%BMB_shelf_exposed_cold,     BMB%wBMB_shelf_exposed_cold    )
    CALL allocate_shared_dp_0D( BMB%BMB_shelf_exposed_warm,     BMB%wBMB_shelf_exposed_warm    )
    CALL allocate_shared_dp_0D( BMB%subshelf_melt_factor,       BMB%wsubshelf_melt_factor      )
    CALL allocate_shared_dp_0D( BMB%deep_ocean_threshold_depth, BMB%wdeep_ocean_threshold_depth)

    IF (region_name == 'NAM') THEN
      BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_NAM
      BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_NAM
      BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_NAM
      BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_NAM
      BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_NAM
      BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_NAM
      BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_NAM
      BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_NAM
      BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_NAM
      BMB%subshelf_melt_factor       = C%subshelf_melt_factor_NAM
      BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_NAM
    ELSEIF (region_name == 'EAS') THEN
      BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_EAS
      BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_EAS
      BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_EAS
      BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_EAS
      BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_EAS
      BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_EAS
      BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_EAS
      BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_EAS
      BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_EAS
      BMB%subshelf_melt_factor       = C%subshelf_melt_factor_EAS
      BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_EAS
    ELSEIF (region_name == 'GRL') THEN
      BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_GRL
      BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_GRL
      BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_GRL
      BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_GRL
      BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_GRL
      BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_GRL
      BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_GRL
      BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_GRL
      BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_GRL
      BMB%subshelf_melt_factor       = C%subshelf_melt_factor_GRL
      BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_GRL
    ELSEIF (region_name == 'ANT') THEN
      BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_ANT
      BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_ANT
      BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_ANT
      BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_ANT
      BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_ANT
      BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_ANT
      BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_ANT
      BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_ANT
      BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_ANT
      BMB%subshelf_melt_factor       = C%subshelf_melt_factor_ANT
      BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_ANT
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=13)

  END SUBROUTINE initialise_BMB_model_ANICE_legacy

! ===== The Favier et al. (2019) sub-shelf melt parameterisations =====
! =====================================================================

  SUBROUTINE run_BMB_model_Favier2019_linear( mesh, ice, ocean, BMB)
    ! Calculate sub-shelf melt with Favier et al. (2019) linear parameterisation

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_Favier2019_linear'
    INTEGER                                            :: vi
    REAL(dp)                                           :: dT

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate ocean temperature and freezing point at the base of the shelf
    CALL calc_ocean_temperature_at_shelf_base(    mesh, ice, ocean, BMB)
    CALL calc_ocean_freezing_point_at_shelf_base( mesh, ice, ocean, BMB)

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise
      BMB%BMB_shelf( vi) = 0._dp

      IF (ice%mask_shelf_a( vi) == 1) THEN

        ! Temperature forcing
        dT = BMB%T_ocean_base( vi) - BMB%T_ocean_freeze_base( vi)

        ! Favier et al. (2019), Eq. 2
        BMB%BMB_shelf( vi) = -sec_per_year * C%BMB_Favier2019_lin_GammaT * seawater_density * cp_ocean * dT / (ice_density * L_fusion)

      END IF ! IF (ice%mask_shelf_a( vi) == 1) THEN

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_Favier2019_linear

  SUBROUTINE run_BMB_model_Favier2019_quadratic( mesh, ice, ocean, BMB)
    ! Calculate sub-shelf melt with Favier et al. (2019) quadratic parameterisation

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_Favier2019_quadratic'
    INTEGER                                            :: vi
    REAL(dp)                                           :: dT

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate ocean temperature and freezing point at the base of the shelf
    CALL calc_ocean_temperature_at_shelf_base(    mesh, ice, ocean, BMB)
    CALL calc_ocean_freezing_point_at_shelf_base( mesh, ice, ocean, BMB)

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise
      BMB%BMB_shelf( vi) = 0._dp

      IF (ice%mask_shelf_a( vi) == 1) THEN

        ! Temperature forcing
        dT = BMB%T_ocean_base( vi) - BMB%T_ocean_freeze_base( vi)

        ! Favier et al. (2019), Eq. 4
        ! BMB%BMB_shelf( vi) = -sec_per_year * C%BMB_Favier2019_quad_GammaT * (seawater_density * cp_ocean * dT / (ice_density * L_fusion))**2._dp
        ! Altered to allow for negative basal melt (i.e. refreezing) when dT < 0
        BMB%BMB_shelf( vi) = -sec_per_year * C%BMB_Favier2019_quad_GammaT * SIGN(dT,1._dp) * (seawater_density * cp_ocean * dT / (ice_density * L_fusion))**2._dp

      END IF ! IF (ice%mask_shelf_a( vi) == 1) THEN

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_Favier2019_quadratic

  SUBROUTINE run_BMB_model_Favier2019_Mplus( mesh, ice, ocean, BMB)
    ! Calculate sub-shelf melt with Favier et al. (2019) M+ parameterisation

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_Favier2019_Mplus'
    INTEGER                                            :: vi,n_shelf,basin_i
    REAL(dp)                                           :: dT,dT_av

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate ocean temperature and freezing point at the base of the shelf
    CALL calc_ocean_temperature_at_shelf_base(    mesh, ice, ocean, BMB)
    CALL calc_ocean_freezing_point_at_shelf_base( mesh, ice, ocean, BMB)

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise
      BMB%BMB_shelf( vi) = 0._dp

      IF (ice%mask_shelf_a( vi) == 1) THEN

        ! Temperature forcing
        dT = BMB%T_ocean_base( vi) - BMB%T_ocean_freeze_base( vi)

        ! Favier et al. (2019), Eq. 5 (part 1, without the basin-averaged term, that comes later)
        BMB%BMB_shelf( vi) = -sec_per_year * C%BMB_Favier2019_Mplus_GammaT * (seawater_density * cp_ocean / (ice_density * L_fusion))**2._dp * dT

      END IF ! IF (ice%mask_shelf_a( vi) == 1) THEN

    END DO
    CALL sync

    ! Calculate and apply basin-averaged temperature forcing
    DO basin_i = 1, ice%nbasins

      dT_av   = 0._dp
      n_shelf = 0

      DO vi = mesh%vi1, mesh%vi2
        IF (ice%basin_ID( vi) == basin_i .AND. ice%mask_shelf_a( vi) == 1) THEN
          dT = BMB%T_ocean_base( vi) - BMB%T_ocean_freeze_base( vi)
          dT_av   = dT_av   + dT
          n_shelf = n_shelf + 1
        END IF
      END DO

      CALL MPI_ALLREDUCE( MPI_IN_PLACE, dT_av,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_shelf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! Safety
      IF (n_shelf == 0) THEN
        dT_av = 0._dp
      ELSE
        dT_av = dT_av / n_shelf
      END IF

      ! Add last term to the equation
      DO vi = mesh%vi1, mesh%vi2
        IF (ice%basin_ID( vi) == basin_i .AND. ice%mask_shelf_a( vi) == 1) THEN
          BMB%BMB_shelf( vi) = BMB%BMB_shelf( vi) * dT_av
        END IF
      END DO
      CALL sync

    END DO ! DO basin_i = 1, ice%nbasins

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_Favier2019_Mplus

  SUBROUTINE calc_ocean_freezing_point_at_shelf_base( mesh, ice, ocean, BMB)
    ! Calculate the ocean freezing point at the base of the shelf, needed to calculate
    ! basal melt in the different parameterisations from Favier et al. (2019)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_ocean_freezing_point_at_shelf_base'
    INTEGER                                            :: vi
    REAL(dp)                                           :: depth
    REAL(dp)                                           :: S0                   ! Practical salinity [PSU]
    REAL(dp), PARAMETER                                :: lambda1 = -0.0575_dp ! Liquidus slope                [degC PSU^-1] (Favier et al. (2019), Table 2)
    REAL(dp), PARAMETER                                :: lambda2 = 0.0832_dp  ! Liquidus intercept            [degC]        (Favier et al. (2019), Table 2)
    REAL(dp), PARAMETER                                :: lambda3 = 7.59E-4_dp ! Liquidus pressure coefficient [degC m^-1]   (Favier et al. (2019), Table 2)

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise at zero
      BMB%T_ocean_freeze_base( vi) = 0._dp

      IF (ice%mask_shelf_a( vi) == 1) THEN

        ! Calculate depth
        depth = MAX( 0.1_dp, ice%Hi_a( vi) * ice_density / seawater_density)

        ! Find salinity at this depth
        CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S_ocean_corr_ext( vi,:), depth, S0)

        ! Calculate ocean freezing temperature (Favier et al. (2019), Eq. 3) in degrees Celsius
        BMB%T_ocean_freeze_base( vi) = lambda1 * S0 + lambda2 - lambda3 * depth

      END IF ! IF (ice%mask_shelf_a( vi) == 1) THEN

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_ocean_freezing_point_at_shelf_base

  SUBROUTINE initialise_BMB_model_Favier2019( mesh, BMB)
    ! Allocate memory for the data fields of the Favier et al. (2019) shelf BMB parameterisations.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_Favier2019'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Variables
    CALL allocate_shared_dp_1D( mesh%nV, BMB%T_ocean_base,        BMB%wT_ocean_base       )
    CALL allocate_shared_dp_1D( mesh%nV, BMB%T_ocean_freeze_base, BMB%wT_ocean_freeze_base)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=2)

  END SUBROUTINE initialise_BMB_model_Favier2019

! ===== The Lazeroms et al. (2018) quasi-2-D plume parameterisation =====
! =======================================================================

  SUBROUTINE run_BMB_model_Lazeroms2018_plume( mesh, ice, ocean, BMB)
    ! Calculate basal melt using the quasi-2-D plume parameterisation by Lazeroms et al. (2018), following the equations in Appendix A

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_Lazeroms2018_plume'
    INTEGER                                            :: vi
    REAL(dp)                                           :: depth                               ! Ice-shelf base depth ("draft") [m]
    REAL(dp)                                           :: Ta                                  ! Ambient temperature at the ice-shelf base [degC]
    REAL(dp)                                           :: Sa                                  ! Ambient salinity    at the ice-shelf base [PSU]
    REAL(dp)                                           :: Tf_GL                               ! Freezing temperature at effective grounding-line plume source
    REAL(dp)                                           :: alpha                               ! Effective local angle ( = atan(slope))
    REAL(dp)                                           :: sinalpha                            ! sin( alpha) (appears so often that it's better to give it its own variable)
    REAL(dp)                                           :: GammaTS                             ! Effective heat exchange coefficient
    REAL(dp)                                           :: g_alpha                             ! Geometry term
    REAL(dp)                                           :: g_alpha_term1, g_alpha_term2, g_alpha_term3, sqrtCd_GammaTS
    REAL(dp)                                           :: M                                   ! Empirically derived melt-rate scale
    REAL(dp)                                           :: l_geo                               ! Geometry- and temperature-dependent length scale
    REAL(dp)                                           :: Xhat                                ! Dimensionless scaled distance along plume path
    REAL(dp)                                           :: Mhat                                ! Dimensionless melt curve

    ! Constant parameters (Lazeroms et al. (2018), Table 1)
    REAL(dp), PARAMETER                                :: E0              =  3.6E-2_dp        ! Entrainment coefficient             [unitless]
    REAL(dp), PARAMETER                                :: Cd              =  2.5E-3_dp        ! Drag coefficient                    [unitless]
    REAL(dp), PARAMETER                                :: lambda1         = -0.0573_dp        ! Freezing point salinity coefficient [K]
    REAL(dp), PARAMETER                                :: lambda2         =  0.0832_dp        ! Freezing point offset               [K]
    REAL(dp), PARAMETER                                :: lambda3         =  7.61E-4_dp       ! Freezing point depth coefficient    [K m^-1]
    REAL(dp), PARAMETER                                :: M0              =  10._dp           ! Melt-rate parameter                 [m yr^1 degC^-2]
    REAL(dp), PARAMETER                                :: sqrtCd_GammaTS0 =  6.0E-4_dp        ! Heat exchange parameter             [unitless]
    REAL(dp), PARAMETER                                :: gamma1          =  0.545_dp         ! Heat exchange parameter             [unitless]
    REAL(dp), PARAMETER                                :: gamma2          =  3.5E-5_dp        ! Heat exchange parameter             [m^-1]
    REAL(dp), PARAMETER                                :: x0              =  0.56_dp          ! Empirically derived dimensionless scaling factor

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not yet implemented for a triangular mesh!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ! ! Calculate ocean temperature at the base of the shelf
    ! CALL calc_ocean_temperature_at_shelf_base( mesh, ice, ocean, BMB)

    ! DO vi = mesh%vi1, mesh%vi2

    !   ! Initialise
    !   BMB%eff_plume_source_depth( vi) = 0._dp
    !   BMB%eff_basal_slope(        vi) = 0._dp
    !   BMB%BMB_shelf(              vi) = 0._dp

    !   IF (ice%mask_shelf_a( vi) == 1) THEN

    !     ! Find ambient temperature and salinity at the ice-shelf base
    !     IF (C%choice_BMB_shelf_model == 'Lazeroms2018_plume') THEN
    !       ! Use the extrapolated ocean temperature+salinity fields

    !       depth = MAX( 0.1_dp, ice%Hi_a( vi) - ice%Hs_a( vi))   ! Depth is positive when below the sea surface!
    !       CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( vi,:), depth, Ta)
    !       CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S_ocean_corr_ext( vi,:), depth, Sa)

    !     ELSEIF (C%choice_BMB_shelf_model == 'PICOP') THEN
    !       ! Use the results from the PICO ocean box model

    !       Ta = BMB%PICO_T( vi)
    !       Sa = BMB%PICO_S( vi)

    !     ELSE
    !       IF (par%master) WRITE(0,*) 'run_BMB_model_Lazeroms2018_plume - ERROR: finding ambient temperature+salinity only defined for choice_BMB_shelf_model = "Lazeroms2018_plume" or "PICOP"!'
    !       CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    !     END IF

    !     ! Determine the effective plume path to find the effective plume source depth and basal slope for this shelf grid cell
    !     CALL find_effective_plume_path( mesh, ice, BMB, vi, BMB%eff_plume_source_depth( vi), BMB%eff_basal_slope( vi))
    !     alpha = ATAN( BMB%eff_basal_slope( vi))
    !     sinalpha = SIN( alpha)

    !     ! Calculate freezing temperature at effective grounding-line plume source (Lazeroms et al., 2018, Eq. A7)
    !     Tf_GL = lambda1 * Sa + lambda2 + lambda3 * BMB%eff_plume_source_depth( vi)

    !     ! Calculate the effective heat exchange coefficient (Lazeroms et al., 2018, Eq. A8)
    !     GammaTS = C%BMB_Lazeroms2018_GammaT * (gamma1 + gamma2 * ((Ta - Tf_GL) / lambda3) * ((E0 * sinalpha) / (sqrtCD_GammaTS0 + E0 * sinalpha)))
    !     sqrtCd_GammaTS = SQRT( Cd) * GammaTS

    !     ! Calculate the geometrical factor in the melt-rate expression (Lazeroms et al., 2018, Eq. A3)
    !     g_alpha_term1 = SQRT(      sinalpha  / (Cd             + E0 * sinalpha))
    !     g_alpha_term2 = SQRT( sqrtCd_GammaTS / (sqrtCd_GammaTS + E0 * sinalpha))
    !     g_alpha_term3 =     ( E0 * sinalpha  / (sqrtCd_GammaTS + E0 * sinalpha))
    !     g_alpha = g_alpha_term1 * g_alpha_term2 * g_alpha_term3

    !     ! Calculate the empirically derived melt-rate scale (Lazeroms et al., 2018, Eq. A9)
    !     M = M0 * g_alpha * (Ta - Tf_GL)**2

    !     ! Calculate the geometry- and temperature-dependent length scale (Lazeroms et al., 2018, Eq. A10)
    !     l_geo = ((Ta - Tf_GL) / lambda3) * (x0 * sqrtCd_GammaTS + E0 * sinalpha) / (x0 * (sqrtCd_GammaTS + E0 * sinalpha))

    !     ! Calculate the dimensionless scaled distance along plume path (Lazeroms et al., 2018, Eq. A11),
    !     !   and evaluate the dimensionless melt curve at that point
    !     Xhat = MIN(1._dp, MAX(0._dp, (depth - BMB%eff_plume_source_depth( vi)) / l_geo ))
    !     Mhat = Lazeroms2018_dimensionless_melt_curve( Xhat)

    !     ! Finally, calculate the basal melt rate (Lazeroms et al., 2018, Eq. A12)
    !     BMB%BMB_shelf( vi) = -M * Mhat

    !   END IF ! IF (ice%mask_shelf_a( vi) == 1) THEN

    ! END DO
    ! CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_Lazeroms2018_plume

  SUBROUTINE find_effective_plume_path( mesh, ice, BMB, vi, eff_plume_source_depth, eff_basal_slope)
    ! Find the effective plume source depth and basal slope for shelf grid cell [i,j]

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp),                            INTENT(OUT)   :: eff_plume_source_depth
    REAL(dp),                            INTENT(OUT)   :: eff_basal_slope

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_effective_plume_path'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! IF     (C%BMB_Lazeroms2018_find_GL_scheme == 'GL_average') THEN
    !   CALL find_effective_plume_path_GL_average(     mesh, ice, BMB, vi, eff_plume_source_depth, eff_basal_slope)
    ! ELSEIF (C%BMB_Lazeroms2018_find_GL_scheme == 'along_ice_flow') THEN
    !   CALL find_effective_plume_path_along_ice_flow( mesh, ice,      vi, eff_plume_source_depth, eff_basal_slope)
    ! ELSE
    !   IF (par%master) WRITE(0,*) '  ERROR: BMB_Lazeroms2018_find_GL_scheme "', TRIM(C%BMB_Lazeroms2018_find_GL_scheme), '" not implemented in find_effective_plume_path!'
    !   CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ! END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE find_effective_plume_path

  SUBROUTINE find_effective_plume_path_GL_average( mesh, ice, BMB, vi, eff_plume_source_depth, eff_basal_slope)
    ! Find the effective plume source depth and basal slope for shelf grid cell [i,j], following
    ! the approach outlined in Lazeroms et al. (2018), Sect. 2.3
    !
    ! Based loosely on code from IMAU-ICE v1.0 by Heiko Goelzer (2019?)
    !
    ! The straight line between the shelf point and a grounding line point indicates the 1-D path that a meltwater plume from the grounding line may
    ! have taken to reach the upper point. From this path, we determine the local slope beneath the ice shelf needed for the UPP melt parametrization.
    ! The average grounding line depth and average slope are a measure of the combined effect of plumes from multiple directions.
    !
    ! ***NOTE*** Only grounding-line points lying DEEPER than the shelf point
    ! can be a valid physical source of the meltwater plume and, consequently,
    ! valid input for the UPP parametrization. Points lying above depth(i,j)
    ! are ignored. Furthermore, we only consider directions in which the LOCAL
    ! slope at (i,j) is positive.

    ! Based on distance_open_ocean by Bas de Boer:
    ! searches for grounding line points in 16 directions:
    !
    !             12 13 14
    !            11 ---- 15
    !           10 --  -- 16
    !          09 -- ij -- 01
    !           08 --  -- 02
    !            07 ---- 03
    !             06 05 04

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp),                            INTENT(OUT)   :: eff_plume_source_depth
    REAL(dp),                            INTENT(OUT)   :: eff_basal_slope

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_effective_plume_path_GL_average'
    INTEGER                                            :: n
    INTEGER                                            :: dpi,dpj,ip1,jp1,ip2,jp2
    REAL(dp)                                           :: basal_slope
    REAL(dp)                                           :: zb_shelf, dist, zb_dp
    LOGICAL                                            :: reached_end, found_source
    REAL(dp)                                           :: TAF1,TAF2,lambda_GL
    REAL(dp)                                           :: Hb1,Hb2,plume_source_depth
    INTEGER                                            :: n_valid_plumes
    REAL(dp)                                           :: sum_plume_source_depths, sum_basal_slopes

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Initialise
    ! n_valid_plumes          = 0
    ! sum_plume_source_depths = 0._dp
    ! sum_basal_slopes        = 0._dp

    ! ! Calculate the shelf base depth (draft)
    ! zb_shelf = ice%Hs_a( vi) - ice%Hi_a( vi)

    ! ! Investigate all 16 search directions
    ! DO n = 1, 16

    !   ! The search direction vector
    !   dpi = BMB%search_directions( n,1)
    !   dpj = BMB%search_directions( n,2)

    !   ! Initialise the search pointer at the shelf grid cell
    !   ip1 = i
    !   jp1 = j
    !   ip2 = i + dpi
    !   jp2 = j + dpj

    !   ! If the search direction already points out of the model domain, don't bother
    !   IF (ip2 < 1 .OR. ip2 > grid%nx .OR. jp2 < 1 .OR. jp2 > grid%ny) CYCLE

    !   ! Calculate the basal slope in this direction (Lazeroms et al. (2018), Eq. 12)
    !   zb_dp = ice%Hs_a( jp2,ip2) - ice%Hi_a( jp2,ip2)
    !   dist  = SQRT( REAL(dpi,dp)**2 + REAL(dpj,dp)**2) * grid%dx
    !   basal_slope = (zb_shelf - zb_dp) / dist

    !   ! If this slope is negative, this plume is not valid
    !   IF (basal_slope < 0._dp) CYCLE

    !   ! Search in this direction
    !   reached_end  = .FALSE.
    !   found_source = .FALSE.
    !   DO WHILE (.NOT. reached_end)

    !     ! If the pointer exits the model domain, stop the search
    !     IF (ip2 < 1 .OR. ip2 > grid%nx .OR. jp2 < 1 .OR. jp2 > grid%ny) THEN
    !       reached_end  = .TRUE.
    !       found_source = .FALSE.
    !       EXIT
    !     END IF

    !     ! If the pointer encounters open ocean, stop the search
    !     IF (ice%mask_ocean_a( jp2,ip2) == 1 .AND. ice%mask_shelf_a( jp2,ip2) == 0) THEN
    !       reached_end  = .TRUE.
    !       found_source = .FALSE.
    !       EXIT
    !     END IF

    !     ! If the pointer encounters grounded ice, there must be a grounding line here
    !     IF (ice%mask_sheet_a( jp2,ip2) == 1) THEN

    !       reached_end  = .TRUE.
    !       found_source = .TRUE.

    !       ! Interpolate the thickness above flotation to find the sub-grid grounding-line depth ( = plume source depth)
    !       TAF1 = ice%TAF_a( jp1,ip1)
    !       TAF2 = ice%TAF_a( jp2,ip2)
    !       lambda_GL = TAF1 / (TAF1 - TAF2)

    !       Hb1  = ice%Hb_a( jp1,ip1)
    !       Hb2  = ice%Hb_a( jp2,ip2)
    !       plume_source_depth = (1._dp - lambda_GL) * Hb1 + lambda_GL * Hb2

    !       ! If this grounding line is less deep than the shelf, don't use it
    !       IF (plume_source_depth > zb_shelf) THEN
    !         found_source = .FALSE.
    !       END IF

    !     END IF ! IF (ice%mask_sheet_a( jp2,ip2) == 1) THEN

    !     ! If none of these exceptions were triggered, advance the pointer in the search direction
    !     ip1 = ip2
    !     jp1 = jp2
    !     ip2 = ip1 + dpi
    !     jp2 = jp1 + dpj

    !   END DO ! DO WHILE (.NOT. reached_end)

    !   ! If a valid plume source was found, add it to the sum for averaging
    !   IF (found_source) THEN
    !     n_valid_plumes          = n_valid_plumes          + 1
    !     sum_plume_source_depths = sum_plume_source_depths + plume_source_depth
    !     sum_basal_slopes        = sum_basal_slopes        + basal_slope
    !   END IF

    ! END DO ! DO n = 1, 16

    ! ! Define the effective plume source depth and basal slope as the average
    ! ! of those values for all valid plume paths

    ! IF (n_valid_plumes > 0) THEN
    !   eff_plume_source_depth = sum_plume_source_depths / REAL(n_valid_plumes,dp)
    !   eff_basal_slope        = sum_basal_slopes        / REAL(n_valid_plumes,dp)
    ! ELSE
    !   ! Exception for when no valid plume sources were found
    !   eff_plume_source_depth = zb_shelf
    !   eff_basal_slope        = 1E-10_dp   ! Because the melt parameterisation yields NaN for a zero slope
    ! END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE find_effective_plume_path_GL_average

  SUBROUTINE find_effective_plume_path_along_ice_flow( mesh, ice, vi, eff_plume_source_depth, eff_basal_slope)
    ! Find the effective plume source depth and basal slope for shelf grid cell [i,j] by
    ! assuming plumes follow the same horizontal flow field as the ice shelf

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp),                            INTENT(OUT)   :: eff_plume_source_depth
    REAL(dp),                            INTENT(OUT)   :: eff_basal_slope

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_effective_plume_path_along_ice_flow'
    REAL(dp), DIMENSION(2)                             :: t1 ,t2 ,uv
    INTEGER                                            :: nit, nit_max
    REAL(dp)                                           :: u1, v1, TAF1, TAF2, depth1, Hi2, Hs2, depth2, lambda, GLx, GLy
    LOGICAL                                            :: found_GL, got_stuck
    REAL(dp)                                           :: zb_shelf, dist

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Calculate the shelf base depth (draft)
    ! zb_shelf = ice%Hs_a( vi) - ice%Hi_a( vi)

    ! ! Track a tracer upstream along the ice flow field until grounded ice is found
    ! t1     = [grid%x( i), grid%y( j)]
    ! t2     = t1
    ! TAF1   = ice%TAF_a( j,i)
    ! depth1 = zb_shelf

    ! found_GL  = .FALSE.
    ! got_stuck = .FALSE.
    ! nit       = 0
    ! nit_max   = grid%nx + grid%ny
    ! DO WHILE (.NOT. found_GL)

    !   ! Safety
    !   nit = nit + 1
    !   IF (nit > nit_max) THEN
    !     got_stuck = .TRUE.
    !     EXIT
    !     !WRITE(0,*) '  find_effective_plume_path_along_ice_flow - ERROR: tracer got stuck!'
    !     !CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    !   END IF

    !   ! Interpolate velocity field to exact tracer location
    !   u1 = interp_bilin_2D( ice%u_vav_a, grid%x, grid%y, t1( 1), t1( 2))
    !   v1 = interp_bilin_2D( ice%v_vav_a, grid%x, grid%y, t1( 1), t1( 2))

    !   ! Move the tracer upstream
    !   uv = [u1, v1]
    !   uv = uv / NORM2( uv)
    !   t2 = t1 - uv * grid%dx

    !   ! Safety
    !   IF (t2( 1) < grid%xmin .OR. t2( 1) > grid%xmax .OR. t2( 2) < grid%ymin .OR. t2( 2) > grid%ymax) THEN
    !     got_stuck = .TRUE.
    !     EXIT
    !     !WRITE(0,*) '  find_effective_plume_path_along_ice_flow - ERROR: tracer outside domain!'
    !     !CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    !   END IF

    !   ! Interpolate data to new tracer location
    !   TAF2   = interp_bilin_2D( ice%TAF_a, grid%x, grid%y, t2( 1), t2( 2))
    !   Hi2    = interp_bilin_2D( ice%Hi_a,  grid%x, grid%y, t2( 1), t2( 2))
    !   Hs2    = interp_bilin_2D( ice%Hs_a,  grid%x, grid%y, t2( 1), t2( 2))
    !   depth2 = Hs2 - Hi2

    !   ! Check if we've found grounded ice
    !   IF (TAF2 > 0._dp) THEN

    !     found_GL = .TRUE.

    !     ! Find exact GL depth
    !     lambda = TAF1 / (TAF1 - TAF2)
    !     eff_plume_source_depth = lambda * depth2 + (1._dp - lambda) * depth1
    !     GLx                    = lambda * t2( 1) + (1._dp - lambda) * t1( 1)
    !     GLy                    = lambda * t2( 2) + (1._dp - lambda) * t1( 2)

    !     ! Calculate the plume source depth and basal slope
    !     dist  = SQRT( (GLx - grid%x( i))**2 + (GLy - grid%y( j))**2 )
    !     eff_basal_slope = (zb_shelf - eff_plume_source_depth) / dist

    !   ELSE
    !     ! Cycle the tracer

    !     t1     = t2
    !     TAF1   = TAF2
    !     depth1 = depth2

    !   END IF ! IF (TAF2 > 0._dp) THEN

    ! END DO ! DO WHILE (.NOT. found_GL)

    ! ! Exception for when the GL source is less deep than the shelf base
    ! IF (eff_plume_source_depth >= zb_shelf .OR. got_stuck) THEN
    !   eff_plume_source_depth = zb_shelf
    !   eff_basal_slope        = 1E-10_dp   ! Because the melt parameterisation yields NaN for a zero slope
    ! END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE find_effective_plume_path_along_ice_flow

  FUNCTION Lazeroms2018_dimensionless_melt_curve( xhat) RESULT( Mhat)
    ! The dimensionless melt curve from Lazeroms et al. (2018), Appendix A, Eq. A13

    IMPLICIT NONE

    REAL(dp), INTENT(IN)                                 :: xhat                        ! Scaled distance  from grounding line
    REAL(dp)                                             :: Mhat                        ! Scaled melt rate from polynomial fit

    ! ! L variables: polynomial coefficients
    ! REAL(dp), PARAMETER                                  :: p11 =  6.387953795485420E4_dp
    ! REAL(dp), PARAMETER                                  :: p10 = -3.520598035764990E5_dp
    ! REAL(dp), PARAMETER                                  :: p9  =  8.466870335320488E5_dp
    ! REAL(dp), PARAMETER                                  :: p8  = -1.166290429178556E6_dp
    ! REAL(dp), PARAMETER                                  :: p7  =  1.015475347943186E6_dp
    ! REAL(dp), PARAMETER                                  :: p6  = -5.820015295669482E5_dp
    ! REAL(dp), PARAMETER                                  :: p5  =  2.218596970948727E5_dp
    ! REAL(dp), PARAMETER                                  :: p4  = -5.563863123811898E4_dp
    ! REAL(dp), PARAMETER                                  :: p3  =  8.927093637594877E3_dp
    ! REAL(dp), PARAMETER                                  :: p2  = -8.951812433987858E2_dp
    ! REAL(dp), PARAMETER                                  :: p1  =  5.527656234709359E1_dp
    ! REAL(dp), PARAMETER                                  :: p0  =  0.1371330075095435_dp

    ! Mhat = p11 * xhat**11 + &
    !        p10 * xhat**10 + &
    !        p9  * xhat**9  + &
    !        p8  * xhat**8  + &
    !        p7  * xhat**7  + &
    !        p6  * xhat**6  + &
    !        p5  * xhat**5  + &
    !        p4  * xhat**4  + &
    !        p3  * xhat**3  + &
    !        p2  * xhat**2  + &
    !        p1  * xhat + p0

  END FUNCTION Lazeroms2018_dimensionless_melt_curve

  SUBROUTINE initialise_BMB_model_Lazeroms2018_plume( mesh, BMB)
    ! Allocate memory for the data fields of the Favier et al. (2019) shelf BMB parameterisations.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_Lazeroms2018_plume'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not yet implemented for a triangular mesh!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ! ! Variables
    ! CALL allocate_shared_dp_1D(  mesh%nV,      BMB%T_ocean_base,           BMB%wT_ocean_base          )
    ! CALL allocate_shared_dp_1D(  mesh%nV,      BMB%eff_plume_source_depth, BMB%weff_plume_source_depth)
    ! CALL allocate_shared_dp_1D(  mesh%nV,      BMB%eff_basal_slope,        BMB%weff_basal_slope       )
    ! CALL allocate_shared_int_2D( 16,      2,   BMB%search_directions,      BMB%wsearch_directions     )

    ! ! Define the 16 search directions
    ! IF (par%master) THEN
    !   BMB%search_directions(  1,:) = (/  1,  0 /)
    !   BMB%search_directions(  2,:) = (/  2, -1 /)
    !   BMB%search_directions(  3,:) = (/  1, -1 /)
    !   BMB%search_directions(  4,:) = (/  1, -2 /)
    !   BMB%search_directions(  5,:) = (/  0, -1 /)
    !   BMB%search_directions(  6,:) = (/ -1, -2 /)
    !   BMB%search_directions(  7,:) = (/ -1, -1 /)
    !   BMB%search_directions(  8,:) = (/ -2, -1 /)
    !   BMB%search_directions(  9,:) = (/ -1,  0 /)
    !   BMB%search_directions( 10,:) = (/ -2,  1 /)
    !   BMB%search_directions( 11,:) = (/ -1,  1 /)
    !   BMB%search_directions( 12,:) = (/ -1,  2 /)
    !   BMB%search_directions( 13,:) = (/  0,  1 /)
    !   BMB%search_directions( 14,:) = (/  1,  2 /)
    !   BMB%search_directions( 15,:) = (/  1,  1 /)
    !   BMB%search_directions( 16,:) = (/  2,  1 /)
    ! END IF
    ! CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model_Lazeroms2018_plume

! ===== The PICO ocean box model =====
! ====================================

  SUBROUTINE run_BMB_model_PICO( mesh, ice, ocean, BMB)
    ! Calculate basal melt using the PICO ocean box model

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_PICO'
    INTEGER                                            :: basin_i

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not yet implemented for a triangular mesh!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ! ! Assign PICO ocean boxes to all shelf grid cells in all basins
    ! CALL PICO_assign_ocean_boxes( grid, ice, BMB)

    ! ! Run PICO for all basins
    ! DO basin_i = 1, ice%nbasins
    !   CALL run_BMB_model_PICO_basin( grid, ice, ocean, BMB, basin_i)
    ! END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_PICO

  SUBROUTINE PICO_assign_ocean_boxes( mesh, ice, BMB)
    ! Assign PICO ocean boxes to shelf grid cells using the distance-to-grounding-line / distance-to-calving-front
    ! approach outlined in Reese et al. (2018)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'PICO_assign_ocean_boxes'
    INTEGER                                            :: i,j,k
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_GL_D
    INTEGER                                            :: wd_GL_D
    REAL(dp)                                           :: d_max
    INTEGER                                            :: basin_i
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: n_cells_per_box
    LOGICAL                                            :: do_reduce_n_D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Determine number of PICO boxes to be used for each basin
    ! ! ========================================================

    !   CALL allocate_shared_dp_1D(  ice%nbasins, d_GL_D, wd_GL_D)

    !   ! Determine relative distance to grounding line for all shelf grid cells in the entire model domain
    !   DO i = grid%i1, grid%i2
    !   DO j = 1, grid%ny
    !     BMB%PICO_d_GL( j,i) = 0._dp
    !     BMB%PICO_d_IF( j,i) = 0._dp
    !     BMB%PICO_r(    j,i) = 0._dp
    !     IF (ice%mask_shelf_a( j,i) == 1) CALL calc_dGL_dIF_r( grid, ice, BMB, i,j, BMB%PICO_d_GL( j,i), BMB%PICO_d_IF( j,i), BMB%PICO_r( j,i))
    !   END DO
    !   END DO
    !   CALL sync

    !   ! Calculate maximum distance to grounding line within each basin
    !   DO basin_i = 1, ice%nbasins
    !     d_max = 0._dp
    !     DO i = grid%i1, grid%i2
    !     DO j = 1, grid%ny
    !       IF (ice%basin_ID( j,i) == basin_i) THEN
    !         d_max = MAX( d_max, BMB%PICO_d_GL( j,i))
    !       END IF
    !     END DO
    !     END DO
    !     CALL MPI_ALLREDUCE( MPI_IN_PLACE, d_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    !     IF (par%master) d_GL_D( basin_i) = d_max
    !     CALL sync
    !   END DO ! DO basin_i = 1, ice%nbasins

    !   ! Calculate maximum distance to grounding line within the entire model doman
    !   IF (par%master) d_max = MAXVAL( d_GL_D)
    !   CALL MPI_BCAST( d_max, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    !   ! Determine number of PICO boxes for each basin
    !   IF (par%master) THEN
    !     DO basin_i = 1, ice%nbasins
    !       ! Reese et al. (2018), Eq. 9
    !       BMB%PICO_n_D( basin_i) = 1 + INT( SQRT( d_GL_D( basin_i) / d_max) * REAL( C%BMB_PICO_nboxes - 1,dp))
    !     END DO
    !   END IF
    !   CALL sync

    ! ! Assign PICO boxes to all shelf grid cells in all basins
    ! ! =======================================================

    !   ALLOCATE( n_cells_per_box( C%BMB_PICO_nboxes))

    !   DO basin_i = 1, ice%nbasins

    !     ! If necessary, reduce the maximum number of boxes for a basin
    !     ! until every box is assigned at least one grid cell

    !     n_cells_per_box = 0

    !     DO WHILE (.TRUE.)

    !       ! Assign ocean boxes according to Reese et al. (2018), Eq. 11
    !       DO i = grid%i1, grid%i2
    !       DO j = 1, grid%ny
    !         IF (ice%basin_ID( j,i) == basin_i) THEN
    !           BMB%PICO_k( j,i) = 0
    !           IF (ice%mask_shelf_a( j,i) == 1) THEN
    !             DO k = 1, BMB%PICO_n_D( basin_i)
    !                 IF (1._dp - SQRT( REAL(BMB%PICO_n_D( basin_i) - k + 1, dp) / REAL( BMB%PICO_n_D( basin_i), dp)) <= BMB%PICO_r( j,i) .AND. &
    !                     1._dp - SQRT( REAL(BMB%PICO_n_D( basin_i) - k    , dp) / REAL( BMB%PICO_n_D( basin_i), dp)) >= BMB%PICO_r( j,i)) THEN
    !                 BMB%PICO_k( j,i) = k
    !                 n_cells_per_box( k) = n_cells_per_box( k) + 1
    !               END IF
    !             END DO ! DO k = 1, BMB%PICO_n_D( basin_i)
    !           END IF ! IF (ice%mask_shelf_a( j,i) == 1) THEN
    !         END IF ! IF (ice%basin_ID( j,i) == basin_i) THEN
    !       END DO
    !       END DO
    !       CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_cells_per_box, C%BMB_PICO_nboxes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    !       IF (par%master) BMB%PICO_A( basin_i,:) = n_cells_per_box * grid%dx**2
    !       CALL sync

    !       ! If any of the boxes have zero grid cells assigned, reduce the maximum number of boxes and try again
    !       do_reduce_n_D = .FALSE.
    !       DO k = 1, BMB%PICO_n_D( basin_i)
    !         IF (n_cells_per_box( k) == 0) THEN
    !           do_reduce_n_D = .TRUE.
    !         END IF
    !       END DO
    !       IF (do_reduce_n_D) THEN
    !         IF (par%master) BMB%PICO_n_D( basin_i) = BMB%PICO_n_D( basin_i) - 1
    !         CALL sync
    !       ELSE
    !         EXIT
    !       END IF

    !     END DO ! DO WHILE (.TRUE.)

    !   END DO ! DO basin_i = 1, ice%nbasins

    !   ! Clean up after yourself
    !   DEALLOCATE( n_cells_per_box)
    !   CALL deallocate_shared( wd_GL_D)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE PICO_assign_ocean_boxes

  SUBROUTINE run_BMB_model_PICO_basin( mesh, ice, ocean, BMB, basin_i)
    ! Calculate basal melt for ice basin i using the PICO ocean box model

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    INTEGER,                             INTENT(IN)    :: basin_i

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_PICO_basin'
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: Tk0,Sk0
    REAL(dp)                                           :: nu,lambda
    REAL(dp)                                           :: g1,g2,s,Crbsa,Tstar,x,y
    REAL(dp)                                           :: q

    ! Constants (Reese et al. (2018), Table 1)
    REAL(dp), PARAMETER                                :: aa          = -0.0572_dp     ! Salinity coefficient of freezing equation         [degC PSU^-1]
    REAL(dp), PARAMETER                                :: bb          =  0.0788_dp     ! Constant coefficient of freezing equation         [degC]
    REAL(dp), PARAMETER                                :: cc          = 7.77E-8_dp     ! Pressure coefficient of freezing equation         [degC Pa^-1]
    REAL(dp), PARAMETER                                :: alpha       =  7.5E-5_dp     ! Thermal expansion coefficient in EOS              [degC^-1]
    REAL(dp), PARAMETER                                :: beta        =  7.7E-4_dp     ! Salt contraction coefficient in EOS               [PSU^-1]
    REAL(dp), PARAMETER                                :: rhostar     = 1033_dp        ! Reference density in EOS                          [kg m^-3]
    REAL(dp), PARAMETER                                :: C_overturn  = 1.0E6_dp       ! Overturning strength                              [m^6 s^-1 kg^-1]

    ! Add routine to path
    CALL init_routine( routine_name)

    !   ! Initialise
    !   DO i = grid%i1, grid%i2
    !   DO j = 1, grid%ny
    !     IF (ice%basin_ID( j,i) == basin_i) THEN
    !       BMB%PICO_T( j,i) = 0._dp
    !       BMB%PICO_S( j,i) = 0._dp
    !       BMB%PICO_m( j,i) = 0._dp
    !     END IF
    !   END DO
    !   END DO
    !   CALL sync

    !   ! Some intermediary constants (Reese et al. (2018), just after Eq. A2)
    !   nu     = ice_density / seawater_density
    !   lambda = L_fusion / cp_ocean

    !   ! Calculate temperature and salinity in box B0 for this basin
    !   CALL PICO_calc_T0_S0( grid, ice, ocean, basin_i, Tk0, Sk0)

    !   ! Calculate 2-D + box-averaged basal pressures
    !   BMB%PICO_p( :,grid%i1:grid%i2) = ice_density * grav * ice%Hi_a( :,grid%i1:grid%i2)
    !   DO k = 1, BMB%PICO_n_D( basin_i)
    !     CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_p, basin_i, k, BMB%PICO_pk( basin_i, k))
    !   END DO

    ! ! Calculate solution for box 1
    ! ! ============================

    !   DO i = grid%i1, grid%i2
    !   DO j = 1, grid%ny

    !     IF (ice%basin_ID( j,i) == basin_i .AND. BMB%PICO_k( j,i) == 1) THEN

    !       ! Reese et al. (2018), just before Eq. A6
    !       g1 = BMB%PICO_A( basin_i, 1) * C%BMB_PICO_GammaTstar
    !       g2 = g1 / (nu * lambda)
    !       Tstar = aa * Sk0 + bb - cc * BMB%PICO_pk( basin_i, 1) - Tk0

    !       ! Reese et al. (2018), just after Eq. A11
    !       s = Sk0 / (nu * lambda)

    !       ! Intermediary constants
    !       Crbsa = C_overturn * rhostar * (beta * s - alpha)

    !       ! Reese et al. (2018), Eq. A12
    !       x = -g1 / (2._dp * Crbsa) + SQRT( (g1 / (2._dp * Crbsa))**2 - (g1 * Tstar / Crbsa))

    !       ! Reese et al. (2018), Eq. A8
    !       y = Sk0 * x / (nu * lambda)

    !       BMB%PICO_T( j,i) = Tk0 - x
    !       BMB%PICO_S( j,i) = Sk0 - y

    !       ! Reese et al. (2019), Eq. 13
    !       BMB%PICO_m( j,i) = sec_per_year * C%BMB_PICO_GammaTstar / (nu*lambda) * (aa * BMB%PICO_S( j,i) + bb - cc * BMB%PICO_p( j,i) - BMB%PICO_T( j,i))

    !     END IF ! IF (BMB%PICO_k( j,i) == 1) THEN

    !   END DO
    !   END DO
    !   CALL sync

    !   ! Calculate box-averaged values
    !   CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_T, basin_i, 1, BMB%PICO_Tk( basin_i, 1))
    !   CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_S, basin_i, 1, BMB%PICO_Sk( basin_i, 1))
    !   CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_m, basin_i, 1, BMB%PICO_mk( basin_i, 1))

    !   ! Calculate overturning strength (Reese et al. (2018), Eq. A9)
    !   q = C_overturn * rhostar * (beta * (Sk0 - BMB%PICO_Sk( basin_i, 1)) - alpha * (Tk0 - BMB%PICO_Tk( basin_i, 1)))

    ! ! Calculate solutions for subsequent boxes
    ! ! ========================================

    !   DO k = 2, BMB%PICO_n_D( basin_i)

    !     DO i = grid%i1, grid%i2
    !     DO j = 1, grid%ny

    !       IF (ice%basin_ID( j,i) == basin_i .AND. BMB%PICO_k( j,i) == k) THEN

    !         ! Reese et al. (2018), just before Eq. A6
    !         g1 = BMB%PICO_A( basin_i, k) * C%BMB_PICO_GammaTstar
    !         g2 = g1 / (nu * lambda)
    !         !Tstar = aa * Sk0 + bb - cc * BMB%PICO_pk( basin_i, k-1) - BMB%PICO_Tk( basin_i, k-1)
    !         Tstar = aa * BMB%PICO_Sk( basin_i, k-1) + bb - cc * BMB%PICO_pk( basin_i, k) - BMB%PICO_Tk( basin_i, k-1)

    !         ! Reese et al. (2018), Eq. A13
    !         x = -g1 * Tstar / (q + g1 - g2 * aa * BMB%PICO_Sk( basin_i, k-1))

    !         ! Reese et al. (2018), Eq. A8
    !         y = BMB%PICO_Sk( basin_i, k-1) * x / (nu * lambda)

    !         BMB%PICO_T( j,i) = BMB%PICO_Tk( basin_i, k-1) - x
    !         BMB%PICO_S( j,i) = BMB%PICO_Sk( basin_i, k-1) - y

    !         ! Reese et al. (2019), Eq. 13
    !         BMB%PICO_m( j,i) = sec_per_year * C%BMB_PICO_GammaTstar / (nu*lambda) * (aa * BMB%PICO_S( j,i) + bb - cc * BMB%PICO_p( j,i) - BMB%PICO_T( j,i))

    !       END IF ! IF (BMB%PICO_k( j,i) == k) THEN

    !     END DO
    !     END DO
    !     CALL sync

    !     ! Calculate box-averaged values
    !     CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_T, basin_i, k, BMB%PICO_Tk( basin_i, k))
    !     CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_S, basin_i, k, BMB%PICO_Sk( basin_i, k))
    !     CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_m, basin_i, k, BMB%PICO_mk( basin_i, k))

    !   END DO ! DO k = 2, BMB%PICO_n_D( basin_i)

    !   ! Copy melt rates to final data field
    !   DO i = grid%i1, grid%i2
    !   DO j = 1, grid%ny
    !     IF (ice%basin_ID( j,i) == basin_i) BMB%BMB_shelf( j,i) = BMB%PICO_m( j,i)
    !   END DO
    !   END DO
    !   CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_PICO_basin

  SUBROUTINE calc_dGL_dIF_r( mesh, ice, BMB, vi, d_GL, d_IF, r)
    ! For each shelf grid cell, calculate the distance to the grounding line dGL,
    ! the distance to the ice front dIF, and the relative distance r (Reese et al. (2018), Eq. 10)
    !
    ! Determines d_GL and d_IF using the 16-directions search scheme.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp),                            INTENT(OUT)   :: d_GL, d_IF, r

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_dGL_dIF_r'
    INTEGER                                            :: n
    INTEGER                                            :: dpi,dpj,ip1,jp1,ip2,jp2
    REAL(dp)                                           :: dist
    LOGICAL                                            :: reached_end

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Exception for when this grid cell isn't shelf
    ! IF (ice%mask_shelf_a( j,i) == 0) THEN
    !   d_GL = 0._dp
    !   d_IF = 0._dp
    !   r    = 0._dp
    !   CALL finalise_routine( routine_name)
    !   RETURN
    ! END IF

    ! ! Initialise
    ! d_GL = REAL( MAX( grid%ny, grid%nx), dp) * grid%dx
    ! d_IF = d_GL

    ! ! Investigate all 16 search directions
    ! DO n = 1, 16

    !   ! The search direction vector
    !   dpi = BMB%search_directions( n,1)
    !   dpj = BMB%search_directions( n,2)

    !   ! Initialise the search pointer at the shelf grid cell
    !   ip1 = i
    !   jp1 = j
    !   ip2 = i + dpi
    !   jp2 = j + dpj

    !   ! If the search direction already points out of the model domain, don't bother
    !   IF (ip2 < 1 .OR. ip2 > grid%nx .OR. jp2 < 1 .OR. jp2 > grid%ny) CYCLE

    !   ! Search in this direction
    !   reached_end  = .FALSE.

    !   DO WHILE (.NOT. reached_end)

    !     ! If the pointer exits the model domain, stop the search
    !     IF (ip2 < 1 .OR. ip2 > grid%nx .OR. jp2 < 1 .OR. jp2 > grid%ny) THEN
    !       reached_end  = .TRUE.
    !       EXIT
    !     END IF

    !     ! If the pointer encounters open ocean, stop the search and update d_IF
    !     IF (ice%mask_sheet_a( jp2,ip2) == 1) THEN
    !       reached_end  = .TRUE.
    !       dist = SQRT( REAL( ip2 - i,dp)**2 + REAL( jp2 - j,dp)**2) * grid%dx
    !       d_GL = MIN( d_GL, dist)
    !       EXIT
    !     END IF

    !     ! If the pointer encounters open ocean, stop the search and update d_IF
    !     IF (ice%mask_ocean_a( jp2,ip2) == 1 .AND. ice%mask_shelf_a( jp2,ip2) == 0) THEN
    !       reached_end  = .TRUE.
    !       dist = SQRT( REAL( ip2 - i,dp)**2 + REAL( jp2 - j,dp)**2) * grid%dx
    !       d_IF = MIN( d_IF, dist)
    !       EXIT
    !     END IF

    !     ! If none of these exceptions were triggered, advance the pointer in the search direction
    !     ip1 = ip2
    !     jp1 = jp2
    !     ip2 = ip1 + dpi
    !     jp2 = jp1 + dpj

    !   END DO ! DO WHILE (.NOT. reached_end)

    ! END DO ! DO n = 1, 16

    ! ! Reese et al. (2018), Eq. 10
    ! r = d_GL / (d_GL + d_IF)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dGL_dIF_r

  SUBROUTINE PICO_calc_T0_S0( mesh, ice, ocean, basin_i, Tk0, Sk0)
    ! Find temperature and salinity in box B0 (defined as mean ocean-floor value at the calving front)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    INTEGER,                             INTENT(IN)    :: basin_i
    REAL(dp),                            INTENT(OUT)   :: Tk0,Sk0

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'PICO_calc_T0_S0'
    INTEGER                                            :: vi,n
    REAL(dp)                                           :: depth, T_floor, S_floor, depth_max
    INTEGER                                            :: ii,jj

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Average ocean-floor temperature and salinity over this basin's ocean-next-to-floating-ice pixels
    ! n   = 0
    ! Tk0 = 0._dp
    ! Sk0 = 0._dp

    ! DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    ! DO j = 2, grid%ny-1

    !   IF (ice%basin_ID( j,i) == basin_i .AND. ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_ice_a( j,i) == 0) THEN
    !     IF (ice%mask_shelf_a( j-1,i-1) == 1 .OR. &
    !         ice%mask_shelf_a( j-1,i  ) == 1 .OR. &
    !         ice%mask_shelf_a( j-1,i+1) == 1 .OR. &
    !         ice%mask_shelf_a( j  ,i-1) == 1 .OR. &
    !         ice%mask_shelf_a( j  ,i+1) == 1 .OR. &
    !         ice%mask_shelf_a( j+1,i-1) == 1 .OR. &
    !         ice%mask_shelf_a( j+1,i  ) == 1 .OR. &
    !         ice%mask_shelf_a( j+1,i+1) == 1) THEN
    !       ! This pixel is open ocean next to floating ice

    !       ! Find ocean-floor temperature and salinity
    !       depth = MAX( 0.1_dp, -ice%Hb_a( j,i))
    !       CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( :,j,i), depth, T_floor)
    !       CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S_ocean_corr_ext( :,j,i), depth, S_floor)

    !       ! Add to sum
    !       n   = n   + 1
    !       Tk0 = Tk0 + T_floor
    !       Sk0 = Sk0 + S_floor

    !     END IF
    !   END IF

    ! END DO
    ! END DO
    ! CALL sync

    ! ! Combine results from processes
    ! CALL MPI_ALLREDUCE( MPI_IN_PLACE, n,   1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    ! CALL MPI_ALLREDUCE( MPI_IN_PLACE, Tk0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    ! CALL MPI_ALLREDUCE( MPI_IN_PLACE, Sk0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Tk0 = Tk0 / REAL(n,dp)
    ! Sk0 = Sk0 / REAL(n,dp)

    ! ! Safety
    ! IF (n == 0) THEN
    !   ! No ocean-next-to-shelf grid cells could be found within this basin;
    !   ! instead, just take the value from the deepest ocean floor grid cell (which must by definition be ice-covered...)

    !   IF (par%master) THEN

    !     depth_max = 0._dp
    !     ii = 0
    !     jj = 0

    !     DO i = 1, grid%nx
    !     DO j = 1, grid%ny
    !       IF (ice%mask_ocean_a( j,i) == 1) THEN
    !         depth = -ice%Hb_a( j,i)
    !         IF (depth > depth_max) THEN
    !           depth_max = depth
    !           ii = i
    !           jj = j
    !         END IF
    !       END IF
    !     END DO
    !     END DO

    !     IF (ii == 0 .OR. jj == 0) THEN
    !       WRITE(0,*) '  PICO_calc_T0_S0 - ERROR: couldnt find deepest ocean floor grid cell!'
    !       CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    !     END IF

    !     ! Find ocean-floor temperature and salinity
    !     CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( :,jj,ii), depth_max, Tk0)
    !     CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S_ocean_corr_ext( :,jj,ii), depth_max, Sk0)

    !   END IF ! IF (par%master) THEN

    !   CALL MPI_BCAST( Tk0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    !   CALL MPI_BCAST( Sk0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! END IF ! IF (n == 0) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE PICO_calc_T0_S0

  SUBROUTINE PICO_calc_box_average( mesh, ice, BMB, d, basin_i, k, d_av)
    ! Calculate the average d_av of field d over ocean box k in basin i

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp), DIMENSION(:   ),           INTENT(IN)    :: d
    INTEGER,                             INTENT(IN)    :: basin_i
    INTEGER,                             INTENT(IN)    :: k
    REAL(dp),                            INTENT(OUT)   :: d_av

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'PICO_calc_box_average'
    INTEGER                                            :: vi,n
    REAL(dp)                                           :: d_sum

    ! Add routine to path
    CALL init_routine( routine_name)

    ! n     = 0
    ! d_sum = 0._dp

    ! DO vi = mesh%vi1, mesh%vi2
    !   IF (ice%basin_ID( vi) == basin_i .AND. BMB%PICO_k( vi) == k) THEN
    !     n     = n     + 1
    !     d_sum = d_sum + d( vi)
    !   END IF
    ! END DO

    ! CALL MPI_ALLREDUCE( MPI_IN_PLACE, n,     1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    ! CALL MPI_ALLREDUCE( MPI_IN_PLACE, d_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! d_av = d_sum / REAL(n,dp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE PICO_calc_box_average

  SUBROUTINE initialise_BMB_model_PICO( mesh, ice, BMB)
    ! Allocate memory for the data fields of the PICO ocean box model

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_PICO'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not yet implemented for a triangular mesh!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ! ! Variables
    ! CALL allocate_shared_int_2D( 16,          2,                     BMB%search_directions, BMB%wsearch_directions)
    ! CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes,     BMB%PICO_A,            BMB%wPICO_A           )
    ! CALL allocate_shared_int_1D( ice%nbasins,                        BMB%PICO_n_D,          BMB%wPICO_n_D         )

    ! CALL allocate_shared_dp_1D(  mesh%nV,     BMB%PICO_d_GL,         BMB%wPICO_d_GL)
    ! CALL allocate_shared_dp_1D(  mesh%nV,     BMB%PICO_d_IF,         BMB%wPICO_d_IF)
    ! CALL allocate_shared_dp_1D(  mesh%nV,     BMB%PICO_r,            BMB%wPICO_r   )
    ! CALL allocate_shared_int_1D( mesh%nV,     BMB%PICO_k,            BMB%wPICO_k   )

    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_T,  BMB%wPICO_T )
    ! CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_Tk, BMB%wPICO_Tk)
    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_S,  BMB%wPICO_S )
    ! CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_Sk, BMB%wPICO_Sk)
    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_p,  BMB%wPICO_p )
    ! CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_pk, BMB%wPICO_pk)
    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_m,  BMB%wPICO_m )
    ! CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_mk, BMB%wPICO_mk)

    ! ! Define the 16 search directions
    ! IF (par%master) THEN
    !   BMB%search_directions(  1,:) = (/  1,  0 /)
    !   BMB%search_directions(  2,:) = (/  2, -1 /)
    !   BMB%search_directions(  3,:) = (/  1, -1 /)
    !   BMB%search_directions(  4,:) = (/  1, -2 /)
    !   BMB%search_directions(  5,:) = (/  0, -1 /)
    !   BMB%search_directions(  6,:) = (/ -1, -2 /)
    !   BMB%search_directions(  7,:) = (/ -1, -1 /)
    !   BMB%search_directions(  8,:) = (/ -2, -1 /)
    !   BMB%search_directions(  9,:) = (/ -1,  0 /)
    !   BMB%search_directions( 10,:) = (/ -2,  1 /)
    !   BMB%search_directions( 11,:) = (/ -1,  1 /)
    !   BMB%search_directions( 12,:) = (/ -1,  2 /)
    !   BMB%search_directions( 13,:) = (/  0,  1 /)
    !   BMB%search_directions( 14,:) = (/  1,  2 /)
    !   BMB%search_directions( 15,:) = (/  1,  1 /)
    !   BMB%search_directions( 16,:) = (/  2,  1 /)
    ! END IF
    ! CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model_PICO

! ===== The PICOP ocean box + plume model =====
! =============================================

  SUBROUTINE run_BMB_model_PICOP( mesh, ice, ocean, BMB)
    ! Calculate basal melt using the PICOP ocean box + plume model

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_PICOP'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not yet implemented for a triangular mesh!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ! ! First run the PICO ocean box model to determine the temperature and salinity in the cavity
    ! CALL run_BMB_model_PICO( mesh, ice, ocean, BMB)

    ! ! Then run the Lazeroms (2018) plume parameterisation to calculate melt rates
    ! CALL run_BMB_model_Lazeroms2018_plume( mesh, ice, ocean, BMB)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_PICOP

  SUBROUTINE initialise_BMB_model_PICOP( mesh, ice, BMB)
    ! Allocate memory for the data fields of the PICOP ocean box + plume model

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_PICOP'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not yet implemented for a triangular mesh!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ! ! Variables
    ! CALL allocate_shared_dp_1D(  mesh%nV,    BMB%T_ocean_base,           BMB%wT_ocean_base          )
    ! CALL allocate_shared_int_2D( 16,      2, BMB%search_directions,      BMB%wsearch_directions     )
    ! CALL allocate_shared_dp_1D(  mesh%nV,    BMB%eff_plume_source_depth, BMB%weff_plume_source_depth)
    ! CALL allocate_shared_dp_1D(  mesh%nV,    BMB%eff_basal_slope,        BMB%weff_basal_slope       )

    ! CALL allocate_shared_int_2D( 16,          2,                 BMB%search_directions, BMB%wsearch_directions)
    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_d_GL,         BMB%wPICO_d_GL        )
    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_d_IF,         BMB%wPICO_d_IF        )
    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_r,            BMB%wPICO_r           )
    ! CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_A,            BMB%wPICO_A           )
    ! CALL allocate_shared_int_1D( ice%nbasins,                    BMB%PICO_n_D,          BMB%wPICO_n_D         )
    ! CALL allocate_shared_int_1D( mesh%nV,                        BMB%PICO_k,            BMB%wPICO_k           )

    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_T,  BMB%wPICO_T )
    ! CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_Tk, BMB%wPICO_Tk)
    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_S,  BMB%wPICO_S )
    ! CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_Sk, BMB%wPICO_Sk)
    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_p,  BMB%wPICO_p )
    ! CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_pk, BMB%wPICO_pk)
    ! CALL allocate_shared_dp_1D(  mesh%nV,                        BMB%PICO_m,  BMB%wPICO_m )
    ! CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_mk, BMB%wPICO_mk)

    ! ! Define the 16 search directions
    ! IF (par%master) THEN
    !   BMB%search_directions(  1,:) = (/  1,  0 /)
    !   BMB%search_directions(  2,:) = (/  2, -1 /)
    !   BMB%search_directions(  3,:) = (/  1, -1 /)
    !   BMB%search_directions(  4,:) = (/  1, -2 /)
    !   BMB%search_directions(  5,:) = (/  0, -1 /)
    !   BMB%search_directions(  6,:) = (/ -1, -2 /)
    !   BMB%search_directions(  7,:) = (/ -1, -1 /)
    !   BMB%search_directions(  8,:) = (/ -2, -1 /)
    !   BMB%search_directions(  9,:) = (/ -1,  0 /)
    !   BMB%search_directions( 10,:) = (/ -2,  1 /)
    !   BMB%search_directions( 11,:) = (/ -1,  1 /)
    !   BMB%search_directions( 12,:) = (/ -1,  2 /)
    !   BMB%search_directions( 13,:) = (/  0,  1 /)
    !   BMB%search_directions( 14,:) = (/  1,  2 /)
    !   BMB%search_directions( 15,:) = (/  1,  1 /)
    !   BMB%search_directions( 16,:) = (/  2,  1 /)
    ! END IF
    ! CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model_PICOP

! ===== The Bernales et al. (202X) sub-shelf parameterisation =====
! =================================================================

  SUBROUTINE run_BMB_model_Bernales202X( mesh, ice, ocean, BMB, refgeo)
    ! Calculate sub-shelf melt with Bernales et al. (202X) parameterisation

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_Bernales202X'
    INTEGER                                            :: vi
    REAL(dp)                                           :: t_melt, t_force, q_factor
    REAL(dp)                                           :: qsh_g, qsh_w, dist_gl, gl_w
    REAL(dp), PARAMETER                                :: gamma_t = 1.0E-04_dp
    REAL(dp), DIMENSION(:), POINTER                    ::  F_melt
    INTEGER                                            :: wF_melt

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise ice shelf basal melt
    BMB%BMB_shelf = 0._dp

    ! Allocate melt factor
    CALL allocate_shared_dp_1D( mesh%nV, F_melt, wF_melt)

    ! Initialise melt factor
    F_melt =  0.0_dp

    ! === Basal ocean temperatures [K] ===
    ! ====================================

    IF (.NOT. C%do_ocean_inv) THEN
      ! Calculate ocean temperature at the base of the ice shelf
      CALL calc_ocean_temperature_at_shelf_base( mesh, ice, ocean, BMB)
    END IF

    ! Calculate ocean salinity at the base of the ice shelf
    CALL calc_ocean_salinity_at_shelf_base( mesh, ice, ocean, BMB)

    ! === Melt factor ===
    ! ===================

    ! Loop over all ice shelf and ocean points
    DO vi = mesh%vi1, mesh%vi2

      IF (ice%mask_shelf_a( vi) == 1 .OR. &
          ice%mask_ocean_a( vi) == 1) THEN

        ! = Flux through floating ice [weight]
        ! ====================================

        qsh_g = ice%uabs_vav_a( vi) * ice%Hi_a( vi)

        qsh_w = (2.0_dp/pi) * ATAN( qsh_g**2.0_dp / 2.0E5_dp**2.0_dp)

        ! = Distance to grounding line [weight]
        ! =====================================

        CALL calc_distance_to_grounding_line( mesh, ice, vi, dist_gl)

        gl_w = EXP(-dist_gl / 20000._dp)

        ! = Weighed-averaged F_melt [K^-1]
        ! ================================

        F_melt( vi) = 1.0E-01_dp * (1.0_dp - qsh_w * gl_w) + &
                      1.0E-00_dp *           qsh_w * gl_w

        ! F_melt( vi) = 1.0E-01_dp

      END IF

    END DO
    CALL sync

    ! === Basal melt ===
    ! ==================

    ! = Merge all constants into one [m/a K^-1]
    ! =========================================

    q_factor = seawater_density * cp_ocean * gamma_t * sec_per_year / (ice_density * L_fusion)

    ! Loop over all ice shelf and ocean points
    DO vi = mesh%vi1, mesh%vi2

      IF (ice%mask_shelf_a( vi) == 1 .OR. &
          ice%mask_ocean_a( vi) == 1) THEN

        ! = Pressure melting point at the bottom of the ice shelf [K]
        ! ===========================================================

          t_melt = 0.0939_dp - 0.057_dp * BMB%S_ocean_base( vi) - &
                   7.64E-04_dp * ice%Hi_a( vi) * ice_density / seawater_density

        ! = Thermal forcing [K]
        ! =====================

        t_force = BMB%T_ocean_base( vi) - t_melt

        ! = Sub-shelf basal melt [m/a]
        ! ============================

        BMB%BMB_shelf( vi) = q_factor * F_melt( vi) * (-t_force) * ABS(t_force)

      END IF

    END DO
    CALL sync

    ! Clean after yourself
    CALL deallocate_shared( wF_melt)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_Bernales202X

  SUBROUTINE initialise_BMB_model_Bernales202X( mesh, ice, ocean, BMB)
    ! Allocate memory for the data fields of the Bernales et al. (202X) shelf BMB parameterisation.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_Bernales202X'
    INTEGER                                            :: k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Variables
    CALL allocate_shared_dp_1D(  mesh%nV, BMB%T_ocean_base, BMB%wT_ocean_base)
    CALL allocate_shared_dp_1D(  mesh%nV, BMB%S_ocean_base, BMB%wS_ocean_base)
    CALL allocate_shared_int_1D( mesh%nV, BMB%M_ocean_base, BMB%wM_ocean_base)

    ! Initialise ocean temperature and salinity at the base of the ice shelf
    CALL calc_ocean_temperature_at_shelf_base( mesh, ice, ocean, BMB)
    CALL calc_ocean_salinity_at_shelf_base( mesh, ice, ocean, BMB)

    IF (C%do_ocean_inv) THEN
      CALL allocate_shared_dp_1D( mesh%nV, BMB%T_base_ave, BMB%wT_base_ave)
      CALL allocate_shared_dp_2D( mesh%nV, C%ocean_inv_window_size, BMB%T_base_window, BMB%wT_base_window)

      BMB%T_base_ave( mesh%vi1:mesh%vi2) = BMB%T_ocean_base( mesh%vi1:mesh%vi2)

      DO k = 1, C%ocean_inv_window_size
        BMB%T_base_window( mesh%vi1:mesh%vi2, k) = BMB%T_ocean_base( mesh%vi1:mesh%vi2)
      END DO
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=3)

  END SUBROUTINE initialise_BMB_model_Bernales202X

  SUBROUTINE calc_ocean_temperature_at_shelf_base( mesh, ice, ocean, BMB)
    ! Calculate ocean temperature at the base of the shelf by interpolating
    ! the 3-D ocean temperature field in the vertical column

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_ocean_temperature_at_shelf_base'
    INTEGER                                            :: vi
    REAL(dp)                                           :: depth

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise at zero
      BMB%T_ocean_base( vi) = 0._dp

      IF (ice%mask_shelf_a( vi) == 1 .OR. &
          ice%mask_ocean_a( vi) == 1) THEN

        ! Calculate depth for ice shelves or ocean
        depth = MAX( 0.1_dp, ice%Hi_a( vi) * ice_density / seawater_density)

        ! Find ocean temperature at this depth
        CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( vi,:), depth, BMB%T_ocean_base( vi))

      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_ocean_temperature_at_shelf_base

  SUBROUTINE calc_ocean_salinity_at_shelf_base( mesh, ice, ocean, BMB)
    ! Calculate ocean temperature at the base of the shelf by interpolating
    ! the 3-D ocean temperature field in the vertical column

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_ocean_salinity_at_shelf_base'
    INTEGER                                            :: vi
    REAL(dp)                                           :: depth

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise at zero
      BMB%S_ocean_base( vi) = 0._dp

      IF (ice%mask_shelf_a( vi) == 1 .OR. &
          ice%mask_ocean_a( vi) == 1) THEN

        ! Calculate depth for ice shelves or ocean
        depth = MAX( 0.1_dp, ice%Hi_a( vi) * ice_density / seawater_density)

        ! Find ocean temperature at this depth
        CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S_ocean_corr_ext( vi,:), depth, BMB%S_ocean_base( vi))

      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_ocean_salinity_at_shelf_base

  SUBROUTINE calc_distance_to_grounding_line( mesh, ice, vi, dist_gl)
    ! Calculate distance-to-grouding-line of a vertex.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),      INTENT(IN)    :: mesh
    TYPE(type_ice_model), INTENT(IN)    :: ice
    INTEGER,              INTENT(IN)    :: vi
    REAL(dp),             INTENT(OUT)   :: dist_gl

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER       :: routine_name = 'calc_distance_to_grounding_line'
    INTEGER                             :: thetai
    REAL(dp)                            :: theta
    INTEGER,  PARAMETER                 :: ntheta = 16 ! Number of directions we'll look into
    INTEGER,  DIMENSION(:), ALLOCATABLE :: sees_gl_in_dir
    REAL(dp), DIMENSION(:), ALLOCATABLE :: dist_gl_in_dir

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (ice%mask_glf_a( vi) == 1) THEN
      ! Ice shelf next to a grounding line

      dist_gl = 0._dp

    ELSEIF (ice%mask_shelf_a( vi) == 1) THEN
      ! Ice shelf away from grounding line or ice-free ocean

      ! Look in 16 directions
      ALLOCATE( sees_gl_in_dir( ntheta))
      ALLOCATE( dist_gl_in_dir( ntheta))

      DO thetai = 1, ntheta
        theta = (thetai-1) * 2._dp * pi / ntheta
        CALL look_for_grounding_line( mesh, ice, vi, theta, sees_gl_in_dir( thetai), dist_gl_in_dir( thetai))
      END DO

      dist_gl = MINVAL( dist_gl_in_dir)

      DEALLOCATE( sees_gl_in_dir)
      DEALLOCATE( dist_gl_in_dir)

    ELSE
      ! Not ice shelf

      dist_gl = 1.0E20_dp

    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_distance_to_grounding_line

  SUBROUTINE look_for_grounding_line( mesh, ice, vi, theta, sees_gl, dist_gl)
    ! Look outward from vertex vi in direction theta and check if we can "see"
    ! a grounding line without any land or grounded ice in between, return true
    ! and the distance to this groudning line point.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp),                            INTENT(IN)    :: theta
    INTEGER,                             INTENT(OUT)   :: sees_gl
    REAL(dp),                            INTENT(OUT)   :: dist_gl

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'look_for_grounding_line'
    REAL(dp), PARAMETER                                :: max_look_distance = 750000._dp
    REAL(dp), DIMENSION(2)                             :: p, q
    REAL(dp)                                           :: d_min, d, distance_along_line
    INTEGER                                            :: vi_prev, vi_cur, vi_next, ci, vc
    LOGICAL                                            :: Finished

    ! Add routine to path
    CALL init_routine( routine_name)

    sees_gl = 0
    dist_gl = mesh%xmax + mesh%ymax - mesh%xmin - mesh%ymin

    ! End points of trace
    p = mesh%V( vi,:)
    q = [mesh%V( vi,1) + max_look_distance * COS( theta), &
         mesh%V( vi,2) + max_look_distance * SIN( theta)]

    vi_prev = 0
    vi_cur  = vi

    Finished = .FALSE.
    DO WHILE (.NOT. Finished)

      ! Find the next vertex point the trace. Instead of doing a "real"
      ! mesh transect, just take the neighbour vertex closest to q,
      ! which is good enough for this and is MUCH faster.

      vi_next = 0
      d_min   = mesh%xmax + mesh%ymax - mesh%xmin - mesh%ymin
      DO ci = 1, mesh%nC( vi_cur)
        vc = mesh%C( vi_cur,ci)
        p  = mesh%V( vc,:)
        d  = NORM2( q-p)
        IF (d < d_min) THEN
          d_min = d
          vi_next = vc
        END IF
      END DO

      distance_along_line = NORM2( mesh%V( vi_next,:) - mesh%V( vi,:) )

      ! Check if we've finished the trace
      IF (ice%mask_gl_a( vi_next) == 1) THEN
        sees_gl = 1
        dist_gl = distance_along_line
        Finished = .TRUE.
      ELSEIF (ice%mask_shelf_a( vi_next) == 0) THEN
        sees_gl = 0
        Finished = .TRUE.
      ELSEIF (distance_along_line > max_look_distance) THEN
        sees_gl = 0
        Finished = .TRUE.
      ELSEIF (mesh%edge_index( vi_next) > 0) THEN
        sees_gl = 0
        Finished = .TRUE.
      ELSEIF (vi_prev == vi_next) THEN
        sees_gl = 0
        Finished = .TRUE.
      END IF

      ! Move to next vertex
      vi_prev = vi_cur
      vi_cur  = vi_next

    END DO ! WHILE (.NOT. Finished)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE look_for_grounding_line

! ===== Inversion of ice shelf basal melt rates =====
! ===================================================

  SUBROUTINE run_BMB_model_shelf_inversion( mesh, ice, BMB, refgeo)
    ! Invert basal melt using the reference topography

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_shelf_inversion'
    INTEGER                                            :: vi
    REAL(dp)                                           :: h_delta, h_scale

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      h_delta = ice%Hi_a( vi) - refgeo%Hi( vi)

      ! Invert only over shelf or ocean points, identified from
      ! the reference topography (on the mesh), the model mask,
      ! or the sub-grid grounded area fraction
      IF ( is_floating( refgeo%Hi( vi), refgeo%Hb( vi), 0._dp) .OR. &
           ice%mask_shelf_a( vi) == 1 .OR. &
           ice%f_grndx_a( vi) < 1.0_dp) THEN

        IF (refgeo%Hi( vi) > 0._dp) THEN
          h_scale = 1.0_dp/C%BMB_inv_scale_shelf
        ELSE
          h_scale = 1.0_dp/C%BMB_inv_scale_ocean
        END IF

        h_delta = MAX(-1.5_dp, MIN(1.5_dp, h_delta * h_scale))

        ! Further adjust only where the previous value is not improving the result
        IF ( (h_delta > 0._dp .AND. ice%dHi_dt_a( vi) >= 0._dp) .OR. &
             (h_delta < 0._dp .AND. ice%dHi_dt_a( vi) <= 0._dp) ) THEN

          BMB%BMB_shelf( vi) = BMB%BMB_shelf( vi) - 1.72476_dp * TAN(h_delta)
                                                  ! The hardcoded jorjonian constant yeah baby.

        END IF ! else BMB_shelf does not change from previous time step

      END IF ! else the reference is grounded ice sheet, so leave it alone

    END DO
    CALL sync

    ! Limit basal melt
    DO vi = mesh%vi1, mesh%vi2
      BMB%BMB_shelf( vi) = MAX( BMB%BMB_shelf( vi), C%BMB_min)
      BMB%BMB_shelf( vi) = MIN( BMB%BMB_shelf( vi), C%BMB_max)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_shelf_inversion

  SUBROUTINE initialise_BMB_model_inversion( mesh, BMB, restart)
    ! Invert basal melt using the reference topography

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    TYPE(type_restart_data),             INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_inversion'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%is_restart .AND. C%BMB_inv_use_restart_field) THEN

      IF (par%master) WRITE (0,*) '   Initialising ice shelf basal mass balance using data read from restart file...'

      ! Assign field from restart file
      BMB%BMB_shelf( mesh%vi1:mesh%vi2) = restart%BMB_shelf( mesh%vi1:mesh%vi2)

    ELSE
      ! No need to do anything
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model_inversion

! ===== Inversion of ocean temperatures =====
! ===========================================

  SUBROUTINE ocean_temperature_inversion( mesh, grid, ice, BMB, refgeo, time)
    ! Invert basal ocean temps using the reference topography

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),               INTENT(IN)    :: mesh
    TYPE(type_grid),               INTENT(IN)    :: grid
    TYPE(type_ice_model),          INTENT(IN)    :: ice
    TYPE(type_BMB_model),          INTENT(INOUT) :: BMB
    TYPE(type_reference_geometry), INTENT(IN)    :: refgeo
    REAL(dp),                      INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'ocean_temperature_inversion'
    INTEGER                                      :: vi
    REAL(dp)                                     :: h_delta, h_scale, t_scale, a_scale, m_scale
    REAL(dp)                                     :: t_melt
    INTEGER,  DIMENSION(:), POINTER              ::  mask,  mask_filled
    INTEGER                                      :: wmask, wmask_filled
    REAL(dp), DIMENSION(:), POINTER              ::  ocean_smoothed
    INTEGER                                      :: wocean_smoothed

    ! Initialisation
    ! ==============

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (time < C%ocean_inv_t_start) THEN
      ! Nothing to do for now. Just return.
      CALL finalise_routine( routine_name)
      RETURN
    ELSEIF (time > C%ocean_inv_t_end) THEN
      ! Inversion is done. Return.
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Allocate masks for extrapolation
    CALL allocate_shared_int_1D( mesh%nV, mask,        wmask       )
    CALL allocate_shared_int_1D( mesh%nV, mask_filled, wmask_filled)

    ! Allocate smoothed bed roughness field
    CALL allocate_shared_dp_1D( mesh%nV, ocean_smoothed, wocean_smoothed)

    ! Adjustment magnitude
    ! ====================

    ! Default values
    a_scale = C%ocean_inv_hi_scale
    m_scale = .000_dp
    t_scale = 0._dp

    ! Time scale
    IF (C%ocean_inv_t_start < C%ocean_inv_t_end) THEN
      ! Compute how much time has passed since start of inversion
      t_scale = (time - C%ocean_inv_t_start) / (C%ocean_inv_t_end - C%ocean_inv_t_start)
      ! Limit t_scale to [0 1]
      t_scale = MAX( 0._dp, MIN( t_scale, 1._dp))
    END IF

    ! Magnitude decay
    IF (C%do_ocean_inv_decay) THEN
      ! Reduce adjustment amount as time goes on
      a_scale = C%ocean_inv_hi_scale * (1._dp - t_scale) + 0.000 * t_scale
    END IF

    ! Do the inversion
    ! ================

    ! milg

    DO vi = mesh%vi1, mesh%vi2

      h_delta = ice%Hi_a( vi) - refgeo%Hi( vi)

      ! Invert only over non-calving-front shelf vertices
      IF ( ice%mask_shelf_a( vi) == 1 .AND. &
           ( ice%mask_cf_a(  vi) == 0 .OR. &
             ice%mask_glf_a( vi) == 1 ) ) THEN

        ! Add this vertex to mask of inverted ocean temperatures
        BMB%M_ocean_base( vi) = 1

        ! Use this vertex during extrapolation
        mask( vi) = 2

        IF ( h_delta > 0._dp .AND. ice%dHi_dt_a( vi) >= .0_dp ) THEN

          IF (t_scale < .9_dp) THEN
            BMB%T_ocean_base( vi) = BMB%T_ocean_base( vi) + a_scale * (1._dp - EXP(-ABS(h_delta*ice%dHi_dt_a( vi))))
          ELSE
            BMB%T_ocean_base( vi) = BMB%T_ocean_base( vi) + C%ocean_inv_hi_scale * (1._dp - EXP(-ABS(ice%dHi_dt_a( vi))))
          END IF

        ELSEIF ( h_delta < 0._dp .AND. ice%dHi_dt_a( vi) <= .0_dp ) THEN

          IF (t_scale < .9_dp) THEN
            BMB%T_ocean_base( vi) = BMB%T_ocean_base( vi) - a_scale * (1._dp - EXP(-ABS(h_delta*ice%dHi_dt_a( vi))))
          ELSE
            BMB%T_ocean_base( vi) = BMB%T_ocean_base( vi) - C%ocean_inv_hi_scale/20._dp * (1._dp - EXP(-ABS(ice%dHi_dt_a( vi))))
          END IF

        ELSEIF ( h_delta > 0._dp .AND. ice%dHi_dt_a( vi) < .0_dp ) THEN

          IF (t_scale < .9_dp) THEN
            BMB%T_ocean_base( vi) = BMB%T_ocean_base( vi) - m_scale * (1._dp - EXP(-ABS(ice%dHi_dt_a( vi))))
          END IF

        ELSEIF ( h_delta < 0._dp .AND. ice%dHi_dt_a( vi) > .0_dp ) THEN

          IF (t_scale < .9_dp) THEN
            BMB%T_ocean_base( vi) = BMB%T_ocean_base( vi) + m_scale * (1._dp - EXP(-ABS(ice%dHi_dt_a( vi))))
          END IF

        END IF

      ELSE

        ! Remove this vertex from mask of inverted ocean temperatures
        BMB%M_ocean_base( vi) = 0

        ! Not ice shelf: mark it for extrapolation
        mask( vi) = 1

      END IF

    END DO
    CALL sync

    ! Limit basal temperatures
    DO vi = mesh%vi1, mesh%vi2

      IF (ice%mask_shelf_a( vi) == 1) THEN

        t_melt = 0.0939_dp - 0.057_dp * BMB%S_ocean_base( vi) - &
                 7.64E-04_dp * ice%Hi_a( vi) * ice_density / seawater_density

        BMB%T_ocean_base( vi) = MAX( BMB%T_ocean_base( vi), t_melt)
        BMB%T_ocean_base( vi) = MIN( BMB%T_ocean_base( vi),  5._dp)

      END IF

    END DO
    CALL sync

    ! Extrapolate the resulting field
    ! ===============================

    ! Perform the extrapolation
    IF (par%master) THEN
      CALL extrapolate_Gaussian_floodfill_mesh( mesh, mask, BMB%T_ocean_base, 40000._dp, mask_filled)
    END IF
    CALL sync

    ! Smoothing
    ! =========

    IF (C%do_ocean_inv_smooth) THEN
      ! Smooth the resulting field

      ! Store the inverted parameters in a local variable
      ocean_smoothed( mesh%vi1:mesh%vi2) = BMB%T_ocean_base( mesh%vi1:mesh%vi2)
      CALL sync

      ! Smooth the local variable
      CALL smooth_Gaussian_2D( mesh, grid, ocean_smoothed, C%ocean_inv_smooth_r)

      ! Combined the smoothed and raw inverted parameter through a weighed average
      DO vi = mesh%vi1, mesh%vi2
          BMB%T_ocean_base( vi) = (1._dp - C%ocean_inv_smooth_w) * BMB%T_ocean_base( vi) + C%ocean_inv_smooth_w * ocean_smoothed( vi)
      END DO
      CALL sync

    END IF ! (C%do_ocean_inv_smooth)
    CALL sync

    ! Running T_ocean_base average
    ! ============================

    ! Update the running window: drop oldest record and push the rest to the back
    BMB%T_base_window( mesh%vi1:mesh%vi2,2:C%ocean_inv_window_size) = BMB%T_base_window( mesh%vi1:mesh%vi2,1:C%ocean_inv_window_size-1)

    ! Update the running window: add new record to beginning of window
    BMB%T_base_window( mesh%vi1:mesh%vi2,1) = BMB%T_ocean_base( mesh%vi1:mesh%vi2)

    ! Compute running average
    BMB%T_base_ave( mesh%vi1:mesh%vi2) = SUM(BMB%T_base_window( mesh%vi1:mesh%vi2,:),2) / REAL(C%ocean_inv_window_size,dp)

    ! Finalisation
    ! ============

    CALL deallocate_shared( wmask)
    CALL deallocate_shared( wmask_filled)
    CALL deallocate_shared( wocean_smoothed)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ocean_temperature_inversion

! ===== Some generally useful tools =====
! =======================================

  SUBROUTINE calculate_sub_angle_dist_open( mesh, ice, vi, sub_angle, dist_open)
    ! Calculate the "subtended angle" (i.e. how many degrees of its horizon have
    ! a view of open ocean) and distance-to-open-ocean of a vertex.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp),                            INTENT(OUT)   :: sub_angle
    REAL(dp),                            INTENT(OUT)   :: dist_open

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calculate_sub_angle_dist_open'
    INTEGER                                            :: thetai
    REAL(dp)                                           :: theta
    INTEGER,  PARAMETER                                :: ntheta            = 16         ! Number of directions we'll look into
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: sees_ocean_in_dir
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: dist_ocean_in_dir

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Only calculate the angle for shelf vertices
    IF (.NOT. ice%mask_shelf_a( vi) == 1) THEN
      IF (ice%mask_ocean_a( vi) == 1) THEN
        sub_angle = 360._dp
        dist_open = 0._dp
        CALL finalise_routine( routine_name)
        RETURN
      ELSE
        sub_angle = 0._dp
        dist_open = mesh%xmax + mesh%ymax - mesh%xmin - mesh%ymin
        CALL finalise_routine( routine_name)
        RETURN
      END IF
    END IF

    ! Look in 16 directions
    ALLOCATE( sees_ocean_in_dir( ntheta))
    ALLOCATE( dist_ocean_in_dir( ntheta))

    DO thetai = 1, ntheta
      theta = (thetai-1) * 2._dp * pi / ntheta
      CALL look_for_ocean( mesh, ice, vi, theta, sees_ocean_in_dir( thetai), dist_ocean_in_dir( thetai))
    END DO

    sub_angle = SUM(    sees_ocean_in_dir) * 360._dp / ntheta
    dist_open = MINVAL( dist_ocean_in_dir)

    DEALLOCATE( sees_ocean_in_dir)
    DEALLOCATE( dist_ocean_in_dir)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_sub_angle_dist_open

  SUBROUTINE look_for_ocean( mesh, ice, vi, theta, sees_ocean, dist_ocean)
    ! Look outward from vertex vi in direction theta and check if we can "see"
    ! open ocean without any land or grounded ice in between, return true (and the distance to this ocean).

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp),                            INTENT(IN)    :: theta
    INTEGER,                             INTENT(OUT)   :: sees_ocean
    REAL(dp),                            INTENT(OUT)   :: dist_ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'look_for_ocean'
    REAL(dp), PARAMETER                                :: max_look_distance = 750000._dp ! Maximum distance to look before we decide no ocean will ever be found.
    REAL(dp), DIMENSION(2)                             :: p, q
    REAL(dp)                                           :: d_min, d, distance_along_line
    INTEGER                                            :: vi_prev, vi_cur, vi_next, ci, vc
    LOGICAL                                            :: Finished

    ! Add routine to path
    CALL init_routine( routine_name)

    sees_ocean = 0
    dist_ocean = mesh%xmax + mesh%ymax - mesh%xmin - mesh%ymin

    ! End points of trace
    p = mesh%V( vi,:)
    q = [mesh%V( vi,1) + max_look_distance * COS( theta), &
         mesh%V( vi,2) + max_look_distance * SIN( theta)]

    vi_prev = 0
    vi_cur  = vi

    Finished = .FALSE.
    DO WHILE (.NOT. Finished)

      ! Find the next vertex point the trace. Instead of doing a "real"
      ! mesh transect, just take the neighbour vertex closest to q,
      ! which is good enough for this and is MUCH faster.

      vi_next = 0
      d_min   = mesh%xmax + mesh%ymax - mesh%xmin - mesh%ymin
      DO ci = 1, mesh%nC( vi_cur)
        vc = mesh%C( vi_cur,ci)
        p  = mesh%V( vc,:)
        d  = NORM2( q-p)
        IF (d < d_min) THEN
          d_min = d
          vi_next = vc
        END IF
      END DO

      distance_along_line = NORM2( mesh%V( vi_next,:) - mesh%V( vi,:) )

      ! Check if we've finished the trace
      IF (ice%mask_land_a( vi_next) == 1) THEN
        sees_ocean = 0
        Finished = .TRUE.
      ELSEIF (ice%mask_ocean_a( vi_next) == 1 .AND. ice%mask_shelf_a( vi_next) == 0) THEN
        sees_ocean = 1
        dist_ocean = distance_along_line
        Finished = .TRUE.
      ELSEIF (distance_along_line > max_look_distance) THEN
        sees_ocean = 0
        Finished = .TRUE.
      ELSEIF (mesh%edge_index( vi_next) > 0) THEN
        sees_ocean = 0
        Finished = .TRUE.
      ELSEIF (vi_prev == vi_next) THEN
        sees_ocean = 0
        Finished = .TRUE.
      END IF

      ! Move to next vertex
      vi_prev = vi_cur
      vi_cur  = vi_next

    END DO ! DO WHILE (.NOT. Finished)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE look_for_ocean

  SUBROUTINE extrapolate_melt_from_FCMP_to_PMP( mesh, ice, BMB)
    ! All the BMB parameterisations are implicitly run using the FCMP sub-grid scheme
    ! (i.e. they are only applied to grid cells whose centre is floating).
    ! Calculating melt rates for partially-floating-but-grounded-at-the-centre cells (which
    ! is needed in the PMP scheme) is not straightforward and extrapolating values into
    ! the grid cells ala IMAU-ICE requires some light puzzle solving. For now, just take
    ! a weighed average of the surrounding FCMP vertices to get something simple going.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extrapolate_melt_from_FCMP_to_PMP'
    INTEGER                                            :: vi, vc, ci
    INTEGER,  DIMENSION(:    ), POINTER                :: mask_FCMP, mask_PMP
    REAL(dp), DIMENSION(:    ), POINTER                :: BMB_shelf_extra
    INTEGER                                            :: wmask_FCMP, wmask_PMP, wBMB_shelf_extra
    REAL(dp)                                           :: BMB_av, sum_dist

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_1D( mesh%nV, mask_FCMP,       wmask_FCMP      )
    CALL allocate_shared_int_1D( mesh%nV, mask_PMP,        wmask_PMP       )
    CALL allocate_shared_dp_1D(  mesh%nV, BMB_shelf_extra, wBMB_shelf_extra)

    ! Define the two masks
    DO vi = mesh%vi1, mesh%vi2
      mask_FCMP( vi) = 0
      mask_PMP(  vi) = 0
      IF (ice%mask_shelf_a( vi) == 1) mask_FCMP( vi) = 1
      IF (ice%f_grndx_a( vi) < 1._dp .AND. ice%mask_ice_a( vi) == 1) mask_PMP( vi) = 1
    END DO
    CALL sync

    ! Extrapolate melt from the FCMP to the PMP mask by taking the average over all FCMP neighbours
    DO vi = mesh%vi1, mesh%vi2

      ! Initialise
      BMB_shelf_extra( vi) = 0._dp

      IF (mask_PMP( vi) == 1) THEN
        IF (mask_FCMP( vi) == 1) THEN
          ! Simply copy the FCMP value

          BMB_shelf_extra( vi) = BMB%BMB_shelf( vi)

        ELSE
          ! Calculate distance-based weighed-average melt rate over all FCMP neighbours

          sum_dist = 0._dp
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            IF (mask_FCMP( vc) == 1) THEN
              sum_dist = sum_dist + NORM2( mesh%V( vi,:) - mesh%V( vc,:))
            END IF
          END DO

          BMB_av = 0._dp
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            IF (mask_FCMP( vc) == 1 .AND. sum_dist > 0._dp) THEN
              BMB_av = BMB_av + BMB%BMB_shelf( vc) * NORM2( mesh%V( vi,:) - mesh%V( vc,:)) / sum_dist
            END IF
          END DO

          BMB_shelf_extra( vi) = BMB_av

        END IF ! (mask_FCMP( vi) == 1)
      END IF ! (mask_PMP( vi) == 1)

    END DO
    CALL sync

    ! Copy results back to original array
    BMB%BMB_shelf( mesh%vi1:mesh%vi2) = BMB_shelf_extra( mesh%vi1:mesh%vi2)

    ! Clean up after yourself
    CALL deallocate_shared( wmask_FCMP      )
    CALL deallocate_shared( wmask_PMP       )
    CALL deallocate_shared( wBMB_shelf_extra)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extrapolate_melt_from_FCMP_to_PMP

!===== Remapping after mesh update =====
!=======================================

  SUBROUTINE remap_BMB_model( mesh_old, mesh_new, map, BMB)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_BMB_model'
    INTEGER                                            :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV
    int_dummy = map%int_dummy

    ! Reallocate rather than remap; after a mesh update we'll immediately run the BMB model anyway.
    CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%BMB,       BMB%wBMB      )
    CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%BMB_sheet, BMB%wBMB_sheet)

    ! Exception for iterative inversion of ice shelf basal melt rates, to avoid resetting it.
    IF (C%choice_BMB_shelf_model == 'inversion') THEN
      CALL remap_field_dp_2D( mesh_old, mesh_new, map, BMB%BMB_shelf, BMB%wBMB_shelf, 'cons_2nd_order')
    ELSE
      CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%BMB_shelf, BMB%wBMB_shelf)
    END IF

    ! Sheet
    ! =====

    IF     (C%choice_BMB_sheet_model == 'uniform') THEN
      ! Nothing else needs to be done
    ELSE
      CALL crash('unknown choice_BMB_sheet_model "' // TRIM(C%choice_BMB_sheet_model) // '"!')
    END IF

    ! Shelf
    ! =====

    IF     (C%choice_BMB_shelf_model == 'uniform' .OR. &
            C%choice_BMB_shelf_model == 'idealised') THEN

      ! Nothing else needs to be done

    ELSEIF (C%choice_BMB_shelf_model == 'ANICE_legacy') THEN

      CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%sub_angle, BMB%wsub_angle)
      CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%dist_open, BMB%wdist_open)

    ELSEIF (C%choice_BMB_shelf_model == 'Favier2019_lin' .OR. &
            C%choice_BMB_shelf_model == 'Favier2019_quad' .OR. &
            C%choice_BMB_shelf_model == 'Favier2019_Mplus') THEN

      CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%T_ocean_base,        BMB%wT_ocean_base       )
      CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%T_ocean_freeze_base, BMB%wT_ocean_freeze_base)

    ELSEIF (C%choice_BMB_shelf_model == 'Lazeroms2018_plume') THEN

      IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not yet implemented for a triangular mesh!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%T_ocean_base,           BMB%wT_ocean_base          )
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%eff_plume_source_depth, BMB%weff_plume_source_depth)
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%eff_basal_slope,        BMB%weff_basal_slope       )

    ELSEIF (C%choice_BMB_shelf_model == 'PICO') THEN

      IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not yet implemented for a triangular mesh!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_d_GL, BMB%wPICO_d_GL)
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_d_IF, BMB%wPICO_d_IF)
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_r,    BMB%wPICO_r   )
      ! CALL reallocate_shared_int_1D( mesh_new%nV, BMB%PICO_k,    BMB%wPICO_k   )

      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_T,    BMB%wPICO_T )
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_S,    BMB%wPICO_S )
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_p,    BMB%wPICO_p )
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_m,    BMB%wPICO_m )

    ELSEIF (C%choice_BMB_shelf_model == 'PICOP') THEN

      IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not yet implemented for a triangular mesh!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%T_ocean_base,           BMB%wT_ocean_base          )
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%eff_plume_source_depth, BMB%weff_plume_source_depth)
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%eff_basal_slope,        BMB%weff_basal_slope       )

      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_d_GL,              BMB%wPICO_d_GL        )
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_d_IF,              BMB%wPICO_d_IF        )
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_r,                 BMB%wPICO_r           )
      ! CALL reallocate_shared_int_1D( mesh_new%nV, BMB%PICO_k,                 BMB%wPICO_k           )

      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_T,                 BMB%wPICO_T )
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_S,                 BMB%wPICO_S )
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_p,                 BMB%wPICO_p )
      ! CALL reallocate_shared_dp_1D(  mesh_new%nV, BMB%PICO_m,                 BMB%wPICO_m )

    ELSEIF (C%choice_BMB_shelf_model == 'inversion') THEN

      ! Nothing else needs to be done for now. Main stuff was done at the start of this routine.

    ELSEIF (C%choice_BMB_shelf_model == 'Bernales202X') THEN

      ! Exception for iterative inversion of ocean temperatures, to avoid resetting it.
      IF (C%do_ocean_inv) THEN
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, BMB%T_ocean_base,  BMB%wT_ocean_base,  'cons_2nd_order')
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, BMB%T_base_ave,    BMB%wT_base_ave,    'nearest_neighbour')
        CALL remap_field_dp_3D( mesh_old, mesh_new, map, BMB%T_base_window, BMB%wT_base_window, 'nearest_neighbour')
      ELSE
        CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%T_ocean_base, BMB%wT_ocean_base)
      END IF

      CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%S_ocean_base, BMB%wS_ocean_base)

    ELSE
      CALL crash('unknown choice_BMB_shelf_model "' // TRIM(C%choice_BMB_shelf_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_BMB_model

!===== OLD routines ======
!=========================

  ! The old routines. Useful as templates for new stuff.

    ! ! Run the SMB model on the region mesh
    ! SUBROUTINE run_BMB_model_OLD( mesh, ice, climate, BMB, region_name)
    !   ! Calculate mean ocean temperature (saved in "climate") and basal mass balance

    !   IMPLICIT NONE

    !   ! In/output variables
    !   TYPE(type_mesh),                      INTENT(IN)    :: mesh
    !   TYPE(type_ice_model),                 INTENT(IN)    :: ice
    !   TYPE(type_climate_snapshot_regional), INTENT(INOUT) :: climate
    !   TYPE(type_BMB_model),                 INTENT(INOUT) :: BMB
    !   CHARACTER(LEN=3),                     INTENT(IN)    :: region_name

    !   ! Local variables
    !   CHARACTER(LEN=64), PARAMETER                        :: routine_name = 'run_BMB_model'
    !   INTEGER                                             :: n1, n2
    !   INTEGER                                             :: vi
    !   REAL(dp)                                            :: BMB_shelf                             ! Sub-shelf melt rate for non-exposed shelf  [m/year]
    !   REAL(dp)                                            :: BMB_shelf_exposed                     ! Sub-shelf melt rate for exposed shelf      [m/year]
    !   REAL(dp)                                            :: BMB_deepocean                         ! Sub-shelf melt rate for deep-ocean areas   [m/year]
    !   REAL(dp)                                            :: w_ins, w_PD, w_warm, w_cold, w_deep, w_expo, weight
    !   REAL(dp)                                            :: T_freeze                              ! Freezing temperature at the base of the shelf (Celcius)
    !   REAL(dp)                                            :: water_depth
    !   REAL(dp), PARAMETER                                 :: cp0        = 3974._dp                 ! specific heat capacity of the ocean mixed layer (J kg-1 K-1)
    !   REAL(dp), PARAMETER                                 :: gamma_T    = 1.0E-04_dp               ! Thermal exchange velocity (m s-1)

    !   n1 = par%mem%n

    !   ! Exceptions for benchmark experiments
    !   IF (C%do_benchmark_experiment) THEN
    !     IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
    !         C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
    !         C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
    !         C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
    !         C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
    !         C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
    !         C%choice_benchmark_experiment == 'Halfar' .OR. &
    !         C%choice_benchmark_experiment == 'Bueler' .OR. &
    !         C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
    !         C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
    !         C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
    !         C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
    !         C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
    !         C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
    !         C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
    !         C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
    !       BMB%BMB( mesh%vi1:mesh%vi2) = 0._dp
    !       CALL sync
    !       RETURN
    !     ELSE
    !       IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_BMB_model!'
    !       CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    !     END IF
    !   END IF ! IF (C%do_benchmark_experiment) THEN

    !   ! Initialise everything at zero
    !   BMB%BMB(       mesh%vi1:mesh%vi2) = 0._dp
    !   BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
    !   BMB%BMB_shelf( mesh%vi1:mesh%vi2) = 0._dp
    !   BMB%sub_angle( mesh%vi1:mesh%vi2) = 360._dp
    !   BMB%dist_open( mesh%vi1:mesh%vi2) = 0._dp
    !   w_ins                             = 0._dp
    !   weight                            = 0._dp
    !   w_PD                              = 0._dp
    !   w_warm                            = 0._dp
    !   w_cold                            = 0._dp
    !   w_deep                            = 0._dp
    !   w_expo                            = 0._dp
    !   BMB_shelf                         = 0._dp
    !   BMB_shelf_exposed                 = 0._dp
    !   BMB_deepocean                     = 0._dp

    !   ! Find the "subtended angle" and distance-to-open-ocean of all shelf pixels
    !   DO vi = mesh%vi1, mesh%vi2
    !     CALL calculate_sub_angle_dist_open( mesh, ice, vi, BMB%sub_angle( vi), BMB%dist_open( vi))
    !   END DO
    !   CALL sync

    !   ! Find the weight from insolation
    !   IF (region_name == 'NAM' .OR. region_name == 'EAS' .OR. region_name == 'GRL') THEN
    !     w_ins = MAX(0._dp, (climate%Q_TOA_jun_65N - 462.29_dp) / 40._dp)
    !   ELSEIF (region_name == 'ANT') THEN
    !     w_ins = MAX(0._dp, (climate%Q_TOA_jan_80S - 532.19_dp) / 40._dp)
    !   END IF

    !   ! Determine mean ocean temperature and basal melt rates for deep ocean and exposed shelves
    !   ! ========================================================================================

    !   IF (C%choice_ocean_temperature_model == 'fixed') THEN
    !     ! Use present-day values

    !     climate%T_ocean_mean = BMB%T_ocean_mean_PD
    !     BMB_deepocean        = BMB%BMB_deepocean_PD
    !     BMB_shelf_exposed    = BMB%BMB_shelf_exposed_PD

    !   ELSEIF (C%choice_ocean_temperature_model == 'scaled') THEN
    !     ! Scale between config values of mean ocean temperature and basal melt rates for PD, cold, and warm climates.

    !     ! Determine weight for scaling between different ocean temperatures
    !     IF (C%choice_forcing_method == 'CO2_direct') THEN

    !       ! Use the prescribed CO2 record as a glacial index
    !       IF (forcing%CO2_obs > 280._dp) THEN
    !         ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 400 ppmv
    !         weight = 2._dp - MAX( 0._dp,   MIN(1.25_dp, ((400._dp - forcing%CO2_obs) / (400._dp - 280._dp) + (3.00_dp - forcing%d18O_obs) / (3.00_dp - 3.23_dp)) / 2._dp )) + w_ins
    !       ELSE
    !         ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
    !         weight = 1._dp - MAX(-0.25_dp, MIN(1._dp,   ((280._dp - forcing%CO2_obs) / (280._dp - 190._dp) + (3.23_dp - forcing%d18O_obs) / (3.23_dp - 4.95_dp)) / 2._dp )) + w_ins
    !       END IF

    !     ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN

    !       ! Use modelled CO2 as a glacial index
    !       IF (forcing%CO2_obs > 280._dp) THEN
    !         ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 400 ppmv
    !         weight = 2._dp - MAX( 0._dp,   MIN(1.25_dp, ((400._dp - forcing%CO2_mod) / (400._dp - 280._dp) + (3.00_dp - forcing%d18O_obs) / (3.00_dp - 3.23_dp)) / 2._dp )) + w_ins
    !       ELSE
    !         ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
    !         weight = 1._dp - MAX(-0.25_dp, MIN(1._dp,   ((280._dp - forcing%CO2_mod) / (280._dp - 190._dp) + (3.23_dp - forcing%d18O_obs) / (3.23_dp - 4.95_dp)) / 2._dp )) + w_ins
    !       END IF

    !     ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN

    !       ! Use modelled global mean annual temperature change as a glacial index
    !       weight = MAX(0._dp, MIN(2._dp, 1._dp + forcing%dT_glob/12._dp + w_ins))

    !     ELSE ! IF (C%choice_forcing_method == 'CO2_direct') THEN
    !       WRITE(0,*) '  ERROR: forcing method "', TRIM(C%choice_forcing_method), '" not implemented in run_BMB_model!'
    !       CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    !     END IF ! IF (C%choice_forcing_method == 'CO2_direct') THEN

    !     IF (weight < 1._dp) THEN
    !       w_PD   = weight
    !       w_cold = 1._dp - w_PD
    !       w_warm = 0._dp
    !     ELSE
    !       w_PD   = 2._dp - weight
    !       w_warm = 1._dp - w_PD
    !       w_cold = 0._dp
    !     END IF

    !     climate%T_ocean_mean = w_PD * BMB%T_ocean_mean_PD      + w_warm * BMB%T_ocean_mean_warm      + w_cold * BMB%T_ocean_mean_cold
    !     BMB_deepocean        = w_PD * BMB%BMB_deepocean_PD     + w_warm * BMB%BMB_deepocean_warm     + w_cold * BMB%BMB_deepocean_cold
    !     BMB_shelf_exposed    = w_PD * BMB%BMB_shelf_exposed_PD + w_warm * BMB%BMB_shelf_exposed_warm + w_cold * BMB%BMB_shelf_exposed_cold

    !   ELSE ! IF (C%choice_ocean_temperature_model == 'fixed') THEN
    !     WRITE(0,*) '  ERROR: choice_ocean_temperature_model "', TRIM(C%choice_ocean_temperature_model), '" not implemented in run_BMB_model!'
    !     CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    !   END IF ! IF (C%choice_ocean_temperature_model == 'fixed') THEN

    !   ! Use the (interpolated, spatially uniform) ocean temperature and the subtended angle + distance-to-open-ocean
    !   ! to calculate sub-shelf melt rates using the parametrisation from Martin et al., 2011
    !   ! ====================================================================================

    !   DO vi = mesh%vi1, mesh%vi2

    !     IF (ice%mask_shelf_a( vi) == 1) THEN
    !       ! Sub-shelf melt

    !       ! Freezing temperature at the bottom of the ice shelves, scaling with depth below water level
    !       T_freeze = 0.0939_dp - 0.057_dp * 35._dp - 7.64E-04_dp * ice%Hi_a( vi) * ice_density / seawater_density

    !       ! Sub-shelf melt rate for non-exposed shelves (Martin, TC, 2011) - melt values, when T_ocean > T_freeze.
    !       BMB_shelf   = seawater_density * cp0 * sec_per_year * gamma_T * BMB%subshelf_melt_factor * &
    !                  (climate%T_ocean_mean - T_freeze) / (L_fusion * ice_density)

    !     ELSE
    !       BMB_shelf = 0._dp
    !     END IF

    !     IF (ice%mask_shelf_a( vi) == 1 .OR. ice%mask_ocean_a( vi) == 1) THEN

    !       water_depth = ice%SL_a( vi) - ice%Hb_a( vi)
    !       w_deep = MAX(0._dp, MIN(1._dp, (water_depth - BMB%deep_ocean_threshold_depth) / 200._dp))
    !       w_expo = MAX(0._dp, MIN(1._dp, (BMB%sub_angle( vi) - 80._dp)/30._dp)) * EXP(-BMB%dist_open( vi) / 100000._dp)

    !       BMB%BMB_shelf( vi) = (w_deep * BMB_deepocean) + (1._dp - w_deep) * (w_expo * BMB_shelf_exposed + (1._dp - w_expo) * BMB_shelf)

    !     ELSE
    !       BMB%BMB_shelf( vi) = 0._dp
    !     END IF

    !   END DO
    !   CALL sync

    !   ! Add sheet and shelf melt together
    !   BMB%BMB( mesh%vi1:mesh%vi2) = BMB%BMB_sheet( mesh%vi1:mesh%vi2) + BMB%BMB_shelf( mesh%vi1:mesh%vi2)
    !   CALL sync

    !   n2 = par%mem%n
    !   !CALL write_to_memory_log( routine_name, n1, n2)

    !   IF (par%master) WRITE(0,*) 'This subroutine (run_BMB_model) is empty. Feel free to fill it before running the model :)'
    !   CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    ! END SUBROUTINE run_BMB_model_OLD

    ! SUBROUTINE remap_BMB_model_OLD( mesh_old, mesh_new, map, BMB)
    !   ! Remap or reallocate all the data fields

    !   ! In/output variables:
    !   TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    !   TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    !   TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    !   TYPE(type_BMB_model),                INTENT(INOUT) :: BMB

    !   ! Local variables:
    !   INTEGER                                            :: int_dummy

    !   ! To prevent compiler warnings for unused variables
    !   int_dummy = mesh_old%nV
    !   int_dummy = mesh_new%nV
    !   int_dummy = map%int_dummy

    !   ! Reallocate rather than remap; after a mesh update we'll immediately run the BMB model anyway
    !   CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%BMB,       BMB%wBMB      )
    !   CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%BMB_sheet, BMB%wBMB_sheet)
    !   CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%BMB_shelf, BMB%wBMB_shelf)
    !   CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%sub_angle, BMB%wsub_angle)
    !   CALL reallocate_shared_dp_1D( mesh_new%nV, BMB%dist_open, BMB%wdist_open)

    ! END SUBROUTINE remap_BMB_model_OLD

END MODULE BMB_module
