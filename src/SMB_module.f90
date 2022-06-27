MODULE SMB_module

  ! All the routines for calculating the surface mass balance from the climate forcing

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
                                             transpose_dp_2D, transpose_dp_3D
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE data_types_module,               ONLY: type_mesh, type_ice_model, &
                                             type_SMB_model, type_remapping_mesh_mesh, &
                                             type_climate_matrix_regional, &
                                             type_climate_snapshot_regional, type_direct_SMB_forcing_regional, &
                                             type_restart_data, type_grid, type_reference_geometry
  USE forcing_module,                  ONLY: forcing
  USE mesh_mapping_module,             ONLY: remap_field_dp_2D, remap_field_dp_3D, &
                                             calc_remapping_operator_grid2mesh, map_grid2mesh_2D, map_grid2mesh_3D, &
                                             deallocate_remapping_operators_grid2mesh

  IMPLICIT NONE

  REAL(dp), PARAMETER :: albedo_water        = 0.1_dp
  REAL(dp), PARAMETER :: albedo_soil         = 0.2_dp
  REAL(dp), PARAMETER :: albedo_ice          = 0.5_dp
  REAL(dp), PARAMETER :: albedo_snow         = 0.85_dp
  REAL(dp), PARAMETER :: initial_snow_depth  = 0.1_dp

CONTAINS

! == The main routines that should be called from the main ice model/program
! ==========================================================================

  SUBROUTINE run_SMB_model( mesh, ice, climate_matrix, time, SMB, mask_noice)
    ! Run the selected regional SMB model

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_matrix_regional),   INTENT(IN)    :: climate_matrix
    REAL(dp),                             INTENT(IN)    :: time
    TYPE(type_SMB_model),                 INTENT(INOUT) :: SMB
    INTEGER,  DIMENSION(:    ),           INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'run_SMB_model'
    INTEGER                                             :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_SMB_model == 'uniform') THEN
      ! Apply a simple uniform SMB

      DO vi = mesh%vi1, mesh%vi2
        IF (mask_noice( vi) == 0) THEN
          SMB%SMB_year( vi) = C%SMB_uniform
        ELSE
          SMB%SMB_year( vi) = 0._dp
        END IF
      END DO
      CALL sync

    ELSEIF (C%choice_SMB_model == 'idealised') THEN
      ! Apply an idealised SMB parameterisation

      CALL run_SMB_model_idealised( mesh, SMB, time, mask_noice)

    ELSEIF (C%choice_SMB_model == 'IMAU-ITM') THEN
      ! Run the IMAU-ITM SMB model

      CALL run_SMB_model_IMAUITM( mesh, ice, climate_matrix%applied, SMB, mask_noice)

    ELSEIF (C%choice_SMB_model == 'IMAU-ITM_wrongrefreezing') THEN
      ! Run the IMAU-ITM SMB model with the old wrong refreezing parameterisation from ANICE

      CALL run_SMB_model_IMAUITM_wrongrefreezing( mesh, ice, climate_matrix%applied, SMB, mask_noice)

    ELSEIF (C%choice_SMB_model == 'direct_global' .OR. &
            C%choice_SMB_model == 'direct_regional') THEN
      ! Use a directly prescribed global/regional SMB

      CALL run_SMB_model_direct( mesh, climate_matrix%SMB_direct, SMB, time, mask_noice)

    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM( C%choice_SMB_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model
  SUBROUTINE initialise_SMB_model( mesh, ice, SMB, region_name, restart)
    ! Allocate memory for the data fields of the SMB model.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_restart_data),             INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_SMB_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '  Initialising regional SMB model "', TRIM(C%choice_SMB_model), '"...'

    ! Allocate shared memory
    IF     (C%choice_SMB_model == 'uniform' .OR. &
            C%choice_SMB_model == 'idealised' .OR. &
            C%choice_SMB_model == 'direct_global' .OR. &
            C%choice_SMB_model == 'direct_regional') THEN
      ! Only need yearly total SMB in these cases

      CALL allocate_shared_dp_1D( mesh%nV, SMB%SMB_year, SMB%wSMB_year)

    ELSEIF (C%choice_SMB_model == 'IMAU-ITM' .OR. &
            C%choice_SMB_model == 'IMAU-ITM_wrongrefreezing') THEN
      ! Allocate memory and initialise some fields for the IMAU-ITM SMB model

      CALL initialise_SMB_model_IMAU_ITM( mesh, ice, SMB, region_name, restart)

    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM( C%choice_SMB_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=22)

  END SUBROUTINE initialise_SMB_model

! == Idealised SMB parameterisations
! ==================================

  SUBROUTINE run_SMB_model_idealised( mesh, SMB, time, mask_noice)
    ! Run the selected SMB model.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                 INTENT(INOUT) :: SMB
    REAL(dp),                             INTENT(IN)    :: time
    INTEGER,  DIMENSION(:    ),           INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_SMB_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_idealised_SMB == 'EISMINT1_A' .OR. &
            C%choice_idealised_SMB == 'EISMINT1_B' .OR. &
            C%choice_idealised_SMB == 'EISMINT1_C' .OR. &
            C%choice_idealised_SMB == 'EISMINT1_D' .OR. &
            C%choice_idealised_SMB == 'EISMINT1_E' .OR. &
            C%choice_idealised_SMB == 'EISMINT1_F') THEN
      CALL run_SMB_model_idealised_EISMINT1( mesh, SMB, time, mask_noice)
    ELSEIF (C%choice_idealised_SMB == 'Bueler') THEN
      CALL run_SMB_model_idealised_Bueler( mesh, SMB, time, mask_noice)
    ELSEIF (C%choice_idealised_SMB == 'BIVMIP_B') THEN
      CALL run_SMB_model_idealised_BIVMIP_B( mesh, SMB, mask_noice)
    ELSE
      CALL crash('unknown choice_idealised_SMB "' // TRIM( C%choice_idealised_SMB) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model_idealised
  SUBROUTINE run_SMB_model_idealised_EISMINT1( mesh, SMB, time, mask_noice)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    REAL(dp),                            INTENT(IN)    :: time
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_SMB_model_idealised_EISMINT1'
    INTEGER                                            :: vi
    REAL(dp)                                           :: E               ! Radius of circle where accumulation is M_max
    REAL(dp)                                           :: dist            ! distance to centre of circle
    REAL(dp)                                           :: S_b             ! Gradient of accumulation-rate change with horizontal distance
    REAL(dp)                                           :: M_max           ! Maximum accumulation rate

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Default EISMINT configuration
    E         = 450000._dp
    S_b       = 0.01_dp / 1000._dp
    M_max     = 0.5_dp

    IF     (C%choice_idealised_SMB == 'EISMINT1_A') THEN ! Moving margin, steady state
      ! No changes
    ELSEIF (C%choice_idealised_SMB == 'EISMINT1_B') THEN ! Moving margin, 20 kyr
      IF (time < 0._dp) THEN
        ! No changes; first 120 kyr are initialised with EISMINT_1
      ELSE
        E         = 450000._dp + 100000._dp * SIN( 2._dp * pi * time / 20000._dp)
      END IF
    ELSEIF (C%choice_idealised_SMB == 'EISMINT1_C') THEN ! Moving margin, 40 kyr
      IF (time < 0._dp) THEN
        ! No changes; first 120 kyr are initialised with EISMINT_1
      ELSE
        E         = 450000._dp + 100000._dp * SIN( 2._dp * pi * time / 40000._dp)
      END IF
    ELSEIF (C%choice_idealised_SMB == 'EISMINT1_D') THEN ! Fixed margin, steady state
      M_max       = 0.3_dp
      E           = 999000._dp
    ELSEIF (C%choice_idealised_SMB == 'EISMINT1_E') THEN ! Fixed margin, 20 kyr
      IF (time < 0._dp) THEN
        M_max     = 0.3_dp
        E         = 999000._dp
      ELSE
        M_max     = 0.3_dp + 0.2_dp * SIN( 2._dp * pi * time / 20000._dp)
        E         = 999000._dp
      END IF
    ELSEIF (C%choice_idealised_SMB == 'EISMINT1_F') THEN ! Fixed margin, 40 kyr
      IF (time < 0._dp) THEN
        M_max     = 0.3_dp
        E         = 999000._dp
      ELSE
        M_max     = 0.3_dp + 0.2_dp * SIN( 2._dp * pi * time / 40000._dp)
        E         = 999000._dp
      END IF
    END IF

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice( vi) == 0) THEN
        dist = NORM2( mesh%V( vi,:))
        SMB%SMB_year( vi) = MIN( M_max, S_b * (E - dist))
      ELSE
        SMB%SMB_year( vi) = 0._dp
      END IF
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

    END SUBROUTINE run_SMB_model_idealised_EISMINT1
    SUBROUTINE run_SMB_model_idealised_Bueler( mesh, SMB, time, mask_noice)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    REAL(dp),                            INTENT(IN)    :: time
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_SMB_model_idealised_Bueler'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice( vi) == 0) THEN
        SMB%SMB_year( vi) = Bueler_solution_MB( mesh%V(vi,1), mesh%V(vi,2), time)
      ELSE
        SMB%SMB_year( vi) = 0._dp
      END IF
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model_idealised_Bueler
  SUBROUTINE run_SMB_model_idealised_BIVMIP_B( mesh, SMB, mask_noice)
    ! Almost the same as the EISMINT1 moving-margin experiment,
    ! but slightly smaller so the ice lobe doesn't reach the domain border

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_SMB_model_idealised_BIVMIP_B'
    INTEGER                                            :: vi

    REAL(dp)                                           :: E               ! Radius of circle where accumulation is M_max
    REAL(dp)                                           :: dist            ! distance to centre of circle
    REAL(dp)                                           :: S_b             ! Gradient of accumulation-rate change with horizontal distance
    REAL(dp)                                           :: M_max           ! Maximum accumulation rate

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Default EISMINT configuration
    E         = 400000._dp
    S_b       = 0.01_dp / 1000._dp
    M_max     = 0.5_dp

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice( vi) == 0) THEN
        dist = NORM2( mesh%V( vi,:))
        SMB%SMB_year( vi) = MIN( M_max, S_b * (E - dist))
      ELSE
        SMB%SMB_year( vi) = 0._dp
      END IF
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model_idealised_BIVMIP_B

! == The IMAU-ITM SMB model
! =========================

  SUBROUTINE run_SMB_model_IMAUITM( mesh, ice, climate, SMB, mask_noice)
    ! Run the IMAU-ITM SMB model.

    ! NOTE: all the SMB components are in meters of water equivalent;
    !       the end result (SMB and SMB_year) are in meters of ice equivalent.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_snapshot_regional), INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(INOUT) :: SMB
    INTEGER,  DIMENSION(:    ),           INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'run_SMB_model_IMAUITM'
    INTEGER                                             :: vi, m
    INTEGER                                             :: mprev
    REAL(dp)                                            :: snowfrac, liquid_water, sup_imp_wat

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Make sure this routine is called correctly
    IF (.NOT. C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL crash('should only be called when choice_SMB_model == "IMAU-ITM"!')
    END IF

    DO vi = mesh%vi1, mesh%vi2

      ! Background albedo
      SMB%AlbedoSurf( vi) = albedo_soil
      IF ((ice%mask_ocean_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 0) .OR. mask_noice( vi) == 1) SMB%AlbedoSurf( vi) = albedo_water
      IF ( ice%mask_ice_a(   vi) == 1) SMB%AlbedoSurf( vi) = albedo_ice

      DO m = 1, 12  ! Month loop

        mprev = m - 1
        IF (mprev==0) mprev = 12

        SMB%Albedo( vi,m) = MIN(albedo_snow, MAX( SMB%AlbedoSurf( vi), albedo_snow - (albedo_snow - SMB%AlbedoSurf( vi))  * &
                             EXP(-15._dp * SMB%FirnDepth( vi,mprev)) - 0.015_dp * SMB%MeltPreviousYear( vi)))
        IF ((ice%mask_ocean_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 0) .OR. mask_noice( vi) == 1) SMB%Albedo( vi,m) = albedo_water

        ! Determine albation as function af surface temperature and albedo/insolation
        ! according to Bintanja et al. (2002)

        IF (C%do_SMB_IMAUITM_inversion) THEN

          SMB%Melt( vi,m) = MAX(0._dp, ( SMB%C_abl_Ts_inv( vi)         * MAX(0._dp, climate%T2m( vi,m) - T0) + &
                                         SMB%C_abl_Q_inv( vi)          * (1.0_dp - SMB%Albedo( vi,m)) * climate%Q_TOA( vi,m) - &
                                         SMB%C_abl_constant_inv( vi))  * sec_per_year / (L_fusion * 1000._dp * 12._dp))
        ELSE

          SMB%Melt( vi,m) = MAX(0._dp, ( SMB%C_abl_Ts                  * MAX(0._dp, climate%T2m( vi,m) - T0) + &
                                         SMB%C_abl_Q                   * (1.0_dp - SMB%Albedo( vi,m)) * climate%Q_TOA( vi,m) - &
                                         SMB%C_abl_constant)           * sec_per_year / (L_fusion * 1000._dp * 12._dp))
        END IF

        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999),
        ! liquid water content (rain and melt water) and snowdepth

        ! NOTE: commented version is the old ANICE version, supposedly based on "physics" (which we cant check), but
        !       the new version was tuned to RACMO output and produced significantly better snow fractions...

       ! snowfrac = MAX(0._dp, MIN(1._dp, 0.5_dp   * (1 - ATAN((climate%T2m(vi,m) - T0) / 3.5_dp)  / 1.25664_dp)))
        snowfrac = MAX(0._dp, MIN(1._dp, 0.725_dp * (1 - ATAN((climate%T2m( vi,m) - T0) / 5.95_dp) / 1.8566_dp)))

        SMB%Snowfall( vi,m) = climate%Precip( vi,m) *          snowfrac
        SMB%Rainfall( vi,m) = climate%Precip( vi,m) * (1._dp - snowfrac)

        ! Refreezing, according to Janssens & Huybrechts, 2000)
        ! The refreezing (=effective retention) is the minimum value of the amount of super imposed
        ! water and the available liquid water, with a maximum value of the total precipitation.
        ! (see also Huybrechts & de Wolde, 1999)

        ! Add this month's snow accumulation to next month's initial snow depth.
        SMB%AddedFirn( vi,m) = SMB%Snowfall( vi,m) - SMB%Melt( vi,m)
        SMB%FirnDepth( vi,m) = MIN(10._dp, MAX(0._dp, SMB%FirnDepth( vi,mprev) + SMB%AddedFirn( vi,m) ))

      END DO ! DO m = 1, 12

      ! Calculate refrezzing for the whole year, divide equally over the 12 months, then calculate resulting runoff and SMB.
      ! This resolves the problem with refreezing, where liquid water is mostly available in summer
      ! but "refreezing potential" mostly in winter, and there is no proper meltwater retention.

      sup_imp_wat  = SMB%C_refr * MAX(0._dp, T0 - SUM(climate%T2m( vi,:))/12._dp)
      liquid_water = SUM(SMB%Rainfall( vi,:)) + SUM(SMB%Melt( vi,:))

      SMB%Refreezing_year( vi) = MIN( MIN( sup_imp_wat, liquid_water), SUM(climate%Precip( vi,:)))
      IF (ice%mask_ice_a( vi)==0) SMB%Refreezing_year( vi) = 0._dp

      DO m = 1, 12
        SMB%Refreezing( vi,m) = SMB%Refreezing_year( vi) / 12._dp
        SMB%Runoff(     vi,m) = SMB%Melt( vi,m) + SMB%Rainfall( vi,m) - SMB%Refreezing( vi,m)
        SMB%SMB(        vi,m) = SMB%Snowfall( vi,m) + SMB%Refreezing( vi,m) - SMB%Melt( vi,m)
      END DO

      SMB%SMB_year( vi) = SUM(SMB%SMB( vi,:))

      ! Calculate total melt over this year, to be used for determining next year's albedo
      SMB%MeltPreviousYear( vi) = SUM(SMB%Melt( vi,:))

    END DO
    CALL sync

    ! Convert final SMB from water to ice equivalent
    SMB%SMB(      mesh%vi1:mesh%vi2,:) = SMB%SMB(      mesh%vi1:mesh%vi2,:) * 1000._dp / ice_density
    SMB%SMB_year( mesh%vi1:mesh%vi2  ) = SMB%SMB_year( mesh%vi1:mesh%vi2  ) * 1000._dp / ice_density
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( SMB%AlbedoSurf      , 'SMB%AlbedoSurf'      )
    CALL check_for_NaN_dp_2D( SMB%Albedo          , 'SMB%Albedo'          )
    CALL check_for_NaN_dp_2D( SMB%Melt            , 'SMB%Melt'            )
    CALL check_for_NaN_dp_2D( SMB%Snowfall        , 'SMB%Snowfall'        )
    CALL check_for_NaN_dp_2D( SMB%Rainfall        , 'SMB%Rainfall'        )
    CALL check_for_NaN_dp_2D( SMB%Refreezing      , 'SMB%Refreezing'      )
    CALL check_for_NaN_dp_2D( SMB%Runoff          , 'SMB%Runoff'          )
    CALL check_for_NaN_dp_2D( SMB%SMB             , 'SMB%SMB'             )
    CALL check_for_NaN_dp_2D( SMB%AddedFirn       , 'SMB%AddedFirn'       )
    CALL check_for_NaN_dp_2D( SMB%FirnDepth       , 'SMB%FirnDepth'       )
    CALL check_for_NaN_dp_1D( SMB%SMB_year        , 'SMB%SMB_year'        )
    CALL check_for_NaN_dp_1D( SMB%MeltPreviousYear, 'SMB%MeltPreviousYear')
    CALL check_for_NaN_dp_1D( SMB%Albedo_year     , 'SMB%Albedo_year'     )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model_IMAUITM
  SUBROUTINE run_SMB_model_IMAUITM_wrongrefreezing( mesh, ice, climate, SMB, mask_noice)
    ! Run the IMAU-ITM SMB model. Old version, exactly as it was in ANICE2.1 (so with the "wrong" refreezing)

    ! NOTE: all the SMB components and the total are in meters of water equivalent

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_snapshot_regional), INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(INOUT) :: SMB
    INTEGER,  DIMENSION(:    ),           INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'run_SMB_model_IMAUITM_wrongrefreezing'
    INTEGER                                             :: vi, m
    INTEGER                                             :: mprev
    REAL(dp)                                            :: snowfrac, liquid_water, sup_imp_wat

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Make sure this routine is called correctly
    IF (.NOT. C%choice_SMB_model == 'IMAU-ITM_wrongrefreezing') THEN
      CALL crash('should only be called when choice_SMB_model == "IMAU-ITM_wrongrefreezing"!')
    END IF

    DO vi = mesh%vi1, mesh%vi2

      ! "Background" albedo (= surface without any firn, so either ice, soil, or water)
      SMB%AlbedoSurf( vi) = albedo_soil
      IF ((ice%mask_ocean_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 0) .OR. mask_noice( vi) == 1) SMB%AlbedoSurf( vi) = albedo_water
      IF (ice%mask_ice_a(    vi) == 1) SMB%AlbedoSurf( vi) = albedo_ice

      DO m = 1, 12  ! Month loop

        mprev = m - 1
        IF (mprev == 0) mprev = 12

        SMB%Albedo( vi,m) = MIN(albedo_snow, MAX( SMB%AlbedoSurf( vi), albedo_snow - (albedo_snow - SMB%AlbedoSurf( vi))  * &
                             EXP(-15._dp * SMB%FirnDepth( vi,mprev)) - 0.015_dp * SMB%MeltPreviousYear( vi)))
        IF ((ice%mask_ocean_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 0) .OR. mask_noice( vi) == 1) SMB%Albedo( vi,m) = albedo_water

        ! Determine ablation as function af surface temperature and albedo/insolation
        ! according to Bintanja et al. (2002)

        SMB%Melt( vi,m) = MAX(0._dp, ( SMB%C_abl_Ts         * (climate%T2m( vi,m) - T0) + &
                                       SMB%C_abl_Q          * (1.0_dp - SMB%Albedo( vi,m)) * climate%Q_TOA( vi,m) - &
                                       SMB%C_abl_constant)  * sec_per_year / (L_fusion * 1000._dp * 12._dp))

        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999),
        ! liquid water content (rain and melt water) and snowdepth
        snowfrac = MAX(0._dp, MIN(1._dp, 0.5_dp   * (1 - ATAN((climate%T2m( vi,m) - T0) / 3.5_dp)  / 1.25664_dp)))

        SMB%Snowfall( vi,m) = climate%Precip( vi,m) *          snowfrac
        SMB%Rainfall( vi,m) = climate%Precip( vi,m) * (1._dp - snowfrac)

        ! Refreezing, according to Janssens & Huybrechts, 2000)
        ! The refreezing (=effective retention) is the minimum value of the amount of super imposed
        ! water and the available liquid water, with a maximum value of the total precipitation.
        ! (see also Huybrechts & de Wolde, 1999)
        sup_imp_wat  = 0.012_dp * MAX(0._dp, T0 - climate%T2m( vi,m))
        liquid_water = SMB%Rainfall( vi,m) + SMB%Melt( vi,m)
        SMB%Refreezing( vi,m) = MIN( MIN( sup_imp_wat, liquid_water), climate%Precip( vi,m))
        IF (ice%mask_ice_a( vi) == 0 .OR. mask_noice( vi) == 1) SMB%Refreezing( vi,m) = 0._dp

        ! Calculate runoff and total SMB
        SMB%Runoff( vi,m) = SMB%Melt(     vi,m) + SMB%Rainfall(   vi,m) - SMB%Refreezing( vi,m)
        SMB%SMB(    vi,m) = SMB%Snowfall( vi,m) + SMB%Refreezing( vi,m) - SMB%Melt(       vi,m)

        ! Add this month's snow accumulation to next month's initial snow depth.
        SMB%AddedFirn( vi,m) = SMB%Snowfall( vi,m) - SMB%Melt( vi,m)
        SMB%FirnDepth( vi,m) = MIN(10._dp, MAX(0._dp, SMB%FirnDepth( vi,mprev) + SMB%AddedFirn( vi,m) ))

      END DO ! DO m = 1, 12

      ! Calculate total SMB for the entire year
      SMB%SMB_year( vi) = SUM(SMB%SMB( vi,:))

      ! Calculate total melt over this year, to be used for determining next year's albedo
      SMB%MeltPreviousYear( vi) = SUM(SMB%Melt( vi,:))

      ! Calculate yearly mean albedo (diagnostic only)
      SMB%Albedo_year( vi) = SUM(SMB%Albedo( vi,:)) / 12._dp

    END DO
    CALL sync

    ! Convert final SMB from water to ice equivalent
    SMB%SMB(      mesh%vi1:mesh%vi2,:) = SMB%SMB(      mesh%vi1:mesh%vi1,:) * 1000._dp / ice_density
    SMB%SMB_year( mesh%vi1:mesh%vi2  ) = SMB%SMB_year( mesh%vi1:mesh%vi2  ) * 1000._dp / ice_density
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( SMB%AlbedoSurf      , 'SMB%AlbedoSurf'      )
    CALL check_for_NaN_dp_2D( SMB%Albedo          , 'SMB%Albedo'          )
    CALL check_for_NaN_dp_2D( SMB%Melt            , 'SMB%Melt'            )
    CALL check_for_NaN_dp_2D( SMB%Snowfall        , 'SMB%Snowfall'        )
    CALL check_for_NaN_dp_2D( SMB%Rainfall        , 'SMB%Rainfall'        )
    CALL check_for_NaN_dp_2D( SMB%Refreezing      , 'SMB%Refreezing'      )
    CALL check_for_NaN_dp_2D( SMB%Runoff          , 'SMB%Runoff'          )
    CALL check_for_NaN_dp_2D( SMB%SMB             , 'SMB%SMB'             )
    CALL check_for_NaN_dp_2D( SMB%AddedFirn       , 'SMB%AddedFirn'       )
    CALL check_for_NaN_dp_2D( SMB%FirnDepth       , 'SMB%FirnDepth'       )
    CALL check_for_NaN_dp_1D( SMB%SMB_year        , 'SMB%SMB_year'        )
    CALL check_for_NaN_dp_1D( SMB%MeltPreviousYear, 'SMB%MeltPreviousYear')
    CALL check_for_NaN_dp_1D( SMB%Albedo_year     , 'SMB%Albedo_year'     )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model_IMAUITM_wrongrefreezing
  SUBROUTINE initialise_SMB_model_IMAU_ITM( mesh, ice, SMB, region_name, restart)
    ! Allocate memory for the data fields of the SMB model.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_restart_data),             INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_SMB_model_IMAU_ITM'
    INTEGER                                            :: vi
    CHARACTER(LEN=256)                                 :: SMB_IMAUITM_choice_init_firn

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '   Initialising the IMAU-ITM SMB model...'

    ! Data fields
    CALL allocate_shared_dp_1D( mesh%nV,     SMB%AlbedoSurf,       SMB%wAlbedoSurf      )
    CALL allocate_shared_dp_1D( mesh%nV,     SMB%MeltPreviousYear, SMB%wMeltPreviousYear)
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%FirnDepth,        SMB%wFirnDepth       )
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%Rainfall,         SMB%wRainfall        )
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%Snowfall,         SMB%wSnowfall        )
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%AddedFirn,        SMB%wAddedFirn       )
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%Melt,             SMB%wMelt            )
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%Refreezing,       SMB%wRefreezing      )
    CALL allocate_shared_dp_1D( mesh%nV,     SMB%Refreezing_year,  SMB%wRefreezing_year )
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%Runoff,           SMB%wRunoff          )
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%Albedo,           SMB%wAlbedo          )
    CALL allocate_shared_dp_1D( mesh%nV,     SMB%Albedo_year,      SMB%wAlbedo_year     )
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%SMB,              SMB%wSMB             )
    CALL allocate_shared_dp_1D( mesh%nV,     SMB%SMB_year,         SMB%wSMB_year        )

    ! Tuning parameters
    CALL allocate_shared_dp_0D( SMB%C_abl_constant, SMB%wC_abl_constant)
    CALL allocate_shared_dp_0D( SMB%C_abl_Ts,       SMB%wC_abl_Ts      )
    CALL allocate_shared_dp_0D( SMB%C_abl_Q,        SMB%wC_abl_Q       )
    CALL allocate_shared_dp_0D( SMB%C_refr,         SMB%wC_refr        )

    ! Inversion parameters
    CALL allocate_shared_dp_1D( mesh%nV, SMB%C_abl_constant_inv, SMB%wC_abl_constant_inv)
    CALL allocate_shared_dp_1D( mesh%nV, SMB%C_abl_Ts_inv,       SMB%wC_abl_Ts_inv      )
    CALL allocate_shared_dp_1D( mesh%nV, SMB%C_abl_Q_inv,        SMB%wC_abl_Q_inv       )
    CALL allocate_shared_dp_1D( mesh%nV, SMB%C_refr_inv,         SMB%wC_refr_inv        )

    ! Initialise regional scalar parameters to specified values
    IF (par%master) THEN
      IF     (region_name == 'NAM') THEN
        SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_NAM
        SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_NAM
        SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_NAM
        SMB%C_refr                   = C%SMB_IMAUITM_C_refr_NAM
      ELSEIF (region_name == 'EAS') THEN
        SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_EAS
        SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_EAS
        SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_EAS
        SMB%C_refr                   = C%SMB_IMAUITM_C_refr_EAS
      ELSEIF (region_name == 'GRL') THEN
        SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_GRL
        SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_GRL
        SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_GRL
        SMB%C_refr                   = C%SMB_IMAUITM_C_refr_GRL
      ELSEIF (region_name == 'ANT') THEN
        SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_ANT
        SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_ANT
        SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_ANT
        SMB%C_refr                   = C%SMB_IMAUITM_C_refr_ANT
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Initialise regional 1-D inversion parameters to specified values
    IF (par%master) THEN
      IF     (region_name == 'NAM') THEN
        SMB%C_abl_constant_inv       = C%SMB_IMAUITM_C_abl_constant_NAM
        SMB%C_abl_Ts_inv             = C%SMB_IMAUITM_C_abl_Ts_NAM
        SMB%C_abl_Q_inv              = C%SMB_IMAUITM_C_abl_Q_NAM
        SMB%C_refr_inv               = C%SMB_IMAUITM_C_refr_NAM
      ELSEIF (region_name == 'EAS') THEN
        SMB%C_abl_constant_inv       = C%SMB_IMAUITM_C_abl_constant_EAS
        SMB%C_abl_Ts_inv             = C%SMB_IMAUITM_C_abl_Ts_EAS
        SMB%C_abl_Q_inv              = C%SMB_IMAUITM_C_abl_Q_EAS
        SMB%C_refr_inv               = C%SMB_IMAUITM_C_refr_EAS
      ELSEIF (region_name == 'GRL') THEN
        SMB%C_abl_constant_inv       = C%SMB_IMAUITM_C_abl_constant_GRL
        SMB%C_abl_Ts_inv             = C%SMB_IMAUITM_C_abl_Ts_GRL
        SMB%C_abl_Q_inv              = C%SMB_IMAUITM_C_abl_Q_GRL
        SMB%C_refr_inv               = C%SMB_IMAUITM_C_refr_GRL
      ELSEIF (region_name == 'ANT') THEN
        SMB%C_abl_constant_inv       = C%SMB_IMAUITM_C_abl_constant_ANT
        SMB%C_abl_Ts_inv             = C%SMB_IMAUITM_C_abl_Ts_ANT
        SMB%C_abl_Q_inv              = C%SMB_IMAUITM_C_abl_Q_ANT
        SMB%C_refr_inv               = C%SMB_IMAUITM_C_refr_ANT
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Initialisation choice
    IF     (region_name == 'NAM') THEN
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_NAM
    ELSEIF (region_name == 'EAS') THEN
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_EAS
    ELSEIF (region_name == 'GRL') THEN
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_GRL
    ELSEIF (region_name == 'ANT') THEN
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_ANT
    END IF

    ! Initialise the firn layer
    IF     (SMB_IMAUITM_choice_init_firn == 'uniform') THEN
      ! Initialise with a uniform firn layer over the ice sheet

      DO vi = mesh%vi1, mesh%vi2
        IF (ice%Hi_a( vi) > 0._dp) THEN
          SMB%FirnDepth(        vi,:) = C%SMB_IMAUITM_initial_firn_thickness
          SMB%MeltPreviousYear( vi  ) = 0._dp
        ELSE
          SMB%FirnDepth(        vi,:) = 0._dp
          SMB%MeltPreviousYear( vi  ) = 0._dp
        END IF
      END DO
      CALL sync

    ELSEIF (SMB_IMAUITM_choice_init_firn == 'restart') THEN
      ! Initialise with the firn layer of a previous run

      IF (par%master) WRITE (0,*) '    Initialising firn layer using data read from restart file...'

      DO vi = mesh%vi1, mesh%vi2
        IF (ice%Hi_a( vi) > 0._dp) THEN
          SMB%FirnDepth(        vi,:) = restart%FirnDepth(        vi,:)
          SMB%MeltPreviousYear( vi  ) = restart%MeltPreviousYear( vi  )
        ELSE
          SMB%FirnDepth(        vi,:) = 0._dp
          SMB%MeltPreviousYear( vi  ) = 0._dp
        END IF
      END DO
      CALL sync

    ELSE
      CALL crash('unknown SMB_IMAUITM_choice_init_firn "' // TRIM( SMB_IMAUITM_choice_init_firn) // '"!')
    END IF

    ! Initialise albedo
    DO vi = mesh%vi1, mesh%vi2

      ! Background albedo
      IF (ice%Hb_a( vi) < 0._dp) THEN
        SMB%AlbedoSurf( vi) = albedo_water
      ELSE
        SMB%AlbedoSurf( vi) = albedo_soil
      END IF
      IF (ice%Hi_a( vi) > 0._dp) THEN
        SMB%AlbedoSurf(  vi) = albedo_snow
      END IF

      SMB%Albedo( vi,:) = SMB%AlbedoSurf( vi)

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=22)

  END SUBROUTINE initialise_SMB_model_IMAU_ITM

! == Directly prescribed global/regional SMB
! ==========================================

  SUBROUTINE run_SMB_model_direct( mesh, SMB_direct, SMB, time, mask_noice)
    ! Run the selected SMB model: direct global/regional SMB forcing.
    !
    ! NOTE: the whole business of reading the data from the NetCDF file and mapping
    !       it to the model mesh is handled by the climate_module!
    ! NOTE ALSO: since the climate_module routines for the "direct_global" option
    !       already map the results to the model mesh (located in region%climate_matrix%SMB_direct),
    !       in this routine here we can treat "direct_global" and "direct_regional" the same way

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_direct_SMB_forcing_regional), INTENT(IN)    :: SMB_direct
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    REAL(dp),                               INTENT(IN)    :: time
    INTEGER,  DIMENSION(:    ),             INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_SMB_model_direct'
    REAL(dp)                                           :: wt0, wt1
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Interpolate the two timeframes in time
    wt0 = (SMB_direct%t1 - time) / (SMB_direct%t1 - SMB_direct%t0)
    wt1 = 1._dp - wt0

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice( vi) == 0) THEN
        SMB%SMB_year( vi) = (wt0 * SMB_direct%SMB_year0( vi)) + (wt1 * SMB_direct%SMB_year1( vi))
      ELSE
        SMB%SMB_year( vi) = 0._dp
      END IF
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model_direct

! == Remapping after mesh update
! ==============================

  SUBROUTINE remap_SMB_model( mesh_old, mesh_new, map, SMB)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_SMB_model'
    INTEGER                                            :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV
    int_dummy = int_dummy

    IF (C%choice_SMB_model == 'IMAU-ITM') THEN
      ! Firn depth and melt-during-previous-year must be remapped
      CALL remap_field_dp_2D( mesh_old, mesh_new, map, SMB%MeltPreviousYear, SMB%wMeltPreviousYear, 'trilin')
      CALL remap_field_dp_3D( mesh_old, mesh_new, map, SMB%FirnDepth,        SMB%wFirnDepth,        'trilin')

      ! IF (C%do_SMB_IMAUITM_inversion) THEN
        ! Remap inverted IMAU-ITM paramaters so they are not reset after a mesh update
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, SMB%C_abl_constant_inv, SMB%wC_abl_constant_inv, 'trilin')
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, SMB%C_abl_Ts_inv,       SMB%wC_abl_Ts_inv,       'trilin')
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, SMB%C_abl_Q_inv,        SMB%wC_abl_Q_inv,        'trilin')
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, SMB%C_refr_inv,         SMB%wC_refr_inv,         'trilin')
      ! END IF

      ! Reallocate rather than remap; after a mesh update we'll immediately run the SMB model anyway
      CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%Q_TOA,            SMB%wQ_TOA           )
      CALL reallocate_shared_dp_1D( mesh_new%nV,     SMB%AlbedoSurf,       SMB%wAlbedoSurf      )
      CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%Rainfall,         SMB%wRainfall        )
      CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%Snowfall,         SMB%wSnowfall        )
      CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%AddedFirn,        SMB%wAddedFirn       )
      CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%Melt,             SMB%wMelt            )
      CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%Refreezing,       SMB%wRefreezing      )
      CALL reallocate_shared_dp_1D( mesh_new%nV,     SMB%Refreezing_year,  SMB%wRefreezing_year )
      CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%Runoff,           SMB%wRunoff          )
      CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%Albedo,           SMB%wAlbedo          )
      CALL reallocate_shared_dp_1D( mesh_new%nV,     SMB%Albedo_year,      SMB%wAlbedo_year     )
      CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%SMB,              SMB%wSMB             )
    END IF

    CALL reallocate_shared_dp_1D(   mesh_new%nV,     SMB%SMB_year,         SMB%wSMB_year        )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_SMB_model

! == Inversion
! ============

  SUBROUTINE SMB_IMAUITM_inversion( mesh, ice, climate, SMB, refgeo)
    ! Iteratively invert for IMAU-ITM SMB parameters over the entire domain

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_snapshot_regional), INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(INOUT) :: SMB
    TYPE(type_reference_geometry),        INTENT(IN)    :: refgeo

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'SMB_IMAUITM_inversion'
    INTEGER                                            :: vi
    REAL(dp)                                           :: h_scale, h_delta, h_dfrac, new_val

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise ice thickness scaling
    h_scale = 1.0_dp/C%SMB_IMAUITM_inv_scale

    DO vi = mesh%vi1, mesh%vi2

      ! Invert only where the model has ice
      IF (ice%mask_ice_a( vi) == 1) THEN

          ! Compute ice thickness difference
          h_delta = ice%Hi_a( vi) - refgeo%Hi( vi)
          ! To what fraction of the ice thickness this difference corresponds
          h_dfrac = h_delta / MAX(refgeo%Hi( vi), 1._dp)
          ! Scale down the difference to the range [-1.5, 1.5] (in the style of Pollard & DeConto 2012)
          h_delta = MAX(-1.5_dp, MIN(1.5_dp, h_delta * h_scale))

          ! Further adjust only where the current value is not improving the result
          IF ( (h_delta > 0._dp .AND. ice%dHi_dt_a( vi) >= 0._dp) .OR. &
               (h_delta < 0._dp .AND. ice%dHi_dt_a( vi) <= 0._dp) ) THEN

            ! NOTE: For some reason, the values for the C_abl_constant parameter are negative
            ! for more melting and positive for less (final total melt never negative tho). This
            ! is in contrast to the other parameters, which use positve values to have a stronger
            ! melting (or refreezin) effect. Thus, pay attention to the sign of the adjustment.

            ! == Constant ablation
            ! ====================

            ! If any month has temperatures above or close enough to melting point
            IF ( ANY(climate%T2m(vi,:) - T0 > -5._dp) ) THEN

              ! Adjust parameter
              IF (C%SMB_IMAUITM_inv_C_abl_constant_min /= C%SMB_IMAUITM_inv_C_abl_constant_max) THEN
                ! Get value to be adjusted
                new_val = SMB%C_abl_constant_inv( vi)
                ! Adjust value based on ice thickness difference (note the negative sign here)
                new_val = new_val - 1.72476_dp * TAN(h_delta)
                ! Limit adjusted value
                new_val = MAX(new_val, C%SMB_IMAUITM_inv_C_abl_constant_min)
                new_val = MIN(new_val, C%SMB_IMAUITM_inv_C_abl_constant_max)
                ! Assign adjusted value to the SMB parameter
                SMB%C_abl_constant_inv( vi) = new_val
              END IF

            END IF ! ( ANY(climate%T2m(vi,:) - T0 > 0._dp) )

            ! == Temperature-based ablation
            ! =============================

            ! If any month has melting temperatures
            IF ( ANY(climate%T2m(vi,:) - T0 > 0._dp) ) THEN

              ! Adjust parameter
              IF (C%SMB_IMAUITM_inv_C_abl_Ts_min /= C%SMB_IMAUITM_inv_C_abl_Ts_max) THEN
                new_val = SMB%C_abl_Ts_inv( vi)
                new_val = new_val + 1.72476_dp * TAN(h_delta) ! Note the positive sign here
                new_val = MAX(new_val, C%SMB_IMAUITM_inv_C_abl_Ts_min)
                new_val = MIN(new_val, C%SMB_IMAUITM_inv_C_abl_Ts_max)

                SMB%C_abl_Ts_inv( vi) = new_val
              END IF

            END IF ! ( ANY(climate%T2m(vi,:) - T0 > 0._dp) )

            ! == Insolation-based ablation
            ! ============================

            ! If any month has less than maximum albedo
            IF ( ANY(1.0_dp - SMB%Albedo( vi,:) > 0._dp) ) THEN

              ! Adjust parameter
              IF (C%SMB_IMAUITM_inv_C_abl_Q_min /= C%SMB_IMAUITM_inv_C_abl_Q_max) THEN
                new_val = SMB%C_abl_Q_inv( vi)
                new_val = new_val + 0.01724_dp * TAN(h_delta) ! Note the positive sign here
                new_val = MAX(new_val, C%SMB_IMAUITM_inv_C_abl_Q_min)
                new_val = MIN(new_val, C%SMB_IMAUITM_inv_C_abl_Q_max)

                SMB%C_abl_Q_inv( vi) = new_val
              END IF

            END IF ! ( ANY(1.0_dp - SMB%Albedo( vi,:) > 0._dp) )

            ! == Refreezing
            ! =============

            ! If any month has freezing temperatures
            IF ( ANY(climate%T2m(vi,:) - T0 < 0._dp) ) THEN

              ! Adjust parameter
              IF (C%SMB_IMAUITM_inv_C_refr_min /= C%SMB_IMAUITM_inv_C_refr_max) THEN
                new_val = SMB%C_refr_inv( vi)
                new_val = new_val - 0.01724_dp * TAN(h_delta) ! Note the negative sign here
                new_val = MAX(new_val, C%SMB_IMAUITM_inv_C_refr_min)
                new_val = MIN(new_val, C%SMB_IMAUITM_inv_C_refr_max)

                SMB%C_refr_inv( vi) = new_val
              END IF

            END IF ! ( ANY(climate%T2m(vi,:) - T0 < 0._dp) )

          END IF ! else the fit is already improving for some other reason, so leave it alone

      END IF ! else the model does not have ice there, so leave it alone

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE SMB_IMAUITM_inversion

! == Some generally useful tools
! ==============================

  FUNCTION Bueler_solution_MB( x, y, t) RESULT(M)
    ! Describes an ice-sheet at time t (in years) conforming to the Bueler solution
    ! with dome thickness H0 and margin radius R0 at t0, with a surface mass balance
    ! determined by lambda.

    ! Input variables
    REAL(dp), INTENT(IN) :: x       ! x coordinate [m]
    REAL(dp), INTENT(IN) :: y       ! y coordinate [m]
    REAL(dp), INTENT(IN) :: t       ! Time from t0 [years]

    ! Result
    REAL(dp)             :: M

    ! Local variables
    REAL(dp) :: A_flow, rho, g, n, alpha, beta, Gamma, f1, f2, t0, tp, f3, f4, H

    REAL(dp), PARAMETER :: H0     = 3000._dp    ! Ice dome thickness at t=0 [m]
    REAL(dp), PARAMETER :: R0     = 500000._dp  ! Ice margin radius  at t=0 [m]
    REAL(dp), PARAMETER :: lambda = 5.0_dp      ! Mass balance parameter

    A_flow  = 1E-16_dp
    rho     = 910._dp
    g       = 9.81_dp
    n       = 3._dp

    alpha = (2._dp - (n+1._dp)*lambda) / ((5._dp*n)+3._dp)
    beta  = (1._dp + ((2._dp*n)+1._dp)*lambda) / ((5._dp*n)+3._dp)
    Gamma = 2._dp/5._dp * (A_flow/sec_per_year) * (rho * g)**n

    f1 = ((2._dp*n)+1)/(n+1._dp)
    f2 = (R0**(n+1._dp))/(H0**((2._dp*n)+1._dp))
    t0 = (beta / Gamma) * (f1**n) * f2

    !tp = (t * sec_per_year) + t0; % Acutal equation needs t in seconds from zero , but we want to supply t in years from t0
    tp = t * sec_per_year

    f1 = (tp / t0)**(-alpha)
    f2 = (tp / t0)**(-beta)
    f3 = SQRT( (x**2._dp) + (y**2._dp) )/R0
    f4 = MAX(0._dp, 1._dp - (f2*f3)**((n+1._dp)/n))
    H = H0 * f1 * f4**(n/((2._dp*n)+1._dp))

    M = (lambda / tp) * H * sec_per_year

  END FUNCTION Bueler_solution_MB

END MODULE SMB_module
