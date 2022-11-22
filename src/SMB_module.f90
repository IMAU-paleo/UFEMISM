module SMB_module
  ! All the routines for calculating the surface mass balance from the climate forcing

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module, only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module,    only : ice_density, L_fusion, sec_per_year, T0, pi
  use parallel_module,      only : par, sync, ierr, cerr, partition_list
  use utilities_module,     only : check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                   check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  use netcdf_module,        only : debug, inquire_restart_file_SMB, &
                                   read_restart_file_SMB
  use data_types_module,    only : type_mesh, type_ice_model, &
                                   type_SMB_model, type_remapping_mesh_mesh, &
                                   type_climate_matrix_regional, &
                                   type_climate_snapshot_regional, &
                                   type_restart_data
  use forcing_module,       only : forcing
  use mesh_mapping_module,  only : remap_field_dp_2D, remap_field_dp_3D, &
                                   calc_remapping_operator_grid2mesh, map_grid2mesh_2D, map_grid2mesh_3D, &
                                   deallocate_remapping_operators_grid2mesh
  use reallocate_mod,       only : reallocate_bounds

  implicit none

  real(dp), parameter :: albedo_water        = 0.1_dp
  real(dp), parameter :: albedo_soil         = 0.2_dp
  real(dp), parameter :: albedo_ice          = 0.5_dp
  real(dp), parameter :: albedo_snow         = 0.85_dp
  real(dp), parameter :: initial_snow_depth  = 0.1_dp

contains

! ===== Main =====
! ================

  subroutine run_SMB_model( mesh, ice, climate_matrix, time, SMB, mask_noice)
    ! Run the selected regional SMB model

    implicit none

    ! In/output variables
    type(type_mesh),                       intent(in)    :: mesh
    type(type_ice_model),                  intent(in)    :: ice
    type(type_climate_matrix_regional),    intent(in)    :: climate_matrix
    real(dp),                              intent(in)    :: time
    type(type_SMB_model),                  intent(inout) :: SMB
    integer, dimension(mesh%vi1:mesh%vi2), intent(in)    :: mask_noice

    ! Local variables:
    character(len=256), parameter                        :: routine_name = 'run_SMB_model'
    integer                                              :: vi
    real(dp)                                             :: R

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)


    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6') THEN
          
        CALL EISMINT_SMB( mesh, time, SMB)
        CALL finalise_routine( routine_name)
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
        SMB%SMB_year( mesh%vi1:mesh%vi2  ) = 0._dp
        SMB%SMB(      mesh%vi1:mesh%vi2,:) = 0._dp
        CALL finalise_routine( routine_name)
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'Bueler') THEN
        DO vi = mesh%vi1, mesh%vi2
          SMB%SMB_year( vi  ) = Bueler_solution_MB( mesh%V(vi,1), mesh%V(vi,2), time)
          SMB%SMB(      vi,:) = SMB%SMB_year( vi) / 12._dp
        END DO
        CALL finalise_routine( routine_name)
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        SMB%SMB_year( mesh%vi1:mesh%vi2  ) = 0.3_dp
        SMB%SMB(      mesh%vi1:mesh%vi2,:) = 0.3_dp / 12._dp
        CALL finalise_routine( routine_name)
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        ! Similar to EISMINT
        DO vi = mesh%vi1, mesh%vi2
          R = NORM2(mesh%V(vi,:))
          IF (R < 250000._dp) THEN
            SMB%SMB_year( vi  ) = 0.3_dp
          ELSE
            SMB%SMB_year( vi  ) = MAX(-2._dp, 0.3_dp - (R - 250000._dp) / 200000._dp)
          END IF
          SMB%SMB(      vi,:) = SMB%SMB_year( vi  ) / 12._dp
        END DO ! DO vi = mesh%vi1, mesh%vi2
        CALL sync
        CALL finalise_routine( routine_name)
        RETURN
      ELSE
        CALL crash('unknown choice_benchmark_experiment "' // TRIM( C%choice_benchmark_experiment) // '"!')
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN

    ! === Choice of model ===
    ! =======================

    select case (C%choice_SMB_model)

      case ('uniform')
        ! Apply a simple uniform SMB

        do vi = mesh%vi1, mesh%vi2
          if (mask_noice( vi) == 0) then
            SMB%SMB_year( vi) = C%SMB_uniform
          else
            SMB%SMB_year( vi) = 0._dp
          end if
        end do

      case ('IMAU-ITM')
        ! Run the IMAU-ITM SMB model
        call run_SMB_model_IMAUITM( mesh, ice, climate_matrix%applied, SMB, mask_noice)

      case default
        !Unknown case
        call crash('unknown choice_SMB_model "' // &
                    trim(C%choice_SMB_model) // '"!')

    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model

  subroutine initialise_SMB_model( mesh, ice, SMB, region_name)
    ! Allocate memory for the data fields of the SMB model.

    implicit none

    ! In/output variables
    type(type_mesh),      intent(inout) :: mesh
    type(type_ice_model), intent(in)    :: ice
    type(type_SMB_model), intent(inout) :: SMB
    character(len=3),     intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_SMB_model'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write (*,"(3A)") '  Initialising regional SMB model "', &
                          trim(C%choice_SMB_model), '"...'
    end if
    call sync

    ! === Choice of model ===
    ! =======================

    select case (C%choice_SMB_model)

      case ('uniform')
        ! Only need yearly total SMB
        allocate( SMB%SMB_year(mesh%vi1:mesh%vi2))

      case ('IMAU-ITM')
        ! Initialise the IMAU-ITM SMB model
        call initialise_SMB_model_IMAU_ITM( mesh, ice, SMB, region_name)

      case default
        !Unknown case
        call crash('unknown choice_SMB_model "' // &
                    trim(C%choice_SMB_model) // '"!')

    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model

! ===== The IMAU-ITM SMB model =====
! ==================================

  subroutine run_SMB_model_IMAUITM( mesh, ice, climate, SMB, mask_noice)
    ! Run the IMAU-ITM SMB model.

    ! NOTE: all the SMB components are in meters of water equivalent;
    !       the end result (SMB and SMB_year) are in meters of ice equivalent.

    implicit none

    ! In/output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_climate_snapshot_regional),   intent(in)    :: climate
    type(type_SMB_model),                   intent(inout) :: SMB
    integer,  dimension(mesh%vi1:mesh%vi2), intent(in)    :: mask_noice

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'run_SMB_model_IMAUITM'
    integer                                               :: vi, m
    integer                                               :: mprev
    real(dp)                                              :: snowfrac, liquid_water, sup_imp_wat

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! === Background albedo ===
      ! =========================

      ! Default value
      SMB%AlbedoSurf( vi) = albedo_soil

      ! Open water value
      if (ice%mask_ocean_a( vi) == 1 .and. &
         (ice%mask_shelf_a( vi) == 0 .or. mask_noice( vi) == 1)) then
        SMB%AlbedoSurf( vi) = albedo_water
      end if

      ! Ice value
      if ( ice%mask_ice_a( vi) == 1) then
        SMB%AlbedoSurf( vi) = albedo_ice
      end if

      ! === Month loop ===
      ! ==================

      do m = 1, 12  ! Month loop

        ! Set previous month index
        ! ========================

        mprev = m - 1

        ! If January, then December
        if (mprev==0) then
          mprev = 12
        end if

        ! Monthly albedo
        ! ==============

        ! Compute value
        SMB%Albedo( vi,m) = min(albedo_snow, max( SMB%AlbedoSurf( vi), &
                                albedo_snow - (albedo_snow - SMB%AlbedoSurf( vi))  * &
                                exp(-15._dp * SMB%FirnDepth( vi,mprev)) - 0.015_dp * SMB%MeltPreviousYear( vi)))

        ! Reset for open water
        if (ice%mask_ocean_a( vi) == 1 .and. &
           (ice%mask_shelf_a( vi) == 0 .or. mask_noice( vi) == 1)) then
          SMB%Albedo( vi,m) = albedo_water
        end if

        ! Monthly ablation
        ! ================

        ! Determine ablation as function af surface temperature and
        ! albedo/insolation according to Bintanja et al. (2002)

        SMB%Melt( vi,m) = max(0._dp, ( SMB%C_abl_Ts        * max(0._dp, climate%T2m( vi,m) - T0) + &
                                       SMB%C_abl_Q         * (1.0_dp - SMB%Albedo( vi,m)) * climate%Q_TOA( vi,m) - &
                                       SMB%C_abl_constant) * sec_per_year / (L_fusion * 1000._dp * 12._dp))

        ! Monthly accumulation
        ! ====================

        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999),
        ! liquid water content (rain and melt water) and snowdepth
        snowfrac = MAX(0._dp, MIN(1._dp, 0.725_dp * (1 - ATAN((climate%T2m( vi,m) - T0) / 5.95_dp) / 1.8566_dp)))

        SMB%Snowfall( vi,m) = climate%Precip( vi,m) *          snowfrac
        SMB%Rainfall( vi,m) = climate%Precip( vi,m) * (1._dp - snowfrac)

        ! Add this month's snow accumulation to next month's initial snow depth.
        SMB%AddedFirn( vi,m) = SMB%Snowfall( vi,m) - SMB%Melt( vi,m)
        SMB%FirnDepth( vi,m) = MIN(10._dp, MAX(0._dp, SMB%FirnDepth( vi,mprev) + SMB%AddedFirn( vi,m) ))

      end do ! m = 1, 12

      ! === Annual Refreezing ===
      ! =========================

      ! According to Janssens & Huybrechts (1999, 2000)
      ! The refreezing (=effective retention) is the minimum value of the amount
      ! of super imposed water and the available liquid water. Calculated for the
      ! whole year, then divided equally over the 12 months. This is done to account
      ! for the fact that liquid water is mostly available in summer, but "refreezing
      ! potential" is mostly available in winter. Not capped at total precipitation
      ! to account for cases where there is more melt than precipitation (e.g. ANT).

      ! Total refreezing estimation
      sup_imp_wat  = SMB%C_refr * max(0._dp, T0 - sum(climate%T2m( vi,:))/12._dp)

      ! Total amount of available liquid water
      liquid_water = sum(SMB%Rainfall( vi,:)) + sum(SMB%Melt( vi,:))

      ! Effective refreezing
      SMB%Refreezing_year( vi) = min( sup_imp_wat, liquid_water)

      ! Limit it to ice areas only
      if (ice%mask_ice_a( vi)==0) then
        SMB%Refreezing_year( vi) = 0._dp
      end if

      ! === Monthly refreezing, runoff, and SMB ===
      ! ===========================================

      do m = 1, 12
        SMB%Refreezing( vi,m) = SMB%Refreezing_year( vi) / 12._dp
        SMB%Runoff(     vi,m) = SMB%Melt( vi,m) + SMB%Rainfall( vi,m) - SMB%Refreezing( vi,m)
        SMB%SMB(        vi,m) = SMB%Snowfall( vi,m) + SMB%Refreezing( vi,m) - SMB%Melt( vi,m)
      end do

      ! === Annual SMB ===
      ! ==================

      ! Sum of monthly fields
      SMB%SMB_year( vi) = sum(SMB%SMB( vi,:))

      ! === Annual ablation ===
      ! =======================

      ! Calculate total melt over this year, to be used for determining next year's albedo
      SMB%MeltPreviousYear( vi) = sum(SMB%Melt( vi,:))

    end do ! vi = mesh%vi1, mesh%vi2

    ! === Final quantities ===
    ! ========================

    ! Convert final SMB from water to ice equivalent
    SMB%SMB(      mesh%vi1:mesh%vi2,:) = SMB%SMB(      mesh%vi1:mesh%vi2,:) * 1000._dp / ice_density
    SMB%SMB_year( mesh%vi1:mesh%vi2  ) = SMB%SMB_year( mesh%vi1:mesh%vi2  ) * 1000._dp / ice_density

    ! === Finalisation ===
    ! ====================

    ! Safety
    call check_for_NaN_dp_1D( SMB%AlbedoSurf      , 'SMB%AlbedoSurf'      )
    call check_for_NaN_dp_2D( SMB%Albedo          , 'SMB%Albedo'          )
    call check_for_NaN_dp_2D( SMB%Melt            , 'SMB%Melt'            )
    call check_for_NaN_dp_2D( SMB%Snowfall        , 'SMB%Snowfall'        )
    call check_for_NaN_dp_2D( SMB%Rainfall        , 'SMB%Rainfall'        )
    call check_for_NaN_dp_2D( SMB%Refreezing      , 'SMB%Refreezing'      )
    call check_for_NaN_dp_2D( SMB%Runoff          , 'SMB%Runoff'          )
    call check_for_NaN_dp_2D( SMB%SMB             , 'SMB%SMB'             )
    call check_for_NaN_dp_2D( SMB%AddedFirn       , 'SMB%AddedFirn'       )
    call check_for_NaN_dp_2D( SMB%FirnDepth       , 'SMB%FirnDepth'       )
    call check_for_NaN_dp_1D( SMB%SMB_year        , 'SMB%SMB_year'        )
    call check_for_NaN_dp_1D( SMB%MeltPreviousYear, 'SMB%MeltPreviousYear')
    !call check_for_NaN_dp_1D( SMB%Albedo_year     , 'SMB%Albedo_year'     ) ! Unused, undefined

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_IMAUITM


  ! The EISMINT SMB parameterisations
  SUBROUTINE EISMINT_SMB( mesh, time, SMB)
    ! Run the IMAU-ITM SMB model. Based on the one from ANICE.
    
    ! NOTE: all the SMB components and the total are in meters of water equivalent
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    REAL(dp),                            INTENT(IN)    :: time
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'EISMINT_SMB'
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
    
    IF     (C%choice_benchmark_experiment == 'EISMINT_1') THEN ! Moving margin, steady state
      ! No changes
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_2') THEN ! Moving margin, 20 kyr
      IF (time < 0._dp) THEN
        ! No changes; first 120 kyr are initialised with EISMINT_1
      ELSE
        E         = 450000._dp + 100000._dp * SIN( 2._dp * pi * time / 20000._dp)
      END IF
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_3') THEN ! Moving margin, 40 kyr
      IF (time < 0._dp) THEN
        ! No changes; first 120 kyr are initialised with EISMINT_1
      ELSE
        E         = 450000._dp + 100000._dp * SIN( 2._dp * pi * time / 40000._dp)
      END IF
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_4') THEN ! Fixed margin, steady state
      M_max       = 0.3_dp       
      E           = 999000._dp
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_5') THEN ! Fixed margin, 20 kyr
      IF (time < 0._dp) THEN
        M_max     = 0.3_dp
        E         = 999000._dp 
      ELSE
        M_max     = 0.3_dp + 0.2_dp * SIN( 2._dp * pi * time / 20000._dp)
        E         = 999000._dp 
      END IF
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_6') THEN ! Fixed margin, 40 kyr
      IF (time < 0._dp) THEN
        M_max     = 0.3_dp
        E         = 999000._dp 
      ELSE
        M_max     = 0.3_dp + 0.2_dp * SIN( 2._dp * pi * time / 40000._dp)
        E         = 999000._dp 
      END IF
    END IF

    DO vi = mesh%vi1, mesh%vi2
      dist = NORM2( mesh%V( vi,:))
      SMB%SMB_year(vi) = MIN( M_max, S_b * (E - dist))
      SMB%SMB(vi,:) = SMB%SMB_year(vi) / 12._dp
    END DO
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
          
  END SUBROUTINE EISMINT_SMB
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

  subroutine initialise_SMB_model_IMAU_ITM( mesh, ice, SMB, region_name)
    ! Allocate memory for the data fields of the SMB model.

    implicit none

    ! In/output variables
    type(type_mesh),      intent(inout) :: mesh
    type(type_ice_model), intent(in)    :: ice
    type(type_SMB_model), intent(inout) :: SMB
    character(len=3),     intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_SMB_model_IMAU_ITM'
    integer                             :: vi
    character(len=256)                  :: SMB_IMAUITM_choice_init_firn

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Data fields
    allocate(  SMB%AlbedoSurf      (mesh%vi1:mesh%vi2    ))
    SMB%AlbedoSurf = 0
    allocate(  SMB%MeltPreviousYear(mesh%vi1:mesh%vi2    ))
    SMB%MeltPreviousYear = 0
    allocate(  SMB%FirnDepth       (mesh%vi1:mesh%vi2, 12))
    SMB%FirnDepth = 0
    allocate(  SMB%Rainfall        (mesh%vi1:mesh%vi2, 12))
    SMB%Rainfall = 0
    allocate(  SMB%Snowfall        (mesh%vi1:mesh%vi2, 12))
    SMB%Snowfall = 0
    allocate(  SMB%AddedFirn       (mesh%vi1:mesh%vi2, 12))
    SMB%AddedFirn = 0
    allocate(  SMB%Melt            (mesh%vi1:mesh%vi2, 12))
    SMB%Melt = 0
    allocate(  SMB%Refreezing      (mesh%vi1:mesh%vi2, 12))
    SMB%Refreezing = 0
    allocate(  SMB%Refreezing_year (mesh%vi1:mesh%vi2    ))
    SMB%Refreezing_year = 0
    allocate(  SMB%Runoff          (mesh%vi1:mesh%vi2, 12))
    SMB%Runoff = 0
    allocate(  SMB%Albedo          (mesh%vi1:mesh%vi2, 12))
    SMB%Albedo = 0
    !allocate(  SMB%Albedo_year     (mesh%vi1:mesh%vi2    ))
    !SMB%Albedo_year = 0
    allocate(  SMB%SMB             (mesh%vi1:mesh%vi2, 12))
    SMB%SMB = 0
    allocate(  SMB%SMB_year        (mesh%vi1:mesh%vi2    ))
    SMB%SMB_year = 0

    ! === Parameters ===
    ! ==================

    ! Initialise regional scalar parameters to specified values
    if     (region_name == 'NAM') then
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_NAM
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_NAM
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_NAM
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_NAM
    elseif (region_name == 'EAS') then
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_EAS
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_EAS
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_EAS
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_EAS
    elseif (region_name == 'GRL') then
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_GRL
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_GRL
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_GRL
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_GRL
    elseif (region_name == 'ANT') then
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_ANT
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_ANT
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_ANT
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_ANT
    end if

    ! === Firn ===
    ! ============

    ! Initialisation choice
    if     (region_name == 'NAM') then
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_NAM
    elseif (region_name == 'EAS') then
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_EAS
    elseif (region_name == 'GRL') then
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_GRL
    elseif (region_name == 'ANT') then
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_ANT
    end if

    select case (SMB_IMAUITM_choice_init_firn)

      case ('uniform')
        ! Initialise with a uniform firn layer over the ice sheet

        do vi = mesh%vi1, mesh%vi2
          if (ice%Hi_a( vi) > 0._dp) then
            SMB%FirnDepth(        vi,:) = C%SMB_IMAUITM_initial_firn_thickness
            SMB%MeltPreviousYear( vi  ) = 0._dp
          else
            SMB%FirnDepth(        vi,:) = 0._dp
            SMB%MeltPreviousYear( vi  ) = 0._dp
          end if
        end do

      case ('restart')
        ! Initialise with the firn layer of a previous run
        call crash('firn from restart not implemented yet!')

      case default
        ! Unknown case
        call crash('unknown SMB_IMAUITM_choice_init_firn "' // &
                    trim(SMB_IMAUITM_choice_init_firn) // '"!')

    end select

    ! === Albedo ===
    ! ==============

    ! Initialise albedo
    do vi = mesh%vi1, mesh%vi2

      ! Background albedo
      if (ice%Hb_a( vi) < 0._dp) then
        SMB%AlbedoSurf( vi) = albedo_water
      else
        SMB%AlbedoSurf( vi) = albedo_soil
      end if

      if (ice%Hi_a( vi) > 0._dp) then
        SMB%AlbedoSurf(  vi) = albedo_snow
      end if

      SMB%Albedo( vi,:) = SMB%AlbedoSurf( vi)

    end do

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_IMAU_ITM

! ===== Remapping =====
! =====================

  subroutine remap_SMB_model( mesh_old, mesh_new, map, SMB)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    type(type_mesh),                     intent(in)    :: mesh_old
    type(type_mesh),                     intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh),      intent(in)    :: map
    type(type_SMB_model),                intent(inout) :: SMB

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'remap_SMB_model'
    integer                                            :: int_dummy

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV
    int_dummy = int_dummy

    ! === Annual SMB ===
    ! ==================

    ! Reallocate annual SMB, needed by all methods
    call reallocate_bounds( SMB%SMB_year, mesh_new%vi1,mesh_new%vi2)

    ! === IMAU-ITM ===
    ! ================

    ! Reallocate IMAU-ITM stuff
    if (C%choice_SMB_model == 'IMAU-ITM') then
      ! Firn depth and melt-during-previous-year must be remapped
      call remap_field_dp_2D( mesh_old, mesh_new, map, SMB%MeltPreviousYear, 'trilin')
      call remap_field_dp_3D( mesh_old, mesh_new, map, SMB%FirnDepth,        'trilin')

      ! Reallocate rather than remap; after a mesh update we'll immediately run the BMB model anyway
      call reallocate_bounds( SMB%Q_TOA           ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%AlbedoSurf      ,mesh_new%vi1,mesh_new%vi2    )
      call reallocate_bounds( SMB%Rainfall        ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%Snowfall        ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%AddedFirn       ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%Melt            ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%Refreezing      ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%Refreezing_year ,mesh_new%vi1,mesh_new%vi2    )
      call reallocate_bounds( SMB%Runoff          ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%Albedo          ,mesh_new%vi1,mesh_new%vi2, 12)
      !call reallocate_bounds( SMB%Albedo_year     ,mesh_new%vi1,mesh_new%vi2    )
      call reallocate_bounds( SMB%SMB             ,mesh_new%vi1,mesh_new%vi2, 12)
    end if

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model

end module SMB_module
