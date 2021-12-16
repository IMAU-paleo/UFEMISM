MODULE SMB_module

  ! All the routines for calculating the surface mass balance from the climate forcing

  ! Import basic functionality
  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, write_to_memory_log, &
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
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  
  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_subclimate_region, type_init_data_fields, &
                                             type_SMB_model, type_remapping
  USE forcing_module,                  ONLY: forcing
  USE mesh_mapping_module,             ONLY: remap_field_dp, remap_field_dp_monthly

  IMPLICIT NONE
  
  REAL(dp), PARAMETER :: albedo_water        = 0.1_dp
  REAL(dp), PARAMETER :: albedo_soil         = 0.2_dp
  REAL(dp), PARAMETER :: albedo_ice          = 0.5_dp
  REAL(dp), PARAMETER :: albedo_snow         = 0.85_dp
  REAL(dp), PARAMETER :: initial_snow_depth  = 0.1_dp
    
CONTAINS

  ! Run the SMB model on the region mesh
  SUBROUTINE run_SMB_model( mesh, ice, climate, time, SMB, mask_noice)
    ! Run the IMAU-ITM SMB model. Old version, exactly as it was in ANICE2.1 (so with the "wrong" refreezing)
    
    ! NOTE: all the SMB components and the total are in meters of water equivalent
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_subclimate_region),        INTENT(IN)    :: climate
    REAL(dp),                            INTENT(IN)    :: time
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_noice
    
    ! Local variables:
    CHARACTER(LEN=64), PARAMETER                       :: routine_name = 'run_SMB_model'
    INTEGER                                            :: n1, n2
    INTEGER                                            :: vi,m
    INTEGER                                            :: mprev
    REAL(dp)                                           :: snowfrac, liquid_water, sup_imp_wat
    REAL(dp)                                           :: R
    
    n1 = par%mem%n
    
    ! Check if we need to apply any special benchmark experiment SMB
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6') THEN
          
        CALL EISMINT_SMB( mesh, time, SMB)
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'Halfar') THEN
        SMB%SMB_year( mesh%v1:mesh%v2  ) = 0._dp
        SMB%SMB(      mesh%v1:mesh%v2,:) = 0._dp
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'Bueler') THEN
        DO vi = mesh%v1, mesh%v2
          SMB%SMB_year( vi  ) = Bueler_solution_MB( C%halfar_solution_H0, C%halfar_solution_R0, C%bueler_solution_lambda, mesh%V(vi,1), mesh%V(vi,2), time)
          SMB%SMB(      vi,:) = SMB%SMB_year( vi) / 12._dp
        END DO
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        SMB%SMB_year( mesh%v1:mesh%v2  ) = 0.3_dp
        SMB%SMB(      mesh%v1:mesh%v2,:) = 0.3_dp / 12._dp
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        ! Similar to EISMINT
        DO vi = mesh%v1, mesh%v2
          R = NORM2(mesh%V(vi,:))
          IF (R < 250000._dp) THEN
            SMB%SMB_year( vi  ) = 0.3_dp
          ELSE
            SMB%SMB_year( vi  ) = MAX(-2._dp, 0.3_dp - (R - 250000._dp) / 200000._dp)
          END IF
          SMB%SMB(      vi,:) = SMB%SMB_year( vi  ) / 12._dp
        END DO ! DO vi = mesh%v1, mesh%v2
        CALL sync
        RETURN
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_SMB_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Calculate SMB components with IMAU_ITM
    ! ======================================
    
    DO vi = mesh%v1, mesh%v2
      
      ! Background albedo
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
        
        ! Firn cannot accumulate on water, silly!
        IF (ice%mask_ocean_a( vi) == 1) SMB%FirnDepth( vi,m) = 0._dp
    
      END DO ! DO m = 1, 12
      
      ! Calculate total SMB for the entire year
      SMB%SMB_year( vi) = SUM(SMB%SMB( vi,:))
      
      ! Calculate total melt over this year, to be used for determining next year's albedo
      SMB%MeltPreviousYear( vi) = SUM(SMB%Melt( vi,:))
      
      ! Calculate yearly mean albedo
      SMB%Albedo_year( vi) = SUM(SMB%Albedo( vi,:)) / 12._dp
      
    END DO
    CALL sync
    
    n2 = par%mem%n
    !CALL write_to_memory_log( routine_name, n1, n2)
          
  END SUBROUTINE run_SMB_model
  
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
    INTEGER                                            :: vi
    
    REAL(dp)                                           :: E               ! Radius of circle where accumulation is M_max
    REAL(dp)                                           :: dist            ! distance to centre of circle
    REAL(dp)                                           :: S_b             ! Gradient of accumulation-rate change with horizontal distance
    REAL(dp)                                           :: M_max           ! Maximum accumulation rate 
    
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

    DO vi = mesh%v1, mesh%v2
      dist = NORM2( mesh%V( vi,:))
      SMB%SMB_year(vi) = MIN( M_max, S_b * (E - dist))
      SMB%SMB(vi,:) = SMB%SMB_year(vi) / 12._dp
    END DO
    CALL sync
          
  END SUBROUTINE EISMINT_SMB
  FUNCTION Bueler_solution_MB( H0, R0, lambda, x, y, t) RESULT(M)
    ! Describes an ice-sheet at time t (in years) conforming to the Bueler solution
    ! with dome thickness H0 and margin radius R0 at t0, with a surface mass balance
    ! determined by lambda.
    
    ! Input variables
    REAL(dp), INTENT(IN) :: H0      ! Ice dome thickness at t=0 [m]
    REAL(dp), INTENT(IN) :: R0      ! Ice margin radius  at t=0 [m]
    REAL(dp), INTENT(IN) :: lambda  ! Mass balance parameter
    REAL(dp), INTENT(IN) :: x       ! x coordinate [m]
    REAL(dp), INTENT(IN) :: y       ! y coordinate [m]
    REAL(dp), INTENT(IN) :: t       ! Time from t0 [years]
    
    ! Result
    REAL(dp)             :: M
    
    ! Local variables
    REAL(dp) :: A_flow, rho, g, n, alpha, beta, Gamma, f1, f2, t0, tp, f3, f4, H
  
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
  
  ! Initialise the SMB model (allocating shared memory)
  SUBROUTINE initialise_SMB_model( mesh, init, SMB, region_name)
    ! Allocate memory for the data fields of the SMB model.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables
    CHARACTER(LEN=64), PARAMETER                       :: routine_name = 'initialise_SMB_model'
    INTEGER                                            :: n1, n2
    INTEGER                                            :: vi
    
    n1 = par%mem%n
    
    IF (par%master) WRITE (0,*) '  Initialising SMB model...'
    
    ! Allocate memory
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%Q_TOA,            SMB%wQ_TOA           )
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
    
    IF (par%master) THEN
      IF     (region_name == 'NAM') THEN
        SMB%C_abl_constant = C%C_abl_constant_NAM
        SMB%C_abl_Ts       = C%C_abl_Ts_NAM
        SMB%C_abl_Q        = C%C_abl_Q_NAM
        SMB%C_refr         = C%C_refr_NAM
      ELSEIF (region_name == 'EAS') THEN
        SMB%C_abl_constant = C%C_abl_constant_EAS
        SMB%C_abl_Ts       = C%C_abl_Ts_EAS
        SMB%C_abl_Q        = C%C_abl_Q_EAS
        SMB%C_refr         = C%C_refr_EAS
      ELSEIF (region_name == 'GRL') THEN
        SMB%C_abl_constant = C%C_abl_constant_GRL
        SMB%C_abl_Ts       = C%C_abl_Ts_GRL
        SMB%C_abl_Q        = C%C_abl_Q_GRL
        SMB%C_refr         = C%C_refr_GRL
      ELSEIF (region_name == 'ANT') THEN
        SMB%C_abl_constant = C%C_abl_constant_ANT
        SMB%C_abl_Ts       = C%C_abl_Ts_ANT
        SMB%C_abl_Q        = C%C_abl_Q_ANT
        SMB%C_refr         = C%C_refr_ANT
      END IF
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Initialise albedo to background albedo
    DO vi = mesh%v1, mesh%v2
      
      ! Background albedo
      IF (init%Hb(vi) < 0._dp) THEN
        SMB%AlbedoSurf(vi) = albedo_water
      ELSE
        SMB%AlbedoSurf(vi) = albedo_soil
      END IF
      
      IF (init%Hi(vi) > 0._dp) THEN
        SMB%AlbedoSurf(vi) = albedo_snow
        SMB%FirnDepth(vi,:) = initial_snow_depth   
      END IF
      
      SMB%Albedo(      vi,:) = SMB%AlbedoSurf(vi)
      SMB%Albedo_year( vi  ) = SMB%AlbedoSurf(vi)
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! If we're doing a restart, initialise with that
    IF (C%is_restart) THEN
    
      SMB%FirnDepth(        mesh%v1:mesh%v2,:) = init%FirnDepth(        mesh%v1:mesh%v2,:)
      SMB%MeltPreviousYear( mesh%v1:mesh%v2  ) = init%MeltPreviousYear( mesh%v1:mesh%v2  )
      
      DO vi = mesh%v1, mesh%v2
        
        ! Background albedo
        IF (init%Hb(vi) < 0._dp) THEN
          SMB%AlbedoSurf(vi) = albedo_water
        ELSE
          SMB%AlbedoSurf(vi) = albedo_soil
        END IF
        
        IF (SMB%FirnDepth( vi,12) > 0._dp) THEN
          SMB%AlbedoSurf(vi)  = albedo_snow
        END IF
        
        SMB%Albedo(      vi,:) = SMB%AlbedoSurf(vi)
        SMB%Albedo_year( vi  ) = SMB%AlbedoSurf(vi)
        
      END DO
      CALL sync
      
    END IF ! IF (C%is_restart) THEN
    
    n2 = par%mem%n
    CALL write_to_memory_log( routine_name, n1, n2)
  
  END SUBROUTINE initialise_SMB_model  
  SUBROUTINE remap_SMB_model( mesh_old, mesh_new, map, SMB)
    ! Remap or reallocate all the data fields
  
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    
    ! Local variables:
    INTEGER                                            :: int_dummy
    
    int_dummy = mesh_old%nV
    int_dummy = map%trilin%vi( 1,1)
    
    ! Firn depth and melt-during-previous-year must be remapped
    CALL remap_field_dp(         mesh_old, mesh_new, map, SMB%MeltPreviousYear, SMB%wMeltPreviousYear, 'trilin')
    CALL remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%FirnDepth,        SMB%wFirnDepth,        'trilin')
        
    ! Reallocate rather than remap; after a mesh update we'll immediately run the BMB model anyway
    CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%Q_TOA,            SMB%wQ_TOA           )
    CALL reallocate_shared_dp_1D( mesh_new%nV,     SMB%AlbedoSurf,       SMB%wAlbedoSurf      )
   !CALL reallocate_shared_dp_1D( mesh_new%nV,     SMB%MeltPreviousYear, SMB%wMeltPreviousYear)
   !CALL reallocate_shared_dp_2D( mesh_new%nV, 12, SMB%FirnDepth,        SMB%wFirnDepth       )
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
    CALL reallocate_shared_dp_1D( mesh_new%nV,     SMB%SMB_year,         SMB%wSMB_year        )
    
  END SUBROUTINE remap_SMB_model

END MODULE SMB_module
