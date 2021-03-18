MODULE SMB_module

  USE mpi
  USE configuration_module,          ONLY: dp, C
  USE parallel_module,               ONLY: par, sync, &
                                           allocate_shared_int_0D, allocate_shared_dp_0D, &
                                           allocate_shared_int_1D, allocate_shared_dp_1D, &
                                           allocate_shared_int_2D, allocate_shared_dp_2D, &
                                           allocate_shared_int_3D, allocate_shared_dp_3D, &
                                           deallocate_shared
  USE data_types_module,             ONLY: type_mesh, type_ice_model, type_climate_model, type_init_data_fields, type_SMB_model, &
                                           type_remapping
  USE netcdf_module,                 ONLY: debug, write_to_debug_file
  USE forcing_module,                ONLY: forcing
  USE parameters_module,             ONLY: T0, L_fusion
  USE mesh_mapping_module,           ONLY: Remap_field_dp, Remap_field_dp_3D, Remap_field_dp_monthly, &
                                           Reallocate_field_dp, Reallocate_field_dp_3D, Reallocate_field_int

  IMPLICIT NONE
  
  REAL(dp), PARAMETER :: albedo_water        = 0.1_dp
  REAL(dp), PARAMETER :: albedo_soil         = 0.2_dp
  REAL(dp), PARAMETER :: albedo_ice          = 0.5_dp
  REAL(dp), PARAMETER :: albedo_snow         = 0.85_dp
  REAL(dp), PARAMETER :: initial_snow_depth  = 0.1_dp
    
CONTAINS

  ! Run the SMB model on the region mesh
  SUBROUTINE run_SMB_model( mesh, ice, climate, time, SMB)
    ! Run the IMAU-ITM SMB model. Based on the one from ANICE.
    
    ! NOTE: all the SMB components and the total are in meters of water equivalent
    
    USE parameters_module, ONLY: sec_per_year
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_climate_model),            INTENT(IN)    :: climate
    REAL(dp),                            INTENT(IN)    :: time
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: vi, m, ilat_l, ilat_u
    REAL(dp)                                           :: R, wt0, wt1, wlat_l, wlat_u
    INTEGER                                            :: mprev
    REAL(dp)                                           :: snowfrac, liquid_water, sup_imp_wat
    
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
    
    ! Update Q_TOA
    ! ============
    
    ! Calculate time interpolation weights
    wt0 = ((forcing%ins_t1 * 1000._dp) - time) / ((forcing%ins_t1 - forcing%ins_t0) * 1000._dp)
    wt1 = 1._dp - wt0
        
    ! Interpolate on the mesh
    DO vi = mesh%v1, mesh%v2
     
      ilat_l = FLOOR(mesh%lat(vi) + 91)
      ilat_u = ilat_l + 1
      
      wlat_l = forcing%ins_lat(ilat_u) - mesh%lat(vi)
      wlat_u = 1._dp - wlat_l
      
      DO m = 1, 12
        SMB%Q_TOA(vi,m) = (wt0 * wlat_l * forcing%ins_Q_TOA0(ilat_l,m)) + &
                          (wt0 * wlat_u * forcing%ins_Q_TOA0(ilat_u,m)) + &
                          (wt1 * wlat_l * forcing%ins_Q_TOA1(ilat_l,m)) + &
                          (wt1 * wlat_u * forcing%ins_Q_TOA1(ilat_u,m))
      END DO    
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! Calculate SMB components with IMAU_ITM
    ! ======================================
    
    DO vi = mesh%v1, mesh%v2
      
      ! Background albedo
      SMB%AlbedoSurf(vi) = albedo_soil
      IF (ice%mask_ocean(vi)==1) SMB%AlbedoSurf(vi) = albedo_water
      IF (ice%mask_ice(  vi)==1) SMB%AlbedoSurf(vi) = albedo_ice
    
      DO m = 1, 12  ! Month loop
        
        mprev = m - 1
        IF (mprev==0) mprev = 12
        
        SMB%Albedo(vi,m) = MIN(albedo_snow, MAX( SMB%AlbedoSurf(vi), albedo_snow - (albedo_snow - SMB%AlbedoSurf(vi))  * &
                            EXP(-15._dp * SMB%FirnDepth(vi,mprev)) - 0.015_dp * SMB%MeltPreviousYear(vi)))
        IF (ice%mask_ocean(vi)==1) SMB%Albedo(vi,m) = albedo_water
               
        ! Determine albation as function af surface temperature and albedo/insolation
        ! according to Bintanja et al. (2002) 
    
        SMB%Melt(vi,m) = MAX(0._dp, ( C%C_abl_Ts         * (climate%applied%T2m(vi,m) - T0) + &
                                      C%C_abl_Q          * (1.0_dp - SMB%Albedo(vi,m)) * SMB%Q_TOA(vi,m) - &
                                      C%C_abl_constant)  * sec_per_year / (L_fusion * 1000._dp * 12._dp))
                
        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999),
        ! liquid water content (rain and melt water) and snowdepth
    
        ! NOTE: commented version is the old ANICE version, supposedly based on physics (which we cant check), but 
        !       the new version was tuned to RACMO output and produced significantly better snow fractions...
    
  !      snowfrac = MAX(0._dp, MIN(1._dp, 0.5_dp   * (1 - ATAN((climate%applied%T2m(vi,m) - T0) / 3.5_dp)  / 1.25664_dp)))
        snowfrac = MAX(0._dp, MIN(1._dp, 0.725_dp * (1 - ATAN((climate%applied%T2m(vi,m) - T0) / 5.95_dp) / 1.8566_dp)))
    
        SMB%Snowfall(vi,m) = climate%applied%Precip(vi,m) *          snowfrac
        SMB%Rainfall(vi,m) = climate%applied%Precip(vi,m) * (1._dp - snowfrac)
    
        ! Refreezing, according to Janssens & Huybrechts, 2000)
        ! The refreezing (=effective retention) is the minimum value of the amount of super imposed 
        ! water and the available liquid water, with a maximum value of the total precipitation.
        ! (see also Huybrechts & de Wolde, 1999)
    
        ! Add this month's snow accumulation to next month's initial snow depth.
        SMB%AddedFirn(vi,m) = SMB%Snowfall(vi,m) - SMB%Melt(vi,m)
        SMB%FirnDepth(vi,m) = MIN(10._dp, MAX(0._dp, SMB%FirnDepth(vi,mprev) + SMB%AddedFirn(vi,m) ))
    
      END DO  ! DO m = 1, 12
    
      ! Calculate refrezzing for the whole year, divide equally over the 12 months, then calculate resulting runoff and SMB.
      ! This resolves the problem with refreezing, where liquid water is mostly available in summer
      ! but "refreezing potential" mostly in winter, and there is no proper meltwater retention.
      
      sup_imp_wat  = C%C_refr * MAX(0._dp, T0 - SUM(climate%applied%T2m(vi,:))/12._dp)
      liquid_water = SUM(SMB%Rainfall(vi,:)) + SUM(SMB%Melt(vi,:))
      
      SMB%Refreezing_year(vi) = MIN( MIN( sup_imp_wat, liquid_water), SUM(climate%applied%Precip(vi,:)))
      IF (ice%mask_ice(vi)==0) SMB%Refreezing_year(vi) = 0._dp
  
      DO m = 1, 12
        SMB%Refreezing(vi,m) = SMB%Refreezing_year(vi) / 12._dp
        SMB%Runoff(vi,m) = SMB%Melt(vi,m) + SMB%Rainfall(vi,m) - SMB%Refreezing(vi,m)
        SMB%SMB(vi,m) = SMB%Snowfall(vi,m) + SMB%Refreezing(vi,m) - SMB%Melt(vi,m)
      END DO
      
      SMB%SMB_year(vi) = SUM(SMB%SMB(vi,:))
      
      ! Calculate total melt over this year, to be used for determining next year's albedo
      SMB%MeltPreviousYear(vi) = SUM(SMB%Melt(vi,:))
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
          
  END SUBROUTINE run_SMB_model
  
  ! The EISMINT SMB parameterisations
  SUBROUTINE EISMINT_SMB( mesh, time, SMB)
    ! Run the IMAU-ITM SMB model. Based on the one from ANICE.
    
    ! NOTE: all the SMB components and the total are in meters of water equivalent
    
    USE parameters_module, ONLY: pi
    
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
    
    USE parameters_module, ONLY: sec_per_year
    
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
  SUBROUTINE initialise_SMB_model( mesh, init, filetype_init, SMB)
    ! Allocate memory for the data fields of the SMB model.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
    CHARACTER(LEN=4),                    INTENT(IN)    :: filetype_init
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    
    ! Local variables
    INTEGER                                            :: vi
    
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
    CALL allocate_shared_dp_2D( mesh%nV, 12, SMB%SMB,              SMB%wSMB             )
    CALL allocate_shared_dp_1D( mesh%nV,     SMB%SMB_year,         SMB%wSMB_year        )
      
    SMB%Q_TOA(            mesh%v1:mesh%v2, :) = 0._dp
    SMB%AlbedoSurf(       mesh%v1:mesh%v2   ) = 0._dp
    SMB%MeltPreviousYear( mesh%v1:mesh%v2   ) = 0._dp
    SMB%FirnDepth(        mesh%v1:mesh%v2, :) = 0._dp
    SMB%Rainfall(         mesh%v1:mesh%v2, :) = 0._dp
    SMB%Snowfall(         mesh%v1:mesh%v2, :) = 0._dp
    SMB%AddedFirn(        mesh%v1:mesh%v2, :) = 0._dp
    SMB%Melt(             mesh%v1:mesh%v2, :) = 0._dp
    SMB%Refreezing(       mesh%v1:mesh%v2, :) = 0._dp
    SMB%Refreezing_year(  mesh%v1:mesh%v2   ) = 0._dp
    SMB%Runoff(           mesh%v1:mesh%v2, :) = 0._dp
    SMB%Albedo(           mesh%v1:mesh%v2, :) = 0._dp
    SMB%SMB(              mesh%v1:mesh%v2, :) = 0._dp
    SMB%SMB_year(         mesh%v1:mesh%v2   ) = 0._dp
    
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
      
      SMB%Albedo(vi,:) = SMB%AlbedoSurf(vi)
      
    END DO
    
    ! If the model is initialised from a previous run, use Albedo, FirnDepth and MeltPreviousYear from init
    IF (filetype_init == 'mesh') THEN
      SMB%Albedo(           mesh%v1:mesh%v2,:) = init%Albedo(           mesh%v1:mesh%v2,:)
      SMB%FirnDepth(        mesh%v1:mesh%v2,:) = init%FirnDepth(        mesh%v1:mesh%v2,:)
      SMB%MeltPreviousYear( mesh%v1:mesh%v2  ) = init%MeltPreviousYear( mesh%v1:mesh%v2  )
    END IF
  
  END SUBROUTINE initialise_SMB_model  
  SUBROUTINE remap_SMB_model( mesh_old, mesh_new, map, SMB)
    ! Remap or reallocate all the data fields
  
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
        
    ! All SMB fields are remapped, since they are not updated in every timestep.
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%Q_TOA,               SMB%wQ_TOA,               'trilin'        )
    CALL Remap_field_dp(         mesh_old, mesh_new, map, SMB%AlbedoSurf,          SMB%wAlbedoSurf,          'trilin'        )
    CALL Remap_field_dp(         mesh_old, mesh_new, map, SMB%MeltPreviousYear,    SMB%wMeltPreviousYear,    'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%FirnDepth,           SMB%wFirnDepth,           'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%Rainfall,            SMB%wRainfall,            'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%Snowfall,            SMB%wSnowfall,            'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%AddedFirn,           SMB%wAddedFirn,           'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%Melt,                SMB%wMelt,                'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%Refreezing,          SMB%wRefreezing,          'trilin'        )
    CALL Remap_field_dp(         mesh_old, mesh_new, map, SMB%Refreezing_year,     SMB%wRefreezing_year,     'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%Runoff,              SMB%wRunoff,              'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%Albedo,              SMB%wAlbedo,              'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, SMB%SMB,                 SMB%wSMB,                 'trilin'        )
    CALL Remap_field_dp(         mesh_old, mesh_new, map, SMB%SMB_year,            SMB%wSMB_year,            'trilin'        )
    
  END SUBROUTINE remap_SMB_model

END MODULE SMB_module
