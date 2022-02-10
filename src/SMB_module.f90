MODULE SMB_module

  ! All the routines for calculating the surface mass balance from the climate forcing

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parameters_module
  USE petsc_module,                    ONLY: perr
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
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             transpose_dp_2D, transpose_dp_3D
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, inquire_restart_file_SMB, &
                                             read_restart_file_SMB
  
  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_subclimate_region, &
                                             type_SMB_model, type_remapping_mesh_mesh, &
                                             type_climate_matrix_regional, &
                                             type_climate_snapshot_regional, type_direct_SMB_forcing_regional, &
                                             type_restart_data
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
        
      ELSEIF (C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
        SMB%SMB_year( mesh%vi1:mesh%vi2  ) = 0._dp
        SMB%SMB(      mesh%vi1:mesh%vi2,:) = 0._dp
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'Bueler') THEN
        DO vi = mesh%vi1, mesh%vi2
          SMB%SMB_year( vi  ) = Bueler_solution_MB( mesh%V(vi,1), mesh%V(vi,2), time)
          SMB%SMB(      vi,:) = SMB%SMB_year( vi) / 12._dp
        END DO
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        SMB%SMB_year( mesh%vi1:mesh%vi2  ) = 0.3_dp
        SMB%SMB(      mesh%vi1:mesh%vi2,:) = 0.3_dp / 12._dp
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
        RETURN
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_SMB_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Calculate SMB components with IMAU_ITM
    ! ======================================
    
    DO vi = mesh%vi1, mesh%vi2
      
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

    DO vi = mesh%vi1, mesh%vi2
      dist = NORM2( mesh%V( vi,:))
      SMB%SMB_year(vi) = MIN( M_max, S_b * (E - dist))
      SMB%SMB(vi,:) = SMB%SMB_year(vi) / 12._dp
    END DO
    CALL sync
          
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
  
  ! Initialise the SMB model (allocating shared memory)
  SUBROUTINE initialise_SMB_model( mesh, ice, SMB, region_name)
    ! Allocate memory for the data fields of the SMB model.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables
    CHARACTER(LEN=64), PARAMETER                       :: routine_name = 'initialise_SMB_model'
    INTEGER                                            :: n1, n2
    INTEGER                                            :: vi
    
    n1 = par%mem%n
    
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

      CALL initialise_SMB_model_IMAU_ITM( mesh, ice, SMB, region_name)

    ELSE
      IF (par%master) WRITE(0,*) 'initialise_SMB_model - ERROR: unknown choice_SMB_model "', TRIM(C%choice_SMB_model), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    n2 = par%mem%n
    CALL write_to_memory_log( routine_name, n1, n2)
  
  END SUBROUTINE initialise_SMB_model  
  SUBROUTINE remap_SMB_model( mesh_old, mesh_new, map, SMB)
    ! Remap or reallocate all the data fields
  
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    
    ! Local variables:
    INTEGER                                            :: int_dummy
    
    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV
    int_dummy = int_dummy
    
    ! Firn depth and melt-during-previous-year must be remapped
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, SMB%MeltPreviousYear, SMB%wMeltPreviousYear, 'trilin')
    CALL remap_field_dp_3D( mesh_old, mesh_new, map, SMB%FirnDepth,        SMB%wFirnDepth,        'trilin')
        
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

  !=============================
  !=============================

  ! == The main routines that should be called from the main ice model/program
  ! ==========================================================================

  SUBROUTINE run_SMB_model_port( mesh, ice, climate_matrix, time, SMB, mask_noice)
    ! Run the selected SMB model.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_matrix_regional),   INTENT(IN)    :: climate_matrix
    REAL(dp),                             INTENT(IN)    :: time
    TYPE(type_SMB_model),                 INTENT(INOUT) :: SMB
    INTEGER,  DIMENSION(:    ),           INTENT(IN)    :: mask_noice

    ! Local variables:
    INTEGER                                             :: vi

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
      IF (par%master) WRITE(0,*) 'run_SMB_model - ERROR: unknown choice_SMB_model "', TRIM(C%choice_SMB_model), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

  END SUBROUTINE run_SMB_model_port

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

    IF     (C%choice_idealised_SMB == 'EISMINT1_A' .OR. &
            C%choice_idealised_SMB == 'EISMINT1_B' .OR. &
            C%choice_idealised_SMB == 'EISMINT1_C' .OR. &
            C%choice_idealised_SMB == 'EISMINT1_D' .OR. &
            C%choice_idealised_SMB == 'EISMINT1_E' .OR. &
            C%choice_idealised_SMB == 'EISMINT1_F') THEN
      CALL run_SMB_model_idealised_EISMINT1( mesh, SMB, time, mask_noice)
    ELSEIF (C%choice_idealised_SMB == 'Bueler') THEN
      CALL run_SMB_model_idealised_Bueler( mesh, SMB, time, mask_noice)
    ELSE
      IF (par%master) WRITE(0,*) 'run_SMB_model_idealised - ERROR: unknown choice_idealised_SMB "', TRIM(C%choice_idealised_SMB), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

  END SUBROUTINE run_SMB_model_idealised
  SUBROUTINE run_SMB_model_idealised_EISMINT1( mesh, SMB, time, mask_noice)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    REAL(dp),                            INTENT(IN)    :: time
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_noice

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

    END SUBROUTINE run_SMB_model_idealised_EISMINT1
    SUBROUTINE run_SMB_model_idealised_Bueler( mesh, SMB, time, mask_noice)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    REAL(dp),                            INTENT(IN)    :: time
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_noice

    ! Local variables:
    INTEGER                                            :: vi

    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice( vi) == 0) THEN
        SMB%SMB_year( vi) = Bueler_solution_MB( mesh%V(vi,1), mesh%V(vi,2), time)
      ELSE
        SMB%SMB_year( vi) = 0._dp
      END IF
    END DO
    CALL sync

  END SUBROUTINE run_SMB_model_idealised_Bueler

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
    INTEGER                                             :: vi, m
    INTEGER                                             :: mprev
    REAL(dp)                                            :: snowfrac, liquid_water, sup_imp_wat

    ! Make sure this routine is called correctly
    IF (.NOT. C%choice_SMB_model == 'IMAU-ITM') THEN
      IF (par%master) WRITE(0,*) ' run_IMAUITM - ERROR: should only be called when choice_SMB_model == "IMAU-ITM"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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

        SMB%Melt( vi,m) = MAX(0._dp, ( SMB%C_abl_Ts         * (climate%T2m( vi,m) - T0) + &
                                       SMB%C_abl_Q          * (1.0_dp - SMB%Albedo( vi,m)) * climate%Q_TOA( vi,m) - &
                                       SMB%C_abl_constant)  * sec_per_year / (L_fusion * 1000._dp * 12._dp))

        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999),
        ! liquid water content (rain and melt water) and snowdepth

        ! NOTE: commented version is the old ANICE version, supposedly based on "physics" (which we cant check), but
        !       the new version was tuned to RACMO output and produced significantly better snow fractions...

  !      snowfrac = MAX(0._dp, MIN(1._dp, 0.5_dp   * (1 - ATAN((climate%T2m(vi,m) - T0) / 3.5_dp)  / 1.25664_dp)))
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
    CALL check_for_NaN_dp_1D( SMB%AlbedoSurf      , 'SMB%AlbedoSurf'      , 'run_IMAUITM')
    CALL check_for_NaN_dp_2D( SMB%Albedo          , 'SMB%Albedo'          , 'run_IMAUITM')
    CALL check_for_NaN_dp_2D( SMB%Melt            , 'SMB%Melt'            , 'run_IMAUITM')
    CALL check_for_NaN_dp_2D( SMB%Snowfall        , 'SMB%Snowfall'        , 'run_IMAUITM')
    CALL check_for_NaN_dp_2D( SMB%Rainfall        , 'SMB%Rainfall'        , 'run_IMAUITM')
    CALL check_for_NaN_dp_2D( SMB%Refreezing      , 'SMB%Refreezing'      , 'run_IMAUITM')
    CALL check_for_NaN_dp_2D( SMB%Runoff          , 'SMB%Runoff'          , 'run_IMAUITM')
    CALL check_for_NaN_dp_2D( SMB%SMB             , 'SMB%SMB'             , 'run_IMAUITM')
    CALL check_for_NaN_dp_2D( SMB%AddedFirn       , 'SMB%AddedFirn'       , 'run_IMAUITM')
    CALL check_for_NaN_dp_2D( SMB%FirnDepth       , 'SMB%FirnDepth'       , 'run_IMAUITM')
    CALL check_for_NaN_dp_1D( SMB%SMB_year        , 'SMB%SMB_year'        , 'run_IMAUITM')
    CALL check_for_NaN_dp_1D( SMB%MeltPreviousYear, 'SMB%MeltPreviousYear', 'run_IMAUITM')
    CALL check_for_NaN_dp_1D( SMB%Albedo_year     , 'SMB%Albedo_year'     , 'run_IMAUITM')

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
    INTEGER                                             :: vi, m
    INTEGER                                             :: mprev
    REAL(dp)                                            :: snowfrac, liquid_water, sup_imp_wat

    ! Make sure this routine is called correctly
    IF (.NOT. C%choice_SMB_model == 'IMAU-ITM_wrongrefreezing') THEN
      IF (par%master) WRITE(0,*) ' run_IMAUITM_wrongrefreezing - ERROR: should only be called when choice_SMB_model == "IMAU-ITM_wrongrefreezing"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
    CALL check_for_NaN_dp_1D( SMB%AlbedoSurf      , 'SMB%AlbedoSurf'      , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_2D( SMB%Albedo          , 'SMB%Albedo'          , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_2D( SMB%Melt            , 'SMB%Melt'            , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_2D( SMB%Snowfall        , 'SMB%Snowfall'        , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_2D( SMB%Rainfall        , 'SMB%Rainfall'        , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_2D( SMB%Refreezing      , 'SMB%Refreezing'      , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_2D( SMB%Runoff          , 'SMB%Runoff'          , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_2D( SMB%SMB             , 'SMB%SMB'             , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_2D( SMB%AddedFirn       , 'SMB%AddedFirn'       , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_2D( SMB%FirnDepth       , 'SMB%FirnDepth'       , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_1D( SMB%SMB_year        , 'SMB%SMB_year'        , 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_1D( SMB%MeltPreviousYear, 'SMB%MeltPreviousYear', 'run_IMAUITM_wrongrefreezing')
    CALL check_for_NaN_dp_1D( SMB%Albedo_year     , 'SMB%Albedo_year'     , 'run_IMAUITM_wrongrefreezing')

  END SUBROUTINE run_SMB_model_IMAUITM_wrongrefreezing
  SUBROUTINE initialise_SMB_model_IMAU_ITM( mesh, ice, SMB, region_name)
    ! Allocate memory for the data fields of the SMB model.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables
    INTEGER                                            :: vi
    CHARACTER(LEN=256)                                 :: SMB_IMAUITM_choice_init_firn

    IF (par%master) WRITE (0,*) '  Initialising the IMAU-ITM SMB model...'

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

      CALL initialise_IMAU_ITM_firn_restart( mesh, SMB, region_name)

    ELSE
      IF (par%master) WRITE(0,*) 'initialise_SMB_model_IMAU_ITM - ERROR: unknown SMB_IMAUITM_choice_init_firn "', TRIM(SMB_IMAUITM_choice_init_firn), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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

  END SUBROUTINE initialise_SMB_model_IMAU_ITM
  SUBROUTINE initialise_IMAU_ITM_firn_restart( mesh, SMB, region_name)
    ! If this is a restarted run, read the firn depth and meltpreviousyear data from the restart file

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),      INTENT(INOUT) :: mesh
    TYPE(type_SMB_model), INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),     INTENT(IN)    :: region_name

    ! Local variables
    INTEGER                             :: i,j,n
    CHARACTER(LEN=256)                  :: filename_restart
    REAL(dp)                            :: time_to_restart_from
    TYPE(type_restart_data)             :: restart

    ! Assume that SMB and geometry are read from the same restart file
    IF     (region_name == 'NAM') THEN
      filename_restart     = C%filename_refgeo_init_NAM
      time_to_restart_from = C%time_to_restart_from_NAM
    ELSEIF (region_name == 'EAS') THEN
      filename_restart     = C%filename_refgeo_init_EAS
      time_to_restart_from = C%time_to_restart_from_EAS
    ELSEIF (region_name == 'GR:') THEN
      filename_restart     = C%filename_refgeo_init_GRL
      time_to_restart_from = C%time_to_restart_from_GRL
    ELSEIF (region_name == 'ANT') THEN
      filename_restart     = C%filename_refgeo_init_ANT
      time_to_restart_from = C%time_to_restart_from_ANT
    END IF

    ! Inquire if all the required fields are present in the specified NetCDF file,
    ! and determine the dimensions of the memory to be allocated.
    CALL allocate_shared_int_0D( restart%grid%nx, restart%grid%wnx)
    CALL allocate_shared_int_0D( restart%grid%ny, restart%grid%wny)
    CALL allocate_shared_int_0D( restart%nt,      restart%wnt)
    IF (par%master) THEN
      restart%netcdf%filename = filename_restart
      CALL inquire_restart_file_SMB( restart)
    END IF
    CALL sync

    ! Assign range to each processor
    CALL partition_list( restart%grid%nx, par%i, par%n, restart%grid%i1, restart%grid%i2)
    CALL partition_list( restart%grid%ny, par%i, par%n, restart%grid%j1, restart%grid%j2)

    ! Allocate memory for raw data
    CALL allocate_shared_dp_1D( restart%grid%nx, restart%grid%x,    restart%grid%wx)
    CALL allocate_shared_dp_1D( restart%grid%ny, restart%grid%y,    restart%grid%wy)
    CALL allocate_shared_dp_1D( restart%nt,      restart%time,      restart%wtime  )

    CALL allocate_shared_dp_3D( restart%grid%nx, restart%grid%ny, 12, restart%FirnDepth,        restart%wFirnDepth       )
    CALL allocate_shared_dp_2D( restart%grid%nx, restart%grid%ny,     restart%MeltPreviousYear, restart%wMeltPreviousYear)

    ! Read data from input file
    IF (par%master) CALL read_restart_file_SMB( restart, time_to_restart_from)
    CALL sync

    ! Safety
    ! CALL check_for_NaN_dp_3D( restart%FirnDepth,        'restart%FirnDepth',        'initialise_IMAU_ITM_firn_restart')
    ! CALL check_for_NaN_dp_2D( restart%MeltPreviousYear, 'restart%MeltPreviousYear', 'initialise_IMAU_ITM_firn_restart')

    ! ! Since we want data represented as [j,i] internally, transpose the data we just read.
    ! CALL transpose_dp_3D( restart%FirnDepth,        restart%wFirnDepth       )
    ! CALL transpose_dp_2D( restart%MeltPreviousYear, restart%wMeltPreviousYear)

    ! Fill in secondary grid parameters
    CALL allocate_shared_dp_0D( restart%grid%dx,   restart%grid%wdx  )
    CALL allocate_shared_dp_0D( restart%grid%xmin, restart%grid%wxmin)
    CALL allocate_shared_dp_0D( restart%grid%xmax, restart%grid%wxmax)
    CALL allocate_shared_dp_0D( restart%grid%ymin, restart%grid%wymin)
    CALL allocate_shared_dp_0D( restart%grid%ymax, restart%grid%wymax)
    IF (par%master) THEN
      restart%grid%dx   = restart%grid%x( 2) - restart%grid%x( 1)
      restart%grid%xmin = restart%grid%x( 1)
      restart%grid%xmax = restart%grid%x( restart%grid%nx)
      restart%grid%ymin = restart%grid%y( 1)
      restart%grid%ymax = restart%grid%y( restart%grid%ny)
    END IF
    CALL sync

    ! Set up grid-to-vector translation tables
    CALL allocate_shared_int_0D( restart%grid%n, restart%grid%wn)
    IF (par%master) restart%grid%n  = restart%grid%nx * restart%grid%ny
    CALL sync
    CALL allocate_shared_int_2D( restart%grid%nx, restart%grid%ny, restart%grid%ij2n, restart%grid%wij2n)
    CALL allocate_shared_int_2D( restart%grid%n , 2              , restart%grid%n2ij, restart%grid%wn2ij)
    IF (par%master) THEN
      n = 0
      DO i = 1, restart%grid%nx
        IF (MOD(i,2) == 1) THEN
          DO j = 1, restart%grid%ny
            n = n+1
            restart%grid%ij2n( i,j) = n
            restart%grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = restart%grid%ny, 1, -1
            n = n+1
            restart%grid%ij2n( i,j) = n
            restart%grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO
    END IF
    CALL sync

    ! Map (transposed) raw data to the model mesh
    CALL calc_remapping_operator_grid2mesh( restart%grid, mesh)
    CALL map_grid2mesh_3D(restart%grid, mesh, restart%FirnDepth,        SMB%FirnDepth)
    CALL map_grid2mesh_2D(restart%grid, mesh, restart%MeltPreviousYear, SMB%MeltPreviousYear)
    CALL deallocate_remapping_operators_grid2mesh( restart%grid)

    ! Deallocate raw data
    CALL deallocate_shared( restart%grid%wnx         )
    CALL deallocate_shared( restart%grid%wny         )
    CALL deallocate_shared( restart%grid%wx          )
    CALL deallocate_shared( restart%grid%wy          )
    CALL deallocate_shared( restart%grid%wdx         )
    CALL deallocate_shared( restart%grid%wxmin       )
    CALL deallocate_shared( restart%grid%wxmax       )
    CALL deallocate_shared( restart%grid%wymin       )
    CALL deallocate_shared( restart%grid%wymax       )
    CALL deallocate_shared( restart%grid%wn          )
    CALL deallocate_shared( restart%grid%wij2n       )
    CALL deallocate_shared( restart%grid%wn2ij       )
    CALL deallocate_shared( restart%wnt              )
    CALL deallocate_shared( restart%wtime            )
    CALL deallocate_shared( restart%wFirnDepth       )
    CALL deallocate_shared( restart%wMeltPreviousYear)

  END SUBROUTINE initialise_IMAU_ITM_firn_restart

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

    ! Local variables
    REAL(dp)                                           :: wt0, wt1
    INTEGER                                            :: vi

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

  END SUBROUTINE run_SMB_model_direct

END MODULE SMB_module
