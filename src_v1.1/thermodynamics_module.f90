MODULE thermodynamics_module

  USE mpi
  USE configuration_module,          ONLY: dp, C
  USE parallel_module,               ONLY: par, sync, &
                                           allocate_shared_int_0D, allocate_shared_dp_0D, &
                                           allocate_shared_int_1D, allocate_shared_dp_1D, &
                                           allocate_shared_int_2D, allocate_shared_dp_2D, &
                                           allocate_shared_int_3D, allocate_shared_dp_3D, &
                                           deallocate_shared
  USE data_types_module,             ONLY: type_mesh, type_ice_model, type_climate_model, type_SMB_model, type_init_data_fields
  USE netcdf_module,                 ONLY: debug, write_to_debug_file
  USE general_ice_model_data_module, ONLY: ice_physical_properties
  USE ice_dynamics_module,           ONLY: solve_SIA_3D    
  USE parameters_module,             ONLY: ice_density, grav, SMT, L_fusion, T0
  USE zeta_module,                   ONLY: calculate_zeta_derivatives, p_zeta
  USE mesh_derivatives_module,       ONLY: get_upwind_derivative_vertex_3D, apply_Neumann_boundary_3D
  
  IMPLICIT NONE
  
CONTAINS
   
  SUBROUTINE update_ice_temperature( mesh, ice, climate, SMB)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_climate_model),            INTENT(IN)    :: climate
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB

    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: vi, k
    REAL(dp)                                           :: dTi_dx, dTi_dy, internal_heating, f1, f2, f3
    REAL(dp), DIMENSION(2:C%NZ)                        :: alpha
    REAL(dp), DIMENSION(C%NZ)                          :: beta
    REAL(dp), DIMENSION(C%NZ-1)                        :: gamma
    REAL(dp), DIMENSION(C%NZ)                          :: delta
    REAL(dp), DIMENSION(:,:), POINTER                  :: Ti_new
    INTEGER                                            :: wTi_new
    INTEGER                                            :: n_unstable
    
    ! Special cases for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6') THEN
        ! Thermodynamics are included in these experiments
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
              C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
              C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler') THEN
        ! Thermodynamics are not included in these experiments
        RETURN
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in update_ice_temperature!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Allocate and initialise temporary shared memory for the new temperature field
    CALL allocate_shared_dp_2D( mesh%nV, C%nZ, Ti_new, wTi_new)
    Ti_new(mesh%v1:mesh%v2,:) = 0._dp
    CALL sync
    
    ! Calculate the 3D ice velocities and zeta derivatives required for solving the heat equation
    CALL solve_SIA_3D(               mesh, ice)
    CALL bottom_frictional_heating(  mesh, ice)
    CALL calculate_zeta_derivatives( mesh, ice)
    
    ! Set ice surface temperature equal to annual mean 2m air temperature
    DO vi = mesh%v1, mesh%v2
      ice%Ti(vi,1) = MIN( T0, SUM(climate%applied%T2m(vi,:)) / 12._dp)
    END DO
    CALL sync
    
    ! Solve the heat equation for all vertices
    DO vi = mesh%v1, mesh%v2
      
      ! Skip the domain boundary
      IF (mesh%edge_index(vi) > 0) CYCLE
      
      ! Skip ice-less elements
      IF (ice%mask_ice(vi) == 0) THEN
         Ti_new(vi,:) = ice%Ti(vi,1)
         CYCLE
      END IF
      
      ! Ice surface boundary condition
      beta(1)  = 1._dp
      gamma(1) = 0._dp
      delta(1) = ice%Ti(vi,1)
  
      ! Loop over the whole vertical domain but not the surface (k=1) and the bottom (k=NZ):
      DO k = 2, C%NZ-1
      
        CALL get_upwind_derivative_vertex_3D(mesh, ice%U_3D, ice%V_3D, ice%Ti, vi, k, dTi_dx, dTi_dy)
        
        IF (ice%mask_sheet(vi) == 1) THEN
          internal_heating = ((- grav * C%zeta(k)) / ice%Cpi(vi,k)) * ( &
               (p_zeta%a_zeta(k) * ice%U_3D(vi,k-1) + p_zeta%b_zeta(k) * ice%U_3D(vi,k) + p_zeta%c_zeta(k) * ice%U_3D(vi,k+1)) * ice%dHs_dx(vi) + &
               (p_zeta%a_zeta(k) * ice%V_3D(vi,k-1) + p_zeta%b_zeta(k) * ice%V_3D(vi,k) + p_zeta%c_zeta(k) * ice%V_3D(vi,k+1)) * ice%dHs_dy(vi) )
        ELSE
          internal_heating = 0._dp
        END IF

        f1 = (ice%Ki(vi,k) * ice%dzeta_dz(vi)**2) / (ice_density * ice%Cpi(vi,k))

        f2 = ice%dzeta_dt(vi,k) + ice%dzeta_dx(vi,k) * ice%U_3D(vi,k) + ice%dzeta_dy(vi,k) * ice%V_3D(vi,k) + ice%dzeta_dz(vi) * ice%W_3D(vi,k)

        f3 = internal_heating + (ice%U_3D(vi,k) * dTi_dx + ice%V_3D(vi,k) * dTi_dy) - ice%Ti(vi,k) / C%dt_thermo

        alpha(k) = f1 * p_zeta%a_zetazeta(k) - f2 * p_zeta%a_zeta(k)
        beta (k) = f1 * p_zeta%b_zetazeta(k) - f2 * p_zeta%b_zeta(k) - 1._dp / C%dt_thermo
        gamma(k) = f1 * p_zeta%c_zetazeta(k) - f2 * p_zeta%c_zeta(k)
        delta(k) = f3
        
!        ! DENK DROM - only vertical diffusion
!        alpha(k) = (p_zeta%a_zeta(k) * ice%dzeta_dt(vi,k)) - ((p_zeta%a_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2))
!        beta( k) = (p_zeta%b_zeta(k) * ice%dzeta_dt(vi,k)) - ((p_zeta%b_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2)) + (1._dp / C%dt_thermo)
!        gamma(k) = (p_zeta%c_zeta(k) * ice%dzeta_dt(vi,k)) - ((p_zeta%c_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2))
!        delta(k) = ice%Ti(vi,k) / C%dt_thermo
        
!        ! DENK DROM - only vertical diffusion + vertical advection
!        alpha(k) = (p_zeta%a_zeta(k) * (ice%dzeta_dt(vi,k) - ice%W_3D(vi,k) / ice%Hi(vi))) - ((p_zeta%a_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2))
!        beta( k) = (p_zeta%b_zeta(k) * (ice%dzeta_dt(vi,k) - ice%W_3D(vi,k) / ice%Hi(vi))) - ((p_zeta%b_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2)) + (1._dp / C%dt_thermo)
!        gamma(k) = (p_zeta%c_zeta(k) * (ice%dzeta_dt(vi,k) - ice%W_3D(vi,k) / ice%Hi(vi))) - ((p_zeta%c_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2))
!        delta(k) = ice%Ti(vi,k) / C%dt_thermo

      END DO ! DO k = 2, C%NZ-1
 
      IF (ice%mask_shelf(vi) == 1 .OR. ice%mask_gl(vi) == 1) THEN
        ! Set ice bottom temperature to seawater temperature
        alpha(C%NZ) = 0._dp
        beta (C%NZ) = 1._dp
        delta(C%NZ) = SMT
      ELSE
        ! Neumann accoring to GHF
        alpha(C%NZ) = 1._dp
        beta (C%NZ) = -1._dp
        delta(C%NZ) = (C%zeta(C%NZ) - C%zeta(C%NZ-1)) * (ice%GHF( vi) + ice%frictional_heating(vi)) / (ice%dzeta_dz(vi) * ice%Ki(vi,C%NZ)) 
        ! Mixed boundary condition depending on PMP limit
        IF (ice%Ti(vi,C%NZ) >= ice%Ti_pmp(vi,C%NZ)) THEN
          ! Dirichlet at PMP
          alpha(C%NZ) = 0._dp
          beta (C%NZ) = 1._dp
          delta(C%NZ) = ice%Ti_pmp(vi,C%NZ)
        END IF
      END IF ! IF (ice%mask_shelf(vi) == 1 .OR. ice%mask_groundingline(vi) == 1) THEN

      Ti_new(vi,:) = tridiagonal_solve(alpha, beta, gamma, delta, 'thermodynamics_module [temperature]')
      
      ! Make sure ice temperature doesn't exceed pressure melting point
      DO k = 1, C%nZ-1
        Ti_new(vi,k) = MIN(Ti_new(vi,k), ice%Ti_pmp(vi,k))
      END DO
      
      IF (Ti_new(vi,C%NZ) >= ice%Ti_pmp(vi,C%NZ)) THEN
        Ti_new(vi,C%NZ) = MIN( ice%Ti_pmp(vi,C%NZ), ice%Ti(vi,C%NZ-1) - (C%zeta(C%NZ) - C%zeta(C%NZ-1)) * (ice%GHF( vi) + ice%frictional_heating(vi)) / (ice%dzeta_dz(vi) * ice%Ki(vi,C%NZ)))
      END IF
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
        
    CALL apply_Neumann_boundary_3D( mesh, Ti_new, C%nZ)
    CALL sync
    
    ice%Ti(mesh%v1:mesh%v2,:) = Ti_new(mesh%v1:mesh%v2,:)
    CALL sync
    
    CALL deallocate_shared( wTi_new)
    NULLIFY( Ti_new)
    
    ! Safety - to prevent the rare instabilities in the heat equation solver from stopping the entire simulation,
    ! find vertices where instability develops (indicated by an ice temperature below 150K) and replace their
    ! temperature profile with the Robin solution. If too many (>1% of nV) vertices become unstable, throw an error.
            
    n_unstable = 0    
    DO vi = mesh%v1, mesh%v2    
      IF (MINVAL(ice%Ti(vi,:)) < 150._dp) THEN  
            
        CALL replace_Ti_with_robin_solution( ice, climate, SMB, vi)
        
        n_unstable = n_unstable + 1
        
      END IF    
    END DO
    CALL sync
    
    ! Check how many unstable vertices were detected.
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_unstable, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    IF (n_unstable > CEILING(REAL(mesh%nV) / 100._dp)) THEN
      IF (par%master) WRITE(0,*) '   ERROR - thermodynamics:  heat equation solver unstable for more than 1% of vertices!'
      STOP
    END IF

  END SUBROUTINE update_ice_temperature
  
  SUBROUTINE replace_Ti_with_robin_solution( ice, climate, SMB, vi)
    ! This function calculates for one horizontal grid point the temperature profiles
    ! using the surface temperature and the geothermal heat flux as boundary conditions.
    ! See Robin solution in: Cuffey & Paterson 2010, 4th ed, chapter 9, eq. (9.13) - (9.22).
    
    USE parameters_module, ONLY: pi, sec_per_year
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_climate_model),            INTENT(IN)    :: climate
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    INTEGER,                             INTENT(IN)    :: vi

    ! Local variables:
    INTEGER                                            :: k
    REAL(dp)                                           :: Ts
    REAL(dp)                                           :: thermal_length_scale
    REAL(dp)                                           :: distance_above_bed
    REAL(dp)                                           :: erf1
    REAL(dp)                                           :: erf2
    
    REAL(dp)                                           :: thermal_conductivity_robin
    REAL(dp)                                           :: thermal_diffusivity_robin
    REAL(dp)                                           :: bottom_temperature_gradient_robin
    
    REAL(dp), PARAMETER                                :: kappa_0_ice_conductivity     = 9.828_dp                   ! The linear constant in the thermal conductivity of ice [J m^-1 K^-1 s^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    REAL(dp), PARAMETER                                :: kappa_e_ice_conductivity     = 0.0057_dp                  ! The exponent constant in the thermal conductivity of ice [K^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    REAL(dp), PARAMETER                                :: c_0_specific_heat            = 2127.5_dp                  ! The constant in the specific heat capacity of ice [J kg^-1 K^-1], see equation (12.5), Zwinger (2007), Cuffey & Paterson (2010, p. 400)
    REAL(dp), PARAMETER                                :: Claus_Clap_gradient          = 8.7E-04_dp                 ! config variable   ! Clausius Clapeyron gradient [K m^-1]
    
    thermal_conductivity_robin        = kappa_0_ice_conductivity * sec_per_year * EXP(-kappa_e_ice_conductivity * T0)           ! Thermal conductivity            [J m^-1 K^-1 y^-1]
    thermal_diffusivity_robin         = thermal_conductivity_robin / (ice_density * c_0_specific_heat)         ! Thermal diffusivity             [m^2 y^-1]
    bottom_temperature_gradient_robin = - ice%GHF( vi) / thermal_conductivity_robin                    ! Temperature gradient at bedrock
    
    Ts = MIN( T0, SUM(climate%applied%T2m(vi,:)) / 12._dp)
    
    IF (ice%mask_sheet(vi) == 1 ) THEN
    
      IF (SMB%SMB_year(vi) > 0._dp) THEN    
        ! The Robin solution can be used to estimate the subsurface temperature profile in an accumulation area
        
        thermal_length_scale = SQRT(2._dp * thermal_diffusivity_robin * ice%Hi(vi) / SMB%SMB_year(vi))
        DO k = 1, C%nZ
          distance_above_bed = (1._dp - C%zeta(k)) * ice%Hi(vi)
          erf1 = erf( distance_above_bed / thermal_length_scale)
          erf2 = erf( ice%Hi(vi) / thermal_length_scale)
          ice%Ti(vi,k) = Ts + SQRT(pi) / 2._dp * thermal_length_scale * bottom_temperature_gradient_robin * (erf1 - erf2)
        END DO
      
      ELSE
    
        ! Ablation area: use linear temperature profile from Ts to (offset below) T_pmp
        ice%Ti(vi,:) = Ts + ((T0 - Claus_Clap_gradient * ice%Hi(vi)) - Ts) * C%zeta(:)
      
      END IF
      
    ELSEIF( ice%mask_shelf(vi) == 1) THEN
    
      ! Use a linear profile between Ts and seawater temperature:
      ice%Ti(vi,:) = Ts + C%zeta(:) * (SMT - Ts)
      
    ELSE
    
      ! No ice present: use Ts everywhere
      ice%Ti(vi,:) = Ts
      
    END IF

    ! Correct all temperatures above T_pmp:
    DO k = 1, C%NZ
      ice%Ti(vi,k) = MIN( ice%Ti(vi,k), T0 - Claus_Clap_gradient * ice%Hi(vi) * C%zeta(k))
    END DO

  END SUBROUTINE replace_Ti_with_robin_solution  
  
  SUBROUTINE bottom_frictional_heating( mesh, ice)
    ! Calculation of the frictional heating at the bottom due to sliding at the sheet/Gl - bedrock interface.
    
    USE parameters_module,           ONLY: ice_density, grav
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables
    INTEGER                                            :: vi
    REAL(dp)                                           :: beta_base
    
    REAL(dp), PARAMETER                                :: delta_v              = 1E-3_dp       ! Normalisation parameter to prevent errors when velocity is zero
    REAL(dp), PARAMETER                                :: q_plastic            = 0.30_dp       ! Parameter used for basal stress (inverse of m_flow)
    REAL(dp), PARAMETER                                :: u_threshold          = 100._dp       ! scaling of tau_yield to get the correct unit (function of q_plastic)

    ice%frictional_heating( mesh%v1:mesh%v2) = 0._dp
    CALL sync
    
    DO vi = mesh%v1, mesh%v2
      IF (ice%mask_sheet(vi) == 1) THEN      
        beta_base = ice%tau_c_AaAc( vi) * ( (delta_v**2 + ice%U_SSA( vi)**2 + ice%V_SSA( vi)**2)**(0.5_dp * (q_plastic-1._dp)) ) / (u_threshold**q_plastic)
        ice%frictional_heating( vi) = beta_base * (ice%U_SSA( vi)**2 + ice%V_SSA( vi)**2)          
      END IF      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE bottom_frictional_heating
  
  FUNCTION tridiagonal_solve( ldiag, diag, udiag, rhs, string_error_message) RESULT(x)
    ! Lapack tridiagnal solver (in double precision):
    ! Matrix system solver for tridiagonal matrices. 
    ! Used e.g. in solving the ADI scheme. 
    ! ldiag = lower diagonal elements (j,j-1) of the matrix
    ! diag  = diagonal elements (j,j) of the matrix
    ! udiag = upper diagonal elements (j,j+1) of the matrix
    ! rhs   = right hand side of the matrix equation in the ADI scheme
    USE configuration_module, ONLY: dp
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:),            INTENT(IN) :: diag
    REAL(dp), DIMENSION(SIZE(diag)-1), INTENT(IN) :: udiag, ldiag
    REAL(dp), DIMENSION(SIZE(diag)),   INTENT(IN) :: rhs
    CHARACTER(LEN=*),                  INTENT(IN) :: string_error_message

    ! Result variable:
    REAL(dp), DIMENSION(SIZE(diag))               :: x
    
    ! Local variables:     
    INTEGER                                       :: info
    REAL(dp), DIMENSION(SIZE(diag))               :: diag_copy
    REAL(dp), DIMENSION(SIZE(udiag))              :: udiag_copy, ldiag_copy

    ! External subroutines:      
    EXTERNAL DGTSV ! Lapack routine that solves tridiagonal systems (in double precision).

    ! The LAPACK solver will overwrite the rhs with the solution x. Therefore we 
    ! first copy the rhs in the solution vector x:
    x = rhs

    ! The LAPACK solver will change the elements in the matrix, therefore we copy them:
    diag_copy  =  diag
    udiag_copy = udiag
    ldiag_copy = ldiag

    CALL DGTSV(SIZE(diag), 1, ldiag_copy, diag_copy, udiag_copy, x, SIZE(diag), info)
    ! Check if solver was successfull:
    IF(info /= 0) THEN
     WRITE(*, FMT='(3a, i5)') ' In the module ', string_error_message, ' in the function tridiagonal_solve_ant: info=', info
     !WRITE(UNIT=*, FMT='(4E30.12)') ((ldiag_copy(i), diag_copy(i), udiag_copy(i), x(i)), i=1, NZ) ! test print
     STOP ' DGTSV problem with tridiagonal system, --STOPPED'
    END IF
    
  END FUNCTION tridiagonal_solve
  
  SUBROUTINE initialise_ice_temperature( mesh, ice, init, climate)
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
    TYPE(type_climate_model),            INTENT(IN)    :: climate
    
    ! Local variables
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: vi, k
    REAL(dp)                                           :: T_surf_annual
    
    ! If we're doing a restart, initialise with that
    IF (C%is_restart) THEN
      ice%Ti( mesh%v1:mesh%v2,:) = init%Ti( mesh%v1:mesh%v2,:)
      CALL sync
      RETURN
    END IF
    
    ! First set all temperatures to -10C so thermal properties can be determined
    ice%Ti( mesh%v1:mesh%v2,:) = 260._dp
    CALL sync
   
    ! Calculate Ti_pmp
    CALL ice_physical_properties( mesh, ice, C%start_time_of_run)
    
    ! Initialise with a linear profile
    DO vi = mesh%v1, mesh%v2
      T_surf_annual = MIN( SUM(climate%applied%T2m(vi,:))/12._dp, T0)
      IF (ice%Hi(vi) > 0._dp) THEN
        DO k = 1, C%nZ
          ice%Ti(vi,k) = T_surf_annual + C%zeta(k) * (ice%Ti_pmp(vi,C%nZ) - T_surf_annual)
        END DO
      ELSE
        ice%Ti(vi,:) = T_surf_annual
      END IF        
    END DO
    CALL sync   
    
    ! Special cases for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3') THEN
              
        ice%Ti( mesh%v1:mesh%v2,:) = 270._dp
        CALL sync
        
      ELSEIF (C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6') THEN
              
        DO vi = mesh%v1, mesh%v2
          ice%Ti( vi,:) = climate%applied%T2m( vi,1)
        END DO
        CALL sync
              
      ELSEIF (C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler' .OR. &
              C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
              C%choice_benchmark_experiment == 'mesh_generation_test') THEN
              
        ice%Ti( mesh%v1:mesh%v2,:) = 270._dp
        CALL sync
        
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_ice_temperature!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  END SUBROUTINE initialise_ice_temperature  
  
END MODULE thermodynamics_module
