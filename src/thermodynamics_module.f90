MODULE thermodynamics_module

  ! All the routines for calculating the englacial temperature profile

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
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_subclimate_region, type_SMB_model
  USE zeta_module,                     ONLY: calculate_zeta_derivatives, p_zeta
  USE utilities_module,                ONLY: tridiagonal_solve, vertical_average
  USE mesh_operators_module,           ONLY: apply_Neumann_BC_direct_3D, ddx_a_to_a_2D, ddy_a_to_a_2D, move_aca_to_a_2D
  USE mesh_help_functions_module,      ONLY: CROSS2
  
  IMPLICIT NONE
  
CONTAINS
   
! == Run the chosen thermodynamics model
  SUBROUTINE run_thermo_model( mesh, ice, climate, SMB, time, do_solve_heat_equation)
    ! Run the thermodynamics model. If so specified, solve the heat equation;
    ! if not, only prescribe a vertically uniform temperature to newly ice-covered grid cells.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_subclimate_region),         INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    REAL(dp),                             INTENT(IN)    :: time
    LOGICAL,                              INTENT(IN)    :: do_solve_heat_equation

    ! Local variables:
    INTEGER                                             :: vi
    REAL(dp)                                            :: T_surf_annual
    
    IF     (C%choice_thermo_model == 'none') THEN
      ! No need to do anything
      ! NOTE: choice_ice_rheology_model should be set to "uniform"!
    ELSEIF (C%choice_thermo_model == '3D_heat_equation') THEN
      ! Solve the 3-D heat equation
      
      ! NOTE: solved asynchronously from the ice dynamical equations.
      !       Since newly ice-covered pixels won't have a temperature assigned
      !       until the heat equation is solved again, treat these separately every time step.
    
      ! Prescribe a simple temperature profile to newly ice-covered grid cells.
      DO vi = mesh%vi1, mesh%vi2
        
        IF (ice%mask_ice_a( vi) == 1 .AND. ice%mask_ice_a_prev( vi) == 0) THEN
          ! This grid cell is newly ice-covered
          ! If one of its neighbours was already ice-covered, assume the temperature
          ! profile here is equal to the profile from the upstream neighbour (due to advection).
          ! If no neighbours were ice-covered, the new ice must come from accumulation;
          ! just set a simple linear profile instead.
          
!          IF     (ice%u_vav_cx( j  ,i-1) > 0._dp .AND. ice%mask_ice_a_prev( j  ,i-1) == 1) THEN
!            ! Ice probably came from the west
!            ice%Ti_a( :,j,i) = ice%Ti_a( :,j  ,i-1)
!          ELSEIF (ice%u_vav_cx( j  ,i  ) < 0._dp .AND. ice%mask_ice_a_prev( j  ,i+1) == 1) THEN
!            ! Ice probably came from the east
!            ice%Ti_a( :,j,i) = ice%Ti_a( :,j  ,i+1)
!          ELSEIF (ice%v_vav_cy( j-1,i  ) > 0._dp .AND. ice%mask_ice_a_prev( j-1,i  ) == 1) THEN
!            ! Ice probably came from the south
!            ice%Ti_a( :,j,i) = ice%Ti_a( :,j-1,i  )
!          ELSEIF (ice%v_vav_cy( j  ,i  ) < 0._dp .AND. ice%mask_ice_a_prev( j+1,i  ) == 1) THEN
!            ! Ice probably came from the north
!            ice%Ti_a( :,j,i) = ice%Ti_a( :,j+1,i  )
!          ELSE
            ! Ice probably came from surface accumulation; initialise with a vertically uniform temperature.
            
            T_surf_annual = MIN( SUM( climate%T2m( vi,:)) / 12._dp, T0)
            ice%Ti_a( vi,:) = T_surf_annual
            
!          END IF
          
        END IF ! IF (ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a_prev( j,i) == 0) THEN
        
      END DO
      CALL sync
    
      ! Calculate various physical terms
      CALL calc_heat_capacity(          mesh, ice)
      CALL calc_thermal_conductivity(   mesh, ice)
      CALL calc_pressure_melting_point( mesh, ice)
      
      ! If so specified, solve the heat equation
      IF (do_solve_heat_equation) CALL solve_3D_heat_equation( mesh, ice, climate, SMB)
      
      ! Safety
      CALL check_for_NaN_dp_2D( ice%Ti_a, 'ice%Ti_a', 'run_thermo_model')
    
    ELSE
      IF (par%master) WRITE(0,*) 'run_thermo_model - ERROR: unknown choice_thermo_model "', TRIM(C%choice_thermo_model), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Calculate the ice flow factor for the new temperature solution
    CALL calc_ice_rheology( mesh, ice, time)

  END SUBROUTINE run_thermo_model
  
! == Solve the 3-D heat equation
  SUBROUTINE solve_3D_heat_equation( mesh, ice, climate, SMB)
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_subclimate_region),        INTENT(IN)    :: climate
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB

    ! Local variables:
    INTEGER                                            :: vi, k
    REAL(dp)                                           :: u_times_dT_dx_upwind, v_times_dT_dy_upwind, f1, f2, f3
    REAL(dp), DIMENSION(2:C%nz)                        :: alpha
    REAL(dp), DIMENSION(C%nz)                          :: beta
    REAL(dp), DIMENSION(C%nz-1)                        :: gamma
    REAL(dp), DIMENSION(C%nz)                          :: delta
    REAL(dp), DIMENSION(:,:  ), POINTER                :: Ti_new
    INTEGER                                            :: wTi_new
    REAL(dp), DIMENSION(:    ), POINTER                ::  T_ocean_at_shelf_base
    INTEGER                                            :: wT_ocean_at_shelf_base
    INTEGER,  DIMENSION(:    ), POINTER                ::  is_unstable
    INTEGER                                            :: wis_unstable
    INTEGER                                            :: n_unstable
    LOGICAL                                            :: hasnan
    
    ! Allocate shared memory
    CALL allocate_shared_int_1D( mesh%nV,       is_unstable,           wis_unstable          )
    CALL allocate_shared_dp_2D(  mesh%nV, C%nz, Ti_new,                wTi_new               )
    CALL allocate_shared_dp_1D(  mesh%nV,       T_ocean_at_shelf_base, wT_ocean_at_shelf_base)
    CALL sync
    
    ! Calculate zeta derivatives required for solving the heat equation
    CALL calculate_zeta_derivatives( mesh, ice)
    
    ! Calculate heating terms
    CALL calc_internal_heating(   mesh, ice)
    CALL calc_frictional_heating( mesh, ice)
    
    ! Set ice surface temperature equal to annual mean 2m air temperature
    DO vi = mesh%vi1, mesh%vi2
      ice%Ti_a( vi,1) = MIN( T0, SUM( climate%T2m( vi,:)) / 12._dp)
    END DO
    CALL sync
    
    ! Find ocean temperature at the shelf base
    DO vi = mesh%vi1, mesh%vi2
    
      T_ocean_at_shelf_base( vi) = SMT
    
!      IF (ice%mask_shelf_a( j,i) == 1) THEN
!        depth = MAX( 0.1_dp, ice%Hi_a( j,i) - ice%Hs_a( j,i))   ! Depth is positive when below the sea surface!
!        CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( :,j,i), depth, T_ocean_at_shelf_base( j,i))
!      ELSE
!        T_ocean_at_shelf_base( j,i) = 0._dp
!      END IF
!      
!      ! NOTE: ocean data gives temperature in Celsius, thermodynamics wants Kelvin!
!      T_ocean_at_shelf_base( j,i) = T_ocean_at_shelf_base( j,i) + T0
      
    END DO
    CALL sync
    
    ! Solve the heat equation for all vertices
    is_unstable( mesh%vi1:mesh%vi2) = 0
    n_unstable                    = 0
    DO vi = mesh%vi1, mesh%vi2
      
      ! Skip the domain boundary
      IF (mesh%edge_index( vi) > 0) CYCLE
      
      ! Skip ice-less elements
      IF (ice%mask_ice_a( vi) == 0) THEN
         Ti_new( vi,:) = ice%Ti_a( vi,1)
         CYCLE
      END IF
      
      ! Ice surface boundary condition
      beta(  1) = 1._dp
      gamma( 1) = 0._dp
      delta( 1) = ice%Ti_a( vi,1)
  
      ! Loop over the whole vertical domain but not the surface (k=1) and the bottom (k=NZ):
      DO k = 2, C%nz-1
      
        CALL calc_upwind_heat_flux_derivatives_vertex( mesh, ice%u_3D_c, ice%v_3D_c, ice%Ti_a, vi, k, u_times_dT_dx_upwind, v_times_dT_dy_upwind)

        f1 = (ice%Ki_a( vi,k) * ice%dzeta_dz_a( vi)**2) / (ice_density * ice%Cpi_a( vi,k))

        f2 = ice%dzeta_dt_a( vi,k) + ice%dzeta_dx_a( vi,k) * ice%u_3D_a(vi,k) + ice%dzeta_dy_a( vi,k) * ice%v_3D_a( vi,k) + ice%dzeta_dz_a( vi) * ice%w_3D_a( vi,k)
 
        f3 = ice%internal_heating_a( vi,k) + (u_times_dT_dx_upwind + v_times_dT_dy_upwind) - ice%Ti_a( vi,k) / C%dt_thermo

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

      END DO ! DO k = 2, C%nz-1
      
      ! Boundary conditions at the surface: set ice temperature equal to annual mean surface temperature
      beta(  1) = 1._dp
      gamma( 1) = 0._dp
      delta( 1) =  MIN( T0, SUM( climate%T2m( vi,:)) / 12._dp)
      
      ! Boundary conditions at the base
      IF (ice%mask_shelf_a( vi) == 1) THEN
        ! Set ice bottom temperature equal to seawater temperature (limited to the PMP)
        alpha( C%nz) = 0._dp
        beta ( C%nz) = 1._dp
        delta( C%nz) = MIN( T0, MIN( ice%Ti_pmp_a( vi,C%nz), T_ocean_at_shelf_base( vi) ))
      ELSE
        IF (ice%Ti_a( vi,C%nz) >= ice%Ti_pmp_a( vi,C%nz)) THEN
          ! Ice is already at/above pressure melting point; set temperature equal to PMP
          alpha( C%nz) = 0._dp
          beta ( C%nz) = 1._dp
          delta( C%nz) = ice%Ti_pmp_a( vi,C%nz)
        ELSE
          ! Set a Neumann BC so the temperature gradient at the base is equal to basal heating rate (= geothermal + friction)
          alpha( C%nz) = 1._dp
          beta ( C%nz) = -1._dp
          delta( C%nz) = (C%zeta(C%nz) - C%zeta(C%nz-1)) * (ice%GHF_a( vi) + ice%frictional_heating_a( vi)) / (ice%dzeta_dz_a( vi) * ice%Ki_a( vi,C%nz)) 
        END IF
      END IF ! IF (ice%mask_shelf_a( vi) == 1) THEN

      ! Solve the tridiagonal matrix equation representing the heat equation for this grid cell
      Ti_new( vi,:) = tridiagonal_solve( alpha, beta, gamma, delta)
      
      ! Make sure ice temperature doesn't exceed pressure melting point
      DO k = 1, C%nz-1
        Ti_new( vi,k) = MIN( Ti_new( vi,k), ice%Ti_pmp_a( vi,k))
      END DO
      
      IF (Ti_new( vi,C%nz) >= ice%Ti_pmp_a( vi,C%nz)) THEN
        Ti_new( vi,C%nz) = MIN( ice%Ti_pmp_a( vi,C%nz), ice%Ti_a( vi,C%nz-1) - (C%zeta(C%nz) - C%zeta(C%nz-1)) * &
          (ice%GHF_a( vi) + ice%frictional_heating_a( vi)) / (ice%dzeta_dz_a( vi) * ice%Ki_a( vi,C%nz)))
      END IF
    
      ! Mark temperatures below 150 K or NaN as unstable, to be replaced with the Robin solution.
      hasnan = .FALSE.
      DO k = 1, C%nz
        IF (Ti_new( vi,k) /= Ti_new( vi,k)) THEN
          hasnan = .TRUE.
        END IF
      END DO
      IF (MINVAL(Ti_new( vi,:)) < 150._dp .OR. hasnan) THEN
        is_unstable( vi) = 1
        n_unstable  = n_unstable + 1
        !WRITE(0,*) 'instability detected; Hi = ', ice%Hi_a( j,i), ', dHi_dt = ', ice%dHi_dt_a( j,i)
      END IF
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
        
    CALL apply_Neumann_BC_direct_3D( mesh, Ti_new)
    
    ! Cope with instability
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_unstable, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)    
    IF (n_unstable < CEILING( REAL( mesh%nV) / 100._dp)) THEN
      ! Instability is limited to an acceptably small number (< 1%) of grid cells;
      ! replace the temperature profile in those cells with the Robin solution
      
      DO vi = mesh%vi1, mesh%vi2
        IF (is_unstable( vi) == 1) CALL replace_Ti_with_robin_solution( ice, climate, SMB, Ti_new, vi)
      END DO
      CALL sync
      
    ELSE
      ! An unacceptably large number of grid cells was unstable; throw an error.
      
      IF (par%master) WRITE(0,*) '   solve_heat_equation - ERROR:  heat equation solver unstable for more than 1% of grid cells!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      
    END IF
    
    ! Move the new temperature field to the ice data structure
    ice%Ti_a( mesh%vi1:mesh%vi2,:) = Ti_new( mesh%vi1:mesh%vi2,:)
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wTi_new)
    CALL deallocate_shared( wis_unstable)
    CALL deallocate_shared( wT_ocean_at_shelf_base)
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%Ti_a, 'ice%Ti_a', 'solve_heat_equation')

  END SUBROUTINE solve_3D_heat_equation
  
! == Calculate upwind heat flux derivatives
  SUBROUTINE calc_upwind_heat_flux_derivatives_vertex( mesh, u_3D_c, v_3D_c, Ti_a, vi, k, u_times_dT_dx_upwind, v_times_dT_dy_upwind)
    ! Calculate upwind heat flux derivatives at vertex vi, vertical layer k
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_3D_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_3D_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Ti_a
    INTEGER,                             INTENT(IN)    :: vi
    INTEGER,                             INTENT(IN)    :: k
    REAL(dp),                            INTENT(OUT)   :: u_times_dT_dx_upwind
    REAL(dp),                            INTENT(OUT)   :: v_times_dT_dy_upwind
    
    ! Local variables:
    INTEGER                                            :: iti, ti, n1, n2, n3, vib, vic, ti_upwind
    REAL(dp), DIMENSION(2)                             :: u_upwind, ab, ac
    INTEGER                                            :: kk, vj
    
    ! Upwind velocity vector
    u_upwind = [u_3D_c( vi,k), v_3D_c( vi,k)]
    
    ! Find the upwind triangle
    ti_upwind = 0
    DO iti = 1, mesh%niTri( vi)
      
      ! Triangle ti is spanned counter-clockwise by vertices [vi,vib,vic]
      ti  = mesh%iTri( vi,iti)
      vib = 0
      vic = 0
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        n3 = n2 + 1
        IF (n3 == 4) n3 = 1
        
        IF (mesh%Tri( ti,n1) == vi) THEN
          vib = mesh%Tri( ti,n2)
          vic = mesh%Tri( ti,n3)
          EXIT
        END IF
      END DO
      
      ! Check if the upwind velocity vector points into this triangle
      ab = mesh%V( vib,:) - mesh%V( vi,:)
      ac = mesh%V( vic,:) - mesh%V( vi,:)
      
      IF (CROSS2( ab, u_upwind) >= 0._dp .AND. CROSS2( u_upwind, ac) >= 0._dp) THEN
        ti_upwind = ti
        EXIT
      END IF
      
    END DO ! DO iti = 1, mesh%niTri( vi)
    
    ! Safety
    IF (ti_upwind == 0) THEN
      WRITE(0,*) 'calc_upwind_heat_flux_derivatives_vertex - ERROR: couldnt find upwind triangle!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Calculate derivatives on this triangle
    
    u_times_dT_dx_upwind = 0._dp
    DO kk = mesh%M_ddx_a_b%ptr( ti_upwind), mesh%M_ddx_a_b%ptr( ti_upwind+1)-1
      vj = mesh%M_ddx_a_b%index( kk)
      u_times_dT_dx_upwind = u_times_dT_dx_upwind + Ti_a( vj,k) * mesh%M_ddx_a_b%val( kk)
    END DO
    
    v_times_dT_dy_upwind = 0._dp
    DO kk = mesh%M_ddy_a_b%ptr( ti_upwind), mesh%M_ddy_a_b%ptr( ti_upwind+1)-1
      vj = mesh%M_ddy_a_b%index( kk)
      v_times_dT_dy_upwind = v_times_dT_dy_upwind + Ti_a( vj,k) * mesh%M_ddy_a_b%val( kk)
    END DO
    
  END SUBROUTINE calc_upwind_heat_flux_derivatives_vertex
  
! == The Robin temperature solution
  SUBROUTINE replace_Ti_with_robin_solution( ice, climate, SMB, Ti, vi)
    ! This function calculates for one horizontal grid point the temperature profiles
    ! using the surface temperature and the geothermal heat flux as boundary conditions.
    ! See Robin solution in: Cuffey & Paterson 2010, 4th ed, chapter 9, eq. (9.13) - (9.22).
    
    USE parameters_module, ONLY: pi, sec_per_year
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_subclimate_region),        INTENT(IN)    :: climate
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: Ti
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
    bottom_temperature_gradient_robin = - ice%GHF_a( vi) / thermal_conductivity_robin                    ! Temperature gradient at bedrock
    
    Ts = MIN( T0, SUM(climate%T2m( vi,:)) / 12._dp)
    
    IF (ice%mask_sheet_a( vi) == 1 ) THEN
    
      IF (SMB%SMB_year( vi) > 0._dp) THEN    
        ! The Robin solution can be used to estimate the subsurface temperature profile in an accumulation area
        
        thermal_length_scale = SQRT(2._dp * thermal_diffusivity_robin * ice%Hi_a( vi) / SMB%SMB_year( vi))
        DO k = 1, C%nz
          distance_above_bed = (1._dp - C%zeta(k)) * ice%Hi_a( vi)
          erf1 = erf( distance_above_bed / thermal_length_scale)
          erf2 = erf( ice%Hi_a( vi) / thermal_length_scale)
          Ti( vi,k) = Ts + SQRT(pi) / 2._dp * thermal_length_scale * bottom_temperature_gradient_robin * (erf1 - erf2)
        END DO
      
      ELSE
    
        ! Ablation area: use linear temperature profile from Ts to (offset below) T_pmp
        Ti( vi,:) = Ts + ((T0 - Claus_Clap_gradient * ice%Hi_a( vi)) - Ts) * C%zeta(:)
      
      END IF
      
    ELSEIF( ice%mask_shelf_a(vi) == 1) THEN
    
      ! Use a linear profile between Ts and seawater temperature:
      Ti( vi,:) = Ts + C%zeta(:) * (SMT - Ts)
      
    ELSE
    
      ! No ice present: use Ts everywhere
      Ti( vi,:) = Ts
      
    END IF

    ! Correct all temperatures above T_pmp:
    DO k = 1, C%nz
      Ti( vi,k) = MIN( Ti( vi,k), T0 - Claus_Clap_gradient * ice%Hi_a( vi) * C%zeta(k))
    END DO

  END SUBROUTINE replace_Ti_with_robin_solution  
  
! == Calculate various physical terms
  SUBROUTINE calc_internal_heating( mesh, ice)
    ! Calculate internal heating due to deformation
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, k
    REAL(dp), DIMENSION(:    ), POINTER                ::  dHs_dx_a,  dHs_dy_a
    INTEGER                                            :: wdHs_dx_a, wdHs_dy_a
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, dHs_dx_a, wdHs_dx_a)
    CALL allocate_shared_dp_1D( mesh%nV, dHs_dy_a, wdHs_dy_a)
    
    ! Calculate surface slopes
    CALL ddx_a_to_a_2D( mesh, ice%Hs_a, dHs_dx_a)
    CALL ddy_a_to_a_2D( mesh, ice%Hs_a, dHs_dy_a)
    
    ! Calculate internal heating
    DO vi = mesh%vi1, mesh%vi2
      
      ice%internal_heating_a( vi,:) = 0._dp
      
      IF (mesh%edge_index( vi) > 0) CYCLE ! Skip the domain boundary
      IF (ice%mask_ice_a( vi) == 0) CYCLE ! Skip ice-less elements
  
      ! Loop over the whole vertical domain but not the surface (k=1) and the bottom (k=NZ):
      DO k = 2, C%nz-1
        ice%internal_heating_a( vi,k) = ((- grav * C%zeta(k)) / ice%Cpi_a( vi,k)) * ( &
             (p_zeta%a_zeta(k) * ice%u_3D_a( vi,k-1) + p_zeta%b_zeta(k) * ice%u_3D_a( vi,k) + p_zeta%c_zeta(k) * ice%u_3D_a( vi,k+1)) * dHs_dx_a( vi) + &
             (p_zeta%a_zeta(k) * ice%v_3D_a( vi,k-1) + p_zeta%b_zeta(k) * ice%v_3D_a( vi,k) + p_zeta%c_zeta(k) * ice%v_3D_a( vi,k+1)) * dHs_dy_a( vi) )
      END DO
      
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx_a)
    CALL deallocate_shared( wdHs_dy_a)
    
  END SUBROUTINE calc_internal_heating
  SUBROUTINE calc_frictional_heating( mesh, ice)
    ! Calculate frictional heating at the base due to sliding
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(:    ), POINTER                ::  beta_a
    INTEGER                                            :: wbeta_a
    
    ! Exception for when no sliding can occur
    IF (C%choice_ice_dynamics == 'SIA' .OR. C%choice_sliding_law == 'no_sliding') THEN
      ice%frictional_heating_a( mesh%vi1:mesh%vi2) = 0._dp
      CALL sync
      RETURN
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, beta_a, wbeta_a)
    
    ! Map beta from the combi-mesh to the regular mesh
    CALL move_aca_to_a_2D( mesh, ice%beta_ac, beta_a)
    
    ! Calculate frictional heating
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_sheet_a( vi) == 1) THEN
        ice%frictional_heating_a( vi) = beta_a( vi) * (ice%u_base_a( vi)**2 + ice%u_base_a( vi)**2)
      ELSE
        ice%frictional_heating_a( vi) = 0._dp
      END IF 
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wbeta_a)

  END SUBROUTINE calc_frictional_heating
  SUBROUTINE calc_heat_capacity( mesh, ice)
    ! Calculate the heat capacity of the ice
      
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
  
    ! Local variables:
    INTEGER                                            :: vi
    
    IF     (C%choice_ice_heat_capacity == 'uniform') THEN
      ! Apply a uniform value for the heat capacity
      
      ice%Cpi_a( mesh%vi1:mesh%vi2,:) = C%uniform_ice_heat_capacity
      CALL sync
      
    ELSEIF (C%choice_ice_heat_capacity == 'Pounder1965') THEN
      ! Calculate the heat capacity of ice according to Pounder: The Physics of Ice (1965)
      
      DO vi = mesh%vi1, mesh%vi2
        ice%Cpi_a( vi,:) = 2115.3_dp + 7.79293_dp * (ice%Ti_a( vi,:) - T0)
      END DO
      CALL sync
    
    ELSE
      IF (par%master) WRITE(0,*) 'calc_heat_capacity - ERROR: unknown choice_ice_heat_capacity "', TRIM(C%choice_ice_heat_capacity), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%Cpi_a, 'ice%Cpi_a', 'calc_heat_capacity')
    
  END SUBROUTINE calc_heat_capacity
  SUBROUTINE calc_thermal_conductivity( mesh, ice)
    ! Calculate the thermal conductivity of the ice
      
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
  
    ! Local variables:
    INTEGER                                            :: vi
    
    IF     (C%choice_ice_thermal_conductivity == 'uniform') THEN
      ! Apply a uniform value for the thermal conductivity
      
      ice%Ki_a( mesh%vi1:mesh%vi2,:) = C%uniform_ice_thermal_conductivity
      CALL sync
      
    ELSEIF (C%choice_ice_thermal_conductivity == 'Ritz1987') THEN
      ! Calculate the thermal conductivity of ice according to Ritz (1987) 
      
      DO vi = mesh%vi1, mesh%vi2
        ice%Ki_a( vi,:) = 3.101E+08_dp * EXP(-0.0057_dp * ice%Ti_a( vi,:))
      END DO
      CALL sync
    
    ELSE
      IF (par%master) WRITE(0,*) 'calc_thermal_conductivity - ERROR: unknown choice_ice_thermal_conductivity "', TRIM(C%choice_ice_thermal_conductivity), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%Ki_a, 'ice%Ki_a', 'calc_thermal_conductivity')
    
  END SUBROUTINE calc_thermal_conductivity
  SUBROUTINE calc_pressure_melting_point( mesh, ice)
    ! Calculate the pressure melting point of the ice according to Huybrechts (1992)
      
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
  
    ! Local variables:
    INTEGER                                            :: vi
    
    DO vi = mesh%vi1, mesh%vi2
      ice%Ti_pmp_a( vi,:) = T0 - CC * ice%Hi_a( vi) * C%zeta
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%Ti_pmp_a, 'ice%Ti_pmp_a', 'calc_pressure_melting_point')
    
  END SUBROUTINE calc_pressure_melting_point
  
! == Calculate the  flow factor A in Glen's flow law
  SUBROUTINE calc_ice_rheology( mesh, ice, time)
    ! Calculate the flow factor A in Glen's flow law
      
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: time
  
    ! Local variables:
    INTEGER                                            :: vi,k
    REAL(dp), DIMENSION(C%nZ)                          :: prof
    REAL(dp), PARAMETER                                :: A_low_temp  = 1.14E-05_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: A_high_temp = 5.47E+10_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_low_temp  = 6.0E+04_dp    ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_high_temp = 13.9E+04_dp   ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp)                                           :: A_flow_MISMIP
    
    IF     (C%choice_ice_rheology == 'uniform') THEN
      ! Apply a uniform value for the ice flow factor
      
      ice%A_flow_3D_a( mesh%vi1:mesh%vi2,:) = C%uniform_flow_factor
      CALL sync
      
    ELSEIF (C%choice_ice_rheology == 'Huybrechts1992') THEN
    
      ! Calculate the ice flow factor as a function of the ice temperature according to the Arrhenius relationship (Huybrechts, 1992)
      DO vi = mesh%vi1, mesh%vi2
      
        DO k = 1, C%nz
          IF (ice%mask_ice_a( vi) == 1) THEN
            IF (ice%Ti_a( vi,k) < 263.15_dp) THEN
              ice%A_flow_3D_a( vi,k) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti_a( vi,k)))  
            ELSE
              ice%A_flow_3D_a( vi,k) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti_a( vi,k)))  
            END IF
          ELSE
            IF (C%choice_ice_margin == 'BC') THEN
              ice%A_flow_3D_a( vi,k) = 0._dp
            ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
              ! In the "infinite slab" case, calculate effective viscosity everywhere
              ! (even when there's technically no ice present)
              ice%A_flow_3D_a( vi,k) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * 263.15_dp))  
            ELSE
              IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in calc_effective_viscosity!'
              CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
            END IF
          END IF
        END DO ! DO k = 1, C%nz
         
      END DO
      CALL sync
    
    ELSEIF (C%choice_ice_rheology == 'MISMIP_mod') THEN
      ! The time-dependent, step-wise changing uniform flow factor in the MISMIP_mod experiment
      
      A_flow_MISMIP = 1.0E-16_dp
      IF     (time < 25000._dp) THEN
        A_flow_MISMIP = 1.0E-16_dp
      ELSEIF (time < 50000._dp) THEN
        A_flow_MISMIP = 1.0E-17_dp
      ELSEIF (time < 75000._dp) THEN
        A_flow_MISMIP = 1.0E-16_dp
      END IF
        
      ice%A_flow_3D_a(  mesh%vi1:mesh%vi2,:) = A_flow_MISMIP
      CALL sync
        
    ELSE
      IF (par%master) WRITE(0,*) 'calc_ice_rheology - ERROR: unknown choice_ice_rheology "', TRIM(C%choice_ice_rheology), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
        
    ! Apply the flow enhancement factors
    DO vi = mesh%vi1, mesh%vi2
      IF     (ice%mask_sheet_a( vi) == 1) THEN
        ice%A_flow_3D_a( vi,:) = ice%A_flow_3D_a( vi,:) * C%m_enh_sheet
      ELSEIF (ice%mask_shelf_a( vi) == 1) THEN
        ice%A_flow_3D_a( vi,:) = ice%A_flow_3D_a( vi,:) * C%m_enh_shelf
      END IF
    END DO
    CALL sync

    ! Calculate vertical average
    DO vi = mesh%vi1, mesh%vi2
      prof = ice%A_flow_3D_a( vi,:)
      CALL vertical_average( prof, ice%A_flow_vav_a( vi))
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%A_flow_3D_a , 'ice%A_flow_3D_a' , 'calc_ice_rheology')
    CALL check_for_NaN_dp_1D( ice%A_flow_vav_a, 'ice%A_flow_vav_a', 'calc_ice_rheology')
    
  END SUBROUTINE calc_ice_rheology
  
! == Initialise the englacial ice temperature at the start of a simulation
  SUBROUTINE initialise_ice_temperature( mesh, ice, climate, SMB, region_name)
    ! Initialise the englacial ice temperature at the start of a simulation
      
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_subclimate_region),         INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name
    
    IF (par%master) WRITE (0,*) '  Initialising ice temperature profile "', TRIM(C%choice_initial_ice_temperature), '"...'
    
    IF     (C%choice_initial_ice_temperature == 'uniform') THEN
      ! Simple uniform temperature
      CALL initialise_ice_temperature_uniform( mesh, ice)
    ELSEIF (C%choice_initial_ice_temperature == 'linear') THEN
      ! Simple linear temperature profile
      CALL initialise_ice_temperature_linear( mesh, ice, climate)
    ELSEIF (C%choice_initial_ice_temperature == 'Robin') THEN
      ! Initialise with the Robin solution
      CALL initialise_ice_temperature_Robin( mesh, ice, climate, SMB)
    ELSEIF (C%choice_initial_ice_temperature == 'restart') THEN
      ! Initialise with the temperature field from the provided restart file
      CALL initialise_ice_temperature_restart( mesh, ice, region_name)
    ELSE
      IF (par%master) WRITE(0,*) 'initialise_ice_temperature - ERROR: unknown choice_initial_ice_temperature "', TRIM(C%choice_thermo_model), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE initialise_ice_temperature
  SUBROUTINE initialise_ice_temperature_uniform( mesh, ice)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Simple uniform temperature
      
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables
    INTEGER                                            :: vi
      
    DO vi = mesh%vi1, mesh%vi2
    
      IF (ice%Hi_a( vi) > 0._dp) THEN
        ice%Ti_a( vi,:) = C%uniform_ice_temperature
      ELSE
        ice%Ti_a( vi,:) = 0._dp
      END IF
      
    END DO
    CALL sync
    
  END SUBROUTINE initialise_ice_temperature_uniform
  SUBROUTINE initialise_ice_temperature_linear( mesh, ice, climate)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Simple linear temperature profile
      
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_subclimate_region),         INTENT(IN)    :: climate
    
    ! Local variables
    INTEGER                                             :: vi
    REAL(dp)                                            :: T_surf_annual, T_PMP_base
      
    DO vi = mesh%vi1, mesh%vi2
    
      IF (ice%Hi_a( vi) > 0._dp) THEN
        T_surf_annual = MIN( SUM( climate%T2m( vi,:)) / 12._dp, T0)
        T_PMP_base    = T0 - (ice%Hi_a( vi) * 8.7E-04_dp)
        ice%Ti_a( vi,:) = T_surf_annual - C%zeta * (T_surf_annual - T_PMP_base)
      ELSE
        ice%Ti_a( vi,:) = 0._dp
      END IF
      
    END DO
    CALL sync
    
  END SUBROUTINE initialise_ice_temperature_linear
  SUBROUTINE initialise_ice_temperature_Robin( mesh, ice, climate, SMB)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Initialise with the Robin solution
      
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_subclimate_region),         INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    
    ! Local variables
    INTEGER                                             :: vi
 
    ! Calculate Ti_pmp
    CALL calc_pressure_melting_point( mesh, ice)
    
    ! Initialise with the Robin solution
    DO vi = mesh%vi1, mesh%vi2
      CALL replace_Ti_with_robin_solution( ice, climate, SMB, ice%Ti_a, vi)    
    END DO
    CALL sync
    
  END SUBROUTINE initialise_ice_temperature_Robin
  SUBROUTINE initialise_ice_temperature_restart( mesh, ice, region_name)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Initialise with the temperature field from the provided restart file
      
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    IF (par%master) WRITE(0,*) 'initialise_ice_temperature_restart - FIXME!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    
!    ! Local variables
!    CHARACTER(LEN=256)                                 :: filename_restart
!    REAL(dp)                                           :: time_to_restart_from
!    TYPE(type_restart_data)                            :: restart
!    
!    ! Assume that temperature and geometry are read from the same restart file
!    IF     (region_name == 'NAM') THEN
!      filename_restart     = C%filename_refgeo_init_NAM
!      time_to_restart_from = C%time_to_restart_from_NAM
!    ELSEIF (region_name == 'EAS') THEN
!      filename_restart     = C%filename_refgeo_init_EAS
!      time_to_restart_from = C%time_to_restart_from_EAS
!    ELSEIF (region_name == 'GR:') THEN
!      filename_restart     = C%filename_refgeo_init_GRL
!      time_to_restart_from = C%time_to_restart_from_GRL
!    ELSEIF (region_name == 'ANT') THEN
!      filename_restart     = C%filename_refgeo_init_ANT
!      time_to_restart_from = C%time_to_restart_from_ANT
!    END IF
!    
!    ! Inquire if all the required fields are present in the specified NetCDF file,
!    ! and determine the dimensions of the memory to be allocated.
!    CALL allocate_shared_int_0D( restart%nx, restart%wnx)
!    CALL allocate_shared_int_0D( restart%ny, restart%wny)
!    CALL allocate_shared_int_0D( restart%nz, restart%wnz)
!    CALL allocate_shared_int_0D( restart%nt, restart%wnt)
!    IF (par%master) THEN
!      restart%netcdf%filename = filename_restart
!      CALL inquire_restart_file_temperature( restart)
!    END IF
!    CALL sync
!    
!    ! Allocate memory for raw data
!    CALL allocate_shared_dp_1D( restart%nx, restart%x,    restart%wx   )
!    CALL allocate_shared_dp_1D( restart%ny, restart%y,    restart%wy   )
!    CALL allocate_shared_dp_1D( restart%nz, restart%zeta, restart%wzeta)
!    CALL allocate_shared_dp_1D( restart%nt, restart%time, restart%wtime)
!    
!    CALL allocate_shared_dp_3D( restart%nx, restart%ny, restart%nz, restart%Ti,               restart%wTi              )
!  
!    ! Read data from input file
!    IF (par%master) CALL read_restart_file_temperature( restart, time_to_restart_from)
!    CALL sync
!    
!    ! Safety
!    CALL check_for_NaN_dp_3D( restart%Ti, 'restart%Ti', 'initialise_ice_temperature')
!    
!    ! Since we want data represented as [j,i] internally, transpose the data we just read.
!    CALL transpose_dp_3D( restart%Ti, restart%wTi)
!    
!    ! Map (transposed) raw data to the model grid
!    CALL map_square_to_square_cons_2nd_order_3D( restart%nx, restart%ny, restart%x, restart%y, grid%nx, grid%ny, grid%x, grid%y, restart%Ti, ice%Ti_a)
!    
!    ! Deallocate raw data
!    CALL deallocate_shared( restart%wnx              )
!    CALL deallocate_shared( restart%wny              )
!    CALL deallocate_shared( restart%wnz              )
!    CALL deallocate_shared( restart%wnt              )
!    CALL deallocate_shared( restart%wx               )
!    CALL deallocate_shared( restart%wy               )
!    CALL deallocate_shared( restart%wzeta            )
!    CALL deallocate_shared( restart%wtime            )
!    CALL deallocate_shared( restart%wTi              )
    
  END SUBROUTINE initialise_ice_temperature_restart
  
END MODULE thermodynamics_module
