MODULE ice_dynamics_module
  ! The actual ice dynamics model, based on the mixed SIA/SSA approach from ANICE.

  USE mpi
  USE configuration_module,        ONLY: dp, C
  USE parameters_module,           ONLY: seawater_density, ice_density, T0
  USE zeta_module,                 ONLY: vertical_integrate, vertical_average
  USE parallel_module,             ONLY: par, sync, &
                                         allocate_shared_memory_int_0D, allocate_shared_memory_dp_0D, &
                                         allocate_shared_memory_int_1D, allocate_shared_memory_dp_1D, &
                                         allocate_shared_memory_int_2D, allocate_shared_memory_dp_2D, &
                                         allocate_shared_memory_int_3D, allocate_shared_memory_dp_3D, &
                                         deallocate_shared_memory                                               
  USE data_types_module,           ONLY: type_mesh, type_ice_model, type_init_data_fields, type_SMB_model, type_BMB_model, type_model_region, &
                                         type_climate_model
  USE mesh_derivatives_module,     ONLY: GetMeshDerivatives, GetMeshCurvatures, GetMeshCurvaturesVertex, ApplyNeumannBoundary, &
                                         ApplyNeumannBoundary_3D, GetMeshDerivativesVertex_3D, GetNeighbourFunctions_masked, GetMeshDerivatives_masked_vertex_3D
  USE mesh_ArakawaC_module,        ONLY: GetAcMeshDerivatives, MapAaToAc, MapAaToAc_3D, MapAcToAa, MapAcToAa_3D
  USE general_ice_model_data_module, ONLY: ice_physical_properties, basal_yield_stress

  IMPLICIT NONE
  
CONTAINS  
  ! Update ice thickness from ice dynamic changes and SMB
  SUBROUTINE calculate_ice_thickness_change(ice, mesh, SMB, BMB, dt)
    ! Use the total ice velocities to update the ice thickness
    
    USE parameters_module, ONLY: ice_density
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp),                            INTENT(IN)    :: dt
    
    ! Local variables
    INTEGER                                            :: aci, vi, vj, cii, ci, cji, cj
    REAL(dp)                                           :: Upar
    REAL(dp)                                           :: dVi, Vi_out, Vi_in, Vi_available, rescale_factor
    REAL(dp), DIMENSION(mesh%nV)                       :: Vi_SMB
            
    ! Initialise at zero
    ice%dVi_in(  mesh%v1:mesh%v2, :) = 0._dp
    ice%dVi_out( mesh%v1:mesh%v2, :) = 0._dp
    CALL sync
    
    Vi_in          = 0._dp
    Vi_available   = 0._dp
    rescale_factor = 0._dp
    Vi_SMB         = 0._dp
        
    ! Calculate ice fluxes across all Aa vertex connections
    ! based on ice velocities calculated on Ac mesh
    ! =============================================
    
    DO aci = mesh%c1, mesh%c2
    
      ! The two Aa vertices connected by the Ac vertex
      vi = mesh%Aci(aci,1)
      vj = mesh%Aci(aci,2)
      
      ! Find their own respective connectivity, for storing the ice flux
      ci = 0
      cj = 0
      DO cii = 1, mesh%nC(vi)
        IF (mesh%C(vi,cii)==vj) THEN
          ci = cii
          EXIT
        END IF
      END DO
      DO cji = 1, mesh%nC(vj)
        IF (mesh%C(vj,cji)==vi) THEN
          cj = cji
          EXIT
        END IF
      END DO
            
      ! Calculate ice volume per year moving along connection from vi to vj as the product of:
      ! - width          (m   - determined by distance between adjacent triangle circumcenters)
      ! - ice thickness  (m   - at flux origin (i when Upar>0, else j))
      ! - ice velocity   (m/y - calculated at midpoint, using surface slope along connection)
      ! - time step      (y)
      Upar = ice%Up_ac(aci) + ice%Up_SSA_ac(aci)
      IF (Upar > 0._dp) THEN
        dVi = ice%Hi(vi) * Upar * mesh%Cw(vi,ci) * dt ! m3
      ELSE
        dVi = ice%Hi(vj) * Upar * mesh%Cw(vi,ci) * dt ! m3
      END IF
     ! dVi = (ice%Hi(vi) + ice%Hi(vj)) * 0.5_dp * Upar * mesh%Cw(vi,ci) * dt ! m3
      
      ! Keep track of ice fluxes across individual connections, to correct for
      ! negative ice thicknesses if necessary         
      ice%dVi_in( vi, ci) = -dVi ! m3
      ice%dVi_in( vj, cj) =  dVi ! m3
           
    END DO ! DO aci = mesh%c1, mesh%c2
    CALL sync
           
    
    ! Correct outfluxes for possible resulting negative ice thicknesses
    ! =================================================================
    
    Vi_SMB(mesh%v1:mesh%v2) = (SMB%SMB_year(mesh%v1:mesh%v2) + BMB%BMB(mesh%v1:mesh%v2))  * mesh%A(mesh%v1:mesh%v2) * dt! * ice_density / 1000._dp     ! m3 ice equivalent
    CALL sync
    
    DO vi = mesh%v1, mesh%v2
    
      ! Check how much ice is available for melting or removing (in m^3)
      Vi_available = mesh%A(vi) * ice%Hi(vi)
      
      dVi = SUM(ice%dVi_in(vi,:))
      
      Vi_in  = 0._dp
      Vi_out = 0._dp
      DO ci = 1, mesh%nC(vi)
        IF (ice%dVi_in(vi,ci) > 0._dp) THEN
          Vi_in  = Vi_in  + ice%dVi_in(vi,ci)
        ELSE
          Vi_out = Vi_out - ice%dVi_in(vi,ci)
        END IF
      END DO
      
      
      rescale_factor = 1._dp      
      
      ! If all the ice already present melts away, there can be no outflux.
      IF (-Vi_SMB(vi) >= Vi_available) THEN      
        ! All ice in this vertex melts, nothing remains to move around. Rescale outfluxes to zero.
        Vi_SMB(vi) = -Vi_available
        rescale_factor = 0._dp
      END IF
            
      ! If the total outflux exceeds the available ice plus SMB plus total influx, rescale outfluxes
      IF (Vi_out > Vi_available + Vi_SMB(vi)) THEN      
        ! Total outflux plus melt exceeds available ice volume. Rescale outfluxes to correct for this.
        rescale_factor = (Vi_available + Vi_SMB(vi)) / Vi_out        
      END IF
      
      ! Rescale ice outfluxes out of vi and into vi's neighbours
      IF (rescale_factor < 1._dp) THEN
        DO ci = 1, mesh%nC(vi)
          vj = mesh%C(vi,ci)
        
          IF (ice%dVi_in(vi,ci) < 0._dp) THEN
            ice%dVi_in(vi,ci) = ice%dVi_in(vi,ci) * rescale_factor
              
            DO cji = 1, mesh%nC(vj)
              IF (mesh%C(vj,cji)==vi) THEN
                ice%dVi_in(vj,cji) = -ice%dVi_in(vi,ci)
                EXIT
              END IF
            END DO
          
          END IF
        END DO
      END IF
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! Calculate change in ice thickness over time at every vertex
    ! ===========================================================
    
    ice%dHi_dt(mesh%v1:mesh%v2) = 0._dp ! m/y    
    DO vi = mesh%v1, mesh%v2
      dVi  = SUM(ice%dVi_in( vi,:))
      ice%dHi_dt(vi) = (dVi + Vi_SMB(vi)) / (mesh%A(vi) * dt)  
    END DO    
    CALL sync
    
    ! Manually set to zero in first timestep.
    IF (dt==0._dp) ice%dHi_dt(mesh%v1:mesh%v2) = 0._dp
    CALL sync
    
    ! Add dynamic ice thickness change and SMB to update the ice thickness
    ! ====================================================================
    
    ice%Hi(mesh%v1:mesh%v2) = ice%Hi(mesh%v1:mesh%v2) + (ice%dHi_dt(mesh%v1:mesh%v2) * dt)

    ! Set ice thickness to zero at the domain boundary
    DO vi = mesh%v1, mesh%v2
      IF (mesh%edge_index(vi) > 0) ice%Hi(vi) = 0._dp
    END DO
    CALL sync
    
    ! MISMIP3D: Neumann boundary everywhere except east side
    IF (C%do_mismip3d_experiment) THEN
      CALL ApplyNeumannBoundary(mesh, ice%Hi)
      DO vi = mesh%v1, mesh%v2
        IF (mesh%edge_index(vi) == 2 .OR. mesh%edge_index(vi) == 3 .OR. mesh%edge_index(vi) == 4) ice%Hi(vi) = 0._dp
      END DO
      CALL sync
    END IF
    
  END SUBROUTINE calculate_ice_thickness_change
  
  ! SIA solver
  SUBROUTINE solve_SIA(ice, mesh)
    ! Calculate ice velocities using the SIA
    
    USE parameters_module, ONLY: n_flow, m_flow, ice_density, seawater_density, grav, A_sliding
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    INTEGER                                            :: ci, k
    REAL(dp)                                           :: D_0
    REAL(dp), DIMENSION(C%NZ)                          :: D_deformation
    REAL(dp)                                           :: D_uv_2D 
    REAL(dp), DIMENSION(C%nZ)                          :: D_SIA_3D_prof    
    REAL(dp), PARAMETER                                :: D_uv_3D_cutoff = -1E5_dp  
        
    ice%D_SIA_ac(   mesh%c1:mesh%c2  ) = 0._dp
    ice%D_SIA_3D_ac(mesh%c1:mesh%c2,:) = 0._dp
    ice%Ux_ac(      mesh%c1:mesh%c2  ) = 0._dp
    ice%Uy_ac(      mesh%c1:mesh%c2  ) = 0._dp
    ice%Up_ac(      mesh%c1:mesh%c2  ) = 0._dp
    ice%Uo_ac(      mesh%c1:mesh%c2  ) = 0._dp
    CALL sync
        
    ! From the ANICE subroutine "calculate_D_uv_3D"
    
    DO ci = mesh%c1, mesh%c2
    
      IF (ice%mask_sheet_ac(ci)==1) THEN     
        
       D_0      = (ice_density * grav * ice%Hi_ac(ci))**n_flow * ((ice%dHs_dp_ac(ci)**2 + ice%dHs_do_ac(ci)**2))**((n_flow - 1._dp) / 2._dp)       
       D_deformation = vertical_integrate(C%m_enh_sia * ice%A_flow_ac(ci,:) * C%zeta**n_flow)       
       D_deformation = 2._dp * ice%Hi_ac(ci) * D_deformation
       
       ice%D_SIA_3D_ac(ci,:) = D_0 * D_deformation

      END IF ! End: Considering only grounded points
      
      ! Check for very large D_uv_3D's, causing very large velocities RESET --DIRTY
      DO k = 1, C%NZ
        IF (ice%D_SIA_3D_ac(ci,k) < D_uv_3D_cutoff) THEN
          ice%D_SIA_3D_ac(ci,k) = D_uv_3D_cutoff          
        END IF        
      END DO

    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
    
    ! From the ANICE subroutine "sheet_velocity_2D"                                                    
    DO ci = mesh%c1, mesh%c2
      IF (ice%mask_sheet_ac(ci)==1) THEN
        D_SIA_3D_prof = ice%D_SIA_3D_ac(ci,:)
        D_uv_2D = vertical_average(D_SIA_3D_prof)
       
        ice%D_SIA_ac(ci) = ice%Hi_ac(ci) * D_uv_2D     ! See equation (8.5)
        ice%Ux_ac(   ci) = D_uv_2D * ice%dHs_dx_ac(ci) ! See equation (7.39)
        ice%Uy_ac(   ci) = D_uv_2D * ice%dHs_dy_ac(ci) ! See equation (7.40)
        ice%Up_ac(   ci) = D_uv_2D * ice%dHs_dp_ac(ci)
        ice%Uo_ac(   ci) = D_uv_2D * ice%dHs_do_ac(ci)        
      END IF      
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
    
    ! Map data to Aa mesh for writing to output (not actually used anywhere...)    
    CALL MapAcToAa(mesh, ice%Ux_ac, ice%U_SIA)
    CALL MapAcToAa(mesh, ice%Uy_ac, ice%V_SIA)
      
  END SUBROUTINE solve_SIA
  
  ! 3D SIA solver (used only for thermodynamics)
  SUBROUTINE solve_SIA_3D(ice, mesh)
    
    USE parameters_module, ONLY: n_flow, m_flow, ice_density, seawater_density, grav, A_sliding
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
      
    ! Local variables:    
    REAL(dp)                                           :: D_0
    REAL(dp), DIMENSION(C%nZ)                          :: D_deformation
    REAL(dp), DIMENSION(C%nZ)                          :: D_SIA_3D
    REAL(dp), PARAMETER                                :: D_uv_3D_cutoff = -1E5_dp 
    REAL(dp)                                           :: dHb_dx
    REAL(dp)                                           :: dHb_dy
    INTEGER                                            :: vi, k
    REAL(dp)                                           :: dUdx_k, dUdy_k, dUdx_kp1, dUdy_kp1
    REAL(dp)                                           :: dVdx_k, dVdy_k, dVdx_kp1, dVdy_kp1
    REAL(dp)                                           :: w1, w2, w3, w4
                                                     
    ! Setting the shelf values to zero:
    ice%U_3D( mesh%v1:mesh%v2,:) = 0._dp
    ice%V_3D( mesh%v1:mesh%v2,:) = 0._dp
    ice%W_3D( mesh%v1:mesh%v2,:) = 0._dp
    CALL sync

    ! Direct calculation of U and V for the ice sheet:
    DO vi = mesh%v1, mesh%v2
    
      ice%U_3D(vi,:) = ice%U_SSA(vi)
      ice%V_3D(vi,:) = ice%V_SSA(vi)
    
      IF (ice%mask(vi) == C%type_shelf) CYCLE
      IF (ice%Hi(vi) == 0._dp) CYCLE
        
      D_0           = (ice_density * grav * ice%Hi(vi))**n_flow * ((ice%dHs_dx(vi)**2 + ice%dHs_dy(vi)**2))**((n_flow - 1._dp) / 2._dp)
      D_deformation = vertical_integrate(C%m_enh_sia * ice%A_flow(vi,:) * C%zeta**n_flow)
      D_deformation = 2._dp * ice%Hi(vi) * D_deformation
      D_SIA_3D      = MAX(D_0 * D_deformation, D_uv_3D_cutoff)
       
      ! To be completely consistent, only add the SSA basal velocities when calculated:
      ice%U_3D(vi,:)    = D_SIA_3D * ice%dHs_dx(vi) + ice%U_SSA(vi)
      ice%V_3D(vi,:)    = D_SIA_3D * ice%dHs_dy(vi) + ice%V_SSA(vi)
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
        
    CALL ApplyNeumannBoundary_3D( mesh, ice%U_3D, C%nZ)
    CALL ApplyNeumannBoundary_3D( mesh, ice%V_3D, C%nZ)
    
    ! Calculate masked neighbour functions for one-sided surface slopes over the grounding line
    CALL GetNeighbourFunctions_masked(mesh, ice%mask_sheet)
    CALL sync

    DO vi = mesh%v1, mesh%v2
    
      IF (mesh%edge_index(vi) > 0) CYCLE
      IF (ice%mask_sheet(vi) == 0) CYCLE
      
      dHb_dx = ice%dHs_dx(vi) - ice%dHi_dx(vi)
      dHb_dy = ice%dHs_dy(vi) - ice%dHi_dy(vi)   
      
      ice%W_3D(vi,C%NZ) = ice%dHb_dt(vi) + ice%U_3D(vi,C%NZ) * dHb_dx + ice%V_3D(vi,C%NZ) * dHb_dy 
                           
      ! The integrant is calculated half way the layer of integration at k+1/2. This integrant is multiplied with the layer thickness and added to the integral
      ! of all layers below, giving the integral up to and including this layer:
      DO k = C%NZ - 1, 1, -1
!        IF (ice%mask(vi) == C%type_sheet) THEN
          CALL GetMeshDerivativesVertex_3D( mesh, ice%U_3D, dUdx_k,   dUdy_k,   vi, k)
          CALL GetMeshDerivativesVertex_3D( mesh, ice%U_3D, dUdx_kp1, dUdy_kp1, vi, k+1)
          CALL GetMeshDerivativesVertex_3D( mesh, ice%V_3D, dVdx_k,   dVdy_k,   vi, k)
          CALL GetMeshDerivativesVertex_3D( mesh, ice%V_3D, dVdx_kp1, dVdy_kp1, vi, k+1)                         
!        ELSE
!          ! One-sided differencing over the grounding line
!          CALL GetMeshDerivatives_masked_vertex_3D( mesh, ice%U_3D, dUdx_k,   dUdy_k,   vi, k)
!          CALL GetMeshDerivatives_masked_vertex_3D( mesh, ice%U_3D, dUdx_kp1, dUdy_kp1, vi, k+1)
!          CALL GetMeshDerivatives_masked_vertex_3D( mesh, ice%V_3D, dVdx_k,   dVdy_k,   vi, k)
!          CALL GetMeshDerivatives_masked_vertex_3D( mesh, ice%V_3D, dVdx_kp1, dVdy_kp1, vi, k+1)
!        END IF
        w1 = (dUdx_k + dUdx_kp1) / 2._dp
        w2 = (dVdy_k + dVdy_kp1) / 2._dp   
        
        w3              = ((ice%dHs_dx(vi) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * ice%dHi_dx(vi)) / MAX(0.1_dp, ice%Hi(vi))) *   &
                          ((ice%U_3D(vi,k+1) - ice%U_3D(vi,k)) / (C%zeta(k+1) - C%zeta(k)))
        w4              = ((ice%dHs_dy(vi) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * ice%dHi_dy(vi)) / MAX(0.1_dp, ice%Hi(vi))) *   &
                          ((ice%V_3D(vi,k+1) - ice%V_3D(vi,k)) / (C%zeta(k+1) - C%zeta(k)))

        ice%W_3D(vi,k) = ice%W_3D(vi,k+1) - ice%Hi(vi) * (w1 + w2 + w3 + w4) * (C%zeta(k+1) - C%zeta(k))
        
      END DO

    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    CALL ApplyNeumannBoundary_3D( mesh, ice%W_3D, C%nZ)
    CALL sync

  END SUBROUTINE solve_SIA_3D
  
  ! SSA solver
  SUBROUTINE solve_SSA(ice, mesh)
    ! Calculate ice velocities using the SSA
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    IF     (C%choice_sliding_law == 'Coulomb') THEN
      CALL solve_SSA_Coulomb_free(ice, mesh)
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      CALL solve_SSA_Weertman_free(ice, mesh)
    ELSE
      IF (par%master) WRITE(0,*) ' ERROR: choice_sliding_law can only be "Coulomb" or "Weertman"!'
      STOP
    END IF
    
  END SUBROUTINE solve_SSA
  SUBROUTINE solve_SSA_Coulomb_free(ice, mesh)
    ! Calculate ice velocities using the SSA
    
    USE parameters_module, ONLY : n_flow, q_plastic, u_threshold, ice_density, grav, epsilon_sq_0, delta_v
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    REAL(dp), PARAMETER                                :: residual_epsilon = 2.5_dp         ! Stop criterium. The residual should drop below this value everywhere
    INTEGER,  PARAMETER                                :: outer_loop_n = 6
    INTEGER,  PARAMETER                                :: inner_loop_n = 10000
    REAL(dp), PARAMETER                                :: omega = 1.1_dp
    INTEGER                                            :: outer_loop_i
    INTEGER                                            :: inner_loop_i, nit
    
    INTEGER                                            :: vi, ci, vc, vj
    REAL(dp)                                           :: Uxxi, Uxyi, Uyyi, Vxxi, Vxyi, Vyyi
    REAL(dp)                                           :: R_max, eu_min
    REAL(dp)                                           :: sumUc, sumVc
    INTEGER                                            :: i, ierr, status(MPI_STATUS_SIZE)
    REAL(dp)                                           :: Ux, Uy, U
    
    ! If we're doing one of the EISMINT experiments, don't solve the SSA.
    IF (C%do_eismint_experiment) THEN
      ice%U_SSA(mesh%v1:mesh%v2) = 0._dp
      ice%V_SSA(mesh%v1:mesh%v2) = 0._dp
      CALL sync
      RETURN
    END IF
    
    ! If there's no grounded ice anywhere, don't solve the SSA.
    IF (SUM(ice%mask_sheet)==0) RETURN
    
    ! Calculate the basal yield stress
    CALL basal_yield_stress(ice, mesh)
    
    ! Solve the SSA
    ! Same approach as an ANICE: an outer loop where ice viscosity is updated,
    ! and an inner loop that uses Successive Over-Relaxation (SOR) to solve the two coupled PDE's for the given viscosity.
                  
    DO outer_loop_i = 1, outer_loop_n
    
      !IF (par%master) WRITE(0,'(A,I1)') ' Outer loop ', outer_loop_i
      
      ! Calculate the terms that don't change in the inner loop: the non-linear term nu,
      ! the centre coefficients eu_i & ev_i, and the right hand sides rhs_x & rhs_y
      
      ! nu requires the spatial derivatives of U and V, easier to use the existing routines for those
      CALL GetMeshDerivatives( mesh, ice%U_SSA, ice%dUdx, ice%dUdy)
      CALL GetMeshDerivatives( mesh, ice%V_SSA, ice%dVdx, ice%dVdy)  
      
      DO vi = mesh%v1, mesh%v2

        ice%nu(vi) = (C%m_enh_ssa * 0.5_dp * ice%A_flow_mean(vi))**(-1._dp / n_flow) * (ice%dUdx(vi)**2 + ice%dVdy(vi)**2 + ice%dUdx(vi) * ice%dVdy(vi) &
                        + 0.25_dp * (ice%dUdy(vi) + ice%dVdx(vi))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow))

        ! Determine beta_base, for basal stress (beta in eqn (27) B&B, 2009), here equivalent to the sliding parameter (A**-1/m)
        ice%beta_base(vi) = ice%tau_yield(vi) * ( (delta_v**2 + ice%U_SSA(vi)**2 + ice%V_SSA(vi)**2)**(0.5_dp * (q_plastic-1._dp)) ) & 
                         / (u_threshold**q_plastic) 

        ! The right hand side of the PDE, the gravitational driving stress
        ice%rhs_x(vi) = (ice_density * grav * ice%dHs_dx_shelf(vi)) / ice%nu(vi)
        ice%rhs_y(vi) = (ice_density * grav * ice%dHs_dy_shelf(vi)) / ice%nu(vi)
       
        ! The coefficients eu_i and ev_i, which don't change during the inner loop iteration
        ice%eu_i(vi) = (4._dp * mesh%Nxx(vi,mesh%nC(vi)+1) + mesh%Nyy(vi,mesh%nC(vi)+1)) - ice%beta_base(vi) / ice%nu(vi)
        ice%ev_i(vi) = (4._dp * mesh%Nyy(vi,mesh%nC(vi)+1) + mesh%Nxx(vi,mesh%nC(vi)+1)) - ice%beta_base(vi) / ice%nu(vi)

      END DO ! DO vi = mesh%v1, mesh%v2
      CALL sync      
      
      ! Inner loop: use SOR to solve U and V for the given right-hand side of the PDE's
      
      nit = 0
      DO inner_loop_i = 1, inner_loop_n
        nit = nit+1
        
        ! Go through one iteration of the SOR scheme
        ! ==========================================
        
        R_max  = 0._dp
        eu_min = 1.0E10_dp
        
        DO vi = mesh%v1, mesh%v2
          IF (mesh%edge_index(vi)>0) CYCLE ! Don't update edge vertices, those will be treated by the boundary conditions
      
          ! Calculate Uxx, Uyy, etc. with partially updated values pf U and V according to Gauss-Seidler
          CALL GetMeshCurvaturesVertex( mesh, ice%U_SSA, Uxxi, Uxyi, Uyyi, vi)
          CALL GetMeshCurvaturesVertex( mesh, ice%V_SSA, Vxxi, Vxyi, Vyyi, vi)
          
          ! The sum terms in the equation
          sumUc = 0._dp
          sumVc = 0._dp                    
          DO ci = 1, mesh%nC(vi)
            vc = mesh%C(vi,ci)
            sumUc = sumUc + ice%U_SSA(vc) * (4._dp * mesh%Nxx(vi,ci) + mesh%Nyy(vi,ci))
            sumVc = sumVc + ice%V_SSA(vc) * (4._dp * mesh%Nyy(vi,ci) + mesh%Nxx(vi,ci))
          END DO
          
          ! Calculate the rest terms
          ice%Rx(vi) = sumUc + (ice%eu_i(vi) * ice%U_SSA(vi)) + (3._dp * Vxyi) - ice%rhs_x(vi)
          ice%Ry(vi) = sumVc + (ice%ev_i(vi) * ice%V_SSA(vi)) + (3._dp * Uxyi) - ice%rhs_y(vi)
          
          ! Update U and V with SOR
          IF (ice%eu_i(vi) == 0._dp) THEN
            ice%U_SSA(vi) = ice%U_SSA(vi)
          ELSE 
            !ice%U_SSA(vi) = ice%U_SSA(vi) - mesh%omega(vi) * Rx / ice%eu_i(vi)
            ice%U_SSA(vi) = ice%U_SSA(vi) - omega * ice%Rx(vi) / ice%eu_i(vi)
            
            R_max  = MAX( R_max,  ABS(ice%Rx(  vi)))
            eu_min = MIN( eu_min, ABS(ice%eu_i(vi)))
          END IF  
          IF (ice%ev_i(vi) == 0._dp) THEN
            ice%V_SSA(vi) = ice%V_SSA(vi)
          ELSE 
            !ice%V_SSA(vi) = ice%V_SSA(vi) - mesh%omega(vi) * Ry / ice%ev_i(vi)
            ice%V_SSA(vi) = ice%V_SSA(vi) - omega * ice%Ry(vi) / ice%ev_i(vi)
            
            R_max  = MAX( R_max,  ABS(ice%Rx(  vi)))
            eu_min = MIN( eu_min, ABS(ice%eu_i(vi)))
          END IF
        
        END DO ! DO vi = mesh%v1, mesh%v2 
        CALL sync
        
        ! Apply Neumann boundary conditions: make sure Ux, Uy, Vx and Vy are zero at domain boundary
        CALL ApplyNeumannBoundary( mesh, ice%U_SSA)
        CALL ApplyNeumannBoundary( mesh, ice%V_SSA)

        ! Check to see if a stable solution is reached  
        ! ============================================
        
        ! Find maximum R and minimum eu across all processes (use native MPI routines, MUCH faster than using SEND/RECEIVE loops!)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, R_max,  1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, eu_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    
       ! IF (par%master) WRITE(0,'(A,I6,A,E8.2)') '  Inner loop ', inner_loop_i, ': epsilon = ', R_max / eu_min
        
        IF (R_max > 1000._dp) THEN
          WRITE(0,*) '    calculate_ice_velocities_SSA - ERROR: numerical instability in SSA solver! R_max > 1000'
          STOP
        END IF
        
        ! Determine if a stable solution is reached
        IF ( R_max / eu_min < residual_epsilon) THEN
          EXIT        
        ELSE IF(inner_loop_i == inner_loop_n) THEN
          ! When iteratation does not converge set Us and Vs to zero and redo the iteration
         ! IF (par%master) WRITE(0,*) '    calculate_ice_velocities_SSA - WARNING: hard reset of SSA velocities!'
          ice%U_SSA = 0._dp
          ice%V_SSA = 0._dp         
        END IF
                
        ! DENK DROM
       ! EXIT
      
      END DO ! inner loop
    
     ! IF (par%master) WRITE(0,*) '   SSA - inner loop solved in ', nit, ' iterations'
      
      ! DENK DROM
     ! EXIT
      
    END DO ! outer loop
    
    ! Map data to Ac mesh for ice thickness calculation
    CALL MapAaToAc(mesh, ice%U_SSA, ice%Ux_SSA_ac)
    CALL MapAaToAc(mesh, ice%V_SSA, ice%Uy_SSA_ac)
    
    ! Get parallel and orthogonal components of SSA velocities for ice thickness calculation
    DO ci = mesh%c1, mesh%c2
      vi = mesh%Aci(ci,1)
      vj = mesh%Aci(ci,2)
      Ux = mesh%V(vj,1)-mesh%V(vi,1)
      Uy = mesh%V(vj,2)-mesh%V(vi,2)
      U  = SQRT(Ux**2+Uy**2)
      
      ice%Up_SSA_ac(ci) = ice%Ux_SSA_ac(ci) * Ux/U + ice%Uy_SSA_ac(ci) * Uy/U
      ice%Uo_SSA_ac(ci) = ice%Uy_SSA_ac(ci) * Ux/U - ice%Ux_SSA_ac(ci) * Uy/U
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
  
  END SUBROUTINE solve_SSA_Coulomb_free
  SUBROUTINE solve_SSA_Weertman_free(ice, mesh)
    ! Calculate ice velocities using the SSA
    
    USE parameters_module, ONLY : n_flow, q_plastic, u_threshold, ice_density, grav, epsilon_sq_0, delta_v
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    REAL(dp), PARAMETER                                :: residual_epsilon = 2.5_dp         ! Stop criterium. The residual should drop below this value everywhere
    INTEGER,  PARAMETER                                :: outer_loop_n = 6
    INTEGER,  PARAMETER                                :: inner_loop_n = 10000
    REAL(dp), PARAMETER                                :: omega = 1.0_dp
    INTEGER                                            :: outer_loop_i
    INTEGER                                            :: inner_loop_i, nit
    
    INTEGER                                            :: vi, ci, vc, vj
    REAL(dp)                                           :: sliding_term_x, sliding_term_y
    REAL(dp)                                           :: Uxxi, Uxyi, Uyyi, Vxxi, Vxyi, Vyyi
    REAL(dp)                                           :: R_max, eu_min
    REAL(dp)                                           :: sumUc, sumVc
    INTEGER                                            :: i, ierr, status(MPI_STATUS_SIZE)
    REAL(dp)                                           :: Ux, Uy, U
    
    ! If we're doing one of the EISMINT experiments, don't solve the SSA.
    IF (C%do_eismint_experiment) THEN
      ice%U_SSA(mesh%v1:mesh%v2) = 0._dp
      ice%V_SSA(mesh%v1:mesh%v2) = 0._dp
      CALL sync
      RETURN
    END IF
    
    ! If there's no grounded ice anywhere, don't solve the SSA.
    IF (SUM(ice%mask_sheet)==0) RETURN
    
    ! Calculate the basal yield stress
    CALL basal_yield_stress(ice, mesh)
    
    ! Solve the SSA
    ! Same approach as an ANICE: an outer loop where ice viscosity is updated,
    ! and an inner loop that uses Successive Over-Relaxation (SOR) to solve the two coupled PDE's for the given viscosity.
            
    DO outer_loop_i = 1, outer_loop_n
    
      !IF (par%master) WRITE(0,'(A,I1)') ' Outer loop ', outer_loop_i
      
      ! Calculate the terms that don't change in the inner loop: the non-linear term nu,
      ! the centre coefficients eu_i & ev_i, and the right hand sides rhs_x & rhs_y
      
      ! nu requires the spatial derivatives of U and V, easier to use the existing routines for those
      CALL GetMeshDerivatives( mesh, ice%U_SSA, ice%dUdx, ice%dUdy)
      CALL GetMeshDerivatives( mesh, ice%V_SSA, ice%dVdx, ice%dVdy)  
      
      DO vi = mesh%v1, mesh%v2

        ice%nu(vi) = (C%m_enh_ssa * 0.5_dp * ice%A_flow_mean(vi))**(-1._dp / n_flow) * (ice%dUdx(vi)**2 + ice%dVdy(vi)**2 + ice%dUdx(vi) * ice%dVdy(vi) &
                        + 0.25_dp * (ice%dUdy(vi) + ice%dVdx(vi))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow))

        ! The right hand side of the PDE, the gravitational driving stress
        ice%rhs_x(vi) = (ice_density * grav * ice%dHs_dx_shelf(vi)) / ice%nu(vi)
        ice%rhs_y(vi) = (ice_density * grav * ice%dHs_dy_shelf(vi)) / ice%nu(vi)
        
        IF (ice%mask_sheet(i)==1) THEN
          sliding_term_x = (C%C_sliding / (C%sec_per_year ** C%m_sliding)) * ((ABS(ice%U_SSA(i)) + delta_v) ** (C%m_sliding - 1._dp)) / (ice%nu(i) * ice%Hi(i))
          sliding_term_y = (C%C_sliding / (C%sec_per_year ** C%m_sliding)) * ((ABS(ice%V_SSA(i)) + delta_v) ** (C%m_sliding - 1._dp)) / (ice%nu(i) * ice%Hi(i))
        ELSE
          sliding_term_x = 0._dp
          sliding_term_y = 0._dp
        END IF
       
        ! The coefficients eu_i and ev_i, which don't change during the inner loop iteration
        ice%eu_i(vi) = (4._dp * mesh%Nxx(vi,mesh%nC(vi)+1) + mesh%Nyy(vi,mesh%nC(vi)+1)) - sliding_term_x
        ice%ev_i(vi) = (4._dp * mesh%Nyy(vi,mesh%nC(vi)+1) + mesh%Nxx(vi,mesh%nC(vi)+1)) - sliding_term_y

      END DO ! DO vi = mesh%v1, mesh%v2
      CALL sync
      
      
      ! Inner loop: use SOR to solve U and V for the given right-hand side of the PDE's
      nit = 0
      DO inner_loop_i = 1, inner_loop_n
        nit = nit+1
        
        ! Go through one iteration of the SOR scheme
        ! ==========================================
        
        R_max  = 0._dp
        eu_min = 1.0E10_dp
        
        DO vi = mesh%v1, mesh%v2
          IF (mesh%edge_index(vi)>0) CYCLE ! Don't update edge vertices, those will be treated by the boundary conditions
      
          ! Calculate Uxx, Uyy, etc. with partially updated values pf U and V according to Gauss-Seidler
          CALL GetMeshCurvaturesVertex( mesh, ice%U_SSA, Uxxi, Uxyi, Uyyi, vi)
          CALL GetMeshCurvaturesVertex( mesh, ice%V_SSA, Vxxi, Vxyi, Vyyi, vi)
          
          ! The sum terms in the equation
          sumUc = 0._dp
          sumVc = 0._dp                    
          DO ci = 1, mesh%nC(vi)
            vc = mesh%C(vi,ci)
            sumUc = sumUc + ice%U_SSA(vc) * (4._dp * mesh%Nxx(vi,ci) + mesh%Nyy(vi,ci))
            sumVc = sumVc + ice%V_SSA(vc) * (4._dp * mesh%Nyy(vi,ci) + mesh%Nxx(vi,ci))
          END DO
          
          ! Calculate the rest terms
          ice%Rx(vi) = sumUc + (ice%eu_i(vi) * ice%U_SSA(vi)) + (3._dp * Vxyi) - ice%rhs_x(vi)
          ice%Ry(vi) = sumVc + (ice%ev_i(vi) * ice%V_SSA(vi)) + (3._dp * Uxyi) - ice%rhs_y(vi)
          
          ! Update U and V with SOR
          IF (ice%eu_i(vi) == 0._dp) THEN
            ice%U_SSA(vi) = ice%U_SSA(vi)
          ELSE 
            !ice%U_SSA(vi) = ice%U_SSA(vi) - mesh%omega(vi) * Rx / ice%eu_i(vi)
            ice%U_SSA(vi) = ice%U_SSA(vi) - omega * ice%Rx(vi) / ice%eu_i(vi)
            
            R_max  = MAX( R_max,  ABS(ice%Rx(  vi)))
            eu_min = MIN( eu_min, ABS(ice%eu_i(vi)))
          END IF  
          IF (ice%ev_i(vi) == 0._dp) THEN
            ice%V_SSA(vi) = ice%V_SSA(vi)
          ELSE 
            !ice%V_SSA(vi) = ice%V_SSA(vi) - mesh%omega(vi) * Ry / ice%ev_i(vi)
            ice%V_SSA(vi) = ice%V_SSA(vi) - omega * ice%Ry(vi) / ice%ev_i(vi)
            
            R_max  = MAX( R_max,  ABS(ice%Rx(  vi)))
            eu_min = MIN( eu_min, ABS(ice%eu_i(vi)))
          END IF
        
        END DO ! DO vi = mesh%v1, mesh%v2
        CALL sync
        
        ! Apply Neumann boundary conditions: make sure Ux, Uy, Vx and Vy are zero at domain boundary
        CALL ApplyNeumannBoundary( mesh, ice%U_SSA)
        CALL ApplyNeumannBoundary( mesh, ice%V_SSA)

        ! Check to see if a stable solution is reached  
        ! ============================================
        
        ! Find maximum R and minimum eu across all processes (use native MPI routines, MUCH faster than using SEND/RECEIVE loops!)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, R_max,  1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, eu_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    
       ! IF (par%master) WRITE(0,'(A,I6,A,E8.2)') '  Inner loop ', inner_loop_i, ': epsilon = ', R_max / eu_min
        
        IF (R_max > 1000._dp) THEN
          WRITE(0,*) '    calculate_ice_velocities_SSA - ERROR: numerical instability in SSA solver! R_max > 1000'
          STOP
        END IF
        
        ! Determine if a stable solution is reached
        IF ( R_max / eu_min < residual_epsilon) THEN
          EXIT        
        ELSE IF(inner_loop_i == inner_loop_n) THEN
          ! When iteratation does not converge set Us and Vs to zero and redo the iteration
         ! IF (par%master) WRITE(0,*) '    calculate_ice_velocities_SSA - WARNING: hard reset of SSA velocities!'
         ! ice%U_SSA = 0._dp
         ! ice%V_SSA = 0._dp         
        END IF
        
        ! DENK DROM
       ! EXIT
      
      END DO ! inner loop
    
      !IF (par%master) WRITE(0,*) '   SSA - inner loop solved in ', nit, ' iterations - max velocity = ', MAXVAL(SQRT(ice%U_SSA**2 + ice%V_SSA**2))
      
      ! DENK DROM
     ! EXIT
      
    END DO ! outer loop
    
  !  IF (par%master) WRITE(0,*) '   SSA max velocity = ', MAXVAL(SQRT(ice%U_SSA**2 + ice%V_SSA**2))
    
    ! Map data to Ac mesh for ice thickness calculation
    CALL MapAaToAc(mesh, ice%U_SSA, ice%Ux_SSA_ac)
    CALL MapAaToAc(mesh, ice%V_SSA, ice%Uy_SSA_ac)
    
    ! Get parallel and orthogonal components of SSA velocities for ice thickness calculation
    DO ci = mesh%c1, mesh%c2
      vi = mesh%Aci(ci,1)
      vj = mesh%Aci(ci,2)
      Ux = mesh%V(vj,1)-mesh%V(vi,1)
      Uy = mesh%V(vj,2)-mesh%V(vi,2)
      U  = SQRT(Ux**2+Uy**2)
      
      ice%Up_SSA_ac(ci) = ice%Ux_SSA_ac(ci) * Ux/U + ice%Uy_SSA_ac(ci) * Uy/U
      ice%Uo_SSA_ac(ci) = ice%Uy_SSA_ac(ci) * Ux/U - ice%Ux_SSA_ac(ci) * Uy/U
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
  
  END SUBROUTINE solve_SSA_Weertman_free
  
  ! Combine SSA and SIA velocities
  SUBROUTINE combine_SIA_SSA_velocities(ice, mesh)

    ! Input variables:
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    INTEGER                                            :: vi
    
    ! On the Aa mesh, for writing to output
    DO vi = mesh%v1, mesh%v2
      IF (ice%Hi(vi) > 0._dp) THEN
        ice%Us(vi) = ice%U_SIA(vi) + ice%U_SSA(vi)
        ice%Vs(vi) = ice%V_SIA(vi) + ice%V_SSA(vi)
      ELSE
        ice%Us(vi) = 0._dp
        ice%Vs(vi) = 0._dp
      END IF
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
  END SUBROUTINE combine_SIA_SSA_velocities
  
  ! Initialisation
  SUBROUTINE initialise_ice_model(ice, mesh, init, filetype_init)
    ! Allocate shared memory for all the data fields of the ice dynamical module, and
    ! initialise some of them    
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
    CHARACTER(LEN=4),                    INTENT(IN)    :: filetype_init
    
    ! Allocate shared memory
    CALL allocate_ice_model(ice, mesh)
        
    ! Initialise with data from initial file    
    ice%Hi(       mesh%v1:mesh%v2) = init%Hi(mesh%v1:mesh%v2)
    ice%Hb(       mesh%v1:mesh%v2) = init%Hb(mesh%v1:mesh%v2)
    ice%sealevel( mesh%v1:mesh%v2) = 0._dp
    ice%Hs(       mesh%v1:mesh%v2) = MAX(ice%sealevel(mesh%v1:mesh%v2), ice%Hb(mesh%v1:mesh%v2) + ice%Hi(mesh%v1:mesh%v2))
    CALL sync
    
    IF (filetype_init == 'cart') THEN
      ice%U_SSA( mesh%v1:mesh%v2) = 0._dp
      ice%V_SSA( mesh%v1:mesh%v2) = 0._dp
    ELSE
      ice%U_SSA( mesh%v1:mesh%v2) = init%U_SSA( mesh%v1:mesh%v2)
      ice%V_SSA( mesh%v1:mesh%v2) = init%V_SSA( mesh%v1:mesh%v2)
    END IF
    CALL sync
    
  END SUBROUTINE initialise_ice_model
  SUBROUTINE allocate_ice_model(ice, mesh)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.    
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
      
    ! Allocate memory
    ! ===============
    
    ! Basic data - ice thickness, bedrock height, surface height, mask, 3D ice velocities and temperature
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Hi,                  ice%wHi                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Hb,                  ice%wHb                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Hs,                  ice%wHs                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%sealevel,            ice%wsealevel           )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Us,                  ice%wUs                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Vs,                  ice%wVs                 )
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%Ti,                  ice%wTi                 )  
    
    ! Different masks
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask,                ice%wmask               )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_land,           ice%wmask_land          )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_ocean,          ice%wmask_ocean         )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_lake,           ice%wmask_lake          )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_ice,            ice%wmask_ice           )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_sheet,          ice%wmask_sheet         )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_shelf,          ice%wmask_shelf         )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_coast,          ice%wmask_coast         )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_margin,         ice%wmask_margin        )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_groundingline,  ice%wmask_groundingline )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_calvingfront,   ice%wmask_calvingfront  )
    
    ! Ice physical properties
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%Ti_pmp,              ice%wTi_pmp             )
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%A_flow,              ice%wA_flow             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%A_flow_mean,         ice%wA_flow_mean        )
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%Cpi,                 ice%wCpi                ) 
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%Ki,                  ice%wKi                 ) 
    
    ! Spatial derivatives and curvatures
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHi_dx,              ice%wdHi_dx             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHi_dy,              ice%wdHi_dy             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHi_dt,              ice%wdHi_dt             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHb_dx,              ice%wdHb_dx             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHb_dy,              ice%wdHb_dy             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHb_dt,              ice%wdHb_dt             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHs_dx,              ice%wdHs_dx             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHs_dy,              ice%wdHs_dy             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHs_dt,              ice%wdHs_dt             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHs_dx_shelf,        ice%wdHs_dx_shelf       )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHs_dy_shelf,        ice%wdHs_dy_shelf       )
    
    ! Zeta derivatives
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%dzeta_dt,            ice%wdzeta_dt           )
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%dzeta_dx,            ice%wdzeta_dx           )
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%dzeta_dy,            ice%wdzeta_dy           )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dzeta_dz,            ice%wdzeta_dz           )
    
    ! Ice dynamics - SIA
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%Hi_ac,               ice%wHi_ac              )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%Hb_ac,               ice%wHb_ac              )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%Hs_ac,               ice%wHs_ac              )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%sealevel_ac,         ice%wsealevel_ac        )
    CALL allocate_shared_memory_int_1D( mesh%nAC,       ice%mask_ac,             ice%wmask_ac            )
    CALL allocate_shared_memory_int_1D( mesh%nAC,       ice%mask_ice_ac,         ice%wmask_ice_ac        )
    CALL allocate_shared_memory_int_1D( mesh%nAC,       ice%mask_sheet_ac,       ice%wmask_sheet_ac      )
    CALL allocate_shared_memory_int_1D( mesh%nAC,       ice%mask_shelf_ac,       ice%wmask_shelf_ac      )
    CALL allocate_shared_memory_int_1D( mesh%nAC,       ice%mask_ocean_ac,       ice%wmask_ocean_ac      )
    CALL allocate_shared_memory_int_1D( mesh%nAC,       ice%mask_groundingline_ac, ice%wmask_groundingline_ac)
    CALL allocate_shared_memory_int_1D( mesh%nAC,       ice%mask_calvingfront_ac,  ice%wmask_calvingfront_ac )
    CALL allocate_shared_memory_dp_2D(  mesh%nAC, C%nZ, ice%Ti_ac,               ice%wTi_ac              )
    CALL allocate_shared_memory_dp_2D(  mesh%nAC, C%nZ, ice%Ti_pmp_ac,           ice%wTi_pmp_ac          )
    CALL allocate_shared_memory_dp_2D(  mesh%nAC, C%nZ, ice%A_flow_ac,           ice%wA_flow_ac          )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%A_flow_mean_ac,      ice%wA_flow_mean_ac     )
    CALL allocate_shared_memory_dp_2D(  mesh%nAC, C%nZ, ice%Cpi_ac,              ice%wCpi_ac             )
    CALL allocate_shared_memory_dp_2D(  mesh%nAC, C%nZ, ice%Ki_ac,               ice%wKi_ac              )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%dHs_dp_ac,           ice%wdHs_dp_ac          )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%dHs_do_ac,           ice%wdHs_do_ac          )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%dHs_dx_ac,           ice%wdHs_dx_ac          )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%dHs_dy_ac,           ice%wdHs_dy_ac          )
    CALL allocate_shared_memory_dp_2D(  mesh%nAC, C%nZ, ice%D_SIA_3D_ac,         ice%wD_SIA_3D_ac        )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%D_SIA_ac,            ice%wD_SIA_ac           )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%Up_ac,               ice%wUp_ac              )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%Uo_ac,               ice%wUo_ac              ) 
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%Ux_ac,               ice%wUx_ac              )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,       ice%Uy_ac,               ice%wUy_ac              ) 
    CALL allocate_shared_memory_dp_1D(  mesh%nV,        ice%U_SIA,               ice%wU_SIA              )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,        ice%V_SIA,               ice%wV_SIA              )
    
    ! Ice dynamics - SSA
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%U_SSA,               ice%wU_SSA              )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%V_SSA,               ice%wV_SSA              )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%tau_yield,           ice%wtau_yield          )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dUdx,                ice%wdUdx               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dUdy,                ice%wdUdy               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dVdx,                ice%wdVdx               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dVdy,                ice%wdVdy               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%nu,                  ice%wnu                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%beta_base,           ice%wbeta_base          )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%rhs_x,               ice%wrhs_x              )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%rhs_y,               ice%wrhs_y              )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%eu_i,                ice%weu_i               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%ev_i,                ice%wev_i               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Rx,                  ice%wRx                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Ry,                  ice%wRy                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,      ice%Ux_SSA_ac,           ice%wUx_SSA_ac          )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,      ice%Uy_SSA_ac,           ice%wUy_SSA_ac          )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,      ice%Up_SSA_ac,           ice%wUp_SSA_ac          )
    CALL allocate_shared_memory_dp_1D(  mesh%nAC,      ice%Uo_SSA_ac,           ice%wUo_SSA_ac          )
    
    ! Ice dynamics - ice thickness calculation
    CALL allocate_shared_memory_dp_2D(  mesh%nV, mesh%nconmax,      ice%dVi_in,              ice%wdVi_in             )
    CALL allocate_shared_memory_dp_2D(  mesh%nV, mesh%nconmax,      ice%dVi_out,             ice%wdVi_out            ) 
    
    ! Thermodynamics
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%frictional_heating,  ice%wfrictional_heating )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Fr,                  ice%wFr                 )
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%U_3D,                ice%wU_3D               )
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%V_3D,                ice%wV_3D               )
    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%W_3D,                ice%wW_3D               )
    
    ! Mesh adaptation data
    CALL allocate_shared_memory_dp_1D(  mesh%nV,        ice%surf_curv,          ice%wsurf_curv          )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,        ice%log_velocity,       ice%wlog_velocity       )
    
  END SUBROUTINE allocate_ice_model
  SUBROUTINE initialise_ice_temperature(ice, mesh, init, filetype_init, climate)
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
    CHARACTER(LEN=4),                    INTENT(IN)    :: filetype_init
    TYPE(type_climate_model),            INTENT(IN)    :: climate
    
    ! Local variables
    INTEGER                                            :: vi, k
    REAL(dp)                                           :: T_surf_annual
    
    IF (filetype_init == 'cart') THEN
     
      ! Calculate Ti_pmp
      CALL ice_physical_properties( ice, mesh)
      
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
      
    ELSEIF (filetype_init == 'mesh') THEN
    
      ice%Ti(mesh%v1:mesh%v2,:) = init%Ti(mesh%v1:mesh%v2,:)
      CALL sync
    
    END IF    
    
    ! EISMINT experiments
    IF (C%do_eismint_experiment .AND. C%use_thermodynamics) THEN
      IF     (C%choice_eismint_experiment == 'A' .OR. &
              C%choice_eismint_experiment == 'N' .OR. &
              C%choice_eismint_experiment == 'O') THEN
        ice%Ti(mesh%v1:mesh%v2,:) = 270._dp
        CALL sync
      ELSEIF (C%choice_eismint_experiment == 'P' .OR. &
              C%choice_eismint_experiment == 'Q' .OR. &
              C%choice_eismint_experiment == 'R') THEN
        DO vi = mesh%v1, mesh%v2
          ice%Ti(vi,:) = climate%applied%T2m(vi,1)
        END DO
        CALL sync
      END IF
    END IF
  END SUBROUTINE initialise_ice_temperature  
  
  
END MODULE ice_dynamics_module
