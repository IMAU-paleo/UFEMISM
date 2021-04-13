MODULE ice_dynamics_module
  ! The actual ice dynamics model, based on the mixed SIA/SSA approach from ANICE.

  USE mpi
  USE configuration_module,          ONLY: dp, C
  USE parameters_module,             ONLY: seawater_density, ice_density, T0
  USE zeta_module,                   ONLY: vertical_integrate, vertical_average
  USE parallel_module,               ONLY: par, sync, &
                                           allocate_shared_int_0D, allocate_shared_dp_0D, &
                                           allocate_shared_int_1D, allocate_shared_dp_1D, &
                                           allocate_shared_int_2D, allocate_shared_dp_2D, &
                                           allocate_shared_int_3D, allocate_shared_dp_3D, &
                                           deallocate_shared                                               
  USE data_types_module,             ONLY: type_mesh, type_ice_model, type_init_data_fields, type_SMB_model, type_BMB_model, &
                                           type_remapping, type_PD_data_fields
  USE netcdf_module,                 ONLY: debug, write_to_debug_file
  USE mesh_derivatives_module,       ONLY: get_mesh_derivatives, get_mesh_curvatures_vertex, apply_Neumann_boundary, &
                                           apply_Neumann_boundary_3D, get_mesh_derivatives_vertex_3D
  USE mesh_mapping_module,           ONLY: remap_field_dp, remap_field_dp_3D, remap_field_dp_monthly, &
                                           reallocate_field_dp, reallocate_field_dp_3D, reallocate_field_int
  USE mesh_ArakawaC_module,          ONLY: map_Aa_to_Ac, map_Ac_to_Aa, rotate_xy_to_po, get_mesh_derivatives_AaAc, &
                                           get_mesh_curvatures_vertex_AaAc, apply_Neumann_boundary_AaAc
  USE general_ice_model_data_module, ONLY: is_floating, update_general_ice_model_data
  USE mesh_help_functions_module,    ONLY: write_mesh_to_text_file

  IMPLICIT NONE
  
CONTAINS  

  ! Update ice thickness from ice dynamic changes and SMB
  SUBROUTINE calculate_ice_thickness_change( mesh, ice, SMB, BMB, dt, mask_noice)
    ! Use the total ice velocities to update the ice thickness
    
    USE parameters_module, ONLY: ice_density
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp),                            INTENT(IN)    :: dt
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_noice
    
    ! Local variables
    INTEGER                                            :: cerr, ierr
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
    
    DO aci = mesh%ac1, mesh%ac2
    
      ! The two Aa vertices connected by the Ac vertex
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      
      ! Find their own respective connectivity, for storing the ice flux
      ci = 0
      cj = 0
      DO cii = 1, mesh%nC( vi)
        IF (mesh%C( vi,cii)==vj) THEN
          ci = cii
          EXIT
        END IF
      END DO
      DO cji = 1, mesh%nC( vj)
        IF (mesh%C( vj,cji)==vi) THEN
          cj = cji
          EXIT
        END IF
      END DO
            
      ! Calculate ice volume per year moving along connection from vi to vj as the product of:
      ! - width          (m   - determined by distance between adjacent triangle circumcenters)
      ! - ice thickness  (m   - at flux origin (upwind scheme, because this is really advection)
      ! - ice velocity   (m/y - calculated at midpoint, using surface slope along connection)
      ! - time step      (y)
      Upar = ice%Up_SIA_Ac( aci) + ice%Up_SSA_Ac( aci)
      IF (Upar > 0._dp) THEN
        dVi = ice%Hi( vi) * Upar * mesh%Cw( vi,ci) * dt ! m3
      ELSE
        dVi = ice%Hi( vj) * Upar * mesh%Cw( vi,ci) * dt ! m3
      END IF
      
      ! Keep track of ice fluxes across individual connections, to correct for
      ! negative ice thicknesses if necessary         
      ice%dVi_in( vi, ci) = -dVi ! m3
      ice%dVi_in( vj, cj) =  dVi ! m3
           
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync
           
    
    ! Correct outfluxes for possible resulting negative ice thicknesses
    ! =================================================================
    
    Vi_SMB( mesh%v1:mesh%v2) = (SMB%SMB_year( mesh%v1:mesh%v2) + BMB%BMB( mesh%v1:mesh%v2))  * mesh%A( mesh%v1:mesh%v2) * dt! * ice_density / 1000._dp     ! m3 ice equivalent
    CALL sync
    
    DO vi = mesh%v1, mesh%v2
    
      ! Check how much ice is available for melting or removing (in m^3)
      Vi_available = mesh%A( vi) * ice%Hi( vi)
      
      dVi = SUM( ice%dVi_in( vi,:))
      
      Vi_in  = 0._dp
      Vi_out = 0._dp
      DO ci = 1, mesh%nC( vi)
        IF (ice%dVi_in( vi,ci) > 0._dp) THEN
          Vi_in  = Vi_in  + ice%dVi_in( vi,ci)
        ELSE
          Vi_out = Vi_out - ice%dVi_in( vi,ci)
        END IF
      END DO
      
      rescale_factor = 1._dp      
      
      ! If all the ice already present melts away, there can be no outflux.
      IF (-Vi_SMB( vi) >= Vi_available) THEN      
        ! All ice in this vertex melts, nothing remains to move around. Rescale outfluxes to zero.
        Vi_SMB( vi) = -Vi_available
        rescale_factor = 0._dp
      END IF
            
      ! If the total outflux exceeds the available ice plus SMB plus total influx, rescale outfluxes
      IF (Vi_out > Vi_available + Vi_SMB( vi)) THEN      
        ! Total outflux plus melt exceeds available ice volume. Rescale outfluxes to correct for this.
        rescale_factor = (Vi_available + Vi_SMB( vi)) / Vi_out        
      END IF
      
      ! Rescale ice outfluxes out of vi and into vi's neighbours
      IF (rescale_factor < 1._dp) THEN
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
        
          IF (ice%dVi_in( vi,ci) < 0._dp) THEN
            ice%dVi_in( vi,ci) = ice%dVi_in( vi,ci) * rescale_factor
              
            DO cji = 1, mesh%nC( vj)
              IF (mesh%C( vj,cji) == vi) THEN
                ice%dVi_in( vj,cji) = -ice%dVi_in( vi,ci)
                EXIT
              END IF
            END DO
          
          END IF ! IF (ice%dVi_in( vi,ci) < 0._dp) THEN
        END DO ! DO ci = 1, mesh%nC( vi)
      END IF ! IF (rescale_factor < 1._dp) THEN
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! Calculate change in ice thickness over time at every vertex
    ! ===========================================================
    
    ice%dHi_dt( mesh%v1:mesh%v2) = 0._dp ! m/y    
    DO vi = mesh%v1, mesh%v2
      dVi  = SUM( ice%dVi_in( vi,:))
      ice%dHi_dt( vi) = (dVi + Vi_SMB( vi)) / (mesh%A( vi) * dt)  
    END DO    
    CALL sync
    
    ! Manually set to zero in first timestep.
    IF (dt==0._dp) ice%dHi_dt( mesh%v1:mesh%v2) = 0._dp
    CALL sync
    
    ! Add dynamic ice thickness change and SMB to update the ice thickness
    ! ====================================================================
    
    ice%Hi_prev( mesh%v1:mesh%v2) = ice%Hi( mesh%v1:mesh%v2)
    ice%Hi(      mesh%v1:mesh%v2) = ice%Hi( mesh%v1:mesh%v2) + (ice%dHi_dt( mesh%v1:mesh%v2) * dt)

    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN

          ! Apply boundary conditions: set ice thickness to zero at the domain boundary
          DO vi = mesh%v1, mesh%v2
            IF (mesh%edge_index( vi) > 0) ice%Hi(vi) = 0._dp
          END DO
          CALL sync
        
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
        ! No exception here, as we already exited the routine at the top.
!      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
!        
!        ! Create a nice circular ice shelf
!        DO vi = mesh%v1, mesh%v2
!          IF (SQRT(mesh%V(vi,1)**2+mesh%V(vi,2)**2) > mesh%xmax * 0.95_dp) ice%Hi( vi) = 0._dp
!        END DO
!        CALL sync
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calculate_ice_thickness_change!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE

      ! Apply boundary conditions: set ice thickness to zero at the domain boundary
      DO vi = mesh%v1, mesh%v2
        IF (mesh%edge_index( vi) > 0) ice%Hi(vi) = 0._dp
      END DO
      CALL sync
      
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Remove ice in areas where no ice is allowed (i.e. Greenland in NAM and EAS, and Ellesmere Island in GRL)
    DO vi = mesh%v1, mesh%v2
      IF (mask_noice( vi) == 1) ice%Hi( vi) = 0._dp
    END DO
    CALL sync
    
  END SUBROUTINE calculate_ice_thickness_change
  
  ! SIA solver
  SUBROUTINE solve_SIA( mesh, ice)
    ! Calculate ice velocities using the SIA
    
    USE parameters_module, ONLY: n_flow, ice_density, grav
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    INTEGER                                            :: aci, k
    REAL(dp)                                           :: D_0
    REAL(dp), DIMENSION(C%nZ)                          :: D_deformation
    REAL(dp)                                           :: D_uv_2D 
    REAL(dp), DIMENSION(C%nZ)                          :: D_SIA_3D_prof    
    REAL(dp), PARAMETER                                :: D_uv_3D_cutoff = -1E5_dp  
        
    ice%D_SIA_3D_Ac(mesh%ac1:mesh%ac2,:) = 0._dp
    ice%D_SIA_Ac(   mesh%ac1:mesh%ac2  ) = 0._dp
    ice%D_SIA(      mesh%v1:mesh%v2    ) = 0._dp
    ice%U_SIA(      mesh%v1:mesh%v2    ) = 0._dp
    ice%V_SIA(      mesh%v1:mesh%v2    ) = 0._dp
    ice%Ux_SIA_Ac(  mesh%ac1:mesh%ac2  ) = 0._dp
    ice%Uy_SIA_Ac(  mesh%ac1:mesh%ac2  ) = 0._dp
    ice%Up_SIA_Ac(  mesh%ac1:mesh%ac2  ) = 0._dp
    ice%Uo_SIA_Ac(  mesh%ac1:mesh%ac2  ) = 0._dp
    CALL sync
        
    ! From the ANICE subroutine "calculate_D_uv_3D"
    
    DO aci = mesh%ac1, mesh%ac2
    
      IF (ice%mask_sheet_Ac( aci) == 1) THEN     
        
       D_0      = (ice_density * grav * ice%Hi_Ac( aci))**n_flow * ((ice%dHs_dp_Ac( aci)**2 + ice%dHs_do_Ac( aci)**2))**((n_flow - 1._dp) / 2._dp)       
       D_deformation = vertical_integrate( C%m_enh_sia * ice%A_flow_Ac( aci,:) * C%zeta**n_flow)       
       D_deformation = 2._dp * ice%Hi_Ac( aci) * D_deformation
       
       ice%D_SIA_3D_Ac( aci,:) = D_0 * D_deformation

      END IF ! End: Considering only grounded points
      
      ! Check for very large D_uv_3D's, causing very large veloacities RESET --DIRTY
      DO k = 1, C%nZ
        IF (ice%D_SIA_3D_Ac( aci,k) < D_uv_3D_cutoff) THEN
          ice%D_SIA_3D_Ac( aci,k) = D_uv_3D_cutoff          
        END IF        
      END DO

    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync
    
    ! From the ANICE subroutine "sheet_veloacity_2D"                                                    
    DO aci = mesh%ac1, mesh%ac2
      IF (ice%mask_sheet_Ac( aci) == 1) THEN
        D_SIA_3D_prof = ice%D_SIA_3D_Ac( aci,:)
        D_uv_2D = vertical_average(D_SIA_3D_prof)
       
        ice%D_SIA_Ac(  aci) = ice%Hi_Ac( aci) * D_uv_2D     ! See equation (8.5)
        ice%Ux_SIA_Ac( aci) = D_uv_2D * ice%dHs_dx_Ac( aci) ! See equation (7.39)
        ice%Uy_SIA_Ac( aci) = D_uv_2D * ice%dHs_dy_Ac( aci) ! See equation (7.40)
        ice%Up_SIA_Ac( aci) = D_uv_2D * ice%dHs_dp_Ac( aci)
        ice%Uo_SIA_Ac( aci) = D_uv_2D * ice%dHs_do_Ac( aci)        
      END IF      
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync
    
    ! Map data to Aa mesh for writing to output (not actually used anywhere...)    
    CALL map_Ac_to_Aa( mesh, ice%Ux_SIA_Ac, ice%U_SIA)
    CALL map_Ac_to_Aa( mesh, ice%Uy_SIA_Ac, ice%V_SIA)
    CALL map_Ac_to_Aa( mesh, ice%D_SIA_Ac,  ice%D_SIA)
      
  END SUBROUTINE solve_SIA
  
  ! 3D SIA solver (used only for thermodynamics)
  SUBROUTINE solve_SIA_3D( mesh, ice)
    
    USE parameters_module, ONLY: n_flow, ice_density, grav
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
      
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
    
      ice%U_3D( vi,:) = ice%U_SSA( vi)
      ice%V_3D( vi,:) = ice%V_SSA( vi)
    
      IF (ice%mask_shelf( vi) == 1) CYCLE
      IF (ice%Hi( vi) == 0._dp) CYCLE
        
      D_0           = (ice_density * grav * ice%Hi( vi))**n_flow * ((ice%dHs_dx( vi)**2 + ice%dHs_dy( vi)**2))**((n_flow - 1._dp) / 2._dp)
      D_deformation = vertical_integrate( C%m_enh_sia * ice%A_flow( vi,:) * C%zeta**n_flow)
      D_deformation = 2._dp * ice%Hi( vi) * D_deformation
      D_SIA_3D      = MAX(D_0 * D_deformation, D_uv_3D_cutoff)
       
      ! To be completely consistent, only add the SSA basal velocities when calculated:
      ice%U_3D( vi,:)    = D_SIA_3D * ice%dHs_dx( vi) + ice%U_SSA( vi)
      ice%V_3D( vi,:)    = D_SIA_3D * ice%dHs_dy( vi) + ice%V_SSA( vi)
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
        
    CALL apply_Neumann_boundary_3D( mesh, ice%U_3D, C%nZ)
    CALL apply_Neumann_boundary_3D( mesh, ice%V_3D, C%nZ)

    DO vi = mesh%v1, mesh%v2
    
      IF (mesh%edge_index( vi) > 0) CYCLE
      IF (ice%mask_sheet( vi) == 0) CYCLE
      
      dHb_dx = ice%dHs_dx( vi) - ice%dHi_dx( vi)
      dHb_dy = ice%dHs_dy( vi) - ice%dHi_dy( vi)   
      
      ice%W_3D( vi,C%nZ) = ice%dHb_dt( vi) + ice%U_3D( vi,C%nZ) * dHb_dx + ice%V_3D( vi,C%nZ) * dHb_dy 
                           
      ! The integrant is calculated half way the layer of integration at k+1/2. This integrant is multiplied with the layer thickness and added to the integral
      ! of all layers below, giving the integral up to and including this layer:
      DO k = C%nZ - 1, 1, -1
      
        CALL get_mesh_derivatives_vertex_3D( mesh, ice%U_3D, dUdx_k,   dUdy_k,   vi, k)
        CALL get_mesh_derivatives_vertex_3D( mesh, ice%U_3D, dUdx_kp1, dUdy_kp1, vi, k+1)
        CALL get_mesh_derivatives_vertex_3D( mesh, ice%V_3D, dVdx_k,   dVdy_k,   vi, k)
        CALL get_mesh_derivatives_vertex_3D( mesh, ice%V_3D, dVdx_kp1, dVdy_kp1, vi, k+1)
        
        w1 = (dUdx_k + dUdx_kp1) / 2._dp
        w2 = (dVdy_k + dVdy_kp1) / 2._dp   
        
        w3              = ((ice%dHs_dx( vi) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * ice%dHi_dx( vi)) / MAX(0.1_dp, ice%Hi( vi))) *   &
                          ((ice%U_3D( vi,k+1) - ice%U_3D( vi,k)) / (C%zeta(k+1) - C%zeta(k)))
        w4              = ((ice%dHs_dy( vi) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * ice%dHi_dy( vi)) / MAX(0.1_dp, ice%Hi( vi))) *   &
                          ((ice%V_3D( vi,k+1) - ice%V_3D( vi,k)) / (C%zeta(k+1) - C%zeta(k)))

        ice%W_3D(vi,k) = ice%W_3D(vi,k+1) - ice%Hi(vi) * (w1 + w2 + w3 + w4) * (C%zeta(k+1) - C%zeta(k))
        
      END DO ! DO k = C%nZ - 1, 1, -1

    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    CALL apply_Neumann_boundary_3D( mesh, ice%W_3D, C%nZ)

  END SUBROUTINE solve_SIA_3D
  
  ! SSA solver
  SUBROUTINE solve_SSA( mesh, ice)
    ! Calculate ice velocities using the SSA
    
    USE parameters_module, ONLY : n_flow, ice_density, grav, sec_per_year
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: cerr,ierr
    LOGICAL                                            :: set_SSA_velocities_to_zero
    LOGICAL                                            :: has_converged
    INTEGER                                            :: viscosity_iteration_i
    REAL(dp)                                           :: sum_DN_sq, sum_N_sq, RN
    LOGICAL                                            :: did_reset_before, did_reset_now
    
    did_reset_before = .FALSE.
    did_reset_now    = .FALSE.
        
    ! Check if we really need to, or are even able to, solve the SSA.
    set_SSA_velocities_to_zero = .FALSE.
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler') THEN
        ! The SSA is not solved in these experiments, set velocities to zero.
        set_SSA_velocities_to_zero = .TRUE.
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
              C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        ! The SSA is solved in these experiments.
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in solve_SSA!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! If there's no grounded ice anywhere, don't solve the SSA.
    IF (SUM( ice%mask_sheet) == 0) set_SSA_velocities_to_zero = .TRUE.
                   
    IF (set_SSA_velocities_to_zero) THEN 
      ice%U_SSA_AaAc( mesh%a1:mesh%a2  ) = 0._dp
      ice%V_SSA_AaAc( mesh%a1:mesh%a2  ) = 0._dp
      ice%U_SSA(      mesh%v1:mesh%v2  ) = 0._dp
      ice%V_SSA(      mesh%v1:mesh%v2  ) = 0._dp
      ice%Ux_SSA_Ac(  mesh%ac1:mesh%ac2) = 0._dp
      ice%Uy_SSA_Ac(  mesh%ac1:mesh%ac2) = 0._dp
      ice%Up_SSA_Ac(  mesh%ac1:mesh%ac2) = 0._dp
      ice%Uo_SSA_Ac(  mesh%ac1:mesh%ac2) = 0._dp
      CALL sync
      RETURN
    END IF
    
    ! Calculate the basal yield stress
    CALL basal_yield_stress( mesh, ice)
    
    ! Calculate the semi-analytical grounding-line flux solution
    IF (C%use_analytical_GL_flux) CALL calculate_GL_flux( mesh, ice)
    
    ! Solve the SSA
    ! Same approach as an ANICE: an outer loop where ice viscosity is updated,
    ! and an inner loop that uses Successive Over-Relaxation (SOR) to solve the two coupled PDE's for the given viscosity.
    
    ! Gather the Aa and Ac fields into a single variables
    ice%Hi_AaAc(                   mesh%v1 :        mesh%v2 ) = ice%Hi(              mesh%v1 :mesh%v2 )
    ice%Hb_AaAc(                   mesh%v1 :        mesh%v2 ) = ice%Hb(              mesh%v1 :mesh%v2 )
    ice%SL_AaAc(                   mesh%v1 :        mesh%v2 ) = ice%SL(              mesh%v1 :mesh%v2 )
    ice%dHs_dx_shelf_AaAc(         mesh%v1 :        mesh%v2 ) = ice%dHs_dx_shelf(    mesh%v1 :mesh%v2 )
    ice%dHs_dy_shelf_AaAc(         mesh%v1 :        mesh%v2 ) = ice%dHs_dy_shelf(    mesh%v1 :mesh%v2 )
    ice%A_flow_mean_AaAc(          mesh%v1 :        mesh%v2 ) = ice%A_flow_mean(     mesh%v1 :mesh%v2 )
    ice%U_SSA_AaAc(                mesh%v1 :        mesh%v2 ) = ice%U_SSA(           mesh%v1 :mesh%v2 )
    ice%V_SSA_AaAc(                mesh%v1 :        mesh%v2 ) = ice%V_SSA(           mesh%v1 :mesh%v2 )
    
    ice%Hi_AaAc(           mesh%nV+mesh%ac1:mesh%nV+mesh%ac2) = ice%Hi_Ac(           mesh%ac1:mesh%ac2)
    ice%Hb_AaAc(           mesh%nV+mesh%ac1:mesh%nV+mesh%ac2) = ice%Hb_Ac(           mesh%ac1:mesh%ac2)
    ice%SL_AaAc(           mesh%nV+mesh%ac1:mesh%nV+mesh%ac2) = ice%SL_Ac(           mesh%ac1:mesh%ac2)
    ice%dHs_dx_shelf_AaAc( mesh%nV+mesh%ac1:mesh%nV+mesh%ac2) = ice%dHs_dx_shelf_Ac( mesh%ac1:mesh%ac2)
    ice%dHs_dy_shelf_AaAc( mesh%nV+mesh%ac1:mesh%nV+mesh%ac2) = ice%dHs_dy_shelf_Ac( mesh%ac1:mesh%ac2)
    ice%A_flow_mean_AaAc(  mesh%nV+mesh%ac1:mesh%nV+mesh%ac2) = ice%A_flow_mean_Ac(  mesh%ac1:mesh%ac2)
    ice%U_SSA_AaAc(        mesh%nV+mesh%ac1:mesh%nV+mesh%ac2) = ice%Ux_SSA_Ac(       mesh%ac1:mesh%ac2)
    ice%V_SSA_AaAc(        mesh%nV+mesh%ac1:mesh%nV+mesh%ac2) = ice%Uy_SSA_Ac(       mesh%ac1:mesh%ac2)
    CALL sync
        
    ! The viscosity iteration
    has_converged = .FALSE.
    viscosity_iteration_i  = 0
    !DO outer_loop_i = 1, outer_loop_n
    viscosity_iteration: DO WHILE ((.NOT. has_converged) .AND. (viscosity_iteration_i < C%SSA_max_outer_loops))
      viscosity_iteration_i = viscosity_iteration_i + 1
      
      ! Store the previous solution, so we can check for convergence later
      ice%N_AaAc_prev( mesh%a1:mesh%a2)  = ice%N_AaAc( mesh%a1:mesh%a2)
      CALL sync
      
      ! Update the effective viscosity eta and the product term N = eta * H
      CALL SSA_effective_viscosity( mesh, ice)
      
      ! Check if a stable solution has been reached
      sum_DN_sq = SUM( (ice%N_AaAc(      mesh%a1:mesh%a2) - &
                        ice%N_AaAc_prev( mesh%a1:mesh%a2))**2)
      sum_N_sq  = SUM(  ice%N_AaAc(      mesh%a1:mesh%a2)**2)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, sum_DN_sq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, sum_N_sq,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      RN = SQRT( sum_DN_sq / sum_N_sq)
      
      !IF (par%master) WRITE(0,'(A,I5,A,E10.4)') '   SSA - viscosity iteration ', viscosity_iteration_i, ', RN = ', RN
      
      IF (RN < C%SSA_RN_tol) THEN
        has_converged = .TRUE.
        EXIT viscosity_iteration
      END IF
      
      ! Calculate the sliding term with the velocity field from the previous viscosity iteration
      CALL SSA_sliding_term( mesh, ice)
    
      ! Solve the linearised SSA using SOR
      CALL solve_SSA_linearised( mesh, ice, did_reset_now)
      
      ! Check if velocities were reset due to instability. If this happens twice, crash.
      IF (did_reset_now) THEN
        IF (.NOT. did_reset_before) THEN
          did_reset_before = .TRUE.
        ELSE
          IF (par%master) WRITE(0,*) 'solve_SSA - ERROR: SSA remains unstable after resetting velocities to zero!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END IF
      
      ! DENK DROM
     ! EXIT
      
    END DO viscosity_iteration ! DO WHILE ((.NOT. has_converged) .AND. (viscosity_iteration_i < C%SSA_max_outer_loops))
    
    ! Map data back to separate Aa and Ac meshes
    ice%U_SSA(     mesh%v1 :mesh%v2 ) = ice%U_SSA_AaAc(         mesh%v1 :        mesh%v2 )
    ice%V_SSA(     mesh%v1 :mesh%v2 ) = ice%V_SSA_AaAc(         mesh%v1 :        mesh%v2 )
    ice%Ux_SSA_Ac( mesh%ac1:mesh%ac2) = ice%U_SSA_AaAc( mesh%nV+mesh%ac1:mesh%nV+mesh%ac2)
    ice%Uy_SSA_Ac( mesh%ac1:mesh%ac2) = ice%V_SSA_AaAc( mesh%nV+mesh%ac1:mesh%nV+mesh%ac2)
    CALL sync
    
    ! Get parallel and orthogonal components of the Ac velocity field for ice thickness update
    CALL rotate_xy_to_po( mesh, ice%Ux_SSA_Ac, ice%Uy_SSA_Ac, ice%Up_SSA_Ac, ice%Uo_SSA_Ac)
    
  END SUBROUTINE solve_SSA
  SUBROUTINE solve_SSA_linearised( mesh, ice, did_reset)
    ! Calculate ice velocities using the SSA
    
    USE parameters_module, ONLY : n_flow, ice_density, grav, sec_per_year
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    LOGICAL,                             INTENT(OUT)   :: did_reset
    
    ! Local variables:
    INTEGER                                            :: cerr,ierr
    INTEGER                                            :: ai, ci, ac
    LOGICAL                                            :: is_edge
    REAL(dp)                                           :: Uxxi, Uxyi, Uyyi, Vxxi, Vxyi, Vyyi
    REAL(dp)                                           :: sumUc, sumVc
    INTEGER                                            :: fci, fcvi1, fcvi2, fcvi
    LOGICAL                                            :: has_converged
    INTEGER                                            :: inner_loop_i
    REAL(dp)                                           :: max_residual_UV
    
    cerr = 0
    ierr = 0
      
    ! Calculate the right-hand sides of the PDE's
    DO ai = mesh%a1, mesh%a2
      ice%RHSx_AaAc( ai) = ice_density * grav * ice%dHs_dx_shelf_AaAc( ai) / ice%eta_AaAc( ai)
      ice%RHSy_AaAc( ai) = ice_density * grav * ice%dHs_dy_shelf_AaAc( ai) / ice%eta_AaAc( ai)
    END DO ! DO ai = mesh%a1, mesh%a2
    
    ! Calculate the centre coefficients
    DO ai = mesh%a1, mesh%a2
      IF (.NOT. is_floating( ice%Hi_AaAc( ai), ice%Hb_AaAc( ai), ice%SL_AaAc( ai))) THEN
        ice%eu_i_AaAc( ai) = (4._dp * mesh%Nxx_AaAc( ai,mesh%nCAaAc( ai)+1) + mesh%Nyy_AaAc( ai,mesh%nCAaAc( ai)+1)) - ice%S_AaAc( ai) / (MAX(0.1_dp,ice%Hi_AaAc( ai)) * ice%eta_AaAc( ai))
        ice%ev_i_AaAc( ai) = (4._dp * mesh%Nyy_AaAc( ai,mesh%nCAaAc( ai)+1) + mesh%Nxx_AaAc( ai,mesh%nCAaAc( ai)+1)) - ice%S_AaAc( ai) / (MAX(0.1_dp,ice%Hi_AaAc( ai)) * ice%eta_AaAc( ai))
      ELSE
        ice%eu_i_AaAc( ai) = (4._dp * mesh%Nxx_AaAc( ai,mesh%nCAaAc( ai)+1) + mesh%Nyy_AaAc( ai,mesh%nCAaAc( ai)+1))
        ice%ev_i_AaAc( ai) = (4._dp * mesh%Nyy_AaAc( ai,mesh%nCAaAc( ai)+1) + mesh%Nxx_AaAc( ai,mesh%nCAaAc( ai)+1))
      END IF
    END DO ! DO ai = mesh%a1, mesh%a2
    CALL sync
      
    ! Use SOR to solve U and V for the given right-hand side of the PDE's
    has_converged = .FALSE.
    inner_loop_i  = 0
    DO WHILE ((.NOT. has_converged) .AND. (inner_loop_i < C%SSA_max_inner_loops))
      inner_loop_i = inner_loop_i + 1
      
      max_residual_UV = 0._dp
      
      ! Go through one iteration of the SOR scheme, using a five-colour partitioning
      DO fci = 1, 5
      
        fcvi1 = mesh%colour_v1( fci)
        fcvi2 = mesh%colour_v2( fci)
      
        DO fcvi = fcvi1, fcvi2
        
          ai = mesh%colour_vi( fcvi, fci)
        
          ! Don't update edge vertices, those will be treated by the boundary conditions
          is_edge = .FALSE.
          IF (ai <= mesh%nV) THEN
            IF (mesh%edge_index(    ai        )>0) is_edge = .TRUE.
          ELSE
            IF (mesh%edge_index_Ac( ai-mesh%nV)>0) is_edge = .TRUE.
          END IF
          IF (is_edge) CYCLE
          
          ! If we're using the analytical GL flux solution, values at the grounding line (on the Ac mesh)
          ! will be prescribed by that solution, so we skip those vertices
          IF (C%use_analytical_GL_flux) THEN
            IF (ai > mesh%nV) THEN
              IF (ice%mask_gl_Ac( ai-mesh%nV) == 1) CYCLE
            END IF
          END IF
      
          ! Calculate Uxy, Vxy with partially updated values pf U and V according to Gauss-Seidler
          CALL get_mesh_curvatures_vertex_AaAc( mesh, ice%U_SSA_AaAc, Uxxi, Uxyi, Uyyi, ai)
          CALL get_mesh_curvatures_vertex_AaAc( mesh, ice%V_SSA_AaAc, Vxxi, Vxyi, Vyyi, ai)
          
          ! The sum terms in the equation
          sumUc = 0._dp
          sumVc = 0._dp
          DO ci = 1, mesh%nCAaAc( ai)
            ac = mesh%CAaAc( ai,ci)
            sumUc = sumUc + ice%U_SSA_AaAc( ac) * (4._dp * mesh%Nxx_AaAc( ai,ci) + mesh%Nyy_AaAc( ai,ci))
            sumVc = sumVc + ice%V_SSA_AaAc( ac) * (4._dp * mesh%Nyy_AaAc( ai,ci) + mesh%Nxx_AaAc( ai,ci))
          END DO
          
          ! Calculate the left-hand side of the equation, using partially updated values
          ice%LHSx_AaAc( ai) = sumUc + (3._dp * Vxyi) + (ice%eu_i_AaAc( ai) * ice%U_SSA_AaAc( ai))
          ice%LHSy_AaAc( ai) = sumVc + (3._dp * Uxyi) + (ice%ev_i_AaAc( ai) * ice%V_SSA_AaAc( ai))
          
          ! Calculate the residuals
          ice%resU_AaAc( ai) = (ice%LHSx_AaAc( ai) - ice%RHSx_AaAc( ai)) / ice%eu_i_AaAc( ai)
          ice%resV_AaAc( ai) = (ice%LHSy_AaAc( ai) - ice%RHSy_AaAc( ai)) / ice%ev_i_AaAc( ai)
        
          max_residual_UV = MAX( max_residual_UV, ABS(ice%resU_AaAc( ai)))
          max_residual_UV = MAX( max_residual_UV, ABS(ice%resV_AaAc( ai)))
          
          ! Update velocities
          ice%U_SSA_AaAc( ai) = ice%U_SSA_AaAc( ai) - C%SSA_SOR_omega * ice%resU_AaAc( ai)
          ice%V_SSA_AaAc( ai) = ice%V_SSA_AaAc( ai) - C%SSA_SOR_omega * ice%resV_AaAc( ai)
        
        END DO ! DO fcvi = fcvi1, fcvi2
        CALL sync
      
      END DO ! DO fci = 1, 5
      
      !IF (par%master) WRITE(0,'(A,I5,A,E10.4)') '   SSA -   SOR iteration ', inner_loop_i, ', max_residual_UV = ', max_residual_UV
      
      ! Apply Neumann boundary conditions
      CALL apply_Neumann_boundary_AaAc( mesh, ice%U_SSA_AaAc)
      CALL apply_Neumann_boundary_AaAc( mesh, ice%V_SSA_AaAc)
      
      ! Check if we've reached a stable solution
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, max_residual_UV, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      
      !IF (par%master) WRITE(0,*) ' SSA -  inner loop ', inner_loop_i, ': largest residual = ', max_residual_UV
      IF (max_residual_UV < C%SSA_max_residual_UV) THEN
        did_reset     = .FALSE.
        has_converged = .TRUE.
      ELSEIF (max_residual_UV > 1E6_dp) THEN
        ! Instability detected. Reset velocities to zero. If that doesn't work, crash (do this from the outer loop)
        ice%U_SSA_AaAc( mesh%a1:mesh%a2) = 0._dp
        ice%V_SSA_AaAc( mesh%a1:mesh%a2) = 0._dp
        did_reset     = .TRUE.
        has_converged = .TRUE.
      ELSEIF (inner_loop_i == C%SSA_max_inner_loops) THEN
        WRITE(0,*) ' WARNING - SSA SOR solver doesnt converge!'
        !IF (par%master) WRITE(0,*) ' ERROR - SSA SOR solver doesnt converge!'
        !CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    
    END DO
    CALL sync
    
  END SUBROUTINE solve_SSA_linearised
  SUBROUTINE SSA_effective_viscosity( mesh, ice)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradient of N in the SSA
    
    USE parameters_module, ONLY: n_flow
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: ai
    
    REAL(dp), PARAMETER                                :: epsilon_sq_0 = 1E-12_dp
      
    ! Calculate derivatives of U and V
    CALL get_mesh_derivatives_AaAc( mesh, ice%U_SSA_AaAc, ice%dU_SSA_dx_AaAc, ice%dU_SSA_dy_AaAc)
    CALL get_mesh_derivatives_AaAc( mesh, ice%V_SSA_AaAc, ice%dV_SSA_dx_AaAc, ice%dV_SSA_dy_AaAc)
          
    ! The effective viscosity eta on both grids, with a normalisation term (following Bueler & Brown, 2009) to prevent divide-by-zero errors,
    ! and the product term N = eta * H
    DO ai = mesh%a1, mesh%a2
      ice%eta_AaAc( ai) = (C%m_enh_ssa * 0.5_dp * ice%A_flow_mean_AaAc( ai))**(-1._dp / n_flow) * (ice%dU_SSA_dx_AaAc( ai)**2 + ice%dV_SSA_dy_AaAc( ai)**2 + &
        ice%dU_SSA_dx_AaAc( ai) * ice%dV_SSA_dy_AaAc( ai) + 0.25_dp * (ice%dU_SSA_dy_AaAc( ai) + ice%dV_SSA_dx_AaAc( ai))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow))
    !  ice%eta_AaAc( ai) = 0.5_dp * (C%m_enh_ssa * ice%A_flow_mean_AaAc( ai))**(-1._dp / n_flow) * (ice%dU_SSA_dx_AaAc( ai)**2 + ice%dV_SSA_dy_AaAc( ai)**2 + &
    !    ice%dU_SSA_dx_AaAc( ai) * ice%dV_SSA_dy_AaAc( ai) + 0.25_dp * (ice%dU_SSA_dy_AaAc( ai) + ice%dV_SSA_dx_AaAc( ai))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow))
      ice%N_AaAc(  ai) = ice%eta_AaAc( ai) * MAX(0.1_dp, ice%Hi_AaAc( ai))
    END DO ! DO ai = mesh%a1, mesh%a2
    CALL sync
    
  END SUBROUTINE SSA_effective_viscosity
  SUBROUTINE SSA_sliding_term( mesh, ice)
    ! Calculate the sliding term S = tau_c / norm([U,V]) in the SSA
    
    USE parameters_module, ONLY: sec_per_year
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: ai
    REAL(dp), PARAMETER                                :: delta_v = 1E-3_dp                    ! Normalisation parameter to prevent errors when velocity is zero
    REAL(dp), PARAMETER                                :: q_plastic            = 0.30_dp       ! Parameter used for basal stress (inverse of m_flow)
    REAL(dp), PARAMETER                                :: u_threshold          = 100._dp       ! scaling of tau_yield to get the correct unit (function of q_plastic)
      
    ! Calculate the velocity-dependent sliding terms (depending on the chosen sliding law)
    IF     (C%choice_sliding_law == 'Coulomb_regularised') THEN
    
      DO ai = mesh%a1, mesh%a2
        ! The sliding term S, with a normalisation term (following Bueler & Brown, 2009) to prevent divide-by-zero errors.
        ice%S_AaAc( ai) = ice%tau_c_AaAc( ai) * ( (delta_v**2 + ice%U_SSA_AaAc( ai)**2 + ice%V_SSA_AaAc( ai)**2)**(0.5_dp * (q_plastic-1._dp)) ) / (u_threshold**q_plastic)
      END DO ! DO ai = mesh%a1, mesh%a2
      CALL sync
      
!    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
!    
!      DO ai = mesh%a1, mesh%a2
!      
!        is_sheet = .FALSE.
!        IF (ai <= mesh%nV) THEN
!          IF (ice%mask_sheet(    ai        ) == 1) is_sheet = .TRUE.
!        ELSE
!          IF (ice%mask_sheet_Ac( ai-mesh%nV) == 1) is_sheet = .TRUE.
!        END IF
!        
!        IF (is_sheet) THEN
!          ice%sliding_term_x_AaAc( ai) = (C%C_sliding / (sec_per_year ** C%m_sliding)) * ((ABS(ice%U_SSA_AaAc( ai)) + delta_v) ** (C%m_sliding - 1._dp)) / (ice%eta_AaAc( ai) * ice%Hi_AaAc( ai))
!          ice%sliding_term_y_AaAc( ai) = (C%C_sliding / (sec_per_year ** C%m_sliding)) * ((ABS(ice%V_SSA_AaAc( ai)) + delta_v) ** (C%m_sliding - 1._dp)) / (ice%eta_AaAc( ai) * ice%Hi_AaAc( ai))
!        ELSE
!          ice%sliding_term_x_AaAc( ai) = 0._dp
!          ice%sliding_term_y_AaAc( ai) = 0._dp
!        END IF
!        
!      END DO ! DO vi = mesh%v1, mesh%v2
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: sliding law "', TRIM(C%choice_sliding_law), '" not implemented in solve_SSA!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF ! IF (C%choice_sliding_law == 'Coulomb') THEN
    
  END SUBROUTINE SSA_sliding_term
  SUBROUTINE basal_yield_stress( mesh, ice)
    ! Calculate the basal yield stress, which is used to determine the basal stress, used 
    ! for sliding. Using the parameterisations given by Martin et al. (TCD: 2010) 
    ! As used in the PISM-PIK model
    
    USE parameters_module, ONLY : pi, ice_density, grav
    
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! local variables
    REAL(dp)                                            :: lambda_p              ! scaling of the pore water pressure
    REAL(dp)                                            :: pore_water_pressure   ! pore water pressure [Pa]
    
    REAL(dp), PARAMETER                                 :: pf1   = -1000._dp
    REAL(dp), PARAMETER                                 :: pf2   = 0._dp
    REAL(dp), PARAMETER                                 :: p_min = 5._dp
    REAL(dp), PARAMETER                                 :: p_max = 20._dp 
    
    INTEGER                                             :: ai, vi, aci
    
    DO vi = mesh%v1, mesh%v2
    
      ai = vi
    
      ! the pore water pressure is scaled with a bedrock height dependend parameterisation
      ! Equation (13) in Martin et al. (2011)
      lambda_p = MAX(0._dp, MIN(1._dp, (1._dp - (ice%Hb( vi) - ice%SL( vi)) / 1000._dp)))
      
      ! The pore water pressure, equation (12) in Martin et al. (2011)
      pore_water_pressure = 0.96_dp * ice_density * grav * MAX(0.1_dp, ice%Hi( vi)) * lambda_p 
      
      ! The friction angle, used for the yield stress, equation (11) in Martin et al. (2011)
      ice%phi_fric_AaAc( ai) = MAX(p_min, MIN(p_max, (p_min + (p_max - p_min) * (1._dp + (ice%Hb( vi) - pf2) / (pf2 - pf1))) ))
      
      ! The yield stress, equation (10) in Martin et al. (2011)
      ice%tau_c_AaAc( ai) = TAN((pi / 180._dp) * ice%phi_fric_AaAc( ai)) * (ice_density * grav * MAX(0.1_dp, ice%Hi( vi)) - pore_water_pressure)
      
    END DO
    CALL sync
    
    DO aci = mesh%ac1, mesh%ac2
    
      ai = aci + mesh%nV
    
      ! the pore water pressure is scaled with a bedrock height dependend parameterisation
      ! Equation (13) in Martin et al. (2011)
      lambda_p = MAX(0._dp, MIN(1._dp, (1._dp - (ice%Hb_Ac( aci) - ice%SL_Ac( aci)) / 1000._dp)))
      
      ! The pore water pressure, equation (12) in Martin et al. (2011)
      pore_water_pressure = 0.96_dp * ice_density * grav * MAX(0.1_dp, ice%Hi_Ac( aci)) * lambda_p 
      
      ! The friction angle, used for the yield stress, equation (11) in Martin et al. (2011)
      ice%phi_fric_AaAc( ai) = MAX(p_min, MIN(p_max, (p_min + (p_max - p_min) * (1._dp + (ice%Hb_Ac( aci) - pf2) / (pf2 - pf1))) ))
      
      ! The yield stress, equation (10) in Martin et al. (2011)
      ice%tau_c_AaAc( ai) = TAN((pi / 180._dp) * ice%phi_fric_AaAc( ai)) * (ice_density * grav * MAX(0.1_dp, ice%Hi_Ac( aci)) - pore_water_pressure)
      
    END DO
    CALL sync
   
  END SUBROUTINE basal_yield_stress
  
  ! Some routines for applying the semi-analytical GL flux solution
  ! as a boundary condition to the SSA solver
  SUBROUTINE calculate_GL_flux( mesh, ice)
    ! Use the semi-analytical solutions by Schoof (2007, in case of a Weertman-type sliding law),
    ! or by Tsai et al. (2012, in case of a Coulomb-type sliding law) to determine
    ! the grounding-line flux.
    
    USE parameters_module, ONLY : pi, n_flow, ice_density, seawater_density, grav, sec_per_year
      
    IMPLICIT NONE
    
    ! In/output variables: 
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: aci, vi, vj
    REAL(dp)                                           :: TAFi, TAFj, lambda_GL, Hi_GL, phi_fric_GL, A_flow_GL
    REAL(dp)                                           :: exponent_Schoof, C_sliding, factor_Schoof, factor_Tsai
    REAL(dp), PARAMETER                                :: Q0 = 0.61_dp
    REAL(dp)                                           :: Fx, Fy, F
    REAL(dp)                                           :: Dx, Dy, D
    
    ! Initialise
    ice%Qabs_GL_Ac(     mesh%ac1:mesh%ac2  ) = 0._dp
    ice%Qp_GL_Ac(       mesh%ac1:mesh%ac2  ) = 0._dp
    CALL sync
      
    ! Some constant terms for the Schoof grounding line flux
    exponent_Schoof  = (C%m_sliding + n_flow + 3._dp) / (C%m_sliding + 1._dp)
    C_sliding        = C%C_sliding / (sec_per_year**C%m_sliding)
    
    DO aci = 1, mesh%nAc
      IF (ice%mask_gl_Ac( aci) == 0) CYCLE

      ! The two vertices spanning this Ac connection
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)

      ! Interpolate ice thickness and till friction angle      
      TAFi = ice%Hi( vi) - ((ice%SL( vi) - ice%Hb( vi)) * (seawater_density / ice_density))
      TAFj = ice%Hi( vj) - ((ice%SL( vj) - ice%Hb( vj)) * (seawater_density / ice_density))
      lambda_GL = TAFi / (TAFi - TAFj)

      ! Interpolate ice thickness to the grounding-line position
      Hi_GL       = (ice%Hi(            vi  ) * (1._dp - lambda_GL)) + (ice%Hi(            vj  ) * lambda_GL)
      
      ! Take the friction angle and flow factor only from the grounded side
      IF (ice%mask_sheet( vi) == 1) THEN
        phi_fric_GL = ice%phi_fric_AaAc( vi)
        A_flow_GL   = ice%A_flow_mean(   vi)
      ELSE
        phi_fric_GL = ice%phi_fric_AaAc( vj)
        A_flow_GL   = ice%A_flow_mean(   vj)
      END IF

      ! The constant factor in the Schoof (2007) grounding line flux solution:
      factor_Schoof = ( A_flow_GL &
                       * ((ice_density * grav) ** (n_flow + 1._dp)) &
                       * ((1._dp - (ice_density / seawater_density)) ** n_flow) &
                       / ((4._dp**n_flow) * C_sliding) ) &
                       ** (1._dp / (REAL(C%m_sliding,dp) + 1._dp))

      ! The constant factor in the Tsai et al. (2015) grounding line flux solution:
      factor_Tsai   = ( 8._dp * Q0 * A_flow_GL &
                       * ((ice_density * grav) ** n_flow) &
                       * ((1._dp - (ice_density / seawater_density)) ** (REAL(n_flow,dp)-1._dp)) &
                       / (4._dp**REAL(n_flow,dp)) )

      ! Calculate grounding line flux
      IF (C%choice_sliding_law == 'Coulomb_regularised') THEN
        ice%Qabs_GL_Ac( aci) = factor_Tsai   * (Hi_GL**(REAL(n_flow,dp) + 2._dp)) / TAN(phi_fric_GL * (pi / 180._dp))
!      ELSEIF     (C%choice_sliding_law == 'Weertman') THEN
!        ice%Qabs_GL_Ac( aci) = factor_Schoof * Hi_GL**exponent_Schoof
      ELSE
        IF (par%master) WRITE(0,*) 'sliding law "'//TRIM(C%choice_sliding_law)//'" not implemented in calculate_GL_flux!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      ! The flux is outward-perpendicular to the grounding line, and therefore antiparallel
      ! to the gradient of the thickness-above-flotation
      Fx = -(ice%dHi_dx_Ac( aci) - ((ice%dSL_dx_Ac( aci) - ice%dHb_dx_Ac( aci)) * (seawater_density / ice_density)))
      Fy = -(ice%dHi_dy_Ac( aci) - ((ice%dSL_dy_Ac( aci) - ice%dHb_dy_Ac( aci)) * (seawater_density / ice_density)))
      F  = NORM2([Fx,Fy])
      Fx = Fx / F
      Fy = Fy / F

      ! Analytical ice velocities (antiparallel to TAF gradient)
      ice%Ux_SSA_Ac( aci) = ice%Qabs_GL_Ac( aci) * Fx / Hi_GL
      ice%Uy_SSA_Ac( aci) = ice%Qabs_GL_Ac( aci) * Fy / Hi_GL
      
      ! Rotate flux-perpendicular-to-the-grounding-line to get flux-parallel-to-the-Ac-connection (i.e. flux-from-vi-to-vj
      Dx = mesh%V( vj,1) - mesh%V( vi,1)
      Dy = mesh%V( vj,2) - mesh%V( vi,2)
      D  = NORM2([Dx,Dy])
      Dx = Dx / D
      Dy = Dy / D
      
      ice%Qp_GL_Ac( aci) = ice%Qabs_GL_Ac( aci) * (Dx*Fx + Dy*Fy)

    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync
    
  END SUBROUTINE calculate_GL_flux
  
  ! Administration: allocation, initialisation, and remapping
  SUBROUTINE initialise_ice_model( mesh, ice, init)
    ! Allocate shared memory for all the data fields of the ice dynamical module, and
    ! initialise some of them    
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    
    IF (par%master) WRITE(0,*) '  Initialising ice dynamics model...'
    
    ! Allocate shared memory
    CALL allocate_ice_model( mesh, ice)
        
    ! Initialise with data from initial file (works for both "raw" initial data and model restarts)
    ice%Hi( mesh%v1:mesh%v2) = init%Hi( mesh%v1:mesh%v2)
    ice%Hb( mesh%v1:mesh%v2) = init%Hb( mesh%v1:mesh%v2)
    ice%Hs( mesh%v1:mesh%v2) = MAX(ice%SL( mesh%v1:mesh%v2), ice%Hb( mesh%v1:mesh%v2) + ice%Hi( mesh%v1:mesh%v2))
    CALL sync
    
    ! Initialise masks and slopes
    CALL update_general_ice_model_data( mesh, ice, C%start_time_of_run)
    
    ! Geothermal heat flux
    IF     (C%choice_geothermal_heat_flux == 'constant') THEN
      ice%GHF( mesh%v1:mesh%v2) = C%constant_geothermal_heat_flux
    ELSEIF (C%choice_geothermal_heat_flux == 'spatial') THEN
      CALL map_geothermal_heat_flux_to_mesh( mesh, ice)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_geothermal_heat_flux "', TRIM(C%choice_geothermal_heat_flux), '" not implemented in initialise_ice_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! If we're doing a restart, initialise SSA velocities too
    IF (C%is_restart) THEN
      ice%U_SSA( mesh%v1:mesh%v2) = init%U_SSA( mesh%v1:mesh%v2)
      ice%V_SSA( mesh%v1:mesh%v2) = init%V_SSA( mesh%v1:mesh%v2)
    END IF
    
  END SUBROUTINE initialise_ice_model
  SUBROUTINE map_geothermal_heat_flux_to_mesh( mesh, ice)
  
    USE data_types_module,          ONLY: type_remapping_latlon2mesh
    USE forcing_module,             ONLY: forcing
    USE mesh_mapping_module,        ONLY: create_remapping_arrays_glob_mesh, map_latlon2mesh_2D, deallocate_remapping_arrays_glob_mesh
    USE mesh_help_functions_module, ONLY: partition_list
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    TYPE(type_remapping_latlon2mesh)                   :: map
    
    ! Calculate mapping arrays
    CALL create_remapping_arrays_glob_mesh( mesh, forcing%grid_ghf, map)
    
    ! Map global climate data to the mesh
    CALL map_latlon2mesh_2D( mesh, map, forcing%ghf_ghf, ice%GHF)
    
    ! Deallocate mapping arrays
    CALL deallocate_remapping_arrays_glob_mesh( map)
    
  END SUBROUTINE map_geothermal_heat_flux_to_mesh
  SUBROUTINE allocate_ice_model( mesh, ice)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.    
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
      
    ! Allocate memory
    ! ===============
    
    ! Basic data - ice thickness, bedrock height, surface height, mask, 3D ice velocities and temperature
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%Hi,                     ice%wHi                   )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Hi_Ac,                  ice%wHi_Ac                )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%Hi_prev,                ice%wHi_prev              )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%Hb,                     ice%wHb                   )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Hb_Ac,                  ice%wHb_Ac                )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHb,                    ice%wdHb                  )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%Hs,                     ice%wHs                   )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Hs_Ac,                  ice%wHs_Ac                )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%SL,                     ice%wSL                   )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%SL_Ac,                  ice%wSL_Ac                )
    CALL allocate_shared_dp_2D(  mesh%nV,   C%nZ, ice%Ti,                     ice%wTi                   ) 
    CALL allocate_shared_dp_2D(  mesh%nAc,  C%nZ, ice%Ti_Ac,                  ice%wTi_Ac                )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%U_SIA,                  ice%wU_SIA                )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%V_SIA,                  ice%wV_SIA                )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Up_SIA_Ac,              ice%wUp_SIA_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Uo_SIA_Ac,              ice%wUo_SIA_Ac            ) 
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Ux_SIA_Ac,              ice%wUx_SIA_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Uy_SIA_Ac,              ice%wUy_SIA_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%U_SSA,                  ice%wU_SSA                )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%V_SSA,                  ice%wV_SSA                )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Up_SSA_Ac,              ice%wUp_SSA_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Uo_SSA_Ac,              ice%wUo_SSA_Ac            ) 
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Ux_SSA_Ac,              ice%wUx_SSA_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Uy_SSA_Ac,              ice%wUy_SSA_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%U_surf,                 ice%wU_surf               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%V_surf,                 ice%wV_surf               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%U_base,                 ice%wU_base               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%V_base,                 ice%wV_base               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%U_vav,                  ice%wU_vav                )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%V_vav,                  ice%wV_vav                )
    
    ! Different masks
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask_land,              ice%wmask_land            )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_land_Ac,           ice%wmask_land_Ac         )
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask_ocean,             ice%wmask_ocean           )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_ocean_Ac,          ice%wmask_ocean_Ac        )
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask_lake,              ice%wmask_lake            )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_lake_Ac,           ice%wmask_lake_Ac         )
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask_ice,               ice%wmask_ice             )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_ice_Ac,            ice%wmask_ice_Ac          )
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask_sheet,             ice%wmask_sheet           )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_sheet_Ac,          ice%wmask_sheet_Ac        )
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask_shelf,             ice%wmask_shelf           )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_shelf_Ac,          ice%wmask_shelf_Ac        )
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask_coast,             ice%wmask_coast           )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_coast_Ac,          ice%wmask_coast_Ac        )
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask_margin,            ice%wmask_margin          )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_margin_Ac,         ice%wmask_margin_Ac       )
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask_gl,                ice%wmask_gl              )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_gl_Ac,             ice%wmask_gl_Ac           )
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask_cf,                ice%wmask_cf              )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_cf_Ac,             ice%wmask_cf_Ac           )
    CALL allocate_shared_int_1D( mesh%nV,         ice%mask,                   ice%wmask                 )
    CALL allocate_shared_int_1D( mesh%nAc,        ice%mask_Ac,                ice%wmask_Ac              )
    
    ! Ice physical properties
    CALL allocate_shared_dp_2D(  mesh%nV,   C%nZ, ice%A_flow,                 ice%wA_flow               )
    CALL allocate_shared_dp_2D(  mesh%nAc,  C%nZ, ice%A_flow_Ac,              ice%wA_flow_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%A_flow_mean,            ice%wA_flow_mean          )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%A_flow_mean_Ac,         ice%wA_flow_mean_Ac       )
    CALL allocate_shared_dp_2D(  mesh%nV,   C%nZ, ice%Ti_pmp,                 ice%wTi_pmp               )
    CALL allocate_shared_dp_2D(  mesh%nV,   C%nZ, ice%Cpi,                    ice%wCpi                  ) 
    CALL allocate_shared_dp_2D(  mesh%nV,   C%nZ, ice%Ki,                     ice%wKi                   ) 
    
    ! Spatial derivatives and curvatures
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHi_dx,                 ice%wdHi_dx               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHi_dy,                 ice%wdHi_dy               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHi_dt,                 ice%wdHi_dt               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHb_dx,                 ice%wdHb_dx               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHb_dy,                 ice%wdHb_dy               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHb_dt,                 ice%wdHb_dt               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHs_dx,                 ice%wdHs_dx               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHs_dy,                 ice%wdHs_dy               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHs_dt,                 ice%wdHs_dt               )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHs_dx_shelf,           ice%wdHs_dx_shelf         )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dHs_dy_shelf,           ice%wdHs_dy_shelf         )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHs_dx_shelf_Ac,        ice%wdHs_dx_shelf_Ac      )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHs_dy_shelf_Ac,        ice%wdHs_dy_shelf_Ac      )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHi_dp_Ac,              ice%wdHi_dp_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHi_do_Ac,              ice%wdHi_do_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHi_dx_Ac,              ice%wdHi_dx_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHi_dy_Ac,              ice%wdHi_dy_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHb_dp_Ac,              ice%wdHb_dp_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHb_do_Ac,              ice%wdHb_do_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHb_dx_Ac,              ice%wdHb_dx_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHb_dy_Ac,              ice%wdHb_dy_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHs_dp_Ac,              ice%wdHs_dp_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHs_do_Ac,              ice%wdHs_do_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHs_dx_Ac,              ice%wdHs_dx_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dHs_dy_Ac,              ice%wdHs_dy_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dSL_dp_Ac,              ice%wdSL_dp_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dSL_do_Ac,              ice%wdSL_do_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dSL_dx_Ac,              ice%wdSL_dx_Ac            )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%dSL_dy_Ac,              ice%wdSL_dy_Ac            )
    
    ! Zeta derivatives
    CALL allocate_shared_dp_2D(  mesh%nV,   C%nZ, ice%dzeta_dt,               ice%wdzeta_dt             )
    CALL allocate_shared_dp_2D(  mesh%nV,   C%nZ, ice%dzeta_dx,               ice%wdzeta_dx             )
    CALL allocate_shared_dp_2D(  mesh%nV,   C%nZ, ice%dzeta_dy,               ice%wdzeta_dy             )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%dzeta_dz,               ice%wdzeta_dz             )
    
    ! Ice dynamics - SIA
    CALL allocate_shared_dp_2D(  mesh%nAc,  C%nZ, ice%D_SIA_3D_Ac,            ice%wD_SIA_3D_Ac          )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%D_SIA_Ac,               ice%wD_SIA_Ac             )
    CALL allocate_shared_dp_1D(  mesh%nV,         ice%D_SIA,                  ice%wD_SIA                )
    
    ! Ice dynamics - SSA
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%Hi_AaAc,                ice%wHi_AaAc              )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%Hb_AaAc,                ice%wHb_AaAc              )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%SL_AaAc,                ice%wSL_AaAc              )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%dHs_dx_shelf_AaAc,      ice%wdHs_dx_shelf_AaAc    )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%dHs_dy_shelf_AaAc,      ice%wdHs_dy_shelf_AaAc    )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%A_flow_mean_AaAc,       ice%wA_flow_mean_AaAc     )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%phi_fric_AaAc,          ice%wphi_fric_AaAc        )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%tau_c_AaAc,             ice%wtau_c_AaAc           )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%dU_SSA_dx_AaAc,         ice%wdU_SSA_dx_AaAc       )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%dU_SSA_dy_AaAc,         ice%wdU_SSA_dy_AaAc       )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%dV_SSA_dx_AaAc,         ice%wdV_SSA_dx_AaAc       )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%dV_SSA_dy_AaAc,         ice%wdV_SSA_dy_AaAc       )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%eta_AaAc,               ice%weta_AaAc             )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%N_AaAc,                 ice%wN_AaAc               )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%N_AaAc_prev,            ice%wN_AaAc_prev          )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%S_AaAc,                 ice%wS_AaAc               )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%eu_i_AaAc,              ice%weu_i_AaAc            )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%ev_i_AaAc,              ice%wev_i_AaAc            )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%RHSx_AaAc,              ice%wRHSx_AaAc            )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%RHSy_AaAc,              ice%wRHSy_AaAc            )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%LHSx_AaAc,              ice%wLHSx_AaAc            )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%LHSy_AaAc,              ice%wLHSy_AaAc            )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%resU_AaAc,              ice%wresU_AaAc            )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%resV_AaAc,              ice%wresV_AaAc            )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%U_SSA_AaAc,             ice%wU_SSA_AaAc           )
    CALL allocate_shared_dp_1D(  mesh%nVAaAc,     ice%V_SSA_AaAc,             ice%wV_SSA_AaAc           )
    
    ! Ice dynamics - SSA - semi-analytical GL flux
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Qabs_GL_Ac,             ice%wQabs_GL_Ac           )
    CALL allocate_shared_dp_1D(  mesh%nAc,        ice%Qp_GL_Ac,               ice%wQp_GL_Ac             )
    
    ! Ice dynamics - ice thickness calculation
    CALL allocate_shared_dp_2D(  mesh%nV, mesh%nC_mem,   ice%dVi_in,         ice%wdVi_in               )
    CALL allocate_shared_dp_2D(  mesh%nV, mesh%nC_mem,   ice%dVi_out,        ice%wdVi_out              ) 
    
    ! Thermodynamics
    CALL allocate_shared_dp_1D(  mesh%nV,        ice%frictional_heating,     ice%wfrictional_heating   )
    CALL allocate_shared_dp_1D(  mesh%nV,        ice%GHF,                    ice%wGHF                  )
    CALL allocate_shared_dp_2D(  mesh%nV,  C%nZ, ice%U_3D,                   ice%wU_3D                 )
    CALL allocate_shared_dp_2D(  mesh%nV,  C%nZ, ice%V_3D,                   ice%wV_3D                 )
    CALL allocate_shared_dp_2D(  mesh%nV,  C%nZ, ice%w_3D,                   ice%wW_3D                 )
    
    ! Mesh adaptation data
    CALL allocate_shared_dp_1D(  mesh%nV,        ice%surf_curv,              ice%wsurf_curv            )
    CALL allocate_shared_dp_1D(  mesh%nV,        ice%log_velocity,           ice%wlog_velocity         )
    
  END SUBROUTINE allocate_ice_model
  SUBROUTINE remap_ice_model( mesh_old, mesh_new, map, ice, PD)
    ! Remap or reallocate all the data fields
  
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_PD_data_fields),           INTENT(IN)    :: PD
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
        
    ! The only fields that actually need to be mapped. The rest only needs memory reallocation.
    CALL remap_field_dp(         mesh_old, mesh_new, map, ice%Hi,                  ice%wHi,                  'cons_1st_order')
    CALL remap_field_dp(         mesh_old, mesh_new, map, ice%Hi_prev,             ice%wHi_prev,             'cons_1st_order')
    CALL remap_field_dp_3D(      mesh_old, mesh_new, map, ice%Ti,                  ice%wTi,                  'cons_1st_order')
!    CALL remap_field_dp(         mesh_old, mesh_new, map, ice%U_SSA,               ice%wU_SSA,               'trilin'        )
!    CALL remap_field_dp(         mesh_old, mesh_new, map, ice%V_SSA,               ice%wV_SSA,               'trilin'        )
    
    ! Remap bedrock change and add up to PD to prevent accumulation of numerical diffusion
    CALL remap_field_dp(         mesh_old, mesh_new, map, ice%dHb,                 ice%wdHb,                 'cons_2nd_order')
    CALL reallocate_field_dp(    mesh_new%nV,  ice%Hb,                     ice%wHb                       )
    ice%Hb( mesh_new%v1:mesh_new%v2) = PD%Hb( mesh_new%v1:mesh_new%v2) + ice%dHb( mesh_new%v1:mesh_new%v2)
      
!    ! Slightly reduce SSA velocities to prevent numerical instability after updating the mesh
!    ice%U_SSA( mesh_new%v1:mesh_new%v2) = ice%U_SSA( mesh_new%v1:mesh_new%v2) * 0.9_dp
!    ice%V_SSA( mesh_new%v1:mesh_new%v2) = ice%V_SSA( mesh_new%v1:mesh_new%v2) * 0.9_dp 
    
    ! Geothermal heat flux
    CALL reallocate_field_dp( mesh_new%nV, ice%GHF, ice%wGHF)
    IF     (C%choice_geothermal_heat_flux == 'constant') THEN
      ice%GHF( mesh_new%v1:mesh_new%v2) = C%constant_geothermal_heat_flux
    ELSEIF (C%choice_geothermal_heat_flux == 'spatial') THEN
      CALL map_geothermal_heat_flux_to_mesh( mesh_new, ice)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_geothermal_heat_flux "', TRIM(C%choice_geothermal_heat_flux), '" not implemented in remap_ice_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF       
    
    ! Simple memory reallocation for all the rest
    ! ===========================================
    
    ! Basic data - ice thickness, bedrock & surface elevation, sea level (geoid elevation), englacial temperature, and ice velocities
!   CALL reallocate_field_dp(    mesh_new%nV,  ice%Hi,                     ice%wHi                       )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Hi_Ac,                  ice%wHi_Ac                    )
!   CALL reallocate_field_dp(    mesh_new%nV,  ice%Hi_prev,                ice%wHi_prev                  )
!    CALL reallocate_field_dp(    mesh_new%nV,  ice%Hb,                     ice%wHb                       )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Hb_Ac,                  ice%wHb_Ac                    )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%Hs,                     ice%wHs                       )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Hs_Ac,                  ice%wHs_Ac                    )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%SL,                     ice%wSL                       )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%SL_Ac,                  ice%wSL_Ac                    )
!   CALL reallocate_field_dp_3D( mesh_new%nV,  ice%Ti,                     ice%wTi                 , C%nZ)
    CALL reallocate_field_dp_3D( mesh_new%nAc, ice%Ti_Ac,                  ice%wTi_Ac              , C%nZ)
    CALL reallocate_field_dp(    mesh_new%nV,  ice%U_SIA,                  ice%wU_SIA                    )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%V_SIA,                  ice%wV_SIA                    )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Ux_SIA_Ac,              ice%wUx_SIA_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Uy_SIA_Ac,              ice%wUy_SIA_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Up_SIA_Ac,              ice%wUp_SIA_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Uo_SIA_Ac,              ice%wUo_SIA_Ac                )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%U_SSA,                  ice%wU_SSA                    )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%V_SSA,                  ice%wV_SSA                    )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Ux_SSA_Ac,              ice%wUx_SSA_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Uy_SSA_Ac,              ice%wUy_SSA_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Up_SSA_Ac,              ice%wUp_SSA_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Uo_SSA_Ac,              ice%wUo_SSA_Ac                )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%U_surf,                 ice%wU_surf                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%V_surf,                 ice%wV_surf                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%U_base,                 ice%wU_base                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%V_base,                 ice%wV_base                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%U_vav,                  ice%wU_vav                    )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%V_vav,                  ice%wV_vav                    )
    
    CALL map_Aa_to_Ac( mesh_new, ice%U_SSA, ice%Ux_SSA_Ac)
    CALL map_Aa_to_Ac( mesh_new, ice%V_SSA, ice%Uy_SSA_Ac)
     
    ! Different masks
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask_land,              ice%wmask_land                ) 
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_land_Ac,           ice%wmask_land_Ac             ) 
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask_ocean,             ice%wmask_ocean               ) 
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_ocean_Ac,          ice%wmask_ocean_Ac            ) 
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask_lake,              ice%wmask_lake                )  
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_lake_Ac,           ice%wmask_lake_Ac             ) 
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask_ice,               ice%wmask_ice                 ) 
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_ice_Ac,            ice%wmask_ice_Ac              ) 
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask_sheet,             ice%wmask_sheet               )  
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_sheet_Ac,          ice%wmask_sheet_Ac            )
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask_shelf,             ice%wmask_shelf               )
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_shelf_Ac,          ice%wmask_shelf_Ac            ) 
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask_coast,             ice%wmask_coast               ) 
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_coast_Ac,          ice%wmask_coast_Ac            )   
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask_margin,            ice%wmask_margin              )
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_margin_Ac,         ice%wmask_margin_Ac           )    
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask_gl,                ice%wmask_gl                  ) 
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_gl_Ac,             ice%wmask_gl_Ac               )  
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask_cf,                ice%wmask_cf                  )
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_cf_Ac,             ice%wmask_cf_Ac               )
    CALL reallocate_field_int(   mesh_new%nV,  ice%mask,                   ice%wmask                     ) 
    CALL reallocate_field_int(   mesh_new%nAc, ice%mask_Ac,                ice%wmask_Ac                  )
    
    ! Ice physical properties
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%A_flow,                 ice%wA_flow             , C%nZ)
    CALL reallocate_field_dp_3D( mesh_new%nAc, ice%A_flow_Ac,              ice%wA_flow_Ac          , C%nZ)
    CALL reallocate_field_dp(    mesh_new%nV,  ice%A_flow_mean,            ice%wA_flow_mean              )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%A_flow_mean_Ac,         ice%wA_flow_mean_Ac           )
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%Ti_pmp,                 ice%wTi_pmp             , C%nZ)
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%Cpi,                    ice%wCpi                , C%nZ)
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%Ki,                     ice%wKi                 , C%nZ)
    
    ! Spatial derivatives
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHi_dx,                 ice%wdHi_dx                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHi_dy,                 ice%wdHi_dy                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHi_dt,                 ice%wdHi_dt                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHb_dx,                 ice%wdHb_dx                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHb_dy,                 ice%wdHb_dy                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHb_dt,                 ice%wdHb_dt                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHs_dx,                 ice%wdHs_dx                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHs_dy,                 ice%wdHs_dy                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHs_dt,                 ice%wdHs_dt                   )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHs_dx_shelf,           ice%wdHs_dx_shelf             )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dHs_dy_shelf,           ice%wdHs_dy_shelf             )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHs_dx_shelf_Ac,        ice%wdHs_dx_shelf_Ac          )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHs_dy_shelf_Ac,        ice%wdHs_dy_shelf_Ac          )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHi_dx_Ac,              ice%wdHi_dx_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHi_dy_Ac,              ice%wdHi_dy_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHi_dp_Ac,              ice%wdHi_dp_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHi_do_Ac,              ice%wdHi_do_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHb_dx_Ac,              ice%wdHb_dx_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHb_dy_Ac,              ice%wdHb_dy_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHb_dp_Ac,              ice%wdHb_dp_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHb_do_Ac,              ice%wdHb_do_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHs_dx_Ac,              ice%wdHs_dx_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHs_dy_Ac,              ice%wdHs_dy_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHs_dp_Ac,              ice%wdHs_dp_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dHs_do_Ac,              ice%wdHs_do_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dSL_dx_Ac,              ice%wdSL_dx_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dSL_dy_Ac,              ice%wdSL_dy_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dSL_dp_Ac,              ice%wdSL_dp_Ac                )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%dSL_do_Ac,              ice%wdSL_do_Ac                )
    
    ! Zeta derivatives
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%dzeta_dt,               ice%wdzeta_dt           , C%nZ)
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%dzeta_dx,               ice%wdzeta_dx           , C%nZ)
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%dzeta_dy,               ice%wdzeta_dy           , C%nZ)
    CALL reallocate_field_dp(    mesh_new%nV,  ice%dzeta_dz,               ice%wdzeta_dz                 )
    
    ! Ice dynamics - SIA
    CALL reallocate_field_dp_3D( mesh_new%nAc, ice%D_SIA_3D_Ac,            ice%wD_SIA_3D_Ac        , C%nZ)
    CALL reallocate_field_dp(    mesh_new%nAc, ice%D_SIA_Ac,               ice%wD_SIA_Ac                 )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%D_SIA,                  ice%wD_SIA                    )
    
    ! Ice dynamics - SSA
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%Hi_AaAc,                ice%wHi_AaAc                  )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%Hb_AaAc,                ice%wHb_AaAc                  )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%SL_AaAc,                ice%wSL_AaAc                  )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%dHs_dx_shelf_AaAc,      ice%wdHs_dx_shelf_AaAc        )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%dHs_dy_shelf_AaAc,      ice%wdHs_dy_shelf_AaAc        )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%A_flow_mean_AaAc,       ice%wA_flow_mean_AaAc         )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%phi_fric_AaAc,          ice%wphi_fric_AaAc            )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%tau_c_AaAc,             ice%wtau_c_AaAc               )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%dU_SSA_dx_AaAc,         ice%wdU_SSA_dx_AaAc           )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%dU_SSA_dy_AaAc,         ice%wdU_SSA_dy_AaAc           )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%dV_SSA_dx_AaAc,         ice%wdV_SSA_dx_AaAc           )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%dV_SSA_dy_AaAc,         ice%wdV_SSA_dy_AaAc           )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%eta_AaAc,               ice%weta_AaAc                 )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%N_AaAc,                 ice%wN_AaAc                   )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%N_AaAc_prev,            ice%wN_AaAc_prev              )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%S_AaAc,                 ice%wS_AaAc                   )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%eu_i_AaAc,              ice%weu_i_AaAc                )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%ev_i_AaAc,              ice%wev_i_AaAc                )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%RHSx_AaAc,              ice%wRHSx_AaAc                )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%RHSy_AaAc,              ice%wRHSy_AaAc                )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%LHSx_AaAc,              ice%wLHSx_AaAc                )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%LHSy_AaAc,              ice%wLHSy_AaAc                )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%resU_AaAc,              ice%wresU_AaAc                )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%resV_AaAc,              ice%wresV_AaAc                )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%U_SSA_AaAc,             ice%wU_SSA_AaAc               )
    CALL reallocate_field_dp(    mesh_new%nVAaAc,  ice%V_SSA_AaAc,             ice%wV_SSA_AaAc               )
    
    ! Ice dynamics- SSA - semi-analytical GL flux
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Qabs_GL_Ac,             ice%wQabs_GL_Ac               )
    CALL reallocate_field_dp(    mesh_new%nAc, ice%Qp_GL_Ac,               ice%wQp_GL_Ac                 )
    
    ! Ice dynamics - ice thickness calculation
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%dVi_in,                 ice%wdVi_in,                  mesh_new%nC_mem)
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%dVi_out,                ice%wdVi_out,                 mesh_new%nC_mem)
    
    ! Thermodynamics
    CALL reallocate_field_dp(    mesh_new%nV,  ice%frictional_heating,     ice%wfrictional_heating       )
    !CALL reallocate_field_dp(    mesh_new%nV,  ice%GHF,                    ice%wGHF                      )
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%U_3D,                   ice%wU_3D,                C%nZ)
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%V_3D,                   ice%wV_3D,                C%nZ)
    CALL reallocate_field_dp_3D( mesh_new%nV,  ice%W_3D,                   ice%wW_3D,                C%nZ)
    
    ! Mesh adaptation data
    CALL reallocate_field_dp(    mesh_new%nV,  ice%surf_curv,              ice%wsurf_curv                )
    CALL reallocate_field_dp(    mesh_new%nV,  ice%log_velocity,           ice%wlog_velocity             )
    
  END SUBROUTINE remap_ice_model  
  
END MODULE ice_dynamics_module
