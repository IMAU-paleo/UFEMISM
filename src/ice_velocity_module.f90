MODULE ice_velocity_module

  ! Contains all the routines needed to calculate instantaneous ice velocities for the
  ! current modelled ice-sheet geometry.

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
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_sparse_matrix_CSR
  USE mesh_operators_module,           ONLY: map_a_to_c_2D, ddx_a_to_c_2D, ddy_a_to_c_2D, map_a_to_c_3D, &
                                             move_aca_to_a_2D, move_acc_to_c_2D, ddx_a_to_a_2D, ddy_a_to_a_2D, &
                                             ddx_c_to_a_3D, ddy_c_to_a_3D, map_a_to_ac_2D, map_a_to_ac_3D, &
                                             ddx_a_to_ac_2D, ddy_a_to_ac_2D, ddx_ac_to_ac_2D, ddy_ac_to_ac_2D, &
                                             move_ac_to_acu_2D, move_ac_to_acv_2D, move_acu_to_ac_2D, move_acv_to_ac_2D
  USE utilities_module,                ONLY: vertical_integration_from_bottom_to_zeta, vertical_average, &
                                             vertical_integrate, multiply_matrix_matrix_CSR, convert_vector_to_diag_CSR, &
                                             add_matrix_matrix_CSR, multiply_matrix_rows_with_vector, &
                                             deallocate_matrix_CSR, overwrite_rows_CSR, solve_matrix_equation_CSR, &
                                             allocate_matrix_CSR_dist, finalise_matrix_CSR_dist
  USE basal_conditions_and_sliding_module, ONLY: calc_basal_conditions, calc_sliding_law
  USE netcdf_module,                   ONLY: write_CSR_matrix_to_NetCDF

  IMPLICIT NONE
  
CONTAINS
  
  ! The different velocity equation solvers that can be called from run_ice_model
  SUBROUTINE solve_SIA(  mesh, ice)
    ! Calculate ice velocities using the SIA.
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
      
    ! Local variables:
    INTEGER                                            :: aci, vi, vj
    REAL(dp), DIMENSION(:    ), POINTER                ::  Hi_c,  dHs_dx_c,  dHs_dy_c
    INTEGER                                            :: wHi_c, wdHs_dx_c, wdHs_dy_c
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  A_flow_3D_c
    INTEGER                                            :: wA_flow_3D_c
    REAL(dp)                                           :: D_0
    REAL(dp), DIMENSION(C%nZ)                          :: D_prof, D_deformation
    REAL(dp), DIMENSION(C%nZ)                          :: D_SIA_3D
    REAL(dp), PARAMETER                                :: D_uv_3D_cutoff = -1E5_dp
    
    ! Safety
    IF (.NOT. (C%choice_ice_dynamics == 'SIA' .OR. C%choice_ice_dynamics == 'SIA/SSA')) THEN
      IF (par%master) WRITE(0,*) 'ERROR - solve_SIA should only be called when choice_ice_dynamics is set to SIA or SIA/SSA!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nAc,       Hi_c       , wHi_c       )
    CALL allocate_shared_dp_1D( mesh%nAc,       dHs_dx_c   , wdHs_dx_c   )
    CALL allocate_shared_dp_1D( mesh%nAc,       dHs_dy_c   , wdHs_dy_c   )
    CALL allocate_shared_dp_2D( mesh%nAc, C%nz, A_flow_3D_c, wA_flow_3D_c)
    
    ! Get ice thickness, surface slopes, and ice flow factor on the staggered mesh
    CALL map_a_to_c_2D( mesh, ice%Hi_a       , Hi_c       )
    CALL ddx_a_to_c_2D( mesh, ice%Hs_a       , dHs_dx_c   )
    CALL ddy_a_to_c_2D( mesh, ice%Hs_a       , dHs_dy_c   )
    CALL map_a_to_c_3D( mesh, ice%A_flow_3D_a, A_flow_3D_c)

    ! Calculate 3D horizontal velocities (on the c-part of the ac-grid)
    DO aci = mesh%ci1, mesh%ci2
      
      vi   = mesh%Aci( aci,1)
      vj   = mesh%Aci( aci,2)
      
      IF (ice%mask_sheet_a( vi) == 0 .AND. ice%mask_sheet_a( vj) == 0) THEN
        ! No ice at either regular vertex; set velocity to zero
  
        ice%u_3D_SIA_c( aci,:) = 0._dp
        ice%v_3D_SIA_c( aci,:) = 0._dp
      
      ELSE
        
        D_0           = (ice_density * grav * Hi_c( aci))**C%n_flow * ((dHs_dx_c( aci)**2 + dHs_dy_c( aci)**2))**((C%n_flow - 1._dp) / 2._dp)
        D_prof        = A_flow_3D_c( aci,:) * C%zeta**C%n_flow
        CALL vertical_integration_from_bottom_to_zeta( D_prof, D_deformation)
        D_deformation = 2._dp * Hi_c( aci) * D_deformation
        D_SIA_3D      = MAX( D_0 * D_deformation, D_uv_3D_cutoff)
       
        ice%u_3D_SIA_c( aci,:) = D_SIA_3D * dHs_dx_c( aci)
        ice%v_3D_SIA_c( aci,:) = D_SIA_3D * dHs_dy_c( aci)
      
      END IF
      
    END DO
    CALL sync
    
    ! Calculate secondary velocities (surface, base, etc.)
    CALL calc_secondary_velocities( mesh, ice)
    
    ! Clean up after yourself
    CALL deallocate_shared( wHi_c       )
    CALL deallocate_shared( wdHs_dx_c   )
    CALL deallocate_shared( wdHs_dy_c   )
    CALL deallocate_shared( wA_flow_3D_c)

  END SUBROUTINE solve_SIA
  SUBROUTINE solve_SSA(  mesh, ice)
    ! Calculate ice velocities using the SSA
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    IF (C%include_SSADIVA_crossterms) THEN
      IF (par%master) WRITE(0,*) 'solve_SSA - ERROR: full SSA solver not yet implemented!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSE
      CALL solve_SSA_sans(  mesh, ice)
    END IF
    
    ! Calculate secondary velocities (surface, base, etc.)
    CALL calc_secondary_velocities( mesh, ice)
    
  END SUBROUTINE solve_SSA
  SUBROUTINE solve_SSA_sans(  mesh, ice)
    ! Calculate ice velocities using the SSA sans the "cross-terms"
    ! 
    ! The full SSA reads:
    ! 
    !   d/dx[ 2N (2*du/dx + dv/dy)] + d/dy[ N (du/dy + dv/dx)] - beta*u = -taud_x
    !   d/dy[ 2N (2*dv/dy + du/dx)] + d/dx[ N (dv/dx + du/dy)] - beta*v = -taud_y
    !
    ! Using the chain rule, this expands to:
    !
    !   4*N*d2u/dx2 + 4*dN/dx*du/dx + 2*N*d2v/dxdy + 2*dN/dx*dv/dy + N*d2u/dy2 + dN/dy*du/dy + N*d2v/dxdy + dN/dy*dv/dx - beta*u = -taud_x
    !   4*N*d2v/dy2 + 4*dN/dy*dv/dy + 2*N*d2u/dxdy + 2*dN/dy*du/dx + N*d2v/dx2 + dN/dx*dv/dx + N*d2u/dxdy + dN/dx*du/dy - beta*v = -taud_y
    !
    ! By neglecting all the "cross-terms" involving gradients of N, this simplifies to:
    !
    !   4*N*d2u/dx2 + N*d2u/dy2 + 3*N*d2v/dxdy - beta*u = -taud_x
    !   4*N*d2v/dy2 + N*d2v/dx2 + 3*N*d2u/dxdy - beta*v = -taud_y
    ! 
    ! Lastly, the left and right sides of the equations can be divided by N to yield:
    ! 
    ! 4*d2u/dx2 + d2u/dy2 + 3*d2v/dxdy - beta*u/N = -taud_x/N
    ! 4*d2v/dy2 + d2v/dx2 + 3*d2u/dxdy - beta*v/N = -taud_y/N
    ! 
    ! By happy coincidence, this equation is much easier to solve, especially because (for unclear reasons)
    ! it doesn't require N to be defined on a staggered grid relative to u,v in order to achieve a stable solution.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    LOGICAL                                            :: set_velocities_to_zero
    LOGICAL                                            :: has_converged
    INTEGER                                            :: viscosity_iteration_i
    REAL(dp)                                           :: resid_UV
    REAL(dp)                                           :: umax_analytical, tauc_analytical
    REAL(dp), DIMENSION(:    ), POINTER                ::  u_ac_prev,  v_ac_prev
    INTEGER                                            :: wu_ac_prev, wv_ac_prev
    
  ! =========
  ! == Safety
    
    ! Check that this routine is called correctly
    IF (.NOT. (C%choice_ice_dynamics == 'SSA' .OR. C%choice_ice_dynamics == 'SIA/SSA')) THEN
      IF (par%master) WRITE(0,*) 'ERROR - solve_SSA should only be called when choice_ice_dynamics is set to SSA or SIA/SSA!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! If there's no grounded ice anywhere, don't bother
    set_velocities_to_zero = .FALSE.
    IF (SUM( ice%mask_sheet_a) == 0) set_velocities_to_zero = .TRUE.
    
    ! If we're prescribing no sliding, set velocities to zero
    IF (C%choice_sliding_law == 'no_sliding') set_velocities_to_zero = .TRUE.
    
    IF (set_velocities_to_zero) THEN
      ice%u_SSA_ac( mesh%avi1:mesh%avi2) = 0._dp
      ice%v_SSA_ac( mesh%avi1:mesh%avi2) = 0._dp
      CALL sync
      CALL calc_secondary_velocities( mesh, ice)
      RETURN
    END IF
    
  ! == Safety - end
  ! ===============
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nVAaAc, u_ac_prev, wu_ac_prev)
    CALL allocate_shared_dp_1D( mesh%nVAaAc, v_ac_prev, wv_ac_prev)
    
    ! Calculate the driving stresses taudx, taudy
    CALL calculate_driving_stress( mesh, ice)

    ! Calculate the basal yield stress tau_c
    CALL calc_basal_conditions( mesh, ice)
    
!    ! Determine sub-mesh grounded fractions for scaling the basal friction
!    CALL determine_grounded_fractions( mesh, ice)
!    
!    ! Find analytical solution for the SSA icestream experiment (used only to print numerical error to screen)
!    CALL SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, ice%A_flow_vav_a( 1,1), 0._dp, umax_analytical, tauc_analytical)
    
!    ! DENK DROM
!    IF (par%master) THEN
!      debug%dp_2D_ac_01 = ice%taudx_ac
!      debug%dp_2D_ac_02 = ice%taudy_ac
!      CALL write_to_debug_file
!    END IF
!    CALL sync
    
    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .FALSE.
    viscosity_iteration: DO WHILE (.NOT. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1
      
      ! Calculate the effective viscosity and the product term N = eta * H
      CALL calc_effective_viscosity_ac( mesh, ice, ice%u_SSA_ac, ice%v_SSA_ac)
    
!    ! DENK DROM
!    IF (par%master) THEN
!      debug%dp_2D_ac_03 = ice%visc_eff_int_ac
!      debug%dp_2D_ac_04 = ice%N_ac
!      CALL write_to_debug_file
!    END IF
!    CALL sync
            
      ! Calculate the sliding term beta
      CALL calc_sliding_term_beta( mesh, ice, ice%u_SSA_ac, ice%v_SSA_ac)
    
!    ! DENK DROM
!    IF (par%master) THEN
!      debug%dp_2D_ac_05 = ice%beta_ac
!      CALL write_to_debug_file
!    END IF
!    CALL sync
    
!      ! Apply the sub-mesh grounded fraction
!      DO i = mesh%i1, mesh%i2
!      DO j = 1, mesh%ny
!        IF (i < mesh%nx) ice%beta_eff_cx( j,i) = ice%beta_eff_cx( j,i) * ice%f_grnd_cx( j,i)**2
!        IF (j < mesh%ny) ice%beta_eff_cy( j,i) = ice%beta_eff_cy( j,i) * ice%f_grnd_cy( j,i)**2
!      END DO
!      END DO
!      CALL sync
      
      ! Store the previous solution so we can check for convergence later
      u_ac_prev( mesh%avi1:mesh%avi2) = ice%u_SSA_ac( mesh%avi1:mesh%avi2)
      v_ac_prev( mesh%avi1:mesh%avi2) = ice%v_SSA_ac( mesh%avi1:mesh%avi2)
      CALL sync
      
      ! Solve the linearised SSA
      CALL solve_SSA_sans_linearised( mesh, ice, ice%u_SSA_ac, ice%v_SSA_ac)
    
!    ! DENK DROM
!    IF (par%master) THEN
!      debug%dp_2D_ac_06 = ice%u_SSA_ac
!      debug%dp_2D_ac_07 = ice%v_SSA_ac
!      CALL write_to_debug_file
!    END IF
!    CALL sync
      
      ! Apply velocity limits (both overflow and underflow) for improved stability
      CALL apply_velocity_limits( mesh, ice%u_SSA_ac, ice%v_SSA_ac)
      
      ! "relax" subsequent viscosity iterations for improved stability
      CALL relax_DIVA_visc_iterations( mesh, u_ac_prev, v_ac_prev, ice%u_SSA_ac, ice%v_SSA_ac, C%DIVA_visc_it_relax)
      
      ! Check if the viscosity iteration has converged
      CALL calc_visc_iter_UV_resid( mesh, u_ac_prev, v_ac_prev, ice%u_SSA_ac, ice%v_SSA_ac, resid_UV)
      IF (par%master) WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ': resid_UV = ', resid_UV, ', u = [', MINVAL(ice%u_SSA_ac), ' - ', MAXVAL(ice%u_SSA_ac), ']'
      
      IF (par%master .AND. C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'SSA_icestream') &
        WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ': err = ', ABS(1._dp - MAXVAL(ice%u_SSA_ac) / umax_analytical), ': resid_UV = ', resid_UV

      has_converged = .FALSE.
      IF     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) THEN
        has_converged = .TRUE.
      ELSEIF (viscosity_iteration_i >= C%DIVA_visc_it_nit) THEN
        has_converged = .TRUE.
      END IF
    
!    ! DENK DROM
!    IF (par%master) THEN
!      debug%dp_2D_ac_08 = ice%u_SSA_ac
!      debug%dp_2D_ac_09 = ice%v_SSA_ac
!      CALL write_to_debug_file
!    END IF
!    CALL sync
      
    END DO viscosity_iteration
    
    ! Clean up after yourself
    CALL deallocate_shared( wu_ac_prev)
    CALL deallocate_shared( wv_ac_prev)
    
  END SUBROUTINE solve_SSA_sans
  
  ! Calculate "secondary" velocities (surface, basal, vertically averaged, on the A-mesh, etc.)
  SUBROUTINE calc_secondary_velocities( mesh, ice)
    ! Calculate "secondary" velocities (surface, basal, vertically averaged, on the A-mesh, etc.)
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: aci, vi, k
    REAL(dp), DIMENSION(C%nz)                          :: prof
    
    IF (C%choice_ice_dynamics == 'SIA') THEN
      ! No SSA or sliding, just the SIA
      
      ! Set 3D velocity field equal to SIA answer
      ice%u_3D_c( mesh%ci1:mesh%ci2,:) = ice%u_3D_SIA_c( mesh%ci1:mesh%ci2,:)
      ice%v_3D_c( mesh%ci1:mesh%ci2,:) = ice%v_3D_SIA_c( mesh%ci1:mesh%ci2,:)
      CALL sync
      
      ! Basal velocity is zero
      ice%u_base_c( mesh%ci1:mesh%ci2) = 0._dp
      ice%v_base_c( mesh%ci1:mesh%ci2) = 0._dp
      CALL sync
    
      ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
      CALL calc_3D_vertical_velocities( mesh, ice)
      
      ! Get surface velocity from the 3D fields
      ice%u_surf_c( mesh%ci1:mesh%ci2) = ice%u_3D_c( mesh%ci1:mesh%ci2,1)
      ice%v_surf_c( mesh%ci1:mesh%ci2) = ice%v_3D_c( mesh%ci1:mesh%ci2,1)
      CALL sync
      
      ! Get vertically averaged velocities
      DO aci = mesh%ci1, mesh%ci2
        prof = ice%u_3D_c( aci,:)
        CALL vertical_average( prof, ice%u_vav_c( aci))
        prof = ice%v_3D_c( aci,:)
        CALL vertical_average( prof, ice%v_vav_c( aci))
      END DO
      CALL sync
      
      ! Map velocity components to the a-grid (as the SIA solver only defines them on the c-grid)
      CALL map_velocity_from_c_to_a_3D( mesh, ice, ice%u_3D_c,   ice%u_3D_a  )
      CALL map_velocity_from_c_to_a_3D( mesh, ice, ice%v_3D_c,   ice%v_3D_a  )
      CALL map_velocity_from_c_to_a_2D( mesh, ice, ice%u_vav_c,  ice%u_vav_a )
      CALL map_velocity_from_c_to_a_2D( mesh, ice, ice%v_vav_c,  ice%v_vav_a )
      CALL map_velocity_from_c_to_a_2D( mesh, ice, ice%u_surf_c, ice%u_surf_a)
      CALL map_velocity_from_c_to_a_2D( mesh, ice, ice%v_surf_c, ice%v_surf_a)
      
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      ! No SIA, just the SSA
      
      ! Set basal velocity equal to SSA answer
      CALL move_aca_to_a_2D( mesh, ice%u_SSA_ac, ice%u_base_a)
      CALL move_aca_to_a_2D( mesh, ice%v_SSA_ac, ice%v_base_a)
      CALL move_acc_to_c_2D( mesh, ice%u_SSA_ac, ice%u_base_c)
      CALL move_acc_to_c_2D( mesh, ice%v_SSA_ac, ice%v_base_c)
      
      ! No vertical variations in velocity
      DO vi = mesh%vi1, mesh%vi2
        ice%u_vav_a(  vi   ) = ice%u_base_a( vi)
        ice%v_vav_a(  vi   ) = ice%v_base_a( vi)
        ice%u_surf_a( vi   ) = ice%u_base_a( vi)
        ice%v_surf_a( vi   ) = ice%v_base_a( vi)
        ice%u_3D_a(   vi ,:) = ice%u_base_a( vi)
        ice%v_3D_a(   vi ,:) = ice%v_base_a( vi)
      END DO
      DO aci = mesh%ci1, mesh%ci2
        ice%u_vav_c(  aci  ) = ice%u_base_c( aci)
        ice%v_vav_c(  aci  ) = ice%v_base_c( aci)
        ice%u_surf_c( aci  ) = ice%u_base_c( aci)
        ice%v_surf_c( aci  ) = ice%v_base_c( aci)
        ice%u_3D_c(   aci,:) = ice%u_base_c( aci)
        ice%v_3D_c(   aci,:) = ice%v_base_c( aci)
      END DO
      CALL sync
      
      ! No vertical velocity anywhere
      ice%v_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
      CALL sync
      
    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      ! Hybrid SIA/SSA
      
      ! Set basal velocity equal to SSA answer
      CALL move_acc_to_c_2D( mesh, ice%u_SSA_ac, ice%u_base_c)
      CALL move_acc_to_c_2D( mesh, ice%v_SSA_ac, ice%v_base_c)
      
      ! Add the two fields together
      DO k = 1, C%nz
        ice%u_3D_c( mesh%ci1:mesh%ci2,k) = ice%u_3D_SIA_c( mesh%ci1:mesh%ci2,k) + ice%u_base_c( mesh%ci1:mesh%ci2)
        ice%v_3D_c( mesh%ci1:mesh%ci2,k) = ice%v_3D_SIA_c( mesh%ci1:mesh%ci2,k) + ice%v_base_c( mesh%ci1:mesh%ci2)
      END DO
      CALL sync
      
      ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
      CALL calc_3D_vertical_velocities( mesh, ice)
      
      ! Get surface velocity from the 3D fields
      ice%u_surf_c( mesh%ci1:mesh%ci2) = ice%u_3D_c( mesh%ci1:mesh%ci2,1)
      ice%v_surf_c( mesh%ci1:mesh%ci2) = ice%v_3D_c( mesh%ci1:mesh%ci2,1)
      CALL sync
      
      ! Get vertically averaged velocities
      DO aci = mesh%ci1, mesh%ci2
        prof = ice%u_3D_c( aci,:)
        CALL vertical_average( prof, ice%u_vav_c( aci))
        prof = ice%v_3D_c( aci,:)
        CALL vertical_average( prof, ice%v_vav_c( aci))
      END DO
      CALL sync
      
      ! Map velocities to the a-grid
      CALL move_aca_to_a_2D( mesh, ice%u_SSA_ac, ice%u_base_a)
      CALL move_aca_to_a_2D( mesh, ice%v_SSA_ac, ice%v_base_a)
      CALL map_velocity_from_c_to_a_3D( mesh, ice, ice%u_3D_SIA_c, ice%u_3D_a)
      CALL map_velocity_from_c_to_a_3D( mesh, ice, ice%v_3D_SIA_c, ice%v_3D_a)
      
      ! Add SIA and SSA contributions together
      DO k = 1, C%nz
        ice%u_3D_a( mesh%vi1:mesh%vi2,k) = ice%u_3D_a( mesh%vi1:mesh%vi2,k) + ice%u_base_a( mesh%vi1:mesh%vi2)
        ice%v_3D_a( mesh%vi1:mesh%vi2,k) = ice%v_3D_a( mesh%vi1:mesh%vi2,k) + ice%v_base_a( mesh%vi1:mesh%vi2)
      END DO
      CALL sync
      
      ! Get surface velocity from the 3D fields
      ice%u_surf_a( mesh%vi1:mesh%vi2) = ice%u_3D_a( mesh%vi1:mesh%vi2,1)
      ice%v_surf_a( mesh%vi1:mesh%vi2) = ice%v_3D_a( mesh%vi1:mesh%vi2,1)
      CALL sync
      
      ! Get vertically averaged velocities
      DO vi = mesh%vi1, mesh%vi2
        prof = ice%u_3D_a( vi,:)
        CALL vertical_average( prof, ice%u_vav_a( vi))
        prof = ice%v_3D_a( vi,:)
        CALL vertical_average( prof, ice%v_vav_a( vi))
      END DO
      CALL sync
    
    ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      ! DIVA
      
      WRITE(0,*) 'calc_secondary_velocities - DIVA: FIXME!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      
!      CALL calc_3D_vertical_velocities( mesh, ice)
!      
!      ! Get surface velocity from the 3D fields
!      ice%u_surf_c( mesh%ci1:mesh%ci2) = ice%u_3D_c( mesh%ci1:mesh%ci2,1)
!      ice%v_surf_c( mesh%ci1:mesh%ci2) = ice%v_3D_c( mesh%ci1:mesh%ci2,1)
!      CALL sync
      
    ELSE ! IF (C%choice_ice_dynamics == 'SIA')
    
      IF (par%master) WRITE(0,*) 'calc_secondary_velocities - ERROR: unknown choice_ice_dynamics "', C%choice_ice_dynamics, '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      
    END IF ! IF (C%choice_ice_dynamics == 'SIA')
    
    ! Get absolute velocities on the a-grid (only used for writing to output)
    DO vi = mesh%vi1, mesh%vi2
      ice%uabs_vav_a(  vi) = SQRT( ice%u_vav_a(  vi)**2 + ice%v_vav_a(  vi)**2)
      ice%uabs_surf_a( vi) = SQRT( ice%u_surf_a( vi)**2 + ice%v_surf_a( vi)**2)
      ice%uabs_base_a( vi) = SQRT( ice%u_base_a( vi)**2 + ice%v_base_a( vi)**2)
    END DO
    CALL sync
    
  END SUBROUTINE calc_secondary_velocities
  SUBROUTINE calc_3D_vertical_velocities( mesh, ice)
    ! Use simple conservation of mass to calculate the vertical velocity w_3D
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: vi, k
    REAL(dp), DIMENSION(:    ), POINTER                ::  dHi_dx_a,  dHi_dy_a,  dHs_dx_a,  dHs_dy_a
    INTEGER                                            :: wdHi_dx_a, wdHi_dy_a, wdHs_dx_a, wdHs_dy_a
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  du_dx_3D_a,  dv_dy_3D_a
    INTEGER                                            :: wdu_dx_3D_a, wdv_dy_3D_a
    REAL(dp)                                           :: dHbase_dx, dHbase_dy
    REAL(dp)                                           :: du_dx_k,   dv_dy_k
    REAL(dp)                                           :: du_dx_kp1, dv_dy_kp1
    REAL(dp)                                           :: w1, w2, w3, w4
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV,       dHi_dx_a  , wdHi_dx_a  )
    CALL allocate_shared_dp_1D( mesh%nV,       dHi_dy_a  , wdHi_dy_a  )
    CALL allocate_shared_dp_1D( mesh%nV,       dHs_dx_a  , wdHs_dx_a  )
    CALL allocate_shared_dp_1D( mesh%nV,       dHs_dy_a  , wdHs_dy_a  )
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, du_dx_3D_a, wdu_dx_3D_a)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, dv_dy_3D_a, wdv_dy_3D_a)
    
    ! Calculate surface & ice thickness gradients
    CALL ddx_a_to_a_2D( mesh, ice%Hs_a  , dHs_dx_a  )
    CALL ddy_a_to_a_2D( mesh, ice%Hs_a  , dHs_dy_a  )
    CALL ddx_a_to_a_2D( mesh, ice%Hi_a  , dHi_dx_a  )
    CALL ddy_a_to_a_2D( mesh, ice%Hi_a  , dHi_dy_a  )
    
    ! Calculate strain rates
    CALL ddx_c_to_a_3D( mesh, ice%u_3D_c, du_dx_3D_a)
    CALL ddy_c_to_a_3D( mesh, ice%v_3D_c, dv_dy_3D_a)
    
    ! Integrate over the incompressibility condition to find the vertical velocity profiles
    DO vi = mesh%vi1, mesh%vi2
    
      IF (ice%mask_ice_a( vi) == 0) CYCLE
      
      ! Calculate the ice basal surface slope (not the same as bedrock slope when ice is floating!)
      dHbase_dx = dHs_dx_a( vi) - dHi_dx_a( vi)
      dHbase_dy = dHs_dy_a( vi) - dHi_dy_a( vi)   
      
      ! Calculate the vertical velocity at the ice base
      IF (ice%mask_sheet_a( vi) == 1) THEN
        ice%w_3D_a( vi,C%nz) = (ice%u_3D_a( vi,C%nz) * dHbase_dx) + (ice%v_3D_a( vi,C%nz) * dHbase_dy) + ice%dHb_dt_a( vi)
      ELSE
        ice%w_3D_a( vi,C%nz) = (ice%u_3D_a( vi,C%nz) * dHbase_dx) + (ice%v_3D_a( vi,C%nz) * dHbase_dy) ! Should this include the ice thinning rate?
      END IF
                           
      ! The integrant is calculated half way the layer of integration at k+1/2. This integrant is multiplied with the layer thickness and added to the integral
      ! of all layers below, giving the integral up to and including this layer:
      DO k = C%nz - 1, 1, -1
        
        du_dx_k   = du_dx_3D_a( vi,k  )
        du_dx_kp1 = du_dx_3D_a( vi,k+1)
        dv_dy_k   = dv_dy_3D_a( vi,k  )
        dv_dy_kp1 = dv_dy_3D_a( vi,k+1)
        
        w1 = (du_dx_k + du_dx_kp1) / 2._dp
        w2 = (dv_dy_k + dv_dy_kp1) / 2._dp   
        
        w3 = ((dHs_dx_a( vi) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * dHi_dx_a( vi)) / MAX(0.1_dp, ice%Hi_a( vi))) *   &
             ((ice%u_3D_a( vi,k+1) - ice%u_3D_a( vi,k)) / (C%zeta(k+1) - C%zeta(k)))
        w4 = ((dHs_dy_a( vi) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * dHi_dy_a( vi)) / MAX(0.1_dp, ice%Hi_a( vi))) *   &
             ((ice%v_3D_a( vi,k+1) - ice%v_3D_a( vi,k)) / (C%zeta(k+1) - C%zeta(k)))

        ice%w_3D_a( vi,k) = ice%w_3D_a( vi,k+1) - ice%Hi_a( vi) * (w1 + w2 + w3 + w4) * (C%zeta(k+1) - C%zeta(k))
        
      END DO ! DO k = C%nZ - 1, 1, -1

    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdHi_dx_a  )
    CALL deallocate_shared( wdHi_dy_a  )
    CALL deallocate_shared( wdHs_dx_a  )
    CALL deallocate_shared( wdHs_dy_a  )
    CALL deallocate_shared( wdu_dx_3D_a)
    CALL deallocate_shared( wdv_dy_3D_a)
    
  END SUBROUTINE calc_3D_vertical_velocities
  
  ! Calculate some physical terms (basal yield stress, effective viscosity, etc.)
  SUBROUTINE calculate_driving_stress( mesh, ice)
    ! Calculate the driving stress taud (in both x and y) on the combi-mesh uv-field
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER                ::  dHs_dx_ac,  dHs_dy_ac
    INTEGER                                            :: wdHs_dx_ac, wdHs_dy_ac
    INTEGER                                            :: avi
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nVAaAc, dHs_dx_ac, wdHs_dx_ac)
    CALL allocate_shared_dp_1D( mesh%nVAaAc, dHs_dy_ac, wdHs_dy_ac)
    
    ! Map ice thickness to the aca-grid
    CALL map_a_to_ac_2D( mesh, ice%Hi_a, ice%Hi_ac)
    
    ! Calculate surface slopes on the aca-grid
    CALL ddx_a_to_ac_2D( mesh, ice%Hs_a, dHs_dx_ac)
    CALL ddy_a_to_ac_2D( mesh, ice%Hs_a, dHs_dy_ac)
    
    ! Calculate driving stresses on the aca-grid
    DO avi =  mesh%avi1, mesh%avi2
      ice%taudx_ac( avi) = -ice_density * grav * ice%Hi_ac( avi) * dHs_dx_ac( avi)
      ice%taudy_ac( avi) = -ice_density * grav * ice%Hi_ac( avi) * dHs_dy_ac( avi)
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx_ac)
    CALL deallocate_shared( wdHs_dy_ac)
    
  END SUBROUTINE calculate_driving_stress
  SUBROUTINE calc_effective_viscosity_ac( mesh, ice, u_ac, v_ac)
    ! Calculate 3D effective viscosity following Lipscomb et al. (2019), Eq. 2

    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_ac
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_ac
    
    ! Local variables
    INTEGER                                            :: avi,k
    REAL(dp)                                           :: eps_sq
    REAL(dp), PARAMETER                                :: epsilon_sq_0 = 1E-15_dp   ! Normalisation term so that zero velocity gives non-zero viscosity
    REAL(dp), DIMENSION(C%nz)                          :: prof
    
    ! Map ice thickness and ice flow factor to the aca-grid
    CALL map_a_to_ac_2D( mesh, ice%Hi_a       , ice%Hi_ac       )
    CALL map_a_to_ac_3D( mesh, ice%A_flow_3D_a, ice%A_flow_3D_ac)

    ! Calculate effective strain components from horizontal stretching on the regular mesh
    CALL ddx_ac_to_ac_2D( mesh, u_ac, ice%du_dx_ac)
    CALL ddy_ac_to_ac_2D( mesh, u_ac, ice%du_dy_ac)
    CALL ddx_ac_to_ac_2D( mesh, v_ac, ice%dv_dx_ac)
    CALL ddy_ac_to_ac_2D( mesh, v_ac, ice%dv_dy_ac)

    DO avi = mesh%avi1, mesh%avi2

      DO k = 1, C%nz
    
        ! Calculate the total effective strain rate from L19, Eq. 21 
        eps_sq = ice%du_dx_ac( avi)**2 + &
                 ice%dv_dy_ac( avi)**2 + &
                 ice%du_dx_ac( avi) * ice%dv_dy_ac( avi) + &
                 0.25_dp * (ice%du_dy_ac( avi) + ice%dv_dx_ac( avi))**2 + &  ! + 0.25_dp * (du_dz_a**2 + dv_dz_a**2) + &
                 epsilon_sq_0
        
        ! Calculate effective viscosity on ab-nodes
        ice%visc_eff_3D_ac( avi,k) = 0.5_dp * ice%A_flow_3D_ac( avi,k)**(-1._dp/C%n_flow) * (eps_sq)**((1._dp - C%n_flow)/(2._dp*C%n_flow))

      END DO ! DO k = 1, C%nz
      
      ! Vertical integral
      prof = ice%visc_eff_3D_ac( avi,:)
      CALL vertical_integrate( prof, ice%visc_eff_int_ac( avi))
      
      ! Product term N = eta * H
      ice%N_ac( avi) = ice%visc_eff_int_ac( avi) * MAX( 0.1_dp, ice%Hi_ac( avi))

    END DO  
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%visc_eff_3D_ac,  'ice%visc_eff_3D_ac' , 'calc_effective_viscosity')
    CALL check_for_NaN_dp_1D( ice%visc_eff_int_ac, 'ice%visc_eff_int_ac', 'calc_effective_viscosity')
    CALL check_for_NaN_dp_1D( ice%N_ac,            'ice%N_ac'           , 'calc_effective_viscosity')

  END SUBROUTINE calc_effective_viscosity_ac
  SUBROUTINE calc_sliding_term_beta( mesh, ice, u_ac, v_ac)
    ! Calculate the sliding term beta in the SSA/DIVA using the specified sliding law
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_ac
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_ac
    
    ! Local variables:
    INTEGER                                            :: avi
    
    ! Calculate the basal friction coefficients beta on the a-grid
    CALL calc_sliding_law( mesh, ice, u_ac, v_ac, ice%beta_ac)
    
    ! Limit beta to improve stability
    DO avi = mesh%avi1, mesh%avi2
      ice%beta_ac( avi) = MIN( C%DIVA_beta_max, ice%beta_ac( avi))
    END DO
    CALL sync
    
!    ! LEGACY - Apply the flotation mask (this is how we did it before we introduced the PISM/CISM-style sub-grid grounded fraction)
!    IF (.NOT. C%do_GL_subgrid_friction) THEN
!      DO i = grid%i1, grid%i2
!      DO j = 1, grid%ny
!        IF (ice%mask_ocean_a( j,i) == 1) ice%beta_a( j,i) = 0._dp
!      END DO
!      END DO
!      CALL sync
!    END IF
    
    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_ac, 'ice%beta_ac', 'calc_sliding_term_beta')
    
  END SUBROUTINE calc_sliding_term_beta
  
  ! Routines for solving the "linearised" SSA
  SUBROUTINE solve_SSA_sans_linearised( mesh, ice, u_ac, v_ac)
    ! Solve the "linearised" version of the SSA (i.e. assuming viscosity and basal stress are
    ! constant rather than functions of velocity).
    ! 
    ! The full SSA reads:
    ! 
    !   d/dx[ 2N (2*du/dx + dv/dy)] + d/dy[ N (du/dy + dv/dx)] - beta*u = -taud_x
    !   d/dy[ 2N (2*dv/dy + du/dx)] + d/dx[ N (dv/dx + du/dy)] - beta*v = -taud_y
    !
    ! Using the chain rule, this expands to:
    !
    !   4*N*d2u/dx2 + 4*dN/dx*du/dx + 2*N*d2v/dxdy + 2*dN/dx*dv/dy + N*d2u/dy2 + dN/dy*du/dy + N*d2v/dxdy + dN/dy*dv/dx - beta*u = -taud_x
    !   4*N*d2v/dy2 + 4*dN/dy*dv/dy + 2*N*d2u/dxdy + 2*dN/dy*du/dx + N*d2v/dx2 + dN/dx*dv/dx + N*d2u/dxdy + dN/dx*du/dy - beta*v = -taud_y
    !
    ! By neglecting all the "cross-terms" involving gradients of N, this simplifies to:
    !
    !   4*N*d2u/dx2 + N*d2u/dy2 + 3*N*d2v/dxdy - beta*u = -taud_x
    !   4*N*d2v/dy2 + N*d2v/dx2 + 3*N*d2u/dxdy - beta*v = -taud_y
    ! 
    ! Lastly, the left and right sides of the equations can be divided by N to yield:
    ! 
    ! 4*d2u/dx2 + d2u/dy2 + 3*d2v/dxdy - beta*u/N = -taud_x/N
    ! 4*d2v/dy2 + d2v/dx2 + 3*d2u/dxdy - beta*v/N = -taud_y/N
    ! 
    ! By happy coincidence, this equation is much easier to solve, especially because (for unclear reasons)
    ! it doesn't require N to be defined on a staggered grid relative to u,v in order to achieve a stable solution.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: u_ac
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: v_ac

    ! Local variables
    INTEGER                                            :: avi, auvi
    REAL(dp), DIMENSION(:    ), POINTER                ::  beta_over_N_ac,  beta_over_N_acu,  beta_over_N_acv,  beta_over_N_acuv
    INTEGER                                            :: wbeta_over_N_ac, wbeta_over_N_acu, wbeta_over_N_acv, wbeta_over_N_acuv
    TYPE(type_sparse_matrix_CSR)                       :: mbeta_over_N_acuv
    REAL(dp), DIMENSION(:    ), POINTER                ::  taudx_over_N_ac,  taudy_over_N_ac,  taudx_over_N_acu,  taudy_over_N_acv,  b_acuv
    INTEGER                                            :: wtaudx_over_N_ac, wtaudy_over_N_ac, wtaudx_over_N_acu, wtaudy_over_N_acv, wb_acuv
    TYPE(type_sparse_matrix_CSR)                       :: A, BC
    REAL(dp), DIMENSION(:    ), POINTER                ::  u_acu,  v_acv,  uv_acuv
    INTEGER                                            :: wu_acu, wv_acv, wuv_acuv
    
  ! Construct the big matrix A
  ! ==========================
    
    ! Calculate beta / N on the acuv-grid
    CALL allocate_shared_dp_1D(   mesh%nVAaAc, beta_over_N_ac  , wbeta_over_N_ac  )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, beta_over_N_acu , wbeta_over_N_acu )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, beta_over_N_acv , wbeta_over_N_acv )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, beta_over_N_acuv, wbeta_over_N_acuv)
    
    DO avi = mesh%avi1, mesh%avi2
      beta_over_N_ac( avi) = ice%beta_ac( avi) / ice%N_ac( avi)
    END DO
    CALL sync
    
    CALL move_ac_to_acu_2D( mesh, beta_over_N_ac, beta_over_N_acu)
    CALL move_ac_to_acv_2D( mesh, beta_over_N_ac, beta_over_N_acv)
    CALL deallocate_shared( wbeta_over_N_ac)
    
    DO auvi = mesh%auvi1, mesh%auvi2
      beta_over_N_acuv( auvi) = beta_over_N_acu( auvi) + beta_over_N_acv( auvi)
    END DO
    CALL sync
    
    CALL deallocate_shared( wbeta_over_N_acu)
    CALL deallocate_shared( wbeta_over_N_acv)
    
    ! Convert sliding term beta from vector to diagonal matrix
    CALL convert_vector_to_diag_CSR( beta_over_N_acuv, mbeta_over_N_acuv)
    CALL deallocate_shared( wbeta_over_N_acuv)
    
    ! Subtract beta / N from the unchanging part of the matrix to get A
    CALL add_matrix_matrix_CSR( ice%M_SSA_sans_no_beta_ac, mbeta_over_N_acuv, A, alpha = 1._dp, beta = -1._dp)
    CALL deallocate_matrix_CSR( mbeta_over_N_acuv)
    
  ! Construct the right-hand side, initial guess, and boundary conditions
  ! =====================================================================
    
    ! Fill in the right-hand side b = -taudx,y / N
    
    CALL allocate_shared_dp_1D(   mesh%nVAaAc, taudx_over_N_ac , wtaudx_over_N_ac )
    CALL allocate_shared_dp_1D(   mesh%nVAaAc, taudy_over_N_ac , wtaudy_over_N_ac )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, taudx_over_N_acu, wtaudx_over_N_acu)
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, taudy_over_N_acv, wtaudy_over_N_acv)
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, b_acuv          , wb_acuv          )
    
    DO avi = mesh%avi1, mesh%avi2
      taudx_over_N_ac( avi) = ice%taudx_ac( avi) / ice%N_ac( avi)
      taudy_over_N_ac( avi) = ice%taudy_ac( avi) / ice%N_ac( avi)
    END DO
    CALL sync
    
    CALL move_ac_to_acu_2D( mesh, taudx_over_N_ac, taudx_over_N_acu)
    CALL move_ac_to_acv_2D( mesh, taudy_over_N_ac, taudy_over_N_acv)
    CALL deallocate_shared( wtaudx_over_N_ac)
    CALL deallocate_shared( wtaudy_over_N_ac)
    
    DO auvi = mesh%auvi1, mesh%auvi2
      b_acuv( auvi) = -1._dp * (taudx_over_N_acu( auvi) + taudy_over_N_acv( auvi))
    END DO
    CALL sync
    
    CALL deallocate_shared( wtaudx_over_N_acu)
    CALL deallocate_shared( wtaudy_over_N_acv)
    
    ! Fill in the initial guess with the current velocity solution
    
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, u_acu  , wu_acu  )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, v_acv  , wv_acv  )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, uv_acuv, wuv_acuv)
    
    CALL move_ac_to_acu_2D( mesh, u_ac, u_acu)
    CALL move_ac_to_acv_2D( mesh, v_ac, v_acv)
    
    DO auvi = mesh%auvi1, mesh%auvi2
      uv_acuv( auvi) = u_acu( auvi) + v_acv( auvi)
    END DO
    CALL sync
    
    CALL deallocate_shared( wu_acu)
    CALL deallocate_shared( wv_acv)
    
    ! Construct the boundary conditions matrix, and adapt the right-hand side and initial guess accordingly
    CALL construct_BC_SSA( mesh, BC, b_acuv, uv_acuv)
    
    ! Add the boundary conditions to A
    CALL overwrite_rows_CSR( A, BC)
    CALL deallocate_matrix_CSR( BC)
    
  ! Solve the matrix equation
  ! =========================
    
    CALL solve_matrix_equation_CSR( A, b_acuv, uv_acuv, &
      C%DIVA_choice_matrix_solver, &
      C%DIVA_SOR_nit             , &
      C%DIVA_SOR_tol             , &
      C%DIVA_SOR_omega           , &
      C%DIVA_PETSc_rtol          , &
      C%DIVA_PETSc_abstol        , &
      colour_v1 = mesh%colour_v1, &
      colour_v2 = mesh%colour_v2, &
      colour_vi = mesh%colour_vi)
    CALL deallocate_matrix_CSR( A)
    CALL deallocate_shared( wb_acuv )
    
  ! Get solution back on the ac-grid
  ! ================================
    
    CALL move_acu_to_ac_2D( mesh, uv_acuv, u_ac)
    CALL move_acv_to_ac_2D( mesh, uv_acuv, v_ac)
    CALL deallocate_shared( wuv_acuv)

  END SUBROUTINE solve_SSA_sans_linearised
  SUBROUTINE construct_BC_SSA( mesh, BC, b_acuv, uv_acuv)
    ! Construct the boundary conditions matrix, and adapt the right-hand side and initial guess accordingly
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: BC
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: b_acuv, uv_acuv
    
    ! Local variables:
    INTEGER                                            :: nnz_BC
    INTEGER                                            :: auvi, avi, vi, aci, edge_index
    
    ! Find maximum number of non-zero entries
    nnz_BC = 0
    DO auvi = mesh%auvi1, mesh%auvi2
    
      IF (MOD(auvi,2) == 1) THEN
        ! u
        avi = (auvi+1)/2
      ELSE
        ! v
        avi = auvi/2
      END IF
    
      IF (avi <= mesh%nV) THEN
        ! vertex
        vi = avi
        edge_index = mesh%edge_index( vi)
      ELSE
        ! edge
        aci = avi - mesh%nV
        edge_index = mesh%edge_index_ac( aci)
      END IF
      
      IF (edge_index > 0) THEN
        ! border vertex/edge
        nnz_BC = nnz_BC + 1 ! + mesh%nCAaAc( avi)
      END IF
      
    END DO ! DO auvi = mesh%auvi1, mesh%auvi2
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nnz_BC, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)   
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( BC, 2*mesh%nVAaAc, 2*mesh%nVAaAc, nnz_BC)
    
    ! Fill in values
    BC%nnz = 0
    BC%ptr = 1
    DO auvi = mesh%auvi1, mesh%auvi2
    
      IF (MOD(auvi,2) == 1) THEN
        ! u
        avi = (auvi+1)/2
      ELSE
        ! v
        avi = auvi/2
      END IF
    
      IF (avi <= mesh%nV) THEN
        ! vertex
        vi = avi
        edge_index = mesh%edge_index( vi)
      ELSE
        ! edge
        aci = avi - mesh%nV
        edge_index = mesh%edge_index_ac( aci)
      END IF
    
      IF     (edge_index == 0) THEN
        ! Free vertex/edge: no boundary conditions apply
        
      ELSE
        ! Border vertex/edge: set velocity to zero
        
        ! Matrix
        BC%nnz  = BC%nnz+1
        BC%index( BC%nnz) = auvi
        BC%val(   BC%nnz) = 1._dp
        
        ! Right-hand side vector
        b_acuv(  auvi) = 0._dp
        uv_acuv( auvi) = 0._dp
        
      END IF
      
      ! Finalise matrix row
      BC%ptr( auvi+1 : 2*mesh%nVAaAc+1) = BC%nnz + 1
      
    END DO ! DO auvi = mesh%auvi1, mesh%auvi2
    CALL sync
    
    ! Finalise matrix
    CALL finalise_matrix_CSR_dist( BC, mesh%auvi1, mesh%auvi2)
    
  END SUBROUTINE construct_BC_SSA
  SUBROUTINE initialise_matrix_SSA_sans( mesh, ice)
    ! Initialise the unchanging part (i.e. without the beta/N term) of
    ! the big matrix equation representing the SSA without the "cross-terms"
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR)                       :: M1
    TYPE(type_sparse_matrix_CSR)                       :: M_d2dx2_acu_acu
    TYPE(type_sparse_matrix_CSR)                       :: M_d2dy2_acu_acu
    TYPE(type_sparse_matrix_CSR)                       :: M_d2dxdy_acv_acu
    TYPE(type_sparse_matrix_CSR)                       :: M_d2dy2_acv_acv
    TYPE(type_sparse_matrix_CSR)                       :: M_d2dx2_acv_acv
    TYPE(type_sparse_matrix_CSR)                       :: M_d2dxdy_acu_acv
    TYPE(type_sparse_matrix_CSR)                       :: Mu
    TYPE(type_sparse_matrix_CSR)                       :: Mv
    
  ! Equation 1
  ! ==========
    
    ! M_d2dx2_acu_acu = M_move_ac_acu * M_d2dx2_ac_ac * M_move_acu_ac
    CALL multiply_matrix_matrix_CSR( mesh%M_d2dx2_ac_ac , mesh%M_move_acu_ac, M1              )
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acu , M1                , M_d2dx2_acu_acu )
    CALL deallocate_matrix_CSR( M1)
    
    ! M_d2dy2_acu_acu = M_move_ac_acu * M_d2dy2_ac_ac * M_move_acu_ac
    CALL multiply_matrix_matrix_CSR( mesh%M_d2dy2_ac_ac , mesh%M_move_acu_ac, M1              )
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acu , M1                , M_d2dy2_acu_acu )
    CALL deallocate_matrix_CSR( M1)
    
    ! M_d2dxdy_acv_acu = M_move_ac_acu * M_d2dxdy_ac_ac * M_move_acv_ac
    CALL multiply_matrix_matrix_CSR( mesh%M_d2dxdy_ac_ac, mesh%M_move_acv_ac, M1              )
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acu , M1                , M_d2dxdy_acv_acu)
    CALL deallocate_matrix_CSR( M1)
    
    ! M1 = 4*M_d2dx2_acu_acu + M_d2dy2_acu_acu
    CALL add_matrix_matrix_CSR( M_d2dx2_acu_acu, M_d2dy2_acu_acu , M1, alpha = 4._dp, beta = 1._dp)
    
    ! Mu = M1 + 3*M_d2dxdy_acv_acu
    CALL add_matrix_matrix_CSR( M1             , M_d2dxdy_acv_acu, Mu, alpha = 1._dp, beta = 3._dp)
    CALL deallocate_matrix_CSR( M_d2dx2_acu_acu )
    CALL deallocate_matrix_CSR( M_d2dy2_acu_acu )
    CALL deallocate_matrix_CSR( M_d2dxdy_acv_acu)
    CALL deallocate_matrix_CSR( M1)
    
  ! Equation 2
  ! ==========
    
    ! M_d2dy2_acv_acv = M_move_ac_acv * M_d2dy2_ac_ac * M_move_acv_ac
    CALL multiply_matrix_matrix_CSR( mesh%M_d2dy2_ac_ac , mesh%M_move_acv_ac, M1              )
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acv , M1                , M_d2dy2_acv_acv )
    CALL deallocate_matrix_CSR( M1)
    
    ! M_d2dx2_acv_acv = M_move_ac_acv * M_d2dx2_ac_ac * M_move_acv_ac
    CALL multiply_matrix_matrix_CSR( mesh%M_d2dx2_ac_ac , mesh%M_move_acv_ac, M1              )
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acv , M1                , M_d2dx2_acv_acv )
    CALL deallocate_matrix_CSR( M1)
    
    ! M_d2dxdy_acu_acv = M_move_ac_acv * M_d2dxdy_ac_ac * M_move_acu_ac
    CALL multiply_matrix_matrix_CSR( mesh%M_d2dxdy_ac_ac, mesh%M_move_acu_ac, M1              )
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acv , M1                , M_d2dxdy_acu_acv)
    CALL deallocate_matrix_CSR( M1)
    
    ! M1 = 4*M_d2dy2_acv_acv + M_d2dx2_acv_acv
    CALL add_matrix_matrix_CSR( M_d2dy2_acv_acv, M_d2dx2_acv_acv , M1, alpha = 4._dp, beta = 1._dp)
    
    ! Mv = M1 + 3*M_d2dxdy_acu_acv
    CALL add_matrix_matrix_CSR( M1             , M_d2dxdy_acu_acv, Mv, alpha = 1._dp, beta = 3._dp)
    CALL deallocate_matrix_CSR( M_d2dy2_acv_acv )
    CALL deallocate_matrix_CSR( M_d2dx2_acv_acv )
    CALL deallocate_matrix_CSR( M_d2dxdy_acu_acv)
    CALL deallocate_matrix_CSR( M1)
    
  ! Add them together
  ! =================
  
    CALL add_matrix_matrix_CSR( Mu, Mv, ice%M_SSA_sans_no_beta_ac)
    CALL deallocate_matrix_CSR( Mu)
    CALL deallocate_matrix_CSR( Mv)
    
  END SUBROUTINE initialise_matrix_SSA_sans
    
  ! Some "administration" routines that help speed up and stabilise the solver
  SUBROUTINE apply_velocity_limits( mesh, u_ac, v_ac)
    ! Apply a velocity limit (for stability)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: u_ac, v_ac

    ! Local variables:
    INTEGER                                            :: avi
    REAL(dp), DIMENSION(:    ), POINTER                ::  uabs_ac
    INTEGER                                            :: wuabs_ac
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nVAaAc, uabs_ac, wuabs_ac)
    
    ! Calculate absolute speed
    DO avi = mesh%avi1, mesh%avi2
      uabs_ac( avi) = SQRT( u_ac( avi)**2 + v_ac( avi)**2)
    END DO
    CALL sync
    
    ! Scale velocities accordingly
    DO avi = mesh%avi1, mesh%avi2
      IF (uabs_ac( avi) > C%DIVA_vel_max) THEN
        u_ac( avi) = u_ac( avi) * C%DIVA_vel_max / uabs_ac( avi)
        v_ac( avi) = v_ac( avi) * C%DIVA_vel_max / uabs_ac( avi)
      END IF
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wuabs_ac)

  END SUBROUTINE apply_velocity_limits
  SUBROUTINE relax_DIVA_visc_iterations( mesh, u_ac_prev, v_ac_prev, u_ac, v_ac, rel)
    ! Relax velocity solution with results of the previous viscosity iteration 
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_ac_prev
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_ac_prev
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: u_ac
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: v_ac
    REAL(dp),                            INTENT(IN)    :: rel
        
    ! Local variables:
    INTEGER                                            :: avi

    DO avi = mesh%avi1, mesh%avi2
      u_ac( avi) = rel * u_ac( avi) + (1._dp - rel) * u_ac_prev( avi)
      v_ac( avi) = rel * v_ac( avi) + (1._dp - rel) * v_ac_prev( avi)
    END DO
    CALL sync

  END SUBROUTINE relax_DIVA_visc_iterations
  SUBROUTINE calc_visc_iter_UV_resid( mesh, u_ac_prev, v_ac_prev, u_ac, v_ac, resid_UV)
    ! Check if the viscosity iteration has converged to a stable solution
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_ac_prev
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_ac_prev
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_ac
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_ac
    REAL(dp),                            INTENT(OUT)   :: resid_UV
    
    ! Local variables:
    INTEGER                                            :: ierr
    INTEGER                                            :: avi,navi
    REAL(dp)                                           :: res1, res2
    REAL(dp), PARAMETER                                :: DIVA_vel_tolerance = 1e-6   ! [m/a] only consider points with velocity above this tolerance limit

    ! Calculate the L2 norm based on velocity solution between previous
    ! and current viscosity iteration (as in Yelmo/SICOPOLIS)
    
    navi = 0
    res1 = 0._dp
    res2 = 0._dp
    
    DO avi = mesh%avi1, mesh%avi2
      
      IF (ABS(u_ac( avi)) > DIVA_vel_tolerance) THEN
        navi = navi + 1
        res1 = res1 + (u_ac( avi) - u_ac_prev( avi))**2._dp
        res2 = res2 + (u_ac( avi) + u_ac_prev( avi))**2._dp
      END IF
    
      IF (ABS(v_ac( avi)) > DIVA_vel_tolerance) THEN
        navi = navi + 1
        res1 = res1 + (v_ac( avi) - v_ac_prev( avi))**2._dp
        res2 = res2 + (v_ac( avi) + v_ac_prev( avi))**2._dp
      END IF
      
    END DO
    
    ! Combine results from all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, navi, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    IF (navi > 0) THEN
      res1 = SQRT( res1)
      res2 = SQRT( res2 )
      res2 = MAX(  res2,1E-8_dp)
      resid_UV = 2._dp * res1 / res2 
    ELSE 
      ! No points available for comparison, set residual equal to zero 
      resid_UV = 0._dp
    END IF

  END SUBROUTINE calc_visc_iter_UV_resid
  
  ! Map velocity components to the a-mesh for writing to output (diagnostic only)
  SUBROUTINE map_velocity_from_c_to_a_2D( mesh, ice, u_c, u_a)
    ! Map velocity components from the staggered c-mesh to the regular a-mesh.
    ! Take the average over those staggered vertices where the opposite regular
    ! vertex has the same ice mask as the current regular vertex.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: u_a
    
    ! Local variables:
    INTEGER                                            :: vi, ci, aci, vj
    REAL(dp)                                           :: uav, w
    
    DO vi = mesh%vi1, mesh%vi2
      
      uav = 0._dp
      w   = 0._dp
      
      DO ci = 1, mesh%nC( vi)
      
        aci = mesh%iAci( vi,ci)
        vj  = mesh%C(    vi,ci)
        
        IF (ice%mask_ice_a( vj) == ice%mask_ice_a( vi)) THEN
          uav = uav + u_c( aci)
          w   = w   + 1._dp
        END IF
        
      END DO
      
      u_a( vi) = uav / w
      
    END DO
    CALL sync
    
  END SUBROUTINE map_velocity_from_c_to_a_2D
  SUBROUTINE map_velocity_from_c_to_a_3D( mesh, ice, u_c, u_a)
    ! Map velocity components from the staggered c-mesh to the regular a-mesh.
    ! Take the average over those staggered vertices where the opposite regular
    ! vertex has the same ice mask as the current regular vertex.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: u_a
    
    ! Local variables:
    INTEGER                                            :: k
    REAL(dp), DIMENSION(:    ), POINTER                ::  u_c_k,  u_a_k
    INTEGER                                            :: wu_c_k, wu_a_k
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nAc, u_c_k, wu_c_k)
    CALL allocate_shared_dp_1D( mesh%nV , u_a_k, wu_a_k)
    
    ! Perform the mapping operations layer-by-layer
    DO k = 1, C%nz
    
      u_c_k( mesh%ci1:mesh%ci2  ) = u_c(   mesh%ci1:mesh%ci2,k)
      CALL sync
      
      CALL map_velocity_from_c_to_a_2D( mesh, ice, u_c_k, u_a_k)
      
      u_a(   mesh%vi1:mesh%vi2,k) = u_a_k( mesh%vi1:mesh%vi2  )
      CALL sync
      
    END DO
    
    ! Clean up after yourself
    CALL deallocate_shared( wu_c_k)
    CALL deallocate_shared( wu_a_k)
    
  END SUBROUTINE map_velocity_from_c_to_a_3D
  
END MODULE ice_velocity_module
