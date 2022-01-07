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
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_sparse_matrix_CSR, type_remapping
  USE mesh_operators_module,           ONLY: map_a_to_c_2D, ddx_a_to_c_2D, ddy_a_to_c_2D, map_a_to_c_3D, &
                                             map_a_to_b_2D, ddx_a_to_b_2D, ddy_a_to_b_2D, map_b_to_c_2D, &
                                             ddx_b_to_a_2D, ddy_b_to_a_2D, map_b_to_a_3D, map_b_to_a_2D, &
                                             ddx_a_to_a_2D, ddy_a_to_a_2D, ddx_c_to_a_3D, ddy_c_to_a_3D, &
                                             map_c_to_a_2D, map_c_to_a_3D, map_a_to_b_3D, map_b_to_c_3D
  USE utilities_module,                ONLY: vertical_integration_from_bottom_to_zeta, vertical_average, &
                                             vertical_integrate, allocate_matrix_CSR_dist, sort_columns_in_CSR_dist, &
                                             finalise_matrix_CSR_dist, solve_matrix_equation_CSR, deallocate_matrix_CSR
  USE basal_conditions_and_sliding_module, ONLY: calc_basal_conditions, calc_sliding_law
  USE netcdf_module,                   ONLY: write_CSR_matrix_to_NetCDF
  USE utilities_module,                ONLY: SSA_Schoof2006_analytical_solution
  USE general_ice_model_data_module,   ONLY: determine_grounded_fractions

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
    
    ! Local variables:
    INTEGER                                            :: ti
    LOGICAL                                            :: set_velocities_to_zero
    LOGICAL                                            :: has_converged
    INTEGER                                            :: viscosity_iteration_i
    REAL(dp)                                           :: resid_UV
    REAL(dp)                                           :: umax_analytical, tauc_analytical
    
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
      ice%u_base_SSA_b( mesh%ti1:mesh%ti2) = 0._dp
      ice%v_base_SSA_b( mesh%ti1:mesh%ti2) = 0._dp
      CALL sync
      CALL calc_secondary_velocities( mesh, ice)
      RETURN
    END IF
    
  ! == Safety - end
  ! ===============
    
    ! Calculate the driving stresses taudx, taudy
    CALL calculate_driving_stress( mesh, ice)

    ! Calculate the basal yield stress tau_c
    CALL calc_basal_conditions( mesh, ice)
    
    ! Determine sub-mesh grounded fractions for scaling the basal friction
    CALL determine_grounded_fractions( mesh, ice)
    
    ! Find analytical solution for the SSA icestream experiment (used only to print numerical error to screen)
    CALL SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, ice%A_flow_vav_a( 1), 0._dp, umax_analytical, tauc_analytical)
    
    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .FALSE.
    viscosity_iteration: DO WHILE (.NOT. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1
      
      ! Calculate the effective viscosity and the product term N = eta * H
      CALL calc_effective_viscosity( mesh, ice, ice%u_base_SSA_b, ice%v_base_SSA_b)
            
      ! Calculate the sliding term beta
      CALL calc_sliding_term_beta( mesh, ice, ice%u_base_SSA_b, ice%v_base_SSA_b)
      
      ! Set beta_eff equal to beta; this turns the DIVA into the SSA
      ice%beta_eff_a( mesh%vi1:mesh%vi2) = ice%beta_a( mesh%vi1:mesh%vi2)
      CALL sync
    
      ! Map beta_eff from the a-grid to the cx/cy-grids
      CALL map_a_to_b_2D( mesh, ice%beta_eff_a, ice%beta_eff_b)
    
      ! Apply the sub-grid grounded fraction
      DO ti = mesh%ti1, mesh%ti2
        ice%beta_eff_b( ti) = ice%beta_eff_b( ti) * ice%f_grnd_b( ti)**2
      END DO
      CALL sync
      
      ! Store the previous solution so we can check for convergence later
      ice%u_prev_b( mesh%ti1:mesh%ti2) = ice%u_base_SSA_b( mesh%ti1:mesh%ti2)
      ice%v_prev_b( mesh%ti1:mesh%ti2) = ice%v_base_SSA_b( mesh%ti1:mesh%ti2)
      CALL sync
      
      ! Solve the linearised SSA
      CALL solve_SSADIVA_linearised( mesh, ice, ice%u_base_SSA_b, ice%v_base_SSA_b)
      
      ! Apply velocity limits (both overflow and underflow) for improved stability
      CALL apply_velocity_limits( mesh, ice%u_base_SSA_b, ice%v_base_SSA_b)
      
      ! "relax" subsequent viscosity iterations for improved stability
      CALL relax_DIVA_visc_iterations( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_base_SSA_b, ice%v_base_SSA_b, C%DIVA_visc_it_relax)
      
      ! Check if the viscosity iteration has converged
      CALL calc_visc_iter_UV_resid( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_base_SSA_b, ice%v_base_SSA_b, resid_UV)
      IF (par%master) WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ': resid_UV = ', resid_UV, ', u = [', MINVAL(ice%u_base_SSA_b), ' - ', MAXVAL(ice%u_base_SSA_b), ']'
      
      IF (par%master .AND. C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'SSA_icestream') &
        WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ': err = ', ABS(1._dp - MAXVAL(ice%u_base_SSA_b) / umax_analytical), ': resid_UV = ', resid_UV

      has_converged = .FALSE.
      IF     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) THEN
        has_converged = .TRUE.
      ELSEIF (viscosity_iteration_i >= C%DIVA_visc_it_nit) THEN
        has_converged = .TRUE.
      END IF
      
    END DO viscosity_iteration
    
    ! Calculate secondary velocities (surface, base, etc.)
    CALL calc_secondary_velocities( mesh, ice)
    
  END SUBROUTINE solve_SSA
  SUBROUTINE solve_DIVA(  mesh, ice)
    ! Calculate ice velocities using the DIVA
      
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
    
  ! =========
  ! == Safety
    
    ! Check that this routine is called correctly
    IF (.NOT. C%choice_ice_dynamics == 'DIVA') THEN
      IF (par%master) WRITE(0,*) 'ERROR - solve_DIVA should only be called when choice_ice_dynamics is set to DIVA!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! If there's no grounded ice anywhere, don't bother
    set_velocities_to_zero = .FALSE.
    IF (SUM( ice%mask_sheet_a) == 0) set_velocities_to_zero = .TRUE.
    
    IF (set_velocities_to_zero) THEN
      ice%u_vav_DIVA_b( mesh%ti1:mesh%ti2) = 0._dp
      ice%v_vav_DIVA_b( mesh%ti1:mesh%ti2) = 0._dp
      CALL sync
      CALL calc_secondary_velocities( mesh, ice)
      RETURN
    END IF
    
  ! == Safety - end
  ! ===============
    
    ! Calculate the driving stresses taudx, taudy
    CALL calculate_driving_stress( mesh, ice)

    ! Calculate the basal yield stress tau_c
    CALL calc_basal_conditions( mesh, ice)
    
    ! Determine sub-mesh grounded fractions for scaling the basal friction
    CALL determine_grounded_fractions( mesh, ice)
    
    ! Find analytical solution for the SSA icestream experiment (used only to print numerical error to screen)
    CALL SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, ice%A_flow_vav_a( 1), 0._dp, umax_analytical, tauc_analytical)
    
    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .FALSE.
    viscosity_iteration: DO WHILE (.NOT. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1
      
      ! Calculate the vertical shear strain rates
      CALL calc_vertical_shear_strain_rates( mesh, ice)
      
      ! Calculate the effective viscosity and the product term N = eta * H
      CALL calc_effective_viscosity( mesh, ice, ice%u_vav_DIVA_b, ice%v_vav_DIVA_b)
            
      ! Calculate the sliding term beta
      CALL calc_sliding_term_beta( mesh, ice, ice%u_vav_DIVA_b, ice%v_vav_DIVA_b)
      
      ! Calculate the F-integral F2
      CALL calc_F_integral( mesh, ice, n = 2._dp)
      
      ! Calculate beta_eff
      CALL calc_beta_eff( mesh, ice)
      
      ! Store the previous solution so we can check for convergence later
      ice%u_prev_b( mesh%ti1:mesh%ti2) = ice%u_vav_DIVA_b( mesh%ti1:mesh%ti2)
      ice%v_prev_b( mesh%ti1:mesh%ti2) = ice%v_vav_DIVA_b( mesh%ti1:mesh%ti2)
      CALL sync
      
      ! Solve the linearised DIVA
      CALL solve_SSADIVA_linearised( mesh, ice, ice%u_vav_DIVA_b, ice%v_vav_DIVA_b)
      
      ! Apply velocity limits (both overflow and underflow) for improved stability
      CALL apply_velocity_limits( mesh, ice%u_vav_DIVA_b, ice%v_vav_DIVA_b)
      
      ! "relax" subsequent viscosity iterations for improved stability
      CALL relax_DIVA_visc_iterations( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_vav_DIVA_b, ice%v_vav_DIVA_b, C%DIVA_visc_it_relax)
      
      ! Check if the viscosity iteration has converged
      CALL calc_visc_iter_UV_resid( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_vav_DIVA_b, ice%v_vav_DIVA_b, resid_UV)
      IF (par%master) WRITE(0,*) '   DIVA - viscosity iteration ', viscosity_iteration_i, ': resid_UV = ', resid_UV, ', u = [', MINVAL(ice%u_vav_DIVA_b), ' - ', MAXVAL(ice%u_vav_DIVA_b), ']'
      
      IF (par%master .AND. C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'SSA_icestream') &
        WRITE(0,*) '   DIVA - viscosity iteration ', viscosity_iteration_i, ': err = ', ABS(1._dp - MAXVAL(ice%u_vav_DIVA_b) / umax_analytical), ': resid_UV = ', resid_UV

      has_converged = .FALSE.
      IF     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) THEN
        has_converged = .TRUE.
      ELSEIF (viscosity_iteration_i >= C%DIVA_visc_it_nit) THEN
        has_converged = .TRUE.
      END IF

      ! Calculate basal stress 
      CALL calc_basal_stress_DIVA( mesh, ice, ice%u_vav_DIVA_b, ice%v_vav_DIVA_b)
      
    END DO viscosity_iteration
    
    ! Calculate secondary velocities (surface, base, etc.)
    CALL calc_secondary_velocities( mesh, ice)
    
  END SUBROUTINE solve_DIVA
  
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
      
      ! Copy surface velocity from the 3D fields
      ice%u_surf_c( mesh%ci1:mesh%ci2) = ice%u_3D_c( mesh%ci1:mesh%ci2,1)
      ice%v_surf_c( mesh%ci1:mesh%ci2) = ice%v_3D_c( mesh%ci1:mesh%ci2,1)
      CALL sync
      
      ! Calculate vertically averaged velocities
      DO aci = mesh%ci1, mesh%ci2
        prof = ice%u_3D_c( aci,:)
        CALL vertical_average( prof, ice%u_vav_c( aci))
        prof = ice%v_3D_c( aci,:)
        CALL vertical_average( prof, ice%v_vav_c( aci))
      END DO
      CALL sync
      
      ! Map velocity components to the a-grid for writing to output
      CALL map_c_to_a_3D( mesh, ice%u_3D_c,   ice%u_3D_a  )
      CALL map_c_to_a_3D( mesh, ice%v_3D_c,   ice%v_3D_a  )
      CALL map_c_to_a_2D( mesh, ice%u_vav_c,  ice%u_vav_a )
      CALL map_c_to_a_2D( mesh, ice%v_vav_c,  ice%v_vav_a )
      CALL map_c_to_a_2D( mesh, ice%u_surf_c, ice%u_surf_a)
      CALL map_c_to_a_2D( mesh, ice%v_surf_c, ice%v_surf_a)
      
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      ! No SIA, just the SSA
      
    ! == From the b-grid to the c-grid (used in the model)
      
      ! Set basal velocity equal to SSA answer
      CALL map_b_to_c_2D( mesh, ice%u_base_SSA_b, ice%u_base_c)
      CALL map_b_to_c_2D( mesh, ice%v_base_SSA_b, ice%v_base_c)
      
      ! No vertical variations in velocity
      DO aci = mesh%ci1, mesh%ci2
        ice%u_vav_c(  aci   ) = ice%u_base_c( aci)
        ice%v_vav_c(  aci   ) = ice%v_base_c( aci)
        ice%u_surf_c( aci   ) = ice%u_base_c( aci)
        ice%v_surf_c( aci   ) = ice%v_base_c( aci)
        ice%u_3D_c(   aci ,:) = ice%u_base_c( aci)
        ice%v_3D_c(   aci ,:) = ice%v_base_c( aci)
      END DO
      CALL sync
      
    ! == From the b-grid to the a-grid (only for writing to output)
      
      ! Set basal velocity equal to SSA answer
      CALL map_b_to_a_2D( mesh, ice%u_base_SSA_b, ice%u_base_a)
      CALL map_b_to_a_2D( mesh, ice%v_base_SSA_b, ice%v_base_a)
      
      ! No vertical variations in velocity
      DO vi = mesh%vi1, mesh%vi2
        ice%u_vav_a(  vi   ) = ice%u_base_a( vi)
        ice%v_vav_a(  vi   ) = ice%v_base_a( vi)
        ice%u_surf_a( vi   ) = ice%u_base_a( vi)
        ice%v_surf_a( vi   ) = ice%v_base_a( vi)
        ice%u_3D_a(   vi ,:) = ice%u_base_a( vi)
        ice%v_3D_a(   vi ,:) = ice%v_base_a( vi)
      END DO
      CALL sync
      
    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      ! Hybrid SIA/SSA
      
    ! == From the b-grid to the c-grid (used in the model)
      
      ! Set basal velocity equal to SSA answer
      CALL map_b_to_c_2D( mesh, ice%u_base_SSA_b, ice%u_base_c)
      CALL map_b_to_c_2D( mesh, ice%v_base_SSA_b, ice%v_base_c)
      
      ! Set 3-D velocities equal to the SIA solution
      ice%u_3D_c( mesh%ci1:mesh%ci2,:) = ice%u_3D_SIA_c( mesh%ci1:mesh%ci2,:)
      ice%v_3D_c( mesh%ci1:mesh%ci2,:) = ice%v_3D_SIA_c( mesh%ci1:mesh%ci2,:)
      CALL sync
      
      ! Add the SSA contributions to the 3-D velocities
      DO k = 1, C%nz
        ice%u_3D_c( mesh%ci1:mesh%ci2,k) = ice%u_3D_c( mesh%ci1:mesh%ci2,k) + ice%u_base_c( mesh%ci1:mesh%ci2)
        ice%v_3D_c( mesh%ci1:mesh%ci2,k) = ice%v_3D_c( mesh%ci1:mesh%ci2,k) + ice%v_base_c( mesh%ci1:mesh%ci2)
      END DO
      CALL sync
      
      ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
      CALL calc_3D_vertical_velocities( mesh, ice)
      
      ! Copy surface velocity from the 3D fields
      ice%u_surf_c( mesh%ci1:mesh%ci2) = ice%u_3D_c( mesh%ci1:mesh%ci2,1)
      ice%v_surf_c( mesh%ci1:mesh%ci2) = ice%v_3D_c( mesh%ci1:mesh%ci2,1)
      CALL sync
      
      ! Calculate vertically averaged velocities
      DO aci = mesh%ci1, mesh%ci2
        prof = ice%u_3D_c( aci,:)
        CALL vertical_average( prof, ice%u_vav_c( aci))
        prof = ice%v_3D_c( aci,:)
        CALL vertical_average( prof, ice%v_vav_c( aci))
      END DO
      CALL sync
      
    ! == From the b-grid to the a-grid (only for writing to output)
      
      ! Set basal velocity equal to SSA answer
      CALL map_b_to_a_2D( mesh, ice%u_base_SSA_b, ice%u_base_a)
      CALL map_b_to_a_2D( mesh, ice%v_base_SSA_b, ice%v_base_a)
      
      ! Set 3-D velocities equal to the SIA solution
      CALL map_c_to_a_3D( mesh, ice%u_3D_SIA_c, ice%u_3D_a)
      CALL map_c_to_a_3D( mesh, ice%v_3D_SIA_c, ice%v_3D_a)
      
      ! Add the SSA contributions to the 3-D velocities
      DO k = 1, C%nz
        ice%u_3D_a( mesh%vi1:mesh%vi2,k) = ice%u_3D_a( mesh%vi1:mesh%vi2,k) + ice%u_base_a( mesh%vi1:mesh%vi2)
        ice%v_3D_a( mesh%vi1:mesh%vi2,k) = ice%v_3D_a( mesh%vi1:mesh%vi2,k) + ice%v_base_a( mesh%vi1:mesh%vi2)
      END DO
      CALL sync
      
      ! Copy surface velocity from the 3D fields
      ice%u_surf_a( mesh%vi1:mesh%vi2) = ice%u_3D_a( mesh%vi1:mesh%vi2,1)
      ice%v_surf_a( mesh%vi1:mesh%vi2) = ice%v_3D_a( mesh%vi1:mesh%vi2,1)
      CALL sync
      
      ! Calculate vertically averaged velocities
      DO vi = mesh%vi1, mesh%vi2
        prof = ice%u_3D_a( vi,:)
        CALL vertical_average( prof, ice%u_vav_a( vi))
        prof = ice%v_3D_a( vi,:)
        CALL vertical_average( prof, ice%v_vav_a( vi))
      END DO
      CALL sync
    
    ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      ! DIVA
      
      ! Calculate basal velocity from depth-averaged solution and basal stress on the b-grid
      CALL calc_basal_velocities_DIVA( mesh, ice)
      
      ! Calculate 3-D velocity solution from the DIVA
      CALL calc_3D_horizontal_velocities_DIVA( mesh, ice)
      
      ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
      CALL calc_3D_vertical_velocities( mesh, ice)
      
    ! == From the b-grid to the c-grid (used in the model)
      
      ! Map 3-D velocity from the b-grid to the c-grid
      CALL map_b_to_c_3D( mesh, ice%u_3D_DIVA_b, ice%u_3D_c)
      CALL map_b_to_c_3D( mesh, ice%v_3D_DIVA_b, ice%v_3D_c)
      
      ! Map vertically averaged velocity from the b-grid to the c-grid
      CALL map_b_to_c_2D( mesh, ice%u_vav_DIVA_b, ice%u_vav_c)
      CALL map_b_to_c_2D( mesh, ice%v_vav_DIVA_b, ice%v_vav_c)
      
      ! Map basal velocity from the b-grid to the c-grid
      CALL map_b_to_c_2D( mesh, ice%u_base_DIVA_b, ice%u_base_c)
      CALL map_b_to_c_2D( mesh, ice%v_base_DIVA_b, ice%v_base_c)
      
      ! Copy surface velocity from the 3-D fields
      ice%u_surf_c( mesh%ti1:mesh%ti2) = ice%u_3D_c( mesh%ti1:mesh%ti2,1)
      ice%v_surf_c( mesh%ti1:mesh%ti2) = ice%v_3D_c( mesh%ti1:mesh%ti2,1)
      CALL sync
      
    ! == From the b-grid to the a-grid (only for writing to output)
      
      ! Map 3-D velocity from the b-grid to the a-grid
      CALL map_b_to_a_3D( mesh, ice%u_3D_DIVA_b, ice%u_3D_a)
      CALL map_b_to_a_3D( mesh, ice%v_3D_DIVA_b, ice%v_3D_a)
      
      ! Map vertically averaged velocity from the b-grid to the a-grid
      CALL map_b_to_a_2D( mesh, ice%u_vav_DIVA_b, ice%u_vav_a)
      CALL map_b_to_a_2D( mesh, ice%v_vav_DIVA_b, ice%v_vav_a)
      
      ! Map basal velocity from the b-grid to the a-grid
      CALL map_b_to_a_2D( mesh, ice%u_base_DIVA_b, ice%u_base_a)
      CALL map_b_to_a_2D( mesh, ice%v_base_DIVA_b, ice%v_base_a)
      
      ! Copy surface velocity from the 3-D fields
      ice%u_surf_a( mesh%vi1:mesh%vi2) = ice%u_3D_a( mesh%vi1:mesh%vi2,1)
      ice%v_surf_a( mesh%vi1:mesh%vi2) = ice%v_3D_a( mesh%vi1:mesh%vi2,1)
      CALL sync
      
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
    ! Calculate the driving stress taud (in both x and y) on the b-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    
    ! Local variables:
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:    ), POINTER                ::  dHs_dx_b,  dHs_dy_b,  Hi_b
    INTEGER                                            :: wdHs_dx_b, wdHs_dy_b, wHi_b
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, dHs_dx_b, wdHs_dx_b)
    CALL allocate_shared_dp_1D( mesh%nTri, dHs_dy_b, wdHs_dy_b)
    CALL allocate_shared_dp_1D( mesh%nTri, Hi_b    , wHi_b    )
    
    ! Map ice thickness to the b-grid
    CALL map_a_to_b_2D( mesh, ice%Hi_a, Hi_b)
    
    ! Calculate surface slopes on the b-grid
    CALL ddx_a_to_b_2D( mesh, ice%Hs_a, dHs_dx_b)
    CALL ddy_a_to_b_2D( mesh, ice%Hs_a, dHs_dy_b)
    
    ! Calculate driving stresses on the b-grid
    DO ti =  mesh%ti1, mesh%ti2
      ice%taudx_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dx_b( ti)
      ice%taudy_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dy_b( ti)
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx_b)
    CALL deallocate_shared( wdHs_dy_b)
    CALL deallocate_shared( wHi_b    )
    
  END SUBROUTINE calculate_driving_stress
  SUBROUTINE calc_vertical_shear_strain_rates( mesh, ice)
    ! Calculate vertical shear rates (Lipscomb et al. 2019, Eq. 36)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: ti,k
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  visc_eff_3D_b
    INTEGER                                            :: wvisc_eff_3D_b
    REAL(dp), PARAMETER                                :: visc_min = 1E3_dp
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, visc_eff_3D_b, wvisc_eff_3D_b)
    
    ! Map 3-D effective viscosity to the b-grid
    CALL map_a_to_b_3D( mesh, ice%visc_eff_3D_a, visc_eff_3D_b)
    
    ! Calculate vertical shear strain rates
    DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        ice%du_dz_3D_b( ti,k) = (ice%taubx_b( ti) / MAX( visc_min, visc_eff_3D_b( ti,k))) * C%zeta( k)
        ice%dv_dz_3D_b( ti,k) = (ice%tauby_b( ti) / MAX( visc_min, visc_eff_3D_b( ti,k))) * C%zeta( k)
      END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wvisc_eff_3D_b)
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%du_dz_3D_b, 'ice%du_dz_3D_b', 'calc_vertical_shear_strain_rates')
    CALL check_for_NaN_dp_2D( ice%dv_dz_3D_b, 'ice%dv_dz_3D_b', 'calc_vertical_shear_strain_rates')
    
  END SUBROUTINE calc_vertical_shear_strain_rates
  SUBROUTINE calc_effective_viscosity( mesh, ice, u_b, v_b)
    ! Calculate 3D effective viscosity following Lipscomb et al. (2019), Eq. 2

    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b
    
    ! Local variables
    INTEGER                                            :: vi, k
    REAL(dp)                                           :: eps_sq
    REAL(dp), PARAMETER                                :: epsilon_sq_0 = 1E-15_dp   ! Normalisation term so that zero velocity gives non-zero viscosity
    REAL(dp), DIMENSION(C%nz)                          :: prof
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  du_dz_3D_a,  dv_dz_3D_a
    INTEGER                                            :: wdu_dz_3D_a, wdv_dz_3D_a

    ! Calculate effective strain components from horizontal stretching
    CALL ddx_b_to_a_2D( mesh, u_b, ice%du_dx_a)
    CALL ddy_b_to_a_2D( mesh, u_b, ice%du_dy_a)
    CALL ddx_b_to_a_2D( mesh, v_b, ice%dv_dx_a)
    CALL ddy_b_to_a_2D( mesh, v_b, ice%dv_dy_a)
    
    ! Map vertical shear strain rates to the a-grid
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, du_dz_3D_a, wdu_dz_3D_a)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, dv_dz_3D_a, wdv_dz_3D_a)
    CALL map_b_to_a_3D( mesh, ice%du_dz_3D_b, du_dz_3D_a)
    CALL map_b_to_a_3D( mesh, ice%dv_dz_3D_b, dv_dz_3D_a)

    DO vi = mesh%vi1, mesh%vi2

      DO k = 1, C%nz
    
        ! Calculate the total effective strain rate from L19, Eq. 21 
        eps_sq = ice%du_dx_a( vi)**2 + &
                 ice%dv_dy_a( vi)**2 + &
                 ice%du_dx_a( vi) * ice%dv_dy_a( vi) + &
                 0.25_dp * (ice%du_dy_a( vi) + ice%dv_dx_a( vi))**2 + &
                 0.25_dp * (du_dz_3D_a( vi,k)**2 + dv_dz_3D_a( vi,k)**2) + &
                 epsilon_sq_0
        
        ! Calculate effective viscosity on ab-nodes
        ice%visc_eff_3D_a( vi,k) = 0.5_dp * ice%A_flow_3D_a( vi,k)**(-1._dp/C%n_flow) * (eps_sq)**((1._dp - C%n_flow)/(2._dp*C%n_flow))

      END DO ! DO k = 1, C%nz
      
      ! Vertical integral
      prof = ice%visc_eff_3D_a( vi,:)
      CALL vertical_integrate( prof, ice%visc_eff_int_a( vi))
      
      ! Product term N = eta * H
      ice%N_a( vi) = ice%visc_eff_int_a( vi) * MAX( 0.1_dp, ice%Hi_a( vi))

    END DO  
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdu_dz_3D_a)
    CALL deallocate_shared( wdv_dz_3D_a)
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%visc_eff_3D_a,  'ice%visc_eff_3D_a' , 'calc_effective_viscosity')
    CALL check_for_NaN_dp_1D( ice%visc_eff_int_a, 'ice%visc_eff_int_a', 'calc_effective_viscosity')
    CALL check_for_NaN_dp_1D( ice%N_a,            'ice%N_a'           , 'calc_effective_viscosity')

  END SUBROUTINE calc_effective_viscosity
  SUBROUTINE calc_sliding_term_beta( mesh, ice, u_b, v_b)
    ! Calculate the sliding term beta in the SSA/DIVA using the specified sliding law
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b
    
    ! Local variables:
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(:    ), POINTER                ::  u_a,  v_a
    INTEGER                                            :: wu_a, wv_a
      
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, u_a, wu_a)
    CALL allocate_shared_dp_1D( mesh%nV, v_a, wv_a)
  
    ! Get velocities on the a-grid
    CALL map_b_to_a_2D( mesh, u_b, u_a)
    CALL map_b_to_a_2D( mesh, v_b, v_a)
    
    ! Calculate the basal friction coefficients beta on the a-grid
    CALL calc_sliding_law( mesh, ice, u_a, v_a, ice%beta_a)
    
    ! Limit beta to improve stability
    DO vi = mesh%vi1, mesh%vi2
      ice%beta_a( vi) = MIN( C%DIVA_beta_max, ice%beta_a( vi))
    END DO
    CALL sync
    
    ! LEGACY - Apply the flotation mask (this is how we did it before we introduced the PISM/CISM-style sub-grid grounded fraction)
    IF (.NOT. C%do_GL_subgrid_friction) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (ice%mask_ocean_a( vi) == 1) ice%beta_a( vi) = 0._dp
      END DO
      CALL sync
    END IF
        
    ! Clean up after yourself
    CALl deallocate_shared( wu_a)
    CALl deallocate_shared( wv_a)
    
    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_a, 'ice%beta_a', 'calc_sliding_term_beta')
    
  END SUBROUTINE calc_sliding_term_beta
  SUBROUTINE calc_F_integral( mesh, ice, n)
    ! Calculate the integral F2 in Lipscomb et al. (2019), Eq. 30
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: vi
    REAL(dp)                                           :: F_int_min
    REAL(dp), PARAMETER                                :: visc_min = 1E5_dp
    REAL(dp), DIMENSION(C%nz)                          :: prof

    ! Set a lower limit for F2 to improve numerical stability
    prof = (1._dp / visc_min) * C%zeta**n
    CALL vertical_integrate( prof, F_int_min)

    DO vi = mesh%vi1, mesh%vi2
      
      prof = (ice%Hi_a( vi) / ice%visc_eff_3D_a( vi,:)) * C%zeta**n
      CALL vertical_integrate( prof, ice%F2_a( vi))
      ice%F2_a( vi) = MAX( ice%F2_a( vi), F_int_min)
      
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_1D( ice%F2_a, 'ice%F2_a', 'calc_F_integral')

  END SUBROUTINE calc_F_integral
  SUBROUTINE calc_beta_eff( mesh, ice)
    ! Calculate the "effective basal friction" beta_eff, used in the DIVA
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: vi, ti

    ! Calculate beta_eff on the a-grid
    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! No basal sliding allowed, impose beta_eff derived from viscosity 

      DO vi = mesh%vi1, mesh%vi2
        ! Lipscomb et al., 2019, Eq. 35
        ice%beta_eff_a( vi) = 1._dp / ice%F2_a( vi)
      END DO
      CALL sync

    ELSE
    
      DO vi = mesh%vi1, mesh%vi2
        ! Lipscomb et al., 2019, Eq. 33
        ice%beta_eff_a( vi) = ice%beta_a( vi) / (1._dp + ice%beta_a( vi) * ice%F2_a( vi))
      END DO
      CALL sync

    END IF
    
    ! Map beta_eff from the a-grid to the b-grid
    CALL map_a_to_b_2D( mesh, ice%beta_eff_a, ice%beta_eff_b)
    
    ! Apply the sub-grid grounded fraction
    IF (C%do_GL_subgrid_friction) THEN
      DO ti = mesh%ti1, mesh%ti2
        ice%beta_eff_b( ti) = ice%beta_eff_b( ti) * ice%f_grnd_b( ti)**2
      END DO
      CALL sync
    END IF
    
    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_eff_a, 'ice%beta_eff_a', 'calc_beta_eff')
    CALL check_for_NaN_dp_1D( ice%beta_eff_b, 'ice%beta_eff_b', 'calc_beta_eff')

  END SUBROUTINE calc_beta_eff
  SUBROUTINE calc_basal_stress_DIVA( mesh, ice, u_b, v_b)
    ! Calculate the basal stress resulting from sliding (friction times velocity)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b, v_b
    
    ! Local variables:
    INTEGER                                            :: ti

    DO ti = mesh%ti1, mesh%ti2
      ice%taubx_b( ti) = ice%beta_eff_b( ti) * u_b( ti)
      ice%tauby_b( ti) = ice%beta_eff_b( ti) * v_b( ti)
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_1D( ice%taubx_b, 'ice%taubx_b', 'calc_basal_stress_DIVA')
    CALL check_for_NaN_dp_1D( ice%tauby_b, 'ice%tauby_b', 'calc_basal_stress_DIVA')

  END SUBROUTINE calc_basal_stress_DIVA
  SUBROUTINE calc_basal_velocities_DIVA( mesh, ice)
    ! Calculate basal sliding following Goldberg (2011), Eq. 34
    ! (or it can also be obtained from L19, Eq. 32 given ub*beta=taub)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:    ), POINTER                ::  F2_b
    INTEGER                                            :: wF2_b

    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! Set basal velocities to zero 
      ! (this comes out naturally more or less with beta_eff set as above, 
      !  but ensuring basal velocity is zero adds stability)
      ice%u_base_DIVA_b( mesh%ti1:mesh%ti2) = 0._dp
      ice%v_base_DIVA_b( mesh%ti1:mesh%ti2) = 0._dp
      CALL sync
      RETURN
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, F2_b    , wF2_b    )
    
    ! Map F2 to the b-grid
    CALL map_a_to_b_2D( mesh, ice%F2_a, F2_b)
    
    ! Calculate basal velocities on the b-grid
    DO ti = mesh%ti1, mesh%ti2
      ice%u_base_DIVA_b( ti) = ice%u_vav_DIVA_b( ti) - ice%taubx_b( ti) * F2_b( ti)
      ice%v_base_DIVA_b( ti) = ice%v_vav_DIVA_b( ti) - ice%tauby_b( ti) * F2_b( ti)
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wF2_b    )
    
    ! Safety
    CALL check_for_NaN_dp_1D( ice%u_base_DIVA_b, 'ice%u_base_DIVA_b', 'calc_basal_velocities_DIVA')
    CALL check_for_NaN_dp_1D( ice%v_base_DIVA_b, 'ice%v_base_DIVA_b', 'calc_basal_velocities_DIVA')
        
  END SUBROUTINE calc_basal_velocities_DIVA
  SUBROUTINE calc_3D_horizontal_velocities_DIVA( mesh, ice)
    ! Calculate the 3D horizontal velocity field (following Lipscomb et al., 2019, Eq. 29)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:    ), POINTER                ::  Hi_b
    INTEGER                                            :: wHi_b
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  visc_eff_3D_b,  F1_3D_b
    INTEGER                                            :: wvisc_eff_3D_b, wF1_3D_b
    REAL(dp), DIMENSION( C%nz)                         :: prof, F1_3D
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri,       Hi_b         , wHi_b         )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, visc_eff_3D_b, wvisc_eff_3D_b)
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, F1_3D_b      , wF1_3D_b      )
    
    ! Map ice thickness and 3-D effective viscosity to the b-grid
    CALL map_a_to_b_2D( mesh, ice%Hi_a         , Hi_b         )
    CALL map_a_to_b_3D( mesh, ice%visc_eff_3D_a, visc_eff_3D_b)
    
    ! Calculate F1_3D on the b-grid
    DO ti = mesh%ti1, mesh%ti2
      prof = (-Hi_b( ti) / visc_eff_3D_b( ti,:)) * C%zeta
      CALL vertical_integration_from_bottom_to_zeta( prof, F1_3D)
      F1_3D_b( ti,:) = F1_3D
    END DO
    CALL sync

    ! Calculate 3D horizontal velocity components
    DO ti = mesh%ti1, mesh%ti2
      ice%u_3D_DIVA_b( ti,:) = ice%u_base_DIVA_b( ti) + ice%taubx_b( ti) * F1_3D_b( ti,:)
      ice%v_3D_DIVA_b( ti,:) = ice%v_base_DIVA_b( ti) + ice%tauby_b( ti) * F1_3D_b( ti,:)
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%u_3D_DIVA_b, 'ice%u_3D_DIVA_b', 'calc_3D_horizontal_velocities_DIVA')
    CALL check_for_NaN_dp_2D( ice%v_3D_DIVA_b, 'ice%v_3D_DIVA_b', 'calc_3D_horizontal_velocities_DIVA')

  END SUBROUTINE calc_3D_horizontal_velocities_DIVA
  
  ! Routines for solving the "linearised" SSA/DIVA
  SUBROUTINE solve_SSADIVA_linearised( mesh, ice, u_b, v_b)
    ! Solve the "linearised" version of the SSA (i.e. assuming viscosity and basal stress are
    ! constant rather than functions of velocity) on the b-grid.
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
    ! Optionally, this can be simplified by neglecting all the "cross-terms" involving gradients of N:
    !
    !   4*N*d2u/dx2 + N*d2u/dy2 + 3*N*d2v/dxdy - beta*u = -taud_x
    !   4*N*d2v/dy2 + N*d2v/dx2 + 3*N*d2u/dxdy - beta*v = -taud_y
    ! 
    ! The left and right sides of these equations can then be divided by N to yield:
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
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: u_b
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: v_b

    ! Local variables
    INTEGER                                            :: tiuv1, tiuv2, tiuv, ti, tiu, tiv, edge_index, n, vi
    REAL(dp), DIMENSION(:    ), POINTER                ::  N_b,  dN_dx_b,  dN_dy_b
    INTEGER                                            :: wN_b, wdN_dx_b, wdN_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                ::  b_buv,  uv_buv
    INTEGER                                            :: wb_buv, wuv_buv
    INTEGER                                            :: nrows, ncols, nnz_per_row_max, nnz_max
    TYPE(type_sparse_matrix_CSR)                       :: A
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D(   mesh%nTri, N_b    , wN_b    )
    CALL allocate_shared_dp_1D(   mesh%nTri, dN_dx_b, wdN_dx_b)
    CALL allocate_shared_dp_1D(   mesh%nTri, dN_dy_b, wdN_dy_b)
    CALL allocate_shared_dp_1D( 2*mesh%nTri, b_buv  , wb_buv  )
    CALL allocate_shared_dp_1D( 2*mesh%nTri, uv_buv , wuv_buv )
    
    ! Calculate N, dN/dx, and dN/dy on the b-grid
    CALL map_a_to_b_2D( mesh, ice%N_a, N_b    )
    CALL ddx_a_to_b_2D( mesh, ice%N_a, dN_dx_b)
    CALL ddy_a_to_b_2D( mesh, ice%N_a, dN_dy_b)
    
  ! Construct the big matrix A
  ! ==========================
    
    ! Allocate shared memory for A
    ncols           = 2*mesh%nTri    ! from
    nrows           = 2*mesh%nTri    ! to
    nnz_per_row_max = 20
    
    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( A, nrows, ncols, nnz_max)
    
    ! Fill in matrix rows, right-hand side, and initial guess
    CALL partition_list( 2*mesh%nTri, par%i, par%n, tiuv1, tiuv2)
    DO tiuv = tiuv1, tiuv2
    
      ! Apply boundary conditions to all triangles that touch the domain boundary
      IF (MOD( tiuv,2) == 1) THEN
        ! u
        ti = (tiuv+1) / 2
      ELSE
        ! v
        ti = tiuv / 2
      END IF
      edge_index = 0
      DO n = 1, 3
        vi = mesh%Tri( ti,n)
        IF (mesh%edge_index( vi) > 0) edge_index = MAX( edge_index, mesh%edge_index( vi))
      END DO
      
      ! Fill matrix coefficients
      IF (MOD( tiuv,2) == 1) THEN
        ! u
      
        IF (edge_index == 0) THEN
          ! Free triangle: fill in matrix row for the SSA/DIVA
          IF (C%include_SSADIVA_crossterms) THEN
            CALL list_DIVA_matrix_coefficients_eq_1_free( mesh, ice, u_b, ti, N_b, dN_dx_b, dN_dy_b, A, b_buv, uv_buv)
          ELSE
            CALL list_DIVA_matrix_coefficients_eq_1_free_sans( mesh, ice, u_b, ti, N_b, A, b_buv, uv_buv)
          END IF
        ELSE
          ! Border triangle: apply boundary conditions
          CALL list_DIVA_matrix_coefficients_eq_1_boundary( mesh, u_b, ti, edge_index, A, b_buv, uv_buv)
        END IF
        
      ELSE ! IF (MOD( tiuv,2) == 1) THEN
        ! v
      
        IF (edge_index == 0) THEN
          ! Free triangle: fill in matrix row for the SSA/DIVA
          IF (C%include_SSADIVA_crossterms) THEN
            CALL list_DIVA_matrix_coefficients_eq_2_free( mesh, ice, v_b, ti, N_b, dN_dx_b, dN_dy_b, A, b_buv, uv_buv)
          ELSE
            CALL list_DIVA_matrix_coefficients_eq_2_free_sans( mesh, ice, v_b, ti, N_b, A, b_buv, uv_buv)
          END IF
        ELSE
          ! Border triangle: apply boundary conditions
          CALL list_DIVA_matrix_coefficients_eq_2_boundary( mesh, v_b, ti, edge_index, A, b_buv, uv_buv)
        END IF
        
      END IF ! IF (MOD( tiuv,2) == 1) THEN
      
    END DO ! DO tiuv = tiuv1, tiuv2
    CALL sync
    
    ! Combine results from the different processes
    CALL sort_columns_in_CSR_dist( A, tiuv1, tiuv2)
    CALL finalise_matrix_CSR_dist( A, tiuv1, tiuv2)
    
  ! Solve the matrix equation
  ! =========================
    
    CALL solve_matrix_equation_CSR( A, b_buv, uv_buv, &
      C%DIVA_choice_matrix_solver, &
      C%DIVA_SOR_nit             , &
      C%DIVA_SOR_tol             , &
      C%DIVA_SOR_omega           , &
      C%DIVA_PETSc_rtol          , &
      C%DIVA_PETSc_abstol)
    
  ! Get solution back on the b-grid
  ! ================================
    
    DO ti = mesh%ti1, mesh%ti2
      
      tiu = 2*ti - 1
      tiv = 2*ti
      
      u_b( ti) = uv_buv( tiu)
      v_b( ti) = uv_buv( tiv)
      
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_matrix_CSR( A)
    CALL deallocate_shared( wN_b    )
    CALL deallocate_shared( wdN_dx_b)
    CALL deallocate_shared( wdN_dy_b)
    CALL deallocate_shared( wb_buv  )
    CALL deallocate_shared( wuv_buv )

  END SUBROUTINE solve_SSADIVA_linearised
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_free(      mesh, ice, u_b, ti, N_b, dN_dx_b, dN_dy_b, A, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    !
    !   4*N*d2u/dx2 + 4*dN/dx*du/dx + 2*N*d2v/dxdy + 2*dN/dx*dv/dy + N*d2u/dy2 + dN/dy*du/dy + N*d2v/dxdy + dN/dy*dv/dx - beta*u = -taud_x
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b
    INTEGER,                             INTENT(IN)    :: ti
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: N_b, dN_dx_b, dN_dy_b
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: b_buv, uv_buv
    
    ! Local variables:
    INTEGER                                            :: tiu, tiv, k1, k2, k, tj, tju, tjv
    
    tiu = 2*ti - 1
    tiv = 2*ti
    
    k1 = mesh%M2_ddx_b_b%ptr( ti)
    k2 = mesh%m2_ddx_b_b%ptr( ti+1) - 1
    
    DO k = k1, k2
      
      tj  = mesh%M2_ddx_b_b%index( k)
      tju = 2*tj - 1
      tjv = 2*tj
      
      ! u-part
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tju
      A%val(   A%nnz) = 4._dp * N_b(     ti) * mesh%M2_d2dx2_b_b%val( k) + &
                        4._dp * dN_dx_b( ti) * mesh%M2_ddx_b_b%val(   k) + &
                        1._dp * N_b(     ti) * mesh%M2_d2dy2_b_b%val( k) + &
                        1._dp * dN_dy_b( ti) * mesh%M2_ddy_b_b%val(   k)
      
      IF (tju == tiu) THEN
        A%val( A%nnz) = A%val( A%nnz) - ice%beta_eff_b( ti)
      END IF
      
      ! v-part
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tjv
      A%val(   A%nnz) = 3._dp * N_b(     ti) * mesh%M2_d2dxdy_b_b%val( k) + &
                        2._dp * dN_dx_b( ti) * mesh%M2_ddy_b_b%val(    k) + &
                        1._dp * dN_dy_b( ti) * mesh%M2_ddx_b_b%val(    k)
      
    END DO
    
    ! Finalise this matrix row
    A%ptr( tiu+1 : A%m+1) = A%nnz+1
    
    ! Right-hand side and initial guess
    b_buv(  tiu) = -ice%taudx_b( ti)
    uv_buv( tiu) = u_b( ti)
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_free
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_free(      mesh, ice, v_b, ti, N_b, dN_dx_b, dN_dy_b, A, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    !
    !   4*N*d2v/dy2 + 4*dN/dy*dv/dy + 2*N*d2u/dxdy + 2*dN/dy*du/dx + N*d2v/dx2 + dN/dx*dv/dx + N*d2u/dxdy + dN/dx*du/dy - beta*v = -taud_y
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b
    INTEGER,                             INTENT(IN)    :: ti
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: N_b, dN_dx_b, dN_dy_b
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: b_buv, uv_buv
    
    ! Local variables:
    INTEGER                                            :: tiu, tiv, k1, k2, k, tj, tju, tjv
    
    tiu = 2*ti - 1
    tiv = 2*ti
    
    k1 = mesh%M2_ddx_b_b%ptr( ti)
    k2 = mesh%m2_ddx_b_b%ptr( ti+1) - 1
    
    DO k = k1, k2
      
      tj  = mesh%M2_ddx_b_b%index( k)
      tju = 2*tj - 1
      tjv = 2*tj
      
      ! v-part
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tjv
      A%val(   A%nnz) = 4._dp * N_b(     ti) * mesh%M2_d2dy2_b_b%val( k) + &
                        4._dp * dN_dy_b( ti) * mesh%M2_ddy_b_b%val(   k) + &
                        1._dp * N_b(     ti) * mesh%M2_d2dx2_b_b%val( k) + &
                        1._dp * dN_dx_b( ti) * mesh%M2_ddx_b_b%val(   k)
      
      IF (tjv == tiv) THEN
        A%val( A%nnz) = A%val( A%nnz) - ice%beta_eff_b( ti)
      END IF
      
      ! u-part
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tju
      A%val(   A%nnz) = 3._dp * N_b(     ti) * mesh%M2_d2dxdy_b_b%val( k) + &
                        2._dp * dN_dy_b( ti) * mesh%M2_ddx_b_b%val(    k) + &
                        1._dp * dN_dx_b( ti) * mesh%M2_ddy_b_b%val(    k)
      
    END DO
    
    ! Finalise this matrix row
    A%ptr( tiv+1 : A%m+1) = A%nnz+1
    
    ! Right-hand side and initial guess
    b_buv(  tiv) = -ice%taudy_b( ti)
    uv_buv( tiv) = v_b( ti)
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_free
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_free_sans( mesh, ice, u_b, ti, N_b, A, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    ! 
    ! 4*d2u/dx2 + d2u/dy2 + 3*d2v/dxdy - beta*u/N = -taud_x/N
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b
    INTEGER,                             INTENT(IN)    :: ti
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: N_b
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: b_buv, uv_buv
    
    ! Local variables:
    INTEGER                                            :: tiu, tiv, k1, k2, k, tj, tju, tjv
    
    tiu = 2*ti - 1
    tiv = 2*ti
    
    k1 = mesh%M2_ddx_b_b%ptr( ti)
    k2 = mesh%m2_ddx_b_b%ptr( ti+1) - 1
    
    DO k = k1, k2
      
      tj  = mesh%M2_ddx_b_b%index( k)
      tju = 2*tj - 1
      tjv = 2*tj
      
      ! u-part
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tju
      A%val(   A%nnz) = 4._dp * mesh%M2_d2dx2_b_b%val( k) + mesh%M2_d2dy2_b_b%val( k)
      
      IF (tju == tiu) THEN
        A%val( A%nnz) = A%val( A%nnz) - (ice%beta_eff_b( ti) / N_b( ti))
      END IF
      
      ! v-part
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tjv
      A%val(   A%nnz) = 3._dp * mesh%M2_d2dxdy_b_b%val( k)
      
    END DO
    
    ! Finalise this matrix row
    A%ptr( tiu+1 : A%m+1) = A%nnz+1
    
    ! Right-hand side and initial guess
    b_buv(  tiu) = -ice%taudx_b( ti) / N_b( ti)
    uv_buv( tiu) = u_b( ti)
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_free_sans
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_free_sans( mesh, ice, v_b, ti, N_b, A, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    ! 
    ! 4*d2v/dy2 + d2v/dx2 + 3*d2u/dxdy - beta*v/N = -taud_y/N
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b
    INTEGER,                             INTENT(IN)    :: ti
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: N_b
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: b_buv, uv_buv
    
    ! Local variables:
    INTEGER                                            :: tiu, tiv, k1, k2, k, tj, tju, tjv
    
    tiu = 2*ti - 1
    tiv = 2*ti
    
    k1 = mesh%M2_ddx_b_b%ptr( ti)
    k2 = mesh%m2_ddx_b_b%ptr( ti+1) - 1
    
    DO k = k1, k2
      
      tj  = mesh%M2_ddx_b_b%index( k)
      tju = 2*tj - 1
      tjv = 2*tj
      
      ! v-part
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tjv
      A%val(   A%nnz) = 4._dp * mesh%M2_d2dy2_b_b%val( k) + mesh%M2_d2dx2_b_b%val( k)
      
      IF (tjv == tiv) THEN
        A%val( A%nnz) = A%val( A%nnz) - (ice%beta_eff_b( ti) / N_b( ti))
      END IF
      
      ! u-part
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tju
      A%val(   A%nnz) = 3._dp * mesh%M2_d2dxdy_b_b%val( k)
      
    END DO
    
    ! Finalise this matrix row
    A%ptr( tiv+1 : A%m+1) = A%nnz+1
    
    ! Right-hand side and initial guess
    b_buv(  tiv) = -ice%taudy_b( ti) / N_b( ti)
    uv_buv( tiv) = v_b( ti)
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_free_sans
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_boundary(  mesh, u_b, ti, edge_index, A, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b
    INTEGER,                             INTENT(IN)    :: ti, edge_index
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: b_buv, uv_buv
    
    ! Local variables:
    INTEGER                                            :: tiu, tiv
    CHARACTER(LEN=256)                                 :: BC
    INTEGER                                            :: tj, tju, tjv, n, tti
    
    tiu = 2*ti - 1
    tiv = 2*ti
   
    ! Determine what kind of boundary conditions to apply 
    IF     (edge_index == 8 .OR. edge_index == 1 .OR. edge_index == 2) THEN
      ! North
      IF     (C%DIVA_boundary_BC_u_north == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_north == 'infinite') THEN
        BC = 'infinite'
      ELSE
        WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_1_boundary - ERROR: unknown DIVA_boundary_BC_u_north "', TRIM(C%DIVA_boundary_BC_u_north), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSEIF (edge_index == 3) THEN
      ! East
      IF     (C%DIVA_boundary_BC_u_east == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_east == 'infinite') THEN
        BC = 'infinite'
      ELSE
        WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_1_boundary - ERROR: unknown DIVA_boundary_BC_u_east "', TRIM(C%DIVA_boundary_BC_u_east), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSEIF (edge_index == 4 .OR. edge_index == 5 .OR. edge_index == 6) THEN
      ! South
      IF     (C%DIVA_boundary_BC_u_south == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_south == 'infinite') THEN
        BC = 'infinite'
      ELSE
        WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_1_boundary - ERROR: unknown DIVA_boundary_BC_u_south "', TRIM(C%DIVA_boundary_BC_u_south), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSEIF (edge_index == 7) THEN
      ! West
      IF     (C%DIVA_boundary_BC_u_west == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_west == 'infinite') THEN
        BC = 'infinite'
      ELSE
        WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_1_boundary - ERROR: unknown DIVA_boundary_BC_u_west "', TRIM(C%DIVA_boundary_BC_u_west), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSE
      WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_1_boundary - ERROR: invalid edge_index ', edge_index, ' at ti = ', ti
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    IF (BC == 'zero') THEN
      ! Let u = 0 at this domain boundary
      
      ! Matrix
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tiu
      A%val(   A%nnz) = 1._dp
    
      ! Finalise this matrix row
      A%ptr( tiu+1 : A%m+1) = A%nnz+1
      
      ! Right-hand side and initial guess
      b_buv(  tiu) = 0._dp
      uv_buv( tiu) = 0._dp
      
    ELSEIF (BC == 'infinite') THEN
      ! Let du/dx = 0 at this domain boundary
      
      ! Matrix
      
      ! Find number of neighbouring triangles
      n = 0
      DO tti = 1, 3
        tj = mesh%TriC( ti,tti)
        IF (tj == 0) CYCLE
        n = n+1
      END DO
      
      ! The triangle itself
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tiu
      A%val(   A%nnz) = -1._dp
      
      ! Neighbouring triangles
      DO tti = 1, 3
        tj = mesh%TriC( ti,tti)
        IF (tj == 0) CYCLE
        tju = 2*tj - 1
        tjv = 2*tj
        A%nnz =  A%nnz + 1
        A%index( A%nnz) = tju
        A%val(   A%nnz) = 1._dp / REAL( n,dp)
      END DO
    
      ! Finalise this matrix row
      A%ptr( tiu+1 : A%m+1) = A%nnz+1
      
      ! Right-hand side and initial guess
      b_buv(  tiu) = 0._dp
      uv_buv( tiu) = u_b( ti)
      
    ELSE
      WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_1_boundary - ERROR: invalid BC = "', BC, '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_boundary
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_boundary(  mesh, v_b, ti, edge_index, A, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b
    INTEGER,                             INTENT(IN)    :: ti, edge_index
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: b_buv, uv_buv
    
    ! Local variables:
    INTEGER                                            :: tiu, tiv
    CHARACTER(LEN=256)                                 :: BC
    INTEGER                                            :: tj, tju, tjv, n, tti
    
    tiu = 2*ti - 1
    tiv = 2*ti
   
    ! Determine what kind of boundary conditions to apply 
    IF     (edge_index == 8 .OR. edge_index == 1 .OR. edge_index == 2) THEN
      ! North
      IF     (C%DIVA_boundary_BC_v_north == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_north == 'infinite') THEN
        BC = 'infinite'
      ELSE
        WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_2_boundary - ERROR: unknown DIVA_boundary_BC_v_north "', TRIM(C%DIVA_boundary_BC_v_north), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSEIF (edge_index == 3) THEN
      ! East
      IF     (C%DIVA_boundary_BC_v_east == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_east == 'infinite') THEN
        BC = 'infinite'
      ELSE
        WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_2_boundary - ERROR: unknown DIVA_boundary_BC_v_east "', TRIM(C%DIVA_boundary_BC_v_east), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSEIF (edge_index == 4 .OR. edge_index == 5 .OR. edge_index == 6) THEN
      ! South
      IF     (C%DIVA_boundary_BC_v_south == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_south == 'infinite') THEN
        BC = 'infinite'
      ELSE
        WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_2_boundary - ERROR: unknown DIVA_boundary_BC_v_south "', TRIM(C%DIVA_boundary_BC_v_south), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSEIF (edge_index == 7) THEN
      ! West
      IF     (C%DIVA_boundary_BC_v_west == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_west == 'infinite') THEN
        BC = 'infinite'
      ELSE
        WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_2_boundary - ERROR: unknown DIVA_boundary_BC_v_west "', TRIM(C%DIVA_boundary_BC_v_west), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSE
      WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_2_boundary - ERROR: invalid edge_index ', edge_index, ' at ti = ', ti
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    IF (BC == 'zero') THEN
      ! Let u = 0 at this domain boundary
      
      ! Matrix
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tiv
      A%val(   A%nnz) = 1._dp
    
      ! Finalise this matrix row
      A%ptr( tiv+1 : A%m+1) = A%nnz+1
      
      ! Right-hand side and initial guess
      b_buv(  tiv) = 0._dp
      uv_buv( tiv) = 0._dp
      
    ELSEIF (BC == 'infinite') THEN
      ! Let du/dx = 0 at this domain boundary
      
      ! Matrix
      
      ! Find number of neighbouring triangles
      n = 0
      DO tti = 1, 3
        tj = mesh%TriC( ti,tti)
        IF (tj == 0) CYCLE
        n = n+1
      END DO
      
      ! The triangle itself
      A%nnz =  A%nnz + 1
      A%index( A%nnz) = tiv
      A%val(   A%nnz) = -1._dp
      
      ! Neighbouring triangles
      DO tti = 1, 3
        tj = mesh%TriC( ti,tti)
        IF (tj == 0) CYCLE
        tju = 2*tj - 1
        tjv = 2*tj
        A%nnz =  A%nnz + 1
        A%index( A%nnz) = tjv
        A%val(   A%nnz) = 1._dp / REAL( n,dp)
      END DO
    
      ! Finalise this matrix row
      A%ptr( tiv+1 : A%m+1) = A%nnz+1
      
      ! Right-hand side and initial guess
      b_buv(  tiv) = 0._dp
      uv_buv( tiv) = v_b( ti)
      
    ELSE
      WRITE(0,*) 'list_DIVA_matrix_coefficients_eq_2_boundary - ERROR: invalid BC = "', BC, '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_boundary
    
  ! Some "administration" routines that help speed up and stabilise the SSA/DIVA solver
  SUBROUTINE apply_velocity_limits( mesh, u_b, v_b)
    ! Apply a velocity limit (for stability)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: u_b, v_b

    ! Local variables:
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:    ), POINTER                ::  uabs_b
    INTEGER                                            :: wuabs_b
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, uabs_b, wuabs_b)
    
    ! Calculate absolute speed
    DO ti = mesh%ti1, mesh%ti2
      uabs_b( ti) = SQRT( u_b( ti)**2 + v_b( ti)**2)
    END DO
    CALL sync
    
    ! Scale velocities accordingly
    DO ti = mesh%ti1, mesh%ti2
      IF (uabs_b( ti) > C%DIVA_vel_max) THEN
        u_b( ti) = u_b( ti) * C%DIVA_vel_max / uabs_b( ti)
        v_b( ti) = v_b( ti) * C%DIVA_vel_max / uabs_b( ti)
      END IF
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wuabs_b)

  END SUBROUTINE apply_velocity_limits
  SUBROUTINE relax_DIVA_visc_iterations( mesh, u_prev_b, v_prev_b, u_b, v_b, rel)
    ! Relax velocity solution with results of the previous viscosity iteration 
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_prev_b
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_prev_b
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: u_b
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: v_b
    REAL(dp),                            INTENT(IN)    :: rel
        
    ! Local variables:
    INTEGER                                            :: ti

    DO ti = mesh%ti1, mesh%ti2
      u_b( ti) = rel * u_b( ti) + (1._dp - rel) * u_prev_b( ti)
      v_b( ti) = rel * v_b( ti) + (1._dp - rel) * v_prev_b( ti)
    END DO
    CALL sync

  END SUBROUTINE relax_DIVA_visc_iterations
  SUBROUTINE calc_visc_iter_UV_resid( mesh, u_prev_b, v_prev_b, u_b, v_b, resid_UV)
    ! Check if the viscosity iteration has converged to a stable solution
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_prev_b
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_prev_b
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: u_b
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: v_b
    REAL(dp),                            INTENT(OUT)   :: resid_UV
    
    ! Local variables:
    INTEGER                                            :: ierr
    INTEGER                                            :: ti, nti
    REAL(dp)                                           :: res1, res2
    REAL(dp), PARAMETER                                :: DIVA_vel_tolerance = 1e-6   ! [m/a] only consider points with velocity above this tolerance limit

    ! Calculate the L2 norm based on velocity solution between previous
    ! and current viscosity iteration (as in Yelmo/SICOPOLIS)
    
    nti  = 0
    res1 = 0._dp
    res2 = 0._dp
    
    DO ti = mesh%ti1, mesh%ti2
      
      IF (ABS(u_b( ti)) > DIVA_vel_tolerance) THEN
        nti = nti + 1
        res1 = res1 + (u_b( ti) - u_prev_b( ti))**2._dp
        res2 = res2 + (u_b( ti) + u_prev_b( ti))**2._dp
      END IF
    
      IF (ABS(v_b( ti)) > DIVA_vel_tolerance) THEN
        nti = nti + 1
        res1 = res1 + (v_b( ti) - v_prev_b( ti))**2._dp
        res2 = res2 + (v_b( ti) + v_prev_b( ti))**2._dp
      END IF
      
    END DO
    
    ! Combine results from all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nti, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    IF (nti > 0) THEN
      res1 = SQRT( res1)
      res2 = SQRT( res2 )
      res2 = MAX(  res2,1E-8_dp)
      resid_UV = 2._dp * res1 / res2 
    ELSE 
      ! No points available for comparison, set residual equal to zero 
      resid_UV = 0._dp
    END IF

  END SUBROUTINE calc_visc_iter_UV_resid
  
  ! Initialise/remap data fields for the velocity solver(s)
  SUBROUTINE initialise_velocity_solver( mesh, ice)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    IF (C%choice_ice_dynamics == 'SIA' .OR. C%choice_ice_dynamics == 'SIA/SSA') THEN
      ! Data fields for the SIA
      
      CALL allocate_shared_dp_2D(   mesh%nAc , C%nz       , ice%u_3D_SIA_c            , ice%wu_3D_SIA_c           )
      CALL allocate_shared_dp_2D(   mesh%nAc , C%nz       , ice%v_3D_SIA_c            , ice%wv_3D_SIA_c           )
      
    END IF
      
    IF (C%choice_ice_dynamics == 'SSA' .OR. C%choice_ice_dynamics == 'SIA/SSA' .OR. C%choice_ice_dynamics == 'DIVA') THEN
      ! Data fields for the SSA / DIVA
      
      IF (C%choice_ice_dynamics == 'SSA' .OR. C%choice_ice_dynamics == 'SIA/SSA') THEN
        ! Velocity fields containing the SSA solution on the b-grid
        CALL allocate_shared_dp_1D(   mesh%nTri,              ice%u_base_SSA_b          , ice%wu_base_SSA_b         )
        CALL allocate_shared_dp_1D(   mesh%nTri,              ice%v_base_SSA_b          , ice%wv_base_SSA_b         )
      ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
        ! Velocity fields containing the DIVA solution on the b-grid
        CALL allocate_shared_dp_2D(   mesh%nTri, C%nz,        ice%u_3D_DIVA_b           , ice%wu_3D_DIVA_b          )
        CALL allocate_shared_dp_2D(   mesh%nTri, C%nz,        ice%v_3D_DIVA_b           , ice%wv_3D_DIVA_b          )
        CALL allocate_shared_dp_1D(   mesh%nTri,              ice%u_vav_DIVA_b          , ice%wu_vav_DIVA_b         )
        CALL allocate_shared_dp_1D(   mesh%nTri,              ice%v_vav_DIVA_b          , ice%wv_vav_DIVA_b         )
        CALL allocate_shared_dp_1D(   mesh%nTri,              ice%u_base_DIVA_b         , ice%wu_base_DIVA_b        )
        CALL allocate_shared_dp_1D(   mesh%nTri,              ice%v_base_DIVA_b         , ice%wv_base_DIVA_b        )
      END IF
      
      ! Useful data fields for solving the SSA/DIVA
      CALL allocate_shared_dp_1D(   mesh%nTri,              ice%taudx_b               , ice%wtaudx_b              )
      CALL allocate_shared_dp_1D(   mesh%nTri,              ice%taudy_b               , ice%wtaudy_b              )
      CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%du_dx_a               , ice%wdu_dx_a              )
      CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%du_dy_a               , ice%wdu_dy_a              )
      CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dv_dx_a               , ice%wdv_dx_a              )
      CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dv_dy_a               , ice%wdv_dy_a              )
      CALL allocate_shared_dp_2D(   mesh%nTri, C%nz,        ice%du_dz_3D_b            , ice%wdu_dz_3D_b           )
      CALL allocate_shared_dp_2D(   mesh%nTri, C%nz,        ice%dv_dz_3D_b            , ice%wdv_dz_3D_b           )
      CALL allocate_shared_dp_2D(   mesh%nV  , C%nz,        ice%visc_eff_3D_a         , ice%wvisc_eff_3D_a        )
      CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%visc_eff_int_a        , ice%wvisc_eff_int_a       )
      CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%N_a                   , ice%wN_a                  )
      CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%beta_a                , ice%wbeta_a               )
      CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%beta_eff_a            , ice%wbeta_eff_a           )
      CALL allocate_shared_dp_1D(   mesh%nTri,              ice%beta_eff_b            , ice%wbeta_eff_b           )
      CALL allocate_shared_dp_1D(   mesh%nTri,              ice%taubx_b               , ice%wtaubx_b              )
      CALL allocate_shared_dp_1D(   mesh%nTri,              ice%tauby_b               , ice%wtauby_b              )
      CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%F2_a                  , ice%wF2_a                 )
      CALL allocate_shared_dp_1D(   mesh%nTri,              ice%u_prev_b              , ice%wu_prev_b             )
      CALL allocate_shared_dp_1D(   mesh%nTri,              ice%v_prev_b              , ice%wv_prev_b             )
      
    END IF
    
  END SUBROUTINE initialise_velocity_solver
  SUBROUTINE remap_velocity_solver( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
!    ! Local variables:
!    INTEGER                                            :: int_dummy
!    
!    ! To prevent compiler warnings for unused variables
!    int_dummy = mesh_old%nV
!    int_dummy = mesh_new%nV
!    int_dummy = map%trilin%vi( 1,1)
!    
!    IF (C%choice_ice_dynamics == 'SIA' .OR. C%choice_ice_dynamics == 'SIA/SSA') THEN
!      ! Data fields for the SIA
!      
!      CALL reallocate_shared_dp_2D(   mesh_new%nAc     , C%nz       , ice%u_3D_SIA_c            , ice%wu_3D_SIA_c           )
!      CALL reallocate_shared_dp_2D(   mesh_new%nAc     , C%nz       , ice%v_3D_SIA_c            , ice%wv_3D_SIA_c           )
!      
!    END IF
!      
!    IF (C%choice_ice_dynamics == 'SSA' .OR. C%choice_ice_dynamics == 'SIA/SSA' .OR. C%choice_ice_dynamics == 'DIVA') THEN
!      ! Data fields for the SSA / DIVA
!      
!      IF (C%choice_ice_dynamics == 'SSA' .OR. C%choice_ice_dynamics == 'SIA/SSA') THEN
!        ! Velocity fields containing the SSA solution
!        CALL reallocate_shared_dp_1D(   mesh_new%nV      ,              ice%u_SSA_a               , ice%wu_SSA_a              )
!        CALL reallocate_shared_dp_1D(   mesh_new%nV      ,              ice%v_SSA_a               , ice%wv_SSA_a              )
!        CALL reallocate_shared_dp_1D(   mesh_new%nAc     ,              ice%u_SSA_c               , ice%wu_SSA_c              )
!        CALL reallocate_shared_dp_1D(   mesh_new%nAc     ,              ice%v_SSA_c               , ice%wv_SSA_c              )
!      END IF
!      
!      ! Useful data fields for solving the SSA/DIVA
!      CALL reallocate_shared_dp_1D(   mesh_new%nVAaAc  ,              ice%Hi_ac                 , ice%wHi_ac                )
!      CALL reallocate_shared_dp_1D(   mesh_new%nVAaAc  ,              ice%taudx_ac              , ice%wtaudx_ac             )
!      CALL reallocate_shared_dp_1D(   mesh_new%nVAaAc  ,              ice%taudy_ac              , ice%wtaudy_ac             )
!      CALL reallocate_shared_dp_1D(   mesh_new%nTriAaAc,              ice%Hi_bb                 , ice%wHi_bb                )
!      CALL reallocate_shared_dp_1D(   mesh_new%nTriAaAc,              ice%du_dx_bb              , ice%wdu_dx_bb             )
!      CALL reallocate_shared_dp_1D(   mesh_new%nTriAaAc,              ice%du_dy_bb              , ice%wdu_dy_bb             )
!      CALL reallocate_shared_dp_1D(   mesh_new%nTriAaAc,              ice%dv_dx_bb              , ice%wdv_dx_bb             )
!      CALL reallocate_shared_dp_1D(   mesh_new%nTriAaAc,              ice%dv_dy_bb              , ice%wdv_dy_bb             )
!      CALL reallocate_shared_dp_2D(   mesh_new%nVAaAc  , C%nz,        ice%du_dz_3D_ac           , ice%wdu_dz_3D_ac          )
!      CALL reallocate_shared_dp_2D(   mesh_new%nVAaAc  , C%nz,        ice%dv_dz_3D_ac           , ice%wdv_dz_3D_ac          )
!      CALL reallocate_shared_dp_2D(   mesh_new%nVAaAc  , C%nz,        ice%A_flow_3D_ac          , ice%wA_flow_3D_ac         )
!      CALL reallocate_shared_dp_2D(   mesh_new%nTriAaAc, C%nz,        ice%A_flow_3D_bb          , ice%wA_flow_3D_bb         )
!      CALL reallocate_shared_dp_2D(   mesh_new%nTriAaAc, C%nz,        ice%visc_eff_3D_bb        , ice%wvisc_eff_3D_bb       )
!      CALL reallocate_shared_dp_2D(   mesh_new%nVAaAc  , C%nz,        ice%visc_eff_3D_ac        , ice%wvisc_eff_3D_ac       )
!      CALL reallocate_shared_dp_1D(   mesh_new%nTriAaAc,              ice%visc_eff_int_bb       , ice%wvisc_eff_int_bb      )
!      CALL reallocate_shared_dp_1D(   mesh_new%nTriAaAc,              ice%N_bb                  , ice%wN_bb                 )
!      CALL reallocate_shared_dp_1D(   mesh_new%nVAaAc  ,              ice%beta_ac               , ice%wbeta_ac              )
!      CALL reallocate_shared_dp_1D(   mesh_new%nVAaAc  ,              ice%beta_eff_ac           , ice%wbeta_eff_ac          )
!      CALL reallocate_shared_dp_1D(   mesh_new%nVAaAc  ,              ice%taubx_ac              , ice%wtaubx_ac             )
!      CALL reallocate_shared_dp_1D(   mesh_new%nVAaAc  ,              ice%tauby_ac              , ice%wtauby_ac             )
!      CALL reallocate_shared_dp_1D(   mesh_new%nVAaAc  ,              ice%F2_ac                 , ice%wF2_ac                )
!      CALL reallocate_shared_dp_2D(   mesh_new%nVAaAc  , C%nz,        ice%F1_3D_ac              , ice%wF1_3D_ac             )
!      
!    END IF
    
  END SUBROUTINE remap_velocity_solver
  
END MODULE ice_velocity_module
