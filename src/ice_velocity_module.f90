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
  USE data_types_module,               ONLY: type_mesh, type_ice_model
  USE mesh_operators_module,           ONLY: map_a_to_c_2D, map_a_to_c_3D, ddx_a_to_c_2D, ddy_a_to_c_2D, &
                                             ddx_a_to_a_2D, ddy_a_to_a_2D, ddx_c_to_a_3D, ddy_c_to_a_3D
  USE utilities_module,                ONLY: vertical_integration_from_bottom_to_zeta, vertical_average

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

    ! Calculate 3D horizontal velocities
    DO aci = mesh%ci1, mesh%ci2
      
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      
      ! u
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
    ! Calculate ice velocities using the SSA.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
!    ! Local variables:
!    INTEGER                                            :: i,j
!    LOGICAL                                            :: set_velocities_to_zero
!    LOGICAL                                            :: has_converged
!    INTEGER                                            :: viscosity_iteration_i
!    REAL(dp)                                           :: resid_UV
!    REAL(dp)                                           :: umax_analytical, tauc_analytical
!    
!    ! Check that this routine is called correctly
!    IF (.NOT. (C%choice_ice_dynamics == 'SSA' .OR. C%choice_ice_dynamics == 'SIA/SSA')) THEN
!      IF (par%master) WRITE(0,*) 'ERROR - solve_SSA should only be called when choice_ice_dynamics is set to SSA or SIA/SSA!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!    
!    ! If there's no grounded ice anywhere, don't bother
!    set_velocities_to_zero = .FALSE.
!    IF (SUM( ice%mask_sheet_a) == 0) set_velocities_to_zero = .TRUE.
!    
!    ! If we're prescribing no sliding, set velocities to zero
!    IF (C%choice_sliding_law == 'no_sliding') set_velocities_to_zero = .TRUE.
!    
!    IF (set_velocities_to_zero) THEN
!      ice%u_SSA_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 0._dp
!      ice%v_SSA_cy( :,grid%i1:              grid%i2 ) = 0._dp
!      CALL sync
!      RETURN
!    END IF
!    
!    ! Calculate the driving stresses taudx, taudy
!    DO i = grid%i1, grid%i2
!    DO j = 1, grid%ny
!      IF (i < grid%nx) ice%taudx_cx( j,i) = -ice_density * grav * ice%Hi_cx( j,i) * ice%dHs_dx_cx( j,i)
!      IF (j < grid%ny) ice%taudy_cy( j,i) = -ice_density * grav * ice%Hi_cy( j,i) * ice%dHs_dy_cy( j,i)
!    END DO
!    END DO
!    CALL sync
!    
!    ! Calculate the basal yield stress tau_c
!    CALL calc_basal_conditions( grid, ice)
!    
!    ! Determine sub-grid grounded fractions for scaling the basal friction
!    CALL determine_grounded_fractions( grid, ice)
!    
!    ! Find analytical solution for the SSA icestream experiment (used only to print numerical error to screen)
!    CALL SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, ice%A_flow_vav_a( 1,1), 0._dp, umax_analytical, tauc_analytical)
!            
!    ! Initially set error very high 
!    ice%DIVA_err_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 1E5_dp
!    ice%DIVA_err_cy( :,grid%i1:              grid%i2 ) = 1E5_dp
!    CALL sync
!    
!    ! The viscosity iteration
!    viscosity_iteration_i = 0
!    has_converged         = .FALSE.
!    viscosity_iteration: DO WHILE (.NOT. has_converged)
!      viscosity_iteration_i = viscosity_iteration_i + 1
!      
!      ! Calculate the effective viscosity and the product term N = eta * H
!      CALL calc_effective_viscosity( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)
!            
!      ! Calculate the sliding term beta (on both the A and Cx/Cy grids)
!      CALL calc_sliding_term_beta( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)
!      
!      ! Set beta_eff equal to beta; this turns the DIVA into the SSA
!      ice%beta_eff_a( :,grid%i1:grid%i2) = ice%beta_a( :,grid%i1:grid%i2)
!    
!      ! Map beta_eff from the a-grid to the cx/cy-grids
!      CALL map_a_to_cx_2D( grid, ice%beta_eff_a, ice%beta_eff_cx)
!      CALL map_a_to_cy_2D( grid, ice%beta_eff_a, ice%beta_eff_cy)
!    
!      ! Apply the sub-grid grounded fraction
!      DO i = grid%i1, grid%i2
!      DO j = 1, grid%ny
!        IF (i < grid%nx) ice%beta_eff_cx( j,i) = ice%beta_eff_cx( j,i) * ice%f_grnd_cx( j,i)**2
!        IF (j < grid%ny) ice%beta_eff_cy( j,i) = ice%beta_eff_cy( j,i) * ice%f_grnd_cy( j,i)**2
!      END DO
!      END DO
!      CALL sync
!      
!      ! Store the previous solution so we can check for convergence later
!      ice%u_cx_prev( :,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_SSA_cx( :,grid%i1:MIN(grid%nx-1,grid%i2))
!      ice%v_cy_prev( :,grid%i1:              grid%i2 ) = ice%v_SSA_cy( :,grid%i1:              grid%i2 )
!      CALL sync
!      
!      ! Solve the linearised DIVA with the SICOPOLIS solver
!      CALL solve_DIVA_stag_linearised( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)
!
!      ! Apply velocity limits (both overflow and underflow) for improved stability
!      CALL apply_velocity_limits( grid, ice%u_SSA_cx, ice%v_SSA_cy)
!      
!      ! "relax" subsequent viscosity iterations for improved stability
!      CALL relax_DIVA_visc_iterations( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy, C%DIVA_visc_it_relax)
!      
!      ! Check if the viscosity iteration has converged
!      CALL calc_visc_iter_UV_resid( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy, resid_UV)
!      !IF (par%master) WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ': resid_UV = ', resid_UV, ', u = [', MINVAL(ice%u_SSA_cx), ' - ', MAXVAL(ice%u_SSA_cx), ']'
!      
!      IF (par%master .AND. C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'SSA_icestream') &
!        WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ': err = ', ABS(1._dp - MAXVAL(ice%u_SSA_cx) / umax_analytical), ': resid_UV = ', resid_UV
!
!      has_converged = .FALSE.
!      IF     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) THEN
!        has_converged = .TRUE.
!      ELSEIF (viscosity_iteration_i >= C%DIVA_visc_it_nit) THEN
!        has_converged = .TRUE.
!      END IF
!      
!      ! If needed, estimate the error in the velocity fields so we can
!      ! update the SICOPOLIS-style DIVA solving masks (for better efficiency)
!      IF (.NOT. has_converged) CALL estimate_visc_iter_UV_errors( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)
!      
!    END DO viscosity_iteration
!    
!    ! Calculate secondary velocities (surface, base, etc.)
!    CALL calc_secondary_velocities( grid, ice)
    
  END SUBROUTINE solve_SSA
  
  ! Calculate "secondary" velocities (surface, basal, vertically averaged, on the A-grid, etc.)
  SUBROUTINE calc_secondary_velocities( mesh, ice)
    ! Calculate "secondary" velocities (surface, basal, vertically averaged, on the A-grid, etc.)
      
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
      
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      ! No SIA, just the SSA
      
      ! Set basal velocity equal to SSA answer
      ice%u_base_c( mesh%ci1:mesh%ci2) = ice%u_SSA_c( mesh%ci1:mesh%ci2)
      ice%v_base_c( mesh%ci1:mesh%ci2) = ice%v_SSA_c( mesh%ci1:mesh%ci2)
      CALL sync
      
      ! No vertical variations in velocity
      DO aci = mesh%ci1, mesh%ci2
        ice%u_vav_c(  aci  ) = ice%u_SSA_c( aci)
        ice%v_vav_c(  aci  ) = ice%v_SSA_c( aci)
        ice%u_surf_c( aci  ) = ice%u_SSA_c( aci)
        ice%v_surf_c( aci  ) = ice%v_SSA_c( aci)
        ice%u_3D_c(   aci,:) = ice%u_SSA_c( aci)
        ice%v_3D_c(   aci,:) = ice%v_SSA_c( aci)
      END DO
      CALL sync
      
    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      ! Hybrid SIA/SSA
      
      ! Add the two fields together
      DO k = 1, C%nz
        ice%u_3D_c( mesh%ci1:mesh%ci2,k) = ice%u_3D_SIA_c( mesh%ci1:mesh%ci2,k) + ice%u_SSA_c( mesh%ci1:mesh%ci2)
        ice%v_3D_c( mesh%ci1:mesh%ci2,k) = ice%v_3D_SIA_c( mesh%ci1:mesh%ci2,k) + ice%v_SSA_c( mesh%ci1:mesh%ci2)
      END DO
      CALL sync
      
      ! Set basal velocity equal to SSA answer
      ice%u_base_c( mesh%ci1:mesh%ci2) = ice%u_SSA_c( mesh%ci1:mesh%ci2)
      ice%v_base_c( mesh%ci1:mesh%ci2) = ice%v_SSA_c( mesh%ci1:mesh%ci2)
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
    
    ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      ! DIVA
      
      CALL calc_3D_vertical_velocities( mesh, ice)
      
      ! Get surface velocity from the 3D fields
      ice%u_surf_c( mesh%ci1:mesh%ci2) = ice%u_3D_c( mesh%ci1:mesh%ci2,1)
      ice%v_surf_c( mesh%ci1:mesh%ci2) = ice%v_3D_c( mesh%ci1:mesh%ci2,1)
      CALL sync
      
    ELSE ! IF (C%choice_ice_dynamics == 'SIA')
    
      IF (par%master) WRITE(0,*) 'calc_secondary_velocities - ERROR: unknown choice_ice_dynamics "', C%choice_ice_dynamics, '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      
    END IF ! IF (C%choice_ice_dynamics == 'SIA')
      
    ! Get velocity components on the A-grid (only used for writing to output)
    CALL map_velocity_from_c_to_a_3D( mesh, ice, ice%u_3D_c,   ice%u_3D_a  )
    CALL map_velocity_from_c_to_a_3D( mesh, ice, ice%v_3D_c,   ice%v_3D_a  )
    CALL map_velocity_from_c_to_a_2D( mesh, ice, ice%u_vav_c,  ice%u_vav_a )
    CALL map_velocity_from_c_to_a_2D( mesh, ice, ice%v_vav_c,  ice%v_vav_a )
    CALL map_velocity_from_c_to_a_2D( mesh, ice, ice%u_surf_c, ice%u_surf_a)
    CALL map_velocity_from_c_to_a_2D( mesh, ice, ice%v_surf_c, ice%v_surf_a)
    CALL map_velocity_from_c_to_a_2D( mesh, ice, ice%u_base_c, ice%u_base_a)
    CALL map_velocity_from_c_to_a_2D( mesh, ice, ice%v_base_c, ice%v_base_a)
    
    ! Get absolute velocities on the A-grid (only used for writing to output)
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
  
  ! Map velocity components to the a-grid for writing to output (diagnostic only)
  SUBROUTINE map_velocity_from_c_to_a_2D( mesh, ice, u_c, u_a)
    ! Map velocity components from the staggered c-grid to the regular a-grid.
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
    ! Map velocity components from the staggered c-grid to the regular a-grid.
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
