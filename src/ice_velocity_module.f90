module ice_velocity_module

  ! Contains all the routines needed to calculate instantaneous ice velocities for the
  ! current modelled ice-sheet geometry.

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module,                only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module
  use parallel_module,                     only : par, sync, ierr, cerr, partition_list
  use utilities_module,                    only : check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                                  check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                                  vertical_average, vertical_integration_from_bottom_to_zeta, &
                                                  vertical_integrate, SSA_Schoof2006_analytical_solution
  use data_types_module,                   only : type_mesh, type_ice_model, type_sparse_matrix_CSR_dp, &
                                                  type_remapping_mesh_mesh
  use mesh_mapping_module,                 only : remap_field_dp_2D
  use mesh_operators_module,               only : map_a_to_b_2d, map_b_to_a_2d, map_a_to_b_3d, map_b_to_a_3D
  use mesh_operators_module,               only : ddx_a_to_a_2D, ddy_a_to_a_2D, ddx_b_to_a_3D, ddy_b_to_a_3D
  use mesh_operators_module,               only : ddx_a_to_b_2D, ddx_b_to_a_2D, ddy_a_to_b_2D, ddy_b_to_a_2D
  use mesh_operators_module,               only : apply_Neumann_BC_direct_a_2D
  use sparse_matrix_module,                only : allocate_matrix_CSR_dist, finalise_matrix_CSR_dist, &
                                                  solve_matrix_equation_CSR, deallocate_matrix_CSR
  use basal_conditions_and_sliding_module, only : calc_basal_conditions, calc_sliding_law
  use general_ice_model_data_module,       only : determine_grounded_fractions
  use reallocate_mod,                      only : reallocate_bounds
  use mpi_module,                          only : allgather_array

  implicit none

contains

! ===== Velocity solvers =====
! ============================

  ! === SIA ===
  ! ===========

  subroutine solve_SIA(  mesh, ice)
    ! Calculate ice velocities using the SIA.

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)      :: mesh
    type(type_ice_model), intent(inout)   :: ice

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'solve_SIA'
    integer                               :: ti, via, vib, vic
    real(dp), dimension(:), allocatable   :: Hi_b, dHs_dx_b, dHs_dy_b
    real(dp)                              :: D_0, w
    real(dp), dimension(C%nZ)             :: A_flow_3D_b, D_prof, D_deformation
    real(dp), dimension(C%nZ)             :: D_SIA_3D
    real(dp), parameter                   :: D_uv_3D_cutoff = -1E5_dp
    real(dp), dimension(:,:), allocatable :: A_flow_3D_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory
    allocate( Hi_b    ( mesh%ti1:mesh%ti2), source=0._dp)
    allocate( dHs_dx_b( mesh%ti1:mesh%ti2), source=0._dp)
    allocate( dHs_dy_b( mesh%ti1:mesh%ti2), source=0._dp)

    allocate( A_flow_3D_a( 1:mesh%nV, C%nz), source=0._dp)
    A_flow_3D_a(mesh%vi1:mesh%vi2,:) = C%uniform_flow_factor!ice%A_flow_3D_a
    call allgather_array(A_flow_3D_a)

    ! Get ice thickness, surface slopes, and ice flow factor on the b-grid
    call map_a_to_b_2D( mesh, ice%Hi_a, Hi_b    )
    call ddx_a_to_b_2D( mesh, ice%Hs_a, dHs_dx_b)
    call ddy_a_to_b_2D( mesh, ice%Hs_a, dHs_dy_b)

    ! Calculate 3D horizontal velocities on the b-grid
    do ti = mesh%ti1, mesh%ti2

      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      if (ice%mask_sheet_a( via) == 0 .and. &
          ice%mask_sheet_a( vib) == 0 .and. &
          ice%mask_sheet_a( vic) == 0) then
         ! No ice at any of the three triangle corners; set velocity to zero

        ice%u_3D_SIA_b( ti,:) = 0._dp
        ice%v_3D_SIA_b( ti,:) = 0._dp

      else

        ! Calculate staggered ice flow factor
        w           = 0._dp
        A_flow_3D_b = 0._dp
        if (ice%mask_sheet_a( via) == 1) then
          w           = w           + 1._dp
          A_flow_3D_b = A_flow_3D_b + A_flow_3D_a( via,:)
        end if
        if (ice%mask_sheet_a( vib) == 1) then
          w           = w           + 1._dp
          A_flow_3D_b = A_flow_3D_b + A_flow_3D_a( vib,:)
        end if
        if (ice%mask_sheet_a( vic) == 1) then
          w           = w           + 1._dp
          A_flow_3D_b = A_flow_3D_b + A_flow_3D_a( vic,:)
        end if
        A_flow_3D_b = A_flow_3D_b / w

        D_0           = (ice_density * grav * Hi_b( ti))**C%n_flow * ((dHs_dx_b( ti)**2 + dHs_dy_b( ti)**2))**((C%n_flow - 1._dp) / 2._dp)
        D_prof        = A_flow_3D_b * C%zeta**C%n_flow
        call vertical_integration_from_bottom_to_zeta( D_prof, D_deformation)
        D_deformation = 2._dp * Hi_b( ti) * D_deformation
        D_SIA_3D      = MAX( D_0 * D_deformation, D_uv_3D_cutoff)

        ice%u_3D_SIA_b( ti,:) = D_SIA_3D * dHs_dx_b( ti)
        ice%v_3D_SIA_b( ti,:) = D_SIA_3D * dHs_dy_b( ti)

      end if

    end do

    ! Calculate secondary velocities (surface, base, etc.)
    call calc_secondary_velocities( mesh, ice)

    ! Clean up after yourself
    deallocate( Hi_b       )
    deallocate( dHs_dx_b   )
    deallocate( dHs_dy_b   )
    deallocate( A_flow_3D_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SIA

  ! === SSA ===
  ! ===========

  subroutine solve_SSA(  mesh, ice)
    ! Calculate ice velocities using the SSA

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'solve_SSA'
    integer                             :: ti
    logical                             :: set_velocities_to_zero
    logical                             :: has_converged
    integer                             :: viscosity_iteration_i
    real(dp)                            :: resid_UV
    real(dp)                            :: umax_analytical, tauc_analytical

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Default
    set_velocities_to_zero = .false.

    ! If there's no grounded ice anywhere, don't bother
    if (sum( ice%mask_sheet_a) == 0) then
      set_velocities_to_zero = .true.
    end if

    call mpi_allreduce(mpi_in_place, set_velocities_to_zero, 1, &
                       mpi_logical, mpi_lor, mpi_comm_world, ierr)

    ! If we're prescribing no sliding, set velocities to zero
    if (C%choice_sliding_law == 'no_sliding') then
      set_velocities_to_zero = .true.
    end if

    if (set_velocities_to_zero) then
      ice%u_base_SSA_b( mesh%ti1:mesh%ti2) = 0._dp
      ice%v_base_SSA_b( mesh%ti1:mesh%ti2) = 0._dp
      call sync
      call calc_secondary_velocities( mesh, ice)
      call finalise_routine( routine_name)
      return
    end if

    ! === Driving stress ===
    ! ======================

    ! Calculate the driving stresses taudx, taudy
    call calculate_driving_stress( mesh, ice)

    ! === Basal friction ===
    ! ======================

    ! Calculate the basal yield stress tau_c
    call calc_basal_conditions( mesh, ice)

    ! Determine sub-mesh grounded fractions for scaling the basal friction
    call determine_grounded_fractions( mesh, ice)

    ! === Idealised case ===
    ! ======================

    ! Find analytical solution for the SSA icestream experiment (used only to print numerical error to screen)
    call SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, ice%A_flow_vav_a( mesh%vi1), &
                                             0._dp, umax_analytical, tauc_analytical)

    ! === Viscosity ===
    ! =================

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .false.

    do while (.not. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! Calculate the effective viscosity and the product term N = eta * H
      call calc_effective_viscosity( mesh, ice, ice%u_base_SSA_b, ice%v_base_SSA_b)

      ! Calculate the sliding term beta
      call calc_sliding_term_beta( mesh, ice, ice%u_base_SSA_b, ice%v_base_SSA_b)

      ! Set beta_eff equal to beta; this turns the DIVA into the SSA
      ice%beta_eff_a( mesh%vi1:mesh%vi2) = ice%beta_a( mesh%vi1:mesh%vi2)

      ! Map beta_eff from the a-grid to the cx/cy-grids
      call map_a_to_b_2D( mesh, ice%beta_eff_a, ice%beta_eff_b)

      ! Apply the sub-grid grounded fraction
      do ti = mesh%ti1, mesh%ti2
        ice%beta_eff_b( ti) = ice%beta_eff_b( ti) * ice%f_grnd_b( ti)**.5_dp
      end do

      ! Store the previous solution so we can check for convergence later
      ice%u_prev_b( mesh%ti1:mesh%ti2) = ice%u_base_SSA_b( mesh%ti1:mesh%ti2)
      ice%v_prev_b( mesh%ti1:mesh%ti2) = ice%v_base_SSA_b( mesh%ti1:mesh%ti2)

      ! Solve the linearised SSA
      call solve_SSADIVA_linearised( mesh, ice, ice%u_base_SSA_b, ice%v_base_SSA_b)

      ! Apply velocity limits (both overflow and underflow) for improved stability
      call apply_velocity_limits( mesh, ice%u_base_SSA_b, ice%v_base_SSA_b)

      ! "Relax" subsequent viscosity iterations for improved stability
      call relax_DIVA_visc_iterations( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_base_SSA_b, ice%v_base_SSA_b, C%DIVA_visc_it_relax)

      ! Check if the viscosity iteration has converged
      call calc_visc_iter_UV_resid( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_base_SSA_b, ice%v_base_SSA_b, resid_UV)

      ! if (par%master) then
      !   write(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, &
      !              ': resid_UV = ', resid_UV, &
      !              ', u = [', MINVAL(ice%u_base_SSA_b), &
      !              ' - ', MAXVAL(ice%u_base_SSA_b), ']'
      ! end if

      ! if (par%master .and. &
      !     C%choice_refgeo_init_ANT == 'idealised' .and. &
      !     C%choice_refgeo_init_idealised == 'SSA_icestream') then
      !   write(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, &
      !              ': err = ', abs(1._dp - maxval(ice%u_base_SSA_b) / umax_analytical), &
      !              ': resid_UV = ', resid_UV
      ! end if

      has_converged = .false.
      if     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) then
        has_converged = .true.
      elseif (viscosity_iteration_i >= C%DIVA_visc_it_nit) then
        has_converged = .true.
      end if

    end do

    ! Calculate secondary velocities (surface, base, etc.)
    call calc_secondary_velocities( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  END SUBROUTINE solve_SSA

  ! === DIVA ===
  ! ============

  subroutine solve_DIVA(  mesh, ice)
    ! Calculate ice velocities using the DIVA

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'solve_DIVA'
    logical                             :: set_velocities_to_zero
    logical                             :: has_converged
    integer                             :: viscosity_iteration_i
    real(dp)                            :: resid_UV
    real(dp)                            :: umax_analytical, tauc_analytical

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If there's no grounded ice anywhere, don't bother
    set_velocities_to_zero = .false.

    if (sum( ice%mask_sheet_a) == 0) then
      set_velocities_to_zero = .true.
    end if

    call mpi_allreduce(mpi_in_place,set_velocities_to_zero,1,mpi_logical,mpi_lor,mpi_comm_world,ierr)

    if (set_velocities_to_zero) then
      ice%u_vav_b( mesh%ti1:mesh%ti2) = 0._dp
      ice%v_vav_b( mesh%ti1:mesh%ti2) = 0._dp
      call calc_secondary_velocities( mesh, ice)
      call finalise_routine( routine_name)
      return
    end if

    ! Calculate the driving stresses taudx, taudy
    CALL calculate_driving_stress( mesh, ice)

    ! Calculate the basal yield stress tau_c
    CALL calc_basal_conditions( mesh, ice)

    ! Determine sub-mesh grounded fractions for scaling the basal friction
    CALL determine_grounded_fractions( mesh, ice)

    ! Find analytical solution for the SSA icestream experiment (used only to print numerical error to screen)
    CALL SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, ice%A_flow_vav_a( mesh%vi1), 0._dp, umax_analytical, tauc_analytical)

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .false.

    do while (.not. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! Calculate the vertical shear strain rates
      call calc_vertical_shear_strain_rates( mesh, ice)

      ! Calculate the effective viscosity and the product term N = eta * H
      call calc_effective_viscosity( mesh, ice, ice%u_vav_b, ice%v_vav_b)

      ! Calculate the sliding term beta
      call calc_sliding_term_beta( mesh, ice, ice%u_vav_b, ice%v_vav_b)

      ! Calculate the F-integral F2
      call calc_F_integral( mesh, ice, n = 2._dp)

      ! Calculate beta_eff
      call calc_beta_eff( mesh, ice)

      ! Store the previous solution so we can check for convergence later
      ice%u_prev_b( mesh%ti1:mesh%ti2) = ice%u_vav_b( mesh%ti1:mesh%ti2)
      ice%v_prev_b( mesh%ti1:mesh%ti2) = ice%v_vav_b( mesh%ti1:mesh%ti2)

      ! Solve the linearised DIVA
      call solve_SSADIVA_linearised( mesh, ice, ice%u_vav_b, ice%v_vav_b)

      ! Apply velocity limits (both overflow and underflow) for improved stability
      call apply_velocity_limits( mesh, ice%u_vav_b, ice%v_vav_b)

      ! "relax" subsequent viscosity iterations for improved stability
      call relax_DIVA_visc_iterations( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_vav_b, ice%v_vav_b, C%DIVA_visc_it_relax)

      ! Check if the viscosity iteration has converged
      call calc_visc_iter_UV_resid( mesh, ice%u_prev_b, ice%v_prev_b, ice%u_vav_b, ice%v_vav_b, resid_UV)

      if (par%master) then
        WRITE(0,*) '   DIVA - viscosity iteration ', viscosity_iteration_i, &
                   ': resid_UV = ', resid_UV, &
                   ', u = [', MINVAL(ice%u_vav_b), &
                   ' - ', MAXVAL(ice%u_vav_b), ']'
      end if

      if (par%master .and. &
          C%choice_refgeo_init_ANT == 'idealised' .and. &
          C%choice_refgeo_init_idealised == 'SSA_icestream') then
        write(0,*) '   DIVA - viscosity iteration ', viscosity_iteration_i, &
                   ': err = ', ABS(1._dp - MAXVAL(ice%u_vav_b) / umax_analytical), &
                   ': resid_UV = ', resid_UV
      end if

      has_converged = .false.
      if     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) then
        has_converged = .true.
      elseif (viscosity_iteration_i >= C%DIVA_visc_it_nit) then
        has_converged = .true.
      end if

      ! Calculate basal stress
      call calc_basal_stress_DIVA( mesh, ice, ice%u_vav_b, ice%v_vav_b)

    end do

    ! Calculate secondary velocities (surface, base, etc.)
    call calc_secondary_velocities( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_DIVA

  ! Calculate "secondary" velocities (surface, basal, vertically averaged, on the A-mesh, etc.)

  SUBROUTINE calc_secondary_velocities( mesh, ice)
    ! Calculate "secondary" velocities (surface, basal, vertically averaged, on the A-mesh, etc.)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_secondary_velocities'
    INTEGER                                            :: ti, vi, k
    REAL(dp), DIMENSION(C%nz)                          :: prof

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_ice_dynamics == 'SIA') THEN
      ! No SSA or sliding, just the SIA

      ! Set 3D velocity field equal to SIA answer
      ice%u_3D_b( mesh%ti1:mesh%ti2,:) = ice%u_3D_SIA_b( mesh%ti1:mesh%ti2,:)
      ice%v_3D_b( mesh%ti1:mesh%ti2,:) = ice%v_3D_SIA_b( mesh%ti1:mesh%ti2,:)

      ! Basal velocity is zero
      ice%u_base_b( mesh%ti1:mesh%ti2) = 0._dp
      ice%v_base_b( mesh%ti1:mesh%ti2) = 0._dp

      ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
      CALL calc_3D_vertical_velocities( mesh, ice)

      ! Copy surface velocity from the 3D fields
      ice%u_surf_b( mesh%ti1:mesh%ti2) = ice%u_3D_b( mesh%ti1:mesh%ti2,1)
      ice%v_surf_b( mesh%ti1:mesh%ti2) = ice%v_3D_b( mesh%ti1:mesh%ti2,1)

      ! Calculate vertically averaged velocities
      DO ti = mesh%ti1, mesh%ti2
        prof = ice%u_3D_b( ti,:)
        CALL vertical_average( prof, ice%u_vav_b( ti))
        prof = ice%v_3D_b( ti,:)
        CALL vertical_average( prof, ice%v_vav_b( ti))
      END DO

    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      ! No SIA, just the SSA

      ! No vertical variations in velocity; all fields are equal to the SSA answer
      DO ti = mesh%ti1, mesh%ti2
        ice%u_3D_b(   ti,:) = ice%u_base_SSA_b( ti)
        ice%v_3D_b(   ti,:) = ice%v_base_SSA_b( ti)
        ice%u_vav_b(  ti  ) = ice%u_base_SSA_b( ti)
        ice%v_vav_b(  ti  ) = ice%v_base_SSA_b( ti)
        ice%u_surf_b( ti  ) = ice%u_base_SSA_b( ti)
        ice%v_surf_b( ti  ) = ice%v_base_SSA_b( ti)
        ice%u_base_b( ti  ) = ice%u_base_SSA_b( ti)
        ice%v_base_b( ti  ) = ice%v_base_SSA_b( ti)
      END DO

      ! No vertical velocity
      ice%w_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp

    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      ! Hybrid SIA/SSA

      ! Set basal velocity equal to SSA answer
      ice%u_base_b( mesh%ti1:mesh%ti2) = ice%u_base_SSA_b( mesh%ti1:mesh%ti2)
      ice%v_base_b( mesh%ti1:mesh%ti2) = ice%v_base_SSA_b( mesh%ti1:mesh%ti2)

      ! Set 3-D velocities equal to the SIA solution
      ice%u_3D_b( mesh%ti1:mesh%ti2,:) = ice%u_3D_SIA_b( mesh%ti1:mesh%ti2,:)
      ice%v_3D_b( mesh%ti1:mesh%ti2,:) = ice%v_3D_SIA_b( mesh%ti1:mesh%ti2,:)

      ! Add the SSA contributions to the 3-D velocities
      DO k = 1, C%nz
        ice%u_3D_b( mesh%ti1:mesh%ti2,k) = ice%u_3D_b( mesh%ti1:mesh%ti2,k) + ice%u_base_b( mesh%ti1:mesh%ti2)
        ice%v_3D_b( mesh%ti1:mesh%ti2,k) = ice%v_3D_b( mesh%ti1:mesh%ti2,k) + ice%v_base_b( mesh%ti1:mesh%ti2)
      END DO

      ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
      CALL calc_3D_vertical_velocities( mesh, ice)

      ! Copy surface velocity from the 3D fields
      ice%u_surf_b( mesh%ti1:mesh%ti2) = ice%u_3D_b( mesh%ti1:mesh%ti2,1)
      ice%v_surf_b( mesh%ti1:mesh%ti2) = ice%v_3D_b( mesh%ti1:mesh%ti2,1)

      ! Calculate vertically averaged velocities
      DO ti = mesh%ti1, mesh%ti2
        prof = ice%u_3D_b( ti,:)
        CALL vertical_average( prof, ice%u_vav_b( ti))
        prof = ice%v_3D_b( ti,:)
        CALL vertical_average( prof, ice%v_vav_b( ti))
      END DO

    ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      ! DIVA

      ! Calculate basal velocity from depth-averaged solution and basal stress on the b-grid
      CALL calc_basal_velocities_DIVA( mesh, ice)

      ! Calculate 3-D velocity solution from the DIVA
      CALL calc_3D_horizontal_velocities_DIVA( mesh, ice)

      ! Copy surface velocity from the 3D fields
      ice%u_surf_b( mesh%ti1:mesh%ti2) = ice%u_3D_b( mesh%ti1:mesh%ti2,1)
      ice%v_surf_b( mesh%ti1:mesh%ti2) = ice%v_3D_b( mesh%ti1:mesh%ti2,1)

      ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
      CALL calc_3D_vertical_velocities( mesh, ice)

    ELSE
      CALL crash('unknown choice_ice_dynamics "' // TRIM( C%choice_ice_dynamics) // '"!')
    END IF

    ! Map velocity components to the a-grid
    CALL map_velocities_b_to_a_3D( mesh, ice%u_3D_b  , ice%v_3D_b  , ice%u_3D_a  , ice%v_3D_a  )
    CALL map_velocities_b_to_a_2D( mesh, ice%u_vav_b , ice%v_vav_b , ice%u_vav_a , ice%v_vav_a )
    CALL map_velocities_b_to_a_2D( mesh, ice%u_surf_b, ice%v_surf_b, ice%u_surf_a, ice%v_surf_a)
    CALL map_velocities_b_to_a_2D( mesh, ice%u_base_b, ice%v_base_b, ice%u_base_a, ice%v_base_a)

    ! Calculate absolute velocities
    DO ti = mesh%ti1, mesh%ti2
      ice%uabs_vav_b(  ti) = SQRT( ice%u_vav_b(  ti)**2 + ice%v_vav_b(  ti)**2)
      ice%uabs_surf_b( ti) = SQRT( ice%u_surf_b( ti)**2 + ice%v_surf_b( ti)**2)
      ice%uabs_base_b( ti) = SQRT( ice%u_base_b( ti)**2 + ice%v_base_b( ti)**2)
    END DO
    DO vi = mesh%vi1, mesh%vi2
      ice%uabs_vav_a(  vi) = SQRT( ice%u_vav_a(  vi)**2 + ice%v_vav_a(  vi)**2)
      ice%uabs_surf_a( vi) = SQRT( ice%u_surf_a( vi)**2 + ice%v_surf_a( vi)**2)
      ice%uabs_base_a( vi) = SQRT( ice%u_base_a( vi)**2 + ice%v_base_a( vi)**2)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_secondary_velocities
  SUBROUTINE calc_3D_vertical_velocities( mesh, ice)
    ! Use simple conservation of mass to calculate the vertical velocity w_3D

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_3D_vertical_velocities'
    INTEGER                                            :: vi, k
    REAL(dp), DIMENSION(:    ), allocatable                ::  dHi_dx_a,  dHi_dy_a,  dHs_dx_a,  dHs_dy_a
    REAL(dp), DIMENSION(:,:  ), allocatable            ::  du_dx_3D_a,  dv_dy_3D_a
    REAL(dp)                                           :: dHbase_dx, dHbase_dy
    REAL(dp)                                           :: du_dx_k,   dv_dy_k
    REAL(dp)                                           :: du_dx_kp1, dv_dy_kp1
    REAL(dp)                                           :: w1, w2, w3, w4

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( dHi_dx_a  (mesh%vi1:mesh%vi2      ))
    allocate( dHi_dy_a  (mesh%vi1:mesh%vi2      ))
    allocate( dHs_dx_a  (mesh%vi1:mesh%vi2      ))
    allocate( dHs_dy_a  (mesh%vi1:mesh%vi2      ))
    allocate( du_dx_3D_a(mesh%vi1:mesh%vi2, C%nz))
    allocate( dv_dy_3D_a(mesh%vi1:mesh%vi2, C%nz))

    ! Calculate surface & ice thickness gradients
    CALL ddx_a_to_a_2D( mesh, ice%Hs_a  , dHs_dx_a  )
    CALL ddy_a_to_a_2D( mesh, ice%Hs_a  , dHs_dy_a  )
    CALL ddx_a_to_a_2D( mesh, ice%Hi_a  , dHi_dx_a  )
    CALL ddy_a_to_a_2D( mesh, ice%Hi_a  , dHi_dy_a  )

    ! Calculate strain rates
    CALL ddx_b_to_a_3D( mesh, ice%u_3D_b, du_dx_3D_a)
    CALL ddy_b_to_a_3D( mesh, ice%v_3D_b, dv_dy_3D_a)

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

    ! Clean up after yourself
    deallocate( dHi_dx_a  )
    deallocate( dHi_dy_a  )
    deallocate( dHs_dx_a  )
    deallocate( dHs_dy_a  )
    deallocate( du_dx_3D_a)
    deallocate( dv_dy_3D_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_3D_vertical_velocities

  ! Calculate some physical terms (basal yield stress, effective viscosity, etc.)
  SUBROUTINE calculate_driving_stress( mesh, ice)
    ! Calculate the driving stress taud (in both x and y) on the b-grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calculate_driving_stress'
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:    ), POINTER                ::  dHs_dx_b,  dHs_dy_b,  Hi_b
    INTEGER                                            :: wdHs_dx_b, wdHs_dy_b, wHi_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( dHs_dx_b(mesh%ti1:mesh%ti2))
    allocate( dHs_dy_b(mesh%ti1:mesh%ti2))
    allocate( Hi_b    (mesh%ti1:mesh%ti2))

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

    ! Clean up after yourself
    deallocate( dHs_dx_b)
    deallocate( dHs_dy_b)
    deallocate( Hi_b    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_driving_stress
  SUBROUTINE calc_vertical_shear_strain_rates( mesh, ice)
    ! Calculate vertical shear rates (Lipscomb et al. 2019, Eq. 36)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_vertical_shear_strain_rates'
    INTEGER                                            :: ti,k
    REAL(dp), DIMENSION(:,:  ), allocatable            ::  visc_eff_3D_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( visc_eff_3D_b( mesh%ti1:mesh%ti2, C%nz))

    ! Map 3-D effective viscosity to the b-grid
    CALL map_a_to_b_3D( mesh, ice%visc_eff_3D_a, visc_eff_3D_b)

    ! Calculate vertical shear strain rates
    DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        ice%du_dz_3D_b( ti,k) = (ice%taubx_b( ti) / visc_eff_3D_b( ti,k)) * C%zeta( k)
        ice%dv_dz_3D_b( ti,k) = (ice%tauby_b( ti) / visc_eff_3D_b( ti,k)) * C%zeta( k)
      END DO
    END DO

    ! Clean up after yourself
    deallocate( visc_eff_3D_b)

    ! Safety
    CALL check_for_NaN_dp_2D( ice%du_dz_3D_b, 'ice%du_dz_3D_b')
    CALL check_for_NaN_dp_2D( ice%dv_dz_3D_b, 'ice%dv_dz_3D_b')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_shear_strain_rates
  SUBROUTINE calc_effective_viscosity( mesh, ice, u_b, v_b)
    ! Calculate 3D effective viscosity following Lipscomb et al. (2019), Eq. 2

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_effective_viscosity'
    INTEGER                                            :: vi, k
    REAL(dp)                                           :: eps_sq
    REAL(dp), DIMENSION(C%nz)                          :: prof
    REAL(dp), DIMENSION(:,:  ), allocatable            :: du_dz_3D_a,  dv_dz_3D_a
    REAL(dp), DIMENSION(  :  ), allocatable            :: N_a

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate effective strain components from horizontal stretching
    CALL ddx_b_to_a_2D( mesh, u_b, ice%du_dx_a)
    CALL ddy_b_to_a_2D( mesh, u_b, ice%du_dy_a)
    CALL ddx_b_to_a_2D( mesh, v_b, ice%dv_dx_a)
    CALL ddy_b_to_a_2D( mesh, v_b, ice%dv_dy_a)

    ! Map vertical shear strain rates to the a-grid
    allocate( du_dz_3D_a(mesh%vi1:mesh%vi2, C%nz))
    allocate( dv_dz_3D_a(mesh%vi1:mesh%vi2, C%nz))
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
                 C%DIVA_epsilon_sq_0

        ! Calculate effective viscosity
        ice%visc_eff_3D_a( vi,k) = 0.5_dp * C%uniform_flow_factor**(-1._dp/C%n_flow) * (eps_sq)**((1._dp - C%n_flow)/(2._dp*C%n_flow))

        ! Safety
        ice%visc_eff_3D_a( vi,k) = MAX( ice%visc_eff_3D_a( vi,k), C%DIVA_visc_eff_min)

      END DO ! DO k = 1, C%nz

      ! Vertical integral
      prof = ice%visc_eff_3D_a( vi,:)
      CALL vertical_integrate( prof, ice%visc_eff_int_a( vi))

      ! Product term N = eta * H
      ice%N_a( vi) = ice%visc_eff_int_a( vi) * MAX( 0.1_dp, ice%Hi_a( vi))

    END DO

    ! Clean up after yourself
    deallocate( du_dz_3D_a)
    deallocate( dv_dz_3D_a)

    ! Apply Neumann boundary conditions to correct inaccurate solutions at the domain border
    allocate(N_a(mesh%nV))
    N_a(mesh%vi1:mesh%vi2) = ice%N_a
    call allgather_array(N_a)

    CALL apply_Neumann_BC_direct_a_2D( mesh, N_a)

    ice%N_a = N_a(mesh%vi1:mesh%vi2)
    deallocate(N_a)

    ! Safety
    CALL check_for_NaN_dp_2D( ice%visc_eff_3D_a,  'ice%visc_eff_3D_a' )
    CALL check_for_NaN_dp_1D( ice%visc_eff_int_a, 'ice%visc_eff_int_a')
    CALL check_for_NaN_dp_1D( ice%N_a,            'ice%N_a'           )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_term_beta'
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(:    ), allocatable            ::  u_a,  v_a

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( u_a(mesh%vi1:mesh%vi2))
    allocate( v_a(mesh%vi1:mesh%vi2))

    ! Get velocities on the a-grid
    CALL map_b_to_a_2D( mesh, u_b, u_a)
    CALL map_b_to_a_2D( mesh, v_b, v_a)

    ! Calculate the basal friction coefficients beta on the a-grid
    CALL calc_sliding_law( mesh, ice, u_a, v_a, ice%beta_a)

    ! Limit beta to improve stability
    DO vi = mesh%vi1, mesh%vi2
      ice%beta_a( vi) = MIN( C%DIVA_beta_max, ice%beta_a( vi))
    END DO

    ! LEGACY - Apply the flotation mask (this is how we did it before we introduced the PISM/CISM-style sub-grid grounded fraction)
    IF (.NOT. C%do_GL_subgrid_friction) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (ice%mask_ocean_a( vi) == 1) ice%beta_a( vi) = 0._dp
      END DO
    END IF

    ! Clean up after yourself
    deallocate( u_a)
    deallocate( v_a)

    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_a, 'ice%beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_term_beta
  SUBROUTINE calc_F_integral( mesh, ice, n)
    ! Calculate the integral F2 in Lipscomb et al. (2019), Eq. 30

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: n

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_F_integral'
    INTEGER                                            :: vi
    REAL(dp)                                           :: F_int_min
    REAL(dp), DIMENSION(C%nz)                          :: prof

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set a lower limit for F2 to improve numerical stability
    prof = (1._dp / C%DIVA_visc_eff_min) * C%zeta**n
    CALL vertical_integrate( prof, F_int_min)

    DO vi = mesh%vi1, mesh%vi2

      prof = (ice%Hi_a( vi) / ice%visc_eff_3D_a( vi,:)) * C%zeta**n
      CALL vertical_integrate( prof, ice%F2_a( vi))
      ice%F2_a( vi) = MAX( ice%F2_a( vi), F_int_min)

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%F2_a, 'ice%F2_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_F_integral
  SUBROUTINE calc_beta_eff( mesh, ice)
    ! Calculate the "effective basal friction" beta_eff, used in the DIVA

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_beta_eff'
    INTEGER                                            :: vi, ti

    ! Add routine to path
    CALL init_routine( routine_name)

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
    CALL check_for_NaN_dp_1D( ice%beta_eff_a, 'ice%beta_eff_a')
    CALL check_for_NaN_dp_1D( ice%beta_eff_b, 'ice%beta_eff_b')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_beta_eff
  SUBROUTINE calc_basal_stress_DIVA( mesh, ice, u_b, v_b)
    ! Calculate the basal stress resulting from sliding (friction times velocity)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN) :: u_b, v_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_stress_DIVA'
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
      ice%taubx_b( ti) = ice%beta_eff_b( ti) * u_b( ti)
      ice%tauby_b( ti) = ice%beta_eff_b( ti) * v_b( ti)
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%taubx_b, 'ice%taubx_b')
    CALL check_for_NaN_dp_1D( ice%tauby_b, 'ice%tauby_b')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_stress_DIVA
  SUBROUTINE calc_basal_velocities_DIVA( mesh, ice)
    ! Calculate basal sliding following Goldberg (2011), Eq. 34
    ! (or it can also be obtained from L19, Eq. 32 given ub*beta=taub)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_velocities_DIVA'
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:    ), allocatable            ::  F2_b

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! Set basal velocities to zero
      ! (this comes out naturally more or less with beta_eff set as above,
      !  but ensuring basal velocity is zero adds stability)
      ice%u_base_b( mesh%ti1:mesh%ti2) = 0._dp
      ice%v_base_b( mesh%ti1:mesh%ti2) = 0._dp
      CALL sync
      RETURN
    END IF

    ! Allocate shared memory
    allocate( F2_b(mesh%ti1:mesh%ti2))

    ! Map F2 to the b-grid
    CALL map_a_to_b_2D( mesh, ice%F2_a, F2_b)

    ! Calculate basal velocities on the b-grid
    DO ti = mesh%ti1, mesh%ti2
      ice%u_base_b( ti) = ice%u_vav_b( ti) - ice%taubx_b( ti) * F2_b( ti)
      ice%v_base_b( ti) = ice%v_vav_b( ti) - ice%tauby_b( ti) * F2_b( ti)
    END DO

    ! Clean up after yourself
    deallocate( F2_b )

    ! Safety
    CALL check_for_NaN_dp_1D( ice%u_base_b, 'ice%u_base_b')
    CALL check_for_NaN_dp_1D( ice%v_base_b, 'ice%v_base_b')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_velocities_DIVA
  SUBROUTINE calc_3D_horizontal_velocities_DIVA( mesh, ice)
    ! Calculate the 3D horizontal velocity field (following Lipscomb et al., 2019, Eq. 29)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_3D_horizontal_velocities_DIVA'
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:    ), allocatable            ::  Hi_b
    REAL(dp), DIMENSION(:,:  ), allocatable            ::  visc_eff_3D_b,  F1_3D_b
    REAL(dp), DIMENSION( C%nz)                         :: prof, F1_3D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( Hi_b          ( mesh%ti1:mesh%ti2      ))
    allocate( visc_eff_3D_b ( mesh%ti1:mesh%ti2, C%nz))
    allocate( F1_3D_b       ( mesh%ti1:mesh%ti2, C%nz))

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
      ice%u_3D_b( ti,:) = ice%u_base_b( ti) + ice%taubx_b( ti) * F1_3D_b( ti,:)
      ice%v_3D_b( ti,:) = ice%v_base_b( ti) + ice%tauby_b( ti) * F1_3D_b( ti,:)
    END DO

    ! Safety
    CALL check_for_NaN_dp_2D( ice%u_3D_b, 'ice%u_3D_b')
    CALL check_for_NaN_dp_2D( ice%v_3D_b, 'ice%v_3D_b')

    ! Clean up after yourself
    deallocate( Hi_b         )
    deallocate( visc_eff_3D_b)
    deallocate( F1_3D_b      )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

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
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(INOUT) :: u_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(INOUT) :: v_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_SSADIVA_linearised'
    INTEGER                                            :: ti, nu, nv, t21,t22, pti1, pti2
    REAL(dp), DIMENSION(:    ), allocatable            ::  N_b,  dN_dx_b,  dN_dy_b
    REAL(dp), DIMENSION(:    ), allocatable            ::  b_buv,  uv_buv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( N_b    (   mesh%ti1:mesh%ti2))
    allocate( dN_dx_b(   mesh%ti1:mesh%ti2))
    allocate( dN_dy_b(   mesh%ti1:mesh%ti2))
    allocate( b_buv  (          1:mesh%nTri*2))
    allocate( uv_buv (          1:mesh%nTri*2))

    !(dividable by 2) ranges for b_buv,uv_buv
    t21 = 1+2*(mesh%ti1-1)
    t22 =   2*(mesh%ti2)

    ! PETSc ranges for b_buv,uv_buv
    call partition_list(mesh%nTri*2, par%i, par%n, pti1, pti2)

    ! Calculate N, dN/dx, and dN/dy on the b-grid
    CALL map_a_to_b_2D( mesh, ice%N_a, N_b    )
    CALL ddx_a_to_b_2D( mesh, ice%N_a, dN_dx_b)
    CALL ddy_a_to_b_2D( mesh, ice%N_a, dN_dy_b)

  ! Fill in the coefficients of the stiffness matrix
  ! ================================================

    DO ti = mesh%ti1, mesh%ti2

      IF (mesh%Tri_edge_index( ti) == 0) THEN
        ! Free triangle: fill in matrix row for the SSA/DIVA
        IF (C%include_SSADIVA_crossterms) THEN
          CALL calc_DIVA_matrix_coefficients_eq_1_free( mesh, ice, u_b, ti, N_b, dN_dx_b, dN_dy_b, b_buv(t21:t22), uv_buv)
          CALL calc_DIVA_matrix_coefficients_eq_2_free( mesh, ice, u_b, ti, N_b, dN_dx_b, dN_dy_b, b_buv(t21:t22), uv_buv)
        ELSE
          CALL calc_DIVA_matrix_coefficients_eq_1_free_sans( mesh, ice, u_b, ti, N_b, b_buv(t21:t22), uv_buv)
          CALL calc_DIVA_matrix_coefficients_eq_2_free_sans( mesh, ice, u_b, ti, N_b, b_buv(t21:t22), uv_buv)
        END IF
      ELSE
        ! Border triangle: apply boundary conditions
        CALL calc_DIVA_matrix_coefficients_eq_1_boundary( mesh, ice, u_b, ti, b_buv(t21:t22), uv_buv)
        CALL calc_DIVA_matrix_coefficients_eq_2_boundary( mesh, ice, u_b, ti, b_buv(t21:t22), uv_buv)
      END IF

    END DO

  ! Solve the matrix equation
  ! =========================

    ! TODO: The problem: b_buv and uv_buv are divideable by two but, petsc doesnt know that and
    ! will divide 6 into 3:3 for two procs for example, whereas the division should be 2:4 here in that case.
    ! Solution: point-to-point communication of boundaries, or fixing it in petsc (possible, but hard)
    call allgather_array(b_buv,t21,t22)
    call allgather_array(uv_buv,t21,t22)

    CALL solve_matrix_equation_CSR( ice%M_SSADIVA, b_buv(pti1:pti2), uv_buv(pti1:pti2), &
      C%DIVA_choice_matrix_solver, &
      C%DIVA_SOR_nit             , &
      C%DIVA_SOR_tol             , &
      C%DIVA_SOR_omega           , &
      C%DIVA_PETSc_rtol          , &
      C%DIVA_PETSc_abstol)

  ! Get solution back on the b-grid
  ! ================================

    DO ti = mesh%ti1, mesh%ti2

      nu = 2*ti - 1
      nv = 2*ti

      u_b( ti) = uv_buv( nu)
      v_b( ti) = uv_buv( nv)

    END DO

    ! Clean up after yourself
    deallocate( N_b    )
    deallocate( dN_dx_b)
    deallocate( dN_dy_b)
    deallocate( b_buv  )
    deallocate( uv_buv )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_SSADIVA_linearised
  SUBROUTINE calc_DIVA_matrix_coefficients_eq_1_free(      mesh, ice, u_b, ti, N_b, dN_dx_b, dN_dy_b, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    !
    !   4*N*d2u/dx2 + 4*dN/dx*du/dx + 2*N*d2v/dxdy + 2*dN/dx*dv/dy + N*d2u/dy2 + dN/dy*du/dy + N*d2v/dxdy + dN/dy*dv/dx - beta*u = -taud_x

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)          , INTENT(IN)    :: u_b
    INTEGER                                         , INTENT(IN)    :: ti
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)          , INTENT(IN)    :: N_b, dN_dx_b, dN_dy_b
    REAL(dp), DIMENSION(2*(mesh%ti1-1)+1:2*mesh%ti2), INTENT(INOUT) :: b_buv, uv_buv

    ! Local variables:
    INTEGER                                            :: nu, nnz, kk, ku, kv, k, mu

    nu = ice%ti2n_u( ti)

    nnz = (ice%M_SSADIVA%ptr( nu+1) - ice%M_SSADIVA%ptr( nu)) / 2

    DO kk = 1, nnz

      ku = ice%M_SSADIVA%ptr(   nu) + 2*kk - 2
      kv = ice%M_SSADIVA%ptr(   nu) + 2*kk - 1
      k  = mesh%M2_ddx_b_b_CSR%ptr( ti) +   kk - 1

      mu = ice%M_SSADIVA%index(   ku)

      ! u-part
      ice%M_SSADIVA%val( ku) = 4._dp * N_b(     ti) * mesh%M2_d2dx2_b_b_CSR%val( k) + &
                               4._dp * dN_dx_b( ti) * mesh%M2_ddx_b_b_CSR%val(   k) + &
                               1._dp * N_b(     ti) * mesh%M2_d2dy2_b_b_CSR%val( k) + &
                               1._dp * dN_dy_b( ti) * mesh%M2_ddy_b_b_CSR%val(   k)

      ! Sliding term on the diagonal
      IF (mu == nu) THEN
        ice%M_SSADIVA%val( ku) = ice%M_SSADIVA%val( ku) - ice%beta_eff_b( ti)
      END IF

      ! v-part
      ice%M_SSADIVA%val( kv) = 3._dp * N_b(     ti) * mesh%M2_d2dxdy_b_b_CSR%val( k) + &
                               2._dp * dN_dx_b( ti) * mesh%M2_ddy_b_b_CSR%val(    k) + &
                               1._dp * dN_dy_b( ti) * mesh%M2_ddx_b_b_CSR%val(    k)

    END DO

    ! Right-hand side and initial guess
    b_buv(  nu) = -ice%taudx_b( ti)
    uv_buv( nu) = u_b( ti)

  END SUBROUTINE calc_DIVA_matrix_coefficients_eq_1_free
  SUBROUTINE calc_DIVA_matrix_coefficients_eq_2_free(      mesh, ice, v_b, ti, N_b, dN_dx_b, dN_dy_b, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    !
    !   4*N*d2v/dy2 + 4*dN/dy*dv/dy + 2*N*d2u/dxdy + 2*dN/dy*du/dx + N*d2v/dx2 + dN/dx*dv/dx + N*d2u/dxdy + dN/dx*du/dy - beta*v = -taud_y

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)          , INTENT(IN)    :: v_b
    INTEGER                                         , INTENT(IN)    :: ti
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)          , INTENT(IN)    :: N_b, dN_dx_b, dN_dy_b
    REAL(dp), DIMENSION(2*(mesh%ti1-1)+1:2*mesh%ti2), INTENT(INOUT) :: b_buv, uv_buv

    ! Local variables:
    INTEGER                                            :: nv, nnz, kk, ku, kv, k, mv

    nv = ice%ti2n_v( ti)

    nnz = (ice%M_SSADIVA%ptr( nv+1) - ice%M_SSADIVA%ptr( nv)) / 2

    DO kk = 1, nnz

      ku = ice%M_SSADIVA%ptr(   nv) + 2*kk - 2
      kv = ice%M_SSADIVA%ptr(   nv) + 2*kk - 1
      k  = mesh%M2_ddx_b_b_CSR%ptr( ti) +   kk - 1

      mv = ice%M_SSADIVA%index(   kv)

      ! v-part
      ice%M_SSADIVA%val( kv) = 4._dp * N_b(     ti) * mesh%M2_d2dy2_b_b_CSR%val( k) + &
                               4._dp * dN_dy_b( ti) * mesh%M2_ddy_b_b_CSR%val(   k) + &
                               1._dp * N_b(     ti) * mesh%M2_d2dx2_b_b_CSR%val( k) + &
                               1._dp * dN_dx_b( ti) * mesh%M2_ddx_b_b_CSR%val(   k)

      ! Sliding term on the diagonal
      IF (mv == nv) THEN
        ice%M_SSADIVA%val( kv) = ice%M_SSADIVA%val( kv) - ice%beta_eff_b( ti)
      END IF

      ! u-part
      ice%M_SSADIVA%val( ku) = 3._dp * N_b(     ti) * mesh%M2_d2dxdy_b_b_CSR%val( k) + &
                               2._dp * dN_dy_b( ti) * mesh%M2_ddx_b_b_CSR%val(    k) + &
                               1._dp * dN_dx_b( ti) * mesh%M2_ddy_b_b_CSR%val(    k)

    END DO

    ! Right-hand side and initial guess
    b_buv(  nv) = -ice%taudy_b( ti)
    uv_buv( nv) = v_b( ti)

  END SUBROUTINE calc_DIVA_matrix_coefficients_eq_2_free
  SUBROUTINE calc_DIVA_matrix_coefficients_eq_1_free_sans( mesh, ice, u_b, ti, N_b, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    !
    ! 4*d2u/dx2 + d2u/dy2 + 3*d2v/dxdy - beta*u/N = -taud_x/N

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)          , INTENT(IN)    :: u_b
    INTEGER                                         , INTENT(IN)    :: ti
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)          , INTENT(IN)    :: N_b
    REAL(dp), DIMENSION(2*(mesh%ti1-1)+1:2*mesh%ti2), INTENT(INOUT) :: b_buv, uv_buv

    ! Local variables:
    INTEGER                                            :: nu, nnz, kk, ku, kv, k, mu

    nu = ice%ti2n_u( ti)

    nnz = (ice%M_SSADIVA%ptr( nu+1) - ice%M_SSADIVA%ptr( nu)) / 2

    DO kk = 1, nnz

      ku = ice%M_SSADIVA%ptr(   nu) + 2*kk - 2
      kv = ice%M_SSADIVA%ptr(   nu) + 2*kk - 1
      k  = mesh%M2_ddx_b_b_CSR%ptr( ti) +   kk - 1

      mu = ice%M_SSADIVA%index(   ku)

      ! u-part
      ice%M_SSADIVA%val( ku) = 4._dp * mesh%M2_d2dx2_b_b_CSR%val( k) + mesh%M2_d2dy2_b_b_CSR%val( k)

      ! Sliding term on the diagonal
      IF (mu == nu) THEN
        ice%M_SSADIVA%val( ku) = ice%M_SSADIVA%val( ku) - (ice%beta_eff_b( ti) / N_b( ti))
      END IF

      ! v-part
      ice%M_SSADIVA%val( kv) = 3._dp * mesh%M2_d2dxdy_b_b_CSR%val( k)

    END DO

    ! Right-hand side and initial guess
    b_buv(  nu) = -ice%taudx_b( ti) / N_b( ti)
    uv_buv( nu) = u_b( ti)

  END SUBROUTINE calc_DIVA_matrix_coefficients_eq_1_free_sans
  SUBROUTINE calc_DIVA_matrix_coefficients_eq_2_free_sans( mesh, ice, v_b, ti, N_b, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists
    !
    ! 4*d2v/dy2 + d2v/dx2 + 3*d2u/dxdy - beta*v/N = -taud_y/N

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)          , INTENT(IN)    :: v_b
    INTEGER                                         , INTENT(IN)    :: ti
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)          , INTENT(IN)    :: N_b
    REAL(dp), DIMENSION(2*(mesh%ti1-1)+1:2*mesh%ti2), INTENT(INOUT) :: b_buv, uv_buv

    ! Local variables:
    INTEGER                                            :: nv, nnz, kk, ku, kv, k, mv

    nv = ice%ti2n_v( ti)

    nnz = (ice%M_SSADIVA%ptr( nv+1) - ice%M_SSADIVA%ptr( nv)) / 2

    DO kk = 1, nnz

      ku = ice%M_SSADIVA%ptr(   nv) + 2*kk - 2
      kv = ice%M_SSADIVA%ptr(   nv) + 2*kk - 1
      k  = mesh%M2_ddx_b_b_CSR%ptr( ti) +   kk - 1

      mv = ice%M_SSADIVA%index(   kv)

      ! v-part
      ice%M_SSADIVA%val( kv) = 4._dp * mesh%M2_d2dy2_b_b_CSR%val( k) + mesh%M2_d2dx2_b_b_CSR%val( k)

      ! Sliding term on the diagonal
      IF (mv == nv) THEN
        ice%M_SSADIVA%val( kv) = ice%M_SSADIVA%val( kv) - (ice%beta_eff_b( ti) / N_b( ti))
      END IF

      ! u-part
      ice%M_SSADIVA%val( ku) = 3._dp * mesh%M2_d2dxdy_b_b_CSR%val( k)

    END DO

    ! Right-hand side and initial guess
    b_buv(  nv) = -ice%taudy_b( ti) / N_b( ti)
    uv_buv( nv) = v_b( ti)

  END SUBROUTINE calc_DIVA_matrix_coefficients_eq_2_free_sans
  SUBROUTINE calc_DIVA_matrix_coefficients_eq_1_boundary(  mesh, ice, u_b, ti, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)          , INTENT(IN)    :: u_b
    INTEGER                                         , INTENT(IN)    :: ti
    REAL(dp), DIMENSION(2*(mesh%ti1-1)+1:2*mesh%ti2), INTENT(INOUT) :: b_buv, uv_buv

    ! Local variables:
    INTEGER                                            :: edge_index
    CHARACTER(LEN=256)                                 :: BC
    INTEGER                                            :: nu, nnz, kk, ku, k, mu, n_neighbours

    edge_index = mesh%Tri_edge_index( ti)

    ! Determine what kind of boundary conditions to apply
    IF     (edge_index == 8 .OR. edge_index == 1 .OR. edge_index == 2) THEN
      ! North
      IF     (C%DIVA_boundary_BC_u_north == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_north == 'infinite') THEN
        BC = 'infinite'
      ELSE
        CALL crash('unknown DIVA_boundary_BC_u_north "' // TRIM( C%DIVA_boundary_BC_u_north) // '"!')
      END IF
    ELSEIF (edge_index == 3) THEN
      ! East
      IF     (C%DIVA_boundary_BC_u_east == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_east == 'infinite') THEN
        BC = 'infinite'
      ELSE
        CALL crash('unknown DIVA_boundary_BC_u_east "' // TRIM( C%DIVA_boundary_BC_u_east) // '"!')
      END IF
    ELSEIF (edge_index == 4 .OR. edge_index == 5 .OR. edge_index == 6) THEN
      ! South
      IF     (C%DIVA_boundary_BC_u_south == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_south == 'infinite') THEN
        BC = 'infinite'
      ELSE
        CALL crash('unknown DIVA_boundary_BC_u_south "' // TRIM( C%DIVA_boundary_BC_u_south) // '"!')
      END IF
    ELSEIF (edge_index == 7) THEN
      ! West
      IF     (C%DIVA_boundary_BC_u_west == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_west == 'infinite') THEN
        BC = 'infinite'
      ELSE
        CALL crash('unknown DIVA_boundary_BC_u_west "' // TRIM( C%DIVA_boundary_BC_u_west) // '"!')
      END IF
    ELSE
      CALL crash('invaled edge_index = {int_01} at ti = {int_02}!', int_01 = edge_index, int_02 = ti)
    END IF

    IF (BC == 'zero') THEN
      ! Let u = 0 at this domain boundary

      nu = ice%ti2n_u( ti)

      nnz = ice%M_SSADIVA%ptr( nu+1) - ice%M_SSADIVA%ptr( nu)
      n_neighbours = nnz - 1

      DO kk = 1, nnz

        ku = ice%M_SSADIVA%ptr(         nu) + kk - 1
        k  = mesh%M_Neumann_BC_b_CSR%ptr( ti) + kk - 1

        mu = ice%M_SSADIVA%index( ku)

        IF (mu == nu) THEN
          ! Diagonal: the triangle itself
          ice%M_SSADIVA%val( ku) = 1._dp
        ELSE
          ! Off-diagonal: neighbours
          ice%M_SSADIVA%val( ku) = 0._dp
        END IF

      END DO

      ! Right-hand side and initial guess
      b_buv(  nu) = 0._dp
      uv_buv( nu) = u_b( ti)

    ELSEIF (BC == 'infinite') THEN
      ! Let du/dx = 0 at this domain boundary

      nu = ice%ti2n_u( ti)

      nnz = ice%M_SSADIVA%ptr( nu+1) - ice%M_SSADIVA%ptr( nu)
      n_neighbours = nnz - 1

      DO kk = 1, nnz

        ku = ice%M_SSADIVA%ptr(         nu) + kk - 1
        k  = mesh%M_Neumann_BC_b_CSR%ptr( ti) + kk - 1

        mu = ice%M_SSADIVA%index( ku)

        IF (mu == nu) THEN
          ! Diagonal: the triangle itself
          ice%M_SSADIVA%val( ku) = REAL( n_neighbours,dp)
        ELSE
          ! Off-diagonal: neighbours
          ice%M_SSADIVA%val( ku) = -1._dp
        END IF

      END DO

      ! Right-hand side and initial guess
      b_buv(  nu) = 0._dp
      uv_buv( nu) = u_b( ti)

    ELSE
      CALL crash('invalid BC = "' // TRIM( BC) // '"!')
    END IF

  END SUBROUTINE calc_DIVA_matrix_coefficients_eq_1_boundary
  SUBROUTINE calc_DIVA_matrix_coefficients_eq_2_boundary(  mesh, ice, v_b, ti, b_buv, uv_buv)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2)          , INTENT(IN)    :: v_b
    INTEGER                                         , INTENT(IN)    :: ti
    REAL(dp), DIMENSION(2*(mesh%ti1-1)+1:2*mesh%ti2), INTENT(INOUT) :: b_buv, uv_buv

    ! Local variables:
    INTEGER                                            :: edge_index
    CHARACTER(LEN=256)                                 :: BC
    INTEGER                                            :: nv, nnz, kk, kv, k, mv, n_neighbours

    edge_index = mesh%Tri_edge_index( ti)

    ! Determine what kind of boundary conditions to apply
    IF     (edge_index == 8 .OR. edge_index == 1 .OR. edge_index == 2) THEN
      ! North
      IF     (C%DIVA_boundary_BC_v_north == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_north == 'infinite') THEN
        BC = 'infinite'
      ELSE
        CALL crash('unknown DIVA_boundary_BC_v_north "' // TRIM( C%DIVA_boundary_BC_v_north) // '"!')
      END IF
    ELSEIF (edge_index == 3) THEN
      ! East
      IF     (C%DIVA_boundary_BC_v_east == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_east == 'infinite') THEN
        BC = 'infinite'
      ELSE
        CALL crash('unknown DIVA_boundary_BC_v_east "' // TRIM( C%DIVA_boundary_BC_v_east) // '"!')
      END IF
    ELSEIF (edge_index == 4 .OR. edge_index == 5 .OR. edge_index == 6) THEN
      ! South
      IF     (C%DIVA_boundary_BC_v_south == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_south == 'infinite') THEN
        BC = 'infinite'
      ELSE
        CALL crash('unknown DIVA_boundary_BC_v_south "' // TRIM( C%DIVA_boundary_BC_v_south) // '"!')
      END IF
    ELSEIF (edge_index == 7) THEN
      ! West
      IF     (C%DIVA_boundary_BC_v_west == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_west == 'infinite') THEN
        BC = 'infinite'
      ELSE
        CALL crash('unknown DIVA_boundary_BC_v_west "' // TRIM( C%DIVA_boundary_BC_v_west) // '"!')
      END IF
    ELSE
      CALL crash('invalid edge_index = {int_01} at ti = {int_02}', int_01 = edge_index, int_02 = ti)
    END IF

    IF (BC == 'zero') THEN
      ! Let v = 0 at this domain boundary

      nv = ice%ti2n_v( ti)

      nnz = ice%M_SSADIVA%ptr( nv+1) - ice%M_SSADIVA%ptr( nv)
      n_neighbours = nnz - 1

      DO kk = 1, nnz

        kv = ice%M_SSADIVA%ptr(         nv) + kk - 1
        k  = mesh%M_Neumann_BC_b_CSR%ptr( ti) + kk - 1

        mv = ice%M_SSADIVA%index( kv)

        IF (mv == nv) THEN
          ! Diagonal: the triangle itself
          ice%M_SSADIVA%val( kv) = 1._dp
        ELSE
          ! Off-diagonal: neighbours
          ice%M_SSADIVA%val( kv) = 0._dp
        END IF

      END DO

      ! Right-hand side and initial guess
      b_buv(  nv) = 0._dp
      uv_buv( nv) = v_b( ti)

    ELSEIF (BC == 'infinite') THEN
      ! Let dv/dx = 0 at this domain boundary

      nv = ice%ti2n_v( ti)

      nnz = ice%M_SSADIVA%ptr( nv+1) - ice%M_SSADIVA%ptr( nv)
      n_neighbours = nnz - 1

      DO kk = 1, nnz

        kv = ice%M_SSADIVA%ptr(         nv) + kk - 1
        k  = mesh%M_Neumann_BC_b_CSR%ptr( ti) + kk - 1

        mv = ice%M_SSADIVA%index( kv)

        IF (mv == nv) THEN
          ! Diagonal: the triangle itself
          ice%M_SSADIVA%val( kv) = REAL( n_neighbours,dp)
        ELSE
          ! Off-diagonal: neighbours
          ice%M_SSADIVA%val( kv) = -1._dp
        END IF

      END DO

      ! Right-hand side and initial guess
      b_buv(  nv) = 0._dp
      uv_buv( nv) = v_b( ti)

    ELSE
      CALL crash('invalid BC = "' // TRIM(BC) // '"!')
    END IF

  END SUBROUTINE calc_DIVA_matrix_coefficients_eq_2_boundary

  ! Some "administration" routines that help speed up and stabilise the SSA/DIVA solver
  SUBROUTINE apply_velocity_limits( mesh, u_b, v_b)
    ! Apply a velocity limit (for stability)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2),          INTENT(INOUT) :: u_b, v_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_velocity_limits'
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:    ), allocatable            ::  uabs_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    allocate( uabs_b(mesh%ti1:mesh%ti2))

    ! Calculate absolute speed
    DO ti = mesh%ti1, mesh%ti2
      uabs_b( ti) = SQRT( u_b( ti)**2 + v_b( ti)**2)
    END DO

    ! Scale velocities accordingly
    DO ti = mesh%ti1, mesh%ti2
      IF (uabs_b( ti) > C%DIVA_vel_max) THEN
        u_b( ti) = u_b( ti) * C%DIVA_vel_max / uabs_b( ti)
        v_b( ti) = v_b( ti) * C%DIVA_vel_max / uabs_b( ti)
      END IF
    END DO

    ! Clean up after yourself
    deallocate( uabs_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_velocity_limits
  SUBROUTINE relax_DIVA_visc_iterations( mesh, u_prev_b, v_prev_b, u_b, v_b, rel)
    ! Relax velocity solution with results of the previous viscosity iteration

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_prev_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_prev_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(INOUT) :: u_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(INOUT) :: v_b
    REAL(dp),                            INTENT(IN)    :: rel

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'relax_DIVA_visc_iterations'
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
      u_b( ti) = rel * u_b( ti) + (1._dp - rel) * u_prev_b( ti)
      v_b( ti) = rel * v_b( ti) + (1._dp - rel) * v_prev_b( ti)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE relax_DIVA_visc_iterations
  SUBROUTINE calc_visc_iter_UV_resid( mesh, u_prev_b, v_prev_b, u_b, v_b, resid_UV)
    ! Check if the viscosity iteration has converged to a stable solution

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: u_prev_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN)    :: v_prev_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(INOUT) :: u_b
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(INOUT) :: v_b
    REAL(dp),                            INTENT(OUT)   :: resid_UV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_visc_iter_UV_resid'
    INTEGER                                            :: ierr
    INTEGER                                            :: ti, nti
    REAL(dp)                                           :: res1, res2
    REAL(dp), PARAMETER                                :: DIVA_vel_tolerance = 1e-6   ! [m/a] only consider points with velocity above this tolerance limit

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_visc_iter_UV_resid

  ! Map velocity fields from the b-grid to the a- and c-grids
  SUBROUTINE map_velocities_b_to_a_2D( mesh, u_b, v_b, u_a, v_a)
    ! Map velocity fields from the b-grid to the a-grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2), INTENT(IN) :: u_b, v_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT):: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_b_to_a_2D'
    real(dp), dimension(:), allocatable                :: l_u_b, l_v_b !local
    INTEGER                                            :: vi, vti, ti

    ! Add routine to path
    CALL init_routine( routine_name)

    allocate(l_u_b(1:mesh%nTri))
    allocate(l_v_b(1:mesh%nTri))
    l_u_b(mesh%ti1:mesh%ti2) = u_b
    l_v_b(mesh%ti1:mesh%ti2) = v_b
    call allgather_array(l_u_b)
    call allgather_array(l_v_b)

    DO vi = mesh%vi1, mesh%vi2

      u_a( vi) = 0._dp
      v_a( vi) = 0._dp

      DO vti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,vti)
        u_a( vi) = u_a( vi) + (l_u_b( ti) / REAL( mesh%niTri( vi),dp))
        v_a( vi) = v_a( vi) + (l_v_b( ti) / REAL( mesh%niTri( vi),dp))
      END DO

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_b_to_a_2D
  SUBROUTINE map_velocities_b_to_a_3D( mesh, u_b, v_b, u_a, v_a)
    ! Map velocity fields from the b-grid to the a-grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(mesh%ti1:mesh%ti2,C%nz ), INTENT(IN)    :: u_b, v_b
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2,C%nz ), INTENT(OUT)   :: u_a, v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_b_to_a_3D'
    real(dp), dimension(:,:), allocatable              :: l_u_b, l_v_b !local
    INTEGER                                            :: vi, vti, ti

    ! Add routine to path
    CALL init_routine( routine_name)

    allocate(l_u_b(1:mesh%nTri,C%nz))
    allocate(l_v_b(1:mesh%nTri,C%nz))
    l_u_b(mesh%ti1:mesh%ti2,:) = u_b
    l_v_b(mesh%ti1:mesh%ti2,:) = v_b
    call allgather_array(l_u_b)
    call allgather_array(l_v_b)


    DO vi = mesh%vi1, mesh%vi2

      u_a( vi,:) = 0._dp
      v_a( vi,:) = 0._dp

      DO vti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,vti)
        u_a( vi,:) = u_a( vi,:) + (l_u_b( ti,:) / REAL( mesh%niTri( vi),dp))
        v_a( vi,:) = v_a( vi,:) + (l_v_b( ti,:) / REAL( mesh%niTri( vi),dp))
      END DO

    END DO

    deallocate(l_u_b)
    deallocate(l_v_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_b_to_a_3D
  SUBROUTINE map_velocities_b_to_c_2D( mesh, u_b, v_b, u_c, v_c)
    ! Map velocity fields from the b-grid to the c-grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b, v_b
    REAL(dp), DIMENSION(mesh%ci1:mesh%ci2), INTENT(OUT)   :: u_c, v_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_b_to_c_2D'
    INTEGER                                            :: aci, til, tir

    ! Add routine to path
    CALL init_routine( routine_name)

    DO aci = mesh%ci1, mesh%ci2

      IF (mesh%edge_index_ac( aci) == 0) THEN
        ! Free edge

        ! Find the two adjacent triangles
        til = mesh%Aci( aci,5)
        tir = mesh%Aci( aci,6)

        u_c( aci) = (u_b( til) + u_b( tir)) / 2._dp
        v_c( aci) = (v_b( til) + v_b( tir)) / 2._dp

      ELSE
        ! Border edge

        ! Find the two adjacent triangles
        til = mesh%Aci( aci,5)
        tir = mesh%Aci( aci,6)

        IF (tir == 0) THEN
          u_c( aci) = u_b( til)
          v_c( aci) = v_b( til)
        ELSE
          u_c( aci) = u_b( tir)
          v_c( aci) = v_b( tir)
        END IF

      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_b_to_c_2D
#if 0
  SUBROUTINE map_velocities_b_to_c_3D( mesh, u_b, v_b, u_c, v_c)
    ! Map velocity fields from the b-grid to the c-grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_b, v_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: u_c, v_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_b_to_c_3D'
    INTEGER                                            :: aci, vi, vj, vl, vr, vti, ti, n1, n2, n3, til, tir

    ! Add routine to path
    CALL init_routine( routine_name)

    DO aci = mesh%ci1, mesh%ci2

      IF (mesh%edge_index_ac( aci) == 0) THEN
        ! Free edge

        ! Find the two adjacent triangles
        vi = mesh%Aci( aci,1)
        vj = mesh%Aci( aci,2)
        vl = mesh%Aci( aci,3)
        vr = mesh%Aci( aci,4)

        til = 0
        tir = 0
        DO vti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,vti)
          DO n1 = 1, 3
            n2 = n1 + 1
            IF (n2 == 4) n2 = 1
            n3 = n2 + 1
            IF (n3 == 4) n3 = 1
            IF (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj .AND. mesh%Tri( ti,n3) == vl) til = ti
            IF (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vr .AND. mesh%Tri( ti,n3) == vj) tir = ti
          END DO
        END DO

        u_c( aci,:) = (u_b( til,:) + u_b( tir,:)) / 2._dp
        v_c( aci,:) = (v_b( til,:) + v_b( tir,:)) / 2._dp

      ELSE
        ! Border edge

        ! Find the adjacent triangle
        vi = mesh%Aci( aci,1)
        vj = mesh%Aci( aci,2)
        vl = mesh%Aci( aci,3)

        til = 0
        DO vti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,vti)
          DO n1 = 1, 3
            n2 = n1 + 1
            IF (n2 == 4) n2 = 1
            n3 = n2 + 1
            IF (n3 == 4) n3 = 1
            IF ((mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj .AND. mesh%Tri( ti,n3) == vl) .OR. &
                (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vl .AND. mesh%Tri( ti,n3) == vj)) THEN
              til = ti
            END IF
          END DO
        END DO

        u_c( aci,:) = u_b( til,:)
        v_c( aci,:) = v_b( til,:)

      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_b_to_c_3D

#endif

  ! Initialise/remap data fields for the velocity solver(s)
  subroutine initialise_velocity_solver( mesh, ice)
    ! Allocate and initialise data fields for the velocity solver

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_velocity_solver'
    real(dp), dimension(:), allocatable ::  u_ISMIP_HOM
    integer                             :: ti
    real(dp)                            :: umin, umax, x, y

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_ice_dynamics == 'SIA' .or. &
        C%choice_ice_dynamics == 'SIA/SSA') then
      ! Data fields for the SIA

      allocate( ice%u_3D_SIA_b ( mesh%ti1:mesh%ti2, C%nz ))
      allocate( ice%v_3D_SIA_b ( mesh%ti1:mesh%ti2, C%nz ))

    end if

    if (C%choice_ice_dynamics == 'SSA' .or. &
        C%choice_ice_dynamics == 'SIA/SSA' .or. &
        C%choice_ice_dynamics == 'DIVA') then
      ! Data fields for the SSA / DIVA

      if (C%choice_ice_dynamics == 'SSA' .or. &
          C%choice_ice_dynamics == 'SIA/SSA') then
        ! Velocity fields containing the SSA solution on the b-grid
        allocate( ice%u_base_SSA_b ( mesh%ti1:mesh%ti2 ))
        allocate( ice%v_base_SSA_b ( mesh%ti1:mesh%ti2 ))
      end if

      ! Physical terms in the SSA/DIVA
      allocate( ice%taudx_b        ( mesh%ti1:mesh%ti2              ))
      allocate( ice%taudy_b        ( mesh%ti1:mesh%ti2              ))
      allocate( ice%du_dx_a        ( mesh%vi1:mesh%vi2              ))
      allocate( ice%du_dy_a        ( mesh%vi1:mesh%vi2              ))
      allocate( ice%dv_dx_a        ( mesh%vi1:mesh%vi2              ))
      allocate( ice%dv_dy_a        ( mesh%vi1:mesh%vi2              ))
      allocate( ice%du_dz_3D_b     ( mesh%ti1:mesh%ti2  , C%nz      ))
      allocate( ice%dv_dz_3D_b     ( mesh%ti1:mesh%ti2  , C%nz      ))
      allocate( ice%visc_eff_3D_a  ( mesh%vi1:mesh%vi2  , C%nz      ))
      allocate( ice%visc_eff_int_a ( mesh%vi1:mesh%vi2              ))
      allocate( ice%N_a            ( mesh%vi1:mesh%vi2              ))
      allocate( ice%beta_a         ( mesh%vi1:mesh%vi2              ))
      allocate( ice%beta_eff_a     ( mesh%vi1:mesh%vi2              ))
      allocate( ice%beta_eff_b     ( mesh%ti1:mesh%ti2              ))
      allocate( ice%taubx_b        ( mesh%ti1:mesh%ti2              ))
      allocate( ice%tauby_b        ( mesh%ti1:mesh%ti2              ))
      allocate( ice%F2_a           ( mesh%vi1:mesh%vi2              ))
      allocate( ice%u_prev_b       ( mesh%ti1:mesh%ti2              ))
      allocate( ice%v_prev_b       ( mesh%ti1:mesh%ti2              ))

      ! Some administrative stuff to make solving the SSA/DIVA more efficient
      allocate( ice%ti2n_u         ( mesh%ti1:mesh%ti2              ))
      allocate( ice%ti2n_v         ( mesh%ti1:mesh%ti2              ))
      allocate( ice%n2ti_uv     ( 2*(mesh%ti1-1)+1:2*mesh%ti2, 2    ))
      call initialise_matrix_conversion_lists(  mesh, ice)
      call initialise_SSADIVA_stiffness_matrix( mesh, ice)

    end if

    ! Initialise the ISMIP-HOM experiments for faster convergence
    if (C%choice_refgeo_init_ANT == 'idealised' .and. &
       (C%choice_refgeo_init_idealised == 'ISMIP_HOM_A' .or. &
        C%choice_refgeo_init_idealised == 'ISMIP_HOM_B' .or. &
        C%choice_refgeo_init_idealised == 'ISMIP_HOM_C' .or. &
        C%choice_refgeo_init_idealised == 'ISMIP_HOM_D')) then

      ! Allocate shared memory
      allocate( u_ISMIP_HOM (mesh%ti1:mesh%ti2))

      umin = 0._dp
      umax = 0._dp

      ! Calculate an approximation of the solution
      if (C%choice_refgeo_init_idealised == 'ISMIP_HOM_A') then

        if     (C%ISMIP_HOM_L == 160000._dp) then
          umin = 1.6_dp
          umax = 108.84_dp
        elseif (C%ISMIP_HOM_L ==  80000._dp) then
          umin = 1.75_dp
          umax = 95.73_dp
        elseif (C%ISMIP_HOM_L ==  40000._dp) then
          umin = 2.27_dp
          umax = 74.45_dp
        elseif (C%ISMIP_HOM_L ==  20000._dp) then
          umin = 4.49_dp
          umax = 49.99_dp
        elseif (C%ISMIP_HOM_L ==  10000._dp) then
          umin = 11.09_dp
          umax = 32.74_dp
        elseif (C%ISMIP_HOM_L ==   5000._dp) then
          umin = 18.38_dp
          umax = 24.79_dp
        end if

        do ti = mesh%ti1, mesh%ti2
          x = mesh%TriGC( ti,1)
          y = mesh%TriGC( ti,2)
          u_ISMIP_HOM( ti) = umin + (umax - umin) * ( (1._dp - (sin( x) * sin( y))) / 2._dp)**2
        end do

      elseif (C%choice_refgeo_init_idealised == 'ISMIP_HOM_B') then

        if     (C%ISMIP_HOM_L == 160000._dp) then
          umin = 1.57_dp
          umax = 111.41_dp
        elseif (C%ISMIP_HOM_L ==  80000._dp) then
          umin = 1.69_dp
          umax = 100.73_dp
        elseif (C%ISMIP_HOM_L ==  40000._dp) then
          umin = 2.09_dp
          umax = 82.3_dp
        elseif (C%ISMIP_HOM_L ==  20000._dp) then
          umin = 3.92_dp
          umax = 57.84_dp
        elseif (C%ISMIP_HOM_L ==  10000._dp) then
          umin = 10.23_dp
          umax = 35.2_dp
        elseif (C%ISMIP_HOM_L ==   5000._dp) then
          umin = 17.22_dp
          umax = 23.53_dp
        end if

        do ti = mesh%ti1, mesh%ti2
          x = mesh%TriGC( ti,1)
          y = mesh%TriGC( ti,2)
          u_ISMIP_HOM( ti) = umin + (umax - umin) * ( (1._dp - sin( x)) / 2._dp)**2
        end do

      elseif (C%choice_refgeo_init_idealised == 'ISMIP_HOM_C') then

        if     (C%ISMIP_HOM_L == 160000._dp) then
          umin = 8.77_dp
          umax = 143.45_dp
        elseif (C%ISMIP_HOM_L ==  80000._dp) then
          umin = 9.8_dp
          umax = 60.28_dp
        elseif (C%ISMIP_HOM_L ==  40000._dp) then
          umin = 11.84_dp
          umax = 28.57_dp
        elseif (C%ISMIP_HOM_L ==  20000._dp) then
          umin = 14.55_dp
          umax = 18.48_dp
        elseif (C%ISMIP_HOM_L ==  10000._dp) then
          umin = 15.7_dp
          umax = 16.06_dp
        elseif (C%ISMIP_HOM_L ==   5000._dp) then
          umin = 13.38_dp
          umax = 13.51_dp
        end if

        do ti = mesh%ti1, mesh%ti2
          x = mesh%TriGC( ti,1)
          y = mesh%TriGC( ti,2)
          u_ISMIP_HOM( ti) = umin + (umax - umin) * ( (1._dp - (SIN( x) * SIN( y))) / 2._dp)**2
        end do

      elseif (C%choice_refgeo_init_idealised == 'ISMIP_HOM_D') then

        if     (C%ISMIP_HOM_L == 160000._dp) then
          umin = 8.62_dp
          umax = 227.23_dp
        elseif (C%ISMIP_HOM_L ==  80000._dp) then
          umin = 9.65_dp
          umax = 94.79_dp
        elseif (C%ISMIP_HOM_L ==  40000._dp) then
          umin = 12.18_dp
          umax = 40.06_dp
        elseif (C%ISMIP_HOM_L ==  20000._dp) then
          umin = 15.28_dp
          umax = 20.29_dp
        elseif (C%ISMIP_HOM_L ==  10000._dp) then
          umin = 15.93_dp
          umax = 16.25_dp
        elseif (C%ISMIP_HOM_L ==   5000._dp) then
          umin = 14.43_dp
          umax = 14.59_dp
        end if

        do ti = mesh%ti1, mesh%ti2
          x = mesh%TriGC( ti,1)
          y = mesh%TriGC( ti,2)
          u_ISMIP_HOM( ti) = umin + (umax - umin) * ( (1._dp - sin( x)) / 2._dp)**2
        end do

      end if

      ! Initialise velocity fields with the approximation
      if     (C%choice_ice_dynamics == 'SIA/SSA') then
        ice%u_base_SSA_b( mesh%ti1:mesh%ti2) = u_ISMIP_HOM( mesh%ti1:mesh%ti2)
      elseif (C%choice_ice_dynamics == 'DIVA') then
        ice%u_vav_b(      mesh%ti1:mesh%ti2) = u_ISMIP_HOM( mesh%ti1:mesh%ti2)
      end if

      ! Clean up after yourself
      deallocate( u_ISMIP_HOM)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = huge( 1))

  end subroutine initialise_velocity_solver

  subroutine initialise_matrix_conversion_lists( mesh, ice)
    ! Initialise lists for converting triangle indices to stiffness matrix rows and vice versa

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(LEN=256), parameter       :: routine_name = 'initialise_matrix_conversion_lists'
    integer                             :: ti, nu, nv

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2

      nu = 2*ti - 1
      nv = 2*ti

      ice%ti2n_u( ti) = nu
      ice%ti2n_v( ti) = nv

      ice%n2ti_uv( nu,:) = [ti,0 ]
      ice%n2ti_uv( nv,:) = [0 ,ti]

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_matrix_conversion_lists

  subroutine initialise_SSADIVA_stiffness_matrix( mesh, ice)
    ! Initialise the non-zero structure template of the SSA/DIVA stiffness matrix

    implicit none

    ! In- and output variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_SSADIVA_stiffness_matrix'
    integer                             :: ncols, nrows, nnz_per_row_max, nnz_max
    integer                             :: n1, n2, n, ti, k1, k2, k, tj, mu, mv
    integer, dimension(:), allocatable  :: ti2n_u, ti2n_v

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory for A
    ncols           = 2*mesh%nTri    ! from
    nrows           = 2*mesh%nTri    ! to
    nnz_per_row_max = 20

    nnz_max = nrows * nnz_per_row_max
    call allocate_matrix_CSR_dist( ice%M_SSADIVA, nrows, ncols, nnz_max)

    ! Fill in matrix rows, right-hand side, and initial guess
    call partition_list( mesh%nTri, par%i, par%n, n1, n2)

    allocate(ti2n_u( 1:mesh%nTri))
    allocate(ti2n_v( 1:mesh%nTri))
    ti2n_u(mesh%ti1:mesh%ti2) = ice%ti2n_u
    ti2n_v(mesh%ti1:mesh%ti2) = ice%ti2n_v
    call allgather_array(ti2n_u)
    call allgather_array(ti2n_v)

    do n = 1+2*(n1-1), n2*2

      ! Fill matrix coefficients
      if (ice%n2ti_uv( n,1) > 0) then
        ! u

        ti = ice%n2ti_uv( n,1)

        if (mesh%Tri_edge_index( ti) == 0) then
          ! Free triangle: fill in matrix row for the SSA/DIVA

          k1 = mesh%M2_ddx_b_b_CSR%ptr( ti)
          k2 = mesh%M2_ddx_b_b_CSR%ptr( ti+1) - 1

          do k = k1, k2

            tj  = mesh%M2_ddx_b_b_CSR%index( k)
            mu = ti2n_u( tj)
            mv = ti2n_v( tj)

            ! u-part
            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mu

            ! v-part
            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mv

          end do

        else ! IF (mesh%Tri_edge_index( ti) == 0) THEN
          ! Border triangle: apply boundary conditions

          k1 = mesh%M_Neumann_BC_b_CSR%ptr( ti)
          k2 = mesh%M_Neumann_BC_b_CSR%ptr( ti+1) - 1

          ! Matrix
          do k = mesh%M_Neumann_BC_b_CSR%ptr( ti), mesh%M_Neumann_BC_b_CSR%ptr( ti+1) - 1

            tj  = mesh%M_Neumann_BC_b_CSR%index( k)
            mu = ti2n_u( tj)
            mv = ti2n_v( tj)

            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mu

          end do

        end if ! IF (mesh%Tri_edge_index( ti) == 0) THEN

      else ! IF (MOD( n,2) == 1) THEN
        ! v

        ti = ice%n2ti_uv( n,2)

        if (mesh%Tri_edge_index( ti) == 0) then
          ! Free triangle: fill in matrix row for the SSA/DIVA

          k1 = mesh%M2_ddx_b_b_CSR%ptr( ti)
          k2 = mesh%M2_ddx_b_b_CSR%ptr( ti+1) - 1

          do k = k1, k2

            tj  = mesh%M2_ddx_b_b_CSR%index( k)
            mu = ti2n_u( tj)
            mv = ti2n_v( tj)

            ! u-part
            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mu

            ! v-part
            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mv

          end do

        else ! IF (mesh%Tri_edge_index( ti) == 0) THEN
          ! Border triangle: apply boundary conditions

          k1 = mesh%M_Neumann_BC_b_CSR%ptr( ti)
          k2 = mesh%M_Neumann_BC_b_CSR%ptr( ti+1) - 1

          ! Matrix
          do k = mesh%M_Neumann_BC_b_CSR%ptr( ti), mesh%M_Neumann_BC_b_CSR%ptr( ti+1) - 1

            tj  = mesh%M_Neumann_BC_b_CSR%index( k)
            mu = ti2n_u( tj)
            mv = ti2n_v( tj)

            ice%M_SSADIVA%nnz =  ice%M_SSADIVA%nnz + 1
            ice%M_SSADIVA%index( ice%M_SSADIVA%nnz) = mv

          end do

        end if ! IF (mesh%Tri_edge_index( ti) == 0) THEN

      end if ! IF (MOD( n,2) == 1) THEN

      ! Finalise this matrix row
      ice%M_SSADIVA%ptr( n+1 : 2*mesh%nTri+1) = ice%M_SSADIVA%nnz+1

    end do ! DO n = n1, n2

    ! Combine results from the different processes
    call finalise_matrix_CSR_dist( ice%M_SSADIVA, 1+2*(n1-1), n2*2)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = 7)

  end subroutine initialise_SSADIVA_stiffness_matrix

  SUBROUTINE remap_velocities( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_velocities'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_dynamics == 'SIA') THEN
      call crash('todo, implement SIA')
      !CALL remap_velocities_SIA( mesh_old, mesh_new, map, ice)
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      call crash('todo, implement SSA')
      !CALL remap_velocities_SSA( mesh_old, mesh_new, map, ice)
    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      call crash('todo, implement SIA/SSA')
      !CALL remap_velocities_SIASSA( mesh_old, mesh_new, map, ice)
    ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      CALL remap_velocities_DIVA( mesh_old, mesh_new, map, ice)
    ELSE
      CALL crash('unknown choice_ice_dynamics "' // TRIM( C%choice_ice_dynamics) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_velocities

#if 0

  SUBROUTINE remap_velocities_SIA( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_velocities_SIA'
    REAL(dp)                                           :: dp_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dp_dummy = mesh_old%V(   1,1)
    dp_dummy = mesh_new%V(   1,1)
    dp_dummy = map%int_dummy
    dp_dummy = ice%Hi_a( 1)

    ! No need to remap anything, just reallocate
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%u_3D_a                , ice%wu_3D_a               )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%v_3D_a                , ice%wv_3D_a               )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%u_3D_b                , ice%wu_3D_b               )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%v_3D_b                , ice%wv_3D_b               )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%w_3D_a                , ice%ww_3D_a               )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_vav_a               , ice%wu_vav_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_vav_a               , ice%wv_vav_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_vav_b               , ice%wu_vav_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_vav_b               , ice%wv_vav_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_vav_a            , ice%wuabs_vav_a           )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_vav_b            , ice%wuabs_vav_b           )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_surf_a              , ice%wu_surf_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_surf_a              , ice%wv_surf_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_surf_b              , ice%wu_surf_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_surf_b              , ice%wv_surf_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_surf_a           , ice%wuabs_surf_a          )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_surf_b           , ice%wuabs_surf_b          )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_base_a              , ice%wu_base_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_base_a              , ice%wv_base_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_base_b              , ice%wu_base_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_base_b              , ice%wv_base_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_base_a           , ice%wuabs_base_a          )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_base_b           , ice%wuabs_base_b          )

    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%u_3D_SIA_b            , ice%wu_3D_SIA_b           )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%v_3D_SIA_b            , ice%wv_3D_SIA_b           )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_base_SSA_b          , ice%wu_base_SSA_b         )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_base_SSA_b          , ice%wv_base_SSA_b         )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_velocities_SIA
  SUBROUTINE remap_velocities_SSA( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_velocities_SSA'
    REAL(dp), DIMENSION(:    ), POINTER                ::  u_a,  v_a
    INTEGER                                            :: wu_a, wv_a

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Remap SSA velocities
  ! =======================

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh_old%nV, u_a, wu_a)
    CALL allocate_shared_dp_1D( mesh_old%nV, v_a, wv_a)

    ! Map velocities to the a-grid
    CALL map_b_to_a_2D( mesh_old, ice%u_base_SSA_b, u_a)
    CALL map_b_to_a_2D( mesh_old, ice%v_base_SSA_b, v_a)

    ! Remap a-grid velocities
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, u_a, wu_a, 'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, v_a, wv_a, 'cons_2nd_order')

    ! Reallocate b-grid velocities
    CALL reallocate_shared_dp_1D( mesh_new%nTri, ice%u_base_SSA_b, ice%wu_base_SSA_b)
    CALL reallocate_shared_dp_1D( mesh_new%nTri, ice%v_base_SSA_b, ice%wv_base_SSA_b)

    ! Map remapped velocities to the b-grid
    CALL map_a_to_b_2D( mesh_new, u_a, ice%u_base_SSA_b)
    CALL map_a_to_b_2D( mesh_new, u_a, ice%v_base_SSA_b)

    ! Clean up after yourself
    CALL deallocate_shared( wu_a)
    CALL deallocate_shared( wv_a)

  ! == Reallocate everything else
  ! =============================

    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%u_3D_a                , ice%wu_3D_a               )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%v_3D_a                , ice%wv_3D_a               )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%u_3D_b                , ice%wu_3D_b               )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%v_3D_b                , ice%wv_3D_b               )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%w_3D_a                , ice%ww_3D_a               )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_vav_a               , ice%wu_vav_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_vav_a               , ice%wv_vav_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_vav_b               , ice%wu_vav_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_vav_b               , ice%wv_vav_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_vav_a            , ice%wuabs_vav_a           )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_vav_b            , ice%wuabs_vav_b           )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_surf_a              , ice%wu_surf_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_surf_a              , ice%wv_surf_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_surf_b              , ice%wu_surf_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_surf_b              , ice%wv_surf_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_surf_a           , ice%wuabs_surf_a          )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_surf_b           , ice%wuabs_surf_b          )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_base_a              , ice%wu_base_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_base_a              , ice%wv_base_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_base_b              , ice%wu_base_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_base_b              , ice%wv_base_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_base_a           , ice%wuabs_base_a          )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_base_b           , ice%wuabs_base_b          )

   !CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%u_3D_SIA_b            , ice%wu_3D_SIA_b           )
   !CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%v_3D_SIA_b            , ice%wv_3D_SIA_b           )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_base_SSA_b          , ice%wu_base_SSA_b         )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_base_SSA_b          , ice%wv_base_SSA_b         )

    ! Physical terms in the SSA/DIVA
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%taudx_b               , ice%wtaudx_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%taudy_b               , ice%wtaudy_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%du_dx_a               , ice%wdu_dx_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%du_dy_a               , ice%wdu_dy_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dv_dx_a               , ice%wdv_dx_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dv_dy_a               , ice%wdv_dy_a              )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz,            ice%du_dz_3D_b            , ice%wdu_dz_3D_b           )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz,            ice%dv_dz_3D_b            , ice%wdv_dz_3D_b           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz,            ice%visc_eff_3D_a         , ice%wvisc_eff_3D_a        )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%visc_eff_int_a        , ice%wvisc_eff_int_a       )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%N_a                   , ice%wN_a                  )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%beta_a                , ice%wbeta_a               )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%beta_eff_a            , ice%wbeta_eff_a           )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%beta_eff_b            , ice%wbeta_eff_b           )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%taubx_b               , ice%wtaubx_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%tauby_b               , ice%wtauby_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%F2_a                  , ice%wF2_a                 )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_prev_b              , ice%wu_prev_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_prev_b              , ice%wv_prev_b             )

    ! Some administrative stuff to make solving the SSA/DIVA more efficient
    CALL reallocate_shared_int_1D(  mesh_new%nTri,                  ice%ti2n_u                , ice%wti2n_u               )
    CALL reallocate_shared_int_1D(  mesh_new%nTri,                  ice%ti2n_v                , ice%wti2n_v               )
    CALL reallocate_shared_int_2D(2*mesh_new%nTri, 2              , ice%n2ti_uv               , ice%wn2ti_uv              )
    CALL deallocate_matrix_CSR( ice%M_SSADIVA)

    CALL initialise_matrix_conversion_lists(  mesh_new, ice)
    CALL initialise_SSADIVA_stiffness_matrix( mesh_new, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_velocities_SSA
  SUBROUTINE remap_velocities_SIASSA( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_velocities_SIASSA'
    REAL(dp), DIMENSION(:    ), POINTER                ::  u_a,  v_a
    INTEGER                                            :: wu_a, wv_a

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Remap SSA velocities
  ! =======================

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh_old%nV, u_a, wu_a)
    CALL allocate_shared_dp_1D( mesh_old%nV, v_a, wv_a)

    ! Map velocities to the a-grid
    CALL map_b_to_a_2D( mesh_old, ice%u_base_SSA_b, u_a)
    CALL map_b_to_a_2D( mesh_old, ice%v_base_SSA_b, v_a)

    ! Remap a-grid velocities
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, u_a, wu_a, 'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, v_a, wv_a, 'cons_2nd_order')

    ! Reallocate b-grid velocities
    CALL reallocate_shared_dp_1D( mesh_new%nTri, ice%u_base_SSA_b, ice%wu_base_SSA_b)
    CALL reallocate_shared_dp_1D( mesh_new%nTri, ice%v_base_SSA_b, ice%wv_base_SSA_b)

    ! Map remapped velocities to the b-grid
    CALL map_a_to_b_2D( mesh_new, u_a, ice%u_base_SSA_b)
    CALL map_a_to_b_2D( mesh_new, u_a, ice%v_base_SSA_b)

    ! Clean up after yourself
    CALL deallocate_shared( wu_a)
    CALL deallocate_shared( wv_a)

  ! == Reallocate everything else
  ! =============================

    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%u_3D_a                , ice%wu_3D_a               )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%v_3D_a                , ice%wv_3D_a               )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%u_3D_b                , ice%wu_3D_b               )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%v_3D_b                , ice%wv_3D_b               )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%w_3D_a                , ice%ww_3D_a               )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_vav_a               , ice%wu_vav_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_vav_a               , ice%wv_vav_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_vav_b               , ice%wu_vav_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_vav_b               , ice%wv_vav_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_vav_a            , ice%wuabs_vav_a           )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_vav_b            , ice%wuabs_vav_b           )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_surf_a              , ice%wu_surf_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_surf_a              , ice%wv_surf_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_surf_b              , ice%wu_surf_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_surf_b              , ice%wv_surf_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_surf_a           , ice%wuabs_surf_a          )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_surf_b           , ice%wuabs_surf_b          )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_base_a              , ice%wu_base_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_base_a              , ice%wv_base_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_base_b              , ice%wu_base_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_base_b              , ice%wv_base_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_base_a           , ice%wuabs_base_a          )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_base_b           , ice%wuabs_base_b          )

    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%u_3D_SIA_b            , ice%wu_3D_SIA_b           )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%v_3D_SIA_b            , ice%wv_3D_SIA_b           )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_base_SSA_b          , ice%wu_base_SSA_b         )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_base_SSA_b          , ice%wv_base_SSA_b         )

    ! Physical terms in the SSA/DIVA
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%taudx_b               , ice%wtaudx_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%taudy_b               , ice%wtaudy_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%du_dx_a               , ice%wdu_dx_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%du_dy_a               , ice%wdu_dy_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dv_dx_a               , ice%wdv_dx_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dv_dy_a               , ice%wdv_dy_a              )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz,            ice%du_dz_3D_b            , ice%wdu_dz_3D_b           )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz,            ice%dv_dz_3D_b            , ice%wdv_dz_3D_b           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz,            ice%visc_eff_3D_a         , ice%wvisc_eff_3D_a        )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%visc_eff_int_a        , ice%wvisc_eff_int_a       )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%N_a                   , ice%wN_a                  )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%beta_a                , ice%wbeta_a               )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%beta_eff_a            , ice%wbeta_eff_a           )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%beta_eff_b            , ice%wbeta_eff_b           )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%taubx_b               , ice%wtaubx_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%tauby_b               , ice%wtauby_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%F2_a                  , ice%wF2_a                 )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_prev_b              , ice%wu_prev_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_prev_b              , ice%wv_prev_b             )

    ! Some administrative stuff to make solving the SSA/DIVA more efficient
    CALL reallocate_shared_int_1D(  mesh_new%nTri,                  ice%ti2n_u                , ice%wti2n_u               )
    CALL reallocate_shared_int_1D(  mesh_new%nTri,                  ice%ti2n_v                , ice%wti2n_v               )
    CALL reallocate_shared_int_2D(2*mesh_new%nTri, 2              , ice%n2ti_uv               , ice%wn2ti_uv              )
    CALL deallocate_matrix_CSR( ice%M_SSADIVA)

    CALL initialise_matrix_conversion_lists(  mesh_new, ice)
    CALL initialise_SSADIVA_stiffness_matrix( mesh_new, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_velocities_SIASSA
#endif
  SUBROUTINE remap_velocities_DIVA( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_velocities_DIVA'
    REAL(dp), DIMENSION(:    ), allocatable            ::  u_a,  v_a

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Remap SSA velocities
  ! =======================

    ! Allocate shared memory
    allocate( u_a( mesh_old%vi1:mesh_old%vi2))
    allocate( v_a( mesh_old%vi1:mesh_old%vi2))

    ! Map velocities to the a-grid
    CALL map_b_to_a_2D( mesh_old, ice%u_vav_b, u_a)
    CALL map_b_to_a_2D( mesh_old, ice%v_vav_b, v_a)

    ! Remap a-grid velocities
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, u_a, 'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, v_a, 'cons_2nd_order')

    ! Reallocate b-grid velocities
    CALL reallocate_bounds( ice%u_vav_b, mesh_new%ti1, mesh_new%ti2)
    CALL reallocate_bounds( ice%v_vav_b, mesh_new%ti1, mesh_new%ti2)

    ! Map remapped velocities to the b-grid
    CALL map_a_to_b_2D( mesh_new, u_a, ice%u_vav_b)
    CALL map_a_to_b_2D( mesh_new, u_a, ice%v_vav_b)

    ! Clean up after yourself
    deallocate( u_a)
    deallocate( v_a)

  ! == Reallocate everything else
  ! =============================

    call reallocate_bounds( ice%u_3D_a        , mesh_new%vi1, mesh_new%vi2  , C%nz )
    call reallocate_bounds( ice%v_3D_a        , mesh_new%vi1, mesh_new%vi2  , C%nz )
    call reallocate_bounds( ice%u_3D_b        , mesh_new%ti1, mesh_new%ti2  , C%nz )
    call reallocate_bounds( ice%v_3D_b        , mesh_new%ti1, mesh_new%ti2  , C%nz )
    call reallocate_bounds( ice%w_3D_a        , mesh_new%vi1, mesh_new%vi2  , C%nz )

    call reallocate_bounds( ice%u_vav_a       , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%v_vav_a       , mesh_new%vi1, mesh_new%vi2         )
   !call reallocate_bounds( ice%u_vav_b       , mesh_new%ti1, mesh_new%ti2         )
   !call reallocate_bounds( ice%v_vav_b       , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%uabs_vav_a    , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%uabs_vav_b    , mesh_new%ti1, mesh_new%ti2         )

    call reallocate_bounds( ice%u_surf_a      , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%v_surf_a      , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%u_surf_b      , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%v_surf_b      , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%uabs_surf_a   , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%uabs_surf_b   , mesh_new%ti1, mesh_new%ti2         )

    call reallocate_bounds( ice%u_base_a      , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%v_base_a      , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%u_base_b      , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%v_base_b      , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%uabs_base_a   , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%uabs_base_b   , mesh_new%ti1, mesh_new%ti2         )

   !call reallocate_bounds( ice%u_3D_SIA_b    , mesh_new%ti1, mesh_new%ti2  , C%nz )
   !call reallocate_bounds( ice%v_3D_SIA_b    , mesh_new%ti1, mesh_new%ti2  , C%nz )
   !call reallocate_bounds( ice%u_base_SSA_b  , mesh_new%ti1, mesh_new%ti2         )
   !call reallocate_bounds( ice%v_base_SSA_b  , mesh_new%ti1, mesh_new%ti2         )

    ! Physical terms in the SSA/DIVA
    call reallocate_bounds( ice%taudx_b       , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%taudy_b       , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%du_dx_a       , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%du_dy_a       , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%dv_dx_a       , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%dv_dy_a       , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%du_dz_3D_b    , mesh_new%ti1, mesh_new%ti2  , C%nz )
    call reallocate_bounds( ice%dv_dz_3D_b    , mesh_new%ti1, mesh_new%ti2  , C%nz )
    call reallocate_bounds( ice%visc_eff_3D_a , mesh_new%vi1, mesh_new%vi2  , C%nz )
    call reallocate_bounds( ice%visc_eff_int_a, mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%N_a           , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%beta_a        , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%beta_eff_a    , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%beta_eff_b    , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%taubx_b       , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%tauby_b       , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%F2_a          , mesh_new%vi1, mesh_new%vi2         )
    call reallocate_bounds( ice%u_prev_b      , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%v_prev_b      , mesh_new%ti1, mesh_new%ti2         )


    ! Some administrative stuff to make solving the SSA/DIVA more efficient

    call reallocate_bounds( ice%ti2n_u        , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%ti2n_v        , mesh_new%ti1, mesh_new%ti2         )
    call reallocate_bounds( ice%n2ti_uv       ,(mesh_new%ti1-1)*2+1, mesh_new%ti2*2, 2  )
    CALL deallocate_matrix_CSR( ice%M_SSADIVA)


    CALL initialise_matrix_conversion_lists(  mesh_new, ice)
    CALL initialise_SSADIVA_stiffness_matrix( mesh_new, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_velocities_DIVA

END MODULE ice_velocity_module
