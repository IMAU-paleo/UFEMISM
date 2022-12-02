MODULE ice_velocity_DIVA_module

! == Contains all the routines needed to solve the Depth-Integrated Viscosity Approximation

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE petscksp
  USE configuration_module,                ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                        ONLY: perr, petscmat_checksum
  USE parallel_module,                     ONLY: par, sync, ierr, cerr, partition_list, &
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
  USE utilities_module,                    ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                                 check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                                 checksum_dp_1D, checksum_dp_2D, checksum_dp_3D
  USE netcdf_debug_module,                 ONLY: debug, write_to_debug_file, &
                                                 save_variable_as_netcdf_int_1D, save_variable_as_netcdf_dp_1D, &
                                                 save_variable_as_netcdf_int_2D, save_variable_as_netcdf_dp_2D, &
                                                 save_variable_as_netcdf_int_3D, save_variable_as_netcdf_dp_3D, &
                                                 write_PETSc_matrix_to_NetCDF, write_CSR_matrix_to_NetCDF

  ! Import specific functionality
  USE data_types_module,                   ONLY: type_mesh, type_ice_model, type_velocity_solver_DIVA
  USE petsc_module,                        ONLY: vec_double2petsc, multiply_PETSc_matrix_with_vector_1D, double2petscmat, solve_matrix_equation_PETSc
  USE mesh_operators_module,               ONLY: map_a_to_b_2D, map_a_to_b_3D, ddx_a_to_b_2D, ddy_a_to_b_2D, &
                                                 map_b_to_a_2D, map_b_to_a_3D, ddx_b_to_a_2D, ddy_b_to_a_2D
  USE utilities_module,                    ONLY: vertical_integrate, vertical_integration_from_bottom_to_zeta
  USe basal_conditions_and_sliding_module, ONLY: calc_basal_conditions, calc_basal_friction_coefficient
  USE mesh_help_functions_module,          ONLY: find_ti_copy_ISMIP_HOM_periodic
  USE general_ice_model_data_module,       ONLY: determine_grounded_fractions

  IMPLICIT NONE

CONTAINS

! == The main routines, to be called from the ice_velocity_module

  SUBROUTINE initialise_DIVA_solver( mesh, DIVA)
    ! Initialise the DIVA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'initialise_DIVA_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory

    ! Solution
    CALL allocate_shared_dp_2D( mesh%nTri   , C%nz, DIVA%u_3D_b           , DIVA%wu_3D_b           )
    CALL allocate_shared_dp_2D( mesh%nTri   , C%nz, DIVA%v_3D_b           , DIVA%wv_3D_b           )

    ! Intermediate data fields (all in vector form)
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%u_vav_b          , DIVA%wu_vav_b          )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%v_vav_b          , DIVA%wv_vav_b          )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%u_base_b         , DIVA%wu_base_b         )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%v_base_b         , DIVA%wv_base_b         )
    CALL allocate_shared_dp_1D( mesh%nV     ,       DIVA%du_dx_a          , DIVA%wdu_dx_a          )
    CALL allocate_shared_dp_1D( mesh%nV     ,       DIVA%du_dy_a          , DIVA%wdu_dy_a          )
    CALL allocate_shared_dp_1D( mesh%nV     ,       DIVA%dv_dx_a          , DIVA%wdv_dx_a          )
    CALL allocate_shared_dp_1D( mesh%nV     ,       DIVA%dv_dy_a          , DIVA%wdv_dy_a          )
    CALL allocate_shared_dp_2D( mesh%nTri   , C%nz, DIVA%du_dz_3D_b       , DIVA%wdu_dz_3D_b       )
    CALL allocate_shared_dp_2D( mesh%nTri   , C%nz, DIVA%dv_dz_3D_b       , DIVA%wdv_dz_3D_b       )
    CALL allocate_shared_dp_2D( mesh%nV     , C%nz, DIVA%du_dz_3D_a       , DIVA%wdu_dz_3D_a       )
    CALL allocate_shared_dp_2D( mesh%nV     , C%nz, DIVA%dv_dz_3D_a       , DIVA%wdv_dz_3D_a       )
    CALL allocate_shared_dp_2D( mesh%nV     , C%nz, DIVA%eta_3D_a         , DIVA%weta_3D_a         )
    CALL allocate_shared_dp_2D( mesh%nTri   , C%nz, DIVA%eta_3D_b         , DIVA%weta_3D_b         )
    CALL allocate_shared_dp_1D( mesh%nV     ,       DIVA%eta_vint_a       , DIVA%weta_vint_a       )
    CALL allocate_shared_dp_1D( mesh%nV     ,       DIVA%N_a              , DIVA%wN_a              )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%N_b              , DIVA%wN_b              )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%dN_dx_b          , DIVA%wdN_dx_b          )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%dN_dy_b          , DIVA%wdN_dy_b          )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%beta_b_b         , DIVA%wbeta_b_b         )
    CALL allocate_shared_dp_1D( mesh%nV     ,       DIVA%F2_a             , DIVA%wF2_a             )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%F2_b             , DIVA%wF2_b             )
    CALL allocate_shared_dp_1D( mesh%nV     ,       DIVA%beta_eff_a       , DIVA%wbeta_eff_a       )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%beta_eff_b       , DIVA%wbeta_eff_b       )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%taubx_b          , DIVA%wtaubx_b          )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%tauby_b          , DIVA%wtauby_b          )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%taudx_b          , DIVA%wtaudx_b          )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%taudy_b          , DIVA%wtaudy_b          )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%u_vav_b_prev     , DIVA%wu_vav_b_prev     )
    CALL allocate_shared_dp_1D( mesh%nTri   ,       DIVA%v_vav_b_prev     , DIVA%wv_vav_b_prev     )

    ! Parameters for the iterative solver used to solve the matrix equation representing the linearised DIVA
    DIVA%PETSc_rtol      = C%DIVA_PETSc_rtol
    DIVA%PETSc_abstol    = C%DIVA_PETSc_abstol

    ! Load vectors
    CALL allocate_shared_dp_1D( mesh%nnbuv  , DIVA%bb               , DIVA%wbb               )
    CALL allocate_shared_dp_1D( mesh%nnbuv  , DIVA%bb_free          , DIVA%wbb_free          )
    CALL allocate_shared_dp_1D( mesh%nnbuv  , DIVA%bb_BC_west       , DIVA%wbb_BC_west       )
    CALL allocate_shared_dp_1D( mesh%nnbuv  , DIVA%bb_BC_east       , DIVA%wbb_BC_east       )
    CALL allocate_shared_dp_1D( mesh%nnbuv  , DIVA%bb_BC_south      , DIVA%wbb_BC_south      )
    CALL allocate_shared_dp_1D( mesh%nnbuv  , DIVA%bb_BC_north      , DIVA%wbb_BC_north      )
    CALL allocate_shared_dp_1D( mesh%nnbuv  , DIVA%bb_BC_prescr     , DIVA%wbb_BC_prescr     )

    ! Calculate combined mesh operators for efficient calculation of the stiffness matrix
    CALL calc_combined_mesh_operators_DIVA( mesh, DIVA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_DIVA_solver

  SUBROUTINE solve_DIVA( mesh, ice, DIVA, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    ! Calculate ice velocities by solving the Shallow Ice Approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA
    INTEGER,  DIMENSION(:    ),          INTENT(IN)   , OPTIONAL :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(:    ),          INTENT(IN)   , OPTIONAL :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:    ),          INTENT(IN)   , OPTIONAL :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_DIVA'
    INTEGER,  DIMENSION(:    ), POINTER                          :: BC_prescr_mask_b_applied
    REAL(dp), DIMENSION(:    ), POINTER                          :: BC_prescr_u_b_applied
    REAL(dp), DIMENSION(:    ), POINTER                          :: BC_prescr_v_b_applied
    INTEGER                                                      :: wBC_prescr_mask_b_applied, wBC_prescr_u_b_applied, wBC_prescr_v_b_applied
    INTEGER                                                      :: viscosity_iteration_i
    LOGICAL                                                      :: has_converged
    REAL(dp)                                                     :: resid_UV

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If there is no grounded ice, no way to solve the DIVA
    IF (SUM( ice%mask_sheet_a) == 0) THEN
      DIVA%u_3D_b = 0._dp
      DIVA%v_3D_b = 0._dp
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Handle the optional prescribed u,v boundary conditions
    CALL allocate_shared_int_1D( mesh%nTri, BC_prescr_mask_b_applied, wBC_prescr_mask_b_applied)
    CALL allocate_shared_dp_1D(  mesh%nTri, BC_prescr_u_b_applied   , wBC_prescr_u_b_applied   )
    CALL allocate_shared_dp_1D(  mesh%nTri, BC_prescr_v_b_applied   , wBC_prescr_v_b_applied   )
    IF (PRESENT( BC_prescr_mask_b) .OR. PRESENT( BC_prescr_u_b) .OR. PRESENT( BC_prescr_v_b)) THEN
      ! Safety
      IF (.NOT. (PRESENT(BC_prescr_mask_b) .AND. PRESENT(BC_prescr_u_b) .AND. PRESENT(BC_prescr_v_b))) THEN
        CALL crash('need to provide prescribed u,v fields and mask!')
      END IF
      BC_prescr_mask_b_applied( mesh%ti1:mesh%ti2) = BC_prescr_mask_b( mesh%ti1:mesh%ti2)
      BC_prescr_u_b_applied(    mesh%ti1:mesh%ti2) = BC_prescr_u_b(    mesh%ti1:mesh%ti2)
      BC_prescr_v_b_applied(    mesh%ti1:mesh%ti2) = BC_prescr_v_b(    mesh%ti1:mesh%ti2)
    ELSE
      BC_prescr_mask_b_applied( mesh%ti1:mesh%ti2) = 0
      BC_prescr_u_b_applied(    mesh%ti1:mesh%ti2) = 0._dp
      BC_prescr_v_b_applied(    mesh%ti1:mesh%ti2) = 0._dp
    END IF

    ! Calculate boundary conditions mask matrices
    CALL calc_DIVA_BC_mask_matrices( mesh, DIVA, BC_prescr_mask_b_applied)

    ! Calculate the driving stress
    CALL calc_driving_stress( mesh, ice, DIVA)

    ! Calculate the basal yield stress tau_c
    CALL calc_basal_conditions( mesh, ice)

    ! Determine sub-mesh grounded fractions for scaling the basal friction
    CALL determine_grounded_fractions( mesh, ice)

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .FALSE.
    viscosity_iteration: DO WHILE (.NOT. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! Calculate the strain rates for the current velocity solution
      CALL calc_strain_rates( mesh, ice, DIVA)

      ! Calculate the effective viscosity for the current velocity solution
      CALL calc_effective_viscosity( mesh, ice, DIVA)

      ! Calculate the sliding + vertical shearing term beta_eff
      CALL calc_beta_eff( mesh, ice, DIVA)

      ! Solve the linearised DIVA to calculate a new velocity solution
      CALL solve_DIVA_linearised( mesh, DIVA, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Limit velocities for improved stability
      CALL apply_velocity_limits( mesh, DIVA)

      ! Reduce the change between velocity solutions
      CALL relax_viscosity_iterations( mesh, DIVA, C%DIVA_visc_it_relax)

      ! Calculate basal shear stress
      CALL calc_basal_shear_stress( mesh, DIVA)

      ! Calculate basal velocities
      CALL calc_basal_velocities( mesh, ice, DIVA)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      CALL calc_visc_iter_UV_resid( mesh, DIVA, resid_UV)

!      ! DENK DROM
!      IF (par%master) WRITE(0,*) '    DIVA - viscosity iteration ', viscosity_iteration_i, ', u = [', MINVAL( DIVA%u_vav_b), ' - ', MAXVAL( DIVA%u_vav_b), '], resid = ', resid_UV

      ! If the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .FALSE.
      IF     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) THEN
        has_converged = .TRUE.
      ELSEIF (viscosity_iteration_i > C%DIVA_visc_it_nit) THEN
        has_converged = .TRUE.
      END IF

    END DO viscosity_iteration

    ! Calculate the 3D horizontal velocity fields
    CALL calc_3D_horizontal_velocities( mesh, ice, DIVA)

    ! Clean up after yourself
    CALL deallocate_shared( wBC_prescr_mask_b_applied)
    CALL deallocate_shared( wBC_prescr_u_b_applied   )
    CALL deallocate_shared( wBC_prescr_v_b_applied   )
    CALL MatDestroy( DIVA%m_free     , perr)
    CALL MatDestroy( DIVA%m_BC_west  , perr)
    CALL MatDestroy( DIVA%m_BC_east  , perr)
    CALL MatDestroy( DIVA%m_BC_south , perr)
    CALL MatDestroy( DIVA%m_BC_north , perr)
    CALL MatDestroy( DIVA%m_BC_prescr, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_DIVA

  SUBROUTINE remap_DIVA_solver( mesh_old, mesh_new, ice, DIVA)
    ! Remap the DIVA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh_old
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh_new
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'remap_DIVA_solver'
    REAL(dp)                                                     :: dp_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dp_dummy = mesh_old%V( 1,1)
    dp_dummy = mesh_new%V( 1,1)
    dp_dummy = ice%Hi_a( 1)

    ! Reallocate shared memory

    ! Solution
    CALL reallocate_shared_dp_2D( mesh_new%nTri   , C%nz, DIVA%u_3D_b           , DIVA%wu_3D_b           )
    CALL reallocate_shared_dp_2D( mesh_new%nTri   , C%nz, DIVA%v_3D_b           , DIVA%wv_3D_b           )

    ! Intermediate data fields (all in vector form)
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%u_vav_b          , DIVA%wu_vav_b          )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%v_vav_b          , DIVA%wv_vav_b          )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%u_base_b         , DIVA%wu_base_b         )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%v_base_b         , DIVA%wv_base_b         )
    CALL reallocate_shared_dp_1D( mesh_new%nV     ,       DIVA%du_dx_a          , DIVA%wdu_dx_a          )
    CALL reallocate_shared_dp_1D( mesh_new%nV     ,       DIVA%du_dy_a          , DIVA%wdu_dy_a          )
    CALL reallocate_shared_dp_1D( mesh_new%nV     ,       DIVA%dv_dx_a          , DIVA%wdv_dx_a          )
    CALL reallocate_shared_dp_1D( mesh_new%nV     ,       DIVA%dv_dy_a          , DIVA%wdv_dy_a          )
    CALL reallocate_shared_dp_2D( mesh_new%nTri   , C%nz, DIVA%du_dz_3D_b       , DIVA%wdu_dz_3D_b       )
    CALL reallocate_shared_dp_2D( mesh_new%nTri   , C%nz, DIVA%dv_dz_3D_b       , DIVA%wdv_dz_3D_b       )
    CALL reallocate_shared_dp_2D( mesh_new%nV     , C%nz, DIVA%du_dz_3D_a       , DIVA%wdu_dz_3D_a       )
    CALL reallocate_shared_dp_2D( mesh_new%nV     , C%nz, DIVA%dv_dz_3D_a       , DIVA%wdv_dz_3D_a       )
    CALL reallocate_shared_dp_2D( mesh_new%nV     , C%nz, DIVA%eta_3D_a         , DIVA%weta_3D_a         )
    CALL reallocate_shared_dp_2D( mesh_new%nTri   , C%nz, DIVA%eta_3D_b         , DIVA%weta_3D_b         )
    CALL reallocate_shared_dp_1D( mesh_new%nV     ,       DIVA%eta_vint_a       , DIVA%weta_vint_a       )
    CALL reallocate_shared_dp_1D( mesh_new%nV     ,       DIVA%N_a              , DIVA%wN_a              )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%N_b              , DIVA%wN_b              )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%dN_dx_b          , DIVA%wdN_dx_b          )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%dN_dy_b          , DIVA%wdN_dy_b          )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%beta_b_b         , DIVA%wbeta_b_b         )
    CALL reallocate_shared_dp_1D( mesh_new%nV     ,       DIVA%F2_a             , DIVA%wF2_a             )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%F2_b             , DIVA%wF2_b             )
    CALL reallocate_shared_dp_1D( mesh_new%nV     ,       DIVA%beta_eff_a       , DIVA%wbeta_eff_a       )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%beta_eff_b       , DIVA%wbeta_eff_b       )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%taubx_b          , DIVA%wtaubx_b          )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%tauby_b          , DIVA%wtauby_b          )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%taudx_b          , DIVA%wtaudx_b          )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%taudy_b          , DIVA%wtaudy_b          )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%u_vav_b_prev     , DIVA%wu_vav_b_prev     )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,       DIVA%v_vav_b_prev     , DIVA%wv_vav_b_prev     )

    ! Parameters for the iterative solver used to solve the matrix equation representing the linearised DIVA
    DIVA%PETSc_rtol      = C%DIVA_PETSc_rtol
    DIVA%PETSc_abstol    = C%DIVA_PETSc_abstol

    ! Load vectors
    CALL reallocate_shared_dp_1D( mesh_new%nnbuv  , DIVA%bb               , DIVA%wbb               )
    CALL reallocate_shared_dp_1D( mesh_new%nnbuv  , DIVA%bb_free          , DIVA%wbb_free          )
    CALL reallocate_shared_dp_1D( mesh_new%nnbuv  , DIVA%bb_BC_west       , DIVA%wbb_BC_west       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbuv  , DIVA%bb_BC_east       , DIVA%wbb_BC_east       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbuv  , DIVA%bb_BC_south      , DIVA%wbb_BC_south      )
    CALL reallocate_shared_dp_1D( mesh_new%nnbuv  , DIVA%bb_BC_north      , DIVA%wbb_BC_north      )
    CALL reallocate_shared_dp_1D( mesh_new%nnbuv  , DIVA%bb_BC_prescr     , DIVA%wbb_BC_prescr     )

    ! Calculate combined mesh operators for efficient calculation of the stiffness matrix
    CALL MatDestroy( DIVA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b, perr)
    CALL MatDestroy( DIVA%M_4_dudx_p_2_dvdy_buv_b              , perr)
    CALL MatDestroy( DIVA%M_dudy_p_dvdx_buv_b                  , perr)
    CALL MatDestroy( DIVA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b, perr)
    CALL MatDestroy( DIVA%M_4_dvdy_p_2_dudx_buv_b              , perr)
    CALL MatDestroy( DIVA%M_dvdx_p_dudy_buv_b                  , perr)
    CALL calc_combined_mesh_operators_DIVA( mesh_new, DIVA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_DIVA_solver

! == Calculate several intermediate terms in the DIVA

  SUBROUTINE calc_driving_stress( mesh, ice, DIVA)
    ! Calculate the driving stress

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_driving_stress'
    REAL(dp), DIMENSION(:    ), POINTER                          :: Hi_b
    REAL(dp), DIMENSION(:    ), POINTER                          :: dHs_dx_b
    REAL(dp), DIMENSION(:    ), POINTER                          :: dHs_dy_b
    INTEGER                                                      :: wHi_b, wdHs_dx_b, wdHs_dy_b
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, Hi_b    , wHi_b    )
    CALL allocate_shared_dp_1D( mesh%nTri, dHs_dx_b, wdHs_dx_b)
    CALL allocate_shared_dp_1D( mesh%nTri, dHs_dy_b, wdHs_dy_b)

    ! Calculate Hi, dHs/dx, and dHs/dy on the b-grid
    CALL map_a_to_b_2D( mesh, ice%Hi_a, Hi_b    )
    CALL ddx_a_to_b_2D( mesh, ice%Hs_a, dHs_dx_b)
    CALL ddy_a_to_b_2D( mesh, ice%Hs_a, dHs_dy_b)

    ! Calculate the driving stress
    DO ti = mesh%ti1, mesh%ti2
      DIVA%taudx_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dx_b( ti)
      DIVA%taudy_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dy_b( ti)
    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wHi_b    )
    CALL deallocate_shared( wdHs_dx_b)
    CALL deallocate_shared( wdHs_dy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_driving_stress

  SUBROUTINE calc_strain_rates( mesh, ice, DIVA)
    ! Calculate the strain rates

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_strain_rates'
    INTEGER                                                      :: ti,k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Horizontal strain rates are calculated explicitly
    CALL ddx_b_to_a_2D( mesh, DIVA%u_vav_b, DIVA%du_dx_a)
    CALL ddy_b_to_a_2D( mesh, DIVA%u_vav_b, DIVA%du_dy_a)
    CALL ddx_b_to_a_2D( mesh, DIVA%v_vav_b, DIVA%dv_dx_a)
    CALL ddy_b_to_a_2D( mesh, DIVA%v_vav_b, DIVA%dv_dy_a)

    ! Vertical shear strain rates are parameterised (Lipscomb et al., 2019, Eq. 36)
    ! (note that the (s-z)/H in the integral is equal to our zeta coordinate)
    DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        DIVA%du_dz_3D_b( ti,k) = DIVA%taubx_b( ti) * C%zeta( k) / MAX( DIVA%eta_3D_b( ti,k), C%DIVA_visc_eff_min)
        DIVA%dv_dz_3D_b( ti,k) = DIVA%tauby_b( ti) * C%zeta( k) / MAX( DIVA%eta_3D_b( ti,k), C%DIVA_visc_eff_min)
      END DO
    END DO
    CALL sync

    ! Map vertical shear strain rates from the b-grid to the a-grid
    CALL map_b_to_a_3D( mesh, DIVA%du_dz_3D_b, DIVA%du_dz_3D_a)
    CALL map_b_to_a_3D( mesh, DIVA%dv_dz_3D_b, DIVA%dv_dz_3D_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_strain_rates

  SUBROUTINE calc_effective_viscosity( mesh, ice, DIVA)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_effective_viscosity'
    INTEGER                                                      :: vi,k
    REAL(dp)                                                     :: epsilon_sq

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Calculate the effective viscosity eta in 3-D on the a-grid
      DO k = 1, C%nz

        ! Calculate the square of the effective strain rate epsilon
        epsilon_sq = DIVA%du_dx_a( vi)**2 + &
                     DIVA%dv_dy_a( vi)**2 + &
                     DIVA%du_dx_a( vi) * DIVA%dv_dy_a( vi) + &
                     0.25_dp * (DIVA%du_dy_a( vi) + DIVA%dv_dx_a( vi))**2 + &
                     0.25_dp * (DIVA%du_dz_3D_a( vi,k)**2 + DIVA%dv_dz_3D_a( vi,k)**2) + &
                     C%DIVA_epsilon_sq_0

        ! Calculate the effective viscosity eta
        DIVA%eta_3D_a( vi,k) = 0.5_dp * ice%A_flow_3D_a( vi,k)**(-1._dp/  C%n_flow) * (epsilon_sq)**((1._dp - C%n_flow)/(2._dp*C%n_flow))

        ! Safety
        DIVA%eta_3D_a( vi,k) = MAX( DIVA%eta_3D_a( vi,k), C%DIVA_visc_eff_min)

      END DO

      ! Calculate the vertical integral of the effective viscosity
      CALL vertical_integrate( C%zeta, DIVA%eta_3D_a( vi,:), DIVA%eta_vint_a( vi))

      ! Calculate the product term N = eta * H
      DIVA%N_a( vi) = DIVA%eta_vint_a( vi) * MAX( 0.1_dp, ice%Hi_a( vi))

    END DO
    CALL sync

    ! Map 3-D effective viscosity to the b-grid
    CALL map_a_to_b_3D( mesh, DIVA%eta_3D_a, DIVA%eta_3D_b)

    ! Calculate the product term N and its gradients on the b-grid
    CALL map_a_to_b_2D( mesh, DIVA%N_a, DIVA%N_b    )
    CALL ddx_a_to_b_2D( mesh, DIVA%N_a, DIVA%dN_dx_b)
    CALL ddy_a_to_b_2D( mesh, DIVA%N_a, DIVA%dN_dy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_effective_viscosity

  FUNCTION Fn( Hi, Hs, eta_3D, n) RESULT( Fn_val)
    ! The Fn-function (Lipscomb et al. (2019), Eq. 30)

    IMPLICIT NONE

    ! In- and output variables:
    REAL(dp),                            INTENT(IN)              :: Hi
    REAL(dp),                            INTENT(IN)              :: Hs
    REAL(dp), DIMENSION(:    ),          INTENT(IN)              :: eta_3D
    REAL(dp),                            INTENT(IN)              :: n
    REAL(dp)                                                     :: Fn_val

    ! Local variables:
    REAL(dp), DIMENSION( C%nz)                                   :: z

    z = Hs - (Hi * C%zeta)

    CALL vertical_integrate( z, (1._dp / eta_3D) * ((Hs - z) / Hi)**n, Fn_val)

  END FUNCTION Fn

  SUBROUTINE calc_beta_eff( mesh, ice, DIVA)
    ! Calculate the sliding + vertical shearing term beta_eff in the DIVA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_beta_eff'
    REAL(dp)                                                     :: F2_min
    INTEGER                                                      :: vi,ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Calculate F2

    ! Set a lower limit for F2 to improve numerical stability
    CALL vertical_integrate( C%zeta, (1._dp / C%DIVA_visc_eff_min) * C%zeta**2._dp, F2_min)

    DO vi = mesh%vi1, mesh%vi2
      ! Lipscomb et al. (2019), Eq. 30
      DIVA%F2_a( vi) = MAX( F2_min, Fn( ice%Hi_a( vi), ice%Hs_a( vi), DIVA%eta_3D_a( vi,:), 2._dp))
    END DO
    CALL sync

  ! == Calculate beta_eff

    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! No basal sliding allowed; calculate beta_eff according to Lipscomb et al. (2019), Eq. 35

      DO vi = mesh%vi1, mesh%vi2
        DIVA%beta_eff_a( vi) = 1._dp / DIVA%F2_a( vi)
      END DO
      CALL sync

    ELSE
      ! Basal sliding is allowed; calculate beta_eff according to Lipscomb et al. (2019), Eq. 33

      ! Calculate the basal friction coefficient beta_b for the current velocity solution
      CALL calc_basal_friction_coefficient( mesh, ice, DIVA%u_base_b, DIVA%v_base_b)

      DO vi = mesh%vi1, mesh%vi2
        DIVA%beta_eff_a( vi) = ice%beta_b_a( vi) / (1._dp + ice%beta_b_a( vi) * DIVA%F2_a( vi))
      END DO
      CALL sync

    END IF

    ! Map beta_eff from the a-grid to the b-grid
    CALL map_a_to_b_2D( mesh, DIVA%beta_eff_a, DIVA%beta_eff_b)

    ! Apply the sub-grid grounded fraction
    IF (C%do_GL_subgrid_friction) THEN
      DO ti = mesh%ti1, mesh%ti2
        DIVA%beta_eff_b( ti) = DIVA%beta_eff_b( ti) * ice%f_grnd_b( ti)**2
      END DO
      CALL sync
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_beta_eff

  SUBROUTINE calc_basal_shear_stress( mesh, DIVA)
    ! Calculate the basal shear stress in the DIVA
    !
    ! This is approximated by tau_b = beta_b * u_base ~= beta_eff * u_vav

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_basal_stress_DIVA'
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
      DIVA%taubx_b( ti) = DIVA%beta_eff_b( ti) * DIVA%u_vav_b( ti)
      DIVA%tauby_b( ti) = DIVA%beta_eff_b( ti) * DIVA%v_vav_b( ti)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_shear_stress

  SUBROUTINE calc_basal_velocities( mesh, ice, DIVA)
    ! Calculate basal velocities in the DIVA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_basal_velocitiesA'
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Exception for the case of no sliding
    IF (C%choice_sliding_law == 'no_sliding') THEN
      DIVA%u_base_b( mesh%ti1:mesh%ti2) = 0._dp
      DIVA%v_base_b( mesh%ti1:mesh%ti2) = 0._dp
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Map beta_b and F2 to the b-grid
    CALL map_a_to_b_2D( mesh, DIVA%F2_a   , DIVA%F2_b    )
    CALL map_a_to_b_2D( mesh, ice%beta_b_a, DIVA%beta_b_b)

    ! Calculate basal velocities on the b-grid according to Lipscomb et al. (2019), Eq. 32
    DO ti = mesh%ti1, mesh%ti2
      DIVA%u_base_b( ti) = DIVA%u_vav_b( ti) / (1._dp + DIVA%beta_b_b( ti) * DIVA%F2_b( ti))
      DIVA%v_base_b( ti) = DIVA%v_vav_b( ti) / (1._dp + DIVA%beta_b_b( ti) * DIVA%F2_b( ti))
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_velocities

  SUBROUTINE calc_3D_horizontal_velocities( mesh, ice, DIVA)
    ! Calculate the 3D horizontal velocity fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_3D_horizontal_velocities'
    INTEGER                                                      :: ti
    REAL(dp), DIMENSION(:    ), POINTER                          ::  Hi_b,  Hs_b
    INTEGER                                                      :: wHi_b, wHs_b
    REAL(dp), DIMENSION( C%nz)                                   :: z
    REAL(dp), DIMENSION( C%nz)                                   :: int_sminz_over_eta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, Hi_b, wHi_b)
    CALL allocate_shared_dp_1D( mesh%nTri, Hs_b, wHs_b)

    ! Map ice thickness and surface elevation to the b-grid
    CALL map_a_to_b_2D( mesh, ice%Hi_a, Hi_b)
    CALL map_a_to_b_2D( mesh, ice%Hs_a, Hs_b)

    ! Calculate 3D horizontal velocity components according to Lipscomb et al. (2019), Eq. 29
    DO ti = mesh%ti1, mesh%ti2
      z = Hs_b( ti) - Hi_b( ti) * C%zeta
      CALL vertical_integration_from_bottom_to_zeta( z, (Hs_b( ti) - z) / DIVA%eta_3D_b( ti,:), int_sminz_over_eta)
      DIVA%u_3D_b( ti,:) = DIVA%u_base_b( ti) + DIVA%taubx_b( ti) * int_sminz_over_eta
      DIVA%v_3D_b( ti,:) = DIVA%v_base_b( ti) + DIVA%tauby_b( ti) * int_sminz_over_eta
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wHi_b)
    CALL deallocate_shared( wHs_b)

    ! DENK DROM
!     IF (par%master) THEN
!       CALL warning('u_vav_b = [{dp_01} - {dp_02}], v_vav_b = [{dp_03} - {dp_04}]', &
!         dp_01 = MINVAL( DIVA%u_vav_b), dp_02 = MAXVAL( DIVA%u_vav_b), dp_03 = MINVAL( DIVA%v_vav_b), dp_04 = MAXVAL( DIVA%v_vav_b))
!       CALL warning('u_3D_b = [{dp_01} - {dp_02}], v_3D_b = [{dp_03} - {dp_04}]', &
!         dp_01 = MINVAL( DIVA%u_3D_b), dp_02 = MAXVAL( DIVA%u_3D_b), dp_03 = MINVAL( DIVA%v_3D_b), dp_04 = MAXVAL( DIVA%v_3D_b))
!       CALL warning('u_base_b = [{dp_01} - {dp_02}], v_base_b = [{dp_03} - {dp_04}]', &
!         dp_01 = MINVAL( DIVA%u_base_b), dp_02 = MAXVAL( DIVA%u_base_b), dp_03 = MINVAL( DIVA%v_base_b), dp_04 = MAXVAL( DIVA%v_base_b))
!     END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_3D_horizontal_velocities

! == Solve the linearised DIVA

  SUBROUTINE solve_DIVA_linearised( mesh, DIVA, BC_prescr_u_b, BC_prescr_v_b)
    ! Solve the linearised DIVA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)              :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:    ),          INTENT(IN)              :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_DIVA_linearised'
    TYPE(tMat)                                                   ::  AA_masked
    REAL(dp), DIMENSION(:    ), POINTER                          ::  bb_masked
    INTEGER                                                      :: wbb_masked
    INTEGER                                                      :: n1,n2
    REAL(dp), DIMENSION(:    ), POINTER                          ::  uv_buv
    INTEGER                                                      :: wuv_buv
    INTEGER                                                      :: ti,nu,nv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nnbuv, bb_masked          , wbb_masked          )
    CALL allocate_shared_dp_1D( mesh%nnbuv, uv_buv             , wuv_buv             )

    ! Calculate the stiffness matrix and load vector representing the DIVA
    CALL calc_stiffness_matrix_DIVA_free( mesh, DIVA)

    ! Calculate the stiffness matrices and load vectors representing the different boundary conditions
    CALL calc_stiffness_matrix_DIVA_BC_west(   mesh, DIVA)
    CALL calc_stiffness_matrix_DIVA_BC_east(   mesh, DIVA)
    CALL calc_stiffness_matrix_DIVA_BC_south(  mesh, DIVA)
    CALL calc_stiffness_matrix_DIVA_BC_north(  mesh, DIVA)
    CALL calc_stiffness_matrix_DIVA_BC_prescr( mesh, DIVA, BC_prescr_u_b, BC_prescr_v_b)

    ! Use the boundary conditions mask matrices to combine the different stiffness matrices and load vectors

    ! Initialise
    CALL MatDuplicate( DIVA%AA_free, MAT_SHARE_NONZERO_PATTERN, DIVA%AA, perr)
    CALL partition_list( mesh%nnbuv, par%i, par%n, n1, n2)
    DIVA%bb( n1:n2) = 0._dp

    ! Free
    CALL MatMatMult(                           DIVA%m_free, DIVA%AA_free, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( DIVA%m_free, DIVA%bb_free,                                         bb_masked      )
    CALL MatAXPY( DIVA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    DIVA%bb( n1:n2) = DIVA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! West
    CALL MatMatMult(                           DIVA%m_BC_west, DIVA%AA_BC_west, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( DIVA%m_BC_west, DIVA%bb_BC_west,                                         bb_masked      )
    CALL MatAXPY( DIVA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    DIVA%bb( n1:n2) = DIVA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! East
    CALL MatMatMult(                           DIVA%m_BC_east, DIVA%AA_BC_east, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( DIVA%m_BC_east, DIVA%bb_BC_east,                                         bb_masked      )
    CALL MatAXPY( DIVA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    DIVA%bb( n1:n2) = DIVA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! South
    CALL MatMatMult(                           DIVA%m_BC_south, DIVA%AA_BC_south, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( DIVA%m_BC_south, DIVA%bb_BC_south,                                         bb_masked      )
    CALL MatAXPY( DIVA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    DIVA%bb( n1:n2) = DIVA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! North
    CALL MatMatMult(                           DIVA%m_BC_north, DIVA%AA_BC_north, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( DIVA%m_BC_north, DIVA%bb_BC_north,                                         bb_masked      )
    CALL MatAXPY( DIVA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    DIVA%bb( n1:n2) = DIVA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! Prescribed
    CALL MatMatMult(                           DIVA%m_BC_prescr, DIVA%AA_BC_prescr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( DIVA%m_BC_prescr, DIVA%bb_BC_prescr,                                         bb_masked      )
    CALL MatAXPY( DIVA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    DIVA%bb( n1:n2) = DIVA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! Combine the u and v velocities into a single vector;
    ! use the current velocity solution as the initial guess for the new one
    DO ti = mesh%ti1, mesh%ti2
      nu = mesh%tiuv2n( ti,1)
      nv = mesh%tiuv2n( ti,2)
      uv_buv( nu) = DIVA%u_vav_b( ti)
      uv_buv( nv) = DIVA%v_vav_b( ti)
    END DO

    ! Save the previous solution
    DIVA%u_vav_b_prev( mesh%ti1:mesh%ti2) = DIVA%u_vav_b( mesh%ti1:mesh%ti2)
    DIVA%v_vav_b_prev( mesh%ti1:mesh%ti2) = DIVA%v_vav_b( mesh%ti1:mesh%ti2)

    ! Solve the matrix equation
    CALL solve_matrix_equation_PETSc( DIVA%AA, DIVA%bb, uv_buv, DIVA%PETSc_rtol, DIVA%PETSc_abstol)

    ! Get the u and v velocities back in field form
    DO ti = mesh%ti1, mesh%ti2
      nu = mesh%tiuv2n( ti,1)
      nv = mesh%tiuv2n( ti,2)
      DIVA%u_vav_b( ti) = uv_buv( nu)
      DIVA%v_vav_b( ti) = uv_buv( nv)
    END DO

    ! Clean up after yourself
    CALL MatDestroy( DIVA%AA          , perr)
    CALL MatDestroy( DIVA%AA_free     , perr)
    CALL MatDestroy( DIVA%AA_BC_west  , perr)
    CALL MatDestroy( DIVA%AA_BC_east  , perr)
    CALL MatDestroy( DIVA%AA_BC_south , perr)
    CALL MatDestroy( DIVA%AA_BC_north , perr)
    CALL MatDestroy( DIVA%AA_BC_prescr, perr)
    CALL deallocate_shared( wbb_masked)
    CALL deallocate_shared( wuv_buv   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_DIVA_linearised

  SUBROUTINE calc_stiffness_matrix_DIVA_free( mesh, DIVA)
    ! Calculate the stiffness matrix and load vector representing the DIVA
    !
    ! The first equation of the DIVA reads:
    !
    !   d/dx [ 2 N ( 2 du/dx + dv/dy )] + d/dy [ N ( du/dy + dv/dx ) ] - beta_eff u = -taud,x
    !
    ! By applying the product rule, this expands to:
    !
    !   4 N d2u/dx2 + 4 dN/dx du/dx + 2 N d2v/dxdy + 2 dN/dx dv/dy + ...
    !     N d2u/dy2 +   dN/dy du/dy +   N d2v/dxdy +   dN/dy dv/dx - beta_eff u = -taud,x
    !
    ! Combining the terms involving N, dN/dx, and dN/dy yields:
    !
    !    N    ( 4 d2u/dx2 + 3 d2v/dxdy + d2u/dy2 ) + ...
    !   dN/dx ( 4 du/dx   + 2 dv/dy ) + ...
    !   dN/dy (   du/dy   +   dv/dx ) - beta_b u = -taud,x
    !
    ! The gradient operators between brackets are calculated beforehand.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_DIVA_free'
    TYPE(tVec)                                                   :: N_b_vec
    TYPE(tVec)                                                   :: dN_dx_b_vec
    TYPE(tVec)                                                   :: dN_dy_b_vec
    TYPE(tVec)                                                   :: beta_eff_b_vec
    TYPE(tMat)                                                   :: Au1, Au2, Au3, Au4, Au, Au_buv
    TYPE(tMat)                                                   :: Av1, Av2, Av3, Av4, Av, Av_buv
    INTEGER                                                      :: ti,nu,nv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get N, dN/dx, and dN/dy as PETSc vectors
    CALL vec_double2petsc( DIVA%N_b       , N_b_vec       )
    CALL vec_double2petsc( DIVA%dN_dx_b   , dN_dx_b_vec   )
    CALL vec_double2petsc( DIVA%dN_dy_b   , dN_dy_b_vec   )
    CALL vec_double2petsc( DIVA%beta_eff_b, beta_eff_b_vec)

  ! == Calculate stiffness matrix for Eq. 1

    ! Au1 = N    ( 4 d2u/dx2 + 3 d2v/dxdy + d2u/dy2 )
    CALL MatDuplicate( DIVA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b, MAT_COPY_VALUES, Au1, perr)
    CALL MatDiagonalScale( Au1, N_b_vec, PETSC_NULL_VEC, perr)

    ! Au2 = dN/dx ( 4 du/dx   + 2 dv/dy )
    CALL MatDuplicate( DIVA%M_4_dudx_p_2_dvdy_buv_b, MAT_COPY_VALUES, Au2, perr)
    CALL MatDiagonalScale( Au2, dN_dx_b_vec, PETSC_NULL_VEC, perr)

    ! Au3 = dN/dy (   du/dy   +   dv/dx )
    CALL MatDuplicate( DIVA%M_dudy_p_dvdx_buv_b, MAT_COPY_VALUES, Au3, perr)
    CALL MatDiagonalScale( Au3, dN_dy_b_vec, PETSC_NULL_VEC, perr)

    ! Au4 = -betab * u
    CALL MatDuplicate( mesh%M_map_bu_b, MAT_COPY_VALUES, Au4, perr)
    CALL MatDiagonalScale( Au4, beta_eff_b_vec, PETSC_NULL_VEC, perr)

    ! Combine all terms of Eq. 1
    CALL MatDuplicate( Au1, MAT_SHARE_NONZERO_PATTERN, Au, perr)
    CALL MatAXPY( Au,  1._dp, Au1, SAME_NONZERO_PATTERN  , perr)
    CALL MatAXPY( Au,  1._dp, Au2, SAME_NONZERO_PATTERN  , perr)
    CALL MatAXPY( Au,  1._dp, Au3, SAME_NONZERO_PATTERN  , perr)
    CALL MatAXPY( Au, -1._dp, Au4, SUBSET_NONZERO_PATTERN, perr)
    CALL MatMatMult( mesh%M_map_b_bu, Au, MAT_INITIAL_MATRIX, PETSC_DEFAuLT_REAL, Au_buv, perr)

    ! Clean up after yourself

  ! == Calculate stiffness matrix for Eq. 2

    ! Av1 = N    ( 4 d2v/dy2 + 3 d2u/dxdy + d2v/dx2 )
    CALL MatDuplicate( DIVA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b, MAT_COPY_VALUES, Av1, perr)
    CALL MatDiagonalScale( Av1, N_b_vec, PETSC_NULL_VEC, perr)

    ! Av2 = dN/dy ( 4 dv/dy   + 2 du/dx )
    CALL MatDuplicate( DIVA%M_4_dvdy_p_2_dudx_buv_b, MAT_COPY_VALUES, Av2, perr)
    CALL MatDiagonalScale( Av2, dN_dy_b_vec, PETSC_NULL_VEC, perr)

    ! Av3 = dN/dx (   dv/dx   +   du/dy )
    CALL MatDuplicate( DIVA%M_dvdx_p_dudy_buv_b, MAT_COPY_VALUES, Av3, perr)
    CALL MatDiagonalScale( Av3, dN_dx_b_vec, PETSC_NULL_VEC, perr)

    ! Av4 = -betab * v
    CALL MatDuplicate( mesh%M_map_bv_b, MAT_COPY_VALUES, Av4, perr)
    CALL MatDiagonalScale( Av4, beta_eff_b_vec, PETSC_NULL_VEC, perr)

    ! Combine all terms of Eq. 2
    CALL MatDuplicate( Av1, MAT_SHARE_NONZERO_PATTERN, Av, perr)
    CALL MatAXPY( Av,  1._dp, Av1, SAME_NONZERO_PATTERN  , perr)
    CALL MatAXPY( Av,  1._dp, Av2, SAME_NONZERO_PATTERN  , perr)
    CALL MatAXPY( Av,  1._dp, Av3, SAME_NONZERO_PATTERN  , perr)
    CALL MatAXPY( Av, -1._dp, Av4, SUBSET_NONZERO_PATTERN, perr)
    CALL MatMatMult( mesh%M_map_b_bv, Av, MAT_INITIAL_MATRIX, PETSC_DEFAuLT_REAL, Av_buv, perr)

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au_buv, MAT_SHARE_NONZERO_PATTERN, DIVA%AA_free, perr)
    CALL MatAXPY( DIVA%AA_free, 1._dp, Au_buv, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( DIVA%AA_free, 1._dp, Av_buv, UNKNOWN_NONZERO_PATTERN, perr)

  ! == Calculate the load vector

    DO ti = mesh%ti1, mesh%ti2

      nu = mesh%tiuv2n( ti,1)
      nv = mesh%tiuv2n( ti,2)

      DIVA%bb_free( nu) = -DIVA%taudx_b( ti)
      DIVA%bb_free( nv) = -DIVA%taudy_b( ti)

    END DO
    CALL sync

    ! Clean up after yourself
    CALL VecDestroy( N_b_vec       , perr)
    CALL VecDestroy( dN_dx_b_vec   , perr)
    CALL VecDestroy( dN_dy_b_vec   , perr)
    CALL VecDestroy( beta_eff_b_vec, perr)
    CALL MatDestroy( Au1           , perr)
    CALL MatDestroy( Au2           , perr)
    CALL MatDestroy( Au3           , perr)
    CALL MatDestroy( Au4           , perr)
    CALL MatDestroy( Au            , perr)
    CALL MatDestroy( Au_buv        , perr)
    CALL MatDestroy( Av1           , perr)
    CALL MatDestroy( Av2           , perr)
    CALL MatDestroy( Av3           , perr)
    CALL MatDestroy( Av4           , perr)
    CALL MatDestroy( Av            , perr)
    CALL MatDestroy( Av_buv        , perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_DIVA_free

  SUBROUTINE calc_stiffness_matrix_DIVA_BC_west( mesh, DIVA)
    ! Calculate the stiffness matrix and load vector representing boundary conditions to the DIVA on the western domain boundary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_DIVA_BC_west'
    TYPE(tMat)                                                   :: M_dudx, Au
    TYPE(tMat)                                                   :: M_dvdx, Av
    INTEGER                                                      :: ti,nu,nv,ti_copy

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == u

    IF     (C%DIVA_boundary_BC_u_west == 'infinite') THEN
      ! du/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddx_b_b, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx, perr)
      CALL MatMatMult( mesh%M_map_b_bu, M_dudx         , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au    , perr)
      CALL MatDestroy( M_dudx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_west( nu) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_west == 'zero') THEN
      ! u = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bu, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_west( nu) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_west == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_u_west == 'periodic_ISMIP_HOM') THEN
      ! Periodic boundary conditions in the ISMIP-HOM experiments are implemented by
      ! taking advantage of the fact that u(x,y) = u(x+L/2,y+L/2)
      !
      ! Velocities at the boundary can therefore be set equal to the interior value
      ! diagonally across from the boundary point (displaced by [L/2,L/2])

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bu, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_west( nu) = DIVA%u_vav_b( ti_copy)
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_u_west "' // TRIM( C%DIVA_boundary_BC_u_west) // '"!')
    END IF

  ! == v

    IF     (C%DIVA_boundary_BC_v_west == 'infinite') THEN
      ! dv/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddx_b_b, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx, perr)
      CALL MatMatMult( mesh%M_map_b_bv, M_dvdx         , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av    , perr)
      CALL MatDestroy( M_dvdx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_west( nv) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_west == 'zero') THEN
      ! v = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bv, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_west( nv) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_west == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_v_west == 'periodic_ISMIP_HOM') THEN
      ! Do nothing here; instead, velocities are prescribed through the BC_prescr arguments

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bv, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_west( nv) = DIVA%v_vav_b( ti_copy)
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_v_west "' // TRIM( C%DIVA_boundary_BC_v_west) // '"!')
    END IF

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au, MAT_SHARE_NONZERO_PATTERN, DIVA%AA_BC_west, perr)
    CALL MatAXPY( DIVA%AA_BC_west, 1._dp, Au, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( DIVA%AA_BC_west, 1._dp, Av, UNKNOWN_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( Au, perr)
    CALL MatDestroy( Av, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_DIVA_BC_west

  SUBROUTINE calc_stiffness_matrix_DIVA_BC_east( mesh, DIVA)
    ! Calculate the stiffness matrix and load vector representing boundary conditions to the DIVA on the eastern domain boundary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_DIVA_BC_east'
    TYPE(tMat)                                                   :: M_dudx, Au
    TYPE(tMat)                                                   :: M_dvdx, Av
    INTEGER                                                      :: ti,nu,nv,ti_copy

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == u

    IF     (C%DIVA_boundary_BC_u_east == 'infinite') THEN
      ! du/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddx_b_b, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx, perr)
      CALL MatMatMult( mesh%M_map_b_bu, M_dudx         , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au    , perr)
      CALL MatDestroy( M_dudx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_east( nu) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_east == 'zero') THEN
      ! u = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bu, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_east( nu) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_east == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_u_east == 'periodic_ISMIP_HOM') THEN
      ! Periodic boundary conditions in the ISMIP-HOM experiments are implemented by
      ! taking advantage of the fact that u(x,y) = u(x+L/2,y+L/2)
      !
      ! Velocities at the boundary can therefore be set equal to the interior value
      ! diagonally across from the boundary point (displaced by [L/2,L/2])

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bu, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_east( nu) = DIVA%u_vav_b( ti_copy)
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_u_east "' // TRIM( C%DIVA_boundary_BC_u_east) // '"!')
    END IF

  ! == v

    IF     (C%DIVA_boundary_BC_v_east == 'infinite') THEN
      ! dv/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddx_b_b, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx, perr)
      CALL MatMatMult( mesh%M_map_b_bv, M_dvdx         , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av    , perr)
      CALL MatDestroy( M_dvdx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_east( nv) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_east == 'zero') THEN
      ! v = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bv, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_east( nv) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_east == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_v_east == 'periodic_ISMIP_HOM') THEN
      ! Do nothing here; instead, velocities are prescribed through the BC_prescr arguments

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bv, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_east( nv) = DIVA%v_vav_b( ti_copy)
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_v_east "' // TRIM( C%DIVA_boundary_BC_v_east) // '"!')
    END IF

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au, MAT_SHARE_NONZERO_PATTERN, DIVA%AA_BC_east, perr)
    CALL MatAXPY( DIVA%AA_BC_east, 1._dp, Au, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( DIVA%AA_BC_east, 1._dp, Av, UNKNOWN_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( Au, perr)
    CALL MatDestroy( Av, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_DIVA_BC_east

  SUBROUTINE calc_stiffness_matrix_DIVA_BC_south( mesh, DIVA)
    ! Calculate the stiffness matrix and load vector representing boundary conditions to the DIVA on the southern domain boundary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_DIVA_BC_south'
    TYPE(tMat)                                                   :: M_dudx, Au
    TYPE(tMat)                                                   :: M_dvdx, Av
    INTEGER                                                      :: ti,nu,nv,ti_copy

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == u

    IF     (C%DIVA_boundary_BC_u_south == 'infinite') THEN
      ! du/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddy_b_b, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx, perr)
      CALL MatMatMult( mesh%M_map_b_bu, M_dudx         , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au    , perr)
      CALL MatDestroy( M_dudx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_south( nu) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_south == 'zero') THEN
      ! u = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bu, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_south( nu) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_south == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_u_south == 'periodic_ISMIP_HOM') THEN
      ! Periodic boundary conditions in the ISMIP-HOM experiments are implemented by
      ! taking advantage of the fact that u(x,y) = u(x+L/2,y+L/2)
      !
      ! Velocities at the boundary can therefore be set equal to the interior value
      ! diagonally across from the boundary point (displaced by [L/2,L/2])

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bu, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_south( nu) = DIVA%u_vav_b( ti_copy)
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_u_south "' // TRIM( C%DIVA_boundary_BC_u_south) // '"!')
    END IF

  ! == v

    IF     (C%DIVA_boundary_BC_v_south == 'infinite') THEN
      ! dv/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddy_b_b, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx, perr)
      CALL MatMatMult( mesh%M_map_b_bv, M_dvdx         , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av    , perr)
      CALL MatDestroy( M_dvdx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_south( nv) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_south == 'zero') THEN
      ! v = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bv, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_south( nv) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_south == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_v_south == 'periodic_ISMIP_HOM') THEN
      ! Do nothing here; instead, velocities are prescribed through the BC_prescr arguments

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bv, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_south( nv) = DIVA%v_vav_b( ti_copy)
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_v_south "' // TRIM( C%DIVA_boundary_BC_v_south) // '"!')
    END IF

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au, MAT_SHARE_NONZERO_PATTERN, DIVA%AA_BC_south, perr)
    CALL MatAXPY( DIVA%AA_BC_south, 1._dp, Au, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( DIVA%AA_BC_south, 1._dp, Av, UNKNOWN_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( Au, perr)
    CALL MatDestroy( Av, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_DIVA_BC_south

  SUBROUTINE calc_stiffness_matrix_DIVA_BC_north( mesh, DIVA)
    ! Calculate the stiffness matrix and load vector representing boundary conditions to the DIVA on the northern domain boundary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_DIVA_BC_north'
    TYPE(tMat)                                                   :: M_dudx, Au
    TYPE(tMat)                                                   :: M_dvdx, Av
    INTEGER                                                      :: ti,nu,nv,ti_copy

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == u

    IF     (C%DIVA_boundary_BC_u_north == 'infinite') THEN
      ! du/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddy_b_b, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx, perr)
      CALL MatMatMult( mesh%M_map_b_bu, M_dudx         , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au    , perr)
      CALL MatDestroy( M_dudx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_north( nu) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_north == 'zero') THEN
      ! u = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bu, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_north( nu) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_north == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_u_north == 'periodic_ISMIP_HOM') THEN
      ! Periodic boundary conditions in the ISMIP-HOM experiments are implemented by
      ! taking advantage of the fact that u(x,y) = u(x+L/2,y+L/2)
      !
      ! Velocities at the boundary can therefore be set equal to the interior value
      ! diagonally across from the boundary point (displaced by [L/2,L/2])

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bu, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        nu = mesh%tiuv2n( ti,1)
        DIVA%bb_BC_north( nu) = DIVA%u_vav_b( ti_copy)
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_u_north "' // TRIM( C%DIVA_boundary_BC_u_north) // '"!')
    END IF

  ! == v

    IF     (C%DIVA_boundary_BC_v_north == 'infinite') THEN
      ! dv/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddy_b_b, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx, perr)
      CALL MatMatMult( mesh%M_map_b_bv, M_dvdx         , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av    , perr)
      CALL MatDestroy( M_dvdx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_north( nv) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_north == 'zero') THEN
      ! v = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bv, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_north( nv) = 0._dp
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_north == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_v_north == 'periodic_ISMIP_HOM') THEN
      ! Do nothing here; instead, velocities are prescribed through the BC_prescr arguments

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_b_bv, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        nv = mesh%tiuv2n( ti,2)
        DIVA%bb_BC_north( nv) = DIVA%v_vav_b( ti_copy)
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_v_north "' // TRIM( C%DIVA_boundary_BC_v_north) // '"!')
    END IF

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au, MAT_SHARE_NONZERO_PATTERN, DIVA%AA_BC_north, perr)
    CALL MatAXPY( DIVA%AA_BC_north, 1._dp, Au, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( DIVA%AA_BC_north, 1._dp, Av, UNKNOWN_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( Au, perr)
    CALL MatDestroy( Av, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_DIVA_BC_north

  SUBROUTINE calc_stiffness_matrix_DIVA_BC_prescr( mesh, DIVA, BC_prescr_u_b, BC_prescr_v_b)
    ! Calculate the stiffness matrix and load vector representing boundary conditions in the form of directly prescribed velocities

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)              :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:    ),          INTENT(IN)              :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_DIVA_BC_prescr'
    TYPE(tMat)                                                   :: Au, Av
    INTEGER                                                      :: ti,nu,nv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Stiffness matrices
    CALL MatMatMult( mesh%M_map_b_bu, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)
    CALL MatMatMult( mesh%M_map_b_bv, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

    ! Combine Eqs. 1 and 2
    CALL MatDuplicate( Au, MAT_SHARE_NONZERO_PATTERN, DIVA%AA_BC_prescr, perr)
    CALL MatAXPY( DIVA%AA_BC_prescr, 1._dp, Au, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( DIVA%AA_BC_prescr, 1._dp, Av, UNKNOWN_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( Au, perr)
    CALL MatDestroy( Av, perr)

    ! Load vector
    DO ti = mesh%ti1, mesh%ti2
      nu = mesh%tiuv2n( ti,1)
      nv = mesh%tiuv2n( ti,2)
      DIVA%bb_BC_prescr( nu) = BC_prescr_u_b( ti)
      DIVA%bb_BC_prescr( nv) = BC_prescr_v_b( ti)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_DIVA_BC_prescr

  SUBROUTINE calc_DIVA_BC_mask_matrices( mesh, DIVA, BC_prescr_mask_b)
    ! Calculate boundary conditions mask matrices

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA
    INTEGER,  DIMENSION(:    ),          INTENT(IN)              :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_DIVA_BC_mask_matrices'
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_free
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_west
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_east
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_north
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_south
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_prescr
    INTEGER                                                      :: wm_free, wm_BC_west, wm_BC_east, wm_BC_north, wm_BC_south, wm_BC_prescr
    INTEGER                                                      :: ti,nu,nv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nnbuv, m_free     , wm_free     )
    CALL allocate_shared_dp_1D( mesh%nnbuv, m_BC_west  , wm_BC_west  )
    CALL allocate_shared_dp_1D( mesh%nnbuv, m_BC_east  , wm_BC_east  )
    CALL allocate_shared_dp_1D( mesh%nnbuv, m_BC_north , wm_BC_north )
    CALL allocate_shared_dp_1D( mesh%nnbuv, m_BC_south , wm_BC_south )
    CALL allocate_shared_dp_1D( mesh%nnbuv, m_BC_prescr, wm_BC_prescr)

    DO ti = mesh%ti1, mesh%ti2

      nu = mesh%tiuv2n( ti,1)
      nv = mesh%tiuv2n( ti,2)

      ! Initialise everything with zero
      m_free(      nu) = 0._dp
      m_BC_west(   nu) = 0._dp
      m_BC_east(   nu) = 0._dp
      m_BC_south(  nu) = 0._dp
      m_BC_north(  nu) = 0._dp
      m_BC_prescr( nu) = 0._dp

      m_free(      nv) = 0._dp
      m_BC_west(   nv) = 0._dp
      m_BC_east(   nv) = 0._dp
      m_BC_south(  nv) = 0._dp
      m_BC_north(  nv) = 0._dp
      m_BC_prescr( nv) = 0._dp

      IF     (BC_prescr_mask_b( ti) == 1) THEN
        ! Velocities are prescribed at this triangle; this overrides everything else

        m_BC_prescr( nu) = 1._dp
        m_BC_prescr( nv) = 1._dp

      ELSEIF (mesh%Tri_edge_index( ti) == 1 .OR. mesh%Tri_edge_index( ti) == 2) THEN
        ! This triangle lies on the northern domain boundary

        m_BC_north( nu) = 1._dp
        m_BC_north( nv) = 1._dp

      ELSEIF (mesh%Tri_edge_index( ti) == 3 .OR. mesh%Tri_edge_index( ti) == 4) THEN
        ! This triangle lies on the eastern domain boundary

        m_BC_east( nu) = 1._dp
        m_BC_east( nv) = 1._dp

      ELSEIF (mesh%Tri_edge_index( ti) == 5 .OR. mesh%Tri_edge_index( ti) == 6) THEN
        ! This triangle lies on the southern domain boundary

        m_BC_south( nu) = 1._dp
        m_BC_south( nv) = 1._dp

      ELSEIF (mesh%Tri_edge_index( ti) == 7 .OR. mesh%Tri_edge_index( ti) == 8) THEN
        ! This triangle lies on the western domain boundary

        m_BC_west( nu) = 1._dp
        m_BC_west( nv) = 1._dp

      ELSE
        ! No boundary conditions apply; solve the DIVA instead

        m_free( nu) = 1._dp
        m_free( nv) = 1._dp

      END IF

    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync

    ! Convert to diagonal matrices
    CALL double2petscmat( m_free     , DIVA%m_free     )
    CALL double2petscmat( m_BC_west  , DIVA%m_BC_west  )
    CALL double2petscmat( m_BC_east  , DIVA%m_BC_east  )
    CALL double2petscmat( m_BC_south , DIVA%m_BC_south )
    CALL double2petscmat( m_BC_north , DIVA%m_BC_north )
    CALL double2petscmat( m_BC_prescr, DIVA%m_BC_prescr)

    ! Clean up after yourself
    CALL deallocate_shared( wm_free     )
    CALL deallocate_shared( wm_BC_west  )
    CALL deallocate_shared( wm_BC_east  )
    CALL deallocate_shared( wm_BC_north )
    CALL deallocate_shared( wm_BC_south )
    CALL deallocate_shared( wm_BC_prescr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_DIVA_BC_mask_matrices

! == Some useful tools for improving numerical stability of the viscosity iteration

  SUBROUTINE relax_viscosity_iterations( mesh, DIVA, R)
    ! Reduce the change between velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA
    REAL(dp),                            INTENT(IN)              :: R

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'relax_viscosity_iterations'
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
      DIVA%u_vav_b( ti) = (R * DIVA%u_vav_b( ti)) + ((1._dp - R) * DIVA%u_vav_b_prev( ti))
      DIVA%v_vav_b( ti) = (R * DIVA%v_vav_b( ti)) + ((1._dp - R) * DIVA%v_vav_b_prev( ti))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE relax_viscosity_iterations

  SUBROUTINE calc_visc_iter_UV_resid( mesh, DIVA, resid_UV)
    ! Calculate the L2-norm of the two consecutive velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_velocity_solver_DIVA),      INTENT(IN)              :: DIVA
    REAL(dp),                            INTENT(OUT)             :: resid_UV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_visc_iter_UV_resid'
    INTEGER                                                      :: ierr
    INTEGER                                                      :: ti
    REAL(dp)                                                     :: res1, res2

    ! Add routine to path
    CALL init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    DO ti = mesh%ti1, mesh%ti2

      res1 = res1 + (DIVA%u_vav_b( ti) - DIVA%u_vav_b_prev( ti))**2
      res1 = res1 + (DIVA%v_vav_b( ti) - DIVA%v_vav_b_prev( ti))**2

      res2 = res2 + (DIVA%u_vav_b( ti) + DIVA%u_vav_b_prev( ti))**2
      res2 = res2 + (DIVA%v_vav_b( ti) + DIVA%v_vav_b_prev( ti))**2

    END DO

    ! Combine results from all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate residual
    resid_UV = 2._dp * res1 / MAX( res2, 1E-8_dp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_visc_iter_UV_resid

  SUBROUTINE apply_velocity_limits( mesh, DIVA)
    ! Limit velocities for improved stability

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'apply_velocity_limits'
    INTEGER                                                      :: ti
    REAL(dp)                                                     :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2

      ! Calculate absolute speed
      uabs = SQRT( DIVA%u_vav_b( ti)**2 + DIVA%v_vav_b( ti)**2)

      ! Reduce velocities if necessary
      IF (uabs > C%DIVA_vel_max) THEN
        DIVA%u_vav_b( ti) = DIVA%u_vav_b( ti) * C%DIVA_vel_max / uabs
        DIVA%v_vav_b( ti) = DIVA%v_vav_b( ti) * C%DIVA_vel_max / uabs
      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_velocity_limits

! == Calculate combined mesh operators for efficient calculation of the stiffness matrix

  SUBROUTINE calc_combined_mesh_operators_DIVA( mesh, DIVA)
    ! Calculate combined mesh operators for efficient calculation of the stiffness matrix

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_DIVA),     INTENT(INOUT)           :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_combined_mesh_operators_DIVA'
    TYPE(tMat)                                                   :: M_dudx_buv_b
    TYPE(tMat)                                                   :: M_dudy_buv_b
    TYPE(tMat)                                                   :: M_d2udx2_buv_b
    TYPE(tMat)                                                   :: M_d2udxdy_buv_b
    TYPE(tMat)                                                   :: M_d2udy2_buv_b
    TYPE(tMat)                                                   :: M_dvdx_buv_b
    TYPE(tMat)                                                   :: M_dvdy_buv_b
    TYPE(tMat)                                                   :: M_d2vdx2_buv_b
    TYPE(tMat)                                                   :: M_d2vdxdy_buv_b
    TYPE(tMat)                                                   :: M_d2vdy2_buv_b
    TYPE(tMat)                                                   :: M1
    REAL(dp), DIMENSION(mesh%nnb)                                :: d_4
    TYPE(tVec)                                                   :: v_4

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Vector containing only  and 4's
    d_4 = 4._dp
    CALL vec_double2petsc( d_4, v_4)

    ! Partial derivatives of u on the b-grid
    CALL MatMatMult( mesh%M2_ddx_b_b   , mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx_buv_b   , perr)
    CALL MatMatMult( mesh%M2_ddy_b_b   , mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudy_buv_b   , perr)
    CALL MatMatMult( mesh%M2_d2dx2_b_b , mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2udx2_buv_b , perr)
    CALL MatMatMult( mesh%M2_d2dxdy_b_b, mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2udxdy_buv_b, perr)
    CALL MatMatMult( mesh%M2_d2dy2_b_b , mesh%M_map_bu_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2udy2_buv_b , perr)

    ! Partial derivatives of v on the b-grid
    CALL MatMatMult( mesh%M2_ddx_b_b   , mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx_buv_b   , perr)
    CALL MatMatMult( mesh%M2_ddy_b_b   , mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdy_buv_b   , perr)
    CALL MatMatMult( mesh%M2_d2dx2_b_b , mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2vdx2_buv_b , perr)
    CALL MatMatMult( mesh%M2_d2dxdy_b_b, mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2vdxdy_buv_b, perr)
    CALL MatMatMult( mesh%M2_d2dy2_b_b , mesh%M_map_bv_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2vdy2_buv_b , perr)

    ! M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b
    CALL MatDuplicate( M_d2udx2_buv_b, MAT_COPY_VALUES, DIVA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b, perr)
    CALL MatDiagonalScale( DIVA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b, v_4, PETSC_NULL_VEC, perr)
    CALL MatAXPY( DIVA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b, 3._dp, M_d2vdxdy_buv_b, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( DIVA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_buv_b, 1._dp, M_d2udy2_buv_b , DIFFERENT_NONZERO_PATTERN, perr)

    ! M_4_dudx_p_2_dvdy_buv_b
    CALL MatDuplicate( M_dudx_buv_b, MAT_COPY_VALUES, DIVA%M_4_dudx_p_2_dvdy_buv_b, perr)
    CALL MatDiagonalScale( DIVA%M_4_dudx_p_2_dvdy_buv_b, v_4, PETSC_NULL_VEC, perr)
    CALL MatAXPY( DIVA%M_4_dudx_p_2_dvdy_buv_b, 2._dp, M_dvdy_buv_b, DIFFERENT_NONZERO_PATTERN, perr)

    ! M_dudy_p_dvdx_buv_b
    CALL MatDuplicate( M_dudy_buv_b, MAT_COPY_VALUES, DIVA%M_dudy_p_dvdx_buv_b, perr)
    CALL MatAXPY( DIVA%M_dudy_p_dvdx_buv_b, 1._dp, M_dvdx_buv_b, DIFFERENT_NONZERO_PATTERN, perr)

    ! M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b
    CALL MatDuplicate( M_d2vdy2_buv_b, MAT_COPY_VALUES, DIVA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b, perr)
    CALL MatDiagonalScale( DIVA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b, v_4, PETSC_NULL_VEC, perr)
    CALL MatAXPY( DIVA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b, 3._dp, M_d2udxdy_buv_b, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( DIVA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_buv_b, 1._dp, M_d2vdx2_buv_b , DIFFERENT_NONZERO_PATTERN, perr)

    ! M_4_dvdy_p_2_dudx_buv_b
    CALL MatDuplicate( M_dvdy_buv_b, MAT_COPY_VALUES, DIVA%M_4_dvdy_p_2_dudx_buv_b, perr)
    CALL MatDiagonalScale( DIVA%M_4_dvdy_p_2_dudx_buv_b, v_4, PETSC_NULL_VEC, perr)
    CALL MatAXPY( DIVA%M_4_dvdy_p_2_dudx_buv_b, 2._dp, M_dudx_buv_b, DIFFERENT_NONZERO_PATTERN, perr)

    ! M_dvdx_p_dudy_buv_b
    CALL MatDuplicate( M_dvdx_buv_b, MAT_COPY_VALUES, DIVA%M_dvdx_p_dudy_buv_b, perr)
    CALL MatAXPY( DIVA%M_dvdx_p_dudy_buv_b, 1._dp, M_dudy_buv_b, DIFFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( M_dudx_buv_b   , perr)
    CALL MatDestroy( M_dudy_buv_b   , perr)
    CALL MatDestroy( M_d2udx2_buv_b , perr)
    CALL MatDestroy( M_d2udxdy_buv_b, perr)
    CALL MatDestroy( M_d2udy2_buv_b , perr)
    CALL MatDestroy( M_dvdx_buv_b   , perr)
    CALL MatDestroy( M_dvdy_buv_b   , perr)
    CALL MatDestroy( M_d2vdx2_buv_b , perr)
    CALL MatDestroy( M_d2vdxdy_buv_b, perr)
    CALL MatDestroy( M_d2vdy2_buv_b , perr)
    CALL VecDestroy( v_4            , perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_combined_mesh_operators_DIVA

END MODULE ice_velocity_DIVA_module
