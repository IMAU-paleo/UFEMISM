MODULE ice_velocity_BPA_module

! == Contains all the routines needed to solve the Blatter-Pattyn Approximation

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
  USE data_types_module,                   ONLY: type_mesh, type_ice_model, type_velocity_solver_BPA
  USE mesh_operators_module,               ONLY: calc_matrix_operators_x_y_z_3D, ddx_a_to_b_2D, ddy_a_to_b_2D, map_a_to_b_2D, &
                                                 field2vec_b, field2vec_ak
  USE petsc_module,                        ONLY: multiply_PETSc_matrix_with_vector_1D, double2petscmat, vec_double2petsc, &
                                                 solve_matrix_equation_PETSc
  USe basal_conditions_and_sliding_module, ONLY: calc_basal_conditions, calc_basal_friction_coefficient
  USE mesh_help_functions_module,          ONLY: find_ti_copy_ISMIP_HOM_periodic
  USE general_ice_model_data_module,       ONLY: determine_grounded_fractions

  IMPLICIT NONE

CONTAINS

! == The main routines, to be called from the ice_velocity_module

  SUBROUTINE initialise_BPA_solver( mesh, BPA)
    ! Initialise the BPA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'initialise_BPA_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory

    ! Solution
    CALL allocate_shared_dp_2D( mesh%nTri   , C%nz  , BPA%u_3D_b              , BPA%wu_3D_b              )
    CALL allocate_shared_dp_2D( mesh%nTri   , C%nz  , BPA%v_3D_b              , BPA%wv_3D_b              )

    ! Intermediate data fields
    CALL allocate_shared_dp_1D( mesh%nnbk   ,         BPA%u_bk_vec            , BPA%wu_bk_vec            )
    CALL allocate_shared_dp_1D( mesh%nnbk   ,         BPA%v_bk_vec            , BPA%wv_bk_vec            )
    CALL allocate_shared_dp_1D( mesh%nnak   ,         BPA%du_dx_ak_vec        , BPA%wdu_dx_ak_vec        )
    CALL allocate_shared_dp_1D( mesh%nnak   ,         BPA%du_dy_ak_vec        , BPA%wdu_dy_ak_vec        )
    CALL allocate_shared_dp_1D( mesh%nnak   ,         BPA%du_dz_ak_vec        , BPA%wdu_dz_ak_vec        )
    CALL allocate_shared_dp_1D( mesh%nnak   ,         BPA%dv_dx_ak_vec        , BPA%wdv_dx_ak_vec        )
    CALL allocate_shared_dp_1D( mesh%nnak   ,         BPA%dv_dy_ak_vec        , BPA%wdv_dy_ak_vec        )
    CALL allocate_shared_dp_1D( mesh%nnak   ,         BPA%dv_dz_ak_vec        , BPA%wdv_dz_ak_vec        )
    CALL allocate_shared_dp_1D( mesh%nnbks  ,         BPA%du_dx_bks_vec       , BPA%wdu_dx_bks_vec       )
    CALL allocate_shared_dp_1D( mesh%nnbks  ,         BPA%du_dy_bks_vec       , BPA%wdu_dy_bks_vec       )
    CALL allocate_shared_dp_1D( mesh%nnbks  ,         BPA%du_dz_bks_vec       , BPA%wdu_dz_bks_vec       )
    CALL allocate_shared_dp_1D( mesh%nnbks  ,         BPA%dv_dx_bks_vec       , BPA%wdv_dx_bks_vec       )
    CALL allocate_shared_dp_1D( mesh%nnbks  ,         BPA%dv_dy_bks_vec       , BPA%wdv_dy_bks_vec       )
    CALL allocate_shared_dp_1D( mesh%nnbks  ,         BPA%dv_dz_bks_vec       , BPA%wdv_dz_bks_vec       )
    CALL allocate_shared_dp_1D( mesh%nnak   ,         BPA%eta_ak_vec          , BPA%weta_ak_vec          )
    CALL allocate_shared_dp_1D( mesh%nnbks  ,         BPA%eta_bks_vec         , BPA%weta_bks_vec         )
    CALL allocate_shared_dp_1D( mesh%nnbk   ,         BPA%eta_bk_vec          , BPA%weta_bk_vec          )
    CALL allocate_shared_dp_1D( mesh%nnbk   ,         BPA%deta_dx_bk_vec      , BPA%wdeta_dx_bk_vec      )
    CALL allocate_shared_dp_1D( mesh%nnbk   ,         BPA%deta_dy_bk_vec      , BPA%wdeta_dy_bk_vec      )
    CALL allocate_shared_dp_1D( mesh%nTri   ,         BPA%beta_b_b            , BPA%wbeta_b_b            )
    CALL allocate_shared_dp_1D( mesh%nnb    ,         BPA%beta_b_b_vec        , BPA%wbeta_b_b_vec        )
    CALL allocate_shared_dp_1D( mesh%nnb    ,         BPA%taudx_b             , BPA%wtaudx_b             )
    CALL allocate_shared_dp_1D( mesh%nnb    ,         BPA%taudy_b             , BPA%wtaudy_b             )
    CALL allocate_shared_dp_1D( mesh%nnbk   ,         BPA%u_bk_prev_vec       , BPA%wu_bk_prev_vec       )
    CALL allocate_shared_dp_1D( mesh%nnbk   ,         BPA%v_bk_prev_vec       , BPA%wv_bk_prev_vec       )

    ! Parameters for the iterative solver used to solve the matrix equation representing the linearised BPA
    BPA%PETSc_rtol      = C%DIVA_PETSc_rtol
    BPA%PETSc_abstol    = C%DIVA_PETSc_abstol

    ! Load vectors
    CALL allocate_shared_dp_1D( mesh%nnbkuv  , BPA%bb               , BPA%wbb               )
    CALL allocate_shared_dp_1D( mesh%nnbkuv  , BPA%bb_free          , BPA%wbb_free          )
    CALL allocate_shared_dp_1D( mesh%nnbkuv  , BPA%bb_BC_surf       , BPA%wbb_BC_surf       )
    CALL allocate_shared_dp_1D( mesh%nnbkuv  , BPA%bb_BC_base       , BPA%wbb_BC_base       )
    CALL allocate_shared_dp_1D( mesh%nnbkuv  , BPA%bb_BC_west       , BPA%wbb_BC_west       )
    CALL allocate_shared_dp_1D( mesh%nnbkuv  , BPA%bb_BC_east       , BPA%wbb_BC_east       )
    CALL allocate_shared_dp_1D( mesh%nnbkuv  , BPA%bb_BC_south      , BPA%wbb_BC_south      )
    CALL allocate_shared_dp_1D( mesh%nnbkuv  , BPA%bb_BC_north      , BPA%wbb_BC_north      )
    CALL allocate_shared_dp_1D( mesh%nnbkuv  , BPA%bb_BC_prescr     , BPA%wbb_BC_prescr     )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BPA_solver

  SUBROUTINE solve_BPA( mesh, ice, BPA, BC_prescr_mask_3D_b, BC_prescr_u_3D_b, BC_prescr_v_3D_b)
    ! Calculate ice velocities by solving the Shallow Ice Approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)   , OPTIONAL :: BC_prescr_mask_3D_b      ! Mask of triangles where velocity is prescribed
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)   , OPTIONAL :: BC_prescr_u_3D_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)   , OPTIONAL :: BC_prescr_v_3D_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_BPA'
    INTEGER,  DIMENSION(:,:  ), POINTER                          :: BC_prescr_mask_3D_b_applied
    REAL(dp), DIMENSION(:,:  ), POINTER                          :: BC_prescr_u_3D_b_applied
    REAL(dp), DIMENSION(:,:  ), POINTER                          :: BC_prescr_v_3D_b_applied
    INTEGER                                                      :: wBC_prescr_mask_3D_b_applied, wBC_prescr_u_3D_b_applied, wBC_prescr_v_3D_b_applied
    INTEGER                                                      :: viscosity_iteration_i
    LOGICAL                                                      :: has_converged
    REAL(dp)                                                     :: resid_UV

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If there is no grounded ice, no way to solve the BPA
    IF (SUM( ice%mask_sheet_a) == 0) THEN
      BPA%u_3D_b = 0._dp
      BPA%v_3D_b = 0._dp
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Convert velocities from field form to vector form
    CALL convert_velocities_field2vec( mesh, BPA)

    ! Calculate all the 3-D matrix operators on the mesh
    CALL calc_matrix_operators_x_y_z_3D( mesh, ice)

    ! Handle the optional prescribed u,v boundary conditions
    CALL allocate_shared_int_2D( mesh%nTri, C%nz, BC_prescr_mask_3D_b_applied, wBC_prescr_mask_3D_b_applied)
    CALL allocate_shared_dp_2D(  mesh%nTri, C%nz, BC_prescr_u_3D_b_applied   , wBC_prescr_u_3D_b_applied   )
    CALL allocate_shared_dp_2D(  mesh%nTri, C%nz, BC_prescr_v_3D_b_applied   , wBC_prescr_v_3D_b_applied   )
    IF (PRESENT( BC_prescr_mask_3D_b) .OR. PRESENT( BC_prescr_u_3D_b) .OR. PRESENT( BC_prescr_v_3D_b)) THEN
      ! Safety
      IF (.NOT. (PRESENT(BC_prescr_mask_3D_b) .AND. PRESENT(BC_prescr_u_3D_b) .AND. PRESENT(BC_prescr_v_3D_b))) THEN
        CALL crash('need to provide prescribed u,v fields and mask!')
      END IF
      BC_prescr_mask_3D_b_applied( mesh%ti1:mesh%ti2,:) = BC_prescr_mask_3D_b( mesh%ti1:mesh%ti2,:)
      BC_prescr_u_3D_b_applied(    mesh%ti1:mesh%ti2,:) = BC_prescr_u_3D_b(    mesh%ti1:mesh%ti2,:)
      BC_prescr_v_3D_b_applied(    mesh%ti1:mesh%ti2,:) = BC_prescr_v_3D_b(    mesh%ti1:mesh%ti2,:)
    ELSE
      BC_prescr_mask_3D_b_applied( mesh%ti1:mesh%ti2,:) = 0
      BC_prescr_u_3D_b_applied(    mesh%ti1:mesh%ti2,:) = 0._dp
      BC_prescr_v_3D_b_applied(    mesh%ti1:mesh%ti2,:) = 0._dp
    END IF

    ! Calculate boundary conditions mask matrices
    CALL calc_BPA_BC_mask_matrices( mesh, BPA, BC_prescr_mask_3D_b_applied)

    ! Calculate the driving stress
    CALL calc_driving_stress( mesh, ice, BPA)

    ! Calculate the basal yield stress tau_c
    CALL calc_basal_conditions( mesh, ice)

    ! Determine sub-mesh grounded fractions for scaling the basal friction
    CALL determine_grounded_fractions( mesh, ice)

    ! Pre-calculate some combined operator matrices
    CALL calc_combined_mesh_operators_BPA( mesh, BPA)

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .FALSE.
    viscosity_iteration: DO WHILE (.NOT. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! Calculate the strain rates for the current velocity solution
      CALL calc_strain_rates( mesh, ice, BPA)

      ! Calculate the effective viscosity for the current velocity solution
      CALL calc_effective_viscosity( mesh, ice, BPA)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      CALL calc_applied_basal_friction_coefficient( mesh, ice, BPA)

      ! Solve the linearised BPA to calculate a new velocity solution
      CALL solve_BPA_linearised( mesh, ice, BPA, BC_prescr_u_3D_b_applied, BC_prescr_v_3D_b_applied)

      ! Limit velocities for improved stability
      CALL apply_velocity_limits( mesh, BPA)

      ! Reduce the change between velocity solutions
      CALL relax_viscosity_iterations( mesh, BPA, C%DIVA_visc_it_relax)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      CALL calc_visc_iter_UV_resid( mesh, BPA, resid_UV)

!      ! DENK DROM
!      IF (par%master) WRITE(0,*) '    BPA - viscosity iteration ', viscosity_iteration_i, ', u = [', MINVAL( BPA%u_bk_vec), ' - ', MAXVAL( BPA%u_bk_vec), '], resid = ', resid_UV

      ! If the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .FALSE.
      IF     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) THEN
        has_converged = .TRUE.
      ELSEIF (viscosity_iteration_i > C%DIVA_visc_it_nit) THEN
        has_converged = .TRUE.
      END IF

    END DO viscosity_iteration

    ! Convert velocities from vector form to field form
    CALL convert_velocities_vec2field( mesh, BPA)

    ! Clean up after yourself
    CALL deallocate_shared( wBC_prescr_mask_3D_b_applied)
    CALL deallocate_shared( wBC_prescr_u_3D_b_applied   )
    CALL deallocate_shared( wBC_prescr_v_3D_b_applied   )
    CALL MatDestroy( BPA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_bkuv_bk, perr)
    CALL MatDestroy( BPA%M_4_dudx_p_2_dvdy_bkuv_bk              , perr)
    CALL MatDestroy( BPA%M_dudy_p_dvdx_bkuv_bk                  , perr)
    CALL MatDestroy( BPA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_bkuv_bk, perr)
    CALL MatDestroy( BPA%M_4_dvdy_p_2_dudx_bkuv_bk              , perr)
    CALL MatDestroy( BPA%M_dvdx_p_dudy_bkuv_bk                  , perr)
    CALL MatDestroy( BPA%M_2_dudx_p_dvdy_bkuv_bk                , perr)
    CALL MatDestroy( BPA%M_2_dvdy_p_dudx_bkuv_bk                , perr)
    CALL MatDestroy( BPA%m_free                                 , perr)
    CALL MatDestroy( BPA%m_BC_surf                              , perr)
    CALL MatDestroy( BPA%m_BC_base                              , perr)
    CALL MatDestroy( BPA%m_BC_west                              , perr)
    CALL MatDestroy( BPA%m_BC_east                              , perr)
    CALL MatDestroy( BPA%m_BC_south                             , perr)
    CALL MatDestroy( BPA%m_BC_north                             , perr)
    CALL MatDestroy( BPA%m_BC_prescr                            , perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_BPA

  SUBROUTINE remap_BPA_solver( mesh_old, mesh_new, ice, BPA)
    ! Remap the BPA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh_old
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh_new
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'remap_BPA_solver'
    REAL(dp)                                                     :: dp_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dp_dummy = mesh_old%V( 1,1)
    dp_dummy = mesh_new%V( 1,1)
    dp_dummy = ice%Hi_a( 1)

    ! Reallocate shared memory

    ! Solution
    CALL reallocate_shared_dp_2D( mesh_new%nTri   , C%nz  , BPA%u_3D_b              , BPA%wu_3D_b              )
    CALL reallocate_shared_dp_2D( mesh_new%nTri   , C%nz  , BPA%v_3D_b              , BPA%wv_3D_b              )

    ! Intermediate data fields
    CALL reallocate_shared_dp_1D( mesh_new%nnbk   ,         BPA%u_bk_vec            , BPA%wu_bk_vec            )
    CALL reallocate_shared_dp_1D( mesh_new%nnbk   ,         BPA%v_bk_vec            , BPA%wv_bk_vec            )
    CALL reallocate_shared_dp_1D( mesh_new%nnak   ,         BPA%du_dx_ak_vec        , BPA%wdu_dx_ak_vec        )
    CALL reallocate_shared_dp_1D( mesh_new%nnak   ,         BPA%du_dy_ak_vec        , BPA%wdu_dy_ak_vec        )
    CALL reallocate_shared_dp_1D( mesh_new%nnak   ,         BPA%du_dz_ak_vec        , BPA%wdu_dz_ak_vec        )
    CALL reallocate_shared_dp_1D( mesh_new%nnak   ,         BPA%dv_dx_ak_vec        , BPA%wdv_dx_ak_vec        )
    CALL reallocate_shared_dp_1D( mesh_new%nnak   ,         BPA%dv_dy_ak_vec        , BPA%wdv_dy_ak_vec        )
    CALL reallocate_shared_dp_1D( mesh_new%nnak   ,         BPA%dv_dz_ak_vec        , BPA%wdv_dz_ak_vec        )
    CALL reallocate_shared_dp_1D( mesh_new%nnbks  ,         BPA%du_dx_bks_vec       , BPA%wdu_dx_bks_vec       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbks  ,         BPA%du_dy_bks_vec       , BPA%wdu_dy_bks_vec       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbks  ,         BPA%du_dz_bks_vec       , BPA%wdu_dz_bks_vec       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbks  ,         BPA%dv_dx_bks_vec       , BPA%wdv_dx_bks_vec       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbks  ,         BPA%dv_dy_bks_vec       , BPA%wdv_dy_bks_vec       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbks  ,         BPA%dv_dz_bks_vec       , BPA%wdv_dz_bks_vec       )
    CALL reallocate_shared_dp_1D( mesh_new%nnak   ,         BPA%eta_ak_vec          , BPA%weta_ak_vec          )
    CALL reallocate_shared_dp_1D( mesh_new%nnbks  ,         BPA%eta_bks_vec         , BPA%weta_bks_vec         )
    CALL reallocate_shared_dp_1D( mesh_new%nnbk   ,         BPA%eta_bk_vec          , BPA%weta_bk_vec          )
    CALL reallocate_shared_dp_1D( mesh_new%nnbk   ,         BPA%deta_dx_bk_vec      , BPA%wdeta_dx_bk_vec      )
    CALL reallocate_shared_dp_1D( mesh_new%nnbk   ,         BPA%deta_dy_bk_vec      , BPA%wdeta_dy_bk_vec      )
    CALL reallocate_shared_dp_1D( mesh_new%nTri   ,         BPA%beta_b_b            , BPA%wbeta_b_b            )
    CALL reallocate_shared_dp_1D( mesh_new%nnb    ,         BPA%beta_b_b_vec        , BPA%wbeta_b_b_vec        )
    CALL reallocate_shared_dp_1D( mesh_new%nnb    ,         BPA%taudx_b             , BPA%wtaudx_b             )
    CALL reallocate_shared_dp_1D( mesh_new%nnb    ,         BPA%taudy_b             , BPA%wtaudy_b             )
    CALL reallocate_shared_dp_1D( mesh_new%nnbk   ,         BPA%u_bk_prev_vec       , BPA%wu_bk_prev_vec       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbk   ,         BPA%v_bk_prev_vec       , BPA%wv_bk_prev_vec       )

    ! Parameters for the iterative solver used to solve the matrix equation representing the linearised BPA
    BPA%PETSc_rtol      = C%DIVA_PETSc_rtol
    BPA%PETSc_abstol    = C%DIVA_PETSc_abstol

    ! Load vectors
    CALL reallocate_shared_dp_1D( mesh_new%nnbkuv  , BPA%bb               , BPA%wbb               )
    CALL reallocate_shared_dp_1D( mesh_new%nnbkuv  , BPA%bb_free          , BPA%wbb_free          )
    CALL reallocate_shared_dp_1D( mesh_new%nnbkuv  , BPA%bb_BC_surf       , BPA%wbb_BC_surf       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbkuv  , BPA%bb_BC_base       , BPA%wbb_BC_base       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbkuv  , BPA%bb_BC_west       , BPA%wbb_BC_west       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbkuv  , BPA%bb_BC_east       , BPA%wbb_BC_east       )
    CALL reallocate_shared_dp_1D( mesh_new%nnbkuv  , BPA%bb_BC_south      , BPA%wbb_BC_south      )
    CALL reallocate_shared_dp_1D( mesh_new%nnbkuv  , BPA%bb_BC_north      , BPA%wbb_BC_north      )
    CALL reallocate_shared_dp_1D( mesh_new%nnbkuv  , BPA%bb_BC_prescr     , BPA%wbb_BC_prescr     )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_BPA_solver

! == Calculate several intermediate terms in the BPA

  SUBROUTINE calc_driving_stress( mesh, ice, BPA)
    ! Calculate the driving stress

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_driving_stress'
    REAL(dp), DIMENSION(:    ), POINTER                          :: dHs_dx_b
    REAL(dp), DIMENSION(:    ), POINTER                          :: dHs_dy_b
    INTEGER                                                      :: wdHs_dx_b, wdHs_dy_b
    INTEGER                                                      :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, dHs_dx_b, wdHs_dx_b)
    CALL allocate_shared_dp_1D( mesh%nTri, dHs_dy_b, wdHs_dy_b)

    ! Calculate dHs/dx, and dHs/dy on the b-grid
    CALL ddx_a_to_b_2D( mesh, ice%Hs_a, dHs_dx_b)
    CALL ddy_a_to_b_2D( mesh, ice%Hs_a, dHs_dy_b)

    ! Calculate the driving stress
    DO ti = mesh%ti1, mesh%ti2
      BPA%taudx_b( ti) = -ice_density * grav * dHs_dx_b( ti)
      BPA%taudy_b( ti) = -ice_density * grav * dHs_dy_b( ti)
    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx_b)
    CALL deallocate_shared( wdHs_dy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_driving_stress

  SUBROUTINE calc_strain_rates( mesh, ice, BPA)
    ! Calculate the strain rates
    !
    ! The velocities u and v are defined on the bk-grid (triangles, regular vertical)
    !
    ! The horizontal strain rates du/dx, du/dy, dv/dx, dv/dy are calculated on the ak-grid (vertices, regular vertical),
    ! and are then mapped to the bks-grid (triangles, staggered vertical)
    !
    ! The vertical shear strain rates du/dz, dv/dz are calculated on the bks-grid (triangles, staggered vertical),
    ! and are then mapped to the ak-grid (vertices, regular vertical)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_strain_rates'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate horizontal stretch/shear strain rates du/dx, du/dy, dv/dx, and dv/dy on the ak-grid
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_bk_ak, BPA%u_bk_vec, BPA%du_dx_ak_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_bk_ak, BPA%u_bk_vec, BPA%du_dy_ak_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_bk_ak, BPA%v_bk_vec, BPA%dv_dx_ak_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_bk_ak, BPA%v_bk_vec, BPA%dv_dy_ak_vec)

    ! Calculate vertical shear strain rates du/dz and dv/dz on the bks-grid
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddz_bk_bks, BPA%u_bk_vec, BPA%du_dz_bks_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddz_bk_bks, BPA%v_bk_vec, BPA%dv_dz_bks_vec)

    ! Map horizontal stretch/shear strain rates du/dx, du/dy, dv/dx, and dv/dy from the ak-grid to the bks-grid
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_ak_bks, BPA%du_dx_ak_vec, BPA%du_dx_bks_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_ak_bks, BPA%du_dy_ak_vec, BPA%du_dy_bks_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_ak_bks, BPA%dv_dx_ak_vec, BPA%dv_dx_bks_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_ak_bks, BPA%dv_dy_ak_vec, BPA%dv_dy_bks_vec)

    ! Map vertical shear strain rates du/dz and dv/dz from the bks-grid to the ak-grid
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_bks_ak, BPA%du_dz_bks_vec, BPA%du_dz_ak_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_bks_ak, BPA%dv_dz_bks_vec, BPA%dv_dz_ak_vec)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_strain_rates

  SUBROUTINE calc_effective_viscosity( mesh, ice, BPA)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradients of N
    !
    ! The effective viscosity eta is calculated separately on both the ak-grid (vertices, regular vertical)
    ! and on the bks-grid (triangles, staggered vertical), using the strain rates calculated in calc_strain_rates.
    !
    ! eta_bk, deta_dx_bk, and deta_dy_bk are calculated from eta_ak

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_effective_viscosity'
    REAL(dp), DIMENSION(:    ), POINTER                          ::  A_flow_ak_vec,  A_flow_bks_vec
    INTEGER                                                      :: wA_flow_ak_vec, wA_flow_bks_vec
    INTEGER                                                      :: vi,ti,k,ks,n
    REAL(dp)                                                     :: epsilon_sq

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nnak , A_flow_ak_vec , wA_flow_ak_vec )
    CALL allocate_shared_dp_1D( mesh%nnbks, A_flow_bks_vec, wA_flow_bks_vec)

    ! Convert ice flow factor on the ak-grid from field form to vector form
    CALL field2vec_ak( mesh, ice%A_flow_3D_a, A_flow_ak_vec)

    ! Map ice flow factor from the ak-grid to the bks-grid
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_ak_bks, A_flow_ak_vec, A_flow_bks_vec)

  ! == Calculate effective viscosity on the ak-grid

    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, C%nz

      n = mesh%vik2n( vi,k)

      ! Calculate the square of the effective strain rate epsilon
      epsilon_sq = BPA%du_dx_ak_vec( n)**2 + &
                   BPA%dv_dy_ak_vec( n)**2 + &
                   BPA%du_dx_ak_vec( n) * BPA%dv_dy_ak_vec( n) + &
                   0.25_dp * (BPA%du_dy_ak_vec( n) + BPA%dv_dx_ak_vec( n))**2 + &
                   0.25_dp * (BPA%du_dz_ak_vec( n)**2 + BPA%dv_dz_ak_vec( n)**2) + &
                   C%DIVA_epsilon_sq_0

      ! Calculate the effective viscosity eta
      BPA%eta_ak_vec( n) = 0.5_dp * A_flow_ak_vec( n)**(-1._dp/  C%n_flow) * (epsilon_sq)**((1._dp - C%n_flow)/(2._dp*C%n_flow))

      ! Safety
      BPA%eta_ak_vec( n) = MAX( BPA%eta_ak_vec( n), C%DIVA_visc_eff_min)

    END DO
    END DO

  ! == Calculate effective viscosity on the bks-grid

    DO ti = mesh%ti1, mesh%ti2
    DO ks = 1, C%nz-1

      n = mesh%tiks2n( ti,ks)

      ! Calculate the square of the effective strain rate epsilon
      epsilon_sq = BPA%du_dx_bks_vec( n)**2 + &
                   BPA%dv_dy_bks_vec( n)**2 + &
                   BPA%du_dx_bks_vec( n) * BPA%dv_dy_bks_vec( n) + &
                   0.25_dp * (BPA%du_dy_bks_vec( n) + BPA%dv_dx_bks_vec( n))**2 + &
                   0.25_dp * (BPA%du_dz_bks_vec( n)**2 + BPA%dv_dz_bks_vec( n)**2) + &
                   C%DIVA_epsilon_sq_0

      ! Calculate the effective viscosity eta
      BPA%eta_bks_vec( n) = 0.5_dp * A_flow_bks_vec( n)**(-1._dp/  C%n_flow) * (epsilon_sq)**((1._dp - C%n_flow)/(2._dp*C%n_flow))

      ! Safety
      BPA%eta_bks_vec( n) = MAX( BPA%eta_bks_vec( n), C%DIVA_visc_eff_min)

    END DO
    END DO

    ! Calculate eta and its gradients on the bk-grid
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_ak_bk, BPA%eta_ak_vec, BPA%eta_bk_vec    )
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_ak_bk, BPA%eta_ak_vec, BPA%deta_dx_bk_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_ak_bk, BPA%eta_ak_vec, BPA%deta_dy_bk_vec)

    ! Clean up after yourself
    CALL deallocate_shared( wA_flow_ak_vec )
    CALL deallocate_shared( wA_flow_bks_vec)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_effective_viscosity

  SUBROUTINE calc_applied_basal_friction_coefficient( mesh, ice, BPA)
    ! Calculate the applied basal friction coefficient beta_b, i.e. on the b-grid
    ! and scaled with the sub-grid grounded fraction

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(INOUT)           :: ice
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_applied_basal_friction_coefficient'
    REAL(dp), DIMENSION(:    ), POINTER                          ::  u_base_b,  v_base_b
    INTEGER                                                      :: wu_base_b, wv_base_b
    INTEGER                                                      :: ti,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, u_base_b, wu_base_b)
    CALL allocate_shared_dp_1D( mesh%nTri, v_base_b, wv_base_b)

    ! Obtain basal velocities
    DO ti = mesh%ti1, mesh%ti2
      n = mesh%tik2n( ti,C%nz)
      u_base_b( ti) = BPA%u_bk_vec( n)
      v_base_b( ti) = BPA%v_bk_vec( n)
    END DO

    ! Calculate the basal friction coefficient beta_b for the current velocity solution
    CALL calc_basal_friction_coefficient( mesh, ice, u_base_b, v_base_b)

    ! Map basal friction coefficient beta_b to the b-grid
    CALL map_a_to_b_2D( mesh, ice%beta_b_a, BPA%beta_b_b)

    ! Apply the sub-grid grounded fraction
    IF (C%do_GL_subgrid_friction) THEN
      DO ti = mesh%ti1, mesh%ti2
        BPA%beta_b_b( ti) = BPA%beta_b_b( ti) * ice%f_grnd_b( ti)**2
      END DO
      CALL sync
    END IF

    ! Convert from field form to vector form
    CALL field2vec_b( mesh, BPA%beta_b_b, BPA%beta_b_b_vec)

    ! Clean up after yourself
    CALL deallocate_shared( wu_base_b)
    CALL deallocate_shared( wv_base_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_applied_basal_friction_coefficient

! == Solve the linearised BPA

  SUBROUTINE solve_BPA_linearised( mesh, ice, BPA, BC_prescr_u_3D_b, BC_prescr_v_3D_b)
    ! Solve the linearised BPA

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)              :: BC_prescr_u_3D_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)              :: BC_prescr_v_3D_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'solve_BPA_linearised'
    TYPE(tMat)                                                   ::  AA_masked
    REAL(dp), DIMENSION(:    ), POINTER                          ::  bb_masked
    INTEGER                                                      :: wbb_masked
    INTEGER                                                      :: n1,n2
    REAL(dp), DIMENSION(:    ), POINTER                          ::  uv_buv
    INTEGER                                                      :: wuv_buv
    INTEGER                                                      :: ti,k,n,nu,nv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nnbkuv, bb_masked          , wbb_masked          )
    CALL allocate_shared_dp_1D( mesh%nnbkuv, uv_buv             , wuv_buv             )

    ! Calculate the stiffness matrix and load vector representing the BPA
    CALL calc_stiffness_matrix_BPA_free( mesh, BPA)

    ! Calculate the stiffness matrices and load vectors representing the different boundary conditions
    CALL calc_stiffness_matrix_BPA_BC_surf(   mesh, ice, BPA)
    CALL calc_stiffness_matrix_BPA_BC_base(   mesh, ice, BPA)
    CALL calc_stiffness_matrix_BPA_BC_west(   mesh,      BPA)
    CALL calc_stiffness_matrix_BPA_BC_east(   mesh,      BPA)
    CALL calc_stiffness_matrix_BPA_BC_south(  mesh,      BPA)
    CALL calc_stiffness_matrix_BPA_BC_north(  mesh,      BPA)
    CALL calc_stiffness_matrix_BPA_BC_prescr( mesh,      BPA, BC_prescr_u_3D_b, BC_prescr_v_3D_b)

    ! Use the boundary conditions mask matrices to combine the different stiffness matrices and load vectors

    ! Initialise
    CALL MatDuplicate( BPA%AA_free, MAT_SHARE_NONZERO_PATTERN, BPA%AA, perr)
    CALL partition_list( mesh%nnbkuv, par%i, par%n, n1, n2)
    BPA%bb( n1:n2) = 0._dp

    ! Free
    CALL MatMatMult(                           BPA%m_free, BPA%AA_free, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( BPA%m_free, BPA%bb_free,                                         bb_masked      )
    CALL MatAXPY( BPA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    BPA%bb( n1:n2) = BPA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! Surface
    CALL MatMatMult(                           BPA%m_BC_surf, BPA%AA_BC_surf, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( BPA%m_BC_surf, BPA%bb_BC_surf,                                         bb_masked      )
    CALL MatAXPY( BPA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    BPA%bb( n1:n2) = BPA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! Base
    CALL MatMatMult(                           BPA%m_BC_base, BPA%AA_BC_base, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( BPA%m_BC_base, BPA%bb_BC_base,                                         bb_masked      )
    CALL MatAXPY( BPA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    BPA%bb( n1:n2) = BPA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! West
    CALL MatMatMult(                           BPA%m_BC_west, BPA%AA_BC_west, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( BPA%m_BC_west, BPA%bb_BC_west,                                         bb_masked      )
    CALL MatAXPY( BPA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    BPA%bb( n1:n2) = BPA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! East
    CALL MatMatMult(                           BPA%m_BC_east, BPA%AA_BC_east, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( BPA%m_BC_east, BPA%bb_BC_east,                                         bb_masked      )
    CALL MatAXPY( BPA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    BPA%bb( n1:n2) = BPA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! South
    CALL MatMatMult(                           BPA%m_BC_south, BPA%AA_BC_south, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( BPA%m_BC_south, BPA%bb_BC_south,                                         bb_masked      )
    CALL MatAXPY( BPA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    BPA%bb( n1:n2) = BPA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! North
    CALL MatMatMult(                           BPA%m_BC_north, BPA%AA_BC_north, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( BPA%m_BC_north, BPA%bb_BC_north,                                         bb_masked      )
    CALL MatAXPY( BPA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    BPA%bb( n1:n2) = BPA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! Prescribed
    CALL MatMatMult(                           BPA%m_BC_prescr, BPA%AA_BC_prescr, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, AA_masked, perr)
    CALL multiply_PETSc_matrix_with_vector_1D( BPA%m_BC_prescr, BPA%bb_BC_prescr,                                         bb_masked      )
    CALL MatAXPY( BPA%AA, 1._dp, AA_masked, SUBSET_NONZERO_PATTERN, perr)
    BPA%bb( n1:n2) = BPA%bb( n1:n2) + bb_masked( n1:n2)
    CALL MatDestroy( AA_masked, perr)

    ! Combine the u and v velocities into a single vector;
    ! use the current velocity solution as the initial guess for the new one
    DO ti = mesh%ti1, mesh%ti2
    DO k  = 1, C%nz

      n  = mesh%tik2n(   ti,k  )
      nu = mesh%tikuv2n( ti,k,1)
      nv = mesh%tikuv2n( ti,k,2)

      uv_buv( nu) = BPA%u_bk_vec( n)
      uv_buv( nv) = BPA%v_bk_vec( n)

      ! Save the previous solution
      BPA%u_bk_prev_vec( n) = BPA%u_bk_vec( n)
      BPA%v_bk_prev_vec( n) = BPA%v_bk_vec( n)

    END DO
    END DO

    ! Solve the matrix equation
    CALL solve_matrix_equation_PETSc( BPA%AA, BPA%bb, uv_buv, BPA%PETSc_rtol, BPA%PETSc_abstol)

    ! Get the u and v velocities back in field form
    DO ti = mesh%ti1, mesh%ti2
    DO k  = 1, C%nz

      n  = mesh%tik2n(   ti,k  )
      nu = mesh%tikuv2n( ti,k,1)
      nv = mesh%tikuv2n( ti,k,2)

      BPA%u_bk_vec( n) = uv_buv( nu)
      BPA%v_bk_vec( n) = uv_buv( nv)

    END DO
    END DO

    ! Clean up after yourself
    CALL MatDestroy( BPA%AA          , perr)
    CALL MatDestroy( BPA%AA_free     , perr)
    CALL MatDestroy( BPA%AA_BC_surf  , perr)
    CALL MatDestroy( BPA%AA_BC_base  , perr)
    CALL MatDestroy( BPA%AA_BC_west  , perr)
    CALL MatDestroy( BPA%AA_BC_east  , perr)
    CALL MatDestroy( BPA%AA_BC_south , perr)
    CALL MatDestroy( BPA%AA_BC_north , perr)
    CALL MatDestroy( BPA%AA_BC_prescr, perr)
    CALL deallocate_shared( wbb_masked)
    CALL deallocate_shared( wuv_buv   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_BPA_linearised

  SUBROUTINE calc_stiffness_matrix_BPA_free( mesh, BPA)
    ! Calculate the stiffness matrix and load vector representing the BPA
    !
    ! The first equation of the BPA reads:
    !
    !   d/dx [ 2 eta ( 2 du/dx + dv/dy )] + d/dy [ eta ( du/dy + dv/dx ) ] + d/dz [ eta du/dz ] = -taud,x
    !
    ! By applying the product rule, this expands to:
    !
    !   4 eta d2u/dx2 + 4 deta/dx du/dx + 2 eta d2v/dxdy + 2 deta/dx dv/dy + ...
    !     eta d2u/dy2 +   deta/dy du/dy +   eta d2v/dxdy +   deta/dy dv/dx + d/dz [ eta du/dz] = -taud,x
    !
    ! Combining the terms involving N, dN/dx, and dN/dy yields:
    !
    !    eta    ( 4 d2u/dx2 + 3 d2v/dxdy + d2u/dy2 ) + ...
    !   deta/dx ( 4 du/dx   + 2 dv/dy ) + ...
    !   deta/dy (   du/dy   +   dv/dx ) + ...
    !   d/dz [ eta du/dz] = -taud,x
    !
    ! The equations are solved on the bk-grid (triangles, regular vertical). In the last term,
    ! the product rule is not explicitly applied; instead, du/dz and eta are defined on the bks-grid
    ! (triangles, staggered vertical).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_BPA_free'
    TYPE(tMat)                                                   :: M1
    TYPE(tMat)                                                   :: Au1, Au2, Au3, Au4, Au, Au_bkuv
    TYPE(tMat)                                                   :: Av1, Av2, Av3, Av4, Av, Av_bkuv
    TYPE(tVec)                                                   :: V1
    INTEGER                                                      :: ti,k,nu,nv

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Calculate the stiffness matrix for Eq. 1

    ! Au1 = 4 d2u/dx2 + 3 d2v/dxdy + d2u/dy2
    CALL MatDuplicate( BPA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_bkuv_bk, MAT_COPY_VALUES, Au1, perr)
    ! V1 = eta
    CALL vec_double2petsc( BPA%eta_bk_vec, V1)
    ! Au1 = eta ( 4 d2u/dx2 + 3 d2v/dxdy + d2u/dy2 )
    CALL MatDiagonalScale( Au1, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Au2 = 4 du/dx + 2 dv/dy
    CALL MatDuplicate( BPA%M_4_dudx_p_2_dvdy_bkuv_bk, MAT_COPY_VALUES, Au2, perr)
    ! V1 = deta/dx
    CALL vec_double2petsc( BPA%deta_dx_bk_vec, V1)
    ! Au2 = deta/dx ( 4 du/dx + 2 dv/dy )
    CALL MatDiagonalScale( Au2, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Au3 = du/dy + dv/dx
    CALL MatDuplicate( BPA%M_dudy_p_dvdx_bkuv_bk, MAT_COPY_VALUES, Au3, perr)
    ! V1 = deta/dy
    CALL vec_double2petsc( BPA%deta_dy_bk_vec, V1)
    ! Au3 = deta/dy ( du/dy + dv/dx )
    CALL MatDiagonalScale( Au3, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! M1 = du/dz on the bks-grid
    CALL MatMatMult( mesh%M_ddz_bk_bks, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M1, perr)
    ! V1 = eta on the bks-grid
    CALL vec_double2petsc( BPA%eta_bks_vec, V1)
    ! M1 = eta du/dz on the bks-grid
    CALL MatDiagonalScale( M1, V1, PETSC_NULL_VEC, perr)
    ! Au4 = d/dz ( eta du/dz)
    CALL MatMatMult( mesh%M_ddz_bks_bk, M1, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au4, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Au = Au1
    CALL MatDuplicate( Au1, MAT_COPY_VALUES, Au, perr)
    ! Au = Au + Au2 = Au1 + Au2
    CALL MatAXPY( Au, 1._dp, Au2, DIFFERENT_NONZERO_PATTERN, perr)
    ! Au = Au + Au3 = Au1 + Au2 + Au3
    CALL MatAXPY( Au, 1._dp, Au3, DIFFERENT_NONZERO_PATTERN, perr)
    ! Au = Au + Au4 = Au1 + Au2 + Au3 + Au4
    CALL MatAXPY( Au, 1._dp, Au4, DIFFERENT_NONZERO_PATTERN, perr)

    ! Au_bkuv = Au on the bkuv-grid
    CALL MatMatMult( mesh%M_map_bk_bku, Au, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au_bkuv, perr)

  ! == Calculate the stiffness matrix for Eq. 2

    ! Av1 = 4 d2v/dy2 + 3 d2u/dxdy + d2v/dx2
    CALL MatDuplicate( BPA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_bkuv_bk, MAT_COPY_VALUES, Av1, perr)
    ! V1 = eta
    CALL vec_double2petsc( BPA%eta_bk_vec, V1)
    ! Av1 = eta ( 4 d2v/dy2 + 3 d2u/dxdy + d2v/dx2 )
    CALL MatDiagonalScale( Av1, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Av2 = 4 dv/dy + 2 du/dx
    CALL MatDuplicate( BPA%M_4_dvdy_p_2_dudx_bkuv_bk, MAT_COPY_VALUES, Av2, perr)
    ! V1 = deta/dy
    CALL vec_double2petsc( BPA%deta_dy_bk_vec, V1)
    ! Av2 = deta/dy ( 4 dv/dy + 2 du/dx )
    CALL MatDiagonalScale( Av2, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Av3 = dv/dx + du/dy
    CALL MatDuplicate( BPA%M_dvdx_p_dudy_bkuv_bk, MAT_COPY_VALUES, Av3, perr)
    ! V1 = deta/dx
    CALL vec_double2petsc( BPA%deta_dx_bk_vec, V1)
    ! Av3 = deta/dx ( dv/dx + du/dy )
    CALL MatDiagonalScale( Av3, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! M1 = dv/dz on the bks-grid
    CALL MatMatMult( mesh%M_ddz_bk_bks, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M1, perr)
    ! V1 = eta on the bks-grid
    CALL vec_double2petsc( BPA%eta_bks_vec, V1)
    ! M1 = eta dv/dz on the bks-grid
    CALL MatDiagonalScale( M1, V1, PETSC_NULL_VEC, perr)
    ! Av4 = d/dz ( eta dv/dz)
    CALL MatMatMult( mesh%M_ddz_bks_bk, M1, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av4, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Av = Av1
    CALL MatDuplicate( Av1, MAT_COPY_VALUES, Av, perr)
    ! Av = Av + Av2 = Av1 + Av2
    CALL MatAXPY( Av, 1._dp, Av2, DIFFERENT_NONZERO_PATTERN, perr)
    ! Av = Av + Av3 = Av1 + Av2 + Av3
    CALL MatAXPY( Av, 1._dp, Av3, DIFFERENT_NONZERO_PATTERN, perr)
    ! Av = Av + Av4 = Av1 + Av2 + Av3 + Av4
    CALL MatAXPY( Av, 1._dp, Av4, DIFFERENT_NONZERO_PATTERN, perr)

    ! Av_bkuv = Av on the bkuv-grid
    CALL MatMatMult( mesh%M_map_bk_bkv, Av, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av_bkuv, perr)

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au_bkuv, MAT_COPY_VALUES, BPA%AA_free, perr)
    CALL MatAXPY( BPA%AA_free, 1._dp, Av_bkuv, DIFFERENT_NONZERO_PATTERN, perr)

  ! == Load vector

    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz
      nu = mesh%tikuv2n( ti,k,1)
      nv = mesh%tikuv2n( ti,k,2)
      BPA%bb_free( nu) = -BPA%taudx_b( ti)
      BPA%bb_free( nv) = -BPA%taudy_b( ti)
    END DO
    END DO

    ! Clean up after yourself
    CALL MatDestroy( M1     , perr)
    CALL MatDestroy( Au1    , perr)
    CALL MatDestroy( Au2    , perr)
    CALL MatDestroy( Au3    , perr)
    CALL MatDestroy( Au4    , perr)
    CALL MatDestroy( Au     , perr)
    CALL MatDestroy( Au_bkuv, perr)
    CALL MatDestroy( Av1    , perr)
    CALL MatDestroy( Av2    , perr)
    CALL MatDestroy( Av3    , perr)
    CALL MatDestroy( Av4    , perr)
    CALL MatDestroy( Av     , perr)
    CALL MatDestroy( Av_bkuv, perr)
    CALL VecDestroy( V1     , perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_BPA_free

  SUBROUTINE calc_stiffness_matrix_BPA_BC_surf( mesh, ice, BPA)
    ! Calculate the stiffness matrix and load vector representing
    ! surface boundary conditions to the BPA:
    !
    ! 2 dh/dx ( 2 du/dx + dv/dy ) + dh/dy ( du/dy + dv/dx ) - du/dz  = 0
    ! 2 dh/dy ( 2 dv/dy + du/dx ) + dh/dx ( dv/dx + du/dy ) - dv/dz  = 0

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_BPA_BC_surf'
    REAL(dp), DIMENSION(:    ), POINTER                          :: dHs_dx_b
    REAL(dp), DIMENSION(:    ), POINTER                          :: dHs_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                          :: dHs_dx_b_vec
    REAL(dp), DIMENSION(:    ), POINTER                          :: dHs_dy_b_vec
    REAL(dp), DIMENSION(:    ), POINTER                          :: dHs_dx_bk_vec
    REAL(dp), DIMENSION(:    ), POINTER                          :: dHs_dy_bk_vec
    INTEGER                                                      :: wdHs_dx_b, wdHs_dy_b, wdHs_dx_b_vec, wdHs_dy_b_vec, wdHs_dx_bk_vec, wdHs_dy_bk_vec
    INTEGER                                                      :: ti,k,n,nu,nv
    TYPE(tMat)                                                   :: Au1, Au2, Au3, Au, Au_bkuv
    TYPE(tMat)                                                   :: Av1, Av2, Av3, Av, Av_bkuv
    TYPE(tVec)                                                   :: V1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, dHs_dx_b     , wdHs_dx_b     )
    CALL allocate_shared_dp_1D( mesh%nTri, dHs_dy_b     , wdHs_dy_b     )
    CALL allocate_shared_dp_1D( mesh%nnb , dHs_dx_b_vec , wdHs_dx_b_vec )
    CALL allocate_shared_dp_1D( mesh%nnb , dHs_dy_b_vec , wdHs_dy_b_vec )
    CALL allocate_shared_dp_1D( mesh%nnbk, dHs_dx_bk_vec, wdHs_dx_bk_vec)
    CALL allocate_shared_dp_1D( mesh%nnbk, dHs_dy_bk_vec, wdHs_dy_bk_vec)

    ! Calculate surface slopes dHs/dx, dHs/dy on the b-grid
    CALL ddx_a_to_b_2D( mesh, ice%Hs_a, dHs_dx_b)
    CALL ddy_a_to_b_2D( mesh, ice%Hs_a, dHs_dy_b)

    ! Convert surface slopes from field form to vector form
    CALL field2vec_b( mesh, dHs_dx_b, dHs_dx_b_vec)
    CALL field2vec_b( mesh, dHs_dy_b, dHs_dy_b_vec)

    ! Map surface slopes to the surface layer of the bk-grid
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_b_bk_surf, dHs_dx_b_vec, dHs_dx_bk_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_b_bk_surf, dHs_dy_b_vec, dHs_dy_bk_vec)

  ! == Calculate the stiffness matrix for Eq. 1

    ! Au1 = (2 du/dx + dv/dy)
    CALL MatDuplicate( BPA%M_2_dudx_p_dvdy_bkuv_bk, MAT_COPY_VALUES, Au1, perr)
    ! V1 = dHs/dx
    CALL vec_double2petsc( dHs_dx_bk_vec, V1)
    ! Au1 = dHs/dx (2 du/dx + dv/dy)
    CALL MatDiagonalScale( Au1, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Au2 = (du/dy + dv/dx)
    CALL MatDuplicate( BPA%M_dudy_p_dvdx_bkuv_bk, MAT_COPY_VALUES, Au2, perr)
    ! V1 = dHs/dy
    CALL vec_double2petsc( dHs_dy_bk_vec, V1)
    ! Au2 = dHs/dy (du/dy + dv/dx)
    CALL MatDiagonalScale( Au2, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Au3 = du/dz
    CALL MatMatMult( mesh%M_ddz_bk_bk, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au3, perr)

    ! Au = Au1 = dHs/dx (2 du/dx + dv/dy)
    CALL MatDuplicate( Au1, MAT_COPY_VALUES, Au, perr)
    ! Au = Au + Au2 = dHs/dx (2 du/dx + dv/dy) + dHs/dy (du/dy + dv/dx)
    CALL MatAXPY( Au, 1._dp, Au2, DIFFERENT_NONZERO_PATTERN, perr)
    ! Au = Au - Au3 = dHs/dx (2 du/dx + dv/dy) + dHs/dy (du/dy + dv/dx) - du/dz
    CALL MatAXPY( Au, -1._dp, Au3, DIFFERENT_NONZERO_PATTERN, perr)

    ! Au_bkuv = Au on the bkuv-grid
    CALL MatMatMult( mesh%M_map_bk_bku, Au, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au_bkuv, perr)

  ! == Calculate the stiffness matrix for Eq. 2

    ! Av1 = (2 dv/dy + du/dx)
    CALL MatDuplicate( BPA%M_2_dvdy_p_dudx_bkuv_bk, MAT_COPY_VALUES, Av1, perr)
    ! V1 = dHs/dy
    CALL vec_double2petsc( dHs_dy_bk_vec, V1)
    ! Av1 = dHs/dy (2 dv/dy + du/dx)
    CALL MatDiagonalScale( Av1, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Av2 = (dv/dx + du/dy)
    CALL MatDuplicate( BPA%M_dvdx_p_dudy_bkuv_bk, MAT_COPY_VALUES, Av2, perr)
    ! V1 = dHs/dx
    CALL vec_double2petsc( dHs_dx_bk_vec, V1)
    ! Av2 = dHs/dx (dv/dx + du/dy)
    CALL MatDiagonalScale( Av2, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Av3 = dv/dz
    CALL MatMatMult( mesh%M_ddz_bk_bk, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av3, perr)

    ! Av = Av1 = dHs/dy (2 dv/dy + du/dx)
    CALL MatDuplicate( Av1, MAT_COPY_VALUES, Av, perr)
    ! Av = Av + Av2 = dHs/dy (2 dv/dy + du/dx) + dHs/dx (dv/dx + du/dy)
    CALL MatAXPY( Av, 1._dp, Av2, DIFFERENT_NONZERO_PATTERN, perr)
    ! Av = Av - Av3 = dHs/dy (2 dv/dy + du/dx) + dHs/dx (dv/dx + du/dy) - dv/dz
    CALL MatAXPY( Av, -1._dp, Av3, DIFFERENT_NONZERO_PATTERN, perr)

    ! Av_bkuv = Av on the bkuv-grid
    CALL MatMatMult( mesh%M_map_bk_bkv, Av, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av_bkuv, perr)

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au_bkuv, MAT_COPY_VALUES, BPA%AA_BC_surf, perr)
    CALL MatAXPY( BPA%AA_BC_surf, 1._dp, Av_bkuv, DIFFERENT_NONZERO_PATTERN, perr)

  ! == Load vector

    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz
      nu = mesh%tikuv2n( ti,k,1)
      nv = mesh%tikuv2n( ti,k,2)
      BPA%bb_BC_surf( nu) = 0._dp
      BPA%bb_BC_surf( nv) = 0._dp
    END DO
    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx_b     )
    CALL deallocate_shared( wdHs_dy_b     )
    CALL deallocate_shared( wdHs_dx_b_vec )
    CALL deallocate_shared( wdHs_dy_b_vec )
    CALL deallocate_shared( wdHs_dx_bk_vec)
    CALL deallocate_shared( wdHs_dy_bk_vec)
    CALL MatDestroy( Au1    , perr)
    CALL MatDestroy( Au2    , perr)
    CALL MatDestroy( Au3    , perr)
    CALL MatDestroy( Au     , perr)
    CALL MatDestroy( Au_bkuv, perr)
    CALL MatDestroy( Av1    , perr)
    CALL MatDestroy( Av2    , perr)
    CALL MatDestroy( Av3    , perr)
    CALL MatDestroy( Av     , perr)
    CALL MatDestroy( Av_bkuv, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_BPA_BC_surf

  SUBROUTINE calc_stiffness_matrix_BPA_BC_base( mesh, ice, BPA)
    ! Calculate the stiffness matrix and load vector representing
    ! basal boundary conditions to the BPA:
    !
    ! 2 db/dx ( 2 du/dx + dv/dy ) + db/dy ( du/dy + dv/dx ) - du/dz + beta_b/eta u = 0
    ! 2 db/dy ( 2 dv/dy + du/dx ) + db/dx ( dv/dx + du/dy ) - dv/dz + beta_b/eta v = 0

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_ice_model),                INTENT(IN)              :: ice
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_BPA_BC_base'
    REAL(dp), DIMENSION(:    ), POINTER                          :: b_a
    REAL(dp), DIMENSION(:    ), POINTER                          :: db_dx_b
    REAL(dp), DIMENSION(:    ), POINTER                          :: db_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                          :: db_dx_b_vec
    REAL(dp), DIMENSION(:    ), POINTER                          :: db_dy_b_vec
    REAL(dp), DIMENSION(:    ), POINTER                          :: db_dx_bk_vec
    REAL(dp), DIMENSION(:    ), POINTER                          :: db_dy_bk_vec
    INTEGER                                                      :: wb_a, wdb_dx_b, wdb_dy_b, wdb_dx_b_vec, wdb_dy_b_vec, wdb_dx_bk_vec, wdb_dy_bk_vec
    REAL(dp), DIMENSION(:    ), POINTER                          :: beta_b_over_eta_b
    REAL(dp), DIMENSION(:    ), POINTER                          :: beta_b_over_eta_b_vec
    REAL(dp), DIMENSION(:    ), POINTER                          :: beta_b_over_eta_bk_vec
    INTEGER                                                      :: wbeta_b_over_eta_b, wbeta_b_over_eta_b_vec, wbeta_b_over_eta_bk_vec
    INTEGER                                                      :: vi,ti,k,n,nu,nv
    TYPE(tMat)                                                   :: Au1, Au2, Au3, Au4, Au, Au_bkuv
    TYPE(tMat)                                                   :: Av1, Av2, Av3, Av4, Av, Av_bkuv
    TYPE(tVec)                                                   :: V1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Exception for the simple case of no sliding
    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! u_base = v_base = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bku, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au_bkuv, perr)
      CALL MatMatMult( mesh%M_map_bk_bkv, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av_bkuv, perr)
      CALL MatDuplicate( Au_bkuv, MAT_COPY_VALUES, BPA%AA_BC_base, perr)
      CALL MatAXPY( BPA%AA_BC_base, 1._dp, Av_bkuv, DIFFERENT_NONZERO_PATTERN, perr)
      CALL MatDestroy( Au_bkuv, perr)
      CALL MatDestroy( Av_bkuv, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nu = mesh%tikuv2n( ti,k,1)
        nv = mesh%tikuv2n( ti,k,2)
        BPA%bb_BC_base( nu) = 0._dp
        BPA%bb_BC_base( nv) = 0._dp
      END DO
      END DO

      ! Finalise routine path
      CALL finalise_routine( routine_name)

      RETURN

    END IF ! IF (C%choice_sliding_law == 'none') THEN

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV  , b_a                   , wb_a                   )
    CALL allocate_shared_dp_1D( mesh%nTri, db_dx_b               , wdb_dx_b               )
    CALL allocate_shared_dp_1D( mesh%nTri, db_dy_b               , wdb_dy_b               )
    CALL allocate_shared_dp_1D( mesh%nnb , db_dx_b_vec           , wdb_dx_b_vec           )
    CALL allocate_shared_dp_1D( mesh%nnb , db_dy_b_vec           , wdb_dy_b_vec           )
    CALL allocate_shared_dp_1D( mesh%nnbk, db_dx_bk_vec          , wdb_dx_bk_vec          )
    CALL allocate_shared_dp_1D( mesh%nnbk, db_dy_bk_vec          , wdb_dy_bk_vec          )
    CALL allocate_shared_dp_1D( mesh%nnb , beta_b_over_eta_b     , wbeta_b_over_eta_b     )
    CALL allocate_shared_dp_1D( mesh%nnb , beta_b_over_eta_b_vec , wbeta_b_over_eta_b_vec )
    CALL allocate_shared_dp_1D( mesh%nnbk, beta_b_over_eta_bk_vec, wbeta_b_over_eta_bk_vec)

    ! Calculate ice base elevation b = h - H
    DO vi = mesh%vi1, mesh%vi2
      b_a( vi) = ice%Hs_a( vi) - ice%Hi_a( vi)
    END DO

    ! Calculate beta_b / eta
    DO ti = mesh%ti1, mesh%ti2
      n = mesh%tik2n( ti,C%nz)
      beta_b_over_eta_b( ti) = BPA%beta_b_b( ti) / BPA%eta_bk_vec( n)
    END DO

    ! Calculate ice basal slopes db/dx, db/dy on the b-grid
    CALL ddx_a_to_b_2D( mesh, b_a, db_dx_b)
    CALL ddy_a_to_b_2D( mesh, b_a, db_dy_b)

    ! Convert ice basal slopes and beta_b / eta from field form to vector form
    CALL field2vec_b( mesh, db_dx_b          , db_dx_b_vec          )
    CALL field2vec_b( mesh, db_dy_b          , db_dy_b_vec          )
    CALL field2vec_b( mesh, beta_b_over_eta_b, beta_b_over_eta_b_vec)

    ! Map ice basal slopes to the surface layer of the bk-grid
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_b_bk_base, db_dx_b_vec          , db_dx_bk_vec          )
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_b_bk_base, db_dy_b_vec          , db_dy_bk_vec          )
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_b_bk_base, beta_b_over_eta_b_vec, beta_b_over_eta_bk_vec)

  ! == Calculate the stiffness matrix for Eq. 1

    ! Au1 = (2 du/dx + dv/dy)
    CALL MatDuplicate( BPA%M_2_dudx_p_dvdy_bkuv_bk, MAT_COPY_VALUES, Au1, perr)
    ! V1 = db/dx
    CALL vec_double2petsc( db_dx_bk_vec, V1)
    ! Au1 = db/dx (2 du/dx + dv/dy)
    CALL MatDiagonalScale( Au1, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Au2 = (du/dy + dv/dx)
    CALL MatDuplicate( BPA%M_dudy_p_dvdx_bkuv_bk, MAT_COPY_VALUES, Au2, perr)
    ! V1 = db/dy
    CALL vec_double2petsc( db_dy_bk_vec, V1)
    ! Au2 = db/dy (du/dy + dv/dx)
    CALL MatDiagonalScale( Au2, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Au3 = du/dz
    CALL MatMatMult( mesh%M_ddz_bk_bk, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au3, perr)

    ! Au4 = u
    CALL MatDuplicate( mesh%M_map_bku_bk, MAT_COPY_VALUES, Au4, perr)
    ! V1 = beta_b / eta
    CALL vec_double2petsc( beta_b_over_eta_bk_vec, V1)
    ! Au4 = beta_b / eta u
    CALL MatDiagonalScale( Au4, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Au = Au1 = db/dx (2 du/dx + dv/dy)
    CALL MatDuplicate( Au1, MAT_COPY_VALUES, Au, perr)
    ! Au = Au + Au2 = db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx)
    CALL MatAXPY( Au, 1._dp, Au2, DIFFERENT_NONZERO_PATTERN, perr)
    ! Au = Au - Au3 = db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) - du/dz
    CALL MatAXPY( Au, -1._dp, Au3, DIFFERENT_NONZERO_PATTERN, perr)
    ! Au = Au + Au4 = db/dx (2 du/dx + dv/dy) + db/dy (du/dy + dv/dx) - du/dz + beta_b / eta u
    CALL MatAXPY( Au, 1._dp, Au4, DIFFERENT_NONZERO_PATTERN, perr)

    ! Au_bkuv = Au on the bkuv-grid
    CALL MatMatMult( mesh%M_map_bk_bku, Au, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au_bkuv, perr)

  ! == Calculate the stiffness matrix for Eq. 2

    ! Av1 = (2 dv/dy + du/dx)
    CALL MatDuplicate( BPA%M_2_dvdy_p_dudx_bkuv_bk, MAT_COPY_VALUES, Av1, perr)
    ! V1 = db/dy
    CALL vec_double2petsc( db_dy_bk_vec, V1)
    ! Av1 = db/dy (2 dv/dy + du/dx)
    CALL MatDiagonalScale( Av1, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Av2 = (dv/dx + du/dy)
    CALL MatDuplicate( BPA%M_dvdx_p_dudy_bkuv_bk, MAT_COPY_VALUES, Av2, perr)
    ! V1 = db/dx
    CALL vec_double2petsc( db_dx_bk_vec, V1)
    ! Av2 = db/dx (dv/dx + du/dy)
    CALL MatDiagonalScale( Av2, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Av3 = dv/dz
    CALL MatMatMult( mesh%M_ddz_bk_bk, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av3, perr)

    ! Av4 = v
    CALL MatDuplicate( mesh%M_map_bkv_bk, MAT_COPY_VALUES, Av4, perr)
    ! V1 = beta_b / eta
    CALL vec_double2petsc( beta_b_over_eta_bk_vec, V1)
    ! Av4 = beta_b / eta v
    CALL MatDiagonalScale( Av4, V1, PETSC_NULL_VEC, perr)
    ! Clean up after yourself
    CALL VecDestroy( V1, perr)

    ! Av = Av1 = db/dy (2 dv/dy + du/dx)
    CALL MatDuplicate( Av1, MAT_COPY_VALUES, Av, perr)
    ! Av = Av + Av2 = db/dy (2 dv/dy + du/dx) + db/dx (dv/dx + du/dy)
    CALL MatAXPY( Av, 1._dp, Av2, DIFFERENT_NONZERO_PATTERN, perr)
    ! Av = Av - Av3 = db/dy (2 dv/dy + du/dx) + db/dx (dv/dx + du/dy) - dv/dz
    CALL MatAXPY( Av, -1._dp, Av3, DIFFERENT_NONZERO_PATTERN, perr)
    ! Av = Av + Av4 = db/dy (2 dv/dy + du/dx) + db/dx (dv/dx + du/dy) - dv/dz + beta_b / eta v
    CALL MatAXPY( Av, 1._dp, Av4, DIFFERENT_NONZERO_PATTERN, perr)

    ! Av_bkuv = Av on the bkuv-grid
    CALL MatMatMult( mesh%M_map_bk_bkv, Av, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av_bkuv, perr)

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au_bkuv, MAT_COPY_VALUES, BPA%AA_BC_base, perr)
    CALL MatAXPY( BPA%AA_BC_base, 1._dp, Av_bkuv, DIFFERENT_NONZERO_PATTERN, perr)

  ! == Load vector

    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz
      nu = mesh%tikuv2n( ti,k,1)
      nv = mesh%tikuv2n( ti,k,2)
      BPA%bb_BC_base( nu) = 0._dp
      BPA%bb_BC_base( nv) = 0._dp
    END DO
    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wb_a                   )
    CALL deallocate_shared( wdb_dx_b               )
    CALL deallocate_shared( wdb_dy_b               )
    CALL deallocate_shared( wdb_dx_b_vec           )
    CALL deallocate_shared( wdb_dy_b_vec           )
    CALL deallocate_shared( wdb_dx_bk_vec          )
    CALL deallocate_shared( wdb_dy_bk_vec          )
    CALL deallocate_shared( wbeta_b_over_eta_b     )
    CALL deallocate_shared( wbeta_b_over_eta_b_vec )
    CALL deallocate_shared( wbeta_b_over_eta_bk_vec)
    CALL MatDestroy( Au1    , perr)
    CALL MatDestroy( Au2    , perr)
    CALL MatDestroy( Au3    , perr)
    CALL MatDestroy( Au4    , perr)
    CALL MatDestroy( Au     , perr)
    CALL MatDestroy( Au_bkuv, perr)
    CALL MatDestroy( Av1    , perr)
    CALL MatDestroy( Av2    , perr)
    CALL MatDestroy( Av3    , perr)
    CALL MatDestroy( Av4    , perr)
    CALL MatDestroy( Av     , perr)
    CALL MatDestroy( Av_bkuv, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_BPA_BC_base

  SUBROUTINE calc_stiffness_matrix_BPA_BC_west( mesh, BPA)
    ! Calculate the stiffness matrix and load vector representing boundary conditions to the BPA on the western domain boundary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_BPA_BC_west'
    TYPE(tMat)                                                   :: M_dudx, Au
    TYPE(tMat)                                                   :: M_dvdx, Av
    INTEGER                                                      :: ti,k,nu,nv,ti_copy,n_copy

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == u

    IF     (C%DIVA_boundary_BC_u_west == 'infinite') THEN
      ! du/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddx_bk_bk, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx, perr)
      CALL MatMatMult( mesh%M_map_bk_bku, M_dudx           , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au    , perr)
      CALL MatDestroy( M_dudx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nu = mesh%tikuv2n( ti,k,1)
        BPA%bb_BC_west( nu) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_west == 'zero') THEN
      ! u = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bku, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nu = mesh%tikuv2n( ti,k,1)
        BPA%bb_BC_west( nu) = 0._dp
      END DO
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
      CALL MatMatMult( mesh%M_map_bk_bku, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        DO k = 1, C%nz
          nu = mesh%tikuv2n( ti,k,1)
          n_copy = mesh%tik2n( ti_copy,k)
          BPA%bb_BC_west( nu) = BPA%u_bk_vec( n_copy)
        END DO
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_u_west "' // TRIM( C%DIVA_boundary_BC_u_west) // '"!')
    END IF

  ! == v

    IF     (C%DIVA_boundary_BC_v_west == 'infinite') THEN
      ! dv/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddx_bk_bk, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx, perr)
      CALL MatMatMult( mesh%M_map_bk_bkv, M_dvdx           , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av    , perr)
      CALL MatDestroy( M_dvdx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nv = mesh%tikuv2n( ti,k,2)
        BPA%bb_BC_west( nv) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_west == 'zero') THEN
      ! v = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bkv, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nv = mesh%tikuv2n( ti,k,2)
        BPA%bb_BC_west( nv) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_west == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_v_west == 'periodic_ISMIP_HOM') THEN
      ! Do nothing here; instead, velocities are prescribed through the BC_prescr arguments

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bkv, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        DO k = 1, C%nz
          nv = mesh%tikuv2n( ti,k,2)
          n_copy = mesh%tik2n( ti_copy,k)
          BPA%bb_BC_west( nv) = BPA%v_bk_vec( n_copy)
        END DO
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_v_west "' // TRIM( C%DIVA_boundary_BC_v_west) // '"!')
    END IF

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au, MAT_SHARE_NONZERO_PATTERN, BPA%AA_BC_west, perr)
    CALL MatAXPY( BPA%AA_BC_west, 1._dp, Au, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( BPA%AA_BC_west, 1._dp, Av, UNKNOWN_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( Au, perr)
    CALL MatDestroy( Av, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_BPA_BC_west

  SUBROUTINE calc_stiffness_matrix_BPA_BC_east( mesh, BPA)
    ! Calculate the stiffness matrix and load vector representing boundary conditions to the BPA on the eastern domain boundary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_BPA_BC_east'
    TYPE(tMat)                                                   :: M_dudx, Au
    TYPE(tMat)                                                   :: M_dvdx, Av
    INTEGER                                                      :: ti,k,nu,nv,ti_copy,n_copy

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == u

    IF     (C%DIVA_boundary_BC_u_east == 'infinite') THEN
      ! du/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddx_bk_bk, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx, perr)
      CALL MatMatMult( mesh%M_map_bk_bku, M_dudx           , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au    , perr)
      CALL MatDestroy( M_dudx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nu = mesh%tikuv2n( ti,k,1)
        BPA%bb_BC_east( nu) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_east == 'zero') THEN
      ! u = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bku, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nu = mesh%tikuv2n( ti,k,1)
        BPA%bb_BC_east( nu) = 0._dp
      END DO
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
      CALL MatMatMult( mesh%M_map_bk_bku, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        DO k = 1, C%nz
          nu = mesh%tikuv2n( ti,k,1)
          n_copy = mesh%tik2n( ti_copy,k)
          BPA%bb_BC_east( nu) = BPA%u_bk_vec( n_copy)
        END DO
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_u_east "' // TRIM( C%DIVA_boundary_BC_u_east) // '"!')
    END IF

  ! == v

    IF     (C%DIVA_boundary_BC_v_east == 'infinite') THEN
      ! dv/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddx_bk_bk, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx, perr)
      CALL MatMatMult( mesh%M_map_bk_bkv, M_dvdx           , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av    , perr)
      CALL MatDestroy( M_dvdx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nv = mesh%tikuv2n( ti,k,2)
        BPA%bb_BC_east( nv) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_east == 'zero') THEN
      ! v = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bkv, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nv = mesh%tikuv2n( ti,k,2)
        BPA%bb_BC_east( nv) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_east == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_v_east == 'periodic_ISMIP_HOM') THEN
      ! Do nothing here; instead, velocities are prescribed through the BC_prescr arguments

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bkv, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        DO k = 1, C%nz
          nv = mesh%tikuv2n( ti,k,2)
          n_copy = mesh%tik2n( ti_copy,k)
          BPA%bb_BC_east( nv) = BPA%v_bk_vec( n_copy)
        END DO
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_v_east "' // TRIM( C%DIVA_boundary_BC_v_east) // '"!')
    END IF

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au, MAT_SHARE_NONZERO_PATTERN, BPA%AA_BC_east, perr)
    CALL MatAXPY( BPA%AA_BC_east, 1._dp, Au, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( BPA%AA_BC_east, 1._dp, Av, UNKNOWN_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( Au, perr)
    CALL MatDestroy( Av, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_BPA_BC_east

  SUBROUTINE calc_stiffness_matrix_BPA_BC_south( mesh, BPA)
    ! Calculate the stiffness matrix and load vector representing boundary conditions to the BPA on the southern domain boundary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_BPA_BC_south'
    TYPE(tMat)                                                   :: M_dudx, Au
    TYPE(tMat)                                                   :: M_dvdx, Av
    INTEGER                                                      :: ti,k,nu,nv,ti_copy,n_copy

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == u

    IF     (C%DIVA_boundary_BC_u_south == 'infinite') THEN
      ! du/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddy_bk_bk, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx, perr)
      CALL MatMatMult( mesh%M_map_bk_bku, M_dudx           , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au    , perr)
      CALL MatDestroy( M_dudx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nu = mesh%tikuv2n( ti,k,1)
        BPA%bb_BC_south( nu) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_south == 'zero') THEN
      ! u = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bku, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nu = mesh%tikuv2n( ti,k,1)
        BPA%bb_BC_south( nu) = 0._dp
      END DO
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
      CALL MatMatMult( mesh%M_map_bk_bku, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        DO k = 1, C%nz
          nu = mesh%tikuv2n( ti,k,1)
          n_copy = mesh%tik2n( ti_copy,k)
          BPA%bb_BC_south( nu) = BPA%u_bk_vec( n_copy)
        END DO
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_u_south "' // TRIM( C%DIVA_boundary_BC_u_south) // '"!')
    END IF

  ! == v

    IF     (C%DIVA_boundary_BC_v_south == 'infinite') THEN
      ! dv/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddy_bk_bk, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx, perr)
      CALL MatMatMult( mesh%M_map_bk_bkv, M_dvdx           , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av    , perr)
      CALL MatDestroy( M_dvdx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nv = mesh%tikuv2n( ti,k,2)
        BPA%bb_BC_south( nv) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_south == 'zero') THEN
      ! v = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bkv, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nv = mesh%tikuv2n( ti,k,2)
        BPA%bb_BC_south( nv) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_south == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_v_south == 'periodic_ISMIP_HOM') THEN
      ! Do nothing here; instead, velocities are prescribed through the BC_prescr arguments

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bkv, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        DO k = 1, C%nz
          nv = mesh%tikuv2n( ti,k,2)
          n_copy = mesh%tik2n( ti_copy,k)
          BPA%bb_BC_south( nv) = BPA%v_bk_vec( n_copy)
        END DO
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_v_south "' // TRIM( C%DIVA_boundary_BC_v_south) // '"!')
    END IF

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au, MAT_SHARE_NONZERO_PATTERN, BPA%AA_BC_south, perr)
    CALL MatAXPY( BPA%AA_BC_south, 1._dp, Au, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( BPA%AA_BC_south, 1._dp, Av, UNKNOWN_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( Au, perr)
    CALL MatDestroy( Av, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_BPA_BC_south

  SUBROUTINE calc_stiffness_matrix_BPA_BC_north( mesh, BPA)
    ! Calculate the stiffness matrix and load vector representing boundary conditions to the BPA on the northern domain boundary

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_BPA_BC_north'
    TYPE(tMat)                                                   :: M_dudx, Au
    TYPE(tMat)                                                   :: M_dvdx, Av
    INTEGER                                                      :: ti,k,nu,nv,ti_copy,n_copy

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == u

    IF     (C%DIVA_boundary_BC_u_north == 'infinite') THEN
      ! du/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddy_bk_bk, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx, perr)
      CALL MatMatMult( mesh%M_map_bk_bku, M_dudx           , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au    , perr)
      CALL MatDestroy( M_dudx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nu = mesh%tikuv2n( ti,k,1)
        BPA%bb_BC_north( nu) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_u_north == 'zero') THEN
      ! u = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bku, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nu = mesh%tikuv2n( ti,k,1)
        BPA%bb_BC_north( nu) = 0._dp
      END DO
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
      CALL MatMatMult( mesh%M_map_bk_bku, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        DO k = 1, C%nz
          nu = mesh%tikuv2n( ti,k,1)
          n_copy = mesh%tik2n( ti_copy,k)
          BPA%bb_BC_north( nu) = BPA%u_bk_vec( n_copy)
        END DO
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_u_north "' // TRIM( C%DIVA_boundary_BC_u_north) // '"!')
    END IF

  ! == v

    IF     (C%DIVA_boundary_BC_v_north == 'infinite') THEN
      ! dv/dx = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M2_ddy_bk_bk, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx, perr)
      CALL MatMatMult( mesh%M_map_bk_bkv, M_dvdx           , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av    , perr)
      CALL MatDestroy( M_dvdx, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nv = mesh%tikuv2n( ti,k,2)
        BPA%bb_BC_north( nv) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_north == 'zero') THEN
      ! v = 0

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bkv, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        nv = mesh%tikuv2n( ti,k,2)
        BPA%bb_BC_north( nv) = 0._dp
      END DO
      END DO

    ELSEIF (C%DIVA_boundary_BC_v_north == 'periodic') THEN
      ! No idea how to do this on a mesh...

      CALL crash('periodic boundary conditions not implemented in UFEMISM!')

    ELSEIF (C%DIVA_boundary_BC_v_north == 'periodic_ISMIP_HOM') THEN
      ! Do nothing here; instead, velocities are prescribed through the BC_prescr arguments

      ! Stiffness matrix
      CALL MatMatMult( mesh%M_map_bk_bkv, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

      ! Load vector
      DO ti = mesh%ti1, mesh%ti2
        CALL find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy)
        DO k = 1, C%nz
          nv = mesh%tikuv2n( ti,k,2)
          n_copy = mesh%tik2n( ti_copy,k)
          BPA%bb_BC_north( nv) = BPA%v_bk_vec( n_copy)
        END DO
      END DO

    ELSE
      CALL crash('unknown DIVA_boundary_BC_v_north "' // TRIM( C%DIVA_boundary_BC_v_north) // '"!')
    END IF

  ! == Combine Eqs. 1 and 2

    CALL MatDuplicate( Au, MAT_SHARE_NONZERO_PATTERN, BPA%AA_BC_north, perr)
    CALL MatAXPY( BPA%AA_BC_north, 1._dp, Au, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( BPA%AA_BC_north, 1._dp, Av, UNKNOWN_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( Au, perr)
    CALL MatDestroy( Av, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_BPA_BC_north

  SUBROUTINE calc_stiffness_matrix_BPA_BC_prescr( mesh, BPA, BC_prescr_u_3D_b, BC_prescr_v_3D_b)
    ! Calculate the stiffness matrix and load vector representing boundary conditions in the form of directly prescribed velocities

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)              :: BC_prescr_u_3D_b         ! Prescribed velocities in the x-direction
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)              :: BC_prescr_v_3D_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_stiffness_matrix_BPA_BC_prescr'
    TYPE(tMat)                                                   :: Au, Av
    INTEGER                                                      :: ti,k,nu,nv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Stiffness matrices
    CALL MatMatMult( mesh%M_map_bk_bku, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Au, perr)
    CALL MatMatMult( mesh%M_map_bk_bkv, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, Av, perr)

    ! Combine Eqs. 1 and 2
    CALL MatDuplicate( Au, MAT_SHARE_NONZERO_PATTERN, BPA%AA_BC_prescr, perr)
    CALL MatAXPY( BPA%AA_BC_prescr, 1._dp, Au, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( BPA%AA_BC_prescr, 1._dp, Av, UNKNOWN_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( Au, perr)
    CALL MatDestroy( Av, perr)

    ! Load vector
    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz
      nu = mesh%tikuv2n( ti,k,1)
      nv = mesh%tikuv2n( ti,k,2)
      BPA%bb_BC_prescr( nu) = BC_prescr_u_3D_b( ti,k)
      BPA%bb_BC_prescr( nv) = BC_prescr_v_3D_b( ti,k)
    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_stiffness_matrix_BPA_BC_prescr

  SUBROUTINE calc_BPA_BC_mask_matrices( mesh, BPA, BC_prescr_mask_3D_b)
    ! Calculate boundary conditions mask matrices

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: BC_prescr_mask_3D_b      ! Mask of 3-D triangles where velocity is prescribed

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_BPA_BC_mask_matrices'
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_free
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_surf
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_base
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_west
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_east
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_north
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_south
    REAL(dp), DIMENSION(:    ), POINTER                          :: m_BC_prescr
    INTEGER                                                      :: wm_free, wm_BC_surf, wm_BC_base, wm_BC_west, wm_BC_east, wm_BC_north, wm_BC_south, wm_BC_prescr
    INTEGER                                                      :: ti,k,nu,nv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nnbkuv, m_free     , wm_free     )
    CALL allocate_shared_dp_1D( mesh%nnbkuv, m_BC_surf  , wm_BC_surf  )
    CALL allocate_shared_dp_1D( mesh%nnbkuv, m_BC_base  , wm_BC_base  )
    CALL allocate_shared_dp_1D( mesh%nnbkuv, m_BC_west  , wm_BC_west  )
    CALL allocate_shared_dp_1D( mesh%nnbkuv, m_BC_east  , wm_BC_east  )
    CALL allocate_shared_dp_1D( mesh%nnbkuv, m_BC_north , wm_BC_north )
    CALL allocate_shared_dp_1D( mesh%nnbkuv, m_BC_south , wm_BC_south )
    CALL allocate_shared_dp_1D( mesh%nnbkuv, m_BC_prescr, wm_BC_prescr)

    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz

      nu = mesh%tikuv2n( ti,k,1)
      nv = mesh%tikuv2n( ti,k,2)

      ! Initialise everything with zero
      m_free(      nu) = 0._dp
      m_BC_surf(   nu) = 0._dp
      m_BC_base(   nu) = 0._dp
      m_BC_west(   nu) = 0._dp
      m_BC_east(   nu) = 0._dp
      m_BC_south(  nu) = 0._dp
      m_BC_north(  nu) = 0._dp
      m_BC_prescr( nu) = 0._dp

      m_free(      nv) = 0._dp
      m_BC_surf(   nv) = 0._dp
      m_BC_base(   nv) = 0._dp
      m_BC_west(   nv) = 0._dp
      m_BC_east(   nv) = 0._dp
      m_BC_south(  nv) = 0._dp
      m_BC_north(  nv) = 0._dp
      m_BC_prescr( nv) = 0._dp

      IF     (BC_prescr_mask_3D_b( ti,k) == 1) THEN
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

      ELSEIF (k == 1) THEN
        ! Apply surface boundary conditions

        m_BC_surf( nu) = 1._dp
        m_BC_surf( nv) = 1._dp

      ELSEIF (k == C%nz) THEN
        ! Apply basal boundary conditions

        m_BC_base( nu) = 1._dp
        m_BC_base( nv) = 1._dp

      ELSE
        ! No boundary conditions apply; solve the BPA instead

        m_free( nu) = 1._dp
        m_free( nv) = 1._dp

      END IF

    END DO ! DO k = 1, C%nz
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync

    ! Convert to diagonal matrices
    CALL double2petscmat( m_free     , BPA%m_free     )
    CALL double2petscmat( m_BC_surf  , BPA%m_BC_surf  )
    CALL double2petscmat( m_BC_base  , BPA%m_BC_base  )
    CALL double2petscmat( m_BC_west  , BPA%m_BC_west  )
    CALL double2petscmat( m_BC_east  , BPA%m_BC_east  )
    CALL double2petscmat( m_BC_south , BPA%m_BC_south )
    CALL double2petscmat( m_BC_north , BPA%m_BC_north )
    CALL double2petscmat( m_BC_prescr, BPA%m_BC_prescr)

    ! Clean up after yourself
    CALL deallocate_shared( wm_free     )
    CALL deallocate_shared( wm_BC_surf  )
    CALL deallocate_shared( wm_BC_base  )
    CALL deallocate_shared( wm_BC_west  )
    CALL deallocate_shared( wm_BC_east  )
    CALL deallocate_shared( wm_BC_north )
    CALL deallocate_shared( wm_BC_south )
    CALL deallocate_shared( wm_BC_prescr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_BPA_BC_mask_matrices

! == Some useful tools for improving numerical stability of the viscosity iteration

  SUBROUTINE relax_viscosity_iterations( mesh, BPA, R)
    ! Reduce the change between velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA
    REAL(dp),                            INTENT(IN)              :: R

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'relax_viscosity_iterations'
    INTEGER                                                      :: ti,k,n

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz
      n = mesh%tik2n( ti,k)
      BPA%u_bk_vec( n) = (R * BPA%u_bk_vec( n)) + ((1._dp - R) * BPA%u_bk_prev_vec( n))
      BPA%v_bk_vec( n) = (R * BPA%v_bk_vec( n)) + ((1._dp - R) * BPA%v_bk_prev_vec( n))
    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE relax_viscosity_iterations

  SUBROUTINE calc_visc_iter_UV_resid( mesh, BPA, resid_UV)
    ! Calculate the L2-norm of the two consecutive velocity solutions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(IN)              :: BPA
    REAL(dp),                            INTENT(OUT)             :: resid_UV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_visc_iter_UV_resid'
    INTEGER                                                      :: ierr
    INTEGER                                                      :: ti,k,n
    REAL(dp)                                                     :: res1, res2

    ! Add routine to path
    CALL init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz

      n = mesh%tik2n( ti,k)

      res1 = res1 + (BPA%u_bk_vec( n) - BPA%u_bk_prev_vec( n))**2
      res1 = res1 + (BPA%v_bk_vec( n) - BPA%v_bk_prev_vec( n))**2

      res2 = res2 + (BPA%u_bk_vec( n) + BPA%u_bk_prev_vec( n))**2
      res2 = res2 + (BPA%v_bk_vec( n) + BPA%v_bk_prev_vec( n))**2

    END DO
    END DO

    ! Combine results from all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate residual
    resid_UV = 2._dp * res1 / MAX( res2, 1E-8_dp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_visc_iter_UV_resid

  SUBROUTINE apply_velocity_limits( mesh, BPA)
    ! Limit velocities for improved stability

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'apply_velocity_limits'
    INTEGER                                                      :: ti,k,n
    REAL(dp)                                                     :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz

      n = mesh%tik2n( ti,k)

      ! Calculate absolute speed
      uabs = SQRT( BPA%u_bk_vec( n)**2 + BPA%v_bk_vec( n)**2)

      ! Reduce velocities if necessary
      IF (uabs > C%DIVA_vel_max) THEN
        BPA%u_bk_vec( n) = BPA%u_bk_vec( n) * C%DIVA_vel_max / uabs
        BPA%v_bk_vec( n) = BPA%v_bk_vec( n) * C%DIVA_vel_max / uabs
      END IF

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_velocity_limits

! == Calculate combined mesh operators for efficient calculation of the stiffness matrix

  SUBROUTINE calc_combined_mesh_operators_BPA( mesh, BPA)
    ! Calculate combined mesh operators for efficient calculation of the stiffness matrix
    !
    ! NOTE: unlike the SSA and the DIVA, these operators depend on ice geometry,
    !       so they need to be recalculated in every time step (but not in every
    !       viscosity iteration, which is why its convenient to store them separately)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_combined_mesh_operators_BPA'
    TYPE(tMat)                                                   :: M_dudx_bkuv_bk
    TYPE(tMat)                                                   :: M_dudy_bkuv_bk
    TYPE(tMat)                                                   :: M_d2udx2_bkuv_bk
    TYPE(tMat)                                                   :: M_d2udxdy_bkuv_bk
    TYPE(tMat)                                                   :: M_d2udy2_bkuv_bk
    TYPE(tMat)                                                   :: M_dvdx_bkuv_bk
    TYPE(tMat)                                                   :: M_dvdy_bkuv_bk
    TYPE(tMat)                                                   :: M_d2vdx2_bkuv_bk
    TYPE(tMat)                                                   :: M_d2vdxdy_bkuv_bk
    TYPE(tMat)                                                   :: M_d2vdy2_bkuv_bk

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Partial derivatives of u and v on the bk-grid
    CALL MatMatMult( mesh%M2_ddx_bk_bk   , mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudx_bkuv_bk   , perr)
    CALL MatMatMult( mesh%M2_ddy_bk_bk   , mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dudy_bkuv_bk   , perr)
    CALL MatMatMult( mesh%M2_d2dx2_bk_bk , mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2udx2_bkuv_bk , perr)
    CALL MatMatMult( mesh%M2_d2dxdy_bk_bk, mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2udxdy_bkuv_bk, perr)
    CALL MatMatMult( mesh%M2_d2dy2_bk_bk , mesh%M_map_bku_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2udy2_bkuv_bk , perr)

    CALL MatMatMult( mesh%M2_ddx_bk_bk   , mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdx_bkuv_bk   , perr)
    CALL MatMatMult( mesh%M2_ddy_bk_bk   , mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_dvdy_bkuv_bk   , perr)
    CALL MatMatMult( mesh%M2_d2dx2_bk_bk , mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2vdx2_bkuv_bk , perr)
    CALL MatMatMult( mesh%M2_d2dxdy_bk_bk, mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2vdxdy_bkuv_bk, perr)
    CALL MatMatMult( mesh%M2_d2dy2_bk_bk , mesh%M_map_bkv_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_d2vdy2_bkuv_bk , perr)

    ! 4 d2u/dx2 + 3 d2v/dxdy + d2u/dy2
    CALL MatDuplicate( M_d2udy2_bkuv_bk, MAT_COPY_VALUES, BPA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_bkuv_bk, perr)
    CALL MatAXPY( BPA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_bkuv_bk, 4._dp, M_d2udx2_bkuv_bk , DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( BPA%M_4_d2udx2_p_3_d2vdxdy_p_d2udy2_bkuv_bk, 3._dp, M_d2vdxdy_bkuv_bk, DIFFERENT_NONZERO_PATTERN, perr)

    ! 4 du/dx + 2 dv/dy
    CALL MatDuplicate( M_dvdy_bkuv_bk, MAT_COPY_VALUES, BPA%M_4_dudx_p_2_dvdy_bkuv_bk, perr)
    CALL MatAXPY( BPA%M_4_dudx_p_2_dvdy_bkuv_bk, 1._dp, M_dvdy_bkuv_bk , DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( BPA%M_4_dudx_p_2_dvdy_bkuv_bk, 4._dp, M_dudx_bkuv_bk , DIFFERENT_NONZERO_PATTERN, perr)

    ! du/dy + dv/dx
    CALL MatDuplicate( M_dudy_bkuv_bk, MAT_COPY_VALUES, BPA%M_dudy_p_dvdx_bkuv_bk, perr)
    CALL MatAXPY( BPA%M_dudy_p_dvdx_bkuv_bk, 1._dp, M_dvdx_bkuv_bk, DIFFERENT_NONZERO_PATTERN, perr)

    ! 4 d2v/dy2 + 3 d2u/dxdy + d2v/dx2
    CALL MatDuplicate( M_d2vdx2_bkuv_bk, MAT_COPY_VALUES, BPA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_bkuv_bk, perr)
    CALL MatAXPY( BPA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_bkuv_bk, 4._dp, M_d2vdy2_bkuv_bk , DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( BPA%M_4_d2vdy2_p_3_d2udxdy_p_d2vdx2_bkuv_bk, 3._dp, M_d2udxdy_bkuv_bk, DIFFERENT_NONZERO_PATTERN, perr)

    ! 4 dv/dy + 2 du/dx
    CALL MatDuplicate( M_dudx_bkuv_bk, MAT_COPY_VALUES, BPA%M_4_dvdy_p_2_dudx_bkuv_bk, perr)
    CALL MatAXPY( BPA%M_4_dvdy_p_2_dudx_bkuv_bk, 1._dp, M_dudx_bkuv_bk , DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( BPA%M_4_dvdy_p_2_dudx_bkuv_bk, 4._dp, M_dvdy_bkuv_bk , DIFFERENT_NONZERO_PATTERN, perr)

    ! dv/dx + du/dy
    CALL MatDuplicate( M_dvdx_bkuv_bk, MAT_COPY_VALUES, BPA%M_dvdx_p_dudy_bkuv_bk, perr)
    CALL MatAXPY( BPA%M_dvdx_p_dudy_bkuv_bk, 1._dp, M_dudy_bkuv_bk, DIFFERENT_NONZERO_PATTERN, perr)

    ! 2 du/dx + dv/dy
    CALL MatDuplicate( M_dvdy_bkuv_bk, MAT_COPY_VALUES, BPA%M_2_dudx_p_dvdy_bkuv_bk, perr)
    CALL MatAXPY( BPA%M_2_dudx_p_dvdy_bkuv_bk, 2._dp, M_dudx_bkuv_bk, DIFFERENT_NONZERO_PATTERN, perr)

    ! 2 dv/dy + du/dx
    CALL MatDuplicate( M_dudx_bkuv_bk, MAT_COPY_VALUES, BPA%M_2_dvdy_p_dudx_bkuv_bk, perr)
    CALL MatAXPY( BPA%M_2_dvdy_p_dudx_bkuv_bk, 2._dp, M_dvdy_bkuv_bk, DIFFERENT_NONZERO_PATTERN, perr)

    ! Clean up after yourself
    CALL MatDestroy( M_dudx_bkuv_bk   , perr)
    CALL MatDestroy( M_dudy_bkuv_bk   , perr)
    CALL MatDestroy( M_d2udx2_bkuv_bk , perr)
    CALL MatDestroy( M_d2udxdy_bkuv_bk, perr)
    CALL MatDestroy( M_d2udy2_bkuv_bk , perr)
    CALL MatDestroy( M_dvdx_bkuv_bk   , perr)
    CALL MatDestroy( M_dvdy_bkuv_bk   , perr)
    CALL MatDestroy( M_d2vdx2_bkuv_bk , perr)
    CALL MatDestroy( M_d2vdxdy_bkuv_bk, perr)
    CALL MatDestroy( M_d2vdy2_bkuv_bk , perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_combined_mesh_operators_BPA

! == Convert final velocity solution between field form and vector form

  SUBROUTINE convert_velocities_field2vec( mesh, BPA)
    ! Convert velocities from field form to vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'convert_velocities_field2vec'
    INTEGER                                                      :: ti,k,n

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz

      n = mesh%tik2n( ti,k)

      BPA%u_bk_vec( n) = BPA%u_3D_b( ti,k)
      BPA%v_bk_vec( n) = BPA%v_3D_b( ti,k)

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE convert_velocities_field2vec

  SUBROUTINE convert_velocities_vec2field( mesh, BPA)
    ! Convert velocities from vector form to field form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    TYPE(type_velocity_solver_BPA),      INTENT(INOUT)           :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'convert_velocities_vec2field'
    INTEGER                                                      :: ti,k,n

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz

      n = mesh%tik2n( ti,k)

      BPA%u_3D_b( ti,k) = BPA%u_bk_vec( n)
      BPA%v_3D_b( ti,k) = BPA%v_bk_vec( n)

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE convert_velocities_vec2field

END MODULE ice_velocity_BPA_module
