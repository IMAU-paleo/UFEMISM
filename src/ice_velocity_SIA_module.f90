MODULE ice_velocity_SIA_module

! == Contains all the routines needed to solve the Shallow Ice Approximation

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
  USE data_types_module,                   ONLY: type_mesh, type_ice_model, type_velocity_solver_SIA
  USE mesh_operators_module,               ONLY: ddx_a_to_a_2D, ddy_a_to_a_2D, map_a_to_b_2D, map_a_to_b_3D, ddx_a_to_b_2D, ddy_a_to_b_2D
  USE utilities_module,                    ONLY: vertical_integration_from_bottom_to_zeta

  IMPLICIT NONE

CONTAINS

! == The main routines, to be called from the ice_velocity_module

  SUBROUTINE initialise_SIA_solver( mesh, SIA)
    ! Initialise the SIA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_velocity_solver_SIA),      INTENT(INOUT) :: SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_SIA_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory

    ! Solution
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, SIA%u_3D_b    , SIA%wu_3D_b    )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, SIA%v_3D_b    , SIA%wv_3D_b    )
    CALL allocate_shared_dp_2D( mesh%nV  , C%nz, SIA%du_dz_3D_a, SIA%wdu_dz_3D_a)
    CALL allocate_shared_dp_2D( mesh%nV  , C%nz, SIA%dv_dz_3D_a, SIA%wdv_dz_3D_a)

    ! Intermediate data fields
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, SIA%D_3D_b    , SIA%wD_3D_b    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SIA_solver

  SUBROUTINE solve_SIA( mesh, ice, SIA)
    ! Calculate ice velocities by solving the Shallow Ice Approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_velocity_solver_SIA),      INTENT(INOUT) :: SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_SIA'
    REAL(dp), PARAMETER                                :: D_SIA_limit = -1E5_dp
    REAL(dp), DIMENSION(:    ), POINTER                ::  Hi_b,  Hs_b,  dHs_dx_a,  dHs_dy_a,  dHs_dx_b,  dHs_dy_b
    INTEGER                                            :: wHi_b, wHs_b, wdHs_dx_a, wdHs_dy_a, wdHs_dx_b, wdHs_dy_b
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  A_flow_3D_b
    INTEGER                                            :: wA_flow_3D_b
    INTEGER                                            :: vi,ti,k
    REAL(dp)                                           :: abs_grad_Hs
    REAL(dp), DIMENSION(C%nz)                          :: z, int_A_hminzetan

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri,       Hi_b       , wHi_b       )
    CALL allocate_shared_dp_1D( mesh%nTri,       Hs_b       , wHs_b       )
    CALL allocate_shared_dp_1D( mesh%nV  ,       dHs_dx_a   , wdHs_dx_a   )
    CALL allocate_shared_dp_1D( mesh%nV  ,       dHs_dy_a   , wdHs_dy_a   )
    CALL allocate_shared_dp_1D( mesh%nTri,       dHs_dx_b   , wdHs_dx_b   )
    CALL allocate_shared_dp_1D( mesh%nTri,       dHs_dy_b   , wdHs_dy_b   )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, A_flow_3D_b, wA_flow_3D_b)

    ! Calculate ice thickness, surface elevation, surface slopes, and ice flow factor on the b-grid
    CALL map_a_to_b_2D( mesh, ice%Hi_a       , Hi_b       )
    CALL map_a_to_b_2D( mesh, ice%Hs_a       , Hs_b       )
    CALL ddx_a_to_a_2D( mesh, ice%Hs_a       , dHs_dx_a   )
    CALL ddy_a_to_a_2D( mesh, ice%Hs_a       , dHs_dy_a   )
    CALL ddx_a_to_b_2D( mesh, ice%Hs_a       , dHs_dx_b   )
    CALL ddy_a_to_b_2D( mesh, ice%Hs_a       , dHs_dy_b   )
    CALL map_a_to_b_3D( mesh, ice%A_flow_3D_a, A_flow_3D_b)

    ! Calculate velocities and strain rates according to the analytical solution of the SIA:
    ! (see also Bueler and Brown, 2009, Eqs. 12-13)
    !
    !   D( z) = -2 (rho g)^n (abs(grad H))^(n-1) int_b_z( A(T*) (h - zeta)^n ) dzeta
    !   u( z) = dh/dx D( z)
    !   v( z) = dh/dy D( z)
    !
    !   du/dz( z) = -2 (rho g)^n (abs(grad H))^(n-1) A(T*) (h - z)^n dh/dx
    !   dv/dz( z) = -2 (rho g)^n (abs(grad H))^(n-1) A(T*) (h - z)^n dh/dy

    ! Calculate velocities
    DO ti = mesh%ti1, mesh%ti2

      ! Calculate the integral from b to z of (A_flow * (h - zeta)^n) dzeta
      z = Hs_b( ti) - C%zeta * Hi_b( ti)
      CALL vertical_integration_from_bottom_to_zeta( z, A_flow_3D_b( ti,:) * (Hs_b( ti) - z)**C%n_flow, int_A_hminzetan)

      ! Calculate the diffusivity term
      abs_grad_Hs = SQRT( dHs_dx_b( ti)**2 + dHs_dy_b( ti)**2)
      SIA%D_3D_b( ti,:) = -2._dp * (ice_density * grav)**C%n_flow * abs_grad_Hs**(C%n_flow - 1._dp) * int_A_hminzetan

      ! Safety
      SIA%D_3D_b( ti,:) = MAX( D_SIA_limit, SIA%D_3D_b( ti,:))

      ! Calculate the velocities
      SIA%u_3D_b( ti,:) = SIA%D_3D_b( ti,:) * dHs_dx_b( ti)
      SIA%v_3D_b( ti,:) = SIA%D_3D_b( ti,:) * dHs_dy_b( ti)

    END DO

    ! Calculate strain rates
    DO vi = mesh%vi1, mesh%vi2

      abs_grad_Hs = SQRT( dHs_dx_a( vi)**2 + dHs_dy_a( vi)**2)
      z = ice%Hs_a( vi) - C%zeta * ice%Hi_a( vi)

      DO k = 1, C%nz
        SIA%du_dz_3D_a( vi,k) = -2._dp * (ice_density * grav)**C%n_flow * abs_grad_Hs**(C%n_flow - 1._dp) * &
          ice%A_flow_3D_a( vi,k) * (ice%Hs_a( vi) - z( k))**C%n_flow * dHs_dx_a( vi)
        SIA%dv_dz_3D_a( vi,k) = -2._dp * (ice_density * grav)**C%n_flow * abs_grad_Hs**(C%n_flow - 1._dp) * &
          ice%A_flow_3D_a( vi,k) * (ice%Hs_a( vi) - z( k))**C%n_flow * dHs_dy_a( vi)
      END DO

    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wHi_b       )
    CALL deallocate_shared( wHs_b       )
    CALL deallocate_shared( wdHs_dx_a   )
    CALL deallocate_shared( wdHs_dy_a   )
    CALL deallocate_shared( wdHs_dx_b   )
    CALL deallocate_shared( wdHs_dy_b   )
    CALL deallocate_shared( wA_flow_3D_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_SIA

  SUBROUTINE remap_SIA_solver( mesh_old, mesh_new, ice, SIA)
    ! Remap the SIA solver

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_velocity_solver_SIA),      INTENT(INOUT) :: SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_SIA_solver'
    REAL(dp)                                           :: dp_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings
    dp_dummy = mesh_old%V( 1,1)
    dp_dummy = mesh_new%V( 1,1)
    dp_dummy = ice%Hi_a( 1)

    ! Reallocate shared memory

    ! Solution
    CALL reallocate_shared_dp_2D( mesh_new%nTri, C%nz, SIA%u_3D_b    , SIA%wu_3D_b    )
    CALL reallocate_shared_dp_2D( mesh_new%nTri, C%nz, SIA%v_3D_b    , SIA%wv_3D_b    )
    CALL reallocate_shared_dp_2D( mesh_new%nV  , C%nz, SIA%du_dz_3D_a, SIA%wdu_dz_3D_a)
    CALL reallocate_shared_dp_2D( mesh_new%nV  , C%nz, SIA%dv_dz_3D_a, SIA%wdv_dz_3D_a)

    ! Intermediate data fields
    CALL reallocate_shared_dp_2D( mesh_new%nTri, C%nz, SIA%D_3D_b    , SIA%wD_3D_b    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_SIA_solver

END MODULE ice_velocity_SIA_module
