MODULE mesh_operators_module

  ! Routines for calculating and applying mapping and gradient operators on the mesh.

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE petscksp
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr, petscmat_checksum
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
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
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             checksum_dp_1D, checksum_dp_2D, checksum_dp_3D

  ! Import specific functionality
  USE parallel_module,                 ONLY: adapt_shared_int_1D,    adapt_shared_dp_1D, &
                                             adapt_shared_int_2D,    adapt_shared_dp_2D, &
                                             adapt_shared_int_3D,    adapt_shared_dp_3D, &
                                             adapt_shared_bool_1D
  USE data_types_module,               ONLY: type_mesh, type_grid, type_sparse_matrix_CSR_dp, type_ice_model
  USE utilities_module,                ONLY: calc_matrix_inverse_2_by_2, calc_matrix_inverse_3_by_3, calc_matrix_inverse_general, &
                                             extend_group_single_iteration_a, extend_group_single_iteration_b, extend_group_single_iteration_c
  USE petsc_module,                    ONLY: multiply_PETSc_matrix_with_vector_1D, multiply_PETSc_matrix_with_vector_2D, mat_CSR2petsc, &
                                             mat_petsc2CSR, mat_A_eq_B_p_diagC_t_D, mat_A_eq_diag_B_t_A
  USE sparse_matrix_module,            ONLY: allocate_matrix_CSR_dist, add_entry_CSR_dist, finalise_matrix_CSR_dist, &
                                             deallocate_matrix_CSR

  IMPLICIT NONE

CONTAINS

! == Mapping

  ! 2-D
  SUBROUTINE map_a_to_b_2D( mesh, d_a, d_b)
    ! Map a 2-D data field from the a (vertex) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_a_to_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_a_b, d_a, d_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_a_to_b_2D

  SUBROUTINE map_a_to_c_2D( mesh, d_a, d_c)
    ! Map a 2-D data field from the a (vertex) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_a_to_c_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_a_c, d_a, d_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_a_to_c_2D

  SUBROUTINE map_b_to_a_2D( mesh, d_b, d_a)
    ! Map a 2-D data field from the b (triangle) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_b_to_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_b_a, d_b, d_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_b_to_a_2D

  SUBROUTINE map_b_to_c_2D( mesh, d_b, d_c)
    ! Map a 2-D data field from the b (triangle) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_b_to_c_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_b_c, d_b, d_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_b_to_c_2D

  SUBROUTINE map_c_to_a_2D( mesh, d_c, d_a)
    ! Map a 2-D data field from the c (edge) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_c_to_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_c_a, d_c, d_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_c_to_a_2D

  SUBROUTINE map_c_to_b_2D( mesh, d_c, d_b)
    ! Map a 2-D data field from the c (edge) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_c_to_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_c_b, d_c, d_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_c_to_b_2D

  ! 3-D
  SUBROUTINE map_a_to_b_3D( mesh, d_a, d_b)
    ! Map a 3-D data field from the a (vertex) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_a_to_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_a,2) /= SIZE( d_b,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_a_b, d_a, d_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_a_to_b_3D

  SUBROUTINE map_a_to_c_3D( mesh, d_a, d_c)
    ! Map a 3-D data field from the a (vertex) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_a_to_c_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_a,2) /= SIZE( d_c,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_a_c, d_a, d_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_a_to_c_3D

  SUBROUTINE map_b_to_a_3D( mesh, d_b, d_a)
    ! Map a 3-D data field from the b (triangle) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_b_to_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_b,2) /= SIZE( d_a,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_b_a, d_b, d_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_b_to_a_3D

  SUBROUTINE map_b_to_c_3D( mesh, d_b, d_c)
    ! Map a 3-D data field from the b (triangle) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_b_to_c_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_b,2) /= SIZE( d_c,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_b_c, d_b, d_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_b_to_c_3D

  SUBROUTINE map_c_to_a_3D( mesh, d_c, d_a)
    ! Map a 3-D data field from the c (edge) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_c_to_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_c,2) /= SIZE( d_a,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_c_a, d_c, d_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_c_to_a_3D

  SUBROUTINE map_c_to_b_3D( mesh, d_c, d_b)
    ! Map a 3-D data field from the c (edge) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_c_to_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_c,2) /= SIZE( d_b,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_c_b, d_c, d_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_c_to_b_3D

! == d/dx

  ! 2-D
  SUBROUTINE ddx_a_to_a_2D( mesh, d_a, ddx_a)
    ! ddx a 2-D data field from the a (vertex) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_a_to_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_a_a, d_a, ddx_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_a_to_a_2D

  SUBROUTINE ddx_a_to_b_2D( mesh, d_a, ddx_b)
    ! ddx a 2-D data field from the a (vertex) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_a_to_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_a_b, d_a, ddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_a_to_b_2D

  SUBROUTINE ddx_a_to_c_2D( mesh, d_a, ddx_c)
    ! ddx a 2-D data field from the a (vertex) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_a_to_c_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_a_c, d_a, ddx_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_a_to_c_2D

  SUBROUTINE ddx_b_to_a_2D( mesh, d_b, ddx_a)
    ! ddx a 2-D data field from the b (triangle) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_b_to_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_b_a, d_b, ddx_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_b_to_a_2D

  SUBROUTINE ddx_b_to_b_2D( mesh, d_b, ddx_b)
    ! ddx a 2-D data field from the b (triangle) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_b_to_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_b_b, d_b, ddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_b_to_b_2D

  SUBROUTINE ddx_b_to_c_2D( mesh, d_b, ddx_c)
    ! ddx a 2-D data field from the b (triangle) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_b_to_c_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_b_c, d_b, ddx_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_b_to_c_2D

  SUBROUTINE ddx_c_to_a_2D( mesh, d_c, ddx_a)
    ! ddx a 2-D data field from the c (edge) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_c_to_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_c_a, d_c, ddx_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_c_to_a_2D

  SUBROUTINE ddx_c_to_b_2D( mesh, d_c, ddx_b)
    ! ddx a 2-D data field from the c (edge) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_c_to_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_c_b, d_c, ddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_c_to_b_2D

  SUBROUTINE ddx_c_to_c_2D( mesh, d_c, ddx_c)
    ! ddx a 2-D data field from the c (edge) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_c_to_c_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_c_c, d_c, ddx_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_c_to_c_2D

  ! 3-D
  SUBROUTINE ddx_a_to_a_3D( mesh, d_a, ddx_a)
    ! ddx a 3-D data field from the a (vertex) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_a_to_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_a_a, d_a, ddx_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_a_to_a_3D

  SUBROUTINE ddx_a_to_b_3D( mesh, d_a, ddx_b)
    ! ddx a 3-D data field from the a (vertex) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_a_to_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_a_b, d_a, ddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_a_to_b_3D

  SUBROUTINE ddx_a_to_c_3D( mesh, d_a, ddx_c)
    ! ddx a 3-D data field from the a (vertex) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_a_to_c_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_a_c, d_a, ddx_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_a_to_c_3D

  SUBROUTINE ddx_b_to_a_3D( mesh, d_b, ddx_a)
    ! ddx a 3-D data field from the b (triangle) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_b_to_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_b_a, d_b, ddx_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_b_to_a_3D

  SUBROUTINE ddx_b_to_b_3D( mesh, d_b, ddx_b)
    ! ddx a 3-D data field from the b (triangle) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_b_to_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_b_b, d_b, ddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_b_to_b_3D

  SUBROUTINE ddx_b_to_c_3D( mesh, d_b, ddx_c)
    ! ddx a 3-D data field from the b (triangle) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_b_to_c_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_b_c, d_b, ddx_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_b_to_c_3D

  SUBROUTINE ddx_c_to_a_3D( mesh, d_c, ddx_a)
    ! ddx a 3-D data field from the c (edge) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_c_to_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_c_a, d_c, ddx_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_c_to_a_3D

  SUBROUTINE ddx_c_to_b_3D( mesh, d_c, ddx_b)
    ! ddx a 3-D data field from the c (edge) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_c_to_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_c_b, d_c, ddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_c_to_b_3D

  SUBROUTINE ddx_c_to_c_3D( mesh, d_c, ddx_c)
    ! ddx a 3-D data field from the c (edge) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddx_c_to_c_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_c_c, d_c, ddx_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddx_c_to_c_3D

! == d/dy

  ! 2-D
  SUBROUTINE ddy_a_to_a_2D( mesh, d_a, ddy_a)
    ! ddy a 2-D data field from the a (vertex) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_a_to_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_a_a, d_a, ddy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_a_to_a_2D

  SUBROUTINE ddy_a_to_b_2D( mesh, d_a, ddy_b)
    ! ddy a 2-D data field from the a (vertex) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_a_to_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_a_b, d_a, ddy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_a_to_b_2D

  SUBROUTINE ddy_a_to_c_2D( mesh, d_a, ddy_c)
    ! ddy a 2-D data field from the a (vertex) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_a_to_c_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_a_c, d_a, ddy_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_a_to_c_2D

  SUBROUTINE ddy_b_to_a_2D( mesh, d_b, ddy_a)
    ! ddy a 2-D data field from the b (triangle) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_b_to_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_b_a, d_b, ddy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_b_to_a_2D

  SUBROUTINE ddy_b_to_b_2D( mesh, d_b, ddy_b)
    ! ddy a 2-D data field from the b (triangle) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_b_to_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_b_b, d_b, ddy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_b_to_b_2D

  SUBROUTINE ddy_b_to_c_2D( mesh, d_b, ddy_c)
    ! ddy a 2-D data field from the b (triangle) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_b_to_c_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_b_c, d_b, ddy_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_b_to_c_2D

  SUBROUTINE ddy_c_to_a_2D( mesh, d_c, ddy_a)
    ! ddy a 2-D data field from the c (edge) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_c_to_a_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_c_a, d_c, ddy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_c_to_a_2D

  SUBROUTINE ddy_c_to_b_2D( mesh, d_c, ddy_b)
    ! ddy a 2-D data field from the c (edge) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_c_to_b_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_c_b, d_c, ddy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_c_to_b_2D

  SUBROUTINE ddy_c_to_c_2D( mesh, d_c, ddy_c)
    ! ddy a 2-D data field from the c (edge) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_c_to_c_2D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_c_c, d_c, ddy_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_c_to_c_2D

  ! 3-D
  SUBROUTINE ddy_a_to_a_3D( mesh, d_a, ddy_a)
    ! ddy a 3-D data field from the a (vertex) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_a_to_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_a_a, d_a, ddy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_a_to_a_3D

  SUBROUTINE ddy_a_to_b_3D( mesh, d_a, ddy_b)
    ! ddy a 3-D data field from the a (vertex) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_a_to_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_a_b, d_a, ddy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_a_to_b_3D

  SUBROUTINE ddy_a_to_c_3D( mesh, d_a, ddy_c)
    ! ddy a 3-D data field from the a (vertex) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_a_to_c_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_a_c, d_a, ddy_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_a_to_c_3D

  SUBROUTINE ddy_b_to_a_3D( mesh, d_b, ddy_a)
    ! ddy a 3-D data field from the b (triangle) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_b_to_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_b_a, d_b, ddy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_b_to_a_3D

  SUBROUTINE ddy_b_to_b_3D( mesh, d_b, ddy_b)
    ! ddy a 3-D data field from the b (triangle) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_b_to_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_b_b, d_b, ddy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_b_to_b_3D

  SUBROUTINE ddy_b_to_c_3D( mesh, d_b, ddy_c)
    ! ddy a 3-D data field from the b (triangle) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_b_to_c_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_b_c, d_b, ddy_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_b_to_c_3D

  SUBROUTINE ddy_c_to_a_3D( mesh, d_c, ddy_a)
    ! ddy a 3-D data field from the c (edge) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_c_to_a_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_c_a, d_c, ddy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_c_to_a_3D

  SUBROUTINE ddy_c_to_b_3D( mesh, d_c, ddy_b)
    ! ddy a 3-D data field from the c (edge) to the b (triangle) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_c_to_b_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_c_b, d_c, ddy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_c_to_b_3D

  SUBROUTINE ddy_c_to_c_3D( mesh, d_c, ddy_c)
    ! ddy a 3-D data field from the c (edge) to the c (edge) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ddy_c_to_c_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_c_c, d_c, ddy_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ddy_c_to_c_3D

! == d2/dx2, d2/dxdy, d2/dy2 on the a-grid

  SUBROUTINE d2dx2_a_to_a_2D(  mesh, d_a, d2dx2_a)
    ! d2/dx2 a 2-D data field from the a (vertex) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dx2_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'd2dx2_a_to_a_2D'
    REAL(dp), DIMENSION(:    ), POINTER                ::  ddx_b
    INTEGER                                            :: wddx_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dx2_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_b, wddx_b)

    ! ddx_b = M_ddx_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_a_b, d_a, ddx_b)

    ! d2dx2_a = M_ddx_b_a * ddx_b
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_b_a, ddx_b, d2dx2_a)

    ! Clean up after yourself
    CALL deallocate_shared( wddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE d2dx2_a_to_a_2D

  SUBROUTINE d2dxdy_a_to_a_2D( mesh, d_a, d2dxdy_a)
    ! d2/dxdy a 2-D data field from the a (vertex) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dxdy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'd2dxdy_a_to_a_2D'
    REAL(dp), DIMENSION(:    ), POINTER                ::  ddx_b
    INTEGER                                            :: wddx_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dxdy_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_b, wddx_b)

    ! ddx_b = M_ddx_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_a_b, d_a, ddx_b)

    ! d2dxdy_a = M_ddy_b_a * ddx_b
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_b_a, ddx_b, d2dxdy_a)

    ! Clean up after yourself
    CALL deallocate_shared( wddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE d2dxdy_a_to_a_2D

  SUBROUTINE d2dy2_a_to_a_2D(  mesh, d_a, d2dy2_a)
    ! d2/dy2 a 2-D data field from the a (vertex) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dy2_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'd2dy2_a_to_a_2D'
    REAL(dp), DIMENSION(:    ), POINTER                ::  ddy_b
    INTEGER                                            :: wddy_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dy2_a,1) /= mesh%nV) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_b, wddy_b)

    ! ddy_b = M_ddy_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_a_b, d_a, ddy_b)

    ! d2dy2_a = M_ddy_b_a * ddy_b
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_b_a, ddy_b, d2dy2_a)

    ! Clean up after yourself
    CALL deallocate_shared( wddy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE d2dy2_a_to_a_2D

  SUBROUTINE d2dx2_a_to_a_3D(  mesh, d_a, d2dx2_a)
    ! d2/dx2 a 3-D data field from the a (vertex) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dx2_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'd2dx2_a_to_a_3D'
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  ddx_b
    INTEGER                                            :: wddx_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dx2_a,1) /= mesh%nV .OR. SIZE( d_a,2) /= SIZE( d2dx2_a,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nTri, SIZE( d_a,2), ddx_b, wddx_b)

    ! ddx_b = M_ddx_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_a_b, d_a, ddx_b)

    ! d2dx2_a = M_ddx_b_a * ddx_b
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_b_a, ddx_b, d2dx2_a)

    ! Clean up after yourself
    CALL deallocate_shared( wddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE d2dx2_a_to_a_3D

  SUBROUTINE d2dxdy_a_to_a_3D( mesh, d_a, d2dxdy_a)
    ! d2/dxdy a 3-D data field from the a (vertex) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dxdy_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'd2dxdy_a_to_a_3D'
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  ddx_b
    INTEGER                                            :: wddx_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dxdy_a,1) /= mesh%nV .OR. SIZE( d_a,2) /= SIZE( d2dxdy_a,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nTri, SIZE( d_a,2), ddx_b, wddx_b)

    ! ddx_b = M_ddx_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_a_b, d_a, ddx_b)

    ! d2dxdy_a = M_ddy_b_a * ddx_b
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_b_a, ddx_b, d2dxdy_a)

    ! Clean up after yourself
    CALL deallocate_shared( wddx_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE d2dxdy_a_to_a_3D

  SUBROUTINE d2dy2_a_to_a_3D(  mesh, d_a, d2dy2_a)
    ! d2/dy2 a 3-D data field from the a (vertex) to the a (vertex) grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dy2_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'd2dy2_a_to_a_3D'
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  ddy_b
    INTEGER                                            :: wddy_b

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dy2_a,1) /= mesh%nV .OR. SIZE( d_a,2) /= SIZE( d2dy2_a,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nTri, SIZE( d_a,2), ddy_b, wddy_b)

    ! ddy_b = M_ddy_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_a_b, d_a, ddy_b)

    ! d2dy2_a = M_ddy_b_a * ddy_b
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_b_a, ddy_b, d2dy2_a)

    ! Clean up after yourself
    CALL deallocate_shared( wddy_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE d2dy2_a_to_a_3D

! == Calculate the matrix operators for mapping and gradients

  SUBROUTINE calc_matrix_operators_mesh_basic( mesh)
    ! Calculate all the geometry-independent matrix operators (i.e. in [x',y',zeta]-coordinates only)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_basic'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate grid-cell-to-matrix-row translation tables
    CALL calc_grid_cell_to_matrix_row_translation_tables( mesh)

    ! Calculate mapping and gradient operators between all three grids (vertices, triangles, and edges)
    CALL calc_matrix_operators_mesh_a_a( mesh)
    CALL calc_matrix_operators_mesh_a_b( mesh)
    CALL calc_matrix_operators_mesh_a_c( mesh)

    CALL calc_matrix_operators_mesh_b_a( mesh)
    CALL calc_matrix_operators_mesh_b_b( mesh)
    CALL calc_matrix_operators_mesh_b_c( mesh)

    CALL calc_matrix_operators_mesh_c_a( mesh)
    CALL calc_matrix_operators_mesh_c_b( mesh)
    CALL calc_matrix_operators_mesh_c_c( mesh)

    ! Calculate 2nd-order accurate operators on the b-grid (triangles), used for solving the SSA/DIVA/BPA
    CALL calc_matrix_operators_mesh_b_b_2nd_order( mesh)

    ! Calculate matrices representing mapping operators between the scalar and vector b grids
    CALL calc_buv_matrices( mesh)

    ! Calculate all the basic 3-D matrix operators on the mesh
    CALL calc_matrix_operators_mesh_basic_3D( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_basic

  SUBROUTINE calc_matrix_operators_mesh_basic_3D( mesh)
    ! Calculate all the geometry-independent 3-D matrix operators on the mesh

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_basic_3D'

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Extend horizontal operators (map, d/dx, d/dy) in the vertical

    ! ak-ak
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2vi , mesh%vi2n , mesh%n2vi , mesh%vi2n , &
      mesh%n2vik, mesh%vik2n, mesh%n2vik, mesh%vik2n, mesh%M_ddx_a_a, mesh%M_ddxp_ak_ak)
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2vi , mesh%vi2n , mesh%n2vi , mesh%vi2n , &
      mesh%n2vik, mesh%vik2n, mesh%n2vik, mesh%vik2n, mesh%M_ddy_a_a, mesh%M_ddyp_ak_ak)

    ! bk-bk
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2ti , mesh%ti2n , mesh%n2ti , mesh%ti2n , &
      mesh%n2tik, mesh%tik2n, mesh%n2tik, mesh%tik2n, mesh%M_ddx_b_b, mesh%M_ddxp_bk_bk)
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2ti , mesh%ti2n , mesh%n2ti , mesh%ti2n , &
      mesh%n2tik, mesh%tik2n, mesh%n2tik, mesh%tik2n, mesh%M_ddy_b_b, mesh%M_ddyp_bk_bk)

    ! 2nd order
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2ti , mesh%ti2n , mesh%n2ti , mesh%ti2n , &
      mesh%n2tik, mesh%tik2n, mesh%n2tik, mesh%tik2n, mesh%M2_ddx_b_b, mesh%M2_ddxp_bk_bk)
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2ti , mesh%ti2n , mesh%n2ti , mesh%ti2n , &
      mesh%n2tik, mesh%tik2n, mesh%n2tik, mesh%tik2n, mesh%M2_ddy_b_b, mesh%M2_ddyp_bk_bk)
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2ti , mesh%ti2n , mesh%n2ti , mesh%ti2n , &
      mesh%n2tik, mesh%tik2n, mesh%n2tik, mesh%tik2n, mesh%M2_d2dx2_b_b, mesh%M2_d2dxp2_bk_bk)
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2ti , mesh%ti2n , mesh%n2ti , mesh%ti2n , &
      mesh%n2tik, mesh%tik2n, mesh%n2tik, mesh%tik2n, mesh%M2_d2dxdy_b_b, mesh%M2_d2dxpdyp_bk_bk)
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2ti , mesh%ti2n , mesh%n2ti , mesh%ti2n , &
      mesh%n2tik, mesh%tik2n, mesh%n2tik, mesh%tik2n, mesh%M2_d2dy2_b_b, mesh%M2_d2dyp2_bk_bk)

    ! ak-bk
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2ti , mesh%ti2n , mesh%n2vi , mesh%vi2n , &
      mesh%n2tik, mesh%tik2n, mesh%n2vik, mesh%vik2n, mesh%M_map_a_b, mesh%M_map_ak_bk)
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2ti , mesh%ti2n , mesh%n2vi , mesh%vi2n , &
      mesh%n2tik, mesh%tik2n, mesh%n2vik, mesh%vik2n, mesh%M_ddx_a_b, mesh%M_ddxp_ak_bk)
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2ti , mesh%ti2n , mesh%n2vi , mesh%vi2n , &
      mesh%n2tik, mesh%tik2n, mesh%n2vik, mesh%vik2n, mesh%M_ddy_a_b, mesh%M_ddyp_ak_bk)

    ! bk-ak
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2vi , mesh%vi2n , mesh%n2ti , mesh%ti2n , &
      mesh%n2vik, mesh%vik2n, mesh%n2tik, mesh%tik2n, mesh%M_map_b_a, mesh%M_map_bk_ak)
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2vi , mesh%vi2n , mesh%n2ti , mesh%ti2n , &
      mesh%n2vik, mesh%vik2n, mesh%n2tik, mesh%tik2n, mesh%M_ddx_b_a, mesh%M_ddxp_bk_ak)
    CALL extend_matrix_operator_2D_to_3D( &
      mesh%n2vi , mesh%vi2n , mesh%n2ti , mesh%ti2n , &
      mesh%n2vik, mesh%vik2n, mesh%n2tik, mesh%tik2n, mesh%M_ddy_b_a, mesh%M_ddyp_bk_ak)

  ! == Calculate mapping and d/dzeta operators between the regular and staggered vertical grids

    ! Operators in the vertical column (1-D, no sense in recalculating the exact same coefficients for every grid point)
    CALL calc_vertical_operators_reg_1D( C%zeta, mesh%M_ddzeta_k_k_1D, mesh%M_d2dzeta2_k_k_1D, mesh%M_fkp1_m_fk_1D, mesh%M_fk_m_fkm1_1D)
    CALL calc_vertical_operators_stag_1D( C%zeta, C%zeta_stag, mesh%M_map_k_ks_1D, mesh%M_ddzeta_k_ks_1D, mesh%M_map_ks_k_1D, mesh%M_ddzeta_ks_k_1D)

    ! Extrude them horizontally to get them on the 3-D grid(s)

    ! ak-ak
    CALL convert_vertical_operator_from_1D_to_3D_k_k( mesh%M_ddzeta_k_k_1D  , mesh%n2vik, mesh%vik2n,                           mesh%M_ddzeta_ak_ak   )
    CALL convert_vertical_operator_from_1D_to_3D_k_k( mesh%M_d2dzeta2_k_k_1D, mesh%n2vik, mesh%vik2n,                           mesh%M_d2dzeta2_ak_ak )
    ! bk-bk
    CALL convert_vertical_operator_from_1D_to_3D_k_k( mesh%M_ddzeta_k_k_1D  , mesh%n2tik, mesh%tik2n,                           mesh%M_ddzeta_bk_bk   )
    CALL convert_vertical_operator_from_1D_to_3D_k_k( mesh%M_d2dzeta2_k_k_1D, mesh%n2tik, mesh%tik2n,                           mesh%M_d2dzeta2_bk_bk )
    CALL convert_vertical_operator_from_1D_to_3D_k_k( mesh%M_fkp1_m_fk_1D   , mesh%n2tik, mesh%tik2n,                           mesh%M_fkp1_m_fk_bk_bk)
    CALL convert_vertical_operator_from_1D_to_3D_k_k( mesh%M_fk_m_fkm1_1D   , mesh%n2tik, mesh%tik2n,                           mesh%M_fk_m_fkm1_bk_bk)
    ! ak-aks
    CALL convert_vertical_operator_from_1D_to_3D_k_ks( mesh%M_map_k_ks_1D   , mesh%n2vik, mesh%vik2n, mesh%n2viks, mesh%viks2n, mesh%M_map_ak_aks     )
    CALL convert_vertical_operator_from_1D_to_3D_k_ks( mesh%M_ddzeta_k_ks_1D, mesh%n2vik, mesh%vik2n, mesh%n2viks, mesh%viks2n, mesh%M_ddzeta_ak_aks  )
    CALL convert_vertical_operator_from_1D_to_3D_ks_k( mesh%M_map_ks_k_1D   , mesh%n2vik, mesh%vik2n, mesh%n2viks, mesh%viks2n, mesh%M_map_aks_ak     )
    CALL convert_vertical_operator_from_1D_to_3D_ks_k( mesh%M_ddzeta_ks_k_1D, mesh%n2vik, mesh%vik2n, mesh%n2viks, mesh%viks2n, mesh%M_ddzeta_aks_ak  )
    ! bk-bks
    CALL convert_vertical_operator_from_1D_to_3D_k_ks( mesh%M_map_k_ks_1D   , mesh%n2tik, mesh%tik2n, mesh%n2tiks, mesh%tiks2n, mesh%M_map_bk_bks     )
    CALL convert_vertical_operator_from_1D_to_3D_k_ks( mesh%M_ddzeta_k_ks_1D, mesh%n2tik, mesh%tik2n, mesh%n2tiks, mesh%tiks2n, mesh%M_ddzeta_bk_bks  )
    CALL convert_vertical_operator_from_1D_to_3D_ks_k( mesh%M_map_ks_k_1D   , mesh%n2tik, mesh%tik2n, mesh%n2tiks, mesh%tiks2n, mesh%M_map_bks_bk     )
    CALL convert_vertical_operator_from_1D_to_3D_ks_k( mesh%M_ddzeta_ks_k_1D, mesh%n2tik, mesh%tik2n, mesh%n2tiks, mesh%tiks2n, mesh%M_ddzeta_bks_bk  )

    ! Zeta operators in tridiagonal form for efficient use in thermodynamics
    CALL calc_zeta_operators_tridiagonal( mesh)

    ! 3-D cross-terms
    CALL MatMatMult( mesh%M_ddzeta_bk_bk, mesh%M2_ddxp_bk_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, mesh%M2_d2dxpdzeta_bk_bk, perr)
    CALL MatMatMult( mesh%M_ddzeta_bk_bk, mesh%M2_ddyp_bk_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, mesh%M2_d2dypdzeta_bk_bk, perr)

  ! == Calculate some useful mapping operators between the grids used in the BPA solver

    CALL MatMatMult( mesh%M_map_ak_bk , mesh%M_ddzeta_ak_ak, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, mesh%M_ddzeta_ak_bk, perr)
    CALL MatMatMult( mesh%M_map_bk_ak , mesh%M_ddzeta_bk_bk, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, mesh%M_ddzeta_bk_ak, perr)

    CALL MatMatMult( mesh%M_map_bk_ak , mesh%M_map_bks_bk  , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, mesh%M_map_bks_ak, perr)
    CALL MatMatMult( mesh%M_map_bk_bks, mesh%M_map_ak_bk   , MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, mesh%M_map_ak_bks, perr)

  ! == Calculate mapping operators between 2-D grid and 3-D grid surface/base layers

    CALL calc_maps_b_to_bk_surf_base( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_basic_3D

  SUBROUTINE calc_matrix_operators_x_y_z_3D( mesh, ice)
    ! Calculate all the geometry-dependent matrix operators on the mesh (i.e. in [x,y,z]-coordinates)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_x_y_z_3D'
    TYPE(tMat)                                         :: M1
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_ak,  d_bk,  d_bks
    INTEGER                                            :: wd_ak, wd_bk, wd_bks
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dzeta_dx_sq,  dzeta_dy_sq,  dzeta_dx_dzeta_dy
    INTEGER                                            :: wdzeta_dx_sq, wdzeta_dy_sq, wdzeta_dx_dzeta_dy
    INTEGER                                            :: ti,k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nnak , d_ak , wd_ak )
    CALL allocate_shared_dp_1D( mesh%nnbk , d_bk , wd_bk )
    CALL allocate_shared_dp_1D( mesh%nnbks, d_bks, wd_bks)

    ! Calculate zeta gradients
    CALL calc_zeta_gradients( mesh, ice)

  ! == ak-ak

    ! d/dx = d/dxp + dzeta/dx d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddx_ak_ak, perr)
    ! Convert dzeta/dx from field form to vector form
    CALL field2vec_ak( mesh, ice%dzeta_dx_ak, d_ak)
    ! d/dx = d/dxp + dzeta/dx d/dzeta
    CALL mat_A_eq_B_p_diagC_t_D( mesh%M_ddxp_ak_ak, d_ak, mesh%M_ddzeta_ak_ak, mesh%M_ddx_ak_ak)

    ! d/dy = d/dyp + dzeta/dy d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddy_ak_ak, perr)
    ! Convert dzeta/dy from field form to vector form
    CALL field2vec_ak( mesh, ice%dzeta_dy_ak, d_ak)
    ! d/dy = d/dyp + dzeta/dy d/dzeta
    CALL mat_A_eq_B_p_diagC_t_D( mesh%M_ddyp_ak_ak, d_ak, mesh%M_ddzeta_ak_ak, mesh%M_ddy_ak_ak)

    ! d/dz = dzeta/dz * d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddz_ak_ak, perr)
    ! d/dz = d/dzeta
    CALL MatDuplicate( mesh%M_ddzeta_ak_ak, MAT_COPY_VALUES, mesh%M_ddz_ak_ak, perr)
    ! Convert dzeta/dz from field form to vector form
    CALL field2vec_ak( mesh, ice%dzeta_dz_ak, d_ak)
    ! d/dz = dzeta/dz * d/dzeta
    CALL mat_A_eq_diag_B_t_A( mesh%M_ddz_ak_ak, d_ak)

  ! == bk-bk

    ! d/dx = d/dxp + dzeta/dx d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddx_bk_bk, perr)
    ! Convert dzeta/dx from field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dx_bk, d_bk)
    ! d/dx = d/dxp + dzeta/dx d/dzeta
    CALL mat_A_eq_B_p_diagC_t_D( mesh%M_ddxp_bk_bk, d_bk, mesh%M_ddzeta_bk_bk, mesh%M_ddx_bk_bk)

    ! d/dy = d/dyp + dzeta/dy d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddy_bk_bk, perr)
    ! Convert dzeta/dy from field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dy_bk, d_bk)
    ! d/dy = d/dyp + dzeta/dy d/dzeta
    CALL mat_A_eq_B_p_diagC_t_D( mesh%M_ddyp_bk_bk, d_bk, mesh%M_ddzeta_bk_bk, mesh%M_ddy_bk_bk)

    ! d/dz = dzeta/dz * d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddz_bk_bk, perr)
    ! d/dz = d/dzeta
    CALL MatDuplicate( mesh%M_ddzeta_bk_bk, MAT_COPY_VALUES, mesh%M_ddz_bk_bk, perr)
    ! Convert dzeta/dz from field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dz_bk, d_bk)
    ! d/dz = dzeta/dz * d/dzeta
    CALL mat_A_eq_diag_B_t_A( mesh%M_ddz_bk_bk, d_bk)

  ! == ak-bk

    ! d/dx = d/dxp + dzeta/dx d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddx_ak_bk, perr)
    ! Convert dzeta/dx from field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dx_bk, d_bk)
    ! d/dx = d/dxp + dzeta/dx d/dzeta
    CALL mat_A_eq_B_p_diagC_t_D( mesh%M_ddxp_ak_bk, d_bk, mesh%M_ddzeta_ak_bk, mesh%M_ddx_ak_bk)

    ! d/dy = d/dyp + dzeta/dy d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddy_ak_bk, perr)
    ! Convert dzeta/dy from field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dy_bk, d_bk)
    ! d/dy = d/dyp + dzeta/dy d/dzeta
    CALL mat_A_eq_B_p_diagC_t_D( mesh%M_ddyp_ak_bk, d_bk, mesh%M_ddzeta_ak_bk, mesh%M_ddy_ak_bk)

    ! d/dz = dzeta/dz * d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddz_ak_bk, perr)
    ! d/dz = d/dzeta
    CALL MatDuplicate( mesh%M_ddzeta_ak_bk, MAT_COPY_VALUES, mesh%M_ddz_ak_bk, perr)
    ! Convert dzeta/dz from field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dz_bk, d_bk)
    ! d/dz = dzeta/dz * d/dzeta
    CALL mat_A_eq_diag_B_t_A( mesh%M_ddz_ak_bk, d_bk)

  ! == bk-ak

    ! d/dx = d/dxp + dzeta/dx d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddx_bk_ak, perr)
    ! Convert dzeta/dx from field form to vector form
    CALL field2vec_ak( mesh, ice%dzeta_dx_ak, d_ak)
    ! d/dx = d/dxp + dzeta/dx d/dzeta
    CALL mat_A_eq_B_p_diagC_t_D( mesh%M_ddxp_bk_ak, d_ak, mesh%M_ddzeta_bk_ak, mesh%M_ddx_bk_ak)

    ! d/dy = d/dyp + dzeta/dy d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddy_bk_ak, perr)
    ! Convert dzeta/dy from field form to vector form
    CALL field2vec_ak( mesh, ice%dzeta_dy_ak, d_ak)
    ! d/dy = d/dyp + dzeta/dy d/dzeta
    CALL mat_A_eq_B_p_diagC_t_D( mesh%M_ddyp_bk_ak, d_ak, mesh%M_ddzeta_bk_ak, mesh%M_ddy_bk_ak)

    ! d/dz = dzeta/dz * d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddz_bk_ak, perr)
    ! d/dz = d/dzeta
    CALL MatDuplicate( mesh%M_ddzeta_bk_ak, MAT_COPY_VALUES, mesh%M_ddz_bk_ak, perr)
    ! Convert dzeta/dz from field form to vector form
    CALL field2vec_ak( mesh, ice%dzeta_dz_ak, d_ak)
    ! d/dz = dzeta/dz * d/dzeta
    CALL mat_A_eq_diag_B_t_A( mesh%M_ddz_bk_ak, d_ak)

  ! == bk-bks

    ! d/dz = dzeta/dz d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddz_bk_bks, perr)
    ! d/dz = d/dzeta
    CALL MatDuplicate( mesh%M_ddzeta_bk_bks, MAT_COPY_VALUES, mesh%M_ddz_bk_bks, perr)
    ! Convert dzeta/dz from field form to vector form
    CALL field2vec_bks( mesh, ice%dzeta_dz_bks, d_bks)
    ! d/dz = dzeta/dz d/dzeta
    CALL mat_A_eq_diag_B_t_A( mesh%M_ddz_bk_bks, d_bks)

    ! d/dz = dzeta/dz d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M_ddz_bks_bk, perr)
    ! d/dz = d/dzeta
    CALL MatDuplicate( mesh%M_ddzeta_bks_bk, MAT_COPY_VALUES, mesh%M_ddz_bks_bk, perr)
    ! Convert dzeta/dz from field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dz_bk, d_bk)
    ! d/dz = dzeta/dz d/dzeta
    CALL mat_A_eq_diag_B_t_A( mesh%M_ddz_bks_bk, d_bk)

  ! == bk-bk, 2nd order

    ! Calculate (dzeta/dx)^2, (dzeta/dy)^2
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, dzeta_dx_sq      , wdzeta_dx_sq      )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, dzeta_dy_sq      , wdzeta_dy_sq      )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, dzeta_dx_dzeta_dy, wdzeta_dx_dzeta_dy)

    DO ti = mesh%ti1, mesh%ti2
      DO k = 1, C%nz
        dzeta_dx_sq(       ti,k) = ice%dzeta_dx_bk( ti,k)**2
        dzeta_dy_sq(       ti,k) = ice%dzeta_dy_bk( ti,k)**2
        dzeta_dx_dzeta_dy( ti,k) = ice%dzeta_dx_bk( ti,k) * ice%dzeta_dy_bk( ti,k)
      END DO
    END DO

    ! d/dx = d/dxp + dzeta/dx d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M2_ddx_bk_bk, perr)
    ! Convert dzeta/dx from field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dx_bk, d_bk)
    ! d/dx = d/dxp + dzeta/dx d/dzeta
    CALL mat_A_eq_B_p_diagC_t_D( mesh%M2_ddxp_bk_bk, d_bk, mesh%M_ddzeta_bk_bk, mesh%M2_ddx_bk_bk)

    ! d/dy = d/dyp + dzeta/dy d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M2_ddy_bk_bk, perr)
    ! Convert dzeta/dy from field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dy_bk, d_bk)
    ! d/dy = d/dyp + dzeta/dy d/dzeta
    CALL mat_A_eq_B_p_diagC_t_D( mesh%M2_ddyp_bk_bk, d_bk, mesh%M_ddzeta_bk_bk, mesh%M2_ddy_bk_bk)

    ! d2/dx2 = d2/dxp2 + (dzeta/dx)^2 d2/dzeta2  + 2 dzeta/dx d2/dxpdzeta + d2zeta/dx2 d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M2_d2dx2_bk_bk, perr)
    ! d2/dx2 = d2/dxp2 ...
    CALL MatDuplicate( mesh%M2_d2dxp2_bk_bk , MAT_COPY_VALUES, mesh%M2_d2dx2_bk_bk, perr)
    ! M1 = d2/dzeta2
    CALL MatDuplicate( mesh%M_d2dzeta2_bk_bk, MAT_COPY_VALUES, M1                 , perr)
    ! Convert (dzeta/dx)^2 from field form to vector form
    CALL field2vec_bk( mesh, dzeta_dx_sq, d_bk)
    ! M1 = (dzeta/dx)^2 d2/dzeta2
    CALL mat_A_eq_diag_B_t_A( M1, d_bk)
    ! d2/dx2 = d2/dxp2 + (dzeta/dx)^2 d2/dzeta2 ...
    CALL MatAXPY( mesh%M2_d2dx2_bk_bk, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    ! Clean up after yourself
    CALL MatDestroy( M1, perr)
    ! M1 = d2/dxpdzeta
    CALL MatDuplicate( mesh%M2_d2dxpdzeta_bk_bk, MAT_COPY_VALUES, M1                 , perr)
    ! Convert dzeta/dx from field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dx_bk, d_bk)
    ! M1 = dzeta/dx d2/dxpdzeta
    CALL mat_A_eq_diag_B_t_A( M1, d_bk)
    ! d2/dx2 = d2/dxp2 + (dzeta/dx)^2 d2/dzeta2  + 2 dzeta/dx d2/dxpdzeta ...
    CALL MatAXPY( mesh%M2_d2dx2_bk_bk, 2._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    ! Clean up after yourself
    CALL MatDestroy( M1, perr)
    ! M1 = d/dzeta
    CALL MatDuplicate( mesh%M_ddzeta_bk_bk, MAT_COPY_VALUES, M1                 , perr)
    ! Convert d2zeta/dx2 form field form to vector form
    CALL field2vec_bk( mesh, ice%d2zeta_dx2_bk, d_bk)
    ! M1 = d2zeta/dx2 d/dzeta
    CALL mat_A_eq_diag_B_t_A( M1, d_bk)
    ! d2/dx2 = d2/dxp2 + (dzeta/dx)^2 d2/dzeta2  + 2 dzeta/dx d2/dxpdzeta + d2zeta/dx2 d/dzeta
    CALL MatAXPY( mesh%M2_d2dx2_bk_bk, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    ! Clean up after yourself
    CALL MatDestroy( M1, perr)

    ! d2/dy2 = d2/dyp2 + (dzeta/dy)^2 d2/dzeta2  + 2 dzeta/dy d2/dypdzeta + d2zeta/dy2 d/dzeta

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M2_d2dy2_bk_bk, perr)
    ! d2/dy2 = d2/dyp2 ...
    CALL MatDuplicate( mesh%M2_d2dyp2_bk_bk , MAT_COPY_VALUES, mesh%M2_d2dy2_bk_bk, perr)
    ! M1 = d2/dzeta2
    CALL MatDuplicate( mesh%M_d2dzeta2_bk_bk, MAT_COPY_VALUES, M1                 , perr)
    ! Convert (dzeta/dy)^2 from field form to vector form
    CALL field2vec_bk( mesh, dzeta_dy_sq, d_bk)
    ! M1 = (dzeta/dy)^2 d2/dzeta2
    CALL mat_A_eq_diag_B_t_A( M1, d_bk)
    ! d2/dy2 = d2/dyp2 + (dzeta/dy)^2 d2/dzeta2 ...
    CALL MatAXPY( mesh%M2_d2dy2_bk_bk, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    ! Clean up after yourself
    CALL MatDestroy( M1, perr)
    ! M1 = d2/dypdzeta
    CALL MatDuplicate( mesh%M2_d2dypdzeta_bk_bk, MAT_COPY_VALUES, M1                 , perr)
    ! Convert dzeta/dy form field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dy_bk, d_bk)
    ! M1 = dzeta/dy d2/dypdzeta
    CALL mat_A_eq_diag_B_t_A( M1, d_bk)
    ! d2/dy2 = d2/dyp2 + (dzeta/dy)^2 d2/dzeta2  + 2 dzeta/dy d2/dypdzeta ...
    CALL MatAXPY( mesh%M2_d2dy2_bk_bk, 2._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    ! Clean up after yourself
    CALL MatDestroy( M1, perr)
    ! M1 = d/dzeta
    CALL MatDuplicate( mesh%M_ddzeta_bk_bk, MAT_COPY_VALUES, M1                 , perr)
    ! Convert d2zeta/dy2 form field form to vector form
    CALL field2vec_bk( mesh, ice%d2zeta_dy2_bk, d_bk)
    ! M1 = d2zeta/dy2 d/dzeta
    CALL mat_A_eq_diag_B_t_A( M1, d_bk)
    ! d2/dy2 = d2/dyp2 + (dzeta/dy)^2 d2/dzeta2  + 2 dzeta/dy d2/dypdzeta + d2zeta/dy2 d/dzeta
    CALL MatAXPY( mesh%M2_d2dy2_bk_bk, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    ! Clean up after yourself
    CALL MatDestroy( M1, perr)

    ! d2/dxdy = d2/dxpdyp + d2zeta/dxdy d/dzeta + dzeta/dx d2/dypdzeta + dzeta/dy d2/dxpdzeta  + dzeta/dx dzeta/dy d2/dzeta2

    ! Clean up extant version if needed
    CALL MatDestroy( mesh%M2_d2dxdy_bk_bk, perr)
    ! d2/dxdy = d2/dxpdyp ...
    CALL MatDuplicate( mesh%M2_d2dxpdyp_bk_bk , MAT_COPY_VALUES, mesh%M2_d2dxdy_bk_bk, perr)
    ! M1 = d/dzeta
    CALL MatDuplicate( mesh%M_ddzeta_bk_bk, MAT_COPY_VALUES, M1                 , perr)
    ! Convert d2zeta/dxdy form field form to vector form
    CALL field2vec_bk( mesh, ice%d2zeta_dxdy_bk, d_bk)
    ! M1 = d2zeta/dxdy d/dzeta
    CALL mat_A_eq_diag_B_t_A( M1, d_bk)
    ! d2/dxdy = d2/dxpdyp + d2zeta/dxdy d/dzeta ...
    CALL MatAXPY( mesh%M2_d2dxdy_bk_bk, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    ! Clean up after yourself
    CALL MatDestroy( M1, perr)
    ! M1 = d2/dypdzeta
    CALL MatDuplicate( mesh%M2_d2dypdzeta_bk_bk, MAT_COPY_VALUES, M1                 , perr)
    ! Convert dzeta/dx form field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dx_bk, d_bk)
    ! M1 = dzeta/dx d2/dypdzeta
    CALL mat_A_eq_diag_B_t_A( M1, d_bk)
    ! d2/dxdy = d2/dxpdyp + d2zeta/dxdy d/dzeta + dzeta/dx d2/dypdzeta ...
    CALL MatAXPY( mesh%M2_d2dxdy_bk_bk, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    ! Clean up after yourself
    CALL MatDestroy( M1, perr)
    ! M1 = d2/dxpdzeta
    CALL MatDuplicate( mesh%M2_d2dxpdzeta_bk_bk, MAT_COPY_VALUES, M1                 , perr)
    ! Convert dzeta/dy form field form to vector form
    CALL field2vec_bk( mesh, ice%dzeta_dy_bk, d_bk)
    ! M1 = dzeta/dy d2/dxpdzeta
    CALL mat_A_eq_diag_B_t_A( M1, d_bk)
    ! d2/dxdy = d2/dxpdyp + d2zeta/dxdy d/dzeta + dzeta/dx d2/dypdzeta + dzeta/dy d2/dxpdzeta ...
    CALL MatAXPY( mesh%M2_d2dxdy_bk_bk, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    ! Clean up after yourself
    CALL MatDestroy( M1, perr)
    ! M1 = d2/dzeta2
    CALL MatDuplicate( mesh%M_d2dzeta2_bk_bk, MAT_COPY_VALUES, M1                 , perr)
    ! Convert dzeta/dx dzeta/dy form field form to vector form
    CALL field2vec_bk( mesh, dzeta_dx_dzeta_dy, d_bk)
    ! M1 = dzeta/dx dzeta/dy d2/dzeta2
    CALL mat_A_eq_diag_B_t_A( M1, d_bk)
    ! d2/dxdy = d2/dxpdyp + d2zeta/dxdy d/dzeta + dzeta/dx d2/dypdzeta + dzeta/dy d2/dxpdzeta  + dzeta/dx dzeta/dy d2/dzeta2
    CALL MatAXPY( mesh%M2_d2dxdy_bk_bk, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    ! Clean up after yourself
    CALL MatDestroy( M1, perr)

    ! Clean up after yourself
    CALL deallocate_shared( wd_ak             )
    CALL deallocate_shared( wd_bk             )
    CALL deallocate_shared( wd_bks            )
    CALL deallocate_shared( wdzeta_dx_sq      )
    CALL deallocate_shared( wdzeta_dy_sq      )
    CALL deallocate_shared( wdzeta_dx_dzeta_dy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_x_y_z_3D

  SUBROUTINE calc_zeta_gradients( mesh, ice)
    ! Calculate all the gradients of zeta, needed to perform the scaled vertical coordinate transformation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_zeta_gradients'
!    REAL(dp), DIMENSION(:    ), POINTER                :: Hi_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dx_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dx2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dxdy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dy2_a
    INTEGER                                            :: wHi_a, wdHi_dx_a, wdHi_dy_a, wd2Hi_dx2_a, wd2Hi_dxdy_a, wd2Hi_dy2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: Hi_b
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dx_b
    REAL(dp), DIMENSION(:    ), POINTER                :: dHi_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dx2_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dxdy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hi_dy2_b
    REAL(dp), DIMENSION(:    ), POINTER                :: Hs_b
    INTEGER                                            :: wHi_b, wdHi_dx_b, wdHi_dy_b, wd2Hi_dx2_b, wd2Hi_dxdy_b, wd2Hi_dy2_b
!    REAL(dp), DIMENSION(:    ), POINTER                :: Hs_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dx_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dx2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dxdy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dy2_a
    INTEGER                                            :: wHs_a, wdHs_dx_a, wdHs_dy_a, wd2Hs_dx2_a, wd2Hs_dxdy_a, wd2Hs_dy2_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dx_b
    REAL(dp), DIMENSION(:    ), POINTER                :: dHs_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dx2_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dxdy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2Hs_dy2_b
    INTEGER                                            :: wHs_b, wdHs_dx_b, wdHs_dy_b, wd2Hs_dx2_b, wd2Hs_dxdy_b, wd2Hs_dy2_b
    INTEGER                                            :: vi,ti,k,ks
    REAL(dp)                                           :: Hi, zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
!    CALL allocate_shared_dp_1D( mesh%nV  , Hi_a       , wHi_a       )
    CALL allocate_shared_dp_1D( mesh%nV  , dHi_dx_a   , wdHi_dx_a   )
    CALL allocate_shared_dp_1D( mesh%nV  , dHi_dy_a   , wdHi_dy_a   )
    CALL allocate_shared_dp_1D( mesh%nV  , d2Hi_dx2_a , wd2Hi_dx2_a )
    CALL allocate_shared_dp_1D( mesh%nV  , d2Hi_dxdy_a, wd2Hi_dxdy_a)
    CALL allocate_shared_dp_1D( mesh%nV  , d2Hi_dy2_a , wd2Hi_dy2_a )

    CALL allocate_shared_dp_1D( mesh%nTri, Hi_b       , wHi_b       )
    CALL allocate_shared_dp_1D( mesh%nTri, dHi_dx_b   , wdHi_dx_b   )
    CALL allocate_shared_dp_1D( mesh%nTri, dHi_dy_b   , wdHi_dy_b   )
    CALL allocate_shared_dp_1D( mesh%nTri, d2Hi_dx2_b , wd2Hi_dx2_b )
    CALL allocate_shared_dp_1D( mesh%nTri, d2Hi_dxdy_b, wd2Hi_dxdy_b)
    CALL allocate_shared_dp_1D( mesh%nTri, d2Hi_dy2_b , wd2Hi_dy2_b )

!    CALL allocate_shared_dp_1D( mesh%nV  , Hs_a       , wHs_a       )
    CALL allocate_shared_dp_1D( mesh%nV  , dHs_dx_a   , wdHs_dx_a   )
    CALL allocate_shared_dp_1D( mesh%nV  , dHs_dy_a   , wdHs_dy_a   )
    CALL allocate_shared_dp_1D( mesh%nV  , d2Hs_dx2_a , wd2Hs_dx2_a )
    CALL allocate_shared_dp_1D( mesh%nV  , d2Hs_dxdy_a, wd2Hs_dxdy_a)
    CALL allocate_shared_dp_1D( mesh%nV  , d2Hs_dy2_a , wd2Hs_dy2_a )

    CALL allocate_shared_dp_1D( mesh%nTri, Hs_b       , wHs_b       )
    CALL allocate_shared_dp_1D( mesh%nTri, dHs_dx_b   , wdHs_dx_b   )
    CALL allocate_shared_dp_1D( mesh%nTri, dHs_dy_b   , wdHs_dy_b   )
    CALL allocate_shared_dp_1D( mesh%nTri, d2Hs_dx2_b , wd2Hs_dx2_b )
    CALL allocate_shared_dp_1D( mesh%nTri, d2Hs_dxdy_b, wd2Hs_dxdy_b)
    CALL allocate_shared_dp_1D( mesh%nTri, d2Hs_dy2_b , wd2Hs_dy2_b )

  ! Calculate gradients of Hi and Hs on both grids

    !CALL map_a_to_a_2D( mesh, ice%Hi_a, Hi_a       )
    CALL ddx_a_to_a_2D( mesh, ice%Hi_a, dHi_dx_a   )
    CALL ddy_a_to_a_2D( mesh, ice%Hi_a, dHi_dy_a   )
    CALL map_a_to_b_2D( mesh, ice%Hi_a, Hi_b       )
    CALL ddx_a_to_b_2D( mesh, ice%Hi_a, dHi_dx_b   )
    CALL ddy_a_to_b_2D( mesh, ice%Hi_a, dHi_dy_b   )
    CALL ddx_b_to_a_2D( mesh, dHi_dx_b, d2Hi_dx2_a )
    CALL ddy_b_to_a_2D( mesh, dHi_dx_b, d2Hi_dxdy_a)
    CALL ddy_b_to_a_2D( mesh, dHi_dy_b, d2Hi_dy2_a )
    CALL ddx_a_to_b_2D( mesh, dHi_dx_a, d2Hi_dx2_b )
    CALL ddy_a_to_b_2D( mesh, dHi_dx_a, d2Hi_dxdy_b)
    CALL ddy_a_to_b_2D( mesh, dHi_dy_a, d2Hi_dy2_b )

    !CALL map_a_to_a_2D( mesh, ice%Hs_a, Hs_a       )
    CALL ddx_a_to_a_2D( mesh, ice%Hs_a, dHs_dx_a   )
    CALL ddy_a_to_a_2D( mesh, ice%Hs_a, dHs_dy_a   )
    CALL map_a_to_b_2D( mesh, ice%Hs_a, Hs_b       )
    CALL ddx_a_to_b_2D( mesh, ice%Hs_a, dHs_dx_b   )
    CALL ddy_a_to_b_2D( mesh, ice%Hs_a, dHs_dy_b   )
    CALL ddx_b_to_a_2D( mesh, dHs_dx_b, d2Hs_dx2_a )
    CALL ddy_b_to_a_2D( mesh, dHs_dx_b, d2Hs_dxdy_a)
    CALL ddy_b_to_a_2D( mesh, dHs_dy_b, d2Hs_dy2_a )
    CALL ddx_a_to_b_2D( mesh, dHs_dx_a, d2Hs_dx2_b )
    CALL ddy_a_to_b_2D( mesh, dHs_dx_a, d2Hs_dxdy_b)
    CALL ddy_a_to_b_2D( mesh, dHs_dy_a, d2Hs_dy2_b )

    ! Calculate zeta gradients on all grids

    ! ak
    DO vi = mesh%vi1, mesh%vi2

      Hi = MAX( 0.1_dp, ice%Hi_a( vi))

      DO k = 1, C%nz

        zeta = C%zeta( k)

        ice%dzeta_dt_ak(    vi,k) = ( 1._dp / Hi) * (ice%dHs_dt_a( vi) - zeta * ice%dHi_dt_a( vi))

        ice%dzeta_dx_ak(    vi,k) = ( 1._dp / Hi) * (dHs_dx_a( vi) - zeta * dHi_dx_a( vi))
        ice%dzeta_dy_ak(    vi,k) = ( 1._dp / Hi) * (dHs_dy_a( vi) - zeta * dHi_dy_a( vi))
        ice%dzeta_dz_ak(    vi,k) = (-1._dp / Hi)

        ice%d2zeta_dx2_ak(  vi,k) = (dHi_dx_a( vi) * -1._dp / Hi) * ice%dzeta_dx_ak( vi,k) + (1._dp / Hi) * (d2Hs_dx2_a(  vi) - zeta * d2Hi_dx2_a(  vi))
        ice%d2zeta_dxdy_ak( vi,k) = (dHi_dy_a( vi) * -1._dp / Hi) * ice%dzeta_dx_ak( vi,k) + (1._dp / Hi) * (d2Hs_dxdy_a( vi) - zeta * d2Hi_dxdy_a( vi))
        ice%d2zeta_dy2_ak(  vi,k) = (dHi_dy_a( vi) * -1._dp / Hi) * ice%dzeta_dy_ak( vi,k) + (1._dp / Hi) * (d2Hs_dy2_a(  vi) - zeta * d2Hi_dy2_a(  vi))

      END DO
    END DO

    ! bk
    DO ti = mesh%ti1, mesh%ti2

      Hi = MAX( 0.1_dp, Hi_b( ti))

      DO k = 1, C%nz

        zeta = C%zeta( k)

        ice%dzeta_dx_bk(    ti,k) = ( 1._dp / Hi) * (dHs_dx_b( ti) - zeta * dHi_dx_b( ti))
        ice%dzeta_dy_bk(    ti,k) = ( 1._dp / Hi) * (dHs_dy_b( ti) - zeta * dHi_dy_b( ti))
        ice%dzeta_dz_bk(    ti,k) = (-1._dp / Hi)

        ice%d2zeta_dx2_bk(  ti,k) = (dHi_dx_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bk( ti,k) + (1._dp / Hi) * (d2Hs_dx2_b(  ti) - zeta * d2Hi_dx2_b(  ti))
        ice%d2zeta_dxdy_bk( ti,k) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bk( ti,k) + (1._dp / Hi) * (d2Hs_dxdy_b( ti) - zeta * d2Hi_dxdy_b( ti))
        ice%d2zeta_dy2_bk(  ti,k) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dy_bk( ti,k) + (1._dp / Hi) * (d2Hs_dy2_b(  ti) - zeta * d2Hi_dy2_b(  ti))

      END DO
    END DO

    ! bks
    DO ti = mesh%ti1, mesh%ti2

      Hi = MAX( 0.1_dp, Hi_b( ti))

      DO ks = 1, C%nz-1

        zeta = C%zeta_stag( ks)

        ice%dzeta_dx_bks(    ti,ks) = ( 1._dp / Hi) * (dHs_dx_b( ti) - zeta * dHi_dx_b( ti))
        ice%dzeta_dy_bks(    ti,ks) = ( 1._dp / Hi) * (dHs_dy_b( ti) - zeta * dHi_dy_b( ti))
        ice%dzeta_dz_bks(    ti,ks) = (-1._dp / Hi)

        ice%d2zeta_dx2_bks(  ti,ks) = (dHi_dx_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dx2_b(  ti) - zeta * d2Hi_dx2_b(  ti))
        ice%d2zeta_dxdy_bks( ti,ks) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dx_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dxdy_b( ti) - zeta * d2Hi_dxdy_b( ti))
        ice%d2zeta_dy2_bks(  ti,ks) = (dHi_dy_b( ti) * -1._dp / Hi) * ice%dzeta_dy_bks( ti,ks) + (1._dp / Hi) * (d2Hs_dy2_b(  ti) - zeta * d2Hi_dy2_b(  ti))

      END DO
    END DO

    ! Clean up after yourself
!    CALL deallocate_shared( wHi_a       )
    CALL deallocate_shared( wdHi_dx_a   )
    CALL deallocate_shared( wdHi_dy_a   )
    CALL deallocate_shared( wd2Hi_dx2_a )
    CALL deallocate_shared( wd2Hi_dxdy_a)
    CALL deallocate_shared( wd2Hi_dy2_a )
    CALL deallocate_shared( wHi_b       )
    CALL deallocate_shared( wdHi_dx_b   )
    CALL deallocate_shared( wdHi_dy_b   )
    CALL deallocate_shared( wd2Hi_dx2_b )
    CALL deallocate_shared( wd2Hi_dxdy_b)
    CALL deallocate_shared( wd2Hi_dy2_b )
!    CALL deallocate_shared( wHs_a       )
    CALL deallocate_shared( wdHs_dx_a   )
    CALL deallocate_shared( wdHs_dy_a   )
    CALL deallocate_shared( wd2Hs_dx2_a )
    CALL deallocate_shared( wd2Hs_dxdy_a)
    CALL deallocate_shared( wd2Hs_dy2_a )
    CALL deallocate_shared( wHs_b       )
    CALL deallocate_shared( wdHs_dx_b   )
    CALL deallocate_shared( wdHs_dy_b   )
    CALL deallocate_shared( wd2Hs_dx2_b )
    CALL deallocate_shared( wd2Hs_dxdy_b)
    CALL deallocate_shared( wd2Hs_dy2_b )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_zeta_gradients

! == Calculate mapping and gradient operators between the a-, b-, and c-grids

  SUBROUTINE calc_matrix_operators_mesh_a_a( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the a-grid (vertices) and the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_a_a'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_ddx, M_ddy
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: vi
    REAL(dp)                                           :: x, y
    INTEGER                                            :: vj
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp)                                           :: Nfx_i, Nfy_i
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 2

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nV      ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2

      ! The vertex represented by this matrix row
      vi = mesh%n2vi( row,1)
      x  = mesh%V( vi,1)
      y  = mesh%V( vi,2)

      ! Initialise the list of neighbours: just vi itself
      mesh%Vmap        = 0
      mesh%Vmap( vi)   = 2
      mesh%Vstack1     = 0
      mesh%VstackN1    = 1
      mesh%Vstack1( 1) = vi

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (mesh%VstackN1 - 1 < n_neighbours_min)
        CALL extend_group_single_iteration_a( mesh)
      END DO

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, mesh%VstackN1
        IF (n_c == n_neighbours_max) EXIT
        vj = mesh%Vstack1( i)
        IF (vj == vi) CYCLE
        n_c = n_c + 1
        i_c( n_c) = vj
        x_c( n_c) = mesh%V( vj,1)
        y_c( n_c) = mesh%V( vj,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_reg_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      CALL add_entry_CSR_dist( M_ddx, row, row, Nfx_i)
      CALL add_entry_CSR_dist( M_ddy, row, row, Nfy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        vj = i_c( i)
        col = mesh%vi2n( vj)
        CALL add_entry_CSR_dist( M_ddx, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( M_ddy, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_ddx, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddy, row1, row2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_ddx, mesh%M_ddx_a_a)
    CALL mat_CSR2petsc( M_ddy, mesh%M_ddy_a_a)

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    CALL deallocate_matrix_CSR( M_ddx)
    CALL deallocate_matrix_CSR( M_ddy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_a_a

  SUBROUTINE calc_matrix_operators_mesh_a_b( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the a-grid (vertices) and the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_a_b'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map, M_ddx, M_ddy
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, vi
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nV        ! from
    nrows           = mesh%nTri      ! to
    nnz_per_row_est = 3
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2

      ! The vertex represented by this matrix row
      ti = mesh%n2ti( row,1)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Initialise the list of neighbours: the three vertices spanning ti
      mesh%Vmap     = 0
      mesh%Vstack1  = 0
      mesh%VstackN1 = 0
      DO n = 1, 3
        vi = mesh%Tri( ti,n)
        mesh%Vmap( vi) = 2
        mesh%VstackN1 = mesh%VstackN1 + 1
        mesh%Vstack1( mesh%VstackN1) = vi
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (mesh%VstackN1 < n_neighbours_min)
        CALL extend_group_single_iteration_a( mesh)
      END DO

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, mesh%VstackN1
        IF (n_c == n_neighbours_max) EXIT
        vi = mesh%Vstack1( i)
        n_c = n_c + 1
        i_c( n_c) = vi
        x_c( n_c) = mesh%V( vi,1)
        y_c( n_c) = mesh%V( vi,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        vi = i_c( i)
        col = mesh%vi2n( vi)
        CALL add_entry_CSR_dist( M_map, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( M_ddx, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( M_ddy, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_map, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddx, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddy, row1, row2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map, mesh%M_map_a_b)
    CALL mat_CSR2petsc( M_ddx, mesh%M_ddx_a_b)
    CALL mat_CSR2petsc( M_ddy, mesh%M_ddy_a_b)

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    CALL deallocate_matrix_CSR( M_map)
    CALL deallocate_matrix_CSR( M_ddx)
    CALL deallocate_matrix_CSR( M_ddy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_a_b

  SUBROUTINE calc_matrix_operators_mesh_a_c( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the a-grid (vertices) and the c-grid (edges)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_a_c'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map, M_ddx, M_ddy
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: aci
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, vi
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nV        ! from
    nrows           = mesh%nAc       ! to
    nnz_per_row_est = 4
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2

      ! The edge represented by this matrix row
      aci = mesh%n2ci( row,1)
      x   = mesh%VAc( aci,1)
      y   = mesh%VAc( aci,2)

      ! Initialise the list of neighbours: the three or four vertices spanning aci
      mesh%Vmap     = 0
      mesh%Vstack1  = 0
      mesh%VstackN1 = 0
      DO n = 1, 4
        vi = mesh%Aci( aci,n)
        IF (vi == 0) CYCLE
        mesh%Vmap( vi) = 2
        mesh%VstackN1 = mesh%VstackN1 + 1
        mesh%Vstack1( mesh%VstackN1) = vi
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (mesh%VstackN1 < n_neighbours_min)
        CALL extend_group_single_iteration_a( mesh)
      END DO

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, mesh%VstackN1
        IF (n_c == n_neighbours_max) EXIT
        vi = mesh%Vstack1( i)
        n_c = n_c + 1
        i_c( n_c) = vi
        x_c( n_c) = mesh%V( vi,1)
        y_c( n_c) = mesh%V( vi,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        vi = i_c( i)
        col = mesh%vi2n( vi)
        CALL add_entry_CSR_dist( M_map, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( M_ddx, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( M_ddy, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_map, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddx, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddy, row1, row2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map, mesh%M_map_a_c)
    CALL mat_CSR2petsc( M_ddx, mesh%M_ddx_a_c)
    CALL mat_CSR2petsc( M_ddy, mesh%M_ddy_a_c)

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    CALL deallocate_matrix_CSR( M_map)
    CALL deallocate_matrix_CSR( M_ddx)
    CALL deallocate_matrix_CSR( M_ddy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_a_c

  SUBROUTINE calc_matrix_operators_mesh_b_a( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the b-grid (triangles) and the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_b_a'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map, M_ddx, M_ddy
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: vi
    REAL(dp)                                           :: x, y
    INTEGER                                            :: iti, ti
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nTri    ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2

      ! The vertex represented by this matrix row
      vi = mesh%n2vi( row,1)
      x  = mesh%V( vi,1)
      y  = mesh%V( vi,2)

      ! Initialise the list of neighbours: all the triangles surrounding vi
      mesh%Trimap     = 0
      mesh%Tristack1  = 0
      mesh%TristackN1 = 0
      DO iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        mesh%Trimap( ti) = 2
        mesh%TristackN1 = mesh%TristackN1 + 1
        mesh%Tristack1( mesh%TristackN1) = ti
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (mesh%TristackN1 < n_neighbours_min)
        CALL extend_group_single_iteration_b( mesh)
      END DO

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, mesh%TristackN1
        IF (n_c == n_neighbours_max) EXIT
        ti = mesh%Tristack1( i)
        n_c = n_c + 1
        i_c( n_c) = ti
        x_c( n_c) = mesh%TriGC( ti,1)
        y_c( n_c) = mesh%TriGC( ti,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        ti = i_c( i)
        col = mesh%ti2n( ti)
        CALL add_entry_CSR_dist( M_map, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( M_ddx, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( M_ddy, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_map, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddx, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddy, row1, row2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map, mesh%M_map_b_a)
    CALL mat_CSR2petsc( M_ddx, mesh%M_ddx_b_a)
    CALL mat_CSR2petsc( M_ddy, mesh%M_ddy_b_a)

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    CALL deallocate_matrix_CSR( M_map)
    CALL deallocate_matrix_CSR( M_ddx)
    CALL deallocate_matrix_CSR( M_ddy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_b_a

  SUBROUTINE calc_matrix_operators_mesh_b_b( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the b-grid (triangles) and the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_b_b'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_ddx, M_ddy
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, tj
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp)                                           :: Nfx_i, Nfy_i
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 2

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nTri      ! from
    nrows           = mesh%nTri      ! to
    nnz_per_row_est = 3
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2

      ! The vertex represented by this matrix row
      ti = mesh%n2ti( row,1)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Initialise the list of neighbours: just ti itself
      mesh%Trimap     = 0
      mesh%Tristack1  = 0
      mesh%TristackN1 = 0
      mesh%Trimap( ti) = 2
      mesh%TristackN1 = 1
      mesh%Tristack1( 1) = ti

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (mesh%TristackN1 - 1 < n_neighbours_min)
        CALL extend_group_single_iteration_b( mesh)
      END DO

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, mesh%TristackN1
        IF (n_c == n_neighbours_max) EXIT
        tj = mesh%Tristack1( i)
        IF (tj == ti) CYCLE
        n_c = n_c + 1
        i_c( n_c) = tj
        x_c( n_c) = mesh%TriGC( tj,1)
        y_c( n_c) = mesh%TriGC( tj,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_reg_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      CALL add_entry_CSR_dist( M_ddx, row, row, Nfx_i)
      CALL add_entry_CSR_dist( M_ddy, row, row, Nfy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        tj = i_c( i)
        col = mesh%ti2n( tj)
        CALL add_entry_CSR_dist( M_ddx, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( M_ddy, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_ddx, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddy, row1, row2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_ddx, mesh%M_ddx_b_b)
    CALL mat_CSR2petsc( M_ddy, mesh%M_ddy_b_b)

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    CALL deallocate_matrix_CSR( M_ddx)
    CALL deallocate_matrix_CSR( M_ddy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_b_b

  SUBROUTINE calc_matrix_operators_mesh_b_b_2nd_order( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the b-grid (triangles) and the b-grid (triangles)
    !
    ! Use 2nd-order accurate shape functions to calculate d/dx, d/dy, d2/dx2, d2/dxdy, and d2/dy2.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_b_b_2nd_order'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_ddx, M_ddy, M_d2dx2, M_d2dxdy, M_d2dy2
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, tj
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp)                                           :: Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 5

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nTri      ! from
    nrows           = mesh%nTri      ! to
    nnz_per_row_est = 3
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_ddx   , nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy   , nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_d2dx2 , nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_d2dxdy, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_d2dy2 , nrows, ncols, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))
    ALLOCATE( Nfxx_c( n_neighbours_max))
    ALLOCATE( Nfxy_c( n_neighbours_max))
    ALLOCATE( Nfyy_c( n_neighbours_max))

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2

      ! The vertex represented by this matrix row
      ti = mesh%n2ti( row,1)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Initialise the list of neighbours: just ti itself
      mesh%Trimap     = 0
      mesh%Tristack1  = 0
      mesh%TristackN1 = 0
      mesh%Trimap( ti) = 2
      mesh%TristackN1 = 1
      mesh%Tristack1( 1) = ti

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (mesh%TristackN1 - 1 < n_neighbours_min)
        CALL extend_group_single_iteration_b( mesh)
      END DO

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, mesh%TristackN1
        IF (n_c == n_neighbours_max) EXIT
        tj = mesh%Tristack1( i)
        IF (tj == ti) CYCLE
        n_c = n_c + 1
        i_c( n_c) = tj
        x_c( n_c) = mesh%TriGC( tj,1)
        y_c( n_c) = mesh%TriGC( tj,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_reg_2nd_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i, Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c)

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      CALL add_entry_CSR_dist( M_ddx   , row, row, Nfx_i )
      CALL add_entry_CSR_dist( M_ddy   , row, row, Nfy_i )
      CALL add_entry_CSR_dist( M_d2dx2 , row, row, Nfxx_i)
      CALL add_entry_CSR_dist( M_d2dxdy, row, row, Nfxy_i)
      CALL add_entry_CSR_dist( M_d2dy2 , row, row, Nfyy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        tj = i_c( i)
        col = mesh%ti2n( tj)
        CALL add_entry_CSR_dist( M_ddx   , row, col, Nfx_c(  i))
        CALL add_entry_CSR_dist( M_ddy   , row, col, Nfy_c(  i))
        CALL add_entry_CSR_dist( M_d2dx2 , row, col, Nfxx_c( i))
        CALL add_entry_CSR_dist( M_d2dxdy, row, col, Nfxy_c( i))
        CALL add_entry_CSR_dist( M_d2dy2 , row, col, Nfyy_c( i))
      END DO

    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_ddx   , row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddy   , row1, row2)
    CALL finalise_matrix_CSR_dist( M_d2dx2 , row1, row2)
    CALL finalise_matrix_CSR_dist( M_d2dxdy, row1, row2)
    CALL finalise_matrix_CSR_dist( M_d2dy2 , row1, row2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_ddx   , mesh%M2_ddx_b_b   )
    CALL mat_CSR2petsc( M_ddy   , mesh%M2_ddy_b_b   )
    CALL mat_CSR2petsc( M_d2dx2 , mesh%M2_d2dx2_b_b )
    CALL mat_CSR2petsc( M_d2dxdy, mesh%M2_d2dxdy_b_b)
    CALL mat_CSR2petsc( M_d2dy2 , mesh%M2_d2dy2_b_b )

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    DEALLOCATE( Nfxx_c)
    DEALLOCATE( Nfxy_c)
    DEALLOCATE( Nfyy_c)
    CALL deallocate_matrix_CSR( M_ddx   )
    CALL deallocate_matrix_CSR( M_ddy   )
    CALL deallocate_matrix_CSR( M_d2dx2 )
    CALL deallocate_matrix_CSR( M_d2dxdy)
    CALL deallocate_matrix_CSR( M_d2dy2 )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_b_b_2nd_order

  SUBROUTINE calc_matrix_operators_mesh_b_c( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the b-grid (triangles) and the c-grid (edges)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_b_c'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map, M_ddx, M_ddy
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: aci
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, ti
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nTri      ! from
    nrows           = mesh%nAc       ! to
    nnz_per_row_est = 6
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2

      ! The edge represented by this matrix row
      aci = mesh%n2ci( row,1)
      x   = mesh%VAc( aci,1)
      y   = mesh%VAc( aci,2)

      ! Initialise the list of neighbours: the two triangles adjacent to aci
      mesh%Trimap     = 0
      mesh%Tristack1  = 0
      mesh%TristackN1 = 0
      DO n = 5, 6
        ti = mesh%Aci( aci,n)
        IF (ti == 0) CYCLE
        mesh%Trimap( ti) = 2
        mesh%TristackN1 = mesh%TristackN1 + 1
        mesh%Tristack1( mesh%TristackN1 ) = ti
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (mesh%TristackN1 < n_neighbours_min)
        CALL extend_group_single_iteration_b( mesh)
      END DO

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, mesh%TristackN1
        IF (n_c == n_neighbours_max) EXIT
        ti = mesh%Tristack1( i)
        n_c = n_c + 1
        i_c( n_c) = ti
        x_c( n_c) = mesh%TriGC( ti,1)
        y_c( n_c) = mesh%TriGC( ti,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        ti = i_c( i)
        col = mesh%ti2n( ti)
        CALL add_entry_CSR_dist( M_map, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( M_ddx, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( M_ddy, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_map, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddx, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddy, row1, row2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map, mesh%M_map_b_c)
    CALL mat_CSR2petsc( M_ddx, mesh%M_ddx_b_c)
    CALL mat_CSR2petsc( M_ddy, mesh%M_ddy_b_c)

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    CALL deallocate_matrix_CSR( M_map)
    CALL deallocate_matrix_CSR( M_ddx)
    CALL deallocate_matrix_CSR( M_ddy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_b_c

  SUBROUTINE calc_matrix_operators_mesh_c_a( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the c-grid (edges) and the a-grid (vertices)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_c_a'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map, M_ddx, M_ddy
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: vi
    REAL(dp)                                           :: x, y
    INTEGER                                            :: ci, aci
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nAc     ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2

      ! The vertex represented by this matrix row
      vi = mesh%n2vi( row,1)
      x  = mesh%V( vi,1)
      y  = mesh%V( vi,2)

      ! Initialise the list of neighbours: all the edges surrounding vi
      mesh%Emap     = 0
      mesh%Estack1  = 0
      mesh%EstackN1 = 0
      DO ci = 1, mesh%nC( vi)
        aci = mesh%iaci( vi,ci)
        mesh%Emap( aci) = 2
        mesh%EstackN1 = mesh%EstackN1 + 1
        mesh%Estack1( mesh%EstackN1) = aci
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (mesh%EstackN1 < n_neighbours_min)
        CALL extend_group_single_iteration_c( mesh)
      END DO

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, mesh%EstackN1
        IF (n_c == n_neighbours_max) EXIT
        aci = mesh%Estack1( i)
        n_c = n_c + 1
        i_c( n_c) = aci
        x_c( n_c) = mesh%VAc( aci,1)
        y_c( n_c) = mesh%VAc( aci,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        aci = i_c( i)
        col = mesh%ci2n( aci)
        CALL add_entry_CSR_dist( M_map, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( M_ddx, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( M_ddy, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_map, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddx, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddy, row1, row2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map, mesh%M_map_c_a)
    CALL mat_CSR2petsc( M_ddx, mesh%M_ddx_c_a)
    CALL mat_CSR2petsc( M_ddy, mesh%M_ddy_c_a)

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    CALL deallocate_matrix_CSR( M_map)
    CALL deallocate_matrix_CSR( M_ddx)
    CALL deallocate_matrix_CSR( M_ddy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_c_a

  SUBROUTINE calc_matrix_operators_mesh_c_b( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the c-grid (edges) and the b-grid (triangles)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_c_b'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map, M_ddx, M_ddy
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y
    INTEGER                                            :: n, n2, vi, vj, ci, aci
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nf_c, Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 3

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nAc       ! from
    nrows           = mesh%nTri      ! to
    nnz_per_row_est = 3
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nf_c(   n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2

      ! The vertex represented by this matrix row
      ti = mesh%n2ti( row,1)
      x  = mesh%TriGC( ti,1)
      y  = mesh%TriGC( ti,2)

      ! Initialise the list of neighbours: the three edges spanning ti
      mesh%Emap     = 0
      mesh%Estack1  = 0
      mesh%EstackN1 = 0
      DO n = 1, 3
        n2 = n+1
        IF (n2 == 4) n2 = 1
        vi = mesh%Tri( ti,n )
        vj = mesh%Tri( ti,n2)
        aci = 0
        DO ci = 1, mesh%nC( vi)
          IF (mesh%C( vi,ci) == vj) THEN
            aci = mesh%iAci( vi,ci)
            EXIT
          END IF
        END DO
        IF (aci == 0) CALL crash('couldnt find edge connecting vi and vj!')
        mesh%Emap( aci) = 2
        mesh%EstackN1 = mesh%EstackN1 + 1
        mesh%Estack1( mesh%EstackN1) = aci
      END DO

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (mesh%EstackN1 < n_neighbours_min)
        CALL extend_group_single_iteration_c( mesh)
      END DO

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, mesh%EstackN1
        IF (n_c == n_neighbours_max) EXIT
        aci = mesh%Estack1( i)
        n_c = n_c + 1
        i_c( n_c) = aci
        x_c( n_c) = mesh%VAc( aci,1)
        y_c( n_c) = mesh%VAc( aci,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        aci = i_c( i)
        col = mesh%ci2n( aci)
        CALL add_entry_CSR_dist( M_map, row, col, Nf_c(  i))
        CALL add_entry_CSR_dist( M_ddx, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( M_ddy, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_map, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddx, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddy, row1, row2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map, mesh%M_map_c_b)
    CALL mat_CSR2petsc( M_ddx, mesh%M_ddx_c_b)
    CALL mat_CSR2petsc( M_ddy, mesh%M_ddy_c_b)

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nf_c  )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    CALL deallocate_matrix_CSR( M_map)
    CALL deallocate_matrix_CSR( M_ddx)
    CALL deallocate_matrix_CSR( M_ddy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_c_b

  SUBROUTINE calc_matrix_operators_mesh_c_c( mesh)
    ! Calculate mapping, d/dx, and d/dy matrix operators between the c-grid (edges) and the c-grid (edges)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_mesh_c_c'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_ddx, M_ddy
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: aci
    REAL(dp)                                           :: x, y
    INTEGER                                            :: acj
    INTEGER                                            :: n_neighbours_min
    INTEGER                                            :: n_neighbours_max
    INTEGER                                            :: i
    INTEGER                                            :: n_c
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    REAL(dp)                                           :: Nfx_i, Nfy_i
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfx_c, Nfy_c
    INTEGER                                            :: col

    ! Add routine to path
    CALL init_routine( routine_name)

    n_neighbours_min = 2

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nAc     ! from
    nrows           = mesh%nAc     ! to
    nnz_per_row_est = mesh%nC_mem+1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_est_proc)

  ! Calculate shape functions and fill them into the matrices
  ! =========================================================

    ! Allocate memory for map, stack, and shape functions
    n_neighbours_max = mesh%nC_mem**2
    ALLOCATE( i_c(    n_neighbours_max))
    ALLOCATE( x_c(    n_neighbours_max))
    ALLOCATE( y_c(    n_neighbours_max))
    ALLOCATE( Nfx_c(  n_neighbours_max))
    ALLOCATE( Nfy_c(  n_neighbours_max))

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2

      ! The vertex represented by this matrix row
      aci = mesh%n2ci( row,1)
      x   = mesh%VAc( aci,1)
      y   = mesh%VAc( aci,2)

      ! Initialise the list of neighbours: just aci itself
      mesh%Emap     = 0
      mesh%Estack1  = 0
      mesh%EstackN1 = 0
      mesh%Emap( aci) = 2
      mesh%EstackN1 = 1
      mesh%Estack1( 1) = aci

      ! Extend outward until enough neighbours are found to calculate the shape functions
      DO WHILE (mesh%EstackN1 - 1 < n_neighbours_min)
        CALL extend_group_single_iteration_c( mesh)
      END DO

      ! Get the coordinates of the neighbours
      n_c = 0
      DO i = 1, mesh%EstackN1
        IF (n_c == n_neighbours_max) EXIT
        acj = mesh%Estack1( i)
        IF (acj == aci) CYCLE
        n_c = n_c + 1
        i_c( n_c) = acj
        x_c( n_c) = mesh%VAc( acj,1)
        y_c( n_c) = mesh%VAc( acj,2)
      END DO

      ! Calculate shape functions
      CALL calc_shape_functions_2D_reg_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c)

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      CALL add_entry_CSR_dist( M_ddx, row, row, Nfx_i)
      CALL add_entry_CSR_dist( M_ddy, row, row, Nfy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      DO i = 1, n_c
        acj = i_c( i)
        col = mesh%ci2n( acj)
        CALL add_entry_CSR_dist( M_ddx, row, col, Nfx_c( i))
        CALL add_entry_CSR_dist( M_ddy, row, col, Nfy_c( i))
      END DO

    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_ddx, row1, row2)
    CALL finalise_matrix_CSR_dist( M_ddy, row1, row2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_ddx, mesh%M_ddx_c_c)
    CALL mat_CSR2petsc( M_ddy, mesh%M_ddy_c_c)

    ! Clean up after yourself
    DEALLOCATE( i_c   )
    DEALLOCATE( x_c   )
    DEALLOCATE( y_c   )
    DEALLOCATE( Nfx_c )
    DEALLOCATE( Nfy_c )
    CALL deallocate_matrix_CSR( M_ddx)
    CALL deallocate_matrix_CSR( M_ddy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_mesh_c_c

! == Shape functions

  SUBROUTINE calc_shape_functions_1D_reg_2nd_order( x, n_max, n_c, x_c, Nfx_i, Nfxx_i, Nfx_c, Nfxx_c)
    ! Calculate shape functions...
    ! ...in one dimension...
    ! ...on the regular grid (i.e. f is known)...
    ! ...to 2nd-order accuracy.
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x          ! The location where we want to know the gradients
    INTEGER,                             INTENT(IN)    :: n_max      ! The maximum number of surrounding points
    INTEGER,                             INTENT(IN)    :: n_c        ! The number  of     surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c        ! Coordinates of the surrounding points where we know f
    REAL(dp),                            INTENT(OUT)   :: Nfx_i      ! d/dx   shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfxx_i     ! d2/dx2 shape function for the point [x,y]
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      ! d/dx   shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfxx_c     ! d2/dx2 shape functions for the surrounding points

    ! Local variables:
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n_c)                           :: dx, w
    REAL(dp), DIMENSION(2,2)                           :: ATWTWA, M

    ! Safety
    IF (n_c < 2) CALL crash('calc_shape_functions_1D_reg_2nd_order needs at least 2 neighbours!')

    ! Calculate distances relative to x
    DO ci = 1, n_c
      dx( ci) = x_c( ci) - x
    END DO

    ! Calculate the weights w
    DO ci = 1, n_c
      w( ci) = 1._dp / (ABS( dx( ci))**q)
    END DO

    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n_c
      ATWTWA( 1,1) = ATWTWA( 1,1) + w(ci)**2 *       dx( ci)    *       dx( ci)
      ATWTWA( 1,2) = ATWTWA( 1,2) + w(ci)**2 *       dx( ci)    * 1/2 * dx( ci)**2

      ATWTWA( 2,1) = ATWTWA( 2,1) + w(ci)**2 * 1/2 * dx( ci)**2 *       dx( ci)
      ATWTWA( 2,2) = ATWTWA( 2,2) + w(ci)**2 * 1/2 * dx( ci)**2 * 1/2 * dx( ci)**2
    END DO

    ! Invert ATWTWA to find M
    CALL calc_matrix_inverse_2_by_2( ATWTWA, M)

    ! Calculate shape functions
    Nfx_c   = 0._dp
    Nfxx_c  = 0._dp
    DO ci = 1, n_c
      Nfx_c(   ci) = w( ci)**2 * ( &
        (M( 1,1) *        dx( ci)   ) + &
        (M( 1,2) * 1/2  * dx( ci)**2))
      Nfxx_c(  ci) = w( ci)**2 * ( &
        (M( 2,1) *        dx( ci)   ) + &
        (M( 2,2) * 1/2  * dx( ci)**2))
    END DO

    Nfx_i  = -SUM( Nfx_c )
    Nfxx_i = -SUM( Nfxx_c)

  END SUBROUTINE calc_shape_functions_1D_reg_2nd_order

  SUBROUTINE calc_shape_functions_1D_stag_2nd_order( x, n_max, n_c, x_c, Nf_c, Nfx_c)
    ! Calculate shape functions...
    ! ...in one dimension...
    ! ...on the staggered grid (i.e. f is not known)...
    ! ...to 2nd-order accuracy.
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x          ! The location where we want to know the gradients
    INTEGER,                             INTENT(IN)    :: n_max      ! The maximum number of surrounding points
    INTEGER,                             INTENT(IN)    :: n_c        ! The number  of     surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c      ! Coordinates of the surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nf_c       ! map    shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      ! d/dx   shape functions for the surrounding points

    ! Local variables:
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n_c)                           :: dx, w
    REAL(dp), DIMENSION(2,2)                           :: ATWTWA, M

    ! Safety
    IF (n_c < 2) CALL crash('calc_shape_functions_1D_stag_2nd_order needs at least 2 neighbours!')

    ! Calculate distances relative to x
    DO ci = 1, n_c
      dx( ci) = x_c( ci) - x
    END DO

    ! Calculate the weights w
    DO ci = 1, n_c
      w( ci) = 1._dp / (ABS( dx( ci))**q)
    END DO

    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n_c
      ATWTWA( 1,1) = ATWTWA( 1,1) + (w( ci)**2 * 1       * 1      )
      ATWTWA( 1,2) = ATWTWA( 1,2) + (w( ci)**2 * 1       * dx( ci))

      ATWTWA( 2,1) = ATWTWA( 2,1) + (w( ci)**2 * dx( ci) * 1      )
      ATWTWA( 2,2) = ATWTWA( 2,2) + (w( ci)**2 * dx( ci) * dx( ci))
    END DO

    ! Invert ATWTWA to find M
    CALL calc_matrix_inverse_2_by_2( ATWTWA, M)

    ! Calculate shape functions
    Nf_c   = 0._dp
    Nfx_c  = 0._dp
    DO ci = 1, n_c
      Nf_c(   ci) = w( ci)**2 * ( &
        (M( 1,1) *        1         ) + &
        (M( 1,2) *        dx( ci)**2))
      Nfx_c(  ci) = w( ci)**2 * ( &
        (M( 2,1) *        1         ) + &
        (M( 2,2) *        dx( ci)**2))
    END DO

  END SUBROUTINE calc_shape_functions_1D_stag_2nd_order

  SUBROUTINE calc_shape_functions_2D_reg_1st_order( x, y, n_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfx_c, Nfy_c)
    ! Calculate shape functions...
    ! ...in two dimensions...
    ! ...on the regular grid (i.e. f is known)...
    ! ...to 1st-order accuracy.
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x, y       ! The location where we want to know the gradients
    INTEGER,                             INTENT(IN)    :: n_max      ! The maximum number of surrounding points
    INTEGER,                             INTENT(IN)    :: n_c        ! The number  of     surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c, y_c   ! Coordinates of the surrounding points where we know f
    REAL(dp),                            INTENT(OUT)   :: Nfx_i      ! d/dx   shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfy_i      ! d/dy   shape function for the point [x,y]
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      ! d/dx   shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfy_c      ! d/dy   shape functions for the surrounding points

    ! Local variables:
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n_c)                           :: dx, dy, w
    REAL(dp), DIMENSION(2,2)                           :: ATWTWA, M

    ! Safety
    IF (n_c < 2) CALL crash('calc_shape_functions_2D_reg_1st_order needs at least 2 neighbours!')

    ! Calculate distances relative to [x,y]
    DO ci = 1, n_c
      dx( ci) = x_c( ci) - x
      dy( ci) = y_c( ci) - y
    END DO

    ! Calculate the weights w
    DO ci = 1, n_c
      w( ci) = 1._dp / (NORM2( [dx( ci), dy( ci)])**q)
    END DO

    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n_c
      ATWTWA( 1,1) = ATWTWA( 1,1) + w(ci)**2 *       dx( ci)    *       dx( ci)
      ATWTWA( 1,2) = ATWTWA( 1,2) + w(ci)**2 *       dx( ci)    *       dy( ci)

      ATWTWA( 2,1) = ATWTWA( 2,1) + w(ci)**2 *       dy( ci)    *       dx( ci)
      ATWTWA( 2,2) = ATWTWA( 2,2) + w(ci)**2 *       dy( ci)    *       dy( ci)
    END DO

    ! Invert ATWTWA to find M
    CALL calc_matrix_inverse_2_by_2( ATWTWA, M)

    ! Calculate shape functions
    Nfx_c   = 0._dp
    Nfy_c   = 0._dp
    DO ci = 1, n_c
      Nfx_c(   ci) = w( ci)**2 * ( &
        (M( 1,1) *        dx( ci)   ) + &
        (M( 1,2) *        dy( ci)   ))
      Nfy_c(   ci) = w( ci)**2 * ( &
        (M( 2,1) *        dx( ci)   ) + &
        (M( 2,2) *        dy( ci)   ))
    END DO

    Nfx_i  = -SUM( Nfx_c )
    Nfy_i  = -SUM( Nfy_c )

  END SUBROUTINE calc_shape_functions_2D_reg_1st_order

  SUBROUTINE calc_shape_functions_2D_reg_2nd_order( x, y, n_max, n_c, x_c, y_c, Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i, Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c)
    ! Calculate shape functions...
    ! ...in two dimensions...
    ! ...on the regular grid (i.e. f is known)...
    ! ...to 2nd-order accuracy.
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x, y       ! The location where we want to know the gradients
    INTEGER,                             INTENT(IN)    :: n_max      ! The maximum number of surrounding points
    INTEGER,                             INTENT(IN)    :: n_c        ! The number  of     surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c, y_c   ! Coordinates of the surrounding points where we know f
    REAL(dp),                            INTENT(OUT)   :: Nfx_i      ! d/dx    shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfy_i      ! d/dy    shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfxx_i     ! d2/dx2  shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfxy_i     ! d2/dxdy shape function for the point [x,y]
    REAL(dp),                            INTENT(OUT)   :: Nfyy_i     ! d2/dxy2 shape function for the point [x,y]
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      ! d/dx    shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfy_c      ! d/dy    shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfxx_c     ! d2/dx2  shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfxy_c     ! d2/dxdy shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfyy_c     ! d2/dy2  shape functions for the surrounding points

    ! Local variables:
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n_c)                           :: dx, dy, w
    REAL(dp), DIMENSION(5,5)                           :: ATWTWA, M

    ! Safety
    IF (n_c < 5) CALL crash('calc_shape_functions_2D_reg_2nd_order needs at least 2 neighbours!')

    ! Calculate distances relative to [x,y]
    DO ci = 1, n_c
      dx( ci) = x_c( ci) - x
      dy( ci) = y_c( ci) - y
    END DO

    ! Calculate the weights w
    DO ci = 1, n_c
      w( ci) = 1._dp / (NORM2( [dx( ci), dy( ci)])**q)
    END DO

    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n_c
      ATWTWA( 1,1) = ATWTWA( 1,1) + w(ci)**2 *       dx( ci)                 *       dx( ci)
      ATWTWA( 1,2) = ATWTWA( 1,2) + w(ci)**2 *       dx( ci)                 *                    dy( ci)
      ATWTWA( 1,3) = ATWTWA( 1,3) + w(ci)**2 *       dx( ci)                 * 1/2 * dx( ci)**2
      ATWTWA( 1,4) = ATWTWA( 1,4) + w(ci)**2 *       dx( ci)                 *       dx( ci)    * dy( ci)
      ATWTWA( 1,5) = ATWTWA( 1,5) + w(ci)**2 *       dx( ci)                 * 1/2 *              dy( ci)**2

      ATWTWA( 2,1) = ATWTWA( 2,1) + w(ci)**2 *                    dy( ci)    *       dx( ci)
      ATWTWA( 2,2) = ATWTWA( 2,2) + w(ci)**2 *                    dy( ci)    *                    dy( ci)
      ATWTWA( 2,3) = ATWTWA( 2,3) + w(ci)**2 *                    dy( ci)    * 1/2 * dx( ci)**2
      ATWTWA( 2,4) = ATWTWA( 2,4) + w(ci)**2 *                    dy( ci)    *       dx( ci)    * dy( ci)
      ATWTWA( 2,5) = ATWTWA( 2,5) + w(ci)**2 *                    dy( ci)    * 1/2 *              dy( ci)**2

      ATWTWA( 3,1) = ATWTWA( 3,1) + w(ci)**2 * 1/2 * dx( ci)**2              *       dx( ci)
      ATWTWA( 3,2) = ATWTWA( 3,2) + w(ci)**2 * 1/2 * dx( ci)**2              *                    dy( ci)
      ATWTWA( 3,3) = ATWTWA( 3,3) + w(ci)**2 * 1/2 * dx( ci)**2              * 1/2 * dx( ci)**2
      ATWTWA( 3,4) = ATWTWA( 3,4) + w(ci)**2 * 1/2 * dx( ci)**2              *       dx( ci)    * dy( ci)
      ATWTWA( 3,5) = ATWTWA( 3,5) + w(ci)**2 * 1/2 * dx( ci)**2              * 1/2 *              dy( ci)**2

      ATWTWA( 4,1) = ATWTWA( 4,1) + w(ci)**2 *       dx( ci)    * dy( ci)    *       dx( ci)
      ATWTWA( 4,2) = ATWTWA( 4,2) + w(ci)**2 *       dx( ci)    * dy( ci)    *                    dy( ci)
      ATWTWA( 4,3) = ATWTWA( 4,3) + w(ci)**2 *       dx( ci)    * dy( ci)    * 1/2 * dx( ci)**2
      ATWTWA( 4,4) = ATWTWA( 4,4) + w(ci)**2 *       dx( ci)    * dy( ci)    *       dx( ci)    * dy( ci)
      ATWTWA( 4,5) = ATWTWA( 4,5) + w(ci)**2 *       dx( ci)    * dy( ci)    * 1/2 *              dy( ci)**2

      ATWTWA( 5,1) = ATWTWA( 5,1) + w(ci)**2 * 1/2 *              dy( ci)**2 *       dx( ci)
      ATWTWA( 5,2) = ATWTWA( 5,2) + w(ci)**2 * 1/2 *              dy( ci)**2 *                    dy( ci)
      ATWTWA( 5,3) = ATWTWA( 5,3) + w(ci)**2 * 1/2 *              dy( ci)**2 * 1/2 * dx( ci)**2
      ATWTWA( 5,4) = ATWTWA( 5,4) + w(ci)**2 * 1/2 *              dy( ci)**2 *       dx( ci)    * dy( ci)
      ATWTWA( 5,5) = ATWTWA( 5,5) + w(ci)**2 * 1/2 *              dy( ci)**2 * 1/2 *              dy( ci)**2
    END DO

    ! Invert ATWTWA to find M
    CALL calc_matrix_inverse_general( ATWTWA, M)

    ! Calculate shape functions
    Nfx_c   = 0._dp
    Nfy_c   = 0._dp
    Nfxx_c  = 0._dp
    Nfxy_c  = 0._dp
    Nfyy_c  = 0._dp
    DO ci = 1, n_c
      Nfx_c(   ci) = w( ci)**2 * ( &
        (M( 1,1) *       dx( ci)                ) + &
        (M( 1,2) *                    dy( ci)   ) + &
        (M( 1,3) * 1/2 * dx( ci)**2             ) + &
        (M( 1,4) *       dx( ci)    * dy( ci)   ) + &
        (M( 1,5) * 1/2 *              dy( ci)**2))
      Nfy_c(   ci) = w( ci)**2 * ( &
        (M( 2,1) *       dx( ci)                ) + &
        (M( 2,2) *                    dy( ci)   ) + &
        (M( 2,3) * 1/2 * dx( ci)**2             ) + &
        (M( 2,4) *       dx( ci)    * dy( ci)   ) + &
        (M( 2,5) * 1/2 *              dy( ci)**2))
      Nfxx_c(   ci) = w( ci)**2 * ( &
        (M( 3,1) *       dx( ci)                ) + &
        (M( 3,2) *                    dy( ci)   ) + &
        (M( 3,3) * 1/2 * dx( ci)**2             ) + &
        (M( 3,4) *       dx( ci)    * dy( ci)   ) + &
        (M( 3,5) * 1/2 *              dy( ci)**2))
      Nfxy_c(   ci) = w( ci)**2 * ( &
        (M( 4,1) *       dx( ci)                ) + &
        (M( 4,2) *                    dy( ci)   ) + &
        (M( 4,3) * 1/2 * dx( ci)**2             ) + &
        (M( 4,4) *       dx( ci)    * dy( ci)   ) + &
        (M( 4,5) * 1/2 *              dy( ci)**2))
      Nfyy_c(   ci) = w( ci)**2 * ( &
        (M( 5,1) *       dx( ci)                ) + &
        (M( 5,2) *                    dy( ci)   ) + &
        (M( 5,3) * 1/2 * dx( ci)**2             ) + &
        (M( 5,4) *       dx( ci)    * dy( ci)   ) + &
        (M( 5,5) * 1/2 *              dy( ci)**2))
    END DO

    Nfx_i   = -SUM( Nfx_c  )
    Nfy_i   = -SUM( Nfy_c  )
    Nfxx_i  = -SUM( Nfxx_c )
    Nfxy_i  = -SUM( Nfxy_c )
    Nfyy_i  = -SUM( Nfyy_c )

  END SUBROUTINE calc_shape_functions_2D_reg_2nd_order

  SUBROUTINE calc_shape_functions_2D_stag_1st_order( x, y, n_max, n_c, x_c, y_c, Nf_c, Nfx_c, Nfy_c)
    ! Calculate shape functions...
    ! ...in two dimensions...
    ! ...on the staggered grid (i.e. f is not known)...
    ! ...to 1st-order accuracy.
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x, y       ! The location where we want to know the gradients
    INTEGER,                             INTENT(IN)    :: n_max      ! The maximum number of surrounding points
    INTEGER,                             INTENT(IN)    :: n_c        ! The number  of     surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c, y_c   ! Coordinates of the surrounding points where we know f
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nf_c       ! map    shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      ! d/dx   shape functions for the surrounding points
    REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfy_c      ! d/dy   shape functions for the surrounding points

    ! Local variables:
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n_c)                           :: dx, dy, w
    REAL(dp), DIMENSION(3,3)                           :: ATWTWA, M

    ! Safety
    IF (n_c < 3) CALL crash('calc_shape_functions_2D_stag_1st_order needs at least 3 neighbours!')

    ! Calculate distances relative to [x,y]
    DO ci = 1, n_c
      dx( ci) = x_c( ci) - x
      dy( ci) = y_c( ci) - y
    END DO

    ! Calculate the weights w
    DO ci = 1, n_c
      w( ci) = 1._dp / (NORM2( [dx( ci), dy( ci)])**q)
    END DO

    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n_c
      ATWTWA( 1,1) = ATWTWA( 1,1) + (w( ci)**2 * 1._dp   * 1._dp  )
      ATWTWA( 1,2) = ATWTWA( 1,2) + (w( ci)**2 * 1._dp   * dx( ci))
      ATWTWA( 1,3) = ATWTWA( 1,3) + (w( ci)**2 * 1._dp   * dy( ci))

      ATWTWA( 2,1) = ATWTWA( 2,1) + (w( ci)**2 * dx( ci) * 1._dp  )
      ATWTWA( 2,2) = ATWTWA( 2,2) + (w( ci)**2 * dx( ci) * dx( ci))
      ATWTWA( 2,3) = ATWTWA( 2,3) + (w( ci)**2 * dx( ci) * dy( ci))

      ATWTWA( 3,1) = ATWTWA( 3,1) + (w( ci)**2 * dy( ci) * 1._dp  )
      ATWTWA( 3,2) = ATWTWA( 3,2) + (w( ci)**2 * dy( ci) * dx( ci))
      ATWTWA( 3,3) = ATWTWA( 3,3) + (w( ci)**2 * dy( ci) * dy( ci))
    END DO

    ! Invert ATWTWA to find M
    CALL calc_matrix_inverse_3_by_3( ATWTWA, M)

    ! Calculate shape functions
    Nf_c    = 0._dp
    Nfx_c   = 0._dp
    Nfy_c   = 0._dp
    DO ci = 1, n_c
      Nf_c(  ci) = w( ci)**2 * ( &
        (M( 1,1) * 1._dp  ) + &
        (M( 1,2) * dx( ci)) + &
        (M( 1,3) * dy( ci)))
      Nfx_c(  ci) = w( ci)**2 * ( &
        (M( 2,1) * 1._dp  ) + &
        (M( 2,2) * dx( ci)) + &
        (M( 2,3) * dy( ci)))
      Nfy_c(  ci) = w( ci)**2 * ( &
        (M( 3,1) * 1._dp  ) + &
        (M( 3,2) * dx( ci)) + &
        (M( 3,3) * dy( ci)))
    END DO

  END SUBROUTINE calc_shape_functions_2D_stag_1st_order

! == Grid-cell-to-matrix-row translation

  ! Calculate translation tables
  SUBROUTINE calc_grid_cell_to_matrix_row_translation_tables( mesh)
    ! Calculate grid-cell-to-matrix-row translation tables

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_grid_cell_to_matrix_row_translation_tables'
    INTEGER                                            :: vi,ti,ci,k,ks,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Grid sizes
    mesh%nna     = mesh%nV
    mesh%nnauv   = mesh%nV              * 2
    mesh%nnak    = mesh%nV   *  C%nz
    mesh%nnakuv  = mesh%nV   *  C%nz    * 2
    mesh%nnaks   = mesh%nV   * (C%nz-1)
    mesh%nnaksuv = mesh%nV   * (C%nz-1) * 2

    mesh%nnb     = mesh%nTri
    mesh%nnbuv   = mesh%nTri            * 2
    mesh%nnbk    = mesh%nTri *  C%nz
    mesh%nnbkuv  = mesh%nTri *  C%nz    * 2
    mesh%nnbks   = mesh%nTri * (C%nz-1)
    mesh%nnbksuv = mesh%nTri * (C%nz-1) * 2

    mesh%nnc     = mesh%nAc
    mesh%nncuv   = mesh%nAc             * 2
    mesh%nnck    = mesh%nAc  *  C%nz
    mesh%nnckuv  = mesh%nAc  *  C%nz    * 2
    mesh%nncks   = mesh%nAc  * (C%nz-1)
    mesh%nncksuv = mesh%nAc  * (C%nz-1) * 2

    ! Allocate shared memory
    CALL allocate_shared_int_2D( mesh%nna    , 1, mesh%n2vi    , mesh%wn2vi    )
    CALL allocate_shared_int_2D( mesh%nnauv  , 2, mesh%n2viuv  , mesh%wn2viuv  )
    CALL allocate_shared_int_2D( mesh%nnak   , 2, mesh%n2vik   , mesh%wn2vik   )
    CALL allocate_shared_int_2D( mesh%nnakuv , 3, mesh%n2vikuv , mesh%wn2vikuv )
    CALL allocate_shared_int_2D( mesh%nnaks  , 2, mesh%n2viks  , mesh%wn2viks  )
    CALL allocate_shared_int_2D( mesh%nnaksuv, 3, mesh%n2viksuv, mesh%wn2viksuv)

    CALL allocate_shared_int_2D( mesh%nnb    , 1, mesh%n2ti    , mesh%wn2ti    )
    CALL allocate_shared_int_2D( mesh%nnbuv  , 2, mesh%n2tiuv  , mesh%wn2tiuv  )
    CALL allocate_shared_int_2D( mesh%nnbk   , 2, mesh%n2tik   , mesh%wn2tik   )
    CALL allocate_shared_int_2D( mesh%nnbkuv , 3, mesh%n2tikuv , mesh%wn2tikuv )
    CALL allocate_shared_int_2D( mesh%nnbks  , 2, mesh%n2tiks  , mesh%wn2tiks  )
    CALL allocate_shared_int_2D( mesh%nnbksuv, 3, mesh%n2tiksuv, mesh%wn2tiksuv)

    CALL allocate_shared_int_2D( mesh%nnc    , 1, mesh%n2ci    , mesh%wn2ci    )
    CALL allocate_shared_int_2D( mesh%nncuv  , 2, mesh%n2ciuv  , mesh%wn2ciuv  )
    CALL allocate_shared_int_2D( mesh%nnck   , 2, mesh%n2cik   , mesh%wn2cik   )
    CALL allocate_shared_int_2D( mesh%nnckuv , 3, mesh%n2cikuv , mesh%wn2cikuv )
    CALL allocate_shared_int_2D( mesh%nncks  , 2, mesh%n2ciks  , mesh%wn2ciks  )
    CALL allocate_shared_int_2D( mesh%nncksuv, 3, mesh%n2ciksuv, mesh%wn2ciksuv)

    CALL allocate_shared_int_1D( mesh%nV  ,            mesh%vi2n    , mesh%wvi2n    )
    CALL allocate_shared_int_2D( mesh%nV          , 2, mesh%viuv2n  , mesh%wviuv2n  )
    CALL allocate_shared_int_2D( mesh%nV  , C%nz  ,    mesh%vik2n   , mesh%wvik2n   )
    CALL allocate_shared_int_3D( mesh%nV  , C%nz  , 2, mesh%vikuv2n , mesh%wvikuv2n )
    CALL allocate_shared_int_2D( mesh%nV  , C%nz-1,    mesh%viks2n  , mesh%wviks2n  )
    CALL allocate_shared_int_3D( mesh%nV  , C%nz-1, 2, mesh%viksuv2n, mesh%wviksuv2n)

    CALL allocate_shared_int_1D( mesh%nTri,            mesh%ti2n    , mesh%wti2n    )
    CALL allocate_shared_int_2D( mesh%nTri        , 2, mesh%tiuv2n  , mesh%wtiuv2n  )
    CALL allocate_shared_int_2D( mesh%nTri, C%nz  ,    mesh%tik2n   , mesh%wtik2n   )
    CALL allocate_shared_int_3D( mesh%nTri, C%nz  , 2, mesh%tikuv2n , mesh%wtikuv2n )
    CALL allocate_shared_int_2D( mesh%nTri, C%nz-1,    mesh%tiks2n  , mesh%wtiks2n  )
    CALL allocate_shared_int_3D( mesh%nTri, C%nz-1, 2, mesh%tiksuv2n, mesh%wtiksuv2n)

    CALL allocate_shared_int_1D( mesh%nAc ,            mesh%ci2n    , mesh%wci2n    )
    CALL allocate_shared_int_2D( mesh%nAc         , 2, mesh%ciuv2n  , mesh%wciuv2n  )
    CALL allocate_shared_int_2D( mesh%nAc , C%nz  ,    mesh%cik2n   , mesh%wcik2n   )
    CALL allocate_shared_int_3D( mesh%nAc , C%nz  , 2, mesh%cikuv2n , mesh%wcikuv2n )
    CALL allocate_shared_int_2D( mesh%nAc , C%nz-1,    mesh%ciks2n  , mesh%wciks2n  )
    CALL allocate_shared_int_3D( mesh%nAc , C%nz-1, 2, mesh%ciksuv2n, mesh%wciksuv2n)

    ! Let the Master do the work
    IF (par%master) THEN

! == a-grid (vertices)

  ! == 2-D

    ! == scalar

      n = 0
      DO vi = 1, mesh%nV
        n = n+1
        mesh%vi2n( vi) = n
        mesh%n2vi( n,1) = vi
      END DO

    ! == vector

      n = 0
      DO vi = 1, mesh%nV
        DO uv = 1, 2
          n = n+1
          mesh%viuv2n( vi,uv) = n
          mesh%n2viuv( n,1) = vi
          mesh%n2viuv( n,2) = uv
        END DO
      END DO

  ! == 3-D regular

    ! == scalar

      n = 0
      DO vi = 1, mesh%nV
        DO k = 1, C%nz
          n = n+1
          mesh%vik2n( vi,k) = n
          mesh%n2vik( n,1) = vi
          mesh%n2vik( n,2) = k
        END DO
      END DO

    ! == vector

      n = 0
      DO vi = 1, mesh%nV
        DO k = 1, C%nz
          DO uv = 1, 2
            n = n+1
            mesh%vikuv2n( vi,k,uv) = n
            mesh%n2vikuv( n,1) = vi
            mesh%n2vikuv( n,2) = k
            mesh%n2vikuv( n,3) = uv
          END DO
        END DO
      END DO

  ! == 3-D staggered

    ! == scalar

      n = 0
      DO vi = 1, mesh%nV
        DO ks = 1, C%nz-1
          n = n+1
          mesh%viks2n( vi,ks) = n
          mesh%n2viks( n,1) = vi
          mesh%n2viks( n,2) = ks
        END DO
      END DO

    ! == vector

      n = 0
      DO vi = 1, mesh%nV
        DO ks = 1, C%nz-1
          DO uv = 1, 2
            n = n+1
            mesh%viksuv2n( vi,ks,uv) = n
            mesh%n2viksuv( n,1) = vi
            mesh%n2viksuv( n,2) = ks
            mesh%n2viksuv( n,3) = uv
          END DO
        END DO
      END DO

! == b-grid (triangles)

  ! == 2-D

    ! == scalar

      n = 0
      DO ti = 1, mesh%nTri
        n = n+1
        mesh%ti2n( ti) = n
        mesh%n2ti( n,1) = ti
      END DO

    ! == vector

      n = 0
      DO ti = 1, mesh%nTri
        DO uv = 1, 2
          n = n+1
          mesh%tiuv2n( ti,uv) = n
          mesh%n2tiuv( n,1) = ti
          mesh%n2tiuv( n,2) = uv
        END DO
      END DO

  ! == 3-D regular

    ! == scalar

      n = 0
      DO ti = 1, mesh%nTri
        DO k = 1, C%nz
          n = n+1
          mesh%tik2n( ti,k) = n
          mesh%n2tik( n,1) = ti
          mesh%n2tik( n,2) = k
        END DO
      END DO

    ! == vector

      n = 0
      DO ti = 1, mesh%nTri
        DO k = 1, C%nz
          DO uv = 1, 2
            n = n+1
            mesh%tikuv2n( ti,k,uv) = n
            mesh%n2tikuv( n,1) = ti
            mesh%n2tikuv( n,2) = k
            mesh%n2tikuv( n,3) = uv
          END DO
        END DO
      END DO

  ! == 3-D regustaggeredlar

    ! == scalar

      n = 0
      DO ti = 1, mesh%nTri
        DO ks = 1, C%nz-1
          n = n+1
          mesh%tiks2n( ti,ks) = n
          mesh%n2tiks( n,1) = ti
          mesh%n2tiks( n,2) = ks
        END DO
      END DO

    ! == vector

      n = 0
      DO ti = 1, mesh%nTri
        DO ks = 1, C%nz-1
          DO uv = 1, 2
            n = n+1
            mesh%tiksuv2n( ti,ks,uv) = n
            mesh%n2tiksuv( n,1) = ti
            mesh%n2tiksuv( n,2) = ks
            mesh%n2tiksuv( n,3) = uv
          END DO
        END DO
      END DO

! == c-grid (edges)

  ! == 2-D

    ! == scalar

      n = 0
      DO ci = 1, mesh%nAc
        n = n+1
        mesh%ci2n( ci) = n
        mesh%n2ci( n,1) = ci
      END DO

    ! == vector

      n = 0
      DO ci = 1, mesh%nAc
        DO uv = 1, 2
          n = n+1
          mesh%ciuv2n( ci,uv) = n
          mesh%n2ciuv( n,1) = ci
          mesh%n2ciuv( n,2) = uv
        END DO
      END DO

  ! == 3-D regular

    ! == scalar

      n = 0
      DO ci = 1, mesh%nAc
        DO k = 1, C%nz
          n = n+1
          mesh%cik2n( ci,k) = n
          mesh%n2cik( n,1) = ci
          mesh%n2cik( n,2) = k
        END DO
      END DO

    ! == vector

      n = 0
      DO ci = 1, mesh%nAc
        DO k = 1, C%nz
          DO uv = 1, 2
            n = n+1
            mesh%cikuv2n( ci,k,uv) = n
            mesh%n2cikuv( n,1) = ci
            mesh%n2cikuv( n,2) = k
            mesh%n2cikuv( n,3) = uv
          END DO
        END DO
      END DO

  ! == 3-D staggered

    ! == scalar

      n = 0
      DO ci = 1, mesh%nAc
        DO ks = 1, C%nz-1
          n = n+1
          mesh%ciks2n( ci,ks) = n
          mesh%n2ciks( n,1) = ci
          mesh%n2ciks( n,2) = ks
        END DO
      END DO

    ! == vector

      n = 0
      DO ci = 1, mesh%nAc
        DO ks = 1, C%nz-1
          DO uv = 1, 2
            n = n+1
            mesh%ciksuv2n( ci,ks,uv) = n
            mesh%n2ciksuv( n,1) = ci
            mesh%n2ciksuv( n,2) = ks
            mesh%n2ciksuv( n,3) = uv
          END DO
        END DO
      END DO

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_grid_cell_to_matrix_row_translation_tables

  ! Calculate matrices representing mapping operators between the scalar and vector b grids
  SUBROUTINE calc_buv_matrices( mesh)
    ! Calculate matrices representing mapping operators between the scalar and vector b grids

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_buv_matrices'
    TYPE(tMat)                                         :: Au, Av, A_bu_b, A_bv_b, A_bu_bu, A_bu_bv, A_bv_bu, A_bv_bv

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_buv_matrices_buv_b(   mesh)
    CALL calc_buv_matrices_b_buv(   mesh)
    CALL calc_buv_matrices_bkuv_bk( mesh)
    CALL calc_buv_matrices_bk_bkuv( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_buv_matrices

  SUBROUTINE calc_buv_matrices_buv_b( mesh)
    ! Calculate matrices representing mapping operators between the scalar and vector b grids

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_buv_matrices_buv_b'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map_bu_b, M_map_bv_b
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nnbuv     ! from
    nrows           = mesh%nnb       ! to
    nnz_per_row_est = 1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map_bu_b, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_map_bv_b, nrows, ncols, nnz_est_proc)

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2
      ti = mesh%n2ti( row,1)
      CALL add_entry_CSR_dist( M_map_bu_b, row, mesh%tiuv2n( ti,1), 1._dp)
      CALL add_entry_CSR_dist( M_map_bv_b, row, mesh%tiuv2n( ti,2), 1._dp)
    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrix
    CALL finalise_matrix_CSR_dist( M_map_bu_b, row1, row2)
    CALL finalise_matrix_CSR_dist( M_map_bv_b, row1, row2)

    ! Convert matrix from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map_bu_b, mesh%M_map_bu_b)
    CALL mat_CSR2petsc( M_map_bv_b, mesh%M_map_bv_b)

    ! Clean up the Fortran version
    CALL deallocate_matrix_CSR( M_map_bu_b)
    CALL deallocate_matrix_CSR( M_map_bv_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_buv_matrices_buv_b

  SUBROUTINE calc_buv_matrices_bkuv_bk( mesh)
    ! Calculate matrices representing mapping operators between the scalar and vector b grids

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_buv_matrices_bkuv_bk'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map_bku_bk, M_map_bkv_bk
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: ti,k

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nnbkuv     ! from
    nrows           = mesh%nnbk       ! to
    nnz_per_row_est = 1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map_bku_bk, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_map_bkv_bk, nrows, ncols, nnz_est_proc)

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2
      ti = mesh%n2tik( row,1)
      k  = mesh%n2tik( row,2)
      CALL add_entry_CSR_dist( M_map_bku_bk, row, mesh%tikuv2n( ti,k,1), 1._dp)
      CALL add_entry_CSR_dist( M_map_bkv_bk, row, mesh%tikuv2n( ti,k,2), 1._dp)
    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrix
    CALL finalise_matrix_CSR_dist( M_map_bku_bk, row1, row2)
    CALL finalise_matrix_CSR_dist( M_map_bkv_bk, row1, row2)

    ! Convert matrix from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map_bku_bk, mesh%M_map_bku_bk)
    CALL mat_CSR2petsc( M_map_bkv_bk, mesh%M_map_bkv_bk)

    ! Clean up the Fortran version
    CALL deallocate_matrix_CSR( M_map_bku_bk)
    CALL deallocate_matrix_CSR( M_map_bkv_bk)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_buv_matrices_bkuv_bk

  SUBROUTINE calc_buv_matrices_b_buv( mesh)
    ! Calculate matrices representing mapping operators between the scalar and vector b grids

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_buv_matrices_b_buv'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map_b_bu, M_map_b_bv
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: ti, uv

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nnb       ! from
    nrows           = mesh%nnbuv     ! to
    nnz_per_row_est = 1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map_b_bu, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_map_b_bv, nrows, ncols, nnz_est_proc)

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2
      ti = mesh%n2tiuv( row,1)
      uv = mesh%n2tiuv( row,2)
      IF (uv == 1) CALL add_entry_CSR_dist( M_map_b_bu, row, mesh%ti2n( ti), 1._dp)
      IF (uv == 2) CALL add_entry_CSR_dist( M_map_b_bv, row, mesh%ti2n( ti), 1._dp)
    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrix
    CALL finalise_matrix_CSR_dist( M_map_b_bu, row1, row2)
    CALL finalise_matrix_CSR_dist( M_map_b_bv, row1, row2)

    ! Convert matrix from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map_b_bu, mesh%M_map_b_bu)
    CALL mat_CSR2petsc( M_map_b_bv, mesh%M_map_b_bv)

    ! Clean up the Fortran version
    CALL deallocate_matrix_CSR( M_map_b_bu)
    CALL deallocate_matrix_CSR( M_map_b_bv)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_buv_matrices_b_buv

  SUBROUTINE calc_buv_matrices_bk_bkuv( mesh)
    ! Calculate matrices representing mapping operators between the scalar and vector b grids

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_buv_matrices_bk_bkuv'
    INTEGER                                            :: ncols, nrows, nnz_per_row_est, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                    :: M_map_bk_bku, M_map_bk_bkv
    INTEGER                                            :: row1, row2, row
    INTEGER                                            :: ti, k, uv

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the matrices using the native UFEMISM CSR-matrix format
  ! =====================================================================

    ! Matrix size
    ncols           = mesh%nnbk       ! from
    nrows           = mesh%nnbkuv     ! to
    nnz_per_row_est = 1
    nnz_est_proc    = CEILING( REAL( nnz_per_row_est * nrows, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map_bk_bku, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_map_bk_bkv, nrows, ncols, nnz_est_proc)

    ! Parallelisation
    CALL partition_list( nrows, par%i, par%n, row1, row2)

    DO row = row1, row2
      ti = mesh%n2tikuv( row,1)
      k  = mesh%n2tikuv( row,2)
      uv = mesh%n2tikuv( row,3)
      IF (uv == 1) CALL add_entry_CSR_dist( M_map_bk_bku, row, mesh%tik2n( ti,k), 1._dp)
      IF (uv == 2) CALL add_entry_CSR_dist( M_map_bk_bkv, row, mesh%tik2n( ti,k), 1._dp)
    END DO ! DO row = row1, row2
    CALL sync

    ! Assemble matrix
    CALL finalise_matrix_CSR_dist( M_map_bk_bku, row1, row2)
    CALL finalise_matrix_CSR_dist( M_map_bk_bkv, row1, row2)

    ! Convert matrix from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map_bk_bku, mesh%M_map_bk_bku)
    CALL mat_CSR2petsc( M_map_bk_bkv, mesh%M_map_bk_bkv)

    ! Clean up the Fortran version
    CALL deallocate_matrix_CSR( M_map_bk_bku)
    CALL deallocate_matrix_CSR( M_map_bk_bkv)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_buv_matrices_bk_bkuv

  ! Transform data between field form and vector form

  ! a-grid, field form to vector form
  SUBROUTINE field2vec_a(     mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_a'
    INTEGER                                            :: vi,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
      n = mesh%vi2n( vi)
      d_vec( n) = d( vi)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_a

  SUBROUTINE field2vec_ak(    mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_ak'
    INTEGER                                            :: vi,k,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV .OR. SIZE( d,2) /= C%nz) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, C%nz
      n = mesh%vik2n( vi,k)
      d_vec( n) = d( vi,k)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_ak

  SUBROUTINE field2vec_aks(   mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_aks'
    INTEGER                                            :: vi,ks,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV .OR. SIZE( d,2) /= C%nz-1) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
    DO ks = 1, C%nz-1
      n = mesh%viks2n( vi,ks)
      d_vec( n) = d( vi,ks)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_aks

  SUBROUTINE field2vec_auv(   mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_auv'
    INTEGER                                            :: vi,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV .OR. SIZE( d,2) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
    DO uv = 1, 2
      n = mesh%viuv2n( vi,uv)
      d_vec( n) = d( vi,uv)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_auv

  SUBROUTINE field2vec_akuv(  mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_akuv'
    INTEGER                                            :: vi,k,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV .OR. SIZE( d,2) /= C%nz .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, C%nz
    DO uv = 1, 2
      n = mesh%vikuv2n( vi,k,uv)
      d_vec( n) = d( vi,k,uv)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_akuv

  SUBROUTINE field2vec_aksuv( mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_aksuv'
    INTEGER                                            :: vi,ks,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV .OR. SIZE( d,2) /= C%nz-1 .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
    DO ks = 1, C%nz-1
    DO uv = 1, 2
      n = mesh%viksuv2n( vi,ks,uv)
      d_vec( n) = d( vi,ks,uv)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_aksuv

  ! a-grid, vector form to field form
  SUBROUTINE vec2field_a(     mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_a'
    INTEGER                                            :: vi,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
      n = mesh%vi2n( vi)
      d( vi) = d_vec( n)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_a

  SUBROUTINE vec2field_ak(    mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_ak'
    INTEGER                                            :: vi,k,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV .OR. SIZE( d,2) /= C%nz) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, C%nz
      n = mesh%vik2n( vi,k)
      d( vi,k) = d_vec( n)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_ak

  SUBROUTINE vec2field_aks(   mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_aks'
    INTEGER                                            :: vi,ks,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV .OR. SIZE( d,2) /= C%nz-1) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
    DO ks = 1, C%nz-1
      n = mesh%viks2n( vi,ks)
      d( vi,ks) = d_vec( n)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_aks

  SUBROUTINE vec2field_auv(   mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_auv'
    INTEGER                                            :: vi,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV .OR. SIZE( d,2) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
    DO uv = 1, 2
      n = mesh%viuv2n( vi,uv)
      d( vi,uv) = d_vec( n)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_auv

  SUBROUTINE vec2field_akuv(  mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_akuv'
    INTEGER                                            :: vi,k,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV .OR. SIZE( d,2) /= C%nz .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, C%nz
    DO uv = 1, 2
      n = mesh%vikuv2n( vi,k,uv)
      d( vi,k,uv) = d_vec( n)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_akuv

  SUBROUTINE vec2field_aksuv( mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_aksuv'
    INTEGER                                            :: vi,ks,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV .OR. SIZE( d,2) /= C%nz-1 .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO vi = mesh%vi1, mesh%vi2
    DO ks = 1, C%nz-1
    DO uv = 1, 2
      n = mesh%viksuv2n( vi,ks,uv)
      d( vi,ks,uv) = d_vec( n)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_aksuv

  ! b-grid, field form to vector form
  SUBROUTINE field2vec_b(     mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_b'
    INTEGER                                            :: ti,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
      n = mesh%ti2n( ti)
      d_vec( n) = d( ti)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_b

  SUBROUTINE field2vec_bk(    mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_bk'
    INTEGER                                            :: ti,k,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri .OR. SIZE( d,2) /= C%nz) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz
      n = mesh%tik2n( ti,k)
      d_vec( n) = d( ti,k)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_bk

  SUBROUTINE field2vec_bks(   mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_bks'
    INTEGER                                            :: ti,ks,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri .OR. SIZE( d,2) /= C%nz-1) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO ks = 1, C%nz-1
      n = mesh%tiks2n( ti,ks)
      d_vec( n) = d( ti,ks)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_bks

  SUBROUTINE field2vec_buv(   mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_buv'
    INTEGER                                            :: ti,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri .OR. SIZE( d,2) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO uv = 1, 2
      n = mesh%tiuv2n( ti,uv)
      d_vec( n) = d( ti,uv)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_buv

  SUBROUTINE field2vec_bkuv(  mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_bkuv'
    INTEGER                                            :: ti,k,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri .OR. SIZE( d,2) /= C%nz .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz
    DO uv = 1, 2
      n = mesh%tikuv2n( ti,k,uv)
      d_vec( n) = d( ti,k,uv)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_bkuv

  SUBROUTINE field2vec_bksuv( mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_bksuv'
    INTEGER                                            :: ti,ks,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri .OR. SIZE( d,2) /= C%nz-1 .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO ks = 1, C%nz-1
    DO uv = 1, 2
      n = mesh%tiksuv2n( ti,ks,uv)
      d_vec( n) = d( ti,ks,uv)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_bksuv

  ! b-grid, vector form to field form
  SUBROUTINE vec2field_b(     mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_b'
    INTEGER                                            :: ti,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
      n = mesh%ti2n( ti)
      d( ti) = d_vec( n)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_b

  SUBROUTINE vec2field_bk(    mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_bk'
    INTEGER                                            :: ti,k,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri .OR. SIZE( d,2) /= C%nz) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz
      n = mesh%tik2n( ti,k)
      d( ti,k) = d_vec( n)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_bk

  SUBROUTINE vec2field_bks(   mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_bks'
    INTEGER                                            :: ti,ks,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri .OR. SIZE( d,2) /= C%nz-1) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO ks = 1, C%nz-1
      n = mesh%tiks2n( ti,ks)
      d( ti,ks) = d_vec( n)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_bks

  SUBROUTINE vec2field_buv(   mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_buv'
    INTEGER                                            :: ti,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri .OR. SIZE( d,2) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO uv = 1, 2
      n = mesh%tiuv2n( ti,uv)
      d( ti,uv) = d_vec( n)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_buv

  SUBROUTINE vec2field_bkuv(  mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_bkuv'
    INTEGER                                            :: ti,k,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri .OR. SIZE( d,2) /= C%nz .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO k = 1, C%nz
    DO uv = 1, 2
      n = mesh%tikuv2n( ti,k,uv)
      d( ti,k,uv) = d_vec( n)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_bkuv

  SUBROUTINE vec2field_bksuv( mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_bksuv'
    INTEGER                                            :: ti,ks,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nTri .OR. SIZE( d,2) /= C%nz-1 .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO ks = 1, C%nz-1
    DO uv = 1, 2
      n = mesh%tiksuv2n( ti,ks,uv)
      d( ti,ks,uv) = d_vec( n)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_bksuv

  ! c-grid, field form to vector form
  SUBROUTINE field2vec_c(     mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_c'
    INTEGER                                            :: ci,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
      n = mesh%ci2n( ci)
      d_vec( n) = d( ci)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_c

  SUBROUTINE field2vec_ck(    mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_ck'
    INTEGER                                            :: ci,k,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc .OR. SIZE( d,2) /= C%nz) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
    DO k = 1, C%nz
      n = mesh%cik2n( ci,k)
      d_vec( n) = d( ci,k)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_ck

  SUBROUTINE field2vec_cks(   mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_cks'
    INTEGER                                            :: ci,ks,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc .OR. SIZE( d,2) /= C%nz-1) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
    DO ks = 1, C%nz-1
      n = mesh%ciks2n( ci,ks)
      d_vec( n) = d( ci,ks)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_cks

  SUBROUTINE field2vec_cuv(   mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_cuv'
    INTEGER                                            :: ci,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc .OR. SIZE( d,2) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
    DO uv = 1, 2
      n = mesh%ciuv2n( ci,uv)
      d_vec( n) = d( ci,uv)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_cuv

  SUBROUTINE field2vec_ckuv(  mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_ckuv'
    INTEGER                                            :: ci,k,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc .OR. SIZE( d,2) /= C%nz .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
    DO k = 1, C%nz
    DO uv = 1, 2
      n = mesh%cikuv2n( ci,k,uv)
      d_vec( n) = d( ci,k,uv)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_ckuv

  SUBROUTINE field2vec_cksuv( mesh, d, d_vec)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_vec

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'field2vec_cksuv'
    INTEGER                                            :: ci,ks,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc .OR. SIZE( d,2) /= C%nz-1 .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
    DO ks = 1, C%nz-1
    DO uv = 1, 2
      n = mesh%ciksuv2n( ci,ks,uv)
      d_vec( n) = d( ci,ks,uv)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE field2vec_cksuv

  ! c-grid, vector form to field form
  SUBROUTINE vec2field_c(     mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_c'
    INTEGER                                            :: ci,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
      n = mesh%ci2n( ci)
      d( ci) = d_vec( n)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_c

  SUBROUTINE vec2field_ck(    mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_ck'
    INTEGER                                            :: ci,k,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc .OR. SIZE( d,2) /= C%nz) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
    DO k = 1, C%nz
      n = mesh%cik2n( ci,k)
      d( ci,k) = d_vec( n)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_ck

  SUBROUTINE vec2field_cks(   mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_cks'
    INTEGER                                            :: ci,ks,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc .OR. SIZE( d,2) /= C%nz-1) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
    DO ks = 1, C%nz-1
      n = mesh%ciks2n( ci,ks)
      d( ci,ks) = d_vec( n)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_cks

  SUBROUTINE vec2field_cuv(   mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_cuv'
    INTEGER                                            :: ci,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc .OR. SIZE( d,2) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
    DO uv = 1, 2
      n = mesh%ciuv2n( ci,uv)
      d( ci,uv) = d_vec( n)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_cuv

  SUBROUTINE vec2field_ckuv(  mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_ckuv'
    INTEGER                                            :: ci,k,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc .OR. SIZE( d,2) /= C%nz .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
    DO k = 1, C%nz
    DO uv = 1, 2
      n = mesh%cikuv2n( ci,k,uv)
      d( ci,k,uv) = d_vec( n)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_ckuv

  SUBROUTINE vec2field_cksuv( mesh, d_vec, d)
    ! Transform data between field form and vector form

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_vec
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'vec2field_cksuv'
    INTEGER                                            :: ci,ks,uv,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nAc .OR. SIZE( d,2) /= C%nz-1 .OR. SIZE( d,3) /= 2) THEN
      CALL crash('arrays are of the wrong size!')
    END IF

    ! Translate indexing form field form to vector form
    DO ci = mesh%ci1, mesh%ci2
    DO ks = 1, C%nz-1
    DO uv = 1, 2
      n = mesh%ciksuv2n( ci,ks,uv)
      d( ci,ks,uv) = d_vec( n)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE vec2field_cksuv

! == 1-D matrix operators on the zeta grid

  SUBROUTINE convert_vertical_operator_from_1D_to_3D_k_k( M_1D, n2ck, ck2n, M_3D)
    ! Take a matrix operator in the 1-D vertical column, and extrude it
    ! to apply the same operation to every column on the 2-D grid.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)              :: M_1D
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: n2ck
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: ck2n
    TYPE(tMat),                          INTENT(INOUT)           :: M_3D

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'convert_vertical_operator_from_1D_to_3D_k_k'
    INTEGER                                                      :: nc, nz
    INTEGER                                                      :: nnz_per_row_max,row,ii1,ii2,nnz_in_row
    INTEGER                                                      :: ncols, nrows, nnz_max, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                              :: M_3D_CSR
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: cols_1D
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: cols_3D
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: vals
    INTEGER                                                      :: row_3D, row_3D1, row_3D2
    INTEGER                                                      :: cc,k
    INTEGER                                                      :: j,col_1D,kn,col_3D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Number of grid cells in the horizontal and vertical dimensions
    nc = SIZE( ck2n,1)
    nz = SIZE( ck2n,2)

    ! Determine maximum number of non-zeros per row
    nnz_per_row_max = 0
    DO row = 1, M_1D%m
      ii1 = M_1D%ptr( row)
      ii2 = M_1D%ptr( row+1)-1
      nnz_in_row = ii2 + 1 - ii1
      nnz_per_row_max = MAX( nnz_per_row_max, nnz_in_row)
    END DO

    ! Fill 3-D operator matrix

    ncols   = nc * nz
    nrows   = nc * nz
    nnz_max = nrows * nnz_per_row_max
    nnz_est_proc = CEILING( REAL( nnz_max, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_3D_CSR, nrows, ncols, nnz_max)

    ! Shape functions
    ALLOCATE( cols_1D( ncols))
    ALLOCATE( cols_3D( ncols))
    ALLOCATE( vals(    ncols))

    CALL partition_list( nrows, par%i, par%n, row_3D1, row_3D2)
    DO row_3D = row_3D1, row_3D2

      ! This row in the 3-D operator matrix corresponds to this 3-D grid cell
      cc = n2ck( row_3D,1)
      k  = n2ck( row_3D,2)

      ! Read the row for layer ks from the 1-D operator matrix
      ii1 = M_1D%ptr( k)
      ii2 = M_1D%ptr( k+1)-1
      nnz_in_row = ii2 + 1 - ii1
      cols_1D( 1:nnz_in_row) = M_1D%index( ii1: ii2)
      vals(    1:nnz_in_row) = M_1D%val(   ii1: ii2)

      ! Convert column indices from 1-D to 3-D
      DO j = 1, nnz_in_row
        col_1D = cols_1D( j)
        ! This column corresponds to this 1-D regular layer
        kn = col_1D
        ! This 3-D grid cell corresponds to this column in the 3-D matrix operator
        col_3D = ck2n( cc,kn)
        cols_3D( j) = col_3D
      END DO

      ! Write entries to the 3-D matrix operator
      DO j = 1, nnz_in_row
        CALL add_entry_CSR_dist( M_3D_CSR, row_3D, cols_3D( j), vals( j))
      END DO

    END DO

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_3D_CSR, row_3D1, row_3D2)

    ! Convert matrices from Fortran to PETSc type
    CALL mat_CSR2petsc( M_3D_CSR, M_3D)

    ! Clean up after yourself
    DEALLOCATE( cols_1D)
    DEALLOCATE( cols_3D)
    DEALLOCATE( vals   )
    CALL deallocate_matrix_CSR( M_3D_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE convert_vertical_operator_from_1D_to_3D_k_k

  SUBROUTINE convert_vertical_operator_from_1D_to_3D_k_ks( M_1D, n2ck, ck2n, n2cks, cks2n, M_3D)
    ! Take a matrix operator in the 1-D vertical column, and extrude it
    ! to apply the same operation to every column on the 2-D grid.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)              :: M_1D
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: n2ck
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: ck2n
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: n2cks
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: cks2n
    TYPE(tMat),                          INTENT(INOUT)           :: M_3D

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'convert_vertical_operator_from_1D_to_3D_k_ks'
    INTEGER                                                      :: nc, nz
    INTEGER                                                      :: nnz_per_row_max,row,ii1,ii2,nnz_in_row
    INTEGER                                                      :: ncols, nrows, nnz_max, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                              :: M_3D_CSR
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: cols_1D
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: cols_3D
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: vals
    INTEGER                                                      :: row_3D, row_3D1, row_3D2
    INTEGER                                                      :: cc,ks
    INTEGER                                                      :: j,col_1D,k,col_3D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Number of grid cells in the horizontal and vertical dimensions
    nc = size( ck2n,1)
    nz = size( ck2n,2)

    ! Determine maximum number of non-zeros per row
    nnz_per_row_max = 0
    DO row = 1, M_1D%m
      ii1 = M_1D%ptr( row)
      ii2 = M_1D%ptr( row+1)-1
      nnz_in_row = ii2 + 1 - ii1
      nnz_per_row_max = MAX( nnz_per_row_max, nnz_in_row)
    END DO

    ! Fill 3-D operator matrices

    ncols = nc *  nz
    nrows = nc * (nz-1)
    nnz_max = nrows * nnz_per_row_max
    nnz_est_proc = CEILING( REAL( nnz_max, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_3D_CSR, nrows, ncols, nnz_max)

    ALLOCATE( cols_1D( ncols))
    ALLOCATE( cols_3D( ncols))
    ALLOCATE( vals(    ncols))

    CALL partition_list( nrows, par%i, par%n, row_3D1, row_3D2)
    DO row_3D = row_3D1, row_3D2

      ! This row in the 3-D operator matrix corresponds to this 3-D grid cell
      cc = n2cks( row_3D,1)
      ks = n2cks( row_3D,2)

      ! Read the row for layer ks from the 1-D operator matrix
      ii1 = M_1D%ptr( ks)
      ii2 = M_1D%ptr( ks+1)-1
      nnz_in_row = ii2 + 1 - ii1
      cols_1D( 1:nnz_in_row) = M_1D%index(  ii1: ii2)
      vals(    1:nnz_in_row) = M_1D%val(    ii1: ii2)

      ! Convert column indices from 1-D to 3-D
      DO j = 1, nnz_in_row
        col_1D = cols_1D( j)
        ! This column corresponds to this 1-D regular layer
        k = col_1D
        ! This 3-D grid cell corresponds to this column in the 3-D matrix operator
        col_3D = ck2n( cc,k)
        cols_3D( j) = col_3D
      END DO

      ! Write entries to the 3-D matrix operator
      DO j = 1, nnz_in_row
        CALL add_entry_CSR_dist( M_3D_CSR, row_3D, cols_3D( j), vals( j))
      END DO

    END DO

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_3D_CSR, row_3D1, row_3D2)

    ! Convert matrices from Fortran to PETSc type
    CALL mat_CSR2petsc( M_3D_CSR, M_3D)

    ! Clean up after yourself
    DEALLOCATE( cols_1D)
    DEALLOCATE( cols_3D)
    DEALLOCATE( vals   )
    CALL deallocate_matrix_CSR( M_3D_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE convert_vertical_operator_from_1D_to_3D_k_ks

  SUBROUTINE convert_vertical_operator_from_1D_to_3D_ks_k( M_1D, n2ck, ck2n, n2cks, cks2n, M_3D)
    ! Take a matrix operator in the 1-D vertical column, and extrude it
    ! to apply the same operation to every column on the 2-D grid.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)              :: M_1D
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: n2ck
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: ck2n
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: n2cks
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: cks2n
    TYPE(tMat),                          INTENT(INOUT)           :: M_3D

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'convert_vertical_operator_from_1D_to_3D_ks_k'
    INTEGER                                                      :: nc, nz
    INTEGER                                                      :: nnz_per_row_max,row,ii1,ii2,nnz_in_row
    INTEGER                                                      :: ncols, nrows, nnz_max, nnz_est_proc
    TYPE(type_sparse_matrix_CSR_dp)                              :: M_3D_CSR
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: cols_1D
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: cols_3D
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: vals
    INTEGER                                                      :: row_3D, row_3D1, row_3D2
    INTEGER                                                      :: cc,ks
    INTEGER                                                      :: j,col_1D,k,col_3D

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Number of grid cells in the horizontal and vertical dimensions
    nc = size( ck2n,1)
    nz = size( ck2n,2)

    ! Determine maximum number of non-zeros per row
    nnz_per_row_max = 0
    DO row = 1, M_1D%m
      ii1 = M_1D%ptr( row)
      ii2 = M_1D%ptr( row+1)-1
      nnz_in_row = ii2 + 1 - ii1
      nnz_per_row_max = MAX( nnz_per_row_max, nnz_in_row)
    END DO

    ! Fill 3-D operator matrices

    ncols = nc * (nz-1)
    nrows = nc * nz
    nnz_max = nrows * nnz_per_row_max

    CALL allocate_matrix_CSR_dist( M_3D_CSR, nrows, ncols, nnz_max)

    ALLOCATE( cols_1D( ncols))
    ALLOCATE( cols_3D( ncols))
    ALLOCATE( vals(    ncols))

    CALL partition_list( nrows, par%i, par%n, row_3D1, row_3D2)
    DO row_3D = row_3D1, row_3D2

      ! This row in the 3-D operator matrix corresponds to this 3-D grid cell
      cc = n2ck( row_3D,1)
      k  = n2ck( row_3D,2)

      ! Read the row for layer ks from the 1-D operator matrix
      ii1 = M_1D%ptr( k)
      ii2 = M_1D%ptr( k+1)-1
      nnz_in_row = ii2 + 1 - ii1
      cols_1D( 1:nnz_in_row) = M_1D%index( ii1: ii2)
      vals(    1:nnz_in_row) = M_1D%val(   ii1: ii2)

      ! Convert column indices from 1-D to 3-D
      DO j = 1, nnz_in_row
        col_1D = cols_1D( j)
        ! This column corresponds to this 1-D staggered layer
        ks = col_1D
        ! This 3-D grid cell corresponds to this column in the 3-D matrix operator
        col_3D = cks2n( cc,ks)
        cols_3D( j) = col_3D
      END DO

      ! Write entries to the 3-D matrix operator
      DO j = 1, nnz_in_row
        CALL add_entry_CSR_dist( M_3D_CSR, row_3D, cols_3D( j), vals( j))
      END DO

    END DO

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_3D_CSR, row_3D1, row_3D2)

    ! Convert matrices from Fortran to PETSc type
    CALL mat_CSR2petsc( M_3D_CSR, M_3D)

    ! Clean up after yourself
    DEALLOCATE( cols_1D)
    DEALLOCATE( vals   )
    DEALLOCATE( cols_3D)
    CALL deallocate_matrix_CSR( M_3D_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE convert_vertical_operator_from_1D_to_3D_ks_k

  SUBROUTINE calc_vertical_operators_reg_1D( zeta, M_ddzeta_k_k_1D, M_d2dzeta2_k_k_1D, M_fkp1_m_fk_1D, M_fk_m_fkm1_1D)
    ! Calculate mapping and d/dzeta operators in the 1-D vertical column
    !
    ! NOTE: since its well possible that the number of cores running the model
    !       exceeds the number of vertical layers, these matrices are stored in
    !       process-local memory

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)              :: zeta
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: M_ddzeta_k_k_1D
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: M_d2dzeta2_k_k_1D
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: M_fkp1_m_fk_1D
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: M_fk_m_fkm1_1D

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_vertical_operators_reg_1D'
    INTEGER                                                      :: nz, ncols, nrows, nnz_max
    INTEGER                                                      :: k
    REAL(dp)                                                     :: z
    INTEGER                                                      :: n_c
    INTEGER,  DIMENSION(2)                                       :: i_c
    REAL(dp), DIMENSION(2)                                       :: z_c
    INTEGER                                                      :: i,kn
    REAL(dp)                                                     :: Nfz_i, Nfzz_i
    REAL(dp), DIMENSION(2)                                       :: Nfz_c, Nfzz_c

    ! Add routine to path
    CALL init_routine( routine_name)

    nz = SIZE( zeta,1)

    ! k (regular) to k (regular)

    ncols   = nz
    nrows   = nz
    nnz_max = nrows * 3

    CALL allocate_matrix_CSR_dist( M_ddzeta_k_k_1D  , nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_d2dzeta2_k_k_1D, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_fkp1_m_fk_1D   , nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_fk_m_fkm1_1D   , nrows, ncols, nnz_max)

    DO k = 1, nz

      ! Source: layer k
      z = zeta( k)

      ! Neighbours: layers k-1 and k+1
      IF     (k == 1) THEN
        n_c = 2
        i_c = [2, 3]
        z_c = [zeta( 2), zeta( 3)]
      ELSEIF (k == nz) THEN
        n_c = 2
        i_c = [nz-2, nz-1]
        z_c = [zeta( nz-2), zeta( nz-1)]
      ELSE
        n_c = 2
        i_c = [k-1, k+1]
        z_c = [zeta( k-1), zeta( k+1)]
      END IF

      ! Calculate shape functions
      CALL calc_shape_functions_1D_reg_2nd_order( z, n_c, n_c, z_c, Nfz_i, Nfzz_i, Nfz_c, Nfzz_c)

      ! Diagonal element: shape function for the source point
      CALL add_entry_CSR_dist( M_ddzeta_k_k_1D  , k, k, Nfz_i )
      CALL add_entry_CSR_dist( M_d2dzeta2_k_k_1D, k, k, Nfzz_i)

      ! Off-diagonal elements: shape functions for neighbours
      DO i = 1, n_c
        kn = i_c( i)
        CALL add_entry_CSR_dist( M_ddzeta_k_k_1D  , k, kn, Nfz_c(  i))
        CALL add_entry_CSR_dist( M_d2dzeta2_k_k_1D, k, kn, Nfzz_c( i))
      END DO

      ! One-sided differencing
      IF (k > 1) THEN
        CALL add_entry_CSR_dist( M_fk_m_fkm1_1D, k, k-1, -1._dp)
        CALL add_entry_CSR_dist( M_fk_m_fkm1_1D, k, k  ,  1._dp)
      END IF
      IF (k < C%nz) THEN
        CALL add_entry_CSR_dist( M_fkp1_m_fk_1D, k, k  , -1._dp)
        CALL add_entry_CSR_dist( M_fkp1_m_fk_1D, k, k+1,  1._dp)
      END IF

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_operators_reg_1D

  SUBROUTINE calc_vertical_operators_stag_1D( zeta, zeta_stag, M_map_k_ks_1D, M_ddzeta_k_ks_1D, M_map_ks_k_1D, M_ddzeta_ks_k_1D)
    ! Calculate mapping and d/dzeta operators in the 1-D vertical column
    !
    ! NOTE: since its well possible that the number of cores running the model
    !       exceeds the number of vertical layers, these matrices are stored in
    !       process-local memory

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)              :: zeta
    REAL(dp), DIMENSION(:    ),          INTENT(IN)              :: zeta_stag
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: M_map_k_ks_1D
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: M_ddzeta_k_ks_1D
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: M_map_ks_k_1D
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT)           :: M_ddzeta_ks_k_1D

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_vertical_operators_stag_1D'
    INTEGER                                                      :: nz, ncols, nrows, nnz_max
    INTEGER                                                      :: ks
    REAL(dp)                                                     :: dzeta
    INTEGER                                                      :: k, row, col
    INTEGER                                                      ::  k_lo,  k_hi,  ks_lo,  ks_hi
    REAL(dp)                                                     :: wk_lo, wk_hi, wks_lo, wks_hi

    ! Add routine to path
    CALL init_routine( routine_name)

    nz = SIZE( zeta,1)

  ! == k (regular) to ks (staggered)

    ncols = nz
    nrows = nz-1
    nnz_max = nrows * 2

    CALL allocate_matrix_CSR_dist( M_map_k_ks_1D   , nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddzeta_k_ks_1D, nrows, ncols, nnz_max)

    DO ks = 1, nz-1

      ! Indices of neighbouring grid points
      k_lo = ks
      k_hi = ks + 1

      ! Local grid spacing
      dzeta = zeta( k_hi) - zeta( k_lo)

      ! Linear interpolation weights
      wk_lo = (zeta( k_hi) - zeta_stag( ks)) / dzeta
      wk_hi = 1._dp - wk_lo

      ! Left-hand neighbour
      row = ks
      col = k_lo
      CALL add_entry_CSR_dist( M_map_k_ks_1D   , row, col, wk_lo         )
      CALL add_entry_CSR_dist( M_ddzeta_k_ks_1D, row, col, -1._dp / dzeta)

      ! Right-hand neighbour
      row = ks
      col = k_hi
      CALL add_entry_CSR_dist( M_map_k_ks_1D   , row, col, wk_hi         )
      CALL add_entry_CSR_dist( M_ddzeta_k_ks_1D, row, col,  1._dp / dzeta)

    END DO

  ! == ks (staggered) to k (regular)

    ncols = nz-1
    nrows = nz
    nnz_max = nrows * 2

    CALL allocate_matrix_CSR_dist( M_map_ks_k_1D   , nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddzeta_ks_k_1D, nrows, ncols, nnz_max)

    DO k = 1, nz

      ! Indices of neighbouring grid points
      IF     (k == 1) THEN
        ks_lo = 1
        ks_hi = 2
      ELSEIF (k == nz) THEN
        ks_lo = nz - 2
        ks_hi = nz - 1
      ELSE
        ks_lo = k - 1
        ks_hi = k
      END IF

      ! Local grid spacing
      dzeta = zeta_stag( ks_hi) - zeta_stag( ks_lo)

      ! Linear interpolation weights
      wks_lo = (zeta_stag( ks_hi) - zeta( k)) / dzeta
      wks_hi = 1._dp - wks_lo

      ! Left-hand neighbour
      row = k
      col = ks_lo
      CALL add_entry_CSR_dist( M_map_ks_k_1D   , row, col, wks_lo        )
      CALL add_entry_CSR_dist( M_ddzeta_ks_k_1D, row, col, -1._dp / dzeta)

      ! Right-hand neighbour
      row = k
      col = ks_hi
      CALL add_entry_CSR_dist( M_map_ks_k_1D   , row, col, wks_hi        )
      CALL add_entry_CSR_dist( M_ddzeta_ks_k_1D, row, col,  1._dp / dzeta)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_operators_stag_1D

  SUBROUTINE calc_zeta_operators_tridiagonal( mesh)
    ! Calculate zeta operators in tridiagonal form for efficient use in thermodynamics

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_zeta_operators_tridiagonal'
    INTEGER                                                      :: k,i,ii1,ii2,ii,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for the three diagonals of both matrices
    ALLOCATE( mesh%M_ddzeta_k_k_ldiag(   C%nz-1))
    ALLOCATE( mesh%M_ddzeta_k_k_diag(    C%nz  ))
    ALLOCATE( mesh%M_ddzeta_k_k_udiag(   C%nz-1))
    ALLOCATE( mesh%M_d2dzeta2_k_k_ldiag( C%nz-1))
    ALLOCATE( mesh%M_d2dzeta2_k_k_diag(  C%nz  ))
    ALLOCATE( mesh%M_d2dzeta2_k_k_udiag( C%nz-1))

    ! Fill in coefficients

    ! d/dzeta
    DO k = 1, C%nz

      ! Lower diagonal
      IF (k > 1) THEN
        ! Initialise
        mesh%M_ddzeta_k_k_ldiag( k-1) = 0._dp
        ! Find matrix element at [k,k-1]
        i = k
        ii1 = mesh%M_ddzeta_k_k_1D%ptr( i)
        ii2 = mesh%M_ddzeta_k_k_1D%ptr( i+1)-1
        DO ii = ii1, ii2
          j = mesh%M_ddzeta_k_k_1D%index( ii)
          IF (j == k-1) mesh%M_ddzeta_k_k_ldiag( k-1) = mesh%M_ddzeta_k_k_1D%val( ii)
        END DO
      END IF ! IF (k > 1) THEN

      ! Central diagonal
      ! Initialise
      mesh%M_ddzeta_k_k_diag( k) = 0._dp
      ! Find matrix element at [k,k-1]
      i = k
      ii1 = mesh%M_ddzeta_k_k_1D%ptr( i)
      ii2 = mesh%M_ddzeta_k_k_1D%ptr( i+1)-1
      DO ii = ii1, ii2
        j = mesh%M_ddzeta_k_k_1D%index( ii)
        IF (j == k) mesh%M_ddzeta_k_k_diag( k) = mesh%M_ddzeta_k_k_1D%val( ii)
      END DO

      ! Upper diagonal
      IF (k < C%nz) THEN
        ! Initialise
        mesh%M_ddzeta_k_k_udiag( k) = 0._dp
        ! Find matrix element at [k,k-1]
        i = k
        ii1 = mesh%M_ddzeta_k_k_1D%ptr( i)
        ii2 = mesh%M_ddzeta_k_k_1D%ptr( i+1)-1
        DO ii = ii1, ii2
          j = mesh%M_ddzeta_k_k_1D%index( ii)
          IF (j == k+1) mesh%M_ddzeta_k_k_udiag( k) = mesh%M_ddzeta_k_k_1D%val( ii)
        END DO
      END IF ! IF (k > 1) THEN

    END DO ! DO k = 1, C%nz

    ! d2/dzeta2
    DO k = 1, C%nz

      ! Lower diagonal
      IF (k > 1) THEN
        ! Initialise
        mesh%M_d2dzeta2_k_k_ldiag( k-1) = 0._dp
        ! Find matrix element at [k,k-1]
        i = k
        ii1 = mesh%M_d2dzeta2_k_k_1D%ptr( i)
        ii2 = mesh%M_d2dzeta2_k_k_1D%ptr( i+1)-1
        DO ii = ii1, ii2
          j = mesh%M_d2dzeta2_k_k_1D%index( ii)
          IF (j == k-1) mesh%M_d2dzeta2_k_k_ldiag( k-1) = mesh%M_d2dzeta2_k_k_1D%val( ii)
        END DO
      END IF ! IF (k > 1) THEN

      ! Central diagonal
      ! Initialise
      mesh%M_d2dzeta2_k_k_diag( k) = 0._dp
      ! Find matrix element at [k,k-1]
      i = k
      ii1 = mesh%M_d2dzeta2_k_k_1D%ptr( i)
      ii2 = mesh%M_d2dzeta2_k_k_1D%ptr( i+1)-1
      DO ii = ii1, ii2
        j = mesh%M_d2dzeta2_k_k_1D%index( ii)
        IF (j == k) mesh%M_d2dzeta2_k_k_diag( k) = mesh%M_d2dzeta2_k_k_1D%val( ii)
      END DO

      ! Upper diagonal
      IF (k < C%nz) THEN
        ! Initialise
        mesh%M_d2dzeta2_k_k_udiag( k) = 0._dp
        ! Find matrix element at [k,k-1]
        i = k
        ii1 = mesh%M_d2dzeta2_k_k_1D%ptr( i)
        ii2 = mesh%M_d2dzeta2_k_k_1D%ptr( i+1)-1
        DO ii = ii1, ii2
          j = mesh%M_d2dzeta2_k_k_1D%index( ii)
          IF (j == k+1) mesh%M_d2dzeta2_k_k_udiag( k) = mesh%M_d2dzeta2_k_k_1D%val( ii)
        END DO
      END IF ! IF (k > 1) THEN

    END DO ! DO k = 1, C%nz

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_zeta_operators_tridiagonal

! == Extend 2-D matrix to 3-D

  SUBROUTINE extend_matrix_operator_2D_to_3D( row2c, c2row, col2c, c2col, row2ck, ck2row, col2ck, ck2col, M_2D_mat, M_3D_mat)
    ! Take a matrix operator that operators on the 2-D plane, and extrude it vertically so that
    ! it operators on every vertical plane in the zeta grid.
    !
    ! E.g.: take the d/dx operators M_ddx_a_a, and use it to construct M_ddxp_a_a

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: row2c
    INTEGER,  DIMENSION(:    ),          INTENT(IN)              :: c2row
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: col2c
    INTEGER,  DIMENSION(:    ),          INTENT(IN)              :: c2col
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: row2ck
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: ck2row
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: col2ck
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)              :: ck2col
    TYPE(tMat),                          INTENT(IN)              :: M_2D_mat
    TYPE(tMat),                          INTENT(INOUT)           :: M_3D_mat

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'extend_matrix_operator_2D_to_3D'
    TYPE(type_sparse_matrix_CSR_dp)                              :: M_2D, M_3D
    INTEGER                                                      :: nrows_2D, ncols_2D, nnz_2D, row_2D, col_2D
    INTEGER                                                      :: nrows_3D, ncols_3D, nnz_3D, row_3D, col_3D
    INTEGER                                                      :: nz
    INTEGER                                                      :: nnz_est_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                      :: cols_2D, cols_3D
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                      :: vals
    INTEGER                                                      :: row_3D1, row_3D2
    INTEGER                                                      :: cc, k, ii1, ii2, nnz_in_row, j, ccc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Convert to CSR format
    CALL mat_petsc2CSR( M_2D_mat, M_2D)

    ! Size of the 2-D operator
    nrows_2D = M_2D%m
    ncols_2D = M_2D%n
    nnz_2D   = M_2D%nnz

    ! Number of vertical layers
    nz = SIZE( ck2row,2)

    ! Size of the 3-D operator
    nrows_3D = nrows_2D * nz
    ncols_3D = ncols_2D * nz
    nnz_3D   = nnz_2D   * nz

    ! Allocate memory for the 3-D operator
    nnz_est_proc = CEILING( REAL( nnz_3D, dp) / REAL( par%n, dp))
    CALL allocate_matrix_CSR_dist( M_3D, nrows_3D, ncols_3D, nnz_est_proc)

    ! Shape functions
    ALLOCATE( cols_2D( ncols_2D))
    ALLOCATE( vals(    ncols_2D))
    ALLOCATE( cols_3D( ncols_2D))

    ! Fill the 3-D matrix
    CALL partition_list( nrows_3D, par%i, par%n, row_3D1, row_3D2)
    DO row_3D = row_3D1, row_3D2

      ! This row corresponds to this 3-D grid cell
      cc = row2ck( row_3D,1)
      k  = row2ck( row_3D,2)

      ! 2-D grid cell cc corresponds to this row in the 2-D matrix operator
      row_2D = c2row( cc);

      ! Read this row from the 2-D matrix operator
      ii1 = M_2D%ptr( row_2D)
      ii2 = M_2D%ptr( row_2D+1)-1
      nnz_in_row = ii2 + 1 - ii1
      IF (nnz_in_row == 0) CYCLE

      cols_2D( 1:nnz_in_row) = M_2D%index( ii1: ii2)
      vals(    1:nnz_in_row) = M_2D%val(   ii1: ii2)

      ! Convert column indices from 2-D to 3-D
      DO j = 1, nnz_in_row
        col_2D = cols_2D( j)
        ! This column corresponds to this 2-D grid cell
        ccc = col2c( col_2D,1)
        ! This 3-D grid cell corresponds to this column in the 3-D matrix operator
        col_3D = ck2col( ccc,k)
        cols_3D( j) = col_3D
      END DO

      ! Write entries to the 3-D matrix operator
      DO j = 1, nnz_in_row
        CALL add_entry_CSR_dist( M_3D, row_3D, cols_3D( j), vals( j))
      END DO

    END DO

    ! Assemble matrix
    CALL finalise_matrix_CSR_dist( M_3D, row_3D1, row_3D2)

    ! Convert matrix from Fortran to PETSc type
    CALL mat_CSR2petsc( M_3D, M_3D_mat)

    ! Clean up after yourself
    DEALLOCATE( cols_2D)
    DEALLOCATE( vals   )
    DEALLOCATE( cols_3D)
    CALL deallocate_matrix_CSR( M_3D)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extend_matrix_operator_2D_to_3D

! == Calculate mapping operators between 2-D grid and 3-D grid surface/base layers

  SUBROUTINE calc_maps_b_to_bk_surf_base( mesh)
    ! Calculate mapping operators between 2-D grid and 3-D grid surface/base layers
    !
    ! Used for efficiently constructing stiffness matrices describing surface/base
    ! boundary conditions to the Blatter-Pattyn Approximation.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT)           :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'calc_maps_b_to_bk_surf_base'
    TYPE(type_sparse_matrix_CSR_dp)                              :: M_map_b_bk_surf_CSR
    TYPE(type_sparse_matrix_CSR_dp)                              :: M_map_b_bk_base_CSR
    TYPE(type_sparse_matrix_CSR_dp)                              :: M_map_bk_surf_b_CSR
    TYPE(type_sparse_matrix_CSR_dp)                              :: M_map_bk_base_b_CSR
    INTEGER                                                      :: ncols, nrows, nnz_max, nnz_est_proc
    INTEGER                                                      :: row_3D, row_3D1, row_3D2
    INTEGER                                                      :: row_2D, row_2D1, row_2D2
    INTEGER                                                      :: ti,k

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == b to bk

    ncols   = mesh%nnb
    nrows   = mesh%nnbk
    nnz_max = ncols
    nnz_est_proc = CEILING( REAL( nnz_max, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map_b_bk_surf_CSR, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_map_b_bk_base_CSR, nrows, ncols, nnz_est_proc)

    CALL partition_list( nrows, par%i, par%n, row_3D1, row_3D2)
    DO row_3D = row_3D1, row_3D2

      ! This matrix represents this 3-D grid cell
      ti = mesh%n2tik( row_3D,1)
      k  = mesh%n2tik( row_3D,2)

      ! This 2-D grid cell corresponds to this matrix row
      row_2D = mesh%ti2n( ti)

      IF     (k == 1) THEN
        CALL add_entry_CSR_dist( M_map_b_bk_surf_CSR, row_3D, row_2D, 1._dp)
      ELSEIF (k == C%nz) THEN
        CALL add_entry_CSR_dist( M_map_b_bk_base_CSR, row_3D, row_2D, 1._dp)
      END IF

    END DO

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_map_b_bk_surf_CSR, row_3D1, row_3D2)
    CALL finalise_matrix_CSR_dist( M_map_b_bk_base_CSR, row_3D1, row_3D2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map_b_bk_surf_CSR, mesh%M_map_b_bk_surf)
    CALL mat_CSR2petsc( M_map_b_bk_base_CSR, mesh%M_map_b_bk_base)

    ! Clean up after yourself
    CALL deallocate_matrix_CSR( M_map_b_bk_surf_CSR)
    CALL deallocate_matrix_CSR( M_map_b_bk_base_CSR)

  ! == bk to b

    ncols   = mesh%nnbk
    nrows   = mesh%nnb
    nnz_max = ncols
    nnz_est_proc = CEILING( REAL( nnz_max, dp) / REAL( par%n, dp))

    CALL allocate_matrix_CSR_dist( M_map_bk_surf_b_CSR, nrows, ncols, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( M_map_bk_base_b_CSR, nrows, ncols, nnz_est_proc)

    CALL partition_list( nrows, par%i, par%n, row_2D1, row_2D2)
    DO row_2D = row_2D1, row_2D2

      ! This matrix represents this 2-D grid cell
      ti = mesh%n2ti( row_2D,1)

      ! The 3-D surface grid cell corresponds to this matrix row
      row_3D = mesh%tik2n( ti,1)
      CALL add_entry_CSR_dist( M_map_bk_surf_b_CSR, row_2D, row_3D, 1._dp)

      ! The 3-D base grid cell corresponds to this matrix row
      row_3D = mesh%tik2n( ti,C%nz)
      CALL add_entry_CSR_dist( M_map_bk_base_b_CSR, row_2D, row_3D, 1._dp)

    END DO

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( M_map_bk_surf_b_CSR, row_2D1, row_2D2)
    CALL finalise_matrix_CSR_dist( M_map_bk_base_b_CSR, row_2D1, row_2D2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( M_map_bk_surf_b_CSR, mesh%M_map_bk_surf_b)
    CALL mat_CSR2petsc( M_map_bk_base_b_CSR, mesh%M_map_bk_base_b)

    ! Clean up after yourself
    CALL deallocate_matrix_CSR( M_map_bk_surf_b_CSR)
    CALL deallocate_matrix_CSR( M_map_bk_base_b_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_maps_b_to_bk_surf_base

! == Calculate d/dx, d/dy operators on a square grid (used in grid-to-mesh remapping)

  SUBROUTINE calc_matrix_operators_grid( grid, M_ddx, M_ddy)
    ! Calculate matrix operators for partial derivatives on a regular grid (needed for conservative remapping)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(tMat),                          INTENT(INOUT) :: M_ddx, M_ddy

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_matrix_operators_grid'
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER                                            :: n, i, j, n_ip1, n_im1, n_jp1, n_jm1
    REAL(dp)                                           :: v1, v2

    ! Add routine to path
    CALL init_routine( routine_name)

    ncols           = grid%n      ! from
    nrows           = grid%n      ! to
    nnz_per_row_max = 2

    ! Initialise the matrix objects
    CALL MatCreate( PETSC_COMM_WORLD, M_ddx, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddy, perr)

    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( M_ddx, 'mpiaij', perr)
    CALL MatSetType( M_ddy, 'mpiaij', perr)

    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( M_ddx, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddy, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)

    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( M_ddx, perr)
    CALL MatSetFromOptions( M_ddy, perr)

    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( M_ddx, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddy, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)

    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( M_ddx, istart, iend, perr)

    DO n = istart+1, iend ! +1 because PETSc indexes from 0

      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)

      ! d/dx
      IF     (i == 1) THEN
        n_ip1 = grid%ij2n( i+1,j)
        v1 = -1._dp / grid%dx
        v2 =  1._dp / grid%dx
        CALL MatSetValues( M_ddx, 1, n-1, 1, n    -1, v1, INSERT_VALUES, perr)
        CALL MatSetValues( M_ddx, 1, n-1, 1, n_ip1-1, v2, INSERT_VALUES, perr)
      ELSEIF (i == grid%nx) THEN
        n_im1 = grid%ij2n( i-1,j)
        v1 =  1._dp / grid%dx
        v2 = -1._dp / grid%dx
        CALL MatSetValues( M_ddx, 1, n-1, 1, n    -1, v1, INSERT_VALUES, perr)
        CALL MatSetValues( M_ddx, 1, n-1, 1, n_im1-1, v2, INSERT_VALUES, perr)
      ELSE
        n_im1 = grid%ij2n( i-1,j)
        n_ip1 = grid%ij2n( i+1,j)
        v1 = -0.5_dp / grid%dx
        v2 =  0.5_dp / grid%dx
        CALL MatSetValues( M_ddx, 1, n-1, 1, n_im1-1, v1, INSERT_VALUES, perr)
        CALL MatSetValues( M_ddx, 1, n-1, 1, n_ip1-1, v2, INSERT_VALUES, perr)
      END IF

      ! d/dy
      IF     (j == 1) THEN
        n_jp1 = grid%ij2n( i,j+1)
        v1 = -1._dp / grid%dx
        v2 =  1._dp / grid%dx
        CALL MatSetValues( M_ddy, 1, n-1, 1, n    -1, v1, INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, n-1, 1, n_jp1-1, v2, INSERT_VALUES, perr)
      ELSEIF (j == grid%ny) THEN
        n_jm1 = grid%ij2n( i,j-1)
        v1 =  1._dp / grid%dx
        v2 = -1._dp / grid%dx
        CALL MatSetValues( M_ddy, 1, n-1, 1, n    -1, v1, INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, n-1, 1, n_jm1-1, v2, INSERT_VALUES, perr)
      ELSE
        n_jm1 = grid%ij2n( i,j-1)
        n_jp1 = grid%ij2n( i,j+1)
        v1 = -0.5_dp / grid%dx
        v2 =  0.5_dp / grid%dx
        CALL MatSetValues( M_ddy, 1, n-1, 1, n_jm1-1, v1, INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, n-1, 1, n_jp1-1, v2, INSERT_VALUES, perr)
      END IF

    END DO ! DO n = n1, n2
    CALL sync

    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.

    CALL MatAssemblyBegin( M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy, MAT_FINAL_ASSEMBLY, perr)

    CALL MatAssemblyEnd(   M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy, MAT_FINAL_ASSEMBLY, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_matrix_operators_grid

! == Directly apply a Neumann boundary condition to a data field
  SUBROUTINE apply_Neumann_BC_direct_a_2D( mesh, d_a)
    ! Directly apply a Neumann boundary condition to a data field

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_Neumann_BC_direct_a_2D'
    INTEGER                                            :: vi, vvi, vj
    REAL(dp)                                           :: sumd, sumw

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV) THEN
      CALL crash('data field is the wrong size!')
    END IF

    ! First the borders - average over adjacent non-boundary vertices
    DO vi = mesh%vi1, mesh%vi2

      IF (mesh%edge_index( vi) == 0 .OR. &
          mesh%edge_index( vi) == 2 .OR. &
          mesh%edge_index( vi) == 4 .OR. &
          mesh%edge_index( vi) == 6 .OR. &
          mesh%edge_index( vi) == 8) CYCLE

      sumd = 0._dp
      sumw = 0._dp

      DO vvi = 1, mesh%nC( vi)
        vj = mesh%C( vi,vvi)
        IF (mesh%edge_index( vj) > 0) CYCLE
        sumd = sumd + d_a( vj)
        sumw = sumw + 1._dp
      END DO

      IF (sumw > 0._dp) THEN
        d_a( vi) = sumd / sumw
      ELSE
        ! This border vertex has no non-border neighbours, which shouldn't be possible!
        CALL crash('border vertex {int_01} has no non-border neighbours!', int_01 = vi)
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Then the corners - average over all adjacent vertices
    IF (par%master) THEN
      DO vi = 1, 4

        sumd = 0._dp
        sumw = 0._dp

        DO vvi = 1, mesh%nC( vi)
          vj = mesh%C( vi,vvi)
          sumd = sumd + d_a( vj)
          sumw = sumw + 1._dp
        END DO

        d_a( vi) = sumd / sumw

      END DO ! DO vi = mesh%vi1, mesh%vi2
    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_Neumann_BC_direct_a_2D

  SUBROUTINE apply_Neumann_BC_direct_a_3D( mesh, d_a)
    ! Directly apply a Neumann boundary condition to a data field

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_Neumann_BC_direct_a_3D'
    INTEGER                                            :: vi, vvi, vj
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: sumd
    REAL(dp)                                           :: sumw

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV) THEN
      CALL crash('data field is the wrong size!')
    END IF

    ALLOCATE( sumd( SIZE( d_a,2)))

    ! First the borders - average over adjacent non-boundary vertices
    DO vi = mesh%vi1, mesh%vi2

      IF (mesh%edge_index( vi) == 0 .OR. &
          mesh%edge_index( vi) == 2 .OR. &
          mesh%edge_index( vi) == 4 .OR. &
          mesh%edge_index( vi) == 6 .OR. &
          mesh%edge_index( vi) == 8) CYCLE

      sumd = 0._dp
      sumw = 0._dp

      DO vvi = 1, mesh%nC( vi)
        vj = mesh%C( vi,vvi)
        IF (mesh%edge_index( vj) > 0) CYCLE
        sumd = sumd + d_a( vj,:)
        sumw = sumw + 1._dp
      END DO

      IF (sumw > 0._dp) THEN
        d_a( vi,:) = sumd / sumw
      ELSE
        ! This border vertex has no non-border neighbours, which shouldn't be possible!
        CALL crash('border vertex {int_01} has no non-border neighbours!', int_01 = vi)
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Then the corners - average over all adjacent vertices
    IF (par%master) THEN
      DO vi = 1, 4

        sumd = 0._dp
        sumw = 0._dp

        DO vvi = 1, mesh%nC( vi)
          vj = mesh%C( vi,vvi)
          sumd = sumd + d_a( vj,:)
          sumw = sumw + 1._dp
        END DO

        d_a( vi,:) = sumd / sumw

      END DO ! DO vi = mesh%vi1, mesh%vi2
    END IF ! IF (par%master) THEN
    CALL sync

    DEALLOCATE( sumd)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_Neumann_BC_direct_a_3D

  SUBROUTINE apply_Neumann_BC_direct_b_2D( mesh, d_b)
    ! Directly apply a Neumann boundary condition to a data field

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_Neumann_BC_direct_b_2D'
    INTEGER                                            :: ti, tti, tj
    REAL(dp)                                           :: sumd, sumw

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri) THEN
      CALL crash('data field is the wrong size!')
    END IF

    DO ti = mesh%ti1, mesh%ti2

      IF (mesh%Tri_edge_index( ti) == 0) CYCLE

      sumd = 0._dp
      sumw = 0._dp

      DO tti = 1, 3
        tj = mesh%TriC( ti,tti)
        IF (tj == 0) CYCLE
        IF (mesh%Tri_edge_index( tj) > 0) CYCLE
        sumd = sumd + d_b( tj)
        sumw = sumw + 1._dp
      END DO

      IF (sumw > 0._dp) THEN
        d_b( ti) = sumd / sumw
      ELSE
        ! This border triangle has no non-border neighbours, which shouldn't be possible!
        CALL crash('border triangle {int_01} has no non-border neighbours!', int_01 = ti)
      END IF

    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_Neumann_BC_direct_b_2D

  SUBROUTINE apply_Neumann_BC_direct_b_3D( mesh, d_b)
    ! Directly apply a Neumann boundary condition to a data field

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_Neumann_BC_direct_b_3D'
    INTEGER                                            :: ti, tti, tj
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: sumd
    REAL(dp)                                           :: sumw

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri) THEN
      CALL crash('data field is the wrong size!')
    END IF

    ALLOCATE( sumd( SIZE( d_b,2)))

    DO ti = mesh%ti1, mesh%ti2

      IF (mesh%Tri_edge_index( ti) == 0) CYCLE

      sumd = 0._dp
      sumw = 0._dp

      DO tti = 1, 3
        tj = mesh%TriC( ti,tti)
        IF (tj == 0) CYCLE
        IF (mesh%Tri_edge_index( tj) > 0) CYCLE
        sumd = sumd + d_b( tj,:)
        sumw = sumw + 1._dp
      END DO

      IF (sumw > 0._dp) THEN
        d_b( ti,:) = sumd / sumw
      ELSE
        ! This border triangle has no non-border neighbours, which shouldn't be possible!
        CALL crash('border triangle {int_01} has no non-border neighbours!', int_01 = ti)
      END IF

    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync

    DEALLOCATE( sumd)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_Neumann_BC_direct_b_3D

END MODULE mesh_operators_module
