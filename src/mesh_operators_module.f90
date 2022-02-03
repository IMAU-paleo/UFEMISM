MODULE mesh_operators_module

  ! Routines for calculating and applying mapping and gradient operators on the mesh.

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE petscksp
  USE configuration_module,            ONLY: dp, C
  USE parameters_module
  USE petsc_module,                    ONLY: perr
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
  
  ! Import specific functionality
  USE parallel_module,                 ONLY: adapt_shared_int_1D,    adapt_shared_dp_1D, &
                                             adapt_shared_int_2D,    adapt_shared_dp_2D, &
                                             adapt_shared_int_3D,    adapt_shared_dp_3D, &
                                             adapt_shared_bool_1D
  USE data_types_module,               ONLY: type_mesh, type_grid
  USE utilities_module,                ONLY: calc_matrix_inverse_2_by_2, calc_matrix_inverse_3_by_3, calc_matrix_inverse_general 
  USE petsc_module,                    ONLY: multiply_PETSc_matrix_with_vector_1D, multiply_PETSc_matrix_with_vector_2D, mat_petsc2CSR

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
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'map_a_to_b_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_a_b, d_a, d_b)
    
  END SUBROUTINE map_a_to_b_2D
  SUBROUTINE map_a_to_c_2D( mesh, d_a, d_c)
    ! Map a 2-D data field from the a (vertex) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_c
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'map_a_to_c_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_a_c, d_a, d_c)
    
  END SUBROUTINE map_a_to_c_2D
  SUBROUTINE map_b_to_a_2D( mesh, d_b, d_a)
    ! Map a 2-D data field from the b (triangle) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_a
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'map_b_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_b_a, d_b, d_a)
    
  END SUBROUTINE map_b_to_a_2D
  SUBROUTINE map_b_to_c_2D( mesh, d_b, d_c)
    ! Map a 2-D data field from the b (triangle) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_c
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'map_b_to_c_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_b_c, d_b, d_c)
    
  END SUBROUTINE map_b_to_c_2D
  SUBROUTINE map_c_to_a_2D( mesh, d_c, d_a)
    ! Map a 2-D data field from the c (edge) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_a
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'map_c_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_c_a, d_c, d_a)
    
  END SUBROUTINE map_c_to_a_2D
  SUBROUTINE map_c_to_b_2D( mesh, d_c, d_b)
    ! Map a 2-D data field from the c (edge) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_b
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'map_c_to_b_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_map_c_b, d_c, d_b)
    
  END SUBROUTINE map_c_to_b_2D
  
  ! 3-D
  SUBROUTINE map_a_to_b_3D( mesh, d_a, d_b)
    ! Map a 3-D data field from the a (vertex) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_a,2) /= SIZE( d_b,2)) THEN
      IF (par%master) WRITE(0,*) 'map_a_to_b_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_a_b, d_a, d_b)
    
  END SUBROUTINE map_a_to_b_3D
  SUBROUTINE map_a_to_c_3D( mesh, d_a, d_c)
    ! Map a 3-D data field from the a (vertex) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_c
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_a,2) /= SIZE( d_c,2)) THEN
      IF (par%master) WRITE(0,*) 'map_a_to_c_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_a_c, d_a, d_c)
    
  END SUBROUTINE map_a_to_c_3D
  SUBROUTINE map_b_to_a_3D( mesh, d_b, d_a)
    ! Map a 3-D data field from the b (triangle) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_a
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_b,2) /= SIZE( d_a,2)) THEN
      IF (par%master) WRITE(0,*) 'map_b_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_b_a, d_b, d_a)
    
  END SUBROUTINE map_b_to_a_3D
  SUBROUTINE map_b_to_c_3D( mesh, d_b, d_c)
    ! Map a 3-D data field from the b (triangle) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_c
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_b,2) /= SIZE( d_c,2)) THEN
      IF (par%master) WRITE(0,*) 'map_b_to_c_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_b_c, d_b, d_c)
    
  END SUBROUTINE map_b_to_c_3D
  SUBROUTINE map_c_to_a_3D( mesh, d_c, d_a)
    ! Map a 3-D data field from the c (edge) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_a
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_c,2) /= SIZE( d_a,2)) THEN
      IF (par%master) WRITE(0,*) 'map_c_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_c_a, d_c, d_a)
    
  END SUBROUTINE map_c_to_a_3D
  SUBROUTINE map_c_to_b_3D( mesh, d_c, d_b)
    ! Map a 3-D data field from the c (edge) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_b
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_c,2) /= SIZE( d_b,2)) THEN
      IF (par%master) WRITE(0,*) 'map_c_to_b_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_map_c_b, d_c, d_b)
    
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
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddx_a_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_a_a, d_a, ddx_a)
    
  END SUBROUTINE ddx_a_to_a_2D
  SUBROUTINE ddx_a_to_b_2D( mesh, d_a, ddx_b)
    ! ddx a 2-D data field from the a (vertex) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddx_a_to_b_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_a_b, d_a, ddx_b)
    
  END SUBROUTINE ddx_a_to_b_2D
  SUBROUTINE ddx_a_to_c_2D( mesh, d_a, ddx_c)
    ! ddx a 2-D data field from the a (vertex) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_c
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_a_to_c_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_a_c, d_a, ddx_c)
    
  END SUBROUTINE ddx_a_to_c_2D
  SUBROUTINE ddx_b_to_a_2D( mesh, d_b, ddx_a)
    ! ddx a 2-D data field from the b (triangle) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_a
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddx_b_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_b_a, d_b, ddx_a)
    
  END SUBROUTINE ddx_b_to_a_2D
  SUBROUTINE ddx_b_to_b_2D( mesh, d_b, ddx_b)
    ! ddx a 2-D data field from the b (triangle) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_b
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddx_b_to_b_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_b_b, d_b, ddx_b)
    
  END SUBROUTINE ddx_b_to_b_2D
  SUBROUTINE ddx_b_to_c_2D( mesh, d_b, ddx_c)
    ! ddx a 2-D data field from the b (triangle) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_c
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_b_to_c_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_b_c, d_b, ddx_c)
    
  END SUBROUTINE ddx_b_to_c_2D
  SUBROUTINE ddx_c_to_a_2D( mesh, d_c, ddx_a)
    ! ddx a 2-D data field from the c (edge) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_a
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddx_c_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_c_a, d_c, ddx_a)
    
  END SUBROUTINE ddx_c_to_a_2D
  SUBROUTINE ddx_c_to_b_2D( mesh, d_c, ddx_b)
    ! ddx a 2-D data field from the c (edge) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_b
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddx_c_to_b_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_c_b, d_c, ddx_b)
    
  END SUBROUTINE ddx_c_to_b_2D
  SUBROUTINE ddx_c_to_c_2D( mesh, d_c, ddx_c)
    ! ddx a 2-D data field from the c (edge) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_c
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_c_to_c_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_c_c, d_c, ddx_c)
    
  END SUBROUTINE ddx_c_to_c_2D

  ! 3-D
  SUBROUTINE ddx_a_to_a_3D( mesh, d_a, ddx_a)
    ! ddx a 3-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_a
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddx_a_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_a_a, d_a, ddx_a)
    
  END SUBROUTINE ddx_a_to_a_3D
  SUBROUTINE ddx_a_to_b_3D( mesh, d_a, ddx_b)
    ! ddx a 3-D data field from the a (vertex) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddx_a_to_b_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_a_b, d_a, ddx_b)
    
  END SUBROUTINE ddx_a_to_b_3D
  SUBROUTINE ddx_a_to_c_3D( mesh, d_a, ddx_c)
    ! ddx a 3-D data field from the a (vertex) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_c
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_a_to_c_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_a_c, d_a, ddx_c)
    
  END SUBROUTINE ddx_a_to_c_3D
  SUBROUTINE ddx_b_to_a_3D( mesh, d_b, ddx_a)
    ! ddx a 3-D data field from the b (triangle) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_a
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddx_b_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_b_a, d_b, ddx_a)
    
  END SUBROUTINE ddx_b_to_a_3D
  SUBROUTINE ddx_b_to_b_3D( mesh, d_b, ddx_b)
    ! ddx a 3-D data field from the b (triangle) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_b
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddx_b_to_b_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_b_b, d_b, ddx_b)
    
  END SUBROUTINE ddx_b_to_b_3D
  SUBROUTINE ddx_b_to_c_3D( mesh, d_b, ddx_c)
    ! ddx a 3-D data field from the b (triangle) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_c
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_b_to_c_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_b_c, d_b, ddx_c)
    
  END SUBROUTINE ddx_b_to_c_3D
  SUBROUTINE ddx_c_to_a_3D( mesh, d_c, ddx_a)
    ! ddx a 3-D data field from the c (edge) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_a
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddx_c_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_c_a, d_c, ddx_a)
    
  END SUBROUTINE ddx_c_to_a_3D
  SUBROUTINE ddx_c_to_b_3D( mesh, d_c, ddx_b)
    ! ddx a 3-D data field from the c (edge) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_b
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddx_c_to_b_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_c_b, d_c, ddx_b)
    
  END SUBROUTINE ddx_c_to_b_3D
  SUBROUTINE ddx_c_to_c_3D( mesh, d_c, ddx_c)
    ! ddx a 3-D data field from the c (edge) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_c
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddx_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_c_to_c_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_c_c, d_c, ddx_c)
    
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
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddy_a_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_a_a, d_a, ddy_a)
    
  END SUBROUTINE ddy_a_to_a_2D
  SUBROUTINE ddy_a_to_b_2D( mesh, d_a, ddy_b)
    ! ddy a 2-D data field from the a (vertex) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddy_a_to_b_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_a_b, d_a, ddy_b)
    
  END SUBROUTINE ddy_a_to_b_2D
  SUBROUTINE ddy_a_to_c_2D( mesh, d_a, ddy_c)
    ! ddy a 2-D data field from the a (vertex) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_c
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_a_to_c_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_a_c, d_a, ddy_c)
    
  END SUBROUTINE ddy_a_to_c_2D
  SUBROUTINE ddy_b_to_a_2D( mesh, d_b, ddy_a)
    ! ddy a 2-D data field from the b (triangle) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_a
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddy_b_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_b_a, d_b, ddy_a)
    
  END SUBROUTINE ddy_b_to_a_2D
  SUBROUTINE ddy_b_to_b_2D( mesh, d_b, ddy_b)
    ! ddy a 2-D data field from the b (triangle) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_b
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddy_b_to_b_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_b_b, d_b, ddy_b)
    
  END SUBROUTINE ddy_b_to_b_2D
  SUBROUTINE ddy_b_to_c_2D( mesh, d_b, ddy_c)
    ! ddy a 2-D data field from the b (triangle) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_c
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_b_to_c_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_b_c, d_b, ddy_c)
    
  END SUBROUTINE ddy_b_to_c_2D
  SUBROUTINE ddy_c_to_a_2D( mesh, d_c, ddy_a)
    ! ddy a 2-D data field from the c (edge) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_a
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddy_c_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_c_a, d_c, ddy_a)
    
  END SUBROUTINE ddy_c_to_a_2D
  SUBROUTINE ddy_c_to_b_2D( mesh, d_c, ddy_b)
    ! ddy a 2-D data field from the c (edge) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_b
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddy_c_to_b_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_c_b, d_c, ddy_b)
    
  END SUBROUTINE ddy_c_to_b_2D
  SUBROUTINE ddy_c_to_c_2D( mesh, d_c, ddy_c)
    ! ddy a 2-D data field from the c (edge) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_c
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_c_to_c_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_c_c, d_c, ddy_c)
    
  END SUBROUTINE ddy_c_to_c_2D

  ! 3-D
  SUBROUTINE ddy_a_to_a_3D( mesh, d_a, ddy_a)
    ! ddy a 3-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_a
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddy_a_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_a_a, d_a, ddy_a)
    
  END SUBROUTINE ddy_a_to_a_3D
  SUBROUTINE ddy_a_to_b_3D( mesh, d_a, ddy_b)
    ! ddy a 3-D data field from the a (vertex) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddy_a_to_b_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_a_b, d_a, ddy_b)
    
  END SUBROUTINE ddy_a_to_b_3D
  SUBROUTINE ddy_a_to_c_3D( mesh, d_a, ddy_c)
    ! ddy a 3-D data field from the a (vertex) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_c
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_a_to_c_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_a_c, d_a, ddy_c)
    
  END SUBROUTINE ddy_a_to_c_3D
  SUBROUTINE ddy_b_to_a_3D( mesh, d_b, ddy_a)
    ! ddy a 3-D data field from the b (triangle) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_a
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddy_b_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_b_a, d_b, ddy_a)
    
  END SUBROUTINE ddy_b_to_a_3D
  SUBROUTINE ddy_b_to_b_3D( mesh, d_b, ddy_b)
    ! ddy a 3-D data field from the b (triangle) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_b
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddy_b_to_b_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_b_b, d_b, ddy_b)
    
  END SUBROUTINE ddy_b_to_b_3D
  SUBROUTINE ddy_b_to_c_3D( mesh, d_b, ddy_c)
    ! ddy a 3-D data field from the b (triangle) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_c
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_b_to_c_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_b_c, d_b, ddy_c)
    
  END SUBROUTINE ddy_b_to_c_3D
  SUBROUTINE ddy_c_to_a_3D( mesh, d_c, ddy_a)
    ! ddy a 3-D data field from the c (edge) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_a
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'ddy_c_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_c_a, d_c, ddy_a)
    
  END SUBROUTINE ddy_c_to_a_3D
  SUBROUTINE ddy_c_to_b_3D( mesh, d_c, ddy_b)
    ! ddy a 3-D data field from the c (edge) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_b
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'ddy_c_to_b_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_c_b, d_c, ddy_b)
    
  END SUBROUTINE ddy_c_to_b_3D
  SUBROUTINE ddy_c_to_c_3D( mesh, d_c, ddy_c)
    ! ddy a 3-D data field from the c (edge) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_c
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( ddy_c,1) /= mesh%nAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_c_to_c_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_c_c, d_c, ddy_c)
    
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
    REAL(dp), DIMENSION(:    ), POINTER                ::  ddx_b
    INTEGER                                            :: wddx_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dx2_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'd2dx2_a_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_b, wddx_b)
    
    ! ddx_b = M_ddx_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_a_b, d_a, ddx_b)
    
    ! d2dx2_a = M_ddx_b_a * ddx_b
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_b_a, ddx_b, d2dx2_a)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddx_b)
    
  END SUBROUTINE d2dx2_a_to_a_2D
  SUBROUTINE d2dxdy_a_to_a_2D( mesh, d_a, d2dxdy_a)
    ! d2/dxdy a 2-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dxdy_a
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER                ::  ddx_b
    INTEGER                                            :: wddx_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dxdy_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'd2dxdy_a_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_b, wddx_b)
    
    ! ddx_b = M_ddx_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddx_a_b, d_a, ddx_b)
    
    ! d2dxdy_a = M_ddy_b_a * ddx_b
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_b_a, ddx_b, d2dxdy_a)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddx_b)
    
  END SUBROUTINE d2dxdy_a_to_a_2D
  SUBROUTINE d2dy2_a_to_a_2D(  mesh, d_a, d2dy2_a)
    ! d2/dy2 a 2-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dy2_a
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER                ::  ddy_b
    INTEGER                                            :: wddy_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dy2_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'd2dy2_a_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_b, wddy_b)
    
    ! ddy_b = M_ddy_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_a_b, d_a, ddy_b)
    
    ! d2dy2_a = M_ddy_b_a * ddy_b
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddy_b_a, ddy_b, d2dy2_a)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddy_b)
    
  END SUBROUTINE d2dy2_a_to_a_2D

  SUBROUTINE d2dx2_a_to_a_3D(  mesh, d_a, d2dx2_a)
    ! d2/dx2 a 3-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dx2_a
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  ddx_b
    INTEGER                                            :: wddx_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dx2_a,1) /= mesh%nV .OR. SIZE( d_a,2) /= SIZE( d2dx2_a,2)) THEN
      IF (par%master) WRITE(0,*) 'd2dx2_a_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nTri, SIZE( d_a,2), ddx_b, wddx_b)
    
    ! ddx_b = M_ddx_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_a_b, d_a, ddx_b)
    
    ! d2dx2_a = M_ddx_b_a * ddx_b
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_b_a, ddx_b, d2dx2_a)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddx_b)
    
  END SUBROUTINE d2dx2_a_to_a_3D
  SUBROUTINE d2dxdy_a_to_a_3D( mesh, d_a, d2dxdy_a)
    ! d2/dxdy a 3-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dxdy_a
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  ddx_b
    INTEGER                                            :: wddx_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dxdy_a,1) /= mesh%nV .OR. SIZE( d_a,2) /= SIZE( d2dxdy_a,2)) THEN
      IF (par%master) WRITE(0,*) 'd2dxdy_a_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nTri, SIZE( d_a,2), ddx_b, wddx_b)
    
    ! ddx_b = M_ddx_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddx_a_b, d_a, ddx_b)
    
    ! d2dxdy_a = M_ddy_b_a * ddx_b
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_b_a, ddx_b, d2dxdy_a)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddx_b)
    
  END SUBROUTINE d2dxdy_a_to_a_3D
  SUBROUTINE d2dy2_a_to_a_3D(  mesh, d_a, d2dy2_a)
    ! d2/dy2 a 3-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dy2_a
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  ddy_b
    INTEGER                                            :: wddy_b
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dy2_a,1) /= mesh%nV .OR. SIZE( d_a,2) /= SIZE( d2dy2_a,2)) THEN
      IF (par%master) WRITE(0,*) 'd2dy2_a_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nTri, SIZE( d_a,2), ddy_b, wddy_b)
    
    ! ddy_b = M_ddy_a_b * d_a
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_a_b, d_a, ddy_b)
    
    ! d2dy2_a = M_ddy_b_a * ddy_b
    CALL multiply_PETSc_matrix_with_vector_2D( mesh%M_ddy_b_a, ddy_b, d2dy2_a)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddy_b)
    
  END SUBROUTINE d2dy2_a_to_a_3D
  
! ==  Combined buv-grid
  SUBROUTINE map_b_to_buv_2D( mesh, d_b, d_buv)
    ! Map a data field defined on the b-grid (triangles) to the vector b-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_buv
    
    ! Local variables:
    INTEGER                                            :: ti, tiu, tiv
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri .OR. SIZE( d_buv,1) /= 2*mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'map_b_to_buv_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    DO ti = mesh%ti1, mesh%ti2
      tiu = 2*ti - 1
      tiv = 2*ti
      d_buv( tiu:tiv) = d_b( ti)
    END DO
    CALL sync
    
  END SUBROUTINE map_b_to_buv_2D
  
! == Calculate the matrix operators for mapping and gradients

  SUBROUTINE calc_matrix_operators_mesh( mesh)
    ! Calculate all the matrix operators
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Calculate matrix operators from the a (vertex)   to the other grids
    CALL calc_matrix_operators_a_a( mesh,                  mesh%M_ddx_a_a,  mesh%M_ddy_a_a)
    CALL calc_matrix_operators_a_b( mesh, mesh%M_map_a_b,  mesh%M_ddx_a_b,  mesh%M_ddy_a_b)
    CALL calc_matrix_operators_a_c( mesh, mesh%M_map_a_c,  mesh%M_ddx_a_c,  mesh%M_ddy_a_c)
    
    ! Calculate matrix operators from the b (triangle) to the other grids
    CALL calc_matrix_operators_b_a( mesh, mesh%M_map_b_a,  mesh%M_ddx_b_a,  mesh%M_ddy_b_a)
    CALL calc_matrix_operators_b_b( mesh,                  mesh%M_ddx_b_b,  mesh%M_ddy_b_b)
    CALL calc_matrix_operators_b_c( mesh, mesh%M_map_b_c,  mesh%M_ddx_b_c,  mesh%M_ddy_b_c)
    
    ! Calculate matrix operators from the c (edge)     to the other grids
    CALL calc_matrix_operators_c_a( mesh, mesh%M_map_c_a,  mesh%M_ddx_c_a,  mesh%M_ddy_c_a)
    CALL calc_matrix_operators_c_b( mesh, mesh%M_map_c_b,  mesh%M_ddx_c_b,  mesh%M_ddy_c_b)
    CALL calc_matrix_operators_c_c( mesh,                  mesh%M_ddx_c_c,  mesh%M_ddy_c_c)
    
    ! Calculate 2nd-order accurate matrix operators on the b-grid
    CALL calc_matrix_operators_2nd_order_b_b( mesh, &
      mesh%M2_ddx_b_b, &
      mesh%M2_ddy_b_b, &
      mesh%M2_d2dx2_b_b, &
      mesh%M2_d2dxdy_b_b, &
      mesh%M2_d2dy2_b_b)
    
    ! Calculate matrix operator for applying Neumann boundary conditions on border triangles
    CALL calc_matrix_operator_Neumann_BC_b( mesh)
    
    ! Set up the matrix operators in CSR-format (since I can't get the velocity solver working yet using only PETSc functionality)
    CALL mat_petsc2CSR( mesh%M2_ddx_b_b    , mesh%M2_ddx_b_b_CSR    )
    CALL mat_petsc2CSR( mesh%M2_ddy_b_b    , mesh%M2_ddy_b_b_CSR    )
    CALL mat_petsc2CSR( mesh%M2_d2dx2_b_b  , mesh%M2_d2dx2_b_b_CSR  )
    CALL mat_petsc2CSR( mesh%M2_d2dxdy_b_b , mesh%M2_d2dxdy_b_b_CSR )
    CALL mat_petsc2CSR( mesh%M2_d2dy2_b_b  , mesh%M2_d2dy2_b_b_CSR  )
    CALL mat_petsc2CSR( mesh%M_Neumann_BC_b, mesh%M_Neumann_BC_b_CSR)
    
  END SUBROUTINE calc_matrix_operators_mesh
  
  SUBROUTINE calc_matrix_operators_a_a( mesh,        M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(tMat),                          INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: vi, n, ci, vj, i
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc
    
    ! Matrix size
    ncols           = mesh%nV      ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_max = mesh%nC_mem+1
    
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
    
    ! Shape functions for a single node
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    DO vi = istart+1, iend ! +1 because PETSc indexes from 0
      
      ! Source points: the neighbouring vertices
      n = 0
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        n = n+1
        i_c( n) = vj
        x_c( n) = mesh%V( vj,1)
        y_c( n) = mesh%V( vj,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 2) CYCLE
      
      ! Destination point: the vertex
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_reg( x, y, n, x_c, y_c, Nfxi, Nfyi, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      CALL MatSetValues( M_ddx, 1, vi-1, 1, vi-1, Nfxi, INSERT_VALUES, perr)
      CALL MatSetValues( M_ddy, 1, vi-1, 1, vi-1, Nfyi, INSERT_VALUES, perr)
      DO i = 1, n
        CALL MatSetValues( M_ddx, 1, vi-1, 1, i_c( i)-1, Nfxc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, vi-1, 1, i_c( i)-1, Nfyc( i), INSERT_VALUES, perr)
      END DO
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    CALL MatAssemblyEnd(   M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_a_a
  SUBROUTINE calc_matrix_operators_a_b( mesh, M_map, M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the a (vertex) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(tMat),                          INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ti, n, vii, vi, i
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nV      ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 3
    
    ! Initialise the matrix objects
    CALL MatCreate( PETSC_COMM_WORLD, M_map, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddx, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddy, perr)
    
    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( M_map, 'mpiaij', perr)
    CALL MatSetType( M_ddx, 'mpiaij', perr)
    CALL MatSetType( M_ddy, 'mpiaij', perr)
    
    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( M_map, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddx, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddy, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    
    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( M_map, perr)
    CALL MatSetFromOptions( M_ddx, perr)
    CALL MatSetFromOptions( M_ddy, perr)
    
    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( M_map, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddx, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddy, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    
    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( M_ddx, istart, iend, perr)
    
    ! Shape functions for a single node
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    DO ti = istart+1, iend ! +1 because PETSc indexes from 0
      
      ! Source points: the neighbouring vertices
      n = 0
      DO vii = 1, 3
        vi = mesh%Tri( ti, vii)
        n = n+1
        i_c( n) = vi
        x_c( n) = mesh%V( vi,1)
        y_c( n) = mesh%V( vi,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) CYCLE
      
      ! Destination point: the triangle's geometric centre
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        CALL MatSetValues( M_map, 1, ti-1, 1, i_c( i)-1, Nfc(  i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddx, 1, ti-1, 1, i_c( i)-1, Nfxc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, ti-1, 1, i_c( i)-1, Nfyc( i), INSERT_VALUES, perr)
      END DO
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    CALL MatAssemblyEnd(   M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfc )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_a_b
  SUBROUTINE calc_matrix_operators_a_c( mesh, M_map, M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the a (vertex) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(tMat),                          INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: aci, n, acii, vi, i
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nV      ! from
    nrows           = mesh%nAc     ! to
    nnz_per_row_max = 4
    
    ! Initialise the matrix objects
    CALL MatCreate( PETSC_COMM_WORLD, M_map, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddx, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddy, perr)
    
    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( M_map, 'mpiaij', perr)
    CALL MatSetType( M_ddx, 'mpiaij', perr)
    CALL MatSetType( M_ddy, 'mpiaij', perr)
    
    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( M_map, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddx, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddy, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    
    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( M_map, perr)
    CALL MatSetFromOptions( M_ddx, perr)
    CALL MatSetFromOptions( M_ddy, perr)
    
    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( M_map, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddx, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddy, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    
    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( M_ddx, istart, iend, perr)
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    DO aci = istart+1, iend ! +1 because PETSc indexes from 0
      
      ! Source points: the three/four regular vertices surrounding the staggered vertex
      n = 0
      DO acii = 1, 4
        vi = mesh%Aci( aci, acii)
        IF (vi == 0) CYCLE
        n = n+1
        i_c( n) = vi
        x_c( n) = mesh%V( vi,1)
        y_c( n) = mesh%V( vi,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) CYCLE
      
      ! Destination point: the staggered vertex
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        CALL MatSetValues( M_map, 1, aci-1, 1, i_c( i)-1, Nfc(  i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddx, 1, aci-1, 1, i_c( i)-1, Nfxc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, aci-1, 1, i_c( i)-1, Nfyc( i), INSERT_VALUES, perr)
      END DO
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    CALL MatAssemblyEnd(   M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfc )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_a_c
  SUBROUTINE calc_matrix_operators_b_a( mesh, M_map, M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the b (triangle) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(tMat),                          INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: vi, n, iti, ti, i
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_max = mesh%nC_mem
    
    ! Initialise the matrix objects
    CALL MatCreate( PETSC_COMM_WORLD, M_map, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddx, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddy, perr)
    
    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( M_map, 'mpiaij', perr)
    CALL MatSetType( M_ddx, 'mpiaij', perr)
    CALL MatSetType( M_ddy, 'mpiaij', perr)
    
    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( M_map, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddx, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddy, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    
    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( M_map, perr)
    CALL MatSetFromOptions( M_ddx, perr)
    CALL MatSetFromOptions( M_ddy, perr)
    
    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( M_map, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddx, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddy, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    
    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( M_ddx, istart, iend, perr)
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    DO vi = istart+1, iend ! +1 because PETSc indexes from 0
      
      ! Source points: the geometric centres of the triangles surrounding the regular vertex
      n = 0
      DO iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi, iti)
        n = n+1
        i_c( n) = ti
        x_c( n) = mesh%TriGC( ti,1)
        y_c( n) = mesh%TriGC( ti,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) CYCLE
      
      ! Destination point: the regular vertex
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        CALL MatSetValues( M_map, 1, vi-1, 1, i_c( i)-1, Nfc(  i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddx, 1, vi-1, 1, i_c( i)-1, Nfxc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, vi-1, 1, i_c( i)-1, Nfyc( i), INSERT_VALUES, perr)
      END DO
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    CALL MatAssemblyEnd(   M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfc )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_b_a
  SUBROUTINE calc_matrix_operators_b_b( mesh,        M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the b (triangle) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(tMat),                          INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ti, n, tii, tj, i
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 4
    
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
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    DO ti = istart+1, iend ! +1 because PETSc indexes from 0
      
      ! Source points: the up to three adjacent triangles's geometric centres
      n = 0
      DO tii = 1, 3
        tj = mesh%TriC( ti,tii)
        IF (tj == 0) CYCLE
        n = n+1
        i_c( n) = tj
        x_c( n) = mesh%TriGC( tj,1)
        y_c( n) = mesh%TriGC( tj,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 2) CYCLE
      
      ! Destination point: the triangle's geometric centre
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_reg( x, y, n, x_c, y_c, Nfxi, Nfyi, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      CALL MatSetValues( M_ddx, 1, ti-1, 1, ti-1, Nfxi, INSERT_VALUES, perr)
      CALL MatSetValues( M_ddy, 1, ti-1, 1, ti-1, Nfyi, INSERT_VALUES, perr)
      DO i = 1, n
        CALL MatSetValues( M_ddx, 1, ti-1, 1, i_c( i)-1, Nfxc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, ti-1, 1, i_c( i)-1, Nfyc( i), INSERT_VALUES, perr)
      END DO
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    CALL MatAssemblyEnd(   M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_b_b
  SUBROUTINE calc_matrix_operators_b_c( mesh, M_map, M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the b (triangle) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(tMat),                          INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: aci, n, acii, vi, iti, ti, n_match, vii, aciii, i
    LOGICAL                                            :: is_listed
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nAc     ! to
    nnz_per_row_max = 6
    
    ! Initialise the matrix objects
    CALL MatCreate( PETSC_COMM_WORLD, M_map, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddx, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddy, perr)
    
    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( M_map, 'mpiaij', perr)
    CALL MatSetType( M_ddx, 'mpiaij', perr)
    CALL MatSetType( M_ddy, 'mpiaij', perr)
    
    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( M_map, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddx, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddy, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    
    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( M_map, perr)
    CALL MatSetFromOptions( M_ddx, perr)
    CALL MatSetFromOptions( M_ddy, perr)
    
    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( M_map, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddx, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddy, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    
    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( M_ddx, istart, iend, perr)
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    DO aci = istart+1, iend ! +1 because PETSc indexes from 0
      
      ! Source points: the two triangles adjacent to the edge, and the (up to) four other triangles adjoining them
      n = 0
      DO acii = 1, 4
        vi = mesh%Aci( aci,acii)
        IF (vi == 0) CYCLE
        
        DO iti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,iti)
          
          n_match = 0
          DO vii = 1, 3
            DO aciii = 1, 4
              IF (mesh%Tri( ti,vii) == mesh%Aci( aci,aciii)) THEN
                n_match = n_match + 1
              END IF
            END DO
          END DO
          
          IF (n_match >= 2) THEN
            is_listed = .FALSE.
            DO i = 1, n
              IF (i_c( i) == ti) THEN
                is_listed = .TRUE.
                EXIT
              END IF
            END DO
            IF (.NOT. is_listed) THEN
              n = n+1
              i_c( n) = ti
              x_c( n) = mesh%TriGC( ti,1)
              y_c( n) = mesh%TriGC( ti,2)
            END IF
          END IF
          
        END DO ! DO iti = 1, mesh.niTri( vi)
      END DO ! DO acii = 1, 4
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) CYCLE
      
      ! Destination point: the staggered vertex
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        CALL MatSetValues( M_map, 1, aci-1, 1, i_c( i)-1, Nfc(  i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddx, 1, aci-1, 1, i_c( i)-1, Nfxc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, aci-1, 1, i_c( i)-1, Nfyc( i), INSERT_VALUES, perr)
      END DO
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    CALL MatAssemblyEnd(   M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfc )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_b_c
  SUBROUTINE calc_matrix_operators_c_a( mesh, M_map, M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the c (edge) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(tMat),                          INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: vi, n, ci, aci, i
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nAc     ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_max = mesh%nC_mem
    
    ! Initialise the matrix objects
    CALL MatCreate( PETSC_COMM_WORLD, M_map, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddx, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddy, perr)
    
    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( M_map, 'mpiaij', perr)
    CALL MatSetType( M_ddx, 'mpiaij', perr)
    CALL MatSetType( M_ddy, 'mpiaij', perr)
    
    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( M_map, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddx, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddy, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    
    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( M_map, perr)
    CALL MatSetFromOptions( M_ddx, perr)
    CALL MatSetFromOptions( M_ddy, perr)
    
    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( M_map, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddx, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddy, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    
    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( M_ddx, istart, iend, perr)
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    DO vi = istart+1, iend ! +1 because PETSc indexes from 0
      
      ! Source points: the staggered vertices surrounding the regular vertex
      n = 0
      DO ci = 1, mesh%nC( vi)
        aci = mesh%iAci( vi, ci)
        n = n+1
        i_c( n) = aci
        x_c( n) = mesh%VAc( aci,1)
        y_c( n) = mesh%VAc( aci,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) CYCLE
      
      ! Destination point: the regular vertex
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        CALL MatSetValues( M_map, 1, vi-1, 1, i_c( i)-1, Nfc(  i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddx, 1, vi-1, 1, i_c( i)-1, Nfxc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, vi-1, 1, i_c( i)-1, Nfyc( i), INSERT_VALUES, perr)
      END DO
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    CALL MatAssemblyEnd(   M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfc )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_c_a
  SUBROUTINE calc_matrix_operators_c_b( mesh, M_map, M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the c (edge) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(tMat),                          INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ti, n, via, vib, vic, ci, vj, acab, acbc, acca, i
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nAc     ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 3
    
    ! Initialise the matrix objects
    CALL MatCreate( PETSC_COMM_WORLD, M_map, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddx, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddy, perr)
    
    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( M_map, 'mpiaij', perr)
    CALL MatSetType( M_ddx, 'mpiaij', perr)
    CALL MatSetType( M_ddy, 'mpiaij', perr)
    
    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( M_map, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddx, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddy, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    
    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( M_map, perr)
    CALL MatSetFromOptions( M_ddx, perr)
    CALL MatSetFromOptions( M_ddy, perr)
    
    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( M_map, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddx, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddy, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    
    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( M_ddx, istart, iend, perr)
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    DO ti = istart+1, iend ! +1 because PETSc indexes from 0
      
      ! Source points: the three staggered vertices on the triangle
      n = 0
      
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      
      acab = 0
      DO ci = 1, mesh%nC( via)
        vj = mesh%C( via,ci)
        IF (vj == vib) THEN
          acab = mesh%iAci( via,ci)
          EXIT
        END IF
      END DO
      
      n = n+1
      i_c( n) = acab
      x_c( n) = mesh%VAc( acab,1)
      y_c( n) = mesh%VAc( acab,2)
      
      acbc = 0
      DO ci = 1, mesh%nC( vib)
        vj = mesh%C( vib,ci)
        IF (vj == vic) THEN
          acbc = mesh%iAci( vib,ci)
          EXIT
        END IF
      END DO
      
      n = n+1
      i_c( n) = acbc
      x_c( n) = mesh%VAc( acbc,1)
      y_c( n) = mesh%VAc( acbc,2)
      
      acca = 0
      DO ci = 1, mesh%nC( vic)
        vj = mesh%C( vic,ci)
        IF (vj == via) THEN
          acca = mesh%iAci( vic,ci)
          EXIT
        END IF
      END DO
      
      n = n+1
      i_c( n) = acca
      x_c( n) = mesh%VAc( acca,1)
      y_c( n) = mesh%VAc( acca,2)
      
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) CYCLE
      
      ! Destination point: the triangle's geometric centre
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        CALL MatSetValues( M_map, 1, ti-1, 1, i_c( i)-1, Nfc(  i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddx, 1, ti-1, 1, i_c( i)-1, Nfxc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, ti-1, 1, i_c( i)-1, Nfyc( i), INSERT_VALUES, perr)
      END DO
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    CALL MatAssemblyEnd(   M_map, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfc )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_c_b
  SUBROUTINE calc_matrix_operators_c_c( mesh,        M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the c (edge) to the c (edge) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(tMat),                          INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: aci, n, acii, vi, ci, vc, acj, i
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc

    ncols           = mesh%nAc     ! from
    nrows           = mesh%nAc     ! to
    nnz_per_row_max = 5
    
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
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    DO aci = istart+1, iend ! +1 because PETSc indexes from 0
      
      ! Source points: all the staggered vertices lying on the edges of the
      ! two triangles adjoining the staggered vertex
      n = 0
      DO acii = 1, 2
        vi = mesh%Aci( aci,acii)
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          acj = mesh%iAci( vi,ci)
          IF (vc == mesh%Aci( aci,3) .OR. vc == mesh%Aci( aci,4)) THEN
            n = n+1
            i_c( n) = acj
            x_c( n) = mesh%VAc( acj,1)
            y_c( n) = mesh%VAc( acj,2)
          END IF
        END DO
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 2) CYCLE
      
      ! Destination point: the staggered vertex
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_reg( x, y, n, x_c, y_c, Nfxi, Nfyi, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      CALL MatSetValues( M_ddx, 1, aci-1, 1, aci-1, Nfxi, INSERT_VALUES, perr)
      CALL MatSetValues( M_ddy, 1, aci-1, 1, aci-1, Nfyi, INSERT_VALUES, perr)
      DO i = 1, n
        CALL MatSetValues( M_ddx, 1, aci-1, 1, i_c( i)-1, Nfxc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy, 1, aci-1, 1, i_c( i)-1, Nfyc( i), INSERT_VALUES, perr)
      END DO
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    CALL MatAssemblyEnd(   M_ddx, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy, MAT_FINAL_ASSEMBLY, perr)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_c_c
  
  SUBROUTINE calc_matrix_operators_2nd_order_b_b( mesh, M_ddx, M_ddy, M_d2dx2, M_d2dxdy, M_d2dy2)
    ! Calculate the matrix operators representing the d/dx, d/dy, d2/dx2,
    ! d2/dxdy, and d2/dy2 operations from the b (triangle) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(tMat),                          INTENT(INOUT) :: M_ddx, M_ddy, M_d2dx2, M_d2dxdy, M_d2dy2
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ti, n, n1, n_cur, tj, tk, i, ii
    LOGICAL                                            :: is_listed
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi, Nfxxi, Nfxyi, Nfyyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc, Nfxxc, Nfxyc, Nfyyc

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 10
    
    ! Initialise the matrix objects
    CALL MatCreate( PETSC_COMM_WORLD, M_ddx   , perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_ddy   , perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_d2dx2 , perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_d2dxdy, perr)
    CALL MatCreate( PETSC_COMM_WORLD, M_d2dy2 , perr)
    
    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( M_ddx   , 'mpiaij', perr)
    CALL MatSetType( M_ddy   , 'mpiaij', perr)
    CALL MatSetType( M_d2dx2 , 'mpiaij', perr)
    CALL MatSetType( M_d2dxdy, 'mpiaij', perr)
    CALL MatSetType( M_d2dy2 , 'mpiaij', perr)
    
    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( M_ddx   , PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_ddy   , PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_d2dx2 , PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_d2dxdy, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    CALL MatSetSizes( M_d2dy2 , PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    
    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( M_ddx   , perr)
    CALL MatSetFromOptions( M_ddy   , perr)
    CALL MatSetFromOptions( M_d2dx2 , perr)
    CALL MatSetFromOptions( M_d2dxdy, perr)
    CALL MatSetFromOptions( M_d2dy2 , perr)
    
    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( M_ddx   , nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_ddy   , nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_d2dx2 , nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_d2dxdy, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    CALL MatMPIAIJSetPreallocation( M_d2dy2 , nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    
    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( M_ddx, istart, iend, perr)
    
    ALLOCATE( i_c(   nnz_per_row_max))
    ALLOCATE( x_c(   nnz_per_row_max))
    ALLOCATE( y_c(   nnz_per_row_max))
    ALLOCATE( Nfxc(  nnz_per_row_max))
    ALLOCATE( Nfyc(  nnz_per_row_max))
    ALLOCATE( Nfxxc( nnz_per_row_max))
    ALLOCATE( Nfxyc( nnz_per_row_max))
    ALLOCATE( Nfyyc( nnz_per_row_max))

    DO ti = istart+1, iend ! +1 because PETSc indexes from 0
      
      ! Source points: 2-star neighbourhood of triangles
      n = 0
      DO n1 = 1, 3
        tj = mesh%TriC( ti,n1)
        IF (tj == 0) CYCLE
        n = n+1
        i_c( n) = tj
        x_c( n) = mesh%TriGC( tj,1)
        y_c( n) = mesh%TriGC( tj,2)
      END DO
      n_cur = n
      DO i = 1, n_cur
        tj = i_c( i)
        DO n1 = 1, 3
          tk = mesh%TriC( tj,n1)
          IF (tk == 0 ) CYCLE
          IF (tk == ti) CYCLE
          is_listed = .FALSE.
          DO ii = 1, n
            IF (i_c( ii) == tk) THEN
              is_listed = .TRUE.
              EXIT
            END IF
          END DO
          IF (.NOT. is_listed) THEN
            n = n+1
            i_c( n) = tk
            x_c( n) = mesh%TriGC( tk,1)
            y_c( n) = mesh%TriGC( tk,2)
          END IF
        END DO
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 2) CYCLE
      
      ! Destination point: the triangle's geometric centre
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_reg_2nd_order( x, y, n, x_c, y_c, Nfxi, Nfyi, Nfxxi, Nfxyi, Nfyyi, Nfxc, Nfyc, Nfxxc, Nfxyc, Nfyyc)
      
      ! Fill into sparse matrices
      CALL MatSetValues( M_ddx   , 1, ti-1, 1, ti-1, Nfxi , INSERT_VALUES, perr)
      CALL MatSetValues( M_ddy   , 1, ti-1, 1, ti-1, Nfyi , INSERT_VALUES, perr)
      CALL MatSetValues( M_d2dx2 , 1, ti-1, 1, ti-1, Nfxxi, INSERT_VALUES, perr)
      CALL MatSetValues( M_d2dxdy, 1, ti-1, 1, ti-1, Nfxyi, INSERT_VALUES, perr)
      CALL MatSetValues( M_d2dy2 , 1, ti-1, 1, ti-1, Nfyyi, INSERT_VALUES, perr)
      DO i = 1, n
        CALL MatSetValues( M_ddx   , 1, ti-1, 1, i_c( i)-1, Nfxc(  i), INSERT_VALUES, perr)
        CALL MatSetValues( M_ddy   , 1, ti-1, 1, i_c( i)-1, Nfyc(  i), INSERT_VALUES, perr)
        CALL MatSetValues( M_d2dx2 , 1, ti-1, 1, i_c( i)-1, Nfxxc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_d2dxdy, 1, ti-1, 1, i_c( i)-1, Nfxyc( i), INSERT_VALUES, perr)
        CALL MatSetValues( M_d2dy2 , 1, ti-1, 1, i_c( i)-1, Nfyyc( i), INSERT_VALUES, perr)
      END DO
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( M_ddx   , MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_ddy   , MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_d2dx2 , MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_d2dxdy, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( M_d2dy2 , MAT_FINAL_ASSEMBLY, perr)
    
    CALL MatAssemblyEnd(   M_ddx   , MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_ddy   , MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_d2dx2 , MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_d2dxdy, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_d2dy2 , MAT_FINAL_ASSEMBLY, perr)
    CALL sync
    
    ! Clean up after yourself
    DEALLOCATE( i_c  )
    DEALLOCATE( x_c  )
    DEALLOCATE( y_c  )
    DEALLOCATE( Nfxc )
    DEALLOCATE( Nfyc )
    DEALLOCATE( Nfxxc)
    DEALLOCATE( Nfxyc)
    DEALLOCATE( Nfyyc)
    
  END SUBROUTINE calc_matrix_operators_2nd_order_b_b
  SUBROUTINE calc_matrix_operator_Neumann_BC_b( mesh)
    ! Calculate matrix operator for applying Neumann boundary conditions on border triangles
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER                                            :: ti, n, tti, tj

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 4
    
    ! Initialise the matrix objects
    CALL MatCreate( PETSC_COMM_WORLD, mesh%M_Neumann_BC_b, perr)
    
    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( mesh%M_Neumann_BC_b, 'mpiaij', perr)
    
    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( mesh%M_Neumann_BC_b, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)
    
    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( mesh%M_Neumann_BC_b, perr)
    
    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( mesh%M_Neumann_BC_b, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)
    
    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( mesh%M_Neumann_BC_b, istart, iend, perr)
    
    DO ti = istart+1, iend ! +1 because PETSc indexes from 0
      
      IF (mesh%Tri_edge_index( ti) > 0) THEN
        ! This triangle touchers the domain border; apply Neumann boundary conditions
        
        ! In this case, that means that the value of f on ti should be equal to the
        ! mean value of f on the neighbours of ti.
      
        ! Find number of neighbouring triangles
        n = 0
        DO tti = 1, 3
          tj = mesh%TriC( ti,tti)
          IF (tj == 0) CYCLE
          n = n+1
        END DO
        
        ! The triangle itself
        CALL MatSetValues( mesh%M_Neumann_BC_b, 1, ti-1, 1, ti-1, -1._dp, INSERT_VALUES, perr)
        
        ! Neighbouring triangles
        DO tti = 1, 3
          tj = mesh%TriC( ti,tti)
          IF (tj == 0) CYCLE
          CALL MatSetValues( mesh%M_Neumann_BC_b, 1, ti-1, 1, tj-1, 1._dp / REAL( n,dp), INSERT_VALUES, perr)
        END DO
        
      END IF ! IF (is_border_triangle) THEN
    
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.
    
    CALL MatAssemblyBegin( mesh%M_Neumann_BC_b, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   mesh%M_Neumann_BC_b, MAT_FINAL_ASSEMBLY, perr)
    
  END SUBROUTINE calc_matrix_operator_Neumann_BC_b
  
  SUBROUTINE calc_matrix_operators_grid( grid, M_ddx, M_ddy)
    ! Calculate matrix operators for partial derivatives on a regular grid (needed for conservative remapping)
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(tMat),                          INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER                                            :: n, i, j, n_ip1, n_im1, n_jp1, n_jm1
    REAL(dp)                                           :: v1, v2

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
    
  END SUBROUTINE calc_matrix_operators_grid
  
! == Routines for calculating neighbour functions for regular and staggered vertices
  SUBROUTINE calc_neighbour_functions_ls_reg( x, y, n, x_c, y_c, Nfxi, Nfyi, Nfxc, Nfyc)
    ! Calculate neighbour functions for regular vertex V at [x,y],
    ! surrounded by n regular vertices at [x_c, y_c]
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).
      
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x,y
    INTEGER,                             INTENT(IN)    :: n
    REAL(dp), DIMENSION(n    ),          INTENT(IN)    :: x_c, y_c
    REAL(dp),                            INTENT(INOUT) :: Nfxi, Nfyi
    REAL(dp), DIMENSION(n    ),          INTENT(INOUT) :: Nfxc, Nfyc
    
    ! Local variables:
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n    )                         :: dx, dy, w
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    REAL(dp), DIMENSION(2,2)                           :: ATWTWA, M
    
    ! Calculate distances relative to v_c
    DO ci = 1, n
      dx( ci) = x_c( ci) - x
      dy( ci) = y_c( ci) - y
    END DO
    
    ! Calculate the weights w
    DO ci = 1, n
      w( ci) = 1._dp / (NORM2( [dx( ci), dy( ci)])**q)
    END DO
    
    ! The matrix ATWTWA that needs to be inverted
    ATWTWA  = 0._dp
    DO ci = 1, n
      ATWTWA( 1,1) = ATWTWA( 1,1) + (w(ci)**2 * dx(ci)**2           )
      ATWTWA( 1,2) = ATWTWA( 1,2) + (w(ci)**2 * dx(ci)    *  dy(ci) )

      ATWTWA( 2,1) = ATWTWA( 1,2)
      ATWTWA( 2,2) = ATWTWA( 2,2) + (w(ci)**2 * dy(ci)**2           )
    END DO

    ! Invert ATWTWA to find M
    CALL calc_matrix_inverse_2_by_2( ATWTWA, M)

    ! Calculate neighbour functions    
    DO ci = 1, n
      Nfxc( ci) = M( 1,1) * w( ci)**2 * dx( ci) + M( 1,2) * w( ci)**2 * dy( ci)
      Nfyc( ci) = M( 2,1) * w( ci)**2 * dx( ci) + M( 2,2) * w( ci)**2 * dy( ci)
    END DO
    
    Nfxi = -SUM( Nfxc)
    Nfyi = -SUM( Nfyc)
    
  END SUBROUTINE calc_neighbour_functions_ls_reg
  SUBROUTINE calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
    ! Calculate neighbour functions for staggered vertex V at [x,y],
    ! surrounded by n regular vertices at [x_c, y_c]
    !
    ! Based on the least-squares approach from Syrakos et al. (2017), here
    ! adapted to the case of a staggered vertex (i.e. where f is not
    ! defined on the point where we want to know the derivatives).
      
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x,y
    INTEGER,                             INTENT(IN)    :: n
    REAL(dp), DIMENSION(n    ),          INTENT(IN)    :: x_c, y_c
    REAL(dp), DIMENSION(n    ),          INTENT(INOUT) :: Nfc, Nfxc, Nfyc
    
    ! Local variables:
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n    )                         :: dx, dy, w
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    REAL(dp), DIMENSION(3,3)                           :: ATWTWA, M
    
    ! Calculate distances relative to v_c
    DO ci = 1, n
      dx( ci) = x_c( ci) - x
      dy( ci) = y_c( ci) - y
    END DO
    
    ! Calculate the weights w
    DO ci = 1, n
      w( ci) = 1._dp / (NORM2( [dx( ci), dy( ci)])**q)
    END DO
    
    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n
      ATWTWA( 1,1) = ATWTWA( 1,1) + (w(ci)**2                       )
      ATWTWA( 1,2) = ATWTWA( 1,2) + (w(ci)**2 * dx(ci)              )
      ATWTWA( 1,3) = ATWTWA( 1,3) + (w(ci)**2 * dy(ci)              )
      
      ATWTWA( 2,1) = ATWTWA( 1,2)
      ATWTWA( 2,2) = ATWTWA( 2,2) + (w(ci)**2 * dx(ci)**2           )
      ATWTWA( 2,3) = ATWTWA( 2,3) + (w(ci)**2 * dx(ci)    *  dy(ci) )

      ATWTWA( 3,1) = ATWTWA( 1,3)
      ATWTWA( 3,2) = ATWTWA( 2,3)
      ATWTWA( 3,3) = ATWTWA( 3,3) + (w(ci)**2 * dy(ci)**2            )
    END DO

    ! Invert ATWTWA to find M
    CALL calc_matrix_inverse_3_by_3( ATWTWA, M)

    ! Calculate neighbour functions    
    DO ci = 1, n
      Nfc(  ci) = M( 1,1) * w( ci)**2 + M( 1,2) * w( ci)**2 * dx( ci) + M( 1,3) * w( ci)**2 * dy( ci)
      Nfxc( ci) = M( 2,1) * w( ci)**2 + M( 2,2) * w( ci)**2 * dx( ci) + M( 2,3) * w( ci)**2 * dy( ci)
      Nfyc( ci) = M( 3,1) * w( ci)**2 + M( 3,2) * w( ci)**2 * dx( ci) + M( 3,3) * w( ci)**2 * dy( ci)
    END DO
    
  END SUBROUTINE calc_neighbour_functions_ls_stag
  SUBROUTINE calc_neighbour_functions_ls_reg_2nd_order( x, y, n, x_c, y_c, Nfxi, Nfyi, Nfxxi, Nfxyi, Nfyyi, Nfxc, Nfyc, Nfxxc, Nfxyc, Nfyyc)
    ! Calculate neighbour functions for regular vertex V at [x,y],
    ! surrounded by n regular vertices at [x_c, y_c]
    !
    ! Based on the least-squares approach from Syrakos et al. (2017).
      
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x,y
    INTEGER,                             INTENT(IN)    :: n
    REAL(dp), DIMENSION(n    ),          INTENT(IN)    :: x_c, y_c
    REAL(dp),                            INTENT(INOUT) :: Nfxi, Nfyi, Nfxxi, Nfxyi, Nfyyi
    REAL(dp), DIMENSION(n    ),          INTENT(INOUT) :: Nfxc, Nfyc, Nfxxc, Nfxyc, Nfyyc
    
    ! Local variables:
    INTEGER                                            :: ci
    REAL(dp), DIMENSION(n    )                         :: dx, dy, w
    REAL(dp), PARAMETER                                :: q = 1.5_dp
    REAL(dp), DIMENSION(5,5)                           :: ATWTWA, M
    REAL(dp)                                           :: tx, ty, txx, txy, tyy
    
    ! Calculate distances relative to v_c
    DO ci = 1, n
      dx( ci) = x_c( ci) - x
      dy( ci) = y_c( ci) - y
    END DO
    
    ! Calculate the weights w
    DO ci = 1, n
      w( ci) = 1._dp / (NORM2( [dx( ci), dy( ci)])**q)
    END DO
    
    ! The matrix ATWTWA that needs to be inverted
    ATWTWA = 0._dp
    DO ci = 1, n
      ATWTWA( 1,1) = ATWTWA( 1,1) + (       w( ci)**2 * dx( ci)**2             )
      ATWTWA( 1,2) = ATWTWA( 1,2) + (       w( ci)**2 * dx( ci)    * dy( ci)   )
      ATWTWA( 1,3) = ATWTWA( 1,3) + (0.5  * w( ci)**2 * dx( ci)**3             )
      ATWTWA( 1,4) = ATWTWA( 1,4) + (       w( ci)**2 * dx( ci)**2 * dy( ci)   )
      ATWTWA( 1,5) = ATWTWA( 1,5) + (0.5  * w( ci)**2 * dx( ci)    * dy( ci)**2)

      ATWTWA( 2,1) = ATWTWA( 1,2)
      ATWTWA( 2,2) = ATWTWA( 2,2) + (       w( ci)**2              * dy( ci)**2)
      ATWTWA( 2,3) = ATWTWA( 2,3) + (0.5  * w( ci)**2 * dx( ci)**2 * dy( ci)   )
      ATWTWA( 2,4) = ATWTWA( 2,4) + (       w( ci)**2 * dx( ci)    * dy( ci)**2)
      ATWTWA( 2,5) = ATWTWA( 2,5) + (0.5  * w( ci)**2              * dy( ci)**3)

      ATWTWA( 3,1) = ATWTWA( 1,3)
      ATWTWA( 3,2) = ATWTWA( 2,3)
      ATWTWA( 3,3) = ATWTWA( 3,3) + (0.25 * w( ci)**2 * dx( ci)**4             )
      ATWTWA( 3,4) = ATWTWA( 3,4) + (0.5  * w( ci)**2 * dx( ci)**3 * dy( ci)   )
      ATWTWA( 3,5) = ATWTWA( 3,5) + (0.25 * w( ci)**2 * dx( ci)**2 * dy( ci)**2)

      ATWTWA( 4,1) = ATWTWA( 1,4)
      ATWTWA( 4,2) = ATWTWA( 2,4)
      ATWTWA( 4,3) = ATWTWA( 3,4)
      ATWTWA( 4,4) = ATWTWA( 4,4) + (       w( ci)**2 * dx( ci)**2 * dy( ci)**2)
      ATWTWA( 4,5) = ATWTWA( 4,5) + (0.5  * w( ci)**2 * dx( ci)    * dy( ci)**3)

      ATWTWA( 5,1) = ATWTWA( 1,5)
      ATWTWA( 5,2) = ATWTWA( 2,5)
      ATWTWA( 5,3) = ATWTWA( 3,5)
      ATWTWA( 5,4) = ATWTWA( 4,5)
      ATWTWA( 5,5) = ATWTWA( 5,5) + (0.25 * w( ci)**2              * dy( ci)**4)
    END DO

    ! Invert ATWTWA to find M
    CALL calc_matrix_inverse_general( ATWTWA, M)

    ! Calculate neighbour functions    
    DO ci = 1, n
      tx  = (M( 1,1) * dx( ci)) + (M( 1,2) * dy( ci)) + (0.5 * M( 1,3) * dx( ci)**2) + (M( 1,4) * dx( ci) * dy( ci)) + (0.5 * M( 1,5) * dy( ci)**2)
      ty  = (M( 2,1) * dx( ci)) + (M( 2,2) * dy( ci)) + (0.5 * M( 2,3) * dx( ci)**2) + (M( 2,4) * dx( ci) * dy( ci)) + (0.5 * M( 2,5) * dy( ci)**2)
      txx = (M( 3,1) * dx( ci)) + (M( 3,2) * dy( ci)) + (0.5 * M( 3,3) * dx( ci)**2) + (M( 3,4) * dx( ci) * dy( ci)) + (0.5 * M( 3,5) * dy( ci)**2)
      txy = (M( 4,1) * dx( ci)) + (M( 4,2) * dy( ci)) + (0.5 * M( 4,3) * dx( ci)**2) + (M( 4,4) * dx( ci) * dy( ci)) + (0.5 * M( 4,5) * dy( ci)**2)
      tyy = (M( 5,1) * dx( ci)) + (M( 5,2) * dy( ci)) + (0.5 * M( 5,3) * dx( ci)**2) + (M( 5,4) * dx( ci) * dy( ci)) + (0.5 * M( 5,5) * dy( ci)**2)

      Nfxc(  ci) = w( ci)**2 * tx
      Nfyc(  ci) = w( ci)**2 * ty
      Nfxxc( ci) = w( ci)**2 * txx
      Nfxyc( ci) = w( ci)**2 * txy
      Nfyyc( ci) = w( ci)**2 * tyy
    END DO
    
    Nfxi  = -SUM( Nfxc)
    Nfyi  = -SUM( Nfyc)
    Nfxxi = -SUM( Nfxxc)
    Nfxyi = -SUM( Nfxyc)
    Nfyyi = -SUM( Nfyyc)
    
  END SUBROUTINE calc_neighbour_functions_ls_reg_2nd_order
  
! == Directly apply a Neumann boundary condition to a data field
  SUBROUTINE apply_Neumann_BC_direct_a_2D( mesh, d_a)
    ! Directly apply a Neumann boundary condition to a data field
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_a
    
    ! Local variables:
    INTEGER                                            :: vi, vvi, vj
    REAL(dp)                                           :: sumd, sumw
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'apply_Neumann_BC_direct_a_2D - ERROR: data field is of the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
        WRITE(0,*) 'apply_Neumann_BC_direct_a_2D - ERROR: border vertex ', vi, ' has no non-border neighbours!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
    
  END SUBROUTINE apply_Neumann_BC_direct_a_2D
  SUBROUTINE apply_Neumann_BC_direct_a_3D( mesh, d_a)
    ! Directly apply a Neumann boundary condition to a data field
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d_a
    
    ! Local variables:
    INTEGER                                            :: vi, vvi, vj
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: sumd
    REAL(dp)                                           :: sumw
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'apply_Neumann_BC_direct_a_3D - ERROR: data field is of the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
        WRITE(0,*) 'apply_Neumann_BC_direct_a_3D - ERROR: border vertex ', vi, ' has no non-border neighbours!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
    
  END SUBROUTINE apply_Neumann_BC_direct_a_3D
  SUBROUTINE apply_Neumann_BC_direct_b_2D( mesh, d_b)
    ! Directly apply a Neumann boundary condition to a data field
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_b
    
    ! Local variables:
    INTEGER                                            :: ti, tti, tj
    REAL(dp)                                           :: sumd, sumw
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'apply_Neumann_BC_direct_b_2D - ERROR: data field is of the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
        WRITE(0,*) 'apply_Neumann_BC_direct_b_2D - ERROR: border triangle ', ti, ' has no non-border neighbours!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
  END SUBROUTINE apply_Neumann_BC_direct_b_2D
  SUBROUTINE apply_Neumann_BC_direct_b_3D( mesh, d_b)
    ! Directly apply a Neumann boundary condition to a data field
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d_b
    
    ! Local variables:
    INTEGER                                            :: ti, tti, tj
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: sumd
    REAL(dp)                                           :: sumw
    
    ! Safety
    IF (SIZE( d_b,1) /= mesh%nTri) THEN
      IF (par%master) WRITE(0,*) 'apply_Neumann_BC_direct_b_3D - ERROR: data field is of the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
        WRITE(0,*) 'apply_Neumann_BC_direct_b_3D - ERROR: border triangle ', ti, ' has no non-border neighbours!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    DEALLOCATE( sumd)
    
  END SUBROUTINE apply_Neumann_BC_direct_b_3D

END MODULE mesh_operators_module
