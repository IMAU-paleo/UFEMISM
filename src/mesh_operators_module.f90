MODULE mesh_operators_module

  ! Routines for calculating and applying mapping and gradient operators on the mesh.

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
  
  ! Import specific functionality
  USE parallel_module,                 ONLY: adapt_shared_int_1D,    adapt_shared_dp_1D, &
                                             adapt_shared_int_2D,    adapt_shared_dp_2D, &
                                             adapt_shared_int_3D,    adapt_shared_dp_3D, &
                                             adapt_shared_bool_1D
  USE data_types_module,               ONLY: type_mesh, type_sparse_matrix_CSR
  USE utilities_module,                ONLY: allocate_matrix_CSR_dist, extend_matrix_CSR_dist, finalise_matrix_CSR_dist, &
                                             deallocate_matrix_CSR, multiply_matrix_matrix_CSR, multiply_matrix_vector_CSR, &
                                             multiply_matrix_vector_2D_CSR, calc_matrix_inverse_2_by_2, calc_matrix_inverse_3_by_3

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
    CALL multiply_matrix_vector_CSR( mesh%M_map_a_b, d_a, d_b)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_map_a_c, d_a, d_c)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_map_b_a, d_b, d_a)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_map_b_c, d_b, d_c)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_map_c_a, d_c, d_a)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_map_c_b, d_c, d_b)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_a_b, d_a, d_b)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_a_c, d_a, d_c)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_b_a, d_b, d_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_b_c, d_b, d_c)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_c_a, d_c, d_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_c_b, d_c, d_b)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_a_a, d_a, ddx_a)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_a_b, d_a, ddx_b)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_a_c, d_a, ddx_c)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_b_a, d_b, ddx_a)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_b_b, d_b, ddx_b)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_b_c, d_b, ddx_c)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_c_a, d_c, ddx_a)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_c_b, d_c, ddx_b)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_c_c, d_c, ddx_c)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_a_a, d_a, ddx_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_a_b, d_a, ddx_b)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_a_c, d_a, ddx_c)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_b_a, d_b, ddx_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_b_b, d_b, ddx_b)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_b_c, d_b, ddx_c)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_c_a, d_c, ddx_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_c_b, d_c, ddx_b)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_c_c, d_c, ddx_c)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_a_a, d_a, ddy_a)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_a_b, d_a, ddy_b)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_a_c, d_a, ddy_c)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_b_a, d_b, ddy_a)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_b_b, d_b, ddy_b)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_b_c, d_b, ddy_c)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_c_a, d_c, ddy_a)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_c_b, d_c, ddy_b)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_c_c, d_c, ddy_c)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_a_a, d_a, ddy_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_a_b, d_a, ddy_b)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_a_c, d_a, ddy_c)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_b_a, d_b, ddy_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_b_b, d_b, ddy_b)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_b_c, d_b, ddy_c)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_c_a, d_c, ddy_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_c_b, d_c, ddy_b)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_c_c, d_c, ddy_c)
    
  END SUBROUTINE ddy_c_to_c_3D
  
! == d2/dx2, d2/dxdy, d2/dy2

  SUBROUTINE d2dx2_a_to_a_2D( mesh, d_a, d2dx2_a)
    ! d2/dx2 a 2-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dx2_a
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dx2_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'd2dx2_a_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_d2dx2_a_a, d_a, d2dx2_a)
    
  END SUBROUTINE d2dx2_a_to_a_2D
  SUBROUTINE d2dx2_a_to_a_3D( mesh, d_a, d2dx2_a)
    ! d2/dx2 a 2-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dx2_a
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dx2_a,1) /= mesh%nV .OR. SIZE( d_a,2) /= SIZE( d2dx2_a,2)) THEN
      IF (par%master) WRITE(0,*) 'd2dx2_a_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_d2dx2_a_a, d_a, d2dx2_a)
    
  END SUBROUTINE d2dx2_a_to_a_3D
  SUBROUTINE d2dxdy_a_to_a_2D( mesh, d_a, d2dxdy_a)
    ! d2/dxdy a 2-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dxdy_a
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dxdy_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'd2dxdy_a_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_d2dxdy_a_a, d_a, d2dxdy_a)
    
  END SUBROUTINE d2dxdy_a_to_a_2D
  SUBROUTINE d2dxdy_a_to_a_3D( mesh, d_a, d2dxdy_a)
    ! d2/dxdy a 2-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dxdy_a
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dxdy_a,1) /= mesh%nV .OR. SIZE( d_a,2) /= SIZE( d2dxdy_a,2)) THEN
      IF (par%master) WRITE(0,*) 'd2dxdy_a_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_d2dxdy_a_a, d_a, d2dxdy_a)
    
  END SUBROUTINE d2dxdy_a_to_a_3D
  SUBROUTINE d2dy2_a_to_a_2D( mesh, d_a, d2dy2_a)
    ! d2/dy2 a 2-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dy2_a
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dy2_a,1) /= mesh%nV) THEN
      IF (par%master) WRITE(0,*) 'd2dy2_a_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_d2dy2_a_a, d_a, d2dy2_a)
    
  END SUBROUTINE d2dy2_a_to_a_2D
  SUBROUTINE d2dy2_a_to_a_3D( mesh, d_a, d2dy2_a)
    ! d2/dy2 a 2-D data field from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dy2_a
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d2dy2_a,1) /= mesh%nV .OR. SIZE( d_a,2) /= SIZE( d2dy2_a,2)) THEN
      IF (par%master) WRITE(0,*) 'd2dy2_a_to_a_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_d2dy2_a_a, d_a, d2dy2_a)
    
  END SUBROUTINE d2dy2_a_to_a_3D
  
! == Operations between the a/c-grids and the ac-grid

  ! From the a/c-grids to the ac-grid
  SUBROUTINE move_a_to_aca_2D( mesh, d_a, d_aca)
    ! Move data from the a-grid to the a-part of the ac-grid
    !
    ! NOTE: sets c-part to zero!
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_aca
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_aca,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'move_a_to_aca_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_move_a_aca, d_a, d_aca)
    
  END SUBROUTINE move_a_to_aca_2D
  SUBROUTINE move_a_to_aca_3D( mesh, d_a, d_aca)
    ! Move data from the a-grid to the a-part of the ac-grid
    !
    ! NOTE: sets c-part to zero!
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_aca
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_aca,1) /= mesh%nVAaAc .OR. SIZE( d_a,2) /= SIZE( d_aca,2)) THEN
      IF (par%master) WRITE(0,*) 'move_a_to_aca_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_move_a_aca, d_a, d_aca)
    
  END SUBROUTINE move_a_to_aca_3D
  SUBROUTINE move_c_to_acc_2D( mesh, d_c, d_acc)
    ! Move data from the c-grid to the c-part of the ac-grid
    !
    ! NOTE: sets a-part to zero!
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_acc
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_acc,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'move_c_to_acc_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_move_c_acc, d_c, d_acc)
    
  END SUBROUTINE move_c_to_acc_2D
  SUBROUTINE move_c_to_acc_3D( mesh, d_c, d_acc)
    ! Move data from the c-grid to the c-part of the ac-grid
    !
    ! NOTE: sets a-part to zero!
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_acc
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_acc,1) /= mesh%nVAaAc .OR. SIZE( d_c,2) /= SIZE( d_acc,2)) THEN
      IF (par%master) WRITE(0,*) 'move_c_to_acc_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_move_c_acc, d_c, d_acc)
    
  END SUBROUTINE move_c_to_acc_3D
  SUBROUTINE move_a_and_c_to_ac_2D( mesh, d_a, d_c, d_ac)
    ! Move data from the a/c-grids to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_ac
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_aca,  d_acc
    INTEGER                                            :: wd_aca, wd_acc
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'move_a_and_c_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nVAaAc, d_aca, wd_aca)
    CALL allocate_shared_dp_1D( mesh%nVAaAc, d_acc, wd_acc)
    
    ! Move both parts
    CALL move_a_to_aca_2D( mesh, d_a, d_aca)
    CALL move_c_to_acc_2D( mesh, d_c, d_acc)
    
    ! Add them together
    d_ac( mesh%avi1:mesh%avi2) = d_aca( mesh%avi1:mesh%avi2) + &
                                 d_acc( mesh%avi1:mesh%avi2)
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_aca)
    CALL deallocate_shared( wd_acc)
    
  END SUBROUTINE move_a_and_c_to_ac_2D
  SUBROUTINE move_a_and_c_to_ac_3D( mesh, d_a, d_c, d_ac)
    ! Move data from the a/c-grids to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_ac
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_aca,  d_acc
    INTEGER                                            :: wd_aca, wd_acc
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_ac,1) /= mesh%nVAaAc &
        .OR. SIZE( d_a,2) /= SIZE( d_c,2) .OR. SIZE( d_a,2) /= SIZE( d_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'move_a_and_c_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nVAaAc, SIZE( d_a,2), d_aca, wd_aca)
    CALL allocate_shared_dp_2D( mesh%nVAaAc, SIZE( d_a,2), d_acc, wd_acc)
    
    ! Move both parts
    CALL move_a_to_aca_3D( mesh, d_a, d_aca)
    CALL move_c_to_acc_3D( mesh, d_c, d_acc)
    
    ! Add them together
    d_ac( mesh%avi1:mesh%avi2,:) = d_aca( mesh%avi1:mesh%avi2,:) + &
                                   d_acc( mesh%avi1:mesh%avi2,:)
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_aca)
    CALL deallocate_shared( wd_acc)
    
  END SUBROUTINE move_a_and_c_to_ac_3D
  
  SUBROUTINE map_a_to_ac_2D( mesh, d_a, d_ac)
    ! Map data from the a-grid to the ac-grid
    !
    ! NOTE: the a-part is moved, the c-part is mapped
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_ac
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_c
    INTEGER                                            :: wd_c
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'map_a_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nAc   , d_c  , wd_c  )
    
    ! Map data from the a-grid to the c-grid
    CALL map_a_to_c_2D( mesh, d_a, d_c)
    
    ! Move both parts
    CALL move_a_and_c_to_ac_2D( mesh, d_a, d_c, d_ac)
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_c)
    
  END SUBROUTINE map_a_to_ac_2D
  SUBROUTINE map_a_to_ac_3D( mesh, d_a, d_ac)
    ! Map data from the a-grid to the ac-grid
    !
    ! NOTE: the a-part is moved, the c-part is mapped
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_ac
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_c
    INTEGER                                            :: wd_c
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_a,2) /= SIZE( d_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'map_a_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nAc   , SIZE( d_a,2), d_c  , wd_c  )
    
    ! Map data from the a-grid to the c-grid
    CALL map_a_to_c_3D( mesh, d_a, d_c)
    
    ! Move both parts
    CALL move_a_and_c_to_ac_3D( mesh, d_a, d_c, d_ac)
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_c)
    
  END SUBROUTINE map_a_to_ac_3D
  SUBROUTINE map_c_to_ac_2D( mesh, d_c, d_ac)
    ! Map data from the a-grid to the ac-grid
    !
    ! NOTE: the c-part is moved, the a-part is mapped
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_ac
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_a
    INTEGER                                            :: wd_a
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'map_c_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV    , d_a  , wd_a  )
    
    ! Map data from the c-grid to the a-grid
    CALL map_c_to_a_2D( mesh, d_c, d_a)
    
    ! Move both parts
    CALL move_a_and_c_to_ac_2D( mesh, d_a, d_c, d_ac)
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_a)
    
  END SUBROUTINE map_c_to_ac_2D
  SUBROUTINE map_c_to_ac_3D( mesh, d_c, d_ac)
    ! Map data from the a-grid to the ac-grid
    !
    ! NOTE: the c-part is moved, the a-part is mapped
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_ac
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_a
    INTEGER                                            :: wd_a
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_c,2) /= SIZE( d_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'map_c_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV    , SIZE( d_c,2), d_a  , wd_a  )
    
    ! Map data from the a-grid to the c-grid
    CALL map_c_to_a_3D( mesh, d_c, d_a)
    
    ! Move both parts
    CALL move_a_and_c_to_ac_3D( mesh, d_a, d_c, d_ac)
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_a)
    
  END SUBROUTINE map_c_to_ac_3D
  
  SUBROUTINE ddx_a_to_ac_2D( mesh, d_a, ddx_ac)
    ! Calculate d/dx from the a-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_ac
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER                ::  ddx_a,  ddx_c,  ddx_aca,  ddx_acc
    INTEGER                                            :: wddx_a, wddx_c, wddx_aca, wddx_acc
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_a_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV    , ddx_a  , wddx_a  )
    CALL allocate_shared_dp_1D( mesh%nAc   , ddx_c  , wddx_c  )
    CALL allocate_shared_dp_1D( mesh%nVAaAc, ddx_aca, wddx_aca)
    CALL allocate_shared_dp_1D( mesh%nVAaAc, ddx_acc, wddx_acc)
    
    ! Calculate d/dx on the a/c-grids
    CALL ddx_a_to_a_2D( mesh, d_a, ddx_a)
    CALL ddx_a_to_c_2D( mesh, d_a, ddx_c)
    
    ! Move both parts
    CALL move_a_to_aca_2D( mesh, ddx_a, ddx_aca)
    CALL move_c_to_acc_2D( mesh, ddx_c, ddx_acc)
    
    ! Add them together
    ddx_ac( mesh%avi1:mesh%avi2) = ddx_aca( mesh%avi1:mesh%avi2) + &
                                   ddx_acc( mesh%avi1:mesh%avi2)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddx_c  )
    CALL deallocate_shared( wddx_aca)
    CALL deallocate_shared( wddx_acc)
    
  END SUBROUTINE ddx_a_to_ac_2D
  SUBROUTINE ddx_a_to_ac_3D( mesh, d_a, ddx_ac)
    ! Calculate d/dx from the a-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_ac
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  ddx_a,  ddx_c,  ddx_aca,  ddx_acc
    INTEGER                                            :: wddx_a, wddx_c, wddx_aca, wddx_acc
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddx_ac,1) /= mesh%nVAaAc .OR. SIZE( d_a,2) /= SIZE( ddx_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'ddx_a_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV    , SIZE( d_a,2), ddx_a  , wddx_a  )
    CALL allocate_shared_dp_2D( mesh%nAc   , SIZE( d_a,2), ddx_c  , wddx_c  )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, SIZE( d_a,2), ddx_aca, wddx_aca)
    CALL allocate_shared_dp_2D( mesh%nVAaAc, SIZE( d_a,2), ddx_acc, wddx_acc)
    
    ! Calculate d/dx on the a/c-grids
    CALL ddx_a_to_a_3D( mesh, d_a, ddx_a)
    CALL ddx_a_to_c_3D( mesh, d_a, ddx_c)
    
    ! Move both parts
    CALL move_a_to_aca_3D( mesh, ddx_a, ddx_aca)
    CALL move_c_to_acc_3D( mesh, ddx_c, ddx_acc)
    
    ! Add them together
    ddx_ac( mesh%avi1:mesh%avi2,:) = ddx_aca( mesh%avi1:mesh%avi2,:) + &
                                     ddx_acc( mesh%avi1:mesh%avi2,:)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddx_c  )
    CALL deallocate_shared( wddx_aca)
    CALL deallocate_shared( wddx_acc)
    
  END SUBROUTINE ddx_a_to_ac_3D
  
  SUBROUTINE ddy_a_to_ac_2D( mesh, d_a, ddy_ac)
    ! Calculate d/dy from the a-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_ac
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER                ::  ddy_a,  ddy_c,  ddy_aca,  ddy_acc
    INTEGER                                            :: wddy_a, wddy_c, wddy_aca, wddy_acc
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_a_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV    , ddy_a  , wddy_a  )
    CALL allocate_shared_dp_1D( mesh%nAc   , ddy_c  , wddy_c  )
    CALL allocate_shared_dp_1D( mesh%nVAaAc, ddy_aca, wddy_aca)
    CALL allocate_shared_dp_1D( mesh%nVAaAc, ddy_acc, wddy_acc)
    
    ! Calculate d/dy on the a/c-grids
    CALL ddy_a_to_a_2D( mesh, d_a, ddy_a)
    CALL ddy_a_to_c_2D( mesh, d_a, ddy_c)
    
    ! Move both parts
    CALL move_a_to_aca_2D( mesh, ddy_a, ddy_aca)
    CALL move_c_to_acc_2D( mesh, ddy_c, ddy_acc)
    
    ! Add them together
    ddy_ac( mesh%avi1:mesh%avi2) = ddy_aca( mesh%avi1:mesh%avi2) + &
                                   ddy_acc( mesh%avi1:mesh%avi2)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddy_c  )
    CALL deallocate_shared( wddy_aca)
    CALL deallocate_shared( wddy_acc)
    
  END SUBROUTINE ddy_a_to_ac_2D
  SUBROUTINE ddy_a_to_ac_3D( mesh, d_a, ddy_ac)
    ! Calculate d/dy from the a-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_ac
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  ddy_a,  ddy_c,  ddy_aca,  ddy_acc
    INTEGER                                            :: wddy_a, wddy_c, wddy_aca, wddy_acc
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( ddy_ac,1) /= mesh%nVAaAc .OR. SIZE( d_a,2) /= SIZE( ddy_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'ddy_a_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV    , SIZE( d_a,2), ddy_a  , wddy_a  )
    CALL allocate_shared_dp_2D( mesh%nAc   , SIZE( d_a,2), ddy_c  , wddy_c  )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, SIZE( d_a,2), ddy_aca, wddy_aca)
    CALL allocate_shared_dp_2D( mesh%nVAaAc, SIZE( d_a,2), ddy_acc, wddy_acc)
    
    ! Calculate d/dy on the a/c-grids
    CALL ddy_a_to_a_3D( mesh, d_a, ddy_a)
    CALL ddy_a_to_c_3D( mesh, d_a, ddy_c)
    
    ! Move both parts
    CALL move_a_to_aca_3D( mesh, ddy_a, ddy_aca)
    CALL move_c_to_acc_3D( mesh, ddy_c, ddy_acc)
    
    ! Add them together
    ddy_ac( mesh%avi1:mesh%avi2,:) = ddy_aca( mesh%avi1:mesh%avi2,:) + &
                                     ddy_acc( mesh%avi1:mesh%avi2,:)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddy_c  )
    CALL deallocate_shared( wddy_aca)
    CALL deallocate_shared( wddy_acc)
    
  END SUBROUTINE ddy_a_to_ac_3D
  
  ! From the ac-grid to the a/c-grids
  SUBROUTINE move_aca_to_a_2D( mesh, d_aca, d_a)
    ! Move data from the a-part of the ac-grid to the a-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_aca
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_a
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_aca,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'move_aca_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_move_aca_a, d_aca, d_a)
    
  END SUBROUTINE move_aca_to_a_2D
  SUBROUTINE move_aca_to_a_3D( mesh, d_aca, d_a)
    ! Move data from the a-part of the ac-grid to the a-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_aca
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_a
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_aca,1) /= mesh%nVAaAc .OR. SIZE( d_aca,2) /= SIZE( d_a,2)) THEN
      IF (par%master) WRITE(0,*) 'move_aca_to_a_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_move_aca_a, d_aca, d_a)
    
  END SUBROUTINE move_aca_to_a_3D
  SUBROUTINE move_acc_to_c_2D( mesh, d_acc, d_c)
    ! Move data from the c-part of the ac-grid to the c-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_acc
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_c
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_acc,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'move_acc_to_c_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_move_acc_c, d_acc, d_c)
    
  END SUBROUTINE move_acc_to_c_2D
  SUBROUTINE move_acc_to_c_3D( mesh, d_acc, d_c)
    ! Move data from the c-part of the ac-grid to the c-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_acc
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_c
    
    ! Safety
    IF (SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_acc,1) /= mesh%nVAaAc .OR. SIZE( d_acc,2) /= SIZE( d_c,2)) THEN
      IF (par%master) WRITE(0,*) 'move_acc_to_c_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_move_acc_c, d_acc, d_c)
    
  END SUBROUTINE move_acc_to_c_3D
  SUBROUTINE move_ac_to_a_and_c_2D( mesh, d_ac, d_a, d_c)
    ! Move data from the ac-grid to the a/c-grids
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_c
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'move_ac_to_a_and_c_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_move_aca_a, d_ac, d_a)
    CALL multiply_matrix_vector_CSR( mesh%M_move_acc_c, d_ac, d_c)
    
  END SUBROUTINE move_ac_to_a_and_c_2D
  SUBROUTINE move_ac_to_a_and_c_3D( mesh, d_ac, d_a, d_c)
    ! Move data from the ac-grid to the a/c-grids
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_c
    
    ! Safety
    IF (SIZE( d_a,1) /= mesh%nV .OR. SIZE( d_c,1) /= mesh%nAc .OR. SIZE( d_ac,1) /= mesh%nVAaAc .OR. &
        SIZE( d_ac,2) /= SIZE( d_a,2) .OR. SIZE( d_a,2) /= SIZE( d_c,2)) THEN
      IF (par%master) WRITE(0,*) 'move_ac_to_a_and_c_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_move_aca_a, d_ac, d_a)
    CALL multiply_matrix_vector_2D_CSR( mesh%M_move_acc_c, d_ac, d_c)
    
  END SUBROUTINE move_ac_to_a_and_c_3D
  
! == Derivatives on the ac-grid
  SUBROUTINE ddx_ac_to_ac_2D( mesh, d_ac, ddx_ac)
    ! Calculate d/dx from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( ddx_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_ac_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_ac_ac, d_ac, ddx_ac)
    
  END SUBROUTINE ddx_ac_to_ac_2D
  SUBROUTINE ddx_ac_to_ac_3D( mesh, d_ac, ddx_ac)
    ! Calculate d/dx from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( ddx_ac,1) /= mesh%nVAaAc .OR. SIZE( d_ac,2) /= SIZE( ddx_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'ddx_ac_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_ac_ac, d_ac, ddx_ac)
    
  END SUBROUTINE ddx_ac_to_ac_3D
  SUBROUTINE ddy_ac_to_ac_2D( mesh, d_ac, ddy_ac)
    ! Calculate d/dy from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( ddy_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_ac_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_ac_ac, d_ac, ddy_ac)
    
  END SUBROUTINE ddy_ac_to_ac_2D
  SUBROUTINE ddy_ac_to_ac_3D( mesh, d_ac, ddy_ac)
    ! Calculate d/dy from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( ddy_ac,1) /= mesh%nVAaAc .OR. SIZE( d_ac,2) /= SIZE( ddy_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'ddy_ac_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_ac_ac, d_ac, ddy_ac)
    
  END SUBROUTINE ddy_ac_to_ac_3D
  SUBROUTINE d2dx2_ac_to_ac_2D( mesh, d_ac, d2dx2_ac)
    ! Calculate d2/dx2 from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dx2_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d2dx2_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'd2dx2_ac_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_d2dx2_ac_ac, d_ac, d2dx2_ac)
    
  END SUBROUTINE d2dx2_ac_to_ac_2D
  SUBROUTINE d2dx2_ac_to_ac_3D( mesh, d_ac, d2dx2_ac)
    ! Calculate d2/dx2 from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dx2_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d2dx2_ac,1) /= mesh%nVAaAc .OR. SIZE( d_ac,2) /= SIZE( d2dx2_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'd2dx2_ac_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_d2dx2_ac_ac, d_ac, d2dx2_ac)
    
  END SUBROUTINE d2dx2_ac_to_ac_3D
  SUBROUTINE d2dxdy_ac_to_ac_2D( mesh, d_ac, d2dxdy_ac)
    ! Calculate d2/dxdy from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dxdy_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d2dxdy_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'd2dxdy_ac_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_d2dxdy_ac_ac, d_ac, d2dxdy_ac)
    
  END SUBROUTINE d2dxdy_ac_to_ac_2D
  SUBROUTINE d2dxdy_ac_to_ac_3D( mesh, d_ac, d2dxdy_ac)
    ! Calculate d2/dxdy from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dxdy_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d2dxdy_ac,1) /= mesh%nVAaAc .OR. SIZE( d_ac,2) /= SIZE( d2dxdy_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'd2dxdy_ac_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_d2dxdy_ac_ac, d_ac, d2dxdy_ac)
    
  END SUBROUTINE d2dxdy_ac_to_ac_3D
  SUBROUTINE d2dy2_ac_to_ac_2D( mesh, d_ac, d2dy2_ac)
    ! Calculate d2/dy2 from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d2dy2_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d2dy2_ac,1) /= mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'd2dy2_ac_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_d2dy2_ac_ac, d_ac, d2dy2_ac)
    
  END SUBROUTINE d2dy2_ac_to_ac_2D
  SUBROUTINE d2dy2_ac_to_ac_3D( mesh, d_ac, d2dy2_ac)
    ! Calculate d2/dy2 from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d2dy2_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d2dy2_ac,1) /= mesh%nVAaAc .OR. SIZE( d_ac,2) /= SIZE( d2dy2_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'd2dy2_ac_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_d2dy2_ac_ac, d_ac, d2dy2_ac)
    
  END SUBROUTINE d2dy2_ac_to_ac_3D
  
! == Operations between the ac-grid and the combined bb-grid
  SUBROUTINE map_ac_to_bb_2D( mesh, d_ac, d_bb)
    ! Map data from the ac-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_bb
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'map_ac_to_bb_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_map_ac_bb, d_ac, d_bb)
    
  END SUBROUTINE map_ac_to_bb_2D
  SUBROUTINE map_ac_to_bb_3D( mesh, d_ac, d_bb)
    ! Map data from the ac-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_bb
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_ac,2) /= SIZE( d_bb,2)) THEN
      IF (par%master) WRITE(0,*) 'map_ac_to_bb_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_ac_bb, d_ac, d_bb)
    
  END SUBROUTINE map_ac_to_bb_3D
  SUBROUTINE map_bb_to_ac_2D( mesh, d_bb, d_ac)
    ! Map data from the bb-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'map_bb_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_map_bb_ac, d_bb, d_ac)
    
  END SUBROUTINE map_bb_to_ac_2D
  SUBROUTINE map_bb_to_ac_3D( mesh, d_bb, d_ac)
    ! Map data from the bb-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_bb,2) /= SIZE( d_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'map_bb_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_bb_ac, d_bb, d_ac)
    
  END SUBROUTINE map_bb_to_ac_3D
  
  SUBROUTINE ddx_ac_to_bb_2D( mesh, d_ac, ddx_bb)
    ! Calculate d/dx from the ac-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_bb
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( ddx_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_ac_to_bb_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_ac_bb, d_ac, ddx_bb)
    
  END SUBROUTINE ddx_ac_to_bb_2D
  SUBROUTINE ddx_ac_to_bb_3D( mesh, d_ac, ddx_bb)
    ! Calculate d/dx from the ac-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_bb
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( ddx_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_ac,2) /= SIZE( ddx_bb,2)) THEN
      IF (par%master) WRITE(0,*) 'ddx_ac_to_bb_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_ac_bb, d_ac, ddx_bb)
    
  END SUBROUTINE ddx_ac_to_bb_3D
  SUBROUTINE ddx_bb_to_ac_2D( mesh, d_bb, ddx_ac)
    ! Calculate d/dx from the bb-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_ac
    
    ! Safety
    IF (SIZE( ddx_ac,1) /= mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_bb_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_bb_ac, d_bb, ddx_ac)
    
  END SUBROUTINE ddx_bb_to_ac_2D
  SUBROUTINE ddx_bb_to_ac_3D( mesh, d_bb, ddx_ac)
    ! Calculate d/dx from the bb-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_ac
    
    ! Safety
    IF (SIZE( ddx_ac,1) /= mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_bb,2) /= SIZE( ddx_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'ddx_bb_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_bb_ac, d_bb, ddx_ac)
    
  END SUBROUTINE ddx_bb_to_ac_3D
  
  SUBROUTINE ddy_ac_to_bb_2D( mesh, d_ac, ddy_bb)
    ! Calculate d/dy from the ac-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_bb
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( ddy_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_ac_to_bb_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_ac_bb, d_ac, ddy_bb)
    
  END SUBROUTINE ddy_ac_to_bb_2D
  SUBROUTINE ddy_ac_to_bb_3D( mesh, d_ac, ddy_bb)
    ! Calculate d/dy from the ac-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_bb
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( ddy_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_ac,2) /= SIZE( ddy_bb,2)) THEN
      IF (par%master) WRITE(0,*) 'ddy_ac_to_bb_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_ac_bb, d_ac, ddy_bb)
    
  END SUBROUTINE ddy_ac_to_bb_3D
  SUBROUTINE ddy_bb_to_ac_2D( mesh, d_bb, ddy_ac)
    ! Calculate d/dy from the bb-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_ac
    
    ! Safety
    IF (SIZE( ddy_ac,1) /= mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_bb_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_bb_ac, d_bb, ddy_ac)
    
  END SUBROUTINE ddy_bb_to_ac_2D
  SUBROUTINE ddy_bb_to_ac_3D( mesh, d_bb, ddy_ac)
    ! Calculate d/dy from the bb-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_ac
    
    ! Safety
    IF (SIZE( ddy_ac,1) /= mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_bb,2) /= SIZE( ddy_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'ddy_bb_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_bb_ac, d_bb, ddy_ac)
    
  END SUBROUTINE ddy_bb_to_ac_3D
  
! == Operations between the ac-grid and the acu/acv-grids

  ! From the ac-grid to the acu/acv-grids
  SUBROUTINE move_ac_to_acu_2D( mesh, d_ac, d_acu)
    ! Move data from the ac-grid to the u-part of the acuv-grid
    !
    ! NOTE: sets v-part to zero!
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_acu
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_acu,1) /= 2*mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'move_ac_to_acu_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_move_ac_acu, d_ac, d_acu)
    
  END SUBROUTINE move_ac_to_acu_2D
  SUBROUTINE move_ac_to_acu_3D( mesh, d_ac, d_acu)
    ! Move data from the ac-grid to the u-part of the acuv-grid
    !
    ! NOTE: sets v-part to zero!
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_acu
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( d_ac,2) /= SIZE( d_acu,2)) THEN
      IF (par%master) WRITE(0,*) 'move_ac_to_acu_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_move_ac_acu, d_ac, d_acu)
    
  END SUBROUTINE move_ac_to_acu_3D
  SUBROUTINE move_ac_to_acv_2D( mesh, d_ac, d_acv)
    ! Move data from the ac-grid to the v-part of the acuv-grid
    !
    ! NOTE: sets u-part to zero!
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_acv
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_acv,1) /= 2*mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'move_ac_to_acv_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_move_ac_acv, d_ac, d_acv)
    
  END SUBROUTINE move_ac_to_acv_2D
  SUBROUTINE move_ac_to_acv_3D( mesh, d_ac, d_acv)
    ! Move data from the ac-grid to the v-part of the acuv-grid
    !
    ! NOTE: sets u-part to zero!
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_ac
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_acv
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( d_ac,2) /= SIZE( d_acv,2)) THEN
      IF (par%master) WRITE(0,*) 'move_ac_to_acv_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_move_ac_acv, d_ac, d_acv)
    
  END SUBROUTINE move_ac_to_acv_3D
  
  ! From the acu/acv-grids to the ac-grid
  SUBROUTINE move_acu_to_ac_2D( mesh, d_acu, d_ac)
    ! Move data from the u-part of the acuv-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_acu
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_acu,1) /= 2*mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'move_acu_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_move_acu_ac, d_acu, d_ac)
    
  END SUBROUTINE move_acu_to_ac_2D
  SUBROUTINE move_acu_to_ac_3D( mesh, d_acu, d_ac)
    ! Move data from the u-part of the acuv-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_acu
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( d_acu,2) /= SIZE( d_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'move_acu_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_move_acu_ac, d_acu, d_ac)
    
  END SUBROUTINE move_acu_to_ac_3D
  SUBROUTINE move_acv_to_ac_2D( mesh, d_acv, d_ac)
    ! Move data from the v-part of the acuv-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_acv
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_acv,1) /= 2*mesh%nVAaAc) THEN
      IF (par%master) WRITE(0,*) 'move_acv_to_ac_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_move_acv_ac, d_acv, d_ac)
    
  END SUBROUTINE move_acv_to_ac_2D
  SUBROUTINE move_acv_to_ac_3D( mesh, d_acv, d_ac)
    ! Move data from the v-part of the acuv-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_acv
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_ac
    
    ! Safety
    IF (SIZE( d_ac,1) /= mesh%nVAaAc .OR. SIZE( d_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( d_acv,2) /= SIZE( d_ac,2)) THEN
      IF (par%master) WRITE(0,*) 'move_acv_to_ac_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_move_acv_ac, d_acv, d_ac)
    
  END SUBROUTINE move_acv_to_ac_3D
  
! == Operations between the acu/acv-grids and the bb-grid

  ! From the acu/acv-grids to the bb-grid
  SUBROUTINE map_acu_to_bb_2D( mesh, d_acu, d_bb)
    ! Map data from the u-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_acu
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_bb
    
    ! Safety
    IF (SIZE( d_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'map_acu_to_bb_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_map_acu_bb, d_acu, d_bb)
    
  END SUBROUTINE map_acu_to_bb_2D
  SUBROUTINE map_acu_to_bb_3D( mesh, d_acu, d_bb)
    ! Map data from the u-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_acu
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_bb
    
    ! Safety
    IF (SIZE( d_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_acu,2) /= SIZE( d_bb,2)) THEN
      IF (par%master) WRITE(0,*) 'map_acu_to_bb_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_acu_bb, d_acu, d_bb)
    
  END SUBROUTINE map_acu_to_bb_3D
  SUBROUTINE map_acv_to_bb_2D( mesh, d_acv, d_bb)
    ! Map data from the v-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_acv
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_bb
    
    ! Safety
    IF (SIZE( d_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'map_acv_to_bb_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_map_acv_bb, d_acv, d_bb)
    
  END SUBROUTINE map_acv_to_bb_2D
  SUBROUTINE map_acv_to_bb_3D( mesh, d_acv, d_bb)
    ! Map data from the v-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_acv
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_bb
    
    ! Safety
    IF (SIZE( d_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_acv,2) /= SIZE( d_bb,2)) THEN
      IF (par%master) WRITE(0,*) 'map_acv_to_bb_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_acv_bb, d_acv, d_bb)
    
  END SUBROUTINE map_acv_to_bb_3D
  
  SUBROUTINE ddx_acu_to_bb_2D( mesh, d_acu, ddx_bb)
    ! Calculate d/dx from the u-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_acu
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_bb
    
    ! Safety
    IF (SIZE( d_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( ddx_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_acu_to_bb_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_acu_bb, d_acu, ddx_bb)
    
  END SUBROUTINE ddx_acu_to_bb_2D
  SUBROUTINE ddx_acu_to_bb_3D( mesh, d_acu, ddx_bb)
    ! Calculate d/dx from the u-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_acu
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_bb
    
    ! Safety
    IF (SIZE( d_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( ddx_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_acu,2) /= SIZE( ddx_bb,2)) THEN
      IF (par%master) WRITE(0,*) 'ddx_acu_to_bb_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_acu_bb, d_acu, ddx_bb)
    
  END SUBROUTINE ddx_acu_to_bb_3D
  SUBROUTINE ddx_acv_to_bb_2D( mesh, d_acv, ddx_bb)
    ! Calculate d/dx from the v-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_acv
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_bb
    
    ! Safety
    IF (SIZE( d_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( ddx_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_acv_to_bb_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_acv_bb, d_acv, ddx_bb)
    
  END SUBROUTINE ddx_acv_to_bb_2D
  SUBROUTINE ddx_acv_to_bb_3D( mesh, d_acv, ddx_bb)
    ! Calculate d/dx from the v-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_acv
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_bb
    
    ! Safety
    IF (SIZE( d_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( ddx_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_acv,2) /= SIZE( ddx_bb,2)) THEN
      IF (par%master) WRITE(0,*) 'ddx_acv_to_bb_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_acv_bb, d_acv, ddx_bb)
    
  END SUBROUTINE ddx_acv_to_bb_3D
  
  SUBROUTINE ddy_acu_to_bb_2D( mesh, d_acu, ddy_bb)
    ! Calculate d/dy from the u-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_acu
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_bb
    
    ! Safety
    IF (SIZE( d_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( ddy_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_acu_to_bb_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_acu_bb, d_acu, ddy_bb)
    
  END SUBROUTINE ddy_acu_to_bb_2D
  SUBROUTINE ddy_acu_to_bb_3D( mesh, d_acu, ddy_bb)
    ! Calculate d/dy from the u-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_acu
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_bb
    
    ! Safety
    IF (SIZE( d_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( ddy_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_acu,2) /= SIZE( ddy_bb,2)) THEN
      IF (par%master) WRITE(0,*) 'ddy_acu_to_bb_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_acu_bb, d_acu, ddy_bb)
    
  END SUBROUTINE ddy_acu_to_bb_3D
  SUBROUTINE ddy_acv_to_bb_2D( mesh, d_acv, ddy_bb)
    ! Calculate d/dy from the v-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_acv
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_bb
    
    ! Safety
    IF (SIZE( d_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( ddy_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_acv_to_bb_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_acv_bb, d_acv, ddy_bb)
    
  END SUBROUTINE ddy_acv_to_bb_2D
  SUBROUTINE ddy_acv_to_bb_3D( mesh, d_acv, ddy_bb)
    ! Calculate d/dy from the v-part of the acuv-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_acv
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_bb
    
    ! Safety
    IF (SIZE( d_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( ddy_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_acv,2) /= SIZE( ddy_bb,2)) THEN
      IF (par%master) WRITE(0,*) 'ddy_acv_to_bb_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_acv_bb, d_acv, ddy_bb)
    
  END SUBROUTINE ddy_acv_to_bb_3D
  
  ! From the bb-grid to the acu/acv-grids
  SUBROUTINE map_bb_to_acu_2D( mesh, d_bb, d_acu)
    ! Map data from the bb-grid to the u-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_acu
    
    ! Safety
    IF (SIZE( d_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'map_bb_to_acu_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_map_bb_acu, d_bb, d_acu)
    
  END SUBROUTINE map_bb_to_acu_2D
  SUBROUTINE map_bb_to_acu_3D( mesh, d_bb, d_acu)
    ! Map data from the bb-grid to the u-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_acu
    
    ! Safety
    IF (SIZE( d_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_bb,2) /= SIZE( d_acu,2)) THEN
      IF (par%master) WRITE(0,*) 'map_bb_to_acu_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_bb_acu, d_bb, d_acu)
    
  END SUBROUTINE map_bb_to_acu_3D
  SUBROUTINE map_bb_to_acv_2D( mesh, d_bb, d_acv)
    ! Map data from the bb-grid to the v-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_acv
    
    ! Safety
    IF (SIZE( d_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'map_bb_to_acv_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_map_bb_acv, d_bb, d_acv)
    
  END SUBROUTINE map_bb_to_acv_2D
  SUBROUTINE map_bb_to_acv_3D( mesh, d_bb, d_acv)
    ! Map data from the bb-grid to the v-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_acv
    
    ! Safety
    IF (SIZE( d_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_bb,2) /= SIZE( d_acv,2)) THEN
      IF (par%master) WRITE(0,*) 'map_bb_to_acv_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_map_bb_acv, d_bb, d_acv)
    
  END SUBROUTINE map_bb_to_acv_3D
  
  SUBROUTINE ddx_bb_to_acu_2D( mesh, d_bb, ddx_acu)
    ! Calculate d/dx from the bb-grid to the u-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_acu
    
    ! Safety
    IF (SIZE( ddx_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_bb_to_acu_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_bb_acu, d_bb, ddx_acu)
    
  END SUBROUTINE ddx_bb_to_acu_2D
  SUBROUTINE ddx_bb_to_acu_3D( mesh, d_bb, ddx_acu)
    ! Calculate d/dx from the bb-grid to the u-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_acu
    
    ! Safety
    IF (SIZE( ddx_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_bb,2) /= SIZE( ddx_acu,2)) THEN
      IF (par%master) WRITE(0,*) 'ddx_bb_to_acu_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_bb_acu, d_bb, ddx_acu)
    
  END SUBROUTINE ddx_bb_to_acu_3D
  SUBROUTINE ddx_bb_to_acv_2D( mesh, d_bb, ddx_acv)
    ! Calculate d/dx from the bb-grid to the v-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddx_acv
    
    ! Safety
    IF (SIZE( ddx_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddx_bb_to_acv_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_bb_acv, d_bb, ddx_acv)
    
  END SUBROUTINE ddx_bb_to_acv_2D
  SUBROUTINE ddx_bb_to_acv_3D( mesh, d_bb, ddx_acv)
    ! Calculate d/dx from the bb-grid to the v-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddx_acv
    
    ! Safety
    IF (SIZE( ddx_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_bb,2) /= SIZE( ddx_acv,2)) THEN
      IF (par%master) WRITE(0,*) 'ddx_bb_to_acv_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_bb_acv, d_bb, ddx_acv)
    
  END SUBROUTINE ddx_bb_to_acv_3D
  
  SUBROUTINE ddy_bb_to_acu_2D( mesh, d_bb, ddy_acu)
    ! Calculate d/dy from the bb-grid to the u-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_acu
    
    ! Safety
    IF (SIZE( ddy_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_bb_to_acu_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_bb_acu, d_bb, ddy_acu)
    
  END SUBROUTINE ddy_bb_to_acu_2D
  SUBROUTINE ddy_bb_to_acu_3D( mesh, d_bb, ddy_acu)
    ! Calculate d/dy from the bb-grid to the u-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_acu
    
    ! Safety
    IF (SIZE( ddy_acu,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_bb,2) /= SIZE( ddy_acu,2)) THEN
      IF (par%master) WRITE(0,*) 'ddy_bb_to_acu_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_bb_acu, d_bb, ddy_acu)
    
  END SUBROUTINE ddy_bb_to_acu_3D
  SUBROUTINE ddy_bb_to_acv_2D( mesh, d_bb, ddy_acv)
    ! Calculate d/dy from the bb-grid to the v-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: ddy_acv
    
    ! Safety
    IF (SIZE( ddy_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc) THEN
      IF (par%master) WRITE(0,*) 'ddy_bb_to_acv_2D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_bb_acv, d_bb, ddy_acv)
    
  END SUBROUTINE ddy_bb_to_acv_2D
  SUBROUTINE ddy_bb_to_acv_3D( mesh, d_bb, ddy_acv)
    ! Calculate d/dy from the bb-grid to the v-part of the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_bb
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ddy_acv
    
    ! Safety
    IF (SIZE( ddy_acv,1) /= 2*mesh%nVAaAc .OR. SIZE( d_bb,1) /= mesh%nTriAaAc .OR. SIZE( d_bb,2) /= SIZE( ddy_acv,2)) THEN
      IF (par%master) WRITE(0,*) 'ddy_bb_to_acv_3D - ERROR: data fields are the wrong size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Perform the operation as a matrix multiplication
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_bb_acv, d_bb, ddy_acv)
    
  END SUBROUTINE ddy_bb_to_acv_3D
  
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
    
    ! Calculate matrix operators for the second partial derivatives on the a-grid (d2/dx2, d2/dxdy, d2/dy2)
    CALL multiply_matrix_matrix_CSR( mesh%M_ddx_b_a, mesh%M_ddx_a_b, mesh%M_d2dx2_a_a )
    CALL multiply_matrix_matrix_CSR( mesh%M_ddy_b_a, mesh%M_ddx_a_b, mesh%M_d2dxdy_a_a)
    CALL multiply_matrix_matrix_CSR( mesh%M_ddy_b_a, mesh%M_ddy_a_b, mesh%M_d2dy2_a_a )
    
    ! Calculate matrix operators involving the "combined" mesh (vertex+edge ac-grid, triangle-X-4 bb-grid)
    CALL calc_matrix_operators_move_to_combi( mesh)
    
    CALL calc_matrix_operators_ac_ac( mesh,                   mesh%M_ddx_ac_ac, mesh%M_ddy_ac_ac)
    CALL calc_matrix_operators_ac_bb( mesh, mesh%M_map_ac_bb, mesh%M_ddx_ac_bb, mesh%M_ddy_ac_bb)
    CALL calc_matrix_operators_bb_ac( mesh, mesh%M_map_bb_ac, mesh%M_ddx_bb_ac, mesh%M_ddy_bb_ac)
    
    CALL multiply_matrix_matrix_CSR( mesh%M_ddx_bb_ac, mesh%M_ddx_ac_bb, mesh%M_d2dx2_ac_ac )
    CALL multiply_matrix_matrix_CSR( mesh%M_ddy_bb_ac, mesh%M_ddx_ac_bb, mesh%M_d2dxdy_ac_ac)
    CALL multiply_matrix_matrix_CSR( mesh%M_ddy_bb_ac, mesh%M_ddy_ac_bb, mesh%M_d2dy2_ac_ac )
    
    CALL calc_matrix_operators_combi_acuv( mesh)
    
  END SUBROUTINE calc_matrix_operators_mesh
  
  SUBROUTINE calc_matrix_operators_a_a( mesh,        M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: vi, n, ci, vj, i, k
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc

    ncols           = mesh%nV      ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_max = mesh%nC_mem+1
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%vi1:mesh%vi2) = 1

    k = 0
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO vi = mesh%vi1, mesh%vi2
      
      ! Source points: the neighbouring vertices
      n = 0
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
!        IF (mask_valid_in( vj) == 0) CYCLE
        n = n+1
        i_c( n) = vj
        x_c( n) = mesh%V( vj,1)
        y_c( n) = mesh%V( vj,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 2) THEN
!        mask_valid_out( vi) = 0
        CYCLE
      END IF
      
      ! Destination point: the vertex
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_reg( x, y, n, x_c, y_c, Nfxi, Nfyi, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO
      
      k = k+1
      M_ddx%index( k) = vi
      M_ddy%index( k) = vi
      M_ddx%val(   k) = Nfxi
      M_ddy%val(   k) = Nfyi

      ! Finish this row
      M_ddx%ptr( vi+1 : nrows+1) = k+1
      M_ddy%ptr( vi+1 : nrows+1) = k+1
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Number of non-zero elements
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%vi1, mesh%vi2)
    
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
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ti, n, vii, vi, i, k
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nV      ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 3
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%ti1:mesh%ti2) = 1

    k = 0
    M_map%ptr = 1
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO ti = mesh%ti1, mesh%ti2
      
      ! Source points: the neighbouring vertices
      n = 0
      DO vii = 1, 3
        vi = mesh%Tri( ti, vii)
!        IF (mask_valid_in( vi) == 0) CYCLE
        n = n+1
        i_c( n) = vi
        x_c( n) = mesh%V( vi,1)
        y_c( n) = mesh%V( vi,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) THEN
!        mask_valid_out( ti) = 0
        CYCLE
      END IF
      
      ! Destination point: the triangle's geometric centre
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_map%index( k) = i_c(  i)
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_map%val(   k) = Nfc(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO

      ! Finish this row
      M_map%ptr( ti+1 : nrows+1) = k+1
      M_ddx%ptr( ti+1 : nrows+1) = k+1
      M_ddy%ptr( ti+1 : nrows+1) = k+1
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Number of non-zero elements
    M_map%nnz = k
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ti1, mesh%ti2)
    
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
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: aci, n, acii, vi, i, k
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nV      ! from
    nrows           = mesh%nAc     ! to
    nnz_per_row_max = 4
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%ci1:mesh%ci2) = 1

    k = 0
    M_map%ptr = 1
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO aci = mesh%ci1, mesh%ci2
      
      ! Source points: the three/four regular vertices surrounding the staggered vertex
      n = 0
      DO acii = 1, 4
        vi = mesh%Aci( aci, acii)
        IF (vi == 0) CYCLE
!        IF (mask_valid_in( vi) == 0) CYCLE
        n = n+1
        i_c( n) = vi
        x_c( n) = mesh%V( vi,1)
        y_c( n) = mesh%V( vi,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) THEN
!        mask_valid_out( aci) = 0
        CYCLE
      END IF
      
      ! Destination point: the staggered vertex
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_map%index( k) = i_c(  i)
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_map%val(   k) = Nfc(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO

      ! Finish this row
      M_map%ptr( aci+1 : nrows+1) = k+1
      M_ddx%ptr( aci+1 : nrows+1) = k+1
      M_ddy%ptr( aci+1 : nrows+1) = k+1
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync
    
    ! Number of non-zero elements
    M_map%nnz = k
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%ci1, mesh%ci2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ci1, mesh%ci2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ci1, mesh%ci2)
    
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
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: vi, n, iti, ti, i, k
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_max = mesh%nC_mem
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%vi1:mesh%vi2) = 1

    k = 0
    M_map%ptr = 1
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO vi = mesh%vi1, mesh%vi2
      
      ! Source points: the geometric centres of the triangles surrounding the regular vertex
      n = 0
      DO iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi, iti)
!        IF (mask_valid_in( ti) == 0) CYCLE
        n = n+1
        i_c( n) = ti
        x_c( n) = mesh%TriGC( ti,1)
        y_c( n) = mesh%TriGC( ti,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) THEN
!        mask_valid_out( vi) = 0
        CYCLE
      END IF
      
      ! Destination point: the regular vertex
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_map%index( k) = i_c(  i)
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_map%val(   k) = Nfc(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO

      ! Finish this row
      M_map%ptr( vi+1 : nrows+1) = k+1
      M_ddx%ptr( vi+1 : nrows+1) = k+1
      M_ddy%ptr( vi+1 : nrows+1) = k+1
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Number of non-zero elements
    M_map%nnz = k
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%vi1, mesh%vi2)
    
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
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ti, n, tii, tj, i, k
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 4
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%ti1:mesh%ti2) = 1

    k = 0
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO ti = mesh%ti1, mesh%ti2
      
      ! Source points: the up to three adjacent triangles's geometric centres
      n = 0
      DO tii = 1, 3
        tj = mesh%TriC( ti,tii)
        IF (tj == 0) CYCLE
!        IF (mask_valid_in( tj) == 0) CYCLE
        n = n+1
        i_c( n) = tj
        x_c( n) = mesh%TriGC( tj,1)
        y_c( n) = mesh%TriGC( tj,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 2) THEN
!        mask_valid_out( ti) = 0
        CYCLE
      END IF
      
      ! Destination point: the triangle's geometric centre
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_reg( x, y, n, x_c, y_c, Nfxi, Nfyi, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO
      
      k = k+1
      M_ddx%index( k) = ti
      M_ddy%index( k) = ti
      M_ddx%val(   k) = Nfxi
      M_ddy%val(   k) = Nfyi

      ! Finish this row
      M_ddx%ptr( ti+1 : nrows+1) = k+1
      M_ddy%ptr( ti+1 : nrows+1) = k+1
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Number of non-zero elements
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ti1, mesh%ti2)
    
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
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: aci, n, acii, vi, iti, ti, n_match, vii, aciii, i, k
    LOGICAL                                            :: is_listed
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nAc     ! to
    nnz_per_row_max = 6
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%ci1:mesh%ci2) = 1

    k = 0
    M_map%ptr = 1
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO aci = mesh%ci1, mesh%ci2
      
      ! Source points: the two triangles adjacent to the edge, and the (up to) four other triangles adjoining them
      n = 0
      DO acii = 1, 4
        vi = mesh%Aci( aci,acii)
        IF (vi == 0) CYCLE
        
        DO iti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,iti)
          
!          IF (mask_valid_in( ti) == 0) CYCLE
          
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
      IF (n < 3) THEN
!        mask_valid_out( aci) = 0
        CYCLE
      END IF
      
      ! Destination point: the staggered vertex
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_map%index( k) = i_c(  i)
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_map%val(   k) = Nfc(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO

      ! Finish this row
      M_map%ptr( aci+1 : nrows+1) = k+1
      M_ddx%ptr( aci+1 : nrows+1) = k+1
      M_ddy%ptr( aci+1 : nrows+1) = k+1
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync
    
    ! Number of non-zero elements
    M_map%nnz = k
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%ci1, mesh%ci2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ci1, mesh%ci2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ci1, mesh%ci2)
    
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
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: vi, n, ci, aci, i, k
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nAc     ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_max = mesh%nC_mem
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%vi1:mesh%vi2) = 1

    k = 0
    M_map%ptr = 1
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO vi = mesh%vi1, mesh%vi2
      
      ! Source points: the staggered vertices surrounding the regular vertex
      n = 0
      DO ci = 1, mesh%nC( vi)
        aci = mesh%iAci( vi, ci)
!        IF (mask_valid_in( aci) == 0) CYCLE
        n = n+1
        i_c( n) = aci
        x_c( n) = mesh%VAc( aci,1)
        y_c( n) = mesh%VAc( aci,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) THEN
!        mask_valid_out( vi) = 0
        CYCLE
      END IF
      
      ! Destination point: the regular vertex
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_map%index( k) = i_c(  i)
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_map%val(   k) = Nfc(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO

      ! Finish this row
      M_map%ptr( vi+1 : nrows+1) = k+1
      M_ddx%ptr( vi+1 : nrows+1) = k+1
      M_ddy%ptr( vi+1 : nrows+1) = k+1
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Number of non-zero elements
    M_map%nnz = k
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%vi1, mesh%vi2)
    
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
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ti, n, via, vib, vic, ci, vj, acab, acbc, acca, i, k
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nAc     ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 3
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%ti1:mesh%ti2) = 1

    k = 0
    M_map%ptr = 1
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO ti = mesh%ti1, mesh%ti2
      
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
!      IF (mask_valid_in( acab) == 1) THEN
        n = n+1
        i_c( n) = acab
        x_c( n) = mesh%VAc( acab,1)
        y_c( n) = mesh%VAc( acab,2)
!      END IF
      
      acbc = 0
      DO ci = 1, mesh%nC( vib)
        vj = mesh%C( vib,ci)
        IF (vj == vic) THEN
          acbc = mesh%iAci( vib,ci)
          EXIT
        END IF
      END DO
!      IF (mask_valid_in( acbc) == 1) THEN
        n = n+1
        i_c( n) = acbc
        x_c( n) = mesh%VAc( acbc,1)
        y_c( n) = mesh%VAc( acbc,2)
!      END IF
      
      acca = 0
      DO ci = 1, mesh%nC( vic)
        vj = mesh%C( vic,ci)
        IF (vj == via) THEN
          acca = mesh%iAci( vic,ci)
          EXIT
        END IF
      END DO
!      IF (mask_valid_in( acca) == 1) THEN
        n = n+1
        i_c( n) = acca
        x_c( n) = mesh%VAc( acca,1)
        y_c( n) = mesh%VAc( acca,2)
!      END IF
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) THEN
        CYCLE
      END IF
      
      ! Destination point: the triangle's geometric centre
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_map%index( k) = i_c(  i)
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_map%val(   k) = Nfc(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO

      ! Finish this row
      M_map%ptr( ti+1 : nrows+1) = k+1
      M_ddx%ptr( ti+1 : nrows+1) = k+1
      M_ddy%ptr( ti+1 : nrows+1) = k+1
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Number of non-zero elements
    M_map%nnz = k
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ti1, mesh%ti2)
    
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
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: aci, n, acii, vi, ci, vc, acj, i, k
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc

    ncols           = mesh%nAc     ! from
    nrows           = mesh%nAc     ! to
    nnz_per_row_max = 5
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%ci1:mesh%ci2) = 1

    k = 0
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO aci = mesh%ci1, mesh%ci2
      
      ! Source points: all the staggered vertices lying on the edges of the
      ! two triangles adjoining the staggered vertex
      n = 0
      DO acii = 1, 2
        vi = mesh%Aci( aci,acii)
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          acj = mesh%iAci( vi,ci)
!          IF (mask_valid_in( acj) == 0) CYCLE
          IF (vc == mesh%Aci( aci,3) .OR. vc == mesh%Aci( aci,4)) THEN
            n = n+1
            i_c( n) = acj
            x_c( n) = mesh%VAc( acj,1)
            y_c( n) = mesh%VAc( acj,2)
          END IF
        END DO
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 2) THEN
!        mask_valid_out( aci) = 0
        CYCLE
      END IF
      
      ! Destination point: the staggered vertex
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_reg( x, y, n, x_c, y_c, Nfxi, Nfyi, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO
      
      k = k+1
      M_ddx%index( k) = aci
      M_ddy%index( k) = aci
      M_ddx%val(   k) = Nfxi
      M_ddy%val(   k) = Nfyi

      ! Finish this row
      M_ddx%ptr( aci+1 : nrows+1) = k+1
      M_ddy%ptr( aci+1 : nrows+1) = k+1
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync
    
    ! Number of non-zero elements
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ci1, mesh%ci2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ci1, mesh%ci2)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_c_c
  
  SUBROUTINE calc_matrix_operators_move_to_combi( mesh)
    ! Calculate matrix operators representing moving operations
    ! between the a/c-grids and the combined ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:
    INTEGER                                            :: nrows, ncols, nnz_max
    INTEGER                                            :: vi, aci, avi
    
  ! a_aca: from the regular mesh a (vertex) grid to the a-part of the combined mesh (ac-grid)
  ! =========================================================================================
    
    ncols   = mesh%nV       ! from
    nrows   = mesh%nVAaAc   ! to
    nnz_max = mesh%nV
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( mesh%M_move_a_aca, nrows, ncols, nnz_max)
    
    ! Initialise
    mesh%M_move_a_aca%nnz = 0
    mesh%M_move_a_aca%ptr = 1
    
    ! Fill the matrices
    DO avi = mesh%avi1, mesh%avi2
    
      IF (avi <= mesh%nV) THEN
        vi = avi
        mesh%M_move_a_aca%nnz  = mesh%M_move_a_aca%nnz + 1
        mesh%M_move_a_aca%index( mesh%M_move_a_aca%nnz) = vi
        mesh%M_move_a_aca%val(   mesh%M_move_a_aca%nnz) = 1._dp
        mesh%M_move_a_aca%ptr( avi+1 : nrows+1) = mesh%M_move_a_aca%nnz + 1
      END IF
      
    END DO ! DO avi = mesh%avi1, mesh%avi2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( mesh%M_move_a_aca, mesh%avi1, mesh%avi2)
    
  ! aca_a: from the a-part of the combined mesh (ac-grid) to the regular mesh a (vertex) grid
  ! =========================================================================================
    
    ncols   = mesh%nVAaAc   ! from
    nrows   = mesh%nV       ! to
    nnz_max = mesh%nV
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( mesh%M_move_aca_a, nrows, ncols, nnz_max)
    
    ! Initialise
    mesh%M_move_aca_a%nnz = 0
    mesh%M_move_aca_a%ptr = 1
    
    ! Fill the matrices
    DO vi = mesh%vi1, mesh%vi2
    
      avi = vi
      mesh%M_move_aca_a%nnz  = mesh%M_move_aca_a%nnz + 1
      mesh%M_move_aca_a%index( mesh%M_move_aca_a%nnz) = avi
      mesh%M_move_aca_a%val(   mesh%M_move_aca_a%nnz) = 1._dp
      mesh%M_move_aca_a%ptr( vi+1 : nrows+1) = mesh%M_move_aca_a%nnz + 1
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( mesh%M_move_aca_a, mesh%vi1, mesh%vi2)
    
  ! c_acc: from the regular mesh c (edge) grid to the c-part of the combined mesh (ac-grid)
  ! =======================================================================================
    
    ncols   = mesh%nAc      ! from
    nrows   = mesh%nVAaAc   ! to
    nnz_max = mesh%nAc
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( mesh%M_move_c_acc, nrows, ncols, nnz_max)
    
    ! Initialise
    mesh%M_move_c_acc%nnz = 0
    mesh%M_move_c_acc%ptr = 1
    
    ! Fill the matrices
    DO avi = mesh%avi1, mesh%avi2
    
      IF (avi > mesh%nV) THEN
        aci = avi - mesh%nV
        mesh%M_move_c_acc%nnz  = mesh%M_move_c_acc%nnz + 1
        mesh%M_move_c_acc%index( mesh%M_move_c_acc%nnz) = aci
        mesh%M_move_c_acc%val(   mesh%M_move_c_acc%nnz) = 1._dp
        mesh%M_move_c_acc%ptr( avi+1 : nrows+1) = mesh%M_move_c_acc%nnz + 1
      END IF
      
    END DO ! DO avi = mesh%avi1, mesh%avi2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( mesh%M_move_c_acc, mesh%avi1, mesh%avi2)
    
  ! acc_c: from the c-part of the combined mesh (ac-grid) to the regular mesh c (edge) grid
  ! =======================================================================================
    
    ncols   = mesh%nVAaAc   ! from
    nrows   = mesh%nAc      ! to
    nnz_max = mesh%nAc
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( mesh%M_move_acc_c, nrows, ncols, nnz_max)
    
    ! Initialise
    mesh%M_move_acc_c%nnz = 0
    mesh%M_move_acc_c%ptr = 1
    
    ! Fill the matrices
    DO aci = mesh%ci1, mesh%ci2
    
      avi = aci + mesh%nV
      mesh%M_move_acc_c%nnz  = mesh%M_move_acc_c%nnz + 1
      mesh%M_move_acc_c%index( mesh%M_move_acc_c%nnz) = avi
      mesh%M_move_acc_c%val(   mesh%M_move_acc_c%nnz) = 1._dp
      mesh%M_move_acc_c%ptr( vi+1 : nrows+1) = mesh%M_move_acc_c%nnz + 1
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( mesh%M_move_acc_c, mesh%ci1, mesh%ci2)
    
  END SUBROUTINE calc_matrix_operators_move_to_combi
  SUBROUTINE calc_matrix_operators_ac_ac( mesh,        M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the
    ! d/dx and d/dy operations from the ac-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: avi, n, ci, avj, i, k
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc

    ncols           = mesh%nVAaAc      ! from
    nrows           = mesh%nVAaAc      ! to
    nnz_per_row_max = mesh%nC_mem+1
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%vi1:mesh%vi2) = 1

    k = 0
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO avi = mesh%avi1, mesh%avi2
      
      ! Source points: the neighbouring vertices
      n = 0
      DO ci = 1, mesh%nCAaAc( avi)
        avj = mesh%CAaAc( avi,ci)
!        IF (mask_valid_in( vj) == 0) CYCLE
        n = n+1
        i_c( n) = avj
        x_c( n) = mesh%VAaAc( avj,1)
        y_c( n) = mesh%VAaAc( avj,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 2) THEN
!        mask_valid_out( vi) = 0
        CYCLE
      END IF
      
      ! Destination point: the vertex
      x = mesh%VAaAc( avi,1)
      y = mesh%VAaAc( avi,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_reg( x, y, n, x_c, y_c, Nfxi, Nfyi, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO
      
      k = k+1
      M_ddx%index( k) = avi
      M_ddy%index( k) = avi
      M_ddx%val(   k) = Nfxi
      M_ddy%val(   k) = Nfyi

      ! Finish this row
      M_ddx%ptr( avi+1 : nrows+1) = k+1
      M_ddy%ptr( avi+1 : nrows+1) = k+1
      
    END DO ! DO avi = mesh%avi1, mesh%avi2
    CALL sync
    
    ! Number of non-zero elements
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%avi1, mesh%avi2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%avi1, mesh%avi2)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_ac_ac
  SUBROUTINE calc_matrix_operators_ac_bb( mesh, M_map, M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the ac-grid to the bb-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ati, n, vii, avi, i, k
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nVAaAc      ! from
    nrows           = mesh%nTriAaAc    ! to
    nnz_per_row_max = 3
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%ti1:mesh%ti2) = 1

    k = 0
    M_map%ptr = 1
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO ati = mesh%ati1, mesh%ati2
      
      ! Source points: the neighbouring vertices
      n = 0
      DO vii = 1, 3
        avi = mesh%TriAaAc( ati, vii)
!        IF (mask_valid_in( vi) == 0) CYCLE
        n = n+1
        i_c( n) = avi
        x_c( n) = mesh%VAaAc( avi,1)
        y_c( n) = mesh%VAaAc( avi,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) THEN
!        mask_valid_out( ti) = 0
        CYCLE
      END IF
      
      ! Destination point: the triangle's geometric centre
      x = mesh%TriGCAaAc( ati,1)
      y = mesh%TriGCAaAc( ati,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_map%index( k) = i_c(  i)
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_map%val(   k) = Nfc(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO

      ! Finish this row
      M_map%ptr( ati+1 : nrows+1) = k+1
      M_ddx%ptr( ati+1 : nrows+1) = k+1
      M_ddy%ptr( ati+1 : nrows+1) = k+1
      
    END DO ! DO ati = mesh%ati1, mesh%ati2
    CALL sync
    
    ! Number of non-zero elements
    M_map%nnz = k
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%ati1, mesh%ati2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ati1, mesh%ati2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ati1, mesh%ati2)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfc )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_ac_bb
  SUBROUTINE calc_matrix_operators_bb_ac( mesh, M_map, M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the bb-grid to the ac-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: avi, n, iti, ati, i, k
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nTriAaAc    ! from
    nrows           = mesh%nVAaAc      ! to
    nnz_per_row_max = mesh%nC_mem
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max)

!    mask_valid_out( mesh%vi1:mesh%vi2) = 1

    k = 0
    M_map%ptr = 1
    M_ddx%ptr = 1
    M_ddy%ptr = 1

    DO avi = mesh%avi1, mesh%avi2
      
      ! Source points: the geometric centres of the triangles surrounding the regular vertex
      n = 0
      DO iti = 1, mesh%niTriAaAc( avi)
        ati = mesh%iTriAaAc( avi, iti)
!        IF (mask_valid_in( ti) == 0) CYCLE
        n = n+1
        i_c( n) = ati
        x_c( n) = mesh%TriGCAaAc( ati,1)
        y_c( n) = mesh%TriGCAaAc( ati,2)
      END DO
      
      ! If not enough valid source points are found, a gradient cannot be calculated here
      IF (n < 3) THEN
!        mask_valid_out( vi) = 0
        CYCLE
      END IF
      
      ! Destination point: the regular vertex
      x = mesh%VAaAc( avi,1)
      y = mesh%VAaAc( avi,2)
      
      ! Calculate local neighbour functions
      CALL calc_neighbour_functions_ls_stag( x, y, n, x_c, y_c, Nfc, Nfxc, Nfyc)
      
      ! Fill into sparse matrices
      DO i = 1, n
        k = k+1
        M_map%index( k) = i_c(  i)
        M_ddx%index( k) = i_c(  i)
        M_ddy%index( k) = i_c(  i)
        M_map%val(   k) = Nfc(  i)
        M_ddx%val(   k) = Nfxc( i)
        M_ddy%val(   k) = Nfyc( i)
      END DO

      ! Finish this row
      M_map%ptr( avi+1 : nrows+1) = k+1
      M_ddx%ptr( avi+1 : nrows+1) = k+1
      M_ddy%ptr( avi+1 : nrows+1) = k+1
      
    END DO ! DO avi = mesh%avi1, mesh%avi2
    CALL sync
    
    ! Number of non-zero elements
    M_map%nnz = k
    M_ddx%nnz = k
    M_ddy%nnz = k
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%avi1, mesh%avi2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%avi1, mesh%avi2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%avi1, mesh%avi2)
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfc )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_bb_ac
  SUBROUTINE calc_matrix_operators_combi_acuv( mesh)
    ! Calculate matrix operators involving the acuv-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:
    INTEGER                                            :: nrows, ncols, nnz_max
    INTEGER                                            :: avi, auvi
    
  ! ac_acu/ac_acv: from the ac-grid to the u/v-parts of the acuv-grid
  ! =================================================================
    
    ncols   = mesh%nVAaAc   ! from
    nrows   = 2*mesh%nVAaAc ! to
    nnz_max = mesh%nVAaAc
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( mesh%M_move_ac_acu, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( mesh%M_move_ac_acv, nrows, ncols, nnz_max)
    
    ! Initialise
    mesh%M_move_ac_acu%nnz = 0
    mesh%M_move_ac_acv%nnz = 0
    mesh%M_move_ac_acu%ptr = 1
    mesh%M_move_ac_acv%ptr = 1
    
    ! Fill the matrices
    DO auvi = mesh%auvi1, mesh%auvi2
    
      IF     (MOD(auvi,2) == 1) THEN
        ! u-field
        avi = (auvi + 1) / 2
        mesh%M_move_ac_acu%nnz  = mesh%M_move_ac_acu%nnz + 1
        mesh%M_move_ac_acu%index( mesh%M_move_ac_acu%nnz) = avi
        mesh%M_move_ac_acu%val(   mesh%M_move_ac_acu%nnz) = 1._dp
        mesh%M_move_ac_acu%ptr( auvi+1 : nrows+1) = mesh%M_move_ac_acu%nnz + 1
      ELSE
        ! v-field
        avi = auvi / 2
        mesh%M_move_ac_acv%nnz  = mesh%M_move_ac_acv%nnz + 1
        mesh%M_move_ac_acv%index( mesh%M_move_ac_acv%nnz) = avi
        mesh%M_move_ac_acv%val(   mesh%M_move_ac_acv%nnz) = 1._dp
        mesh%M_move_ac_acv%ptr( auvi+1 : nrows+1) = mesh%M_move_ac_acv%nnz + 1
      END IF
      
    END DO ! DO auvi = mesh%auvi1, mesh%auvi2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( mesh%M_move_ac_acu, mesh%auvi1, mesh%auvi2)
    CALL finalise_matrix_CSR_dist( mesh%M_move_ac_acv, mesh%auvi1, mesh%auvi2)
    
  ! acau_aca/acav_aca: from the u/v-parts of the acuv-grid to the ac-grid
  ! =====================================================================
    
    ncols   = 2*mesh%nVAaAc ! from
    nrows   = mesh%nVAaAc   ! to
    nnz_max = mesh%nVAaAc
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( mesh%M_move_acu_ac, nrows, ncols, nnz_max)
    CALL allocate_matrix_CSR_dist( mesh%M_move_acv_ac, nrows, ncols, nnz_max)
    
    ! Initialise
    mesh%M_move_acu_ac%nnz = 0
    mesh%M_move_acv_ac%nnz = 0
    mesh%M_move_acu_ac%ptr = 1
    mesh%M_move_acv_ac%ptr = 1
    
    ! Fill the matrices
    DO avi = mesh%avi1, mesh%avi2
    
      ! u-field
      auvi = 2*avi - 1
      mesh%M_move_acu_ac%nnz  = mesh%M_move_acu_ac%nnz + 1
      mesh%M_move_acu_ac%index( mesh%M_move_acu_ac%nnz) = auvi
      mesh%M_move_acu_ac%val(   mesh%M_move_acu_ac%nnz) = 1._dp
      mesh%M_move_acu_ac%ptr( avi+1 : nrows+1) = mesh%M_move_acu_ac%nnz + 1
      
      ! v-field
      auvi = 2*avi
      mesh%M_move_acv_ac%nnz  = mesh%M_move_acv_ac%nnz + 1
      mesh%M_move_acv_ac%index( mesh%M_move_acv_ac%nnz) = auvi
      mesh%M_move_acv_ac%val(   mesh%M_move_acv_ac%nnz) = 1._dp
      mesh%M_move_acv_ac%ptr( avi+1 : nrows+1) = mesh%M_move_acv_ac%nnz + 1
      
    END DO ! DO avi = mesh%avi1, mesh%avi2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( mesh%M_move_acu_ac, mesh%avi1, mesh%avi2)
    CALL finalise_matrix_CSR_dist( mesh%M_move_acv_ac, mesh%avi1, mesh%avi2)
    
  ! acu_bb/acv_bb: from the u/v-parts of the acuv-grid to the bb-grid
  ! =================================================================
    
    ! These can be defined as product of first moving from acu/acv to ac, then mapping/gradient from ac to bb
    
    CALL multiply_matrix_matrix_CSR( mesh%M_map_ac_bb, mesh%M_move_acu_ac, mesh%M_map_acu_bb)
    CALL multiply_matrix_matrix_CSR( mesh%M_map_ac_bb, mesh%M_move_acv_ac, mesh%M_map_acv_bb)
    CALL multiply_matrix_matrix_CSR( mesh%M_ddx_ac_bb, mesh%M_move_acu_ac, mesh%M_ddx_acu_bb)
    CALL multiply_matrix_matrix_CSR( mesh%M_ddx_ac_bb, mesh%M_move_acv_ac, mesh%M_ddx_acv_bb)
    CALL multiply_matrix_matrix_CSR( mesh%M_ddy_ac_bb, mesh%M_move_acu_ac, mesh%M_ddy_acu_bb)
    CALL multiply_matrix_matrix_CSR( mesh%M_ddy_ac_bb, mesh%M_move_acv_ac, mesh%M_ddy_acv_bb)
    
  ! acu_bb/acv_bb: from the bb-grid to the u/v-parts of the ac-grid
  ! ===========================================================================================
    
    ! These can be defined as product of first mapping/gradient from bb to ac, then moving from ac to acu/acv
    
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acu, mesh%M_map_bb_ac, mesh%M_map_bb_acu)
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acv, mesh%M_map_bb_ac, mesh%M_map_bb_acv)
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acu, mesh%M_ddx_bb_ac, mesh%M_ddx_bb_acu)
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acv, mesh%M_ddx_bb_ac, mesh%M_ddx_bb_acv)
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acu, mesh%M_ddy_bb_ac, mesh%M_ddy_bb_acu)
    CALL multiply_matrix_matrix_CSR( mesh%M_move_ac_acv, mesh%M_ddy_bb_ac, mesh%M_ddy_bb_acv)
    
  END SUBROUTINE calc_matrix_operators_combi_acuv
  
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
  
! == Directly apply a Neumann boundary condition to a data field
  SUBROUTINE apply_Neumann_BC_direct_2D( mesh, d_a)
    ! Directly apply a Neumann boundary condition to a data field
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_a
    
    ! Local variables:
    INTEGER                                            :: vi, ci, vc, n
    
    DO vi = mesh%vi1, mesh%vi2
      IF (mesh%edge_index( vi) == 0) CYCLE
      
      d_a( vi) = 0._dp
      n        = 0
      DO ci = 1, mesh%nC( vi)
        vc = mesh%C( vi,ci)
        IF (mesh%edge_index( vc) == 0) CYCLE
        d_a( vi) = d_a( vi) + d_a( vc)
        n = n + 1
      END DO
      
      IF (n > 0) THEN
        d_a( vi) = d_a( vi) / REAL(n,dp)
      ELSE
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          d_a( vi) = d_a( vi) + d_a( vc) / mesh%nC( vi)
        END DO
      END IF
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
  END SUBROUTINE apply_Neumann_BC_direct_2D
  SUBROUTINE apply_Neumann_BC_direct_3D( mesh, d_a)
    ! Directly apply a Neumann boundary condition to a data field
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d_a
    
    ! Local variables:
    INTEGER                                            :: vi, ci, vc, n
    
    DO vi = mesh%vi1, mesh%vi2
      IF (mesh%edge_index( vi) == 0) CYCLE
      
      d_a( vi,:) = 0._dp
      n        = 0
      DO ci = 1, mesh%nC( vi)
        vc = mesh%C( vi,ci)
        IF (mesh%edge_index( vc) == 0) CYCLE
        d_a( vi,:) = d_a( vi,:) + d_a( vc,:)
        n = n + 1
      END DO
      
      IF (n > 0) THEN
        d_a( vi,:) = d_a( vi,:) / REAL(n,dp)
      ELSE
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          d_a( vi,:) = d_a( vi,:) + d_a( vc,:) / mesh%nC( vi)
        END DO
      END IF
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
  END SUBROUTINE apply_Neumann_BC_direct_3D

END MODULE mesh_operators_module
