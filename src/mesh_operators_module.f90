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
  USE data_types_module,               ONLY: type_mesh, type_sparse_matrix_CSR_dp, type_grid
  USE utilities_module,                ONLY: calc_matrix_inverse_2_by_2, calc_matrix_inverse_3_by_3, calc_matrix_inverse_general 
  USE sparse_matrix_module,            ONLY: allocate_matrix_CSR_dist, add_entry_CSR_dist, extend_matrix_CSR_dist, finalise_matrix_CSR_dist, &
                                             multiply_matrix_vector_CSR, multiply_matrix_vector_2D_CSR, check_CSR, check_CSR_dist

  IMPLICIT NONE
  
  LOGICAL, PARAMETER :: do_check_matrices = .TRUE.

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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_a_b, d_a, ddx_b)
    
    ! d2dx2_a = M_ddx_b_a * ddx_b
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_b_a, ddx_b, d2dx2_a)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_a_b, d_a, ddx_b)
    
    ! d2dxdy_a = M_ddy_b_a * ddx_b
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_b_a, ddx_b, d2dxdy_a)
    
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
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_a_b, d_a, ddy_b)
    
    ! d2dy2_a = M_ddy_b_a * ddy_b
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_b_a, ddy_b, d2dy2_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_a_b, d_a, ddx_b)
    
    ! d2dx2_a = M_ddx_b_a * ddx_b
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_b_a, ddx_b, d2dx2_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddx_a_b, d_a, ddx_b)
    
    ! d2dxdy_a = M_ddy_b_a * ddx_b
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_b_a, ddx_b, d2dxdy_a)
    
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
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_a_b, d_a, ddy_b)
    
    ! d2dy2_a = M_ddy_b_a * ddy_b
    CALL multiply_matrix_vector_2D_CSR( mesh%M_ddy_b_a, ddy_b, d2dy2_a)
    
    ! Clean up after yourself
    CALL deallocate_shared( wddy_b)
    
  END SUBROUTINE d2dy2_a_to_a_3D
  
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
    CALL calc_matrix_operators_2nd_order_b_b( mesh)
    
    ! Calculate matrix operator for applying Neumann boundary conditions on border triangles
    CALL calc_matrix_operator_Neumann_BC_b_b( mesh)
    
  END SUBROUTINE calc_matrix_operators_mesh
  
  SUBROUTINE calc_matrix_operators_a_a( mesh,        M_ddx, M_ddy)
    ! Calculate all the matrix operators representing the mapping,
    ! d/dx and d/dy operations from the a (vertex) to the a (vertex) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: vi, n, ci, vj, i
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc

    ncols           = mesh%nV      ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_max = mesh%nC_mem+1
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max_proc)

    DO vi = mesh%vi1, mesh%vi2
      
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
      CALL add_entry_CSR_dist( M_ddx, vi, vi, Nfxi)
      CALL add_entry_CSR_dist( M_ddy, vi, vi, Nfyi)
      DO i = 1, n
        CALL add_entry_CSR_dist( M_ddx, vi, i_c( i), Nfxc( i))
        CALL add_entry_CSR_dist( M_ddy, vi, i_c( i), Nfyc( i))
      END DO
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( M_ddx, 'mesh%M_ddx_a_a', mesh%vi1, mesh%vi2)
      CALL check_CSR_dist( M_ddy, 'mesh%M_ddy_a_a', mesh%vi1, mesh%vi2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%vi1, mesh%vi2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( M_ddx, 'mesh%M_ddx_a_a')
      CALL check_CSR( M_ddy, 'mesh%M_ddy_a_a')
    END IF
    
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
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ti, n, vii, vi, i
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nV      ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 3
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max_proc)

    DO ti = mesh%ti1, mesh%ti2
      
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
        CALL add_entry_CSR_dist( M_map, ti, i_c( i), Nfc(  i))
        CALL add_entry_CSR_dist( M_ddx, ti, i_c( i), Nfxc( i))
        CALL add_entry_CSR_dist( M_ddy, ti, i_c( i), Nfyc( i))
      END DO
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( M_map, 'mesh%M_map_a_b', mesh%ti1, mesh%ti2)
      CALL check_CSR_dist( M_ddx, 'mesh%M_ddx_a_b', mesh%ti1, mesh%ti2)
      CALL check_CSR_dist( M_ddy, 'mesh%M_ddy_a_b', mesh%ti1, mesh%ti2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ti1, mesh%ti2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( M_map, 'mesh%M_map_a_b')
      CALL check_CSR( M_ddx, 'mesh%M_ddx_a_b')
      CALL check_CSR( M_ddy, 'mesh%M_ddy_a_b')
    END IF
    
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
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: aci, n, acii, vi, i
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nV      ! from
    nrows           = mesh%nAc     ! to
    nnz_per_row_max = 4
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max_proc)

    DO aci = mesh%ci1, mesh%ci2
      
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
        CALL add_entry_CSR_dist( M_map, aci, i_c( i), Nfc(  i))
        CALL add_entry_CSR_dist( M_ddx, aci, i_c( i), Nfxc( i))
        CALL add_entry_CSR_dist( M_ddy, aci, i_c( i), Nfyc( i))
      END DO
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( M_map, 'mesh%M_map_a_c', mesh%ci1, mesh%ci2)
      CALL check_CSR_dist( M_ddx, 'mesh%M_ddx_a_c', mesh%ci1, mesh%ci2)
      CALL check_CSR_dist( M_ddy, 'mesh%M_ddy_a_c', mesh%ci1, mesh%ci2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%ci1, mesh%ci2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ci1, mesh%ci2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ci1, mesh%ci2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( M_map, 'mesh%M_map_a_c')
      CALL check_CSR( M_ddx, 'mesh%M_ddx_a_c')
      CALL check_CSR( M_ddy, 'mesh%M_ddy_a_c')
    END IF
    
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
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: vi, n, iti, ti, i
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_max = mesh%nC_mem
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max_proc)

    DO vi = mesh%vi1, mesh%vi2
      
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
        CALL add_entry_CSR_dist( M_map, vi, i_c( i), Nfc(  i))
        CALL add_entry_CSR_dist( M_ddx, vi, i_c( i), Nfxc( i))
        CALL add_entry_CSR_dist( M_ddy, vi, i_c( i), Nfyc( i))
      END DO
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( M_map, 'mesh%M_map_b_a', mesh%vi1, mesh%vi2)
      CALL check_CSR_dist( M_ddx, 'mesh%M_ddx_b_a', mesh%vi1, mesh%vi2)
      CALL check_CSR_dist( M_ddy, 'mesh%M_ddy_b_a', mesh%vi1, mesh%vi2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%vi1, mesh%vi2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( M_map, 'mesh%M_map_b_a')
      CALL check_CSR( M_ddx, 'mesh%M_ddx_b_a')
      CALL check_CSR( M_ddy, 'mesh%M_ddy_b_a')
    END IF
    
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
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ti, n, tii, tj, i
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 4
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    nnz_max = nrows * nnz_per_row_max
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max_proc)

    DO ti = mesh%ti1, mesh%ti2
      
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
      CALL add_entry_CSR_dist( M_ddx, ti, ti, Nfxi)
      CALL add_entry_CSR_dist( M_ddy, ti, ti, Nfyi)
      DO i = 1, n
        CALL add_entry_CSR_dist( M_ddx, ti, i_c( i), Nfxc( i))
        CALL add_entry_CSR_dist( M_ddy, ti, i_c( i), Nfyc( i))
      END DO
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( M_ddx, 'mesh%M_ddx_b_b', mesh%ti1, mesh%ti2)
      CALL check_CSR_dist( M_ddy, 'mesh%M_ddy_b_b', mesh%ti1, mesh%ti2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ti1, mesh%ti2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( M_ddx, 'mesh%M_ddx_b_b')
      CALL check_CSR( M_ddy, 'mesh%M_ddy_b_b')
    END IF
    
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
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: aci, n, acii, vi, iti, ti, n_match, vii, aciii, i
    LOGICAL                                            :: is_listed
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nAc     ! to
    nnz_per_row_max = 6
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max_proc)

    DO aci = mesh%ci1, mesh%ci2
      
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
        CALL add_entry_CSR_dist( M_map, aci, i_c( i), Nfc(  i))
        CALL add_entry_CSR_dist( M_ddx, aci, i_c( i), Nfxc( i))
        CALL add_entry_CSR_dist( M_ddy, aci, i_c( i), Nfyc( i))
      END DO
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( M_map, 'mesh%M_map_b_c', mesh%ci1, mesh%ci2)
      CALL check_CSR_dist( M_ddx, 'mesh%M_ddx_b_c', mesh%ci1, mesh%ci2)
      CALL check_CSR_dist( M_ddy, 'mesh%M_ddy_b_c', mesh%ci1, mesh%ci2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%ci1, mesh%ci2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ci1, mesh%ci2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ci1, mesh%ci2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( M_map, 'mesh%M_map_b_c')
      CALL check_CSR( M_ddx, 'mesh%M_ddx_b_c')
      CALL check_CSR( M_ddy, 'mesh%M_ddy_b_c')
    END IF
    
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
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: vi, n, ci, aci, i
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nAc     ! from
    nrows           = mesh%nV      ! to
    nnz_per_row_max = mesh%nC_mem
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max_proc)

    DO vi = mesh%vi1, mesh%vi2
      
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
        CALL add_entry_CSR_dist( M_map, vi, i_c( i), Nfc(  i))
        CALL add_entry_CSR_dist( M_ddx, vi, i_c( i), Nfxc( i))
        CALL add_entry_CSR_dist( M_ddy, vi, i_c( i), Nfyc( i))
      END DO
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( M_map, 'mesh%M_map_c_a', mesh%vi1, mesh%vi2)
      CALL check_CSR_dist( M_ddx, 'mesh%M_ddx_c_a', mesh%vi1, mesh%vi2)
      CALL check_CSR_dist( M_ddy, 'mesh%M_ddy_c_a', mesh%vi1, mesh%vi2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%vi1, mesh%vi2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( M_map, 'mesh%M_map_c_a')
      CALL check_CSR( M_ddx, 'mesh%M_ddx_c_a')
      CALL check_CSR( M_ddy, 'mesh%M_ddy_c_a')
    END IF
    
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
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: M_map, M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: ti, n, via, vib, vic, ci, vj, acab, acbc, acca, i
    REAL(dp)                                           :: x, y
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfc, Nfxc, Nfyc

    ncols           = mesh%nAc     ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 3
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfc(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    CALL allocate_matrix_CSR_dist( M_map, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max_proc)

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
        CALL add_entry_CSR_dist( M_map, ti, i_c( i), Nfc(  i))
        CALL add_entry_CSR_dist( M_ddx, ti, i_c( i), Nfxc( i))
        CALL add_entry_CSR_dist( M_ddy, ti, i_c( i), Nfyc( i))
      END DO
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( M_map, 'mesh%M_map_c_b', mesh%ti1, mesh%ti2)
      CALL check_CSR_dist( M_ddx, 'mesh%M_ddx_c_b', mesh%ti1, mesh%ti2)
      CALL check_CSR_dist( M_ddy, 'mesh%M_ddy_c_b', mesh%ti1, mesh%ti2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_map, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ti1, mesh%ti2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( M_map, 'mesh%M_map_c_b')
      CALL check_CSR( M_ddx, 'mesh%M_ddx_c_b')
      CALL check_CSR( M_ddy, 'mesh%M_ddy_c_b')
    END IF
    
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
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: i_c
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: x_c, y_c
    INTEGER                                            :: aci, n, acii, vi, ci, vc, acj, i
    REAL(dp)                                           :: x, y
    REAL(dp)                                           :: Nfxi, Nfyi
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: Nfxc, Nfyc

    ncols           = mesh%nAc     ! from
    nrows           = mesh%nAc     ! to
    nnz_per_row_max = 5
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))
    
    ALLOCATE( i_c(  nnz_per_row_max))
    ALLOCATE( x_c(  nnz_per_row_max))
    ALLOCATE( y_c(  nnz_per_row_max))
    ALLOCATE( Nfxc( nnz_per_row_max))
    ALLOCATE( Nfyc( nnz_per_row_max))

    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max_proc)

    DO aci = mesh%ci1, mesh%ci2
      
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
      CALL add_entry_CSR_dist( M_ddx, aci, aci, Nfxi)
      CALL add_entry_CSR_dist( M_ddy, aci, aci, Nfyi)
      DO i = 1, n
        CALL add_entry_CSR_dist( M_ddx, aci, i_c( i), Nfxc( i))
        CALL add_entry_CSR_dist( M_ddy, aci, i_c( i), Nfyc( i))
      END DO
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( M_ddx, 'mesh%M_ddx_c_c', mesh%ci1, mesh%ci2)
      CALL check_CSR_dist( M_ddy, 'mesh%M_ddy_c_c', mesh%ci1, mesh%ci2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_ddx, mesh%ci1, mesh%ci2)
    CALL finalise_matrix_CSR_dist( M_ddy, mesh%ci1, mesh%ci2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( M_ddx, 'mesh%M_ddx_c_c')
      CALL check_CSR( M_ddy, 'mesh%M_ddy_c_c')
    END IF
    
    ! Clean up after yourself
    DEALLOCATE( i_c )
    DEALLOCATE( x_c )
    DEALLOCATE( y_c )
    DEALLOCATE( Nfxc)
    DEALLOCATE( Nfyc)
    
  END SUBROUTINE calc_matrix_operators_c_c
  
  SUBROUTINE calc_matrix_operators_2nd_order_b_b( mesh)
    ! Calculate the matrix operators representing the d/dx, d/dy, d2/dx2,
    ! d2/dxdy, and d2/dy2 operations from the b (triangle) to the b (triangle) grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
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
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))
    
    ALLOCATE( i_c(   nnz_per_row_max))
    ALLOCATE( x_c(   nnz_per_row_max))
    ALLOCATE( y_c(   nnz_per_row_max))
    ALLOCATE( Nfxc(  nnz_per_row_max))
    ALLOCATE( Nfyc(  nnz_per_row_max))
    ALLOCATE( Nfxxc( nnz_per_row_max))
    ALLOCATE( Nfxyc( nnz_per_row_max))
    ALLOCATE( Nfyyc( nnz_per_row_max))

    CALL allocate_matrix_CSR_dist( mesh%M2_ddx_b_b   , nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_ddy_b_b   , nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_d2dx2_b_b , nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_d2dxdy_b_b, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( mesh%M2_d2dy2_b_b , nrows, ncols, nnz_max_proc)

    DO ti = mesh%ti1, mesh%ti2
      
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
      CALL add_entry_CSR_dist( mesh%M2_ddx_b_b   , ti, ti, Nfxi )
      CALL add_entry_CSR_dist( mesh%M2_ddy_b_b   , ti, ti, Nfyi )
      CALL add_entry_CSR_dist( mesh%M2_d2dx2_b_b , ti, ti, Nfxxi)
      CALL add_entry_CSR_dist( mesh%M2_d2dxdy_b_b, ti, ti, Nfxyi)
      CALL add_entry_CSR_dist( mesh%M2_d2dy2_b_b , ti, ti, Nfyyi)
      DO i = 1, n
        CALL add_entry_CSR_dist( mesh%M2_ddx_b_b   , ti, i_c( i), Nfxc(  i))
        CALL add_entry_CSR_dist( mesh%M2_ddy_b_b   , ti, i_c( i), Nfyc(  i))
        CALL add_entry_CSR_dist( mesh%M2_d2dx2_b_b , ti, i_c( i), Nfxxc( i))
        CALL add_entry_CSR_dist( mesh%M2_d2dxdy_b_b, ti, i_c( i), Nfxyc( i))
        CALL add_entry_CSR_dist( mesh%M2_d2dy2_b_b , ti, i_c( i), Nfyyc( i))
      END DO
      
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( mesh%M2_ddx_b_b   , 'mesh%M2_ddx_b_b'   , mesh%ti1, mesh%ti2)
      CALL check_CSR_dist( mesh%M2_ddy_b_b   , 'mesh%M2_ddy_b_b'   , mesh%ti1, mesh%ti2)
      CALL check_CSR_dist( mesh%M2_d2dx2_b_b , 'mesh%M2_d2dx2_b_b' , mesh%ti1, mesh%ti2)
      CALL check_CSR_dist( mesh%M2_d2dxdy_b_b, 'mesh%M2_d2dxdy_b_b', mesh%ti1, mesh%ti2)
      CALL check_CSR_dist( mesh%M2_d2dy2_b_b , 'mesh%M2_d2dy2_b_b' , mesh%ti1, mesh%ti2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( mesh%M2_ddx_b_b   , mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( mesh%M2_ddy_b_b   , mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( mesh%M2_d2dx2_b_b , mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( mesh%M2_d2dxdy_b_b, mesh%ti1, mesh%ti2)
    CALL finalise_matrix_CSR_dist( mesh%M2_d2dy2_b_b , mesh%ti1, mesh%ti2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( mesh%M2_ddx_b_b   , 'mesh%M2_ddx_b_b'   )
      CALL check_CSR( mesh%M2_ddy_b_b   , 'mesh%M2_ddy_b_b'   )
      CALL check_CSR( mesh%M2_d2dx2_b_b , 'mesh%M2_d2dx2_b_b' )
      CALL check_CSR( mesh%M2_d2dxdy_b_b, 'mesh%M2_d2dxdy_b_b')
      CALL check_CSR( mesh%M2_d2dy2_b_b , 'mesh%M2_d2dy2_b_b' )
    END IF
    
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
  
  SUBROUTINE calc_matrix_operator_Neumann_BC_b_b( mesh)
    ! Calculate matrix operator for applying Neumann boundary conditions on border triangles
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER                                            :: ti, n, tti, tj

  ! Then fill in the stiffness matrix coefficients
  ! ==============================================

    ncols           = mesh%nTri    ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 4
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))

    CALL allocate_matrix_CSR_dist( mesh%M_Neumann_BC_b_b, nrows, ncols, nnz_max_proc)
    
    DO ti = mesh%ti1, mesh%ti2
      
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
        CALL add_entry_CSR_dist( mesh%M_Neumann_BC_b_b, ti, ti, -1._dp)
        
        ! Neighbouring triangles
        DO tti = 1, 3
          tj = mesh%TriC( ti,tti)
          IF (tj == 0) CYCLE
          CALL add_entry_CSR_dist( mesh%M_Neumann_BC_b_b, ti, tj, 1._dp / REAL( n,dp))
        END DO
        
      END IF ! IF (is_border_triangle) THEN
    
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( mesh%M_Neumann_BC_b_b, 'mesh%M_Neumann_BC_b_b', mesh%ti1, mesh%ti2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( mesh%M_Neumann_BC_b_b, mesh%ti1, mesh%ti2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( mesh%M_Neumann_BC_b_b, 'mesh%M_Neumann_BC_b_b')
    END IF
    
  END SUBROUTINE calc_matrix_operator_Neumann_BC_b_b
  
  SUBROUTINE calc_matrix_operators_grid( grid, M_ddx, M_ddy)
    ! Calculate matrix operators for partial derivatives on a regular grid (needed for conservative remapping)
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: M_ddx, M_ddy
    
    ! Local variables:
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc
    INTEGER                                            :: n1, n2, n, i, j, n_ip1, n_im1, n_jp1, n_jm1
    REAL(dp)                                           :: v1, v2

    ncols           = grid%n      ! from
    nrows           = grid%n      ! to
    nnz_per_row_max = 3
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL( nnz_max,dp) / REAL( par%n,dp))

    CALL allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nnz_max_proc)
    CALL allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nnz_max_proc)
    
    CALL partition_list( grid%n, par%i, par%n, n1, n2)
    
    DO n = n1, n2
      
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      
      ! d/dx
      IF     (i == 1) THEN
        n_ip1 = grid%ij2n( i+1,j)
        v1 = -1._dp / grid%dx
        v2 =  1._dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx, n, n    , v1)
        CALL add_entry_CSR_dist( M_ddx, n, n_ip1, v2)
      ELSEIF (i == grid%nx) THEN
        n_im1 = grid%ij2n( i-1,j)
        v1 =  1._dp / grid%dx
        v2 = -1._dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx, n, n    , v1)
        CALL add_entry_CSR_dist( M_ddx, n, n_im1, v2)
      ELSE
        n_im1 = grid%ij2n( i-1,j)
        n_ip1 = grid%ij2n( i+1,j)
        v1 = -0.5_dp / grid%dx
        v2 =  0.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddx, n, n_im1, v1)
        CALL add_entry_CSR_dist( M_ddx, n, n_ip1, v2)
      END IF
      j = grid%n2ij( n,2)
      
      ! d/dy
      IF     (j == 1) THEN
        n_jp1 = grid%ij2n( i,j+1)
        v1 = -1._dp / grid%dx
        v2 =  1._dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy, n, n    , v1)
        CALL add_entry_CSR_dist( M_ddy, n, n_jp1, v2)
      ELSEIF (j == grid%ny) THEN
        n_jm1 = grid%ij2n( i,j-1)
        v1 =  1._dp / grid%dx
        v2 = -1._dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy, n, n    , v1)
        CALL add_entry_CSR_dist( M_ddy, n, n_jm1, v2)
      ELSE
        n_jm1 = grid%ij2n( i,j-1)
        n_jp1 = grid%ij2n( i,j+1)
        v1 = -0.5_dp / grid%dx
        v2 =  0.5_dp / grid%dx
        CALL add_entry_CSR_dist( M_ddy, n, n_jm1, v1)
        CALL add_entry_CSR_dist( M_ddy, n, n_jp1, v2)
      END IF
      
    END DO ! DO n = n1, n2
    CALL sync
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR_dist( M_ddx, 'M_grid_ddx', n1, n2)
      CALL check_CSR_dist( M_ddy, 'M_grid_ddy', n1, n2)
    END IF
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( M_ddx, n1, n2)
    CALL finalise_matrix_CSR_dist( M_ddy, n1, n2)
    
    ! Safety
    IF (do_check_matrices) THEN
      CALL check_CSR( M_ddx, 'M_grid_ddx')
      CALL check_CSR( M_ddy, 'M_grid_ddy')
    END IF
    
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
