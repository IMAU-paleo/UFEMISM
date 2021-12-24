MODULE tests_and_checks_module

  ! A bunch of tests and checks of the CSR matrix operators

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
  USE mesh_operators_module
  USE utilities_module
  USE netcdf_module,                   ONLY: write_CSR_matrix_to_NetCDF

  IMPLICIT NONE

CONTAINS

  SUBROUTINE run_all_matrix_tests( mesh)
    ! Test all the matrix operators
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) '========================================'
    IF (par%master) WRITE(0,*) '=== Testing all the matrix operators ==='
    IF (par%master) WRITE(0,*) '========================================'
    IF (par%master) WRITE(0,*) ''
    
  ! Test basic CSR matrix operations
  ! ================================
    
    CALL test_CSR_matrix_operations
    
  ! Test all the matrix operators
  ! =============================
    
    ! Between the a,b,c-grids
    CALL test_matrix_operators_abc_2D(     mesh)
    CALL test_matrix_operators_abc_3D(     mesh)
    
    ! Between the a/c-grids and the ac-grid
    CALL test_matrix_operators_a_c_ac_2D(  mesh)
    CALL test_matrix_operators_a_c_ac_3D(  mesh)
    
    ! On the ac-grid
    CALL test_matrix_operators_ac_ac_2D(   mesh)
    CALL test_matrix_operators_ac_ac_3D(   mesh)
    
    ! Between the ac-grid and the bb-grid
    CALL test_matrix_operators_ac_bb_2D(   mesh)
    CALL test_matrix_operators_ac_bb_3D(   mesh)
    
    ! Between the ac-grid and the acu/acv-grids
    CALL test_matrix_operators_ac_acuv_2D( mesh)
    CALL test_matrix_operators_ac_acuv_3D( mesh)
  
    ! Between the acu/acv-grids and the bb-grid
    CALL test_matrix_operators_acuv_bb_2D( mesh)
    CALL test_matrix_operators_acuv_bb_3D( mesh)
    
  ! Solve the (modified) Laplace equation as the ultimate test
  ! ==========================================================
  
    CALL solve_modified_Laplace_equation( mesh)
    CALL solve_modified_Laplace_equation_combi( mesh)
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) '================================================='
    IF (par%master) WRITE(0,*) '=== Finished testing all the matrix operators ==='
    IF (par%master) WRITE(0,*) '================================================='
    IF (par%master) WRITE(0,*) ''
    
  END SUBROUTINE run_all_matrix_tests
  
! == Test basic CSR matrix operations
  SUBROUTINE test_CSR_matrix_operations
    ! Testing some basic CSR-formatted matrix operations
    
    IMPLICIT NONE
        
    IF (par%master) WRITE(0,*) '  Testing basic CSR matrix operations...'
    
    CALL test_CSR_matrix_operations_add    
    CALL test_CSR_matrix_operations_overwrite_rows
    CALL test_CSR_matrix_operations_multiply_vector
    CALL test_CSR_matrix_operations_multiply_matrix
    CALL test_CSR_matrix_operations_multiply_matrix_rows_with_vector
    
  END SUBROUTINE test_CSR_matrix_operations
  
  SUBROUTINE test_CSR_matrix_operations_add
    ! Testing some basic CSR-formatted matrix operations
    
    USE mesh_operators_module
    USE utilities_module
    
    IMPLICIT NONE
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR)                       :: AAA, BBB, CCC, CCC_ex
    LOGICAL                                            :: are_identical
    
    IF (par%master) WRITE(0,*) '   Testing matrix addition (C = A+B)...'
    
    ! Define some simple matrices
    CALL allocate_matrix_CSR_shared( AAA   , 64, 3, 96)
    CALL allocate_matrix_CSR_shared( BBB   , 64, 3, 96)
    CALL allocate_matrix_CSR_shared( CCC_ex, 64, 3, 144)
    
    IF (par%master) THEN
AAA%nnz = 96
AAA%ptr   = [1,1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,76,79,82,85,88,91,94,97]
AAA%index = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3]
AAA%val   = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
BBB%nnz = 96
BBB%ptr   = [1,1,2,3,4,6,8,10,13,13,14,15,16,18,20,22,25,25,26,27,28,30,32,34,37,37,38,39,40,42,44,46,49,49,50,51,52,54,56,58,61,61,62,63,64,66,68,70,73,73,74,75,76,78,80,82,85,85,86,87,88,90,92,94,97]
BBB%index = [1,2,3,1,2,1,3,2,3,1,2,3,1,2,3,1,2,1,3,2,3,1,2,3,1,2,3,1,2,1,3,2,3,1,2,3,1,2,3,1,2,1,3,2,3,1,2,3,1,2,3,1,2,1,3,2,3,1,2,3,1,2,3,1,2,1,3,2,3,1,2,3,1,2,3,1,2,1,3,2,3,1,2,3,1,2,3,1,2,1,3,2,3,1,2,3]
BBB%val   = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
CCC_ex%nnz = 144
CCC_ex%ptr   = [1,1,2,3,4,6,8,10,13,14,15,17,19,21,23,26,29,30,32,33,35,37,40,42,45,46,48,50,51,54,56,58,61,63,65,67,70,72,75,78,81,83,85,88,90,93,95,98,101,103,106,108,110,113,116,118,121,124,127,130,133,136,139,142,145]
CCC_ex%index = [1,2,3,1,2,1,3,2,3,1,2,3,1,1,1,2,1,3,1,2,1,3,1,2,3,1,2,3,2,1,2,2,2,3,1,2,1,2,3,2,3,1,2,3,3,1,3,2,3,3,1,2,3,1,3,2,3,1,2,3,1,2,1,2,1,2,1,2,3,1,2,1,2,3,1,2,3,1,2,3,1,3,1,3,1,2,3,1,3,1,2,3,1,3,1,2,3,1,2,3,2,3,1,2,3,2,3,2,3,1,2,3,1,2,3,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3]
CCC_ex%val   = [2,2,2,2,2,2,2,2,2,2,2,2,1,3,1,2,1,2,3,2,3,2,1,2,2,3,2,2,1,2,1,3,1,2,2,3,2,1,2,3,2,2,3,2,1,2,1,2,1,3,2,2,1,2,3,2,3,2,2,3,1,1,3,1,1,3,1,1,2,3,3,3,1,2,1,3,2,3,3,2,1,1,3,1,1,2,1,1,3,3,2,1,3,3,1,2,3,3,2,3,1,1,2,1,1,3,1,1,3,2,3,1,2,1,3,3,3,2,3,3,1,1,1,3,1,1,1,3,1,1,1,3,3,3,1,3,1,3,1,3,3,3,3,3]
    END IF
    CALL sync
    
    CALL write_CSR_matrix_to_NetCDF( AAA,'A.nc')
    CALL write_CSR_matrix_to_NetCDF( BBB,'B.nc')
    CALL write_CSR_matrix_to_NetCDF( CCC_ex,'C_ex.nc')
    
    ! Perform multiplication
    CALL add_matrix_matrix_CSR( AAA, BBB, CCC)
    CALL write_CSR_matrix_to_NetCDF( CCC,'C.nc')
    
    ! Check if operation gave the correct result
    CALL are_identical_matrices_CSR( CCC, CCC_ex, are_identical)
    IF (par%master .AND. (.NOT. are_identical)) THEN
      WRITE(0,*) 'test_CSR_matrix_operations_add - ERROR: CSR matrix addition gave wrong answer!'
      WRITE(0,*) '  C%nnz   = ', CCC%nnz
      WRITE(0,*) '  C%ptr   = ', CCC%ptr
      WRITE(0,*) '  C%index = ', CCC%index
      WRITE(0,*) '  C%val   = ', CCC%val
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_matrix_CSR( AAA)
    CALL deallocate_matrix_CSR( BBB)
    CALL deallocate_matrix_CSR( CCC)
    CALL deallocate_matrix_CSR( CCC_ex)
    
  END SUBROUTINE test_CSR_matrix_operations_add
  SUBROUTINE test_CSR_matrix_operations_overwrite_rows
    ! Testing some basic CSR-formatted matrix operations
    
    USE mesh_operators_module
    USE utilities_module
    
    IMPLICIT NONE
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR)                       :: AAA, BBB
    
    IF (par%master) WRITE(0,*) '   Testing matrix row overwriting (C = A (+O) B)...'
    
    ! Define some simple matrices
    CALL allocate_matrix_CSR_shared( AAA, 4, 4, 8)
    CALL allocate_matrix_CSR_shared( BBB, 4, 4, 8)
    
    IF (par%master) THEN
      ! A
      AAA%nnz   = 8
      AAA%ptr   = [1,3,5,6,9]
      AAA%index = [1,3,1,4,2,1,3,4]
      AAA%val   = [1,2,3,4,5,6,7,8]
      ! B
      BBB%nnz   = 4
      BBB%ptr   = [1,3,3,5,5]
      BBB%index = [1,3,2,4]
      BBB%val   = [15,16,17,18]
    END IF
    CALL sync
    
    ! Perform multiplication
    CALL overwrite_rows_CSR( AAA, BBB)
    
    ! Check if operation gave the correct result
    IF (par%master) THEN
      IF (AAA%nnz == 9 .AND. &
          AAA%ptr(    1) ==  1 .AND. &
          AAA%ptr(    2) ==  3 .AND. &
          AAA%ptr(    3) ==  5 .AND. &
          AAA%ptr(    4) ==  7 .AND. &
          AAA%ptr(    5) == 10 .AND. &
          AAA%index(  1) == 1 .AND. &
          AAA%index(  2) == 3 .AND. &
          AAA%index(  3) == 1 .AND. &
          AAA%index(  4) == 4 .AND. &
          AAA%index(  5) == 2 .AND. &
          AAA%index(  6) == 4 .AND. &
          AAA%index(  7) == 1 .AND. &
          AAA%index(  8) == 3 .AND. &
          AAA%index(  9) == 4 .AND. &
          AAA%val(    1) ==  15._dp .AND. &
          AAA%val(    2) ==  16._dp .AND. &
          AAA%val(    3) ==   3._dp .AND. &
          AAA%val(    4) ==   4._dp .AND. &
          AAA%val(    5) ==  17._dp .AND. &
          AAA%val(    6) ==  18._dp .AND. &
          AAA%val(    7) ==   6._dp .AND. &
          AAA%val(    8) ==   7._dp .AND. &
          AAA%val(    9) ==   8._dp) THEN
      ELSE
        WRITE(0,*) 'test_CSR_matrix_operations_overwrite_rows - ERROR: CSR matrix row overwriting gave wrong answer!'
        WRITE(0,*) '  A%nnz   = ', AAA%nnz
        WRITE(0,*) '  A%ptr   = ', AAA%ptr
        WRITE(0,*) '  A%index = ', AAA%index
        WRITE(0,*) '  A%val   = ', AAA%val
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_matrix_CSR( AAA)
    CALL deallocate_matrix_CSR( BBB)
    
  END SUBROUTINE test_CSR_matrix_operations_overwrite_rows
  SUBROUTINE test_CSR_matrix_operations_multiply_vector
    ! Testing some basic CSR-formatted matrix operations
    
    USE mesh_operators_module
    USE utilities_module
    
    IMPLICIT NONE
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR)                       :: AAA
    REAL(dp), DIMENSION(:    ), POINTER                ::  x,  b
    INTEGER                                            :: wx, wb
    
    IF (par%master) WRITE(0,*) '   Testing matrix-vector multiplication (b = A*x)...'
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( 4, x, wx)
    CALL allocate_shared_dp_1D( 4, b, wb)
    
    ! Define some simple matrices
    CALL allocate_matrix_CSR_shared( AAA, 4, 4, 8)
    
    IF (par%master) THEN
      ! A
      AAA%nnz   = 8
      AAA%ptr   = [1,3,5,6,9]
      AAA%index = [1,3,1,4,2,1,3,4]
      AAA%val   = [1,2,3,4,5,6,7,8]
      ! x
      x = [1,2,3,4]
    END IF
    CALL sync
    
    ! Perform multiplication
    CALL multiply_matrix_vector_CSR( AAA, x, b)
    
    ! Check if operation gave the correct result
    IF (par%master) THEN
      IF (b( 1) == 7._dp .AND. b( 2) == 19._dp .AND. b( 3) == 10._dp .AND. b( 4) == 59._dp) THEN
      ELSE
        WRITE(0,*) 'test_CSR_matrix_operations_multiply_vector - ERROR: CSR matrix-vector multiplication gave wrong answer!'
        WRITE(0,*) '  b = ', b
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_matrix_CSR( AAA)
    CALL deallocate_shared( wx)
    CALL deallocate_shared( wb)
    
  END SUBROUTINE test_CSR_matrix_operations_multiply_vector
  SUBROUTINE test_CSR_matrix_operations_multiply_matrix
    ! Testing some basic CSR-formatted matrix operations
    
    USE mesh_operators_module
    USE utilities_module
    
    IMPLICIT NONE
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR)                       :: AAA, BBB, CCC
    
    IF (par%master) WRITE(0,*) '   Testing matrix-matrix multiplication (C = A*B)...'
    
    ! Define some simple matrices
    CALL allocate_matrix_CSR_shared( AAA, 4, 4, 8)
    CALL allocate_matrix_CSR_shared( BBB, 4, 4, 8)
    
    IF (par%master) THEN
      ! A
      AAA%nnz   = 8
      AAA%ptr   = [1,3,5,6,9]
      AAA%index = [1,3,1,4,2,1,3,4]
      AAA%val   = [1,2,3,4,5,6,7,8]
      ! B
      BBB%nnz   = 8
      BBB%ptr   = [1,3,5,7,9]
      BBB%index = [1,3,3,4,2,3,2,4]
      BBB%val   = [9,10,11,12,13,14,15,16]
    END IF
    CALL sync
    
    ! Perform multiplication
    CALL multiply_matrix_matrix_CSR( AAA, BBB, CCC)
    
    ! Check if operation gave the correct result
    IF (par%master) THEN
      IF (CCC%nnz == 13 .AND. &
          CCC%ptr(    1) ==  1 .AND. &
          CCC%ptr(    2) ==  4 .AND. &
          CCC%ptr(    3) ==  8 .AND. &
          CCC%ptr(    4) == 10 .AND. &
          CCC%ptr(    5) == 14 .AND. &
          CCC%index(  1) == 1 .AND. &
          CCC%index(  2) == 2 .AND. &
          CCC%index(  3) == 3 .AND. &
          CCC%index(  4) == 1 .AND. &
          CCC%index(  5) == 2 .AND. &
          CCC%index(  6) == 3 .AND. &
          CCC%index(  7) == 4 .AND. &
          CCC%index(  8) == 3 .AND. &
          CCC%index(  9) == 4 .AND. &
          CCC%index( 10) == 1 .AND. &
          CCC%index( 11) == 2 .AND. &
          CCC%index( 12) == 3 .AND. &
          CCC%index( 13) == 4 .AND. &
          CCC%val(    1) ==   9._dp .AND. &
          CCC%val(    2) ==  26._dp .AND. &
          CCC%val(    3) ==  38._dp .AND. &
          CCC%val(    4) ==  27._dp .AND. &
          CCC%val(    5) ==  60._dp .AND. &
          CCC%val(    6) ==  30._dp .AND. &
          CCC%val(    7) ==  64._dp .AND. &
          CCC%val(    8) ==  55._dp .AND. &
          CCC%val(    9) ==  60._dp .AND. &
          CCC%val(   10) ==  54._dp .AND. &
          CCC%val(   11) == 211._dp .AND. &
          CCC%val(   12) == 158._dp .AND. &
          CCC%val(   13) == 128._dp) THEN
      ELSE
        WRITE(0,*) 'test_CSR_matrix_operations_multiply_matrix - ERROR: CSR matrix-matrix multiplication gave wrong answer!'
        WRITE(0,*) '  C%nnz   = ', CCC%nnz
        WRITE(0,*) '  C%ptr   = ', CCC%ptr
        WRITE(0,*) '  C%index = ', CCC%index
        WRITE(0,*) '  C%val   = ', CCC%val
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_matrix_CSR( AAA)
    CALL deallocate_matrix_CSR( BBB)
    CALL deallocate_matrix_CSR( CCC)
    
  END SUBROUTINE test_CSR_matrix_operations_multiply_matrix
  SUBROUTINE test_CSR_matrix_operations_multiply_matrix_rows_with_vector
    ! Testing some basic CSR-formatted matrix operations
    
    USE mesh_operators_module
    USE utilities_module
    
    IMPLICIT NONE
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR)                       :: AAA, CCC
    REAL(dp), DIMENSION(:    ), POINTER                ::  BB
    INTEGER                                            :: wBB
    
    IF (par%master) WRITE(0,*) '   Testing multiplication of matrix rows with vector elements (C = DIAG(B) * A)...'
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( 4, BB, wBB)
    
    ! Define some simple matrices
    CALL allocate_matrix_CSR_shared( AAA, 4, 4, 8)
    
    IF (par%master) THEN
      ! A
      AAA%nnz   = 8
      AAA%ptr   = [1,3,5,6,9]
      AAA%index = [1,3,1,4,2,1,3,4]
      AAA%val   = [1,2,3,4,5,6,7,8]
      ! BB
      BB = [1,2,3,4]
    END IF
    CALL sync
    
    ! Perform multiplication
    CALL multiply_matrix_rows_with_vector( AAA, BB, CCC)
    
    ! Check if operation gave the correct result
    IF (par%master) THEN
      IF (CCC%nnz == 8 .AND. &
          CCC%ptr(    1) ==  1 .AND. &
          CCC%ptr(    2) ==  3 .AND. &
          CCC%ptr(    3) ==  5 .AND. &
          CCC%ptr(    4) ==  6 .AND. &
          CCC%ptr(    5) ==  9 .AND. &
          CCC%index(  1) == 1 .AND. &
          CCC%index(  2) == 3 .AND. &
          CCC%index(  3) == 1 .AND. &
          CCC%index(  4) == 4 .AND. &
          CCC%index(  5) == 2 .AND. &
          CCC%index(  6) == 1 .AND. &
          CCC%index(  7) == 3 .AND. &
          CCC%index(  8) == 4 .AND. &
          CCC%val(    1) ==   1._dp .AND. &
          CCC%val(    2) ==   2._dp .AND. &
          CCC%val(    3) ==   6._dp .AND. &
          CCC%val(    4) ==   8._dp .AND. &
          CCC%val(    5) ==  15._dp .AND. &
          CCC%val(    6) ==  24._dp .AND. &
          CCC%val(    7) ==  28._dp .AND. &
          CCC%val(    8) ==  32._dp) THEN
      ELSE
        WRITE(0,*) 'test_CSR_matrix_operations_multiply_matrix - ERROR: CSR matrix rows with vector elements multiplication gave wrong answer!'
        WRITE(0,*) '  C%nnz   = ', CCC%nnz
        WRITE(0,*) '  C%ptr   = ', CCC%ptr
        WRITE(0,*) '  C%index = ', CCC%index
        WRITE(0,*) '  C%val   = ', CCC%val
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_matrix_CSR( AAA)
    CALL deallocate_shared(     wBB)
    CALL deallocate_matrix_CSR( CCC)
    
  END SUBROUTINE test_CSR_matrix_operations_multiply_matrix_rows_with_vector
  
! == Test all the matrix operators (mapping+gradients)

  ! Between the a,b,c-grids
  SUBROUTINE test_matrix_operators_abc_2D( mesh)
    ! Test all the matrix operators on the three basic grids (a,b,c)
    ! 
    ! map_a_to_b_2D
    ! map_a_to_c_2D
    ! map_b_to_a_2D
    ! map_b_to_c_2D
    ! map_c_to_a_2D
    ! map_c_to_b_2D
    ! 
    ! ddx_a_to_a_2D
    ! ddx_a_to_b_2D
    ! ddx_a_to_c_2D
    ! ddx_b_to_a_2D
    ! ddx_b_to_b_2D
    ! ddx_b_to_c_2D
    ! ddx_c_to_a_2D
    ! ddx_c_to_b_2D
    ! ddx_c_to_c_2D
    ! 
    ! ddy_a_to_a_2D
    ! ddy_a_to_b_2D
    ! ddy_a_to_c_2D
    ! ddy_b_to_a_2D
    ! ddy_b_to_b_2D
    ! ddy_b_to_c_2D
    ! ddy_c_to_a_2D
    ! ddy_c_to_b_2D
    ! ddy_c_to_c_2D
    ! 
    ! d2dx2_a_to_a_2D
    ! d2dxdy_a_to_a_2D
    ! d2dy2_a_to_a_2D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, ti, aci
    REAL(dp)                                           :: x, y, d, ddx, ddy, d2dx2, d2dxdy, d2dy2
    
    REAL(dp), DIMENSION(:    ), POINTER                :: d_a_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_c_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_a_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_c_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_a_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_c_ex
    INTEGER :: wd_a_ex, wd_b_ex, wd_c_ex, wddx_a_ex, wddx_b_ex, wddx_c_ex, wddy_a_ex, wddy_b_ex, wddy_c_ex
    
    REAL(dp), DIMENSION(:    ), POINTER                :: d_a_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d_a_to_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d_a_to_c
    REAL(dp), DIMENSION(:    ), POINTER                :: d_b_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d_b_to_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d_b_to_c
    REAL(dp), DIMENSION(:    ), POINTER                :: d_c_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d_c_to_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d_c_to_c
    INTEGER :: wd_a_to_a, wd_a_to_b, wd_a_to_c, wd_b_to_a, wd_b_to_b, wd_b_to_c, wd_c_to_a, wd_c_to_b, wd_c_to_c
    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_a_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_a_to_b
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_a_to_c
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_b_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_b_to_b
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_b_to_c
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_c_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_c_to_b
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_c_to_c
    INTEGER :: wddx_a_to_a, wddx_a_to_b, wddx_a_to_c, wddx_b_to_a, wddx_b_to_b, wddx_b_to_c, wddx_c_to_a, wddx_c_to_b, wddx_c_to_c
    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_a_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_a_to_b
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_a_to_c
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_b_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_b_to_b
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_b_to_c
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_c_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_c_to_b
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_c_to_c
    INTEGER :: wddy_a_to_a, wddy_a_to_b, wddy_a_to_c, wddy_b_to_a, wddy_b_to_b, wddy_b_to_c, wddy_c_to_a, wddy_c_to_b, wddy_c_to_c
    
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dx2_a_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dxdy_a_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dy2_a_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dx2_a_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dxdy_a_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dy2_a_to_a
    INTEGER :: wd2dx2_a_ex, wd2dxdy_a_ex, wd2dy2_a_ex, wd2dx2_a_to_a, wd2dxdy_a_to_a, wd2dy2_a_to_a
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (a/b/c - a/b/c, 2D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
    
    CALL allocate_shared_dp_1D( mesh%nV,   d_a_ex,          wd_a_ex         )
    CALL allocate_shared_dp_1D( mesh%nTri, d_b_ex,          wd_b_ex         )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_c_ex,          wd_c_ex         )
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_a_ex,        wddx_a_ex       )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_b_ex,        wddx_b_ex       )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_c_ex,        wddx_c_ex       )
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_a_ex,        wddy_a_ex       )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_b_ex,        wddy_b_ex       )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_c_ex,        wddy_c_ex       )
    
    CALL allocate_shared_dp_1D( mesh%nV,   d_a_to_a,        wd_a_to_a       )
    CALL allocate_shared_dp_1D( mesh%nTri, d_a_to_b,        wd_a_to_b       )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_a_to_c,        wd_a_to_c       )
    CALL allocate_shared_dp_1D( mesh%nV,   d_b_to_a,        wd_b_to_a       )
    CALL allocate_shared_dp_1D( mesh%nTri, d_b_to_b,        wd_b_to_b       )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_b_to_c,        wd_b_to_c       )
    CALL allocate_shared_dp_1D( mesh%nV,   d_c_to_a,        wd_c_to_a       )
    CALL allocate_shared_dp_1D( mesh%nTri, d_c_to_b,        wd_c_to_b       )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_c_to_c,        wd_c_to_c       )
    
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_a_to_a,      wddx_a_to_a     )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_a_to_b,      wddx_a_to_b     )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_a_to_c,      wddx_a_to_c     )
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_b_to_a,      wddx_b_to_a     )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_b_to_b,      wddx_b_to_b     )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_b_to_c,      wddx_b_to_c     )
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_c_to_a,      wddx_c_to_a     )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_c_to_b,      wddx_c_to_b     )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_c_to_c,      wddx_c_to_c     )
    
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_a_to_a,      wddy_a_to_a     )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_a_to_b,      wddy_a_to_b     )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_a_to_c,      wddy_a_to_c     )
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_b_to_a,      wddy_b_to_a     )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_b_to_b,      wddy_b_to_b     )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_b_to_c,      wddy_b_to_c     )
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_c_to_a,      wddy_c_to_a     )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_c_to_b,      wddy_c_to_b     )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_c_to_c,      wddy_c_to_c     )
    
    CALL allocate_shared_dp_1D( mesh%nV,   d2dx2_a_ex,      wd2dx2_a_ex     )
    CALL allocate_shared_dp_1D( mesh%nV,   d2dxdy_a_ex,     wd2dxdy_a_ex    )
    CALL allocate_shared_dp_1D( mesh%nV,   d2dy2_a_ex,      wd2dy2_a_ex     )
    CALL allocate_shared_dp_1D( mesh%nV,   d2dx2_a_to_a,    wd2dx2_a_to_a   )
    CALL allocate_shared_dp_1D( mesh%nV,   d2dxdy_a_to_a,   wd2dxdy_a_to_a  )
    CALL allocate_shared_dp_1D( mesh%nV,   d2dy2_a_to_a,    wd2dy2_a_to_a   )
    
  ! Calculate exact solutions
  ! =========================
    
    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_a_ex(      vi) = d
      ddx_a_ex(    vi) = ddx
      ddy_a_ex(    vi) = ddy
      d2dx2_a_ex(  vi) = d2dx2
      d2dxdy_a_ex( vi) = d2dxdy
      d2dy2_a_ex(  vi) = d2dy2
    END DO
    DO ti = mesh%ti1, mesh%ti2
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_b_ex(   ti) = d
      ddx_b_ex( ti) = ddx
      ddy_b_ex( ti) = ddy
    END DO
    DO aci = mesh%ci1, mesh%ci2
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_c_ex(   aci) = d
      ddx_c_ex( aci) = ddx
      ddy_c_ex( aci) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
  
    ! Mapping
    d_a_to_a( mesh%vi1:mesh%vi2) = d_a_ex( mesh%vi1:mesh%vi2)
    CALL map_a_to_b_2D( mesh, d_a_ex, d_a_to_b)
    CALL map_a_to_c_2D( mesh, d_a_ex, d_a_to_c)
    
    CALL map_b_to_a_2D( mesh, d_b_ex, d_b_to_a)
    d_b_to_b( mesh%ti1:mesh%ti2) = d_b_ex( mesh%ti1:mesh%ti2)
    CALL map_b_to_c_2D( mesh, d_b_ex, d_b_to_c)
    
    CALL map_c_to_a_2D( mesh, d_c_ex, d_c_to_a)
    CALL map_c_to_b_2D( mesh, d_c_ex, d_c_to_b)
    d_c_to_c( mesh%ci1:mesh%ci2) = d_c_ex( mesh%ci1:mesh%ci2)
    
    ! d/dx
    CALL ddx_a_to_a_2D( mesh, d_a_ex, ddx_a_to_a)
    CALL ddx_a_to_b_2D( mesh, d_a_ex, ddx_a_to_b)
    CALL ddx_a_to_c_2D( mesh, d_a_ex, ddx_a_to_c)
    
    CALL ddx_b_to_a_2D( mesh, d_b_ex, ddx_b_to_a)
    CALL ddx_b_to_b_2D( mesh, d_b_ex, ddx_b_to_b)
    CALL ddx_b_to_c_2D( mesh, d_b_ex, ddx_b_to_c)
    
    CALL ddx_c_to_a_2D( mesh, d_c_ex, ddx_c_to_a)
    CALL ddx_c_to_b_2D( mesh, d_c_ex, ddx_c_to_b)
    CALL ddx_c_to_c_2D( mesh, d_c_ex, ddx_c_to_c)
    
    ! d/dy
    CALL ddy_a_to_a_2D( mesh, d_a_ex, ddy_a_to_a)
    CALL ddy_a_to_b_2D( mesh, d_a_ex, ddy_a_to_b)
    CALL ddy_a_to_c_2D( mesh, d_a_ex, ddy_a_to_c)
    
    CALL ddy_b_to_a_2D( mesh, d_b_ex, ddy_b_to_a)
    CALL ddy_b_to_b_2D( mesh, d_b_ex, ddy_b_to_b)
    CALL ddy_b_to_c_2D( mesh, d_b_ex, ddy_b_to_c)
    
    CALL ddy_c_to_a_2D( mesh, d_c_ex, ddy_c_to_a)
    CALL ddy_c_to_b_2D( mesh, d_c_ex, ddy_c_to_b)
    CALL ddy_c_to_c_2D( mesh, d_c_ex, ddy_c_to_c)
    
    ! second partial derivatives
    CALL d2dx2_a_to_a_2D(  mesh, d_a_ex, d2dx2_a_to_a )
    CALL d2dxdy_a_to_a_2D( mesh, d_a_ex, d2dxdy_a_to_a)
    CALL d2dy2_a_to_a_2D(  mesh, d_a_ex, d2dy2_a_to_a )
    
  ! Write results to debug file
  ! ===========================
  
!    IF (par%master) THEN
!      
!      ! Mapping
!      debug%dp_2D_a_01 = d_a_to_a
!      debug%dp_2D_a_02 = d_b_to_a
!      debug%dp_2D_a_03 = d_c_to_a
!      debug%dp_2D_b_01 = d_a_to_b
!      debug%dp_2D_b_02 = d_b_to_b
!      debug%dp_2D_b_03 = d_c_to_b
!      debug%dp_2D_c_01 = d_a_to_c
!      debug%dp_2D_c_02 = d_b_to_c
!      debug%dp_2D_c_03 = d_c_to_c
!      
!      ! d/dx
!      debug%dp_2D_a_04 = ddx_a_to_a
!      debug%dp_2D_a_05 = ddx_b_to_a
!      debug%dp_2D_a_06 = ddx_c_to_a
!      debug%dp_2D_b_04 = ddx_a_to_b
!      debug%dp_2D_b_05 = ddx_b_to_b
!      debug%dp_2D_b_06 = ddx_c_to_b
!      debug%dp_2D_c_04 = ddx_a_to_c
!      debug%dp_2D_c_05 = ddx_b_to_c
!      debug%dp_2D_c_06 = ddx_c_to_c
!      
!      ! d/dy
!      debug%dp_2D_a_07 = ddy_a_to_a
!      debug%dp_2D_a_08 = ddy_b_to_a
!      debug%dp_2D_a_09 = ddy_c_to_a
!      debug%dp_2D_b_07 = ddy_a_to_b
!      debug%dp_2D_b_08 = ddy_b_to_b
!      debug%dp_2D_b_09 = ddy_c_to_b
!      debug%dp_2D_c_07 = ddy_a_to_c
!      debug%dp_2D_c_08 = ddy_b_to_c
!      debug%dp_2D_c_09 = ddy_c_to_c
!      
!      CALL write_to_debug_file
!      
!    END IF
!    CALL sync
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_a_ex      )
    CALL deallocate_shared( wd_b_ex      )
    CALL deallocate_shared( wd_c_ex      )
    CALL deallocate_shared( wddx_a_ex    )
    CALL deallocate_shared( wddx_b_ex    )
    CALL deallocate_shared( wddx_c_ex    )
    CALL deallocate_shared( wddy_a_ex    )
    CALL deallocate_shared( wddy_b_ex    )
    CALL deallocate_shared( wddy_c_ex    )
    
    CALL deallocate_shared( wd_a_to_a    )
    CALL deallocate_shared( wd_a_to_b    )
    CALL deallocate_shared( wd_a_to_c    )
    CALL deallocate_shared( wd_b_to_a    )
    CALL deallocate_shared( wd_b_to_b    )
    CALL deallocate_shared( wd_b_to_c    )
    CALL deallocate_shared( wd_c_to_a    )
    CALL deallocate_shared( wd_c_to_b    )
    CALL deallocate_shared( wd_c_to_c    )
    
    CALL deallocate_shared( wddx_a_to_a  )
    CALL deallocate_shared( wddx_a_to_b  )
    CALL deallocate_shared( wddx_a_to_c  )
    CALL deallocate_shared( wddx_b_to_a  )
    CALL deallocate_shared( wddx_b_to_b  )
    CALL deallocate_shared( wddx_b_to_c  )
    CALL deallocate_shared( wddx_c_to_a  )
    CALL deallocate_shared( wddx_c_to_b  )
    CALL deallocate_shared( wddx_c_to_c  )
    
    CALL deallocate_shared( wddy_a_to_a  )
    CALL deallocate_shared( wddy_a_to_b  )
    CALL deallocate_shared( wddy_a_to_c  )
    CALL deallocate_shared( wddy_b_to_a  )
    CALL deallocate_shared( wddy_b_to_b  )
    CALL deallocate_shared( wddy_b_to_c  )
    CALL deallocate_shared( wddy_c_to_a  )
    CALL deallocate_shared( wddy_c_to_b  )
    CALL deallocate_shared( wddy_c_to_c  )
    
    CALL deallocate_shared( wd2dx2_a_ex     )
    CALL deallocate_shared( wd2dxdy_a_ex    )
    CALL deallocate_shared( wd2dy2_a_ex     )
    CALL deallocate_shared( wd2dx2_a_to_a   )
    CALL deallocate_shared( wd2dxdy_a_to_a  )
    CALL deallocate_shared( wd2dy2_a_to_a   )
    
  END SUBROUTINE test_matrix_operators_abc_2D
  SUBROUTINE test_matrix_operators_abc_3D( mesh)
    ! Test all the matrix operators on the three basic grids (a,b,c)
    ! 
    ! map_a_to_b_3D
    ! map_a_to_c_3D
    ! map_b_to_a_3D
    ! map_b_to_c_3D
    ! map_c_to_a_3D
    ! map_c_to_b_3D
    ! 
    ! ddx_a_to_a_3D
    ! ddx_a_to_b_3D
    ! ddx_a_to_c_3D
    ! ddx_b_to_a_3D
    ! ddx_b_to_b_3D
    ! ddx_b_to_c_3D
    ! ddx_c_to_a_3D
    ! ddx_c_to_b_3D
    ! ddx_c_to_c_3D
    ! 
    ! ddy_a_to_a_3D
    ! ddy_a_to_b_3D
    ! ddy_a_to_c_3D
    ! ddy_b_to_a_3D
    ! ddy_b_to_b_3D
    ! ddy_b_to_c_3D
    ! ddy_c_to_a_3D
    ! ddy_c_to_b_3D
    ! ddy_c_to_c_3D
    ! 
    ! d3Dx2_a_to_a_3D
    ! d3Dxdy_a_to_a_3D
    ! d3Dy2_a_to_a_3D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, ti, aci
    REAL(dp)                                           :: x, y, d, ddx, ddy, d2dx2, d2dxdy, d2dy2
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_a_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_b_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_c_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_a_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_b_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_c_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_a_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_b_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_c_ex
    INTEGER :: wd_a_ex, wd_b_ex, wd_c_ex, wddx_a_ex, wddx_b_ex, wddx_c_ex, wddy_a_ex, wddy_b_ex, wddy_c_ex
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_a_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_a_to_b
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_a_to_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_b_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_b_to_b
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_b_to_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_c_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_c_to_b
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_c_to_c
    INTEGER :: wd_a_to_a, wd_a_to_b, wd_a_to_c, wd_b_to_a, wd_b_to_b, wd_b_to_c, wd_c_to_a, wd_c_to_b, wd_c_to_c
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_a_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_a_to_b
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_a_to_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_b_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_b_to_b
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_b_to_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_c_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_c_to_b
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_c_to_c
    INTEGER :: wddx_a_to_a, wddx_a_to_b, wddx_a_to_c, wddx_b_to_a, wddx_b_to_b, wddx_b_to_c, wddx_c_to_a, wddx_c_to_b, wddx_c_to_c
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_a_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_a_to_b
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_a_to_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_b_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_b_to_b
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_b_to_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_c_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_c_to_b
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_c_to_c
    INTEGER :: wddy_a_to_a, wddy_a_to_b, wddy_a_to_c, wddy_b_to_a, wddy_b_to_b, wddy_b_to_c, wddy_c_to_a, wddy_c_to_b, wddy_c_to_c
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d2dx2_a_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d2dxdy_a_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d2dy2_a_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d2dx2_a_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d2dxdy_a_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d2dy2_a_to_a
    INTEGER :: wd2dx2_a_ex, wd2dxdy_a_ex, wd2dy2_a_ex, wd2dx2_a_to_a, wd2dxdy_a_to_a, wd2dy2_a_to_a
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (a/b/c - a/b/c, 3D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
    
    ! Exact solutions
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, d_a_ex,          wd_a_ex         )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, d_b_ex,          wd_b_ex         )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, d_c_ex,          wd_c_ex         )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, ddx_a_ex,        wddx_a_ex       )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, ddx_b_ex,        wddx_b_ex       )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, ddx_c_ex,        wddx_c_ex       )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, ddy_a_ex,        wddy_a_ex       )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, ddy_b_ex,        wddy_b_ex       )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, ddy_c_ex,        wddy_c_ex       )
    
    ! Discretised approximations
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, d_a_to_a,        wd_a_to_a       )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, d_a_to_b,        wd_a_to_b       )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, d_a_to_c,        wd_a_to_c       )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, d_b_to_a,        wd_b_to_a       )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, d_b_to_b,        wd_b_to_b       )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, d_b_to_c,        wd_b_to_c       )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, d_c_to_a,        wd_c_to_a       )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, d_c_to_b,        wd_c_to_b       )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, d_c_to_c,        wd_c_to_c       )
    
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, ddx_a_to_a,      wddx_a_to_a     )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, ddx_a_to_b,      wddx_a_to_b     )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, ddx_a_to_c,      wddx_a_to_c     )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, ddx_b_to_a,      wddx_b_to_a     )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, ddx_b_to_b,      wddx_b_to_b     )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, ddx_b_to_c,      wddx_b_to_c     )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, ddx_c_to_a,      wddx_c_to_a     )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, ddx_c_to_b,      wddx_c_to_b     )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, ddx_c_to_c,      wddx_c_to_c     )
    
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, ddy_a_to_a,      wddy_a_to_a     )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, ddy_a_to_b,      wddy_a_to_b     )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, ddy_a_to_c,      wddy_a_to_c     )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, ddy_b_to_a,      wddy_b_to_a     )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, ddy_b_to_b,      wddy_b_to_b     )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, ddy_b_to_c,      wddy_b_to_c     )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, ddy_c_to_a,      wddy_c_to_a     )
    CALL allocate_shared_dp_2D( mesh%nTri, C%nz, ddy_c_to_b,      wddy_c_to_b     )
    CALL allocate_shared_dp_2D( mesh%nAc,  C%nz, ddy_c_to_c,      wddy_c_to_c     )
    
    ! Second partial derivatives
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, d2dx2_a_ex,      wd2dx2_a_ex     )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, d2dxdy_a_ex,     wd2dxdy_a_ex    )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, d2dy2_a_ex,      wd2dy2_a_ex     )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, d2dx2_a_to_a,    wd2dx2_a_to_a   )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, d2dxdy_a_to_a,   wd2dxdy_a_to_a  )
    CALL allocate_shared_dp_2D( mesh%nV,   C%nz, d2dy2_a_to_a,    wd2dy2_a_to_a   )
    
  ! Calculate exact solutions
  ! =========================
    
    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_a_ex(      vi,:) = d
      ddx_a_ex(    vi,:) = ddx
      ddy_a_ex(    vi,:) = ddy
      d2dx2_a_ex(  vi,:) = d2dx2
      d2dxdy_a_ex( vi,:) = d2dxdy
      d2dy2_a_ex(  vi,:) = d2dy2
    END DO
    DO ti = mesh%ti1, mesh%ti2
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_b_ex(   ti,:) = d
      ddx_b_ex( ti,:) = ddx
      ddy_b_ex( ti,:) = ddy
    END DO
    DO aci = mesh%ci1, mesh%ci2
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_c_ex(   aci,:) = d
      ddx_c_ex( aci,:) = ddx
      ddy_c_ex( aci,:) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
  
    ! Mapping
    d_a_to_a( mesh%vi1:mesh%vi2,:) = d_a_ex( mesh%vi1:mesh%vi2,:)
    CALL map_a_to_b_3D( mesh, d_a_ex, d_a_to_b)
    CALL map_a_to_c_3D( mesh, d_a_ex, d_a_to_c)
    
    CALL map_b_to_a_3D( mesh, d_b_ex, d_b_to_a)
    d_b_to_b( mesh%ti1:mesh%ti2,:) = d_b_ex( mesh%ti1:mesh%ti2,:)
    CALL map_b_to_c_3D( mesh, d_b_ex, d_b_to_c)
    
    CALL map_c_to_a_3D( mesh, d_c_ex, d_c_to_a)
    CALL map_c_to_b_3D( mesh, d_c_ex, d_c_to_b)
    d_c_to_c( mesh%ci1:mesh%ci2,:) = d_c_ex( mesh%ci1:mesh%ci2,:)
    
    ! d/dx
    CALL ddx_a_to_a_3D( mesh, d_a_ex, ddx_a_to_a)
    CALL ddx_a_to_b_3D( mesh, d_a_ex, ddx_a_to_b)
    CALL ddx_a_to_c_3D( mesh, d_a_ex, ddx_a_to_c)
    
    CALL ddx_b_to_a_3D( mesh, d_b_ex, ddx_b_to_a)
    CALL ddx_b_to_b_3D( mesh, d_b_ex, ddx_b_to_b)
    CALL ddx_b_to_c_3D( mesh, d_b_ex, ddx_b_to_c)
    
    CALL ddx_c_to_a_3D( mesh, d_c_ex, ddx_c_to_a)
    CALL ddx_c_to_b_3D( mesh, d_c_ex, ddx_c_to_b)
    CALL ddx_c_to_c_3D( mesh, d_c_ex, ddx_c_to_c)
    
    ! d/dy
    CALL ddy_a_to_a_3D( mesh, d_a_ex, ddy_a_to_a)
    CALL ddy_a_to_b_3D( mesh, d_a_ex, ddy_a_to_b)
    CALL ddy_a_to_c_3D( mesh, d_a_ex, ddy_a_to_c)
    
    CALL ddy_b_to_a_3D( mesh, d_b_ex, ddy_b_to_a)
    CALL ddy_b_to_b_3D( mesh, d_b_ex, ddy_b_to_b)
    CALL ddy_b_to_c_3D( mesh, d_b_ex, ddy_b_to_c)
    
    CALL ddy_c_to_a_3D( mesh, d_c_ex, ddy_c_to_a)
    CALL ddy_c_to_b_3D( mesh, d_c_ex, ddy_c_to_b)
    CALL ddy_c_to_c_3D( mesh, d_c_ex, ddy_c_to_c)
    
    ! second partial derivatives
    CALL d2dx2_a_to_a_3D(  mesh, d_a_ex, d2dx2_a_to_a )
    CALL d2dxdy_a_to_a_3D( mesh, d_a_ex, d2dxdy_a_to_a)
    CALL d2dy2_a_to_a_3D(  mesh, d_a_ex, d2dy2_a_to_a )
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_a_ex      )
    CALL deallocate_shared( wd_b_ex      )
    CALL deallocate_shared( wd_c_ex      )
    CALL deallocate_shared( wddx_a_ex    )
    CALL deallocate_shared( wddx_b_ex    )
    CALL deallocate_shared( wddx_c_ex    )
    CALL deallocate_shared( wddy_a_ex    )
    CALL deallocate_shared( wddy_b_ex    )
    CALL deallocate_shared( wddy_c_ex    )
    
    CALL deallocate_shared( wd_a_to_a    )
    CALL deallocate_shared( wd_a_to_b    )
    CALL deallocate_shared( wd_a_to_c    )
    CALL deallocate_shared( wd_b_to_a    )
    CALL deallocate_shared( wd_b_to_b    )
    CALL deallocate_shared( wd_b_to_c    )
    CALL deallocate_shared( wd_c_to_a    )
    CALL deallocate_shared( wd_c_to_b    )
    CALL deallocate_shared( wd_c_to_c    )
    
    CALL deallocate_shared( wddx_a_to_a  )
    CALL deallocate_shared( wddx_a_to_b  )
    CALL deallocate_shared( wddx_a_to_c  )
    CALL deallocate_shared( wddx_b_to_a  )
    CALL deallocate_shared( wddx_b_to_b  )
    CALL deallocate_shared( wddx_b_to_c  )
    CALL deallocate_shared( wddx_c_to_a  )
    CALL deallocate_shared( wddx_c_to_b  )
    CALL deallocate_shared( wddx_c_to_c  )
    
    CALL deallocate_shared( wddy_a_to_a  )
    CALL deallocate_shared( wddy_a_to_b  )
    CALL deallocate_shared( wddy_a_to_c  )
    CALL deallocate_shared( wddy_b_to_a  )
    CALL deallocate_shared( wddy_b_to_b  )
    CALL deallocate_shared( wddy_b_to_c  )
    CALL deallocate_shared( wddy_c_to_a  )
    CALL deallocate_shared( wddy_c_to_b  )
    CALL deallocate_shared( wddy_c_to_c  )
    
    CALL deallocate_shared( wd2dx2_a_ex     )
    CALL deallocate_shared( wd2dxdy_a_ex    )
    CALL deallocate_shared( wd2dy2_a_ex     )
    CALL deallocate_shared( wd2dx2_a_to_a   )
    CALL deallocate_shared( wd2dxdy_a_to_a  )
    CALL deallocate_shared( wd2dy2_a_to_a   )
    
  END SUBROUTINE test_matrix_operators_abc_3D
  
  ! Between the a/c-grids and the ac-grid
  SUBROUTINE test_matrix_operators_a_c_ac_2D( mesh)
    ! Test all the matrix operators between the a/c-grids and the ac-grid
    !
    ! move_a_to_aca_2D
    ! move_c_to_acc_2D
    ! move a_and_c_to_ac_2D
    ! map_a_to_ac_2D
    ! map_c_to_ac_2D
    ! ddx_a_to_ac_2D
    ! ddy_a_to_ac_2D
    ! d_aca_to_a_2D
    ! d_acc_to_c_2D
    ! d_ac_to_a_and_c_2D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, aci, avi
    REAL(dp)                                           :: x, y, d, ddx, ddy, d2dx2, d2dxdy, d2dy2
    
    REAL(dp), DIMENSION(:    ), POINTER                :: d_a_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_c_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ac_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ac_ex
    INTEGER :: wd_a_ex, wd_c_ex, wd_ac_ex, wddx_ac_ex, wddy_ac_ex

    REAL(dp), DIMENSION(:    ), POINTER                :: d_a_to_aca
    REAL(dp), DIMENSION(:    ), POINTER                :: d_c_to_acc
    REAL(dp), DIMENSION(:    ), POINTER                :: d_a_and_c_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: d_a_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: d_c_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_a_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_a_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: d_aca_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acc_to_c
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac_to_a
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac_to_c
    INTEGER :: wd_a_to_aca, wd_c_to_acc, wd_a_and_c_to_ac, wd_a_to_ac, wd_c_to_ac, wddx_a_to_ac, wddy_a_to_ac, wd_aca_to_a, wd_acc_to_c, wd_ac_to_a, wd_ac_to_c
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (a/c - ac, 2D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
  
    CALL allocate_shared_dp_1D( mesh%nV    ,   d_a_ex         , wd_a_ex         )
    CALL allocate_shared_dp_1D( mesh%nAc   ,   d_c_ex         , wd_c_ex         )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d_ac_ex        , wd_ac_ex        )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddx_ac_ex      , wddx_ac_ex      )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddy_ac_ex      , wddy_ac_ex      )
    
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d_a_to_aca     , wd_a_to_aca     )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d_c_to_acc     , wd_c_to_acc     )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d_a_and_c_to_ac, wd_a_and_c_to_ac)
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d_a_to_ac      , wd_a_to_ac      )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d_c_to_ac      , wd_c_to_ac      )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddx_a_to_ac    , wddx_a_to_ac    )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddy_a_to_ac    , wddy_a_to_ac    )
    CALL allocate_shared_dp_1D( mesh%nV    ,   d_aca_to_a     , wd_aca_to_a     )
    CALL allocate_shared_dp_1D( mesh%nAc   ,   d_acc_to_c     , wd_acc_to_c     )
    CALL allocate_shared_dp_1D( mesh%nV    ,   d_ac_to_a      , wd_ac_to_a      )
    CALL allocate_shared_dp_1D( mesh%nAc   ,   d_ac_to_c      , wd_ac_to_c      )
    
  ! Calculate exact solutions
  ! =========================
    
    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_a_ex(      vi) = d
    END DO
    DO aci = mesh%ci1, mesh%ci2
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_c_ex(   aci) = d
    END DO
    DO avi = mesh%avi1, mesh%avi2
      x = mesh%VAaAc( avi,1)
      y = mesh%VAaAc( avi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_ac_ex(   avi) = d
      ddx_ac_ex( avi) = ddx
      ddy_ac_ex( avi) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL move_a_to_aca_2D(      mesh, d_a_ex,         d_a_to_aca                )
    CALL move_c_to_acc_2D(      mesh, d_c_ex,         d_c_to_acc                )
    CALL move_a_and_c_to_ac_2D( mesh, d_a_ex, d_c_ex, d_a_and_c_to_ac           )
    CALL map_a_to_ac_2D(        mesh, d_a_ex,         d_a_to_ac                 )
    CALL map_c_to_ac_2D(        mesh, d_c_ex,         d_c_to_ac                 )
    CALL ddx_a_to_ac_2D(        mesh, d_a_ex,         ddx_a_to_ac               )
    CALL ddy_a_to_ac_2D(        mesh, d_a_ex,         ddy_a_to_ac               )
    CALL move_aca_to_a_2D(      mesh, d_ac_ex,        d_aca_to_a                )
    CALL move_acc_to_c_2D(      mesh, d_ac_ex,        d_acc_to_c                )
    CALL move_ac_to_a_and_c_2D( mesh, d_ac_ex,        d_ac_to_a      , d_ac_to_c)
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_a_ex         )
    CALL deallocate_shared( wd_c_ex         )
    CALL deallocate_shared( wd_ac_ex        )
    CALL deallocate_shared( wddx_ac_ex      )
    CALL deallocate_shared( wddy_ac_ex      )
    
    CALL deallocate_shared( wd_a_to_aca     )
    CALL deallocate_shared( wd_c_to_acc     )
    CALL deallocate_shared( wd_a_and_c_to_ac)
    CALL deallocate_shared( wd_a_to_ac      )
    CALL deallocate_shared( wd_c_to_ac      )
    CALL deallocate_shared( wddx_a_to_ac    )
    CALL deallocate_shared( wddy_a_to_ac    )
    CALL deallocate_shared( wd_aca_to_a     )
    CALL deallocate_shared( wd_acc_to_c     )
    CALL deallocate_shared( wd_ac_to_a      )
    CALL deallocate_shared( wd_ac_to_c      )
    
  END SUBROUTINE test_matrix_operators_a_c_ac_2D
  SUBROUTINE test_matrix_operators_a_c_ac_3D( mesh)
    ! Test all the matrix operators between the a/c-grids and the ac-grid
    !
    ! move_a_to_aca_3D
    ! move_c_to_acc_3D
    ! move a_and_c_to_ac_3D
    ! map_a_to_ac_3D
    ! map_c_to_ac_3D
    ! ddx_a_to_ac_3D
    ! ddy_a_to_ac_3D
    ! d_aca_to_a_3D
    ! d_acc_to_c_3D
    ! d_ac_to_a_and_c_3D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, aci, avi
    REAL(dp)                                           :: x, y, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_a_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_c_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_ac_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_ac_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_ac_ex
    INTEGER :: wd_a_ex, wd_c_ex, wd_ac_ex, wddx_ac_ex, wddy_ac_ex

    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_a_to_aca
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_c_to_acc
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_a_and_c_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_a_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_c_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_a_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_a_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_aca_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_acc_to_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_ac_to_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_ac_to_c
    INTEGER :: wd_a_to_aca, wd_c_to_acc, wd_a_and_c_to_ac, wd_a_to_ac, wd_c_to_ac, wddx_a_to_ac, wddy_a_to_ac, wd_aca_to_a, wd_acc_to_c, wd_ac_to_a, wd_ac_to_c
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (a/c - ac, 3D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
  
    CALL allocate_shared_dp_2D( mesh%nV    , C%nz,   d_a_ex         , wd_a_ex         )
    CALL allocate_shared_dp_2D( mesh%nAc   , C%nz,   d_c_ex         , wd_c_ex         )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz,   d_ac_ex        , wd_ac_ex        )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz,   ddx_ac_ex      , wddx_ac_ex      )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz,   ddy_ac_ex      , wddy_ac_ex      )
    
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz,   d_a_to_aca     , wd_a_to_aca     )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz,   d_c_to_acc     , wd_c_to_acc     )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz,   d_a_and_c_to_ac, wd_a_and_c_to_ac)
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz,   d_a_to_ac      , wd_a_to_ac      )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz,   d_c_to_ac      , wd_c_to_ac      )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz,   ddx_a_to_ac    , wddx_a_to_ac    )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz,   ddy_a_to_ac    , wddy_a_to_ac    )
    CALL allocate_shared_dp_2D( mesh%nV    , C%nz,   d_aca_to_a     , wd_aca_to_a     )
    CALL allocate_shared_dp_2D( mesh%nAc   , C%nz,   d_acc_to_c     , wd_acc_to_c     )
    CALL allocate_shared_dp_2D( mesh%nV    , C%nz,   d_ac_to_a      , wd_ac_to_a      )
    CALL allocate_shared_dp_2D( mesh%nAc   , C%nz,   d_ac_to_c      , wd_ac_to_c      )
    
  ! Calculate exact solutions
  ! =========================
    
    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2)
      d_a_ex(      vi,:) = d
    END DO
    DO aci = mesh%ci1, mesh%ci2
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2)
      d_c_ex(   aci,:) = d
    END DO
    DO avi = mesh%avi1, mesh%avi2
      x = mesh%VAaAc( avi,1)
      y = mesh%VAaAc( avi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2)
      d_ac_ex(   avi,:) = d
      ddx_ac_ex( avi,:) = ddx
      ddy_ac_ex( avi,:) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL move_a_to_aca_3D(      mesh, d_a_ex,         d_a_to_aca                )
    CALL move_c_to_acc_3D(      mesh, d_c_ex,         d_c_to_acc                )
    CALL move_a_and_c_to_ac_3D( mesh, d_a_ex, d_c_ex, d_a_and_c_to_ac           )
    CALL map_a_to_ac_3D(        mesh, d_a_ex,         d_a_to_ac                 )
    CALL map_c_to_ac_3D(        mesh, d_c_ex,         d_c_to_ac                 )
    CALL ddx_a_to_ac_3D(        mesh, d_a_ex,         ddx_a_to_ac               )
    CALL ddy_a_to_ac_3D(        mesh, d_a_ex,         ddy_a_to_ac               )
    CALL move_aca_to_a_3D(      mesh, d_ac_ex,        d_aca_to_a                )
    CALL move_acc_to_c_3D(      mesh, d_ac_ex,        d_acc_to_c                )
    CALL move_ac_to_a_and_c_3D( mesh, d_ac_ex,        d_ac_to_a      , d_ac_to_c)
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_a_ex         )
    CALL deallocate_shared( wd_c_ex         )
    CALL deallocate_shared( wd_ac_ex        )
    CALL deallocate_shared( wddx_ac_ex      )
    CALL deallocate_shared( wddy_ac_ex      )
    
    CALL deallocate_shared( wd_a_to_aca     )
    CALL deallocate_shared( wd_c_to_acc     )
    CALL deallocate_shared( wd_a_and_c_to_ac)
    CALL deallocate_shared( wd_a_to_ac      )
    CALL deallocate_shared( wd_c_to_ac      )
    CALL deallocate_shared( wddx_a_to_ac    )
    CALL deallocate_shared( wddy_a_to_ac    )
    CALL deallocate_shared( wd_aca_to_a     )
    CALL deallocate_shared( wd_acc_to_c     )
    CALL deallocate_shared( wd_ac_to_a      )
    CALL deallocate_shared( wd_ac_to_c      )
    
  END SUBROUTINE test_matrix_operators_a_c_ac_3D
    
  ! On the ac-grid
  SUBROUTINE test_matrix_operators_ac_ac_2D( mesh)
    ! Test all the matrix operators on the ac-grid
    ! 
    ! ddx_ac_to_ac_2D
    ! ddy_ac_to_ac_2D
    ! d2dx2_ac_to_ac_2D
    ! d2dxdy_ac_to_ac_2D
    ! d2dy2_ac_to_ac_2D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: avi
    REAL(dp)                                           :: x, y, d, ddx, ddy, d2dx2, d2dxdy, d2dy2
    
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ac_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ac_ex
    INTEGER :: wd_ac_ex, wddx_ac_ex, wddy_ac_ex

    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ac_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ac_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dx2_ac_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dxdy_ac_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dy2_ac_to_ac
    INTEGER :: wddx_ac_to_ac, wddy_ac_to_ac, wd2dx2_ac_to_ac, wd2dxdy_ac_to_ac, wd2dy2_ac_to_ac
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (ac - ac, 2D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
  
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d_ac_ex        , wd_ac_ex        )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddx_ac_ex      , wddx_ac_ex      )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddy_ac_ex      , wddy_ac_ex      )
    
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddx_ac_to_ac   , wddx_ac_to_ac   )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddy_ac_to_ac   , wddy_ac_to_ac   )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d2dx2_ac_to_ac , wd2dx2_ac_to_ac )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d2dxdy_ac_to_ac, wd2dxdy_ac_to_ac)
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d2dy2_ac_to_ac , wd2dy2_ac_to_ac )
    
  ! Calculate exact solutions
  ! =========================
    
    DO avi = mesh%avi1, mesh%avi2
      x = mesh%VAaAc( avi,1)
      y = mesh%VAaAc( avi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_ac_ex(   avi) = d
      ddx_ac_ex( avi) = ddx
      ddy_ac_ex( avi) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL ddx_ac_to_ac_2D(    mesh, d_ac_ex, ddx_ac_to_ac   )
    CALL ddy_ac_to_ac_2D(    mesh, d_ac_ex, ddy_ac_to_ac   )
    CALL d2dx2_ac_to_ac_2D(  mesh, d_ac_ex, d2dx2_ac_to_ac )
    CALL d2dxdy_ac_to_ac_2D( mesh, d_ac_ex, d2dxdy_ac_to_ac)
    CALL d2dy2_ac_to_ac_2D(  mesh, d_ac_ex, d2dy2_ac_to_ac )
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_ac_ex        )
    CALL deallocate_shared( wddx_ac_ex      )
    CALL deallocate_shared( wddy_ac_ex      )
    
    CALL deallocate_shared( wddx_ac_to_ac   )
    CALL deallocate_shared( wddy_ac_to_ac   )
    CALL deallocate_shared( wd2dx2_ac_to_ac )
    CALL deallocate_shared( wd2dxdy_ac_to_ac)
    CALL deallocate_shared( wd2dy2_ac_to_ac )
    
  END SUBROUTINE test_matrix_operators_ac_ac_2D
  SUBROUTINE test_matrix_operators_ac_ac_3D( mesh)
    ! Test all the matrix operators on the ac-grid
    ! 
    ! ddx_ac_to_ac_3D
    ! ddy_ac_to_ac_3D
    ! d2dx2_ac_to_ac_3D
    ! d2dxdy_ac_to_ac_3D
    ! d2dy2_ac_to_ac_3D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: avi
    REAL(dp)                                           :: x, y, d, ddx, ddy, d2dx2, d2dxdy, d2dy2
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_ac_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_ac_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_ac_ex
    INTEGER :: wd_ac_ex, wddx_ac_ex, wddy_ac_ex

    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_ac_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_ac_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d2dx2_ac_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d2dxdy_ac_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d2dy2_ac_to_ac
    INTEGER :: wddx_ac_to_ac, wddy_ac_to_ac, wd2dx2_ac_to_ac, wd2dxdy_ac_to_ac, wd2dy2_ac_to_ac
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (ac - ac, 3D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
  
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz, d_ac_ex        , wd_ac_ex        )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz, ddx_ac_ex      , wddx_ac_ex      )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz, ddy_ac_ex      , wddy_ac_ex      )
    
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz, ddx_ac_to_ac   , wddx_ac_to_ac   )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz, ddy_ac_to_ac   , wddy_ac_to_ac   )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz, d2dx2_ac_to_ac , wd2dx2_ac_to_ac )
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz, d2dxdy_ac_to_ac, wd2dxdy_ac_to_ac)
    CALL allocate_shared_dp_2D( mesh%nVAaAc, C%nz, d2dy2_ac_to_ac , wd2dy2_ac_to_ac )
    
  ! Calculate exact solutions
  ! =========================
    
    DO avi = mesh%avi1, mesh%avi2
      x = mesh%VAaAc( avi,1)
      y = mesh%VAaAc( avi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_ac_ex(   avi,:) = d
      ddx_ac_ex( avi,:) = ddx
      ddy_ac_ex( avi,:) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL ddx_ac_to_ac_3D(    mesh, d_ac_ex, ddx_ac_to_ac   )
    CALL ddy_ac_to_ac_3D(    mesh, d_ac_ex, ddy_ac_to_ac   )
    CALL d2dx2_ac_to_ac_3D(  mesh, d_ac_ex, d2dx2_ac_to_ac )
    CALL d2dxdy_ac_to_ac_3D( mesh, d_ac_ex, d2dxdy_ac_to_ac)
    CALL d2dy2_ac_to_ac_3D(  mesh, d_ac_ex, d2dy2_ac_to_ac )
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_ac_ex        )
    CALL deallocate_shared( wddx_ac_ex      )
    CALL deallocate_shared( wddy_ac_ex      )
    
    CALL deallocate_shared( wddx_ac_to_ac   )
    CALL deallocate_shared( wddy_ac_to_ac   )
    CALL deallocate_shared( wd2dx2_ac_to_ac )
    CALL deallocate_shared( wd2dxdy_ac_to_ac)
    CALL deallocate_shared( wd2dy2_ac_to_ac )
    
  END SUBROUTINE test_matrix_operators_ac_ac_3D
    
  ! Between the ac-grid and the bb-grid
  SUBROUTINE test_matrix_operators_ac_bb_2D( mesh)
    ! Test all the matrix operators between the ac-grid and the bb-grid
    ! 
    ! map_ac_to_bb_2D
    ! ddx_ac_to_bb_2D
    ! ddy_ac_to_bb_2D
    ! map_bb_to_ac_2D
    ! ddx_bb_to_ac_2D
    ! ddy_bb_to_ac_2D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: avi, ati
    REAL(dp)                                           :: x, y, d, ddx, ddy, d2dx2, d2dxdy, d2dy2
    
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ac_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ac_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_bb_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_bb_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_bb_ex
    INTEGER :: wd_ac_ex, wddx_ac_ex, wddy_ac_ex, wd_bb_ex, wddx_bb_ex, wddy_bb_ex

    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac_to_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ac_to_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ac_to_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: d_bb_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_bb_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_bb_to_ac
    INTEGER :: wd_ac_to_bb, wddx_ac_to_bb, wddy_ac_to_bb, wd_bb_to_ac, wddx_bb_to_ac, wddy_bb_to_ac
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (ac to bb, 2D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
  
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d_ac_ex        , wd_ac_ex        )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddx_ac_ex      , wddx_ac_ex      )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddy_ac_ex      , wddy_ac_ex      )
    CALL allocate_shared_dp_1D( mesh%nTriAaAc, d_bb_ex        , wd_bb_ex        )
    CALL allocate_shared_dp_1D( mesh%nTriAaAc, ddx_bb_ex      , wddx_bb_ex      )
    CALL allocate_shared_dp_1D( mesh%nTriAaAc, ddy_bb_ex      , wddy_bb_ex      )
    
    CALL allocate_shared_dp_1D( mesh%nTriAaAc, d_ac_to_bb     , wd_ac_to_bb     )
    CALL allocate_shared_dp_1D( mesh%nTriAaAc, ddx_ac_to_bb   , wddx_ac_to_bb   )
    CALL allocate_shared_dp_1D( mesh%nTriAaAc, ddy_ac_to_bb   , wddy_ac_to_bb   )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   d_bb_to_ac     , wd_bb_to_ac     )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddx_bb_to_ac   , wddx_bb_to_ac   )
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   ddy_bb_to_ac   , wddy_bb_to_ac   )
    
  ! Calculate exact solutions
  ! =========================
    
    DO avi = mesh%avi1, mesh%avi2
      x = mesh%VAaAc( avi,1)
      y = mesh%VAaAc( avi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_ac_ex(   avi) = d
      ddx_ac_ex( avi) = ddx
      ddy_ac_ex( avi) = ddy
    END DO
    DO ati = mesh%ati1, mesh%ati2
      x = mesh%TriGCAaAc( ati,1)
      y = mesh%TriGCAaAc( ati,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_bb_ex(   ati) = d
      ddx_bb_ex( ati) = ddx
      ddy_bb_ex( ati) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL map_ac_to_bb_2D( mesh, d_ac_ex, d_ac_to_bb  )
    CALL ddx_ac_to_bb_2D( mesh, d_ac_ex, ddx_ac_to_bb)
    CALL ddy_ac_to_bb_2D( mesh, d_ac_ex, ddy_ac_to_bb)
    CALL map_bb_to_ac_2D( mesh, d_bb_ex, d_bb_to_ac  )
    CALL ddx_bb_to_ac_2D( mesh, d_bb_ex, ddx_bb_to_ac)
    CALL ddy_bb_to_ac_2D( mesh, d_bb_ex, ddy_bb_to_ac)
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_ac_ex     )
    CALL deallocate_shared( wddx_ac_ex   )
    CALL deallocate_shared( wddy_ac_ex   )
    CALL deallocate_shared( wd_bb_ex     )
    CALL deallocate_shared( wddx_bb_ex   ) 
    CALL deallocate_shared( wddy_bb_ex   )
      
    CALL deallocate_shared( wd_ac_to_bb  )
    CALL deallocate_shared( wddx_ac_to_bb)
    CALL deallocate_shared( wddy_ac_to_bb)
    CALL deallocate_shared( wd_bb_to_ac  )
    CALL deallocate_shared( wddx_bb_to_ac)
    CALL deallocate_shared( wddy_bb_to_ac)
    
  END SUBROUTINE test_matrix_operators_ac_bb_2D
  SUBROUTINE test_matrix_operators_ac_bb_3D( mesh)
    ! Test all the matrix operators between the ac-grid and the bb-grid
    ! 
    ! map_ac_to_bb_3D
    ! ddx_ac_to_bb_3D
    ! ddy_ac_to_bb_3D
    ! map_bb_to_ac_3D
    ! ddx_bb_to_ac_3D
    ! ddy_bb_to_ac_3D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: avi, ati
    REAL(dp)                                           :: x, y, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_ac_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_ac_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_ac_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_bb_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_bb_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_bb_ex
    INTEGER :: wd_ac_ex, wddx_ac_ex, wddy_ac_ex, wd_bb_ex, wddx_bb_ex, wddy_bb_ex

    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_ac_to_bb
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_ac_to_bb
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_ac_to_bb
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_bb_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_bb_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_bb_to_ac
    INTEGER :: wd_ac_to_bb, wddx_ac_to_bb, wddy_ac_to_bb, wd_bb_to_ac, wddx_bb_to_ac, wddy_bb_to_ac
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (ac to bb, 3D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
  
    CALL allocate_shared_dp_2D( mesh%nVAaAc,   C%nz, d_ac_ex        , wd_ac_ex        )
    CALL allocate_shared_dp_2D( mesh%nVAaAc,   C%nz, ddx_ac_ex      , wddx_ac_ex      )
    CALL allocate_shared_dp_2D( mesh%nVAaAc,   C%nz, ddy_ac_ex      , wddy_ac_ex      )
    CALL allocate_shared_dp_2D( mesh%nTriAaAc, C%nz, d_bb_ex        , wd_bb_ex        )
    CALL allocate_shared_dp_2D( mesh%nTriAaAc, C%nz, ddx_bb_ex      , wddx_bb_ex      )
    CALL allocate_shared_dp_2D( mesh%nTriAaAc, C%nz, ddy_bb_ex      , wddy_bb_ex      )
    
    CALL allocate_shared_dp_2D( mesh%nTriAaAc, C%nz, d_ac_to_bb     , wd_ac_to_bb     )
    CALL allocate_shared_dp_2D( mesh%nTriAaAc, C%nz, ddx_ac_to_bb   , wddx_ac_to_bb   )
    CALL allocate_shared_dp_2D( mesh%nTriAaAc, C%nz, ddy_ac_to_bb   , wddy_ac_to_bb   )
    CALL allocate_shared_dp_2D( mesh%nVAaAc,   C%nz, d_bb_to_ac     , wd_bb_to_ac     )
    CALL allocate_shared_dp_2D( mesh%nVAaAc,   C%nz, ddx_bb_to_ac   , wddx_bb_to_ac   )
    CALL allocate_shared_dp_2D( mesh%nVAaAc,   C%nz, ddy_bb_to_ac   , wddy_bb_to_ac   )
    
  ! Calculate exact solutions
  ! =========================
    
    DO avi = mesh%avi1, mesh%avi2
      x = mesh%VAaAc( avi,1)
      y = mesh%VAaAc( avi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2)
      d_ac_ex(   avi,:) = d
      ddx_ac_ex( avi,:) = ddx
      ddy_ac_ex( avi,:) = ddy
    END DO
    DO ati = mesh%ati1, mesh%ati2
      x = mesh%TriGCAaAc( ati,1)
      y = mesh%TriGCAaAc( ati,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2)
      d_bb_ex(   ati,:) = d
      ddx_bb_ex( ati,:) = ddx
      ddy_bb_ex( ati,:) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL map_ac_to_bb_3D( mesh, d_ac_ex, d_ac_to_bb  )
    CALL ddx_ac_to_bb_3D( mesh, d_ac_ex, ddx_ac_to_bb)
    CALL ddy_ac_to_bb_3D( mesh, d_ac_ex, ddy_ac_to_bb)
    CALL map_bb_to_ac_3D( mesh, d_bb_ex, d_bb_to_ac  )
    CALL ddx_bb_to_ac_3D( mesh, d_bb_ex, ddx_bb_to_ac)
    CALL ddy_bb_to_ac_3D( mesh, d_bb_ex, ddy_bb_to_ac)
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_ac_ex     )
    CALL deallocate_shared( wddx_ac_ex   )
    CALL deallocate_shared( wddy_ac_ex   )
    CALL deallocate_shared( wd_bb_ex     )
    CALL deallocate_shared( wddx_bb_ex   ) 
    CALL deallocate_shared( wddy_bb_ex   )
      
    CALL deallocate_shared( wd_ac_to_bb  )
    CALL deallocate_shared( wddx_ac_to_bb)
    CALL deallocate_shared( wddy_ac_to_bb)
    CALL deallocate_shared( wd_bb_to_ac  )
    CALL deallocate_shared( wddx_bb_to_ac)
    CALL deallocate_shared( wddy_bb_to_ac)
    
  END SUBROUTINE test_matrix_operators_ac_bb_3D
    
  ! Between the ac-grid and the acu/acv-grids
  SUBROUTINE test_matrix_operators_ac_acuv_2D( mesh)
    ! Test all the matrix operators between the ac-grid and the acu/acv-grids
    !
    ! move_ac_to_acu_2D
    ! move_ac_to_acv_2D
    ! move_acu_to_ac_2D
    ! move_acv_to_ac_2D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: avi, auvi
    REAL(dp)                                           :: x, y, d, ddx, ddy, d2dx2, d2dxdy, d2dy2
    
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acu_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acv_ex
    INTEGER :: wd_ac_ex, wd_acu_ex, wd_acv_ex

    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac_to_acu
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac_to_acv
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acu_to_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acv_to_ac
    INTEGER :: wd_ac_to_acu, wd_ac_to_acv, wd_acu_to_ac, wd_acv_to_ac
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (ac - acu/acv, 2D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
  
    CALL allocate_shared_dp_1D(   mesh%nVAaAc, d_ac_ex         , wd_ac_ex         )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, d_acu_ex        , wd_acu_ex        )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, d_acv_ex        , wd_acv_ex        )
    
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, d_ac_to_acu     , wd_ac_to_acu     )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc, d_ac_to_acv     , wd_ac_to_acv     )
    CALL allocate_shared_dp_1D(   mesh%nVAaAc, d_acu_to_ac     , wd_acu_to_ac     )
    CALL allocate_shared_dp_1D(   mesh%nVAaAc, d_acv_to_ac     , wd_acv_to_ac     )
    
  ! Calculate exact solutions
  ! =========================
    
    DO avi = mesh%avi1, mesh%avi2
      x = mesh%VAaAc( avi,1)
      y = mesh%VAaAc( avi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_ac_ex(      avi) = d
    END DO
    DO auvi = mesh%auvi1, mesh%auvi2
      IF (MOD( auvi,2) == 1) THEN
        ! u
        avi = (auvi+1)/2
        x = mesh%VAaAc( avi,1)
        y = mesh%VAaAc( avi,2)
        CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
        d_acu_ex(   auvi) = d
      END IF
    END DO
    DO auvi = mesh%auvi1, mesh%auvi2
      IF (MOD( auvi,2) == 0) THEN
        ! v
        avi = auvi/2
        x = mesh%VAaAc( avi,1)
        y = mesh%VAaAc( avi,2)
        CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
        d_acv_ex(   auvi) = d
      END IF
    END DO
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL move_ac_to_acu_2D( mesh, d_ac_ex , d_ac_to_acu)
    CALL move_ac_to_acv_2D( mesh, d_ac_ex , d_ac_to_acv)
    CALL move_acu_to_ac_2D( mesh, d_acu_ex, d_acu_to_ac)
    CALL move_acv_to_ac_2D( mesh, d_acv_ex, d_acv_to_ac)
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_ac_ex         )
    CALL deallocate_shared( wd_acu_ex        )
    CALL deallocate_shared( wd_acv_ex        )
    
    CALL deallocate_shared( wd_ac_to_acu     )
    CALL deallocate_shared( wd_ac_to_acv     )
    CALL deallocate_shared( wd_acu_to_ac     )
    CALL deallocate_shared( wd_acv_to_ac     )
    
  END SUBROUTINE test_matrix_operators_ac_acuv_2D
  SUBROUTINE test_matrix_operators_ac_acuv_3D( mesh)
    ! Test all the matrix operators between the ac-grid and the acu/acv-grids
    !
    ! move_ac_to_acu_3D
    ! move_ac_to_acv_3D
    ! move_acu_to_ac_3D
    ! move_acv_to_ac_3D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: avi, auvi
    REAL(dp)                                           :: x, y, d, ddx, ddy, d2dx2, d2dxdy, d2dy2
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_ac_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_acu_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_acv_ex
    INTEGER :: wd_ac_ex, wd_acu_ex, wd_acv_ex

    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_ac_to_acu
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_ac_to_acv
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_acu_to_ac
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_acv_to_ac
    INTEGER :: wd_ac_to_acu, wd_ac_to_acv, wd_acu_to_ac, wd_acv_to_ac
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (ac - acu/acv, 3D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
  
    CALL allocate_shared_dp_2D(   mesh%nVAaAc, C%nz, d_ac_ex         , wd_ac_ex         )
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc, C%nz, d_acu_ex        , wd_acu_ex        )
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc, C%nz, d_acv_ex        , wd_acv_ex        )
    
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc, C%nz, d_ac_to_acu     , wd_ac_to_acu     )
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc, C%nz, d_ac_to_acv     , wd_ac_to_acv     )
    CALL allocate_shared_dp_2D(   mesh%nVAaAc, C%nz, d_acu_to_ac     , wd_acu_to_ac     )
    CALL allocate_shared_dp_2D(   mesh%nVAaAc, C%nz, d_acv_to_ac     , wd_acv_to_ac     )
    
  ! Calculate exact solutions
  ! =========================
    
    DO avi = mesh%avi1, mesh%avi2
      x = mesh%VAaAc( avi,1)
      y = mesh%VAaAc( avi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_ac_ex(      avi,:) = d
    END DO
    DO auvi = mesh%auvi1, mesh%auvi2
      IF (MOD( auvi,2) == 1) THEN
        ! u
        avi = (auvi+1)/2
        x = mesh%VAaAc( avi,1)
        y = mesh%VAaAc( avi,2)
        CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
        d_acu_ex(   auvi,:) = d
      END IF
    END DO
    DO auvi = mesh%auvi1, mesh%auvi2
      IF (MOD( auvi,2) == 0) THEN
        ! v
        avi = auvi/2
        x = mesh%VAaAc( avi,1)
        y = mesh%VAaAc( avi,2)
        CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
        d_acv_ex(   auvi,:) = d
      END IF
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL move_ac_to_acu_3D( mesh, d_ac_ex , d_ac_to_acu)
    CALL move_ac_to_acv_3D( mesh, d_ac_ex , d_ac_to_acv)
    CALL move_acu_to_ac_3D( mesh, d_acu_ex, d_acu_to_ac)
    CALL move_acv_to_ac_3D( mesh, d_acv_ex, d_acv_to_ac)
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_ac_ex         )
    CALL deallocate_shared( wd_acu_ex        )
    CALL deallocate_shared( wd_acv_ex        )
    
    CALL deallocate_shared( wd_ac_to_acu     )
    CALL deallocate_shared( wd_ac_to_acv     )
    CALL deallocate_shared( wd_acu_to_ac     )
    CALL deallocate_shared( wd_acv_to_ac     )
    
  END SUBROUTINE test_matrix_operators_ac_acuv_3D
  
  ! Between the acu/acv-grids and the bb-grid
  SUBROUTINE test_matrix_operators_acuv_bb_2D( mesh)
    ! Test all the matrix operators between the acu/acv-grids and the bb-grid
    !
    ! map_acu_to_bb_2D
    ! ddx_acu_to_bb_2D
    ! ddy_acu_to_bb_2D
    ! map_acv_to_bb_2D
    ! ddx_acv_to_bb_2D
    ! ddy_acv_to_bb_2D
    ! 
    ! map_bb_to_acu_2D
    ! ddx_bb_to_acu_2D
    ! ddy_bb_to_acu_2D
    ! map_bb_to_acv_2D
    ! ddx_bb_to_acv_2D
    ! ddy_bb_to_acv_2D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: auvi, avi, ati
    REAL(dp)                                           :: x, y, d, ddx, ddy, d2dx2, d2dxdy, d2dy2
    
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acu_ex    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_acu_ex    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_acu_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acv_ex    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_acv_ex    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_acv_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_bb_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_bb_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_bb_ex
    INTEGER :: wd_acu_ex, wddx_acu_ex, wddy_acu_ex, wd_acv_ex, wddx_acv_ex, wddy_acv_ex, wd_bb_ex, wddx_bb_ex, wddy_bb_ex

    REAL(dp), DIMENSION(:    ), POINTER                :: d_acu_to_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_acu_to_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_acu_to_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acv_to_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_acv_to_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_acv_to_bb
    INTEGER :: wd_acu_to_bb, wddx_acu_to_bb, wddy_acu_to_bb, wd_acv_to_bb, wddx_acv_to_bb, wddy_acv_to_bb
    
    REAL(dp), DIMENSION(:    ), POINTER                :: d_bb_to_acu
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_bb_to_acu
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_bb_to_acu
    REAL(dp), DIMENSION(:    ), POINTER                :: d_bb_to_acv
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_bb_to_acv
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_bb_to_acv
    INTEGER :: wd_bb_to_acu, wddx_bb_to_acu, wddy_bb_to_acu, wd_bb_to_acv, wddx_bb_to_acv, wddy_bb_to_acv
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (acu/acv - bb, 2D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
  
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   d_acu_ex     , wd_acu_ex     )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   ddx_acu_ex   , wddx_acu_ex   )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   ddy_acu_ex   , wddy_acu_ex   )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   d_acv_ex     , wd_acv_ex     )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   ddx_acv_ex   , wddx_acv_ex   )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   ddy_acv_ex   , wddy_acv_ex   )
    CALL allocate_shared_dp_1D(   mesh%nTriAaAc, d_bb_ex      , wd_bb_ex      )
    CALL allocate_shared_dp_1D(   mesh%nTriAaAc, ddx_bb_ex    , wddx_bb_ex    )
    CALL allocate_shared_dp_1D(   mesh%nTriAaAc, ddy_bb_ex    , wddy_bb_ex    )
    
    CALL allocate_shared_dp_1D(   mesh%nTriAaAc, d_acu_to_bb  , wd_acu_to_bb  )
    CALL allocate_shared_dp_1D(   mesh%nTriAaAc, ddx_acu_to_bb, wddx_acu_to_bb)
    CALL allocate_shared_dp_1D(   mesh%nTriAaAc, ddy_acu_to_bb, wddy_acu_to_bb)
    CALL allocate_shared_dp_1D(   mesh%nTriAaAc, d_acv_to_bb  , wd_acv_to_bb  )
    CALL allocate_shared_dp_1D(   mesh%nTriAaAc, ddx_acv_to_bb, wddx_acv_to_bb)
    CALL allocate_shared_dp_1D(   mesh%nTriAaAc, ddy_acv_to_bb, wddy_acv_to_bb)
    
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   d_bb_to_acu  , wd_bb_to_acu  )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   ddx_bb_to_acu, wddx_bb_to_acu)
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   ddy_bb_to_acu, wddy_bb_to_acu)
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   d_bb_to_acv  , wd_bb_to_acv  )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   ddx_bb_to_acv, wddx_bb_to_acv)
    CALL allocate_shared_dp_1D( 2*mesh%nVAaAc,   ddy_bb_to_acv, wddy_bb_to_acv)
    
  ! Calculate exact solutions
  ! =========================
    
    DO auvi = mesh%auvi1, mesh%auvi2
      IF (MOD( auvi,2) == 1) THEN
        ! u
        avi = (auvi+1)/2
        x = mesh%VAaAc( avi,1)
        y = mesh%VAaAc( avi,2)
        CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
        d_acu_ex(   auvi) = d
        ddx_acu_ex( auvi) = ddx
        ddy_acu_ex( auvi) = ddy
      END IF
    END DO
    DO auvi = mesh%auvi1, mesh%auvi2
      IF (MOD( auvi,2) == 0) THEN
        ! v
        avi = auvi/2
        x = mesh%VAaAc( avi,1)
        y = mesh%VAaAc( avi,2)
        CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
        d_acv_ex(   auvi) = d
        ddx_acv_ex( auvi) = ddx
        ddy_acv_ex( auvi) = ddy
      END IF
    END DO
    DO ati = mesh%ati1, mesh%ati2
      x = mesh%TriGCAaAc( ati,1)
      y = mesh%TriGCAaAc( ati,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_bb_ex(   ati) = d
      ddx_bb_ex( ati) = ddx
      ddy_bb_ex( ati) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL map_acu_to_bb_2D( mesh, d_acu_ex, d_acu_to_bb  )
    CALL ddx_acu_to_bb_2D( mesh, d_acu_ex, ddx_acu_to_bb)
    CALL ddy_acu_to_bb_2D( mesh, d_acu_ex, ddy_acu_to_bb)
    CALL map_acv_to_bb_2D( mesh, d_acv_ex, d_acv_to_bb  )
    CALL ddx_acv_to_bb_2D( mesh, d_acv_ex, ddx_acv_to_bb)
    CALL ddy_acv_to_bb_2D( mesh, d_acv_ex, ddy_acv_to_bb)
    
    CALL map_bb_to_acu_2D( mesh, d_bb_ex , d_bb_to_acu  )
    CALL ddx_bb_to_acu_2D( mesh, d_bb_ex , ddx_bb_to_acu)
    CALL ddy_bb_to_acu_2D( mesh, d_bb_ex , ddy_bb_to_acu)
    CALL map_bb_to_acv_2D( mesh, d_bb_ex , d_bb_to_acv  )
    CALL ddx_bb_to_acv_2D( mesh, d_bb_ex , ddx_bb_to_acv)
    CALL ddy_bb_to_acv_2D( mesh, d_bb_ex , ddy_bb_to_acv)
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_acu_ex     )
    CALL deallocate_shared( wddx_acu_ex   )
    CALL deallocate_shared( wddy_acu_ex   )
    CALL deallocate_shared( wd_acv_ex     )
    CALL deallocate_shared( wddx_acv_ex   )
    CALL deallocate_shared( wddy_acv_ex   )
    CALL deallocate_shared( wd_bb_ex      )
    CALL deallocate_shared( wddx_bb_ex    )
    CALL deallocate_shared( wddy_bb_ex    )
    
    CALL deallocate_shared( wd_acu_to_bb  )
    CALL deallocate_shared( wddx_acu_to_bb)
    CALL deallocate_shared( wddy_acu_to_bb)
    CALL deallocate_shared( wd_acv_to_bb  )
    CALL deallocate_shared( wddx_acv_to_bb)
    CALL deallocate_shared( wddy_acv_to_bb)
    
    CALL deallocate_shared( wd_bb_to_acu  )
    CALL deallocate_shared( wddx_bb_to_acu)
    CALL deallocate_shared( wddy_bb_to_acu)
    CALL deallocate_shared( wd_bb_to_acv  )
    CALL deallocate_shared( wddx_bb_to_acv)
    CALL deallocate_shared( wddy_bb_to_acv)
    
  END SUBROUTINE test_matrix_operators_acuv_bb_2D
  SUBROUTINE test_matrix_operators_acuv_bb_3D( mesh)
    ! Test all the matrix operators between the acu/acv-grids and the bb-grid
    !
    ! map_acu_to_bb_3D
    ! ddx_acu_to_bb_3D
    ! ddy_acu_to_bb_3D
    ! map_acv_to_bb_3D
    ! ddx_acv_to_bb_3D
    ! ddy_acv_to_bb_3D
    ! 
    ! map_bb_to_acu_3D
    ! ddx_bb_to_acu_3D
    ! ddy_bb_to_acu_3D
    ! map_bb_to_acv_3D
    ! ddx_bb_to_acv_3D
    ! ddy_bb_to_acv_3D
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: auvi, avi, ati
    REAL(dp)                                           :: x, y, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_acu_ex    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_acu_ex    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_acu_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_acv_ex    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_acv_ex    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_acv_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_bb_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_bb_ex
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_bb_ex
    INTEGER :: wd_acu_ex, wddx_acu_ex, wddy_acu_ex, wd_acv_ex, wddx_acv_ex, wddy_acv_ex, wd_bb_ex, wddx_bb_ex, wddy_bb_ex

    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_acu_to_bb
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_acu_to_bb
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_acu_to_bb
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_acv_to_bb
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_acv_to_bb
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_acv_to_bb
    INTEGER :: wd_acu_to_bb, wddx_acu_to_bb, wddy_acu_to_bb, wd_acv_to_bb, wddx_acv_to_bb, wddy_acv_to_bb
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_bb_to_acu
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_bb_to_acu
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_bb_to_acu
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_bb_to_acv
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddx_bb_to_acv
    REAL(dp), DIMENSION(:,:  ), POINTER                :: ddy_bb_to_acv
    INTEGER :: wd_bb_to_acu, wddx_bb_to_acu, wddy_bb_to_acu, wd_bb_to_acv, wddx_bb_to_acv, wddy_bb_to_acv
        
    IF (par%master) WRITE(0,*) '  Testing the matrix operators (acu/acv - bb, 3D)...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
  
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, d_acu_ex     , wd_acu_ex     )
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, ddx_acu_ex   , wddx_acu_ex   )
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, ddy_acu_ex   , wddy_acu_ex   )
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, d_acv_ex     , wd_acv_ex     )
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, ddx_acv_ex   , wddx_acv_ex   )
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, ddy_acv_ex   , wddy_acv_ex   )
    CALL allocate_shared_dp_2D(   mesh%nTriAaAc, C%nz, d_bb_ex      , wd_bb_ex      )
    CALL allocate_shared_dp_2D(   mesh%nTriAaAc, C%nz, ddx_bb_ex    , wddx_bb_ex    )
    CALL allocate_shared_dp_2D(   mesh%nTriAaAc, C%nz, ddy_bb_ex    , wddy_bb_ex    )
    
    CALL allocate_shared_dp_2D(   mesh%nTriAaAc, C%nz, d_acu_to_bb  , wd_acu_to_bb  )
    CALL allocate_shared_dp_2D(   mesh%nTriAaAc, C%nz, ddx_acu_to_bb, wddx_acu_to_bb)
    CALL allocate_shared_dp_2D(   mesh%nTriAaAc, C%nz, ddy_acu_to_bb, wddy_acu_to_bb)
    CALL allocate_shared_dp_2D(   mesh%nTriAaAc, C%nz, d_acv_to_bb  , wd_acv_to_bb  )
    CALL allocate_shared_dp_2D(   mesh%nTriAaAc, C%nz, ddx_acv_to_bb, wddx_acv_to_bb)
    CALL allocate_shared_dp_2D(   mesh%nTriAaAc, C%nz, ddy_acv_to_bb, wddy_acv_to_bb)
    
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, d_bb_to_acu  , wd_bb_to_acu  )
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, ddx_bb_to_acu, wddx_bb_to_acu)
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, ddy_bb_to_acu, wddy_bb_to_acu)
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, d_bb_to_acv  , wd_bb_to_acv  )
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, ddx_bb_to_acv, wddx_bb_to_acv)
    CALL allocate_shared_dp_2D( 2*mesh%nVAaAc  , C%nz, ddy_bb_to_acv, wddy_bb_to_acv)
    
  ! Calculate exact solutions
  ! =========================
    
    DO auvi = mesh%auvi1, mesh%auvi2
      IF (MOD( auvi,2) == 1) THEN
        ! u
        avi = (auvi+1)/2
        x = mesh%VAaAc( avi,1)
        y = mesh%VAaAc( avi,2)
        CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2)
        d_acu_ex(   auvi,:) = d
        ddx_acu_ex( auvi,:) = ddx
        ddy_acu_ex( auvi,:) = ddy
      END IF
    END DO
    DO auvi = mesh%auvi1, mesh%auvi2
      IF (MOD( auvi,2) == 0) THEN
        ! v
        avi = auvi/2
        x = mesh%VAaAc( avi,1)
        y = mesh%VAaAc( avi,2)
        CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2)
        d_acv_ex(   auvi,:) = d
        ddx_acv_ex( auvi,:) = ddx
        ddy_acv_ex( auvi,:) = ddy
      END IF
    END DO
    DO ati = mesh%ati1, mesh%ati2
      x = mesh%TriGCAaAc( ati,1)
      y = mesh%TriGCAaAc( ati,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d3Dx2, d3Dxdy, d3Dy2)
      d_bb_ex(   ati,:) = d
      ddx_bb_ex( ati,:) = ddx
      ddy_bb_ex( ati,:) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL map_acu_to_bb_3D( mesh, d_acu_ex, d_acu_to_bb  )
    CALL ddx_acu_to_bb_3D( mesh, d_acu_ex, ddx_acu_to_bb)
    CALL ddy_acu_to_bb_3D( mesh, d_acu_ex, ddy_acu_to_bb)
    CALL map_acv_to_bb_3D( mesh, d_acv_ex, d_acv_to_bb  )
    CALL ddx_acv_to_bb_3D( mesh, d_acv_ex, ddx_acv_to_bb)
    CALL ddy_acv_to_bb_3D( mesh, d_acv_ex, ddy_acv_to_bb)
    
    CALL map_bb_to_acu_3D( mesh, d_bb_ex , d_bb_to_acu  )
    CALL ddx_bb_to_acu_3D( mesh, d_bb_ex , ddx_bb_to_acu)
    CALL ddy_bb_to_acu_3D( mesh, d_bb_ex , ddy_bb_to_acu)
    CALL map_bb_to_acv_3D( mesh, d_bb_ex , d_bb_to_acv  )
    CALL ddx_bb_to_acv_3D( mesh, d_bb_ex , ddx_bb_to_acv)
    CALL ddy_bb_to_acv_3D( mesh, d_bb_ex , ddy_bb_to_acv)
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_acu_ex     )
    CALL deallocate_shared( wddx_acu_ex   )
    CALL deallocate_shared( wddy_acu_ex   )
    CALL deallocate_shared( wd_acv_ex     )
    CALL deallocate_shared( wddx_acv_ex   )
    CALL deallocate_shared( wddy_acv_ex   )
    CALL deallocate_shared( wd_bb_ex      )
    CALL deallocate_shared( wddx_bb_ex    )
    CALL deallocate_shared( wddy_bb_ex    )
    
    CALL deallocate_shared( wd_acu_to_bb  )
    CALL deallocate_shared( wddx_acu_to_bb)
    CALL deallocate_shared( wddy_acu_to_bb)
    CALL deallocate_shared( wd_acv_to_bb  )
    CALL deallocate_shared( wddx_acv_to_bb)
    CALL deallocate_shared( wddy_acv_to_bb)
    
    CALL deallocate_shared( wd_bb_to_acu  )
    CALL deallocate_shared( wddx_bb_to_acu)
    CALL deallocate_shared( wddy_bb_to_acu)
    CALL deallocate_shared( wd_bb_to_acv  )
    CALL deallocate_shared( wddx_bb_to_acv)
    CALL deallocate_shared( wddy_bb_to_acv)
    
  END SUBROUTINE test_matrix_operators_acuv_bb_3D
  
  SUBROUTINE test_matrix_operators_test_function( x, y, xmin, xmax, ymin, ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
    ! Simple function for testing all the matrix operators
      
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x
    REAL(dp),                            INTENT(IN)    :: y
    REAL(dp),                            INTENT(IN)    :: xmin
    REAL(dp),                            INTENT(IN)    :: xmax
    REAL(dp),                            INTENT(IN)    :: ymin
    REAL(dp),                            INTENT(IN)    :: ymax
    REAL(dp),                            INTENT(OUT)   :: d
    REAL(dp),                            INTENT(OUT)   :: ddx
    REAL(dp),                            INTENT(OUT)   :: ddy
    REAL(dp),                            INTENT(OUT)   :: d2dx2
    REAL(dp),                            INTENT(OUT)   :: d2dxdy
    REAL(dp),                            INTENT(OUT)   :: d2dy2
    
    ! Local variables:
    REAL(dp)                                           :: xp, yp
    
    xp = 2._dp * pi * (x - xmin) / (xmax - xmin)
    yp = 2._dp * pi * (y - ymin) / (ymax - ymin)
    
    d      =  SIN( xp) * SIN( yp)
    ddx    =  COS( xp) * SIN( yp) *  2._dp * pi / (xmax - xmin)
    ddy    =  SIN( xp) * COS( yp) *  2._dp * pi / (xmax - xmin)
    d2dx2  = -SIN( xp) * SIN( yp) * (2._dp * pi / (xmax - xmin))**2
    d2dxdy =  COS( xp) * COS( yp) * (2._dp * pi / (xmax - xmin))**2
    d2dy2  = -SIN( xp) * SIN( yp) * (2._dp * pi / (xmax - xmin))**2
        
  END SUBROUTINE test_matrix_operators_test_function
  
! == Solve a (modified version of) the Laplace equation on the mesh
  SUBROUTINE solve_modified_Laplace_equation( mesh)
    ! Test the discretisation by solving (a modified version of) the Laplace equation:
    !
    ! d/dx ( N * df/dx) + d/dx ( N * df/dy) = 0
    !
    ! (Modified to include a spatially variable stiffness, making it much more similar to
    ! the SSA/DIVA, and necessitating the use of a staggered grid/mesh to solve it)
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, ti
    REAL(dp), DIMENSION(:    ), POINTER                ::  f_a,  N_b,  b_a
    INTEGER                                            :: wf_a, wN_b, wb_a
    TYPE(type_sparse_matrix_CSR)                       :: M_a, BC_a, mNddx_b
    REAL(dp)                                           :: x, y, xp, yp
    REAL(dp), PARAMETER                                :: amp = 0._dp
    INTEGER                                            :: nnz_BC
        
    IF (par%master) WRITE(0,*) '  Solving the modified Laplace equation...'
    
  ! Solve the equation on the a-grid (vertex)
  ! =========================================
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV,   f_a, wf_a)
    CALL allocate_shared_dp_1D( mesh%nV,   b_a, wb_a)
    CALL allocate_shared_dp_1D( mesh%nTri, N_b, wN_b)
    
    ! Set up the spatially variable stiffness N on the b-grid (triangle)
    DO ti = mesh%ti1, mesh%ti2
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      xp = (x - mesh%xmin) / (mesh%xmax - mesh%xmin)
      yp = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
      N_b( ti) = 1._dp + EXP(amp * SIN( 2._dp * pi * xp) * SIN(2._dp * pi * yp))
    END DO
    CALL sync
    
    ! Create the operator matrix M
    CALL multiply_matrix_rows_with_vector( mesh%M_ddx_a_b, N_b, mNddx_b     )
    CALL multiply_matrix_matrix_CSR(       mesh%M_ddx_b_a,      mNddx_b, M_a)
    
    ! Create the boundary conditions matrix BC
    
    ! Find maximum number of non-zero entries
    nnz_BC = 0
    DO vi = mesh%vi1, mesh%vi2
      IF (mesh%edge_index( vi) > 0) THEN
        nnz_BC = nnz_BC + 1 + mesh%nC( vi)
      END IF
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nnz_BC, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)   
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( BC_a, mesh%nV, mesh%nV, nnz_BC)
    
    ! Fill in values
    BC_a%nnz = 0
    BC_a%ptr = 1
    DO vi = mesh%vi1, mesh%vi2
    
      IF (mesh%edge_index( vi) == 0) THEN
        ! Free vertex: no boundary conditions apply
      
      ELSEIF (mesh%edge_index( vi) == 6 .OR. mesh%edge_index( vi) == 7 .OR. mesh%edge_index( vi) == 8) THEN
        ! West: bump
        
        ! Matrix
        BC_a%nnz  = BC_a%nnz+1
        BC_a%index( BC_a%nnz) = vi
        BC_a%val(   BC_a%nnz) = 1._dp
        
        ! Right-hand side vector
        y = mesh%V( vi,2)
        y = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
        b_a( vi) = 0.5_dp * (1._dp - COS( 2._dp * pi * y))
        
      ELSE
        ! Otherwise: zero
        
        ! Matrix
        BC_a%nnz  = BC_a%nnz+1
        BC_a%index( BC_a%nnz) = vi
        BC_a%val(   BC_a%nnz) = 1._dp
        
        ! Right-hand side vector
        b_a( vi) = 0._dp
        
      END IF
      
      ! Finalise matrix row
      BC_a%ptr( vi+1) = BC_a%nnz + 1
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Finalise matrix
    CALL finalise_matrix_CSR_dist( BC_a, mesh%vi1, mesh%vi2)
    
    ! Add boundary conditions to the operator matrix
    CALL overwrite_rows_CSR( M_a, BC_a)
    
    ! Solve the equation
    CALL solve_matrix_equation_CSR( M_a, b_a, f_a, C%DIVA_choice_matrix_solver, &
      SOR_nit = 5000, SOR_tol = 0.0001_dp, SOR_omega = 1.3_dp, &
      PETSc_rtol = C%DIVA_PETSc_rtol, PETSc_abstol = C%DIVA_PETSc_abstol)
      
    ! Write result to debug file
    IF (par%master) THEN
      debug%dp_2D_a_01 = f_a
      CALL write_to_debug_file
    END IF
    CALL sync
    
  END SUBROUTINE solve_modified_Laplace_equation
  SUBROUTINE solve_modified_Laplace_equation_combi( mesh)
    ! Test the discretisation by solving (a modified version of) the Laplace equation:
    !
    ! d/dx ( N * df/dx) + d/dx ( N * df/dy) = 0
    !
    ! (Modified to include a spatially variable stiffness, making it much more similar to
    ! the SSA/DIVA, and necessitating the use of a staggered grid/mesh to solve it)
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: avi, ati, vi, aci, edge_index
    REAL(dp), DIMENSION(:    ), POINTER                ::  f_ac,  N_bb,  b_ac
    INTEGER                                            :: wf_ac, wN_bb, wb_ac
    TYPE(type_sparse_matrix_CSR)                       :: M_ac, BC_ac, mNddx_bb
    REAL(dp)                                           :: x, y, xp, yp
    REAL(dp), PARAMETER                                :: amp = 0._dp
    INTEGER                                            :: nnz_BC
        
    IF (par%master) WRITE(0,*) '  Solving the modified Laplace equation on the combined ac-grid...'
    
  ! Solve the equation on the a-grid (vertex)
  ! =========================================
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   f_ac, wf_ac)
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   b_ac, wb_ac)
    CALL allocate_shared_dp_1D( mesh%nTriAaAc, N_bb, wN_bb)
    
    ! Set up the spatially variable stiffness N on the bb-grid
    DO ati = mesh%ati1, mesh%ati2
      x = mesh%TriGCAaAc( ati,1)
      y = mesh%TriGCAaAc( ati,2)
      xp = (x - mesh%xmin) / (mesh%xmax - mesh%xmin)
      yp = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
      N_bb( ati) = 1._dp + EXP(amp * SIN( 2._dp * pi * xp) * SIN(2._dp * pi * yp))
    END DO
    CALL sync
    
    ! Create the operator matrix M
    CALL multiply_matrix_rows_with_vector( mesh%M_ddx_ac_bb, N_bb, mNddx_bb      )
    CALL multiply_matrix_matrix_CSR(       mesh%M_ddx_bb_ac,       mNddx_bb, M_ac)
    
    ! Create the boundary conditions matrix BC
    
    ! Find maximum number of non-zero entries
    nnz_BC = 0
    DO avi = mesh%avi1, mesh%avi2
    
      IF (avi <= mesh%nV) THEN
        vi = avi
        edge_index = mesh%edge_index( vi)
      ELSE
        aci = avi - mesh%nV
        edge_index = mesh%edge_index_ac( aci)
      END IF
      
      IF (edge_index > 0) THEN
        nnz_BC = nnz_BC + 1 + mesh%nCAaAc( avi)
      END IF
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nnz_BC, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)   
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( BC_ac, mesh%nVAaAc, mesh%nVAaAc, nnz_BC)
    
    ! Fill in values
    BC_ac%nnz = 0
    BC_ac%ptr = 1
    DO avi = mesh%avi1, mesh%avi2
    
      IF (avi <= mesh%nV) THEN
        vi = avi
        edge_index = mesh%edge_index( vi)
      ELSE
        aci = avi - mesh%nV
        edge_index = mesh%edge_index_ac( aci)
      END IF
    
      IF (edge_index == 0) THEN
        ! Free vertex: no boundary conditions apply
        
      ELSEIF (edge_index == 6 .OR. edge_index == 7 .OR. edge_index == 8) THEN
        ! West: bump
        
        ! Matrix
        BC_ac%nnz  = BC_ac%nnz+1
        BC_ac%index( BC_ac%nnz) = avi
        BC_ac%val(   BC_ac%nnz) = 1._dp
        
        ! Right-hand side vector
        y = mesh%VAaAc( avi,2)
        y = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
        b_ac( avi) = 0.5_dp * (1._dp - COS( 2._dp * pi * y))
        
      ELSE
        ! Otherwise: zero
        
        ! Matrix
        BC_ac%nnz  = BC_ac%nnz+1
        BC_ac%index( BC_ac%nnz) = avi
        BC_ac%val(   BC_ac%nnz) = 1._dp
        
        ! Right-hand side vector
        b_ac( avi) = 0._dp
        
      END IF
      
      ! Finalise matrix row
      BC_ac%ptr( avi+1 : mesh%nVAaAc) = BC_ac%nnz + 1
      
    END DO ! DO avi = mesh%avi1, mesh%avi2
    CALL sync
    
    ! Finalise matrix
    CALL finalise_matrix_CSR_dist( BC_ac, mesh%avi1, mesh%avi2)
    
    ! Add boundary conditions to the operator matrix
    CALL overwrite_rows_CSR( M_ac, BC_ac)
    
    ! Solve the equation
    CALL solve_matrix_equation_CSR( M_ac, b_ac, f_ac, C%DIVA_choice_matrix_solver, &
      SOR_nit = 5000, SOR_tol = 0.0001_dp, SOR_omega = 1.7_dp, &
      PETSc_rtol = C%DIVA_PETSc_rtol, PETSc_abstol = C%DIVA_PETSc_abstol)
      
    ! Write result to debug file
    IF (par%master) THEN
      debug%dp_2D_ac_01 = f_ac
      CALL write_to_debug_file
    END IF
    CALL sync
    
  END SUBROUTINE solve_modified_Laplace_equation_combi

END MODULE tests_and_checks_module
