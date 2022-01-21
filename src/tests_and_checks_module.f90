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
  USE sparse_matrix_module
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
    CALL test_matrix_operators_abc_2D( mesh)
    CALL test_matrix_operators_abc_3D( mesh)
    
    ! 2nd-order accurate matrix operators on the b-grid
    CALL test_matrix_operators_2nd_order_b_to_b_2D( mesh)
    
  ! Solve the (modified) Laplace equation as the ultimate test
  ! ==========================================================
  
    CALL solve_modified_Laplace_equation_b( mesh)
    
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
    CALL test_CSR_matrix_operations_multiply_vector
    CALL test_CSR_matrix_operations_multiply_matrix
    
  END SUBROUTINE test_CSR_matrix_operations
  
  SUBROUTINE test_CSR_matrix_operations_add
    ! Testing some basic CSR-formatted matrix operations
    
    USE mesh_operators_module
    USE utilities_module
    
    IMPLICIT NONE
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR_dp)                    :: AAA, BBB, CCC, CCC_ex
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
    
    ! Perform multiplication
    CALL add_matrix_matrix_CSR( AAA, BBB, CCC)
    
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
  SUBROUTINE test_CSR_matrix_operations_multiply_vector
    ! Testing some basic CSR-formatted matrix operations
    
    USE mesh_operators_module
    USE utilities_module
    
    IMPLICIT NONE
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR_dp)                    :: AAA
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
    TYPE(type_sparse_matrix_CSR_dp)                    :: AAA, BBB, CCC, CCC_ex
    LOGICAL                                            :: are_identical
    
    IF (par%master) WRITE(0,*) '   Testing matrix-matrix multiplication (C = A*B)...'
    
    ! Define some simple matrices
    CALL allocate_matrix_CSR_shared( AAA   , 4, 4, 8 )
    CALL allocate_matrix_CSR_shared( BBB   , 4, 4, 8 )
    CALL allocate_matrix_CSR_shared( CCC_ex, 4, 4, 13)
    
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
      ! C
      CCC_ex%nnz = 13
      CCC_ex%ptr   = [1,4,8,10,14]
      CCC_ex%index = [1,2,3,1,2,3,4,3,4,1,2,3,4]
      CCC_ex%val   = [9,26,38,27,60,30,64,55,60,54,211,158,128]
    END IF
    CALL sync
    
    CALL write_CSR_matrix_to_NetCDF( AAA, 'A.nc')
    CALL write_CSR_matrix_to_NetCDF( BBB, 'B.nc')
    CALL write_CSR_matrix_to_NetCDF( CCC_ex, 'C_ex.nc')
    
    ! Perform multiplication
    CALL multiply_matrix_matrix_CSR( AAA, BBB, CCC)
    CALL write_CSR_matrix_to_NetCDF( CCC, 'C.nc')
    
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
    
  END SUBROUTINE test_CSR_matrix_operations_multiply_matrix
  
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
!    CALl MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    
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
  
  ! 2nd-order accurate matrix operators on the b-grid
  SUBROUTINE test_matrix_operators_2nd_order_b_to_b_2D( mesh)
    ! Test all the 2nd-order accurate matrix operators on the b-grid
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables:
    INTEGER                                            :: ti
    REAL(dp)                                           :: x, y, d, ddx, ddy, d2dx2, d2dxdy, d2dy2
    
    REAL(dp), DIMENSION(:    ), POINTER                :: d_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dx2_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dxdy_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dy2_b_ex
    INTEGER :: wd_b_ex, wddx_b_ex, wddy_b_ex, wd2dx2_b_ex, wd2dxdy_b_ex, wd2dy2_b_ex
    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_b
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dx2_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dxdy_b
    REAL(dp), DIMENSION(:    ), POINTER                :: d2dy2_b
    INTEGER :: wddx_b, wddy_b, wd2dx2_b, wd2dxdy_b, wd2dy2_b
        
    IF (par%master) WRITE(0,*) '  Testing the 2nd-order accurate matrix operators on the b-grid...'
    CALL sync
    
  ! Allocate shared memory
  ! ======================
    
    CALL allocate_shared_dp_1D( mesh%nTri, d_b_ex,          wd_b_ex         )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_b_ex,        wddx_b_ex       )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_b_ex,        wddy_b_ex       )
    CALL allocate_shared_dp_1D( mesh%nTri, d2dx2_b_ex,      wd2dx2_b_ex     )
    CALL allocate_shared_dp_1D( mesh%nTri, d2dxdy_b_ex,     wd2dxdy_b_ex    )
    CALL allocate_shared_dp_1D( mesh%nTri, d2dy2_b_ex,      wd2dy2_b_ex     )
    
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_b,           wddx_b          )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_b,           wddy_b          )
    CALL allocate_shared_dp_1D( mesh%nTri, d2dx2_b,         wd2dx2_b        )
    CALL allocate_shared_dp_1D( mesh%nTri, d2dxdy_b,        wd2dxdy_b       )
    CALL allocate_shared_dp_1D( mesh%nTri, d2dy2_b,         wd2dy2_b        )
    
  ! Calculate exact solutions
  ! =========================
    
    DO ti = mesh%ti1, mesh%ti2
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy, d2dx2, d2dxdy, d2dy2)
      d_b_ex(      ti) = d
      ddx_b_ex(    ti) = ddx
      ddy_b_ex(    ti) = ddy
      d2dx2_b_ex(  ti) = d2dx2
      d2dxdy_b_ex( ti) = d2dxdy
      d2dy2_b_ex(  ti) = d2dy2
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
    
    CALL multiply_matrix_vector_CSR( mesh%M2_ddx_b_b   , d_b_ex, ddx_b   )
    CALL multiply_matrix_vector_CSR( mesh%M2_ddy_b_b   , d_b_ex, ddy_b   )
    CALL multiply_matrix_vector_CSR( mesh%M2_d2dx2_b_b , d_b_ex, d2dx2_b )
    CALL multiply_matrix_vector_CSR( mesh%M2_d2dxdy_b_b, d_b_ex, d2dxdy_b)
    CALL multiply_matrix_vector_CSR( mesh%M2_d2dy2_b_b , d_b_ex, d2dy2_b )
    
  ! Write results to debug file
  ! ===========================
  
!    IF (par%master) THEN
!      
!      ! Exact solutions
!      debug%dp_2D_b_01 = ddx_b_ex
!      debug%dp_2D_b_02 = ddy_b_ex
!      debug%dp_2D_b_03 = d2dx2_b_ex
!      debug%dp_2D_b_04 = d2dxdy_b_ex
!      debug%dp_2D_b_05 = d2dy2_b_ex
!      
!      ! Discretised approximations
!      debug%dp_2D_b_06 = ddx_b
!      debug%dp_2D_b_07 = ddy_b
!      debug%dp_2D_b_08 = d2dx2_b
!      debug%dp_2D_b_09 = d2dxdy_b
!      debug%dp_2D_b_10 = d2dy2_b
!      
!      CALL write_to_debug_file
!      
!    END IF
!    CALL sync
!    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_b_ex     )
    CALL deallocate_shared( wddx_b_ex   )
    CALL deallocate_shared( wddy_b_ex   )
    CALL deallocate_shared( wd2dx2_b_ex )
    CALL deallocate_shared( wd2dxdy_b_ex)
    CALL deallocate_shared( wd2dy2_b_ex )
    
    CALL deallocate_shared( wddx_b      )
    CALL deallocate_shared( wddy_b      )
    CALL deallocate_shared( wd2dx2_b    )
    CALL deallocate_shared( wd2dxdy_b   )
    CALL deallocate_shared( wd2dy2_b    )
    
  END SUBROUTINE test_matrix_operators_2nd_order_b_to_b_2D
  
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
  SUBROUTINE solve_modified_Laplace_equation_b( mesh)
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
    REAL(dp), DIMENSION(:    ), POINTER                ::  N_a,  N_b,  dN_dx_b,  dN_dy_b
    INTEGER                                            :: wN_a, wN_b, wdN_dx_b, wdN_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                ::  x_b,  b_b
    INTEGER                                            :: wx_b, wb_b
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, nnz_max, nnz_max_proc, n, edge_index
    INTEGER                                            :: k1, k2, k
    TYPE(type_sparse_matrix_CSR_dp)                    :: A_b
    REAL(dp)                                           :: x, y, xp, yp
    REAL(dp), PARAMETER                                :: amp = 7._dp
        
    IF (par%master) WRITE(0,*) '  Solving the modified Laplace equation on the b-grid...'
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV  , N_a    , wN_a    )
    CALL allocate_shared_dp_1D( mesh%nTri, N_b    , wN_b    )
    CALL allocate_shared_dp_1D( mesh%nTri, dN_dx_b, wdN_dx_b)
    CALL allocate_shared_dp_1D( mesh%nTri, dN_dy_b, wdN_dy_b)
    CALL allocate_shared_dp_1D( mesh%nTri, x_b    , wx_b    )
    CALL allocate_shared_dp_1D( mesh%nTri, b_b    , wb_b    )
    
    ! Calculate N on the a-grid
    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      xp = (x - mesh%xmin) / (mesh%xmax - mesh%xmin)
      yp = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
      N_a( vi) = 1._dp + EXP(amp * SIN( 2._dp * pi * xp) * SIN(2._dp * pi * yp))
    END DO
    CALL sync
        
    ! Calculate N, dN/dx, and dN/dy on the b-grid
    CALL map_a_to_b_2D( mesh, N_a, N_b    )
    CALL ddx_a_to_b_2D( mesh, N_a, dN_dx_b)
    CALL ddy_a_to_b_2D( mesh, N_a, dN_dy_b)
    
    ! Set up the matrix A
    ncols           = mesh%nTri    ! from
    nrows           = mesh%nTri    ! to
    nnz_per_row_max = 10
    
    nnz_max         = nrows * nnz_per_row_max
    nnz_max_proc    = CEILING( REAL(nnz_max,dp) / REAL(par%n,dp))
    CALL allocate_matrix_CSR_dist( A_b, nrows, ncols, nnz_max_proc)
    
    DO ti = mesh%ti1, mesh%ti2
      
      ! If this triangle touches the domain border, apply boundary conditions
      edge_index = 0
      DO n = 1, 3
        vi = mesh%Tri( ti,n)
        IF (mesh%edge_index( vi) > 0) edge_index = MAX( edge_index, mesh%edge_index( vi))
      END DO
      
      IF (edge_index > 0) THEN
        ! This triangle touches the domain border; apply boundary conditions
        
        IF (edge_index == 6 .OR. edge_index == 7 .OR. edge_index == 8) THEN
          ! Western border: bump
        
          ! Matrix
          A_b%nnz  = A_b%nnz+1
          A_b%index( A_b%nnz) = ti
          A_b%val(   A_b%nnz) = 1._dp
          ! Extend memory if necessary
          IF (A_b%nnz > A_b%nnz_max - 1000) CALL extend_matrix_CSR_dist( A_b, A_b%nnz_max + 1000)
        
          ! Right-hand side and initial guess
          y = mesh%V( vi,2)
          y = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
          b_b( ti) = 0.5_dp * (1._dp - COS( 2._dp * pi * y))
          x_b( ti) = b_b( ti)
        
        ELSE
          ! Other borders: zero
        
          ! Matrix
          A_b%nnz  = A_b%nnz+1
          A_b%index( A_b%nnz) = ti
          A_b%val(   A_b%nnz) = 1._dp
          ! Extend memory if necessary
          IF (A_b%nnz > A_b%nnz_max - 1000) CALL extend_matrix_CSR_dist( A_b, A_b%nnz_max + 1000)
        
          ! Right-hand side and initial guess
          b_b( ti) = 0._dp
          x_b( ti) = 0._dp
          
        END IF
        
      ELSE ! IF (edge_index > 0) THEN
        ! This is a free triangle; fill in matrix coefficients
        
        ! Matrix
        k1 = mesh%M2_ddx_b_b%ptr( ti)
        k2 = mesh%M2_ddx_b_b%ptr( ti+1) - 1
        
        DO k = k1, k2
        
          A_b%nnz  = A_b%nnz+1
          A_b%index( A_b%nnz) = mesh%M2_d2dx2_b_b%index( k)
          A_b%val(   A_b%nnz) = N_b(     ti) * mesh%M2_d2dx2_b_b%val( k) + &
                                dN_dx_b( ti) * mesh%M2_ddx_b_b%val(   k) + &
                                N_b(     ti) * mesh%M2_d2dy2_b_b%val( k) + &
                                dN_dy_b( ti) * mesh%M2_ddy_b_b%val(   k)
          ! Extend memory if necessary
          IF (A_b%nnz > A_b%nnz_max - 1000) CALL extend_matrix_CSR_dist( A_b, A_b%nnz_max + 1000)
          
        END DO
        
        ! Right-hand side and initial guess
        b_b( ti) = 0._dp
        x_b( ti) = 0._dp
        
      END IF ! IF (edge_index > 0) THEN
      
      ! Finalise matrix row
      A_b%ptr( ti+1: A_b%m+1) = A_b%nnz+1
      
    END DO
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( A_b, mesh%ti1, mesh%ti2)
    
    ! Solve the equation
    CALL solve_matrix_equation_CSR( A_b, b_b, x_b, &
      choice_matrix_solver = C%DIVA_choice_matrix_solver, &
      SOR_nit = 5000, SOR_tol = 0.00001_dp, SOR_omega = 1.3_dp, &
      PETSc_rtol = 0.001_dp, PETSc_abstol = 0.001_dp)
      
    ! Write result to debug file
    IF (par%master) THEN
      debug%dp_2D_b_01 = x_b
      CALL write_to_debug_file
    END IF
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_matrix_CSR( A_b)
    CALL deallocate_shared( wN_a)
    CALL deallocate_shared( wN_b)
    CALL deallocate_shared( wdN_dx_b)
    CALL deallocate_shared( wdN_dy_b)
    CALL deallocate_shared( wx_b)
    CALL deallocate_shared( wb_b)
    
  END SUBROUTINE solve_modified_Laplace_equation_b

END MODULE tests_and_checks_module
