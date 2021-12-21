MODULE tests_and_checks_module

  ! The main regional ice-sheet model

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
  USE netcdf_module,                   ONLY: write_CSR_matrix_to_NetCDF

  IMPLICIT NONE

CONTAINS
  
! == Test basic CSR matrix operations
  SUBROUTINE test_CSR_matrix_operations
    ! Testing some basic CSR-formatted matrix operations
    
    IMPLICIT NONE
        
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Testing basic CSR matrix operations...'
    
    CALL test_CSR_matrix_operations_add    
    CALL test_CSR_matrix_operations_overwrite_rows
    CALL test_CSR_matrix_operations_multiply_vector
    CALL test_CSR_matrix_operations_multiply_matrix
        
    IF (par%master) WRITE(0,*) ' Finished testing basic CSR matrix operations!'
    
  END SUBROUTINE test_CSR_matrix_operations
  
  SUBROUTINE test_CSR_matrix_operations_add
    ! Testing some basic CSR-formatted matrix operations
    
    USE mesh_operators_module
    USE utilities_module
    
    IMPLICIT NONE
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR)                       :: AAA, BBB, CCC
    
    IF (par%master) WRITE(0,*) '  Testing matrix addition (C = A+B)...'
    
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
    CALL add_matrix_matrix_CSR( AAA, BBB, CCC)
    
    ! Check if operation gave the correct result
    IF (par%master) THEN
      IF (CCC%nnz == 11 .AND. &
          CCC%ptr(    1) ==  1 .AND. &
          CCC%ptr(    2) ==  3 .AND. &
          CCC%ptr(    3) ==  6 .AND. &
          CCC%ptr(    4) ==  8 .AND. &
          CCC%ptr(    5) == 12 .AND. &
          CCC%index(  1) == 1 .AND. &
          CCC%index(  2) == 3 .AND. &
          CCC%index(  3) == 1 .AND. &
          CCC%index(  4) == 3 .AND. &
          CCC%index(  5) == 4 .AND. &
          CCC%index(  6) == 2 .AND. &
          CCC%index(  7) == 3 .AND. &
          CCC%index(  8) == 1 .AND. &
          CCC%index(  9) == 2 .AND. &
          CCC%index( 10) == 3 .AND. &
          CCC%index( 11) == 4 .AND. &
          CCC%val(    1) ==  10._dp .AND. &
          CCC%val(    2) ==  12._dp .AND. &
          CCC%val(    3) ==   3._dp .AND. &
          CCC%val(    4) ==  11._dp .AND. &
          CCC%val(    5) ==  16._dp .AND. &
          CCC%val(    6) ==  18._dp .AND. &
          CCC%val(    7) ==  14._dp .AND. &
          CCC%val(    8) ==   6._dp .AND. &
          CCC%val(    9) ==  15._dp .AND. &
          CCC%val(   10) ==   7._dp .AND. &
          CCC%val(   11) ==  24._dp) THEN
      ELSE
        WRITE(0,*) 'test_CSR_matrix_operations_add - ERROR: CSR matrix addition gave wrong answer!'
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
    
  END SUBROUTINE test_CSR_matrix_operations_add
  SUBROUTINE test_CSR_matrix_operations_overwrite_rows
    ! Testing some basic CSR-formatted matrix operations
    
    USE mesh_operators_module
    USE utilities_module
    
    IMPLICIT NONE
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR)                       :: AAA, BBB
    
    IF (par%master) WRITE(0,*) '  Testing matrix row overwriting (C = A (+O) B)...'
    
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
    
    IF (par%master) WRITE(0,*) '  Testing matrix-vector multiplication (b = A*x)...'
    
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
    
    IF (par%master) WRITE(0,*) '  Testing matrix-matrix multiplication (C = A*B)...'
    
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
  
! == Test all the matrix operators (mapping+gradients)
  SUBROUTINE test_matrix_operators( mesh)
    ! Test all the matrix operators
    
    USE mesh_operators_module
    USE utilities_module
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, ti, aci, ai, auvi, auvi1, auvi2
    REAL(dp)                                           :: x, y, d, ddx, ddy
    
    ! Exact solutions
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
    
    ! Discretised approximations
    REAL(dp), DIMENSION(:    ), POINTER                :: d_aa
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ab
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ba
    REAL(dp), DIMENSION(:    ), POINTER                :: d_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: d_bc
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ca
    REAL(dp), DIMENSION(:    ), POINTER                :: d_cb
    REAL(dp), DIMENSION(:    ), POINTER                :: d_cc
    INTEGER :: wd_aa, wd_ab, wd_ac, wd_ba, wd_bb, wd_bc, wd_ca, wd_cb, wd_cc
    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_aa
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ab
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ba
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_bc
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ca
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_cb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_cc
    INTEGER :: wddx_aa, wddx_ab, wddx_ac, wddx_ba, wddx_bb, wddx_bc, wddx_ca, wddx_cb, wddx_cc
    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_aa
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ab
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ba
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_bc
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ca
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_cb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_cc
    INTEGER :: wddy_aa, wddy_ab, wddy_ac, wddy_ba, wddy_bb, wddy_bc, wddy_ca, wddy_cb, wddy_cc
    
    ! Data on the combined mesh
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acau
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acav
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acau_acb
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acav_acb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_acau_acb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_acav_acb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_acau_acb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_acav_acb
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acb_acau
    REAL(dp), DIMENSION(:    ), POINTER                :: d_acb_acav
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_acb_acau
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_acb_acav
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_acb_acau
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_acb_acav
    INTEGER :: wd_acau, wd_acav
    INTEGER :: wd_acau_acb, wd_acav_acb, wddx_acau_acb, wddx_acav_acb, wddy_acau_acb, wddy_acav_acb
    INTEGER :: wd_acb_acau, wd_acb_acav, wddx_acb_acau, wddx_acb_acav, wddy_acb_acau, wddy_acb_acav
        
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Testing the matrix operators...'
    CALL sync
    
  ! Write all the matrix operators to files, for double-checking stuff in Matlab
  ! ============================================================================
    
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_a_b, 'M_map_a_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_a_c, 'M_map_a_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_b_a, 'M_map_b_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_b_c, 'M_map_b_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_c_a, 'M_map_c_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_c_b, 'M_map_c_b.nc')
    
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_a_a, 'M_ddx_a_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_a_b, 'M_ddx_a_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_a_c, 'M_ddx_a_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_b_a, 'M_ddx_b_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_b_b, 'M_ddx_b_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_b_c, 'M_ddx_b_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_c_a, 'M_ddx_c_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_c_b, 'M_ddx_c_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_c_c, 'M_ddx_c_c.nc')
    
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_a_a, 'M_ddy_a_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_a_b, 'M_ddy_a_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_a_c, 'M_ddy_a_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_b_a, 'M_ddy_b_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_b_b, 'M_ddy_b_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_b_c, 'M_ddy_b_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_c_a, 'M_ddy_c_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_c_b, 'M_ddy_c_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_c_c, 'M_ddy_c_c.nc')
    
  ! Allocate shared memory
  ! ======================
    
    ! Exact solutions
    CALL allocate_shared_dp_1D( mesh%nV,   d_a_ex,   wd_a_ex  )
    CALL allocate_shared_dp_1D( mesh%nTri, d_b_ex,   wd_b_ex  )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_c_ex,   wd_c_ex  )
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_a_ex, wddx_a_ex)
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_b_ex, wddx_b_ex)
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_c_ex, wddx_c_ex)
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_a_ex, wddy_a_ex)
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_b_ex, wddy_b_ex)
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_c_ex, wddy_c_ex)
    
    ! Discretised approximations
    CALL allocate_shared_dp_1D( mesh%nV,   d_aa,     wd_aa    )
    CALL allocate_shared_dp_1D( mesh%nTri, d_ab,     wd_ab    )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_ac,     wd_ac    )
    CALL allocate_shared_dp_1D( mesh%nV,   d_ba,     wd_ba    )
    CALL allocate_shared_dp_1D( mesh%nTri, d_bb,     wd_bb    )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_bc,     wd_bc    )
    CALL allocate_shared_dp_1D( mesh%nV,   d_ca,     wd_ca    )
    CALL allocate_shared_dp_1D( mesh%nTri, d_cb,     wd_cb    )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_cc,     wd_cc    )
    
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_aa,   wddx_aa  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_ab,   wddx_ab  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_ac,   wddx_ac  )
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_ba,   wddx_ba  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_bb,   wddx_bb  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_bc,   wddx_bc  )
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_ca,   wddx_ca  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_cb,   wddx_cb  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_cc,   wddx_cc  )
    
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_aa,   wddy_aa  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_ab,   wddy_ab  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_ac,   wddy_ac  )
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_ba,   wddy_ba  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_bb,   wddy_bb  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_bc,   wddy_bc  )
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_ca,   wddy_ca  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_cb,   wddy_cb  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_cc,   wddy_cc  )
    
    ! Data on the combined mesh
    CALL allocate_shared_dp_1D( 2*mesh%nVAaac,   d_acau        , wd_acau        )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaac,   d_acav        , wd_acav        )
    CALL allocate_shared_dp_1D(   mesh%nTriAaac, d_acau_acb    , wd_acau_acb    )
    CALL allocate_shared_dp_1D(   mesh%nTriAaac, d_acav_acb    , wd_acav_acb    )
    CALL allocate_shared_dp_1D(   mesh%nTriAaac, ddx_acau_acb  , wddx_acau_acb  )
    CALL allocate_shared_dp_1D(   mesh%nTriAaac, ddx_acav_acb  , wddx_acav_acb  )
    CALL allocate_shared_dp_1D(   mesh%nTriAaac, ddy_acau_acb  , wddy_acau_acb  )
    CALL allocate_shared_dp_1D(   mesh%nTriAaac, ddy_acav_acb  , wddy_acav_acb  )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaac,   d_acb_acau    , wd_acb_acau    )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaac,   d_acb_acav    , wd_acb_acav    )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaac,   ddx_acb_acau  , wddx_acb_acau  )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaac,   ddx_acb_acav  , wddx_acb_acav  )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaac,   ddy_acb_acau  , wddy_acb_acau  )
    CALL allocate_shared_dp_1D( 2*mesh%nVAaac,   ddy_acb_acav  , wddy_acb_acav  )
    
  ! Calculate exact solutions
  ! =========================
    
    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy)
      d_a_ex(   vi) = d
      ddx_a_ex( vi) = ddx
      ddy_a_ex( vi) = ddy
    END DO
    DO ti = mesh%ti1, mesh%ti2
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy)
      d_b_ex(   ti) = d
      ddx_b_ex( ti) = ddx
      ddy_b_ex( ti) = ddy
    END DO
    DO aci = mesh%ci1, mesh%ci2
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      CALL test_matrix_operators_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy)
      d_c_ex(   aci) = d
      ddx_c_ex( aci) = ddx
      ddy_c_ex( aci) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
  
    ! Mapping
    d_aa( mesh%vi1:mesh%vi2) = d_a_ex( mesh%vi1:mesh%vi2)
    CALL map_a_to_b_2D( mesh, d_a_ex, d_ab)
    CALL map_a_to_c_2D( mesh, d_a_ex, d_ac)
    
    CALL map_b_to_a_2D( mesh, d_b_ex, d_ba)
    d_bb( mesh%ti1:mesh%ti2) = d_b_ex( mesh%ti1:mesh%ti2)
    CALL map_b_to_c_2D( mesh, d_b_ex, d_bc)
    
    CALL map_c_to_a_2D( mesh, d_c_ex, d_ca)
    CALL map_c_to_b_2D( mesh, d_c_ex, d_cb)
    d_cc( mesh%ci1:mesh%ci2) = d_c_ex( mesh%ci1:mesh%ci2)
    
    ! d/dx
    CALL ddx_a_to_a_2D( mesh, d_a_ex, ddx_aa)
    CALL ddx_a_to_b_2D( mesh, d_a_ex, ddx_ab)
    CALL ddx_a_to_c_2D( mesh, d_a_ex, ddx_ac)
    
    CALL ddx_b_to_a_2D( mesh, d_b_ex, ddx_ba)
    CALL ddx_b_to_b_2D( mesh, d_b_ex, ddx_bb)
    CALL ddx_b_to_c_2D( mesh, d_b_ex, ddx_bc)
    
    CALL ddx_c_to_a_2D( mesh, d_c_ex, ddx_ca)
    CALL ddx_c_to_b_2D( mesh, d_c_ex, ddx_cb)
    CALL ddx_c_to_c_2D( mesh, d_c_ex, ddx_cc)
    
    ! d/dy
    CALL ddy_a_to_a_2D( mesh, d_a_ex, ddy_aa)
    CALL ddy_a_to_b_2D( mesh, d_a_ex, ddy_ab)
    CALL ddy_a_to_c_2D( mesh, d_a_ex, ddy_ac)
    
    CALL ddy_b_to_a_2D( mesh, d_b_ex, ddy_ba)
    CALL ddy_b_to_b_2D( mesh, d_b_ex, ddy_bb)
    CALL ddy_b_to_c_2D( mesh, d_b_ex, ddy_bc)
    
    CALL ddy_c_to_a_2D( mesh, d_c_ex, ddy_ca)
    CALL ddy_c_to_b_2D( mesh, d_c_ex, ddy_cb)
    CALL ddy_c_to_c_2D( mesh, d_c_ex, ddy_cc)
    
    ! Write results to debug file
    IF (par%master) THEN
    
      ! a-grid (vertex)
      debug%dp_2D_a_01 = d_aa   - d_a_ex
      debug%dp_2D_a_02 = d_ba   - d_a_ex
      debug%dp_2D_a_03 = d_ca   - d_a_ex
      debug%dp_2D_a_04 = ddx_aa - ddx_a_ex
      debug%dp_2D_a_05 = ddx_ba - ddx_a_ex
      debug%dp_2D_a_06 = ddx_ca - ddx_a_ex
      debug%dp_2D_a_07 = ddy_aa - ddy_a_ex
      debug%dp_2D_a_08 = ddy_ba - ddy_a_ex
      debug%dp_2D_a_09 = ddy_ca - ddy_a_ex
    
      ! b-grid (triangle)
      debug%dp_2D_b_01 = d_ab   - d_b_ex
      debug%dp_2D_b_02 = d_bb   - d_b_ex
      debug%dp_2D_b_03 = d_cb   - d_b_ex
      debug%dp_2D_b_04 = ddx_ab - ddx_b_ex
      debug%dp_2D_b_05 = ddx_bb - ddx_b_ex
      debug%dp_2D_b_06 = ddx_cb - ddx_b_ex
      debug%dp_2D_b_07 = ddy_ab - ddy_b_ex
      debug%dp_2D_b_08 = ddy_bb - ddy_b_ex
      debug%dp_2D_b_09 = ddy_cb - ddy_b_ex
    
      ! c-grid (edge)
      debug%dp_2D_c_01 = d_ac   - d_c_ex
      debug%dp_2D_c_02 = d_bc   - d_c_ex
      debug%dp_2D_c_03 = d_cc   - d_c_ex
      debug%dp_2D_c_04 = ddx_ac - ddx_c_ex
      debug%dp_2D_c_05 = ddx_bc - ddx_c_ex
      debug%dp_2D_c_06 = ddx_cc - ddx_c_ex
      debug%dp_2D_c_07 = ddy_ac - ddy_c_ex
      debug%dp_2D_c_08 = ddy_bc - ddy_c_ex
      debug%dp_2D_c_09 = ddy_cc - ddy_c_ex
      
      CALL write_to_debug_file
      
    END IF ! IF (par%master) THEN
    CALL sync
    
  ! Stuff on the combined ACUV-mesh
  ! ===============================
    
    CALL partition_list( 2*mesh%nVAaAc, par%i, par%n, auvi1, auvi2)
  
    ! Map data from regular a/c-grids to combined ac-grid (same data on u- and v-fields)
    CALL map_a_and_c_to_acau( mesh, d_a_ex, d_c_ex, d_acau)
    CALL map_a_and_c_to_acav( mesh, d_a_ex, d_c_ex, d_acav)
    
    ! Calculate mapped data and derivatives on the combi-triangles
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_acau_acb, d_acau, ddx_acau_acb)
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_acav_acb, d_acav, ddx_acav_acb)
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_acau_acb, d_acau, ddy_acau_acb)
    CALL multiply_matrix_vector_CSR( mesh%M_ddy_acav_acb, d_acav, ddy_acav_acb)
    
    ! Map everything back to the combi-vertices for writing to output
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_acb_acau, ddx_acau_acb, ddx_acb_acau)
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_acb_acav, ddx_acav_acb, ddx_acb_acav)
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_acb_acau, ddy_acau_acb, ddy_acb_acau)
    CALL multiply_matrix_vector_CSR( mesh%M_ddx_acb_acav, ddy_acav_acb, ddy_acb_acav)
    
    ! Map data to single-field AC-mesh for writing to debug output
    DO auvi = auvi1, auvi2
      IF (MOD( auvi,2) == 1) THEN
        ! u
        ai = (auvi+1) / 2
        debug%dp_2D_ac_01( ai) = d_acau(       auvi)
        debug%dp_2D_ac_03( ai) = ddx_acb_acau( auvi)
        debug%dp_2D_ac_05( ai) = ddy_acb_acau( auvi)
      ELSE
        ! v
        ai =  auvi    / 2
        debug%dp_2D_ac_02( ai) = -d_acav(      auvi)
        debug%dp_2D_ac_04( ai) = ddx_acb_acav( auvi)
        debug%dp_2D_ac_06( ai) = ddy_acb_acav( auvi)
      END IF
    END DO ! DO auvi = auvi1, auvi2
    CALL sync
    
    IF (par%master) THEN
      
      CALL write_to_debug_file
      
    END IF ! IF (par%master) THEN
    CALL sync
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_a_ex  )
    CALL deallocate_shared( wd_b_ex  )
    CALL deallocate_shared( wd_c_ex  )
    CALL deallocate_shared( wddx_a_ex)
    CALL deallocate_shared( wddx_b_ex)
    CALL deallocate_shared( wddx_c_ex)
    CALL deallocate_shared( wddy_a_ex)
    CALL deallocate_shared( wddy_b_ex)
    CALL deallocate_shared( wddy_c_ex)
    
    CALL deallocate_shared( wd_aa    )
    CALL deallocate_shared( wd_ab    )
    CALL deallocate_shared( wd_ac    )
    CALL deallocate_shared( wd_ba    )
    CALL deallocate_shared( wd_bb    )
    CALL deallocate_shared( wd_bc    )
    CALL deallocate_shared( wd_ca    )
    CALL deallocate_shared( wd_cb    )
    CALL deallocate_shared( wd_cc    )
    
    CALL deallocate_shared( wddx_aa  )
    CALL deallocate_shared( wddx_ab  )
    CALL deallocate_shared( wddx_ac  )
    CALL deallocate_shared( wddx_ba  )
    CALL deallocate_shared( wddx_bb  )
    CALL deallocate_shared( wddx_bc  )
    CALL deallocate_shared( wddx_ca  )
    CALL deallocate_shared( wddx_cb  )
    CALL deallocate_shared( wddx_cc  )
    
    CALL deallocate_shared( wddy_aa  )
    CALL deallocate_shared( wddy_ab  )
    CALL deallocate_shared( wddy_ac  )
    CALL deallocate_shared( wddy_ba  )
    CALL deallocate_shared( wddy_bb  )
    CALL deallocate_shared( wddy_bc  )
    CALL deallocate_shared( wddy_ca  )
    CALL deallocate_shared( wddy_cb  )
    CALL deallocate_shared( wddy_cc  )
    
    CALL deallocate_shared( wd_acau        )
    CALL deallocate_shared( wd_acav        )
    CALL deallocate_shared( wd_acau_acb    )
    CALL deallocate_shared( wd_acav_acb    )
    CALL deallocate_shared( wddx_acau_acb  )
    CALL deallocate_shared( wddx_acav_acb  )
    CALL deallocate_shared( wddy_acau_acb  )
    CALL deallocate_shared( wddy_acav_acb  )
    CALL deallocate_shared( wd_acb_acau    )
    CALL deallocate_shared( wd_acb_acav    )
    CALL deallocate_shared( wddx_acb_acau  )
    CALL deallocate_shared( wddx_acb_acav  )
    CALL deallocate_shared( wddy_acb_acau  )
    CALL deallocate_shared( wddy_acb_acav  )
        
    IF (par%master) WRITE(0,*) ' Finished testing the matrix operators!'
    CALL sync
    
  END SUBROUTINE test_matrix_operators
  SUBROUTINE test_matrix_operators_test_function( x, y, xmin, xmax, ymin, ymax, d, ddx, ddy)
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
    
    ! Local variables:
    REAL(dp)                                           :: xp, yp
    
    xp = 2._dp * pi * (x - xmin) / (xmax - xmin)
    yp = 2._dp * pi * (y - ymin) / (ymax - ymin)
    
    d   = SIN( xp) * SIN( yp)
    ddx = COS( xp) * SIN( yp) * 2._dp * pi / (xmax - xmin)
    ddy = SIN( xp) * COS( yp) * 2._dp * pi / (xmax - xmin)
        
  END SUBROUTINE test_matrix_operators_test_function
  
! == Solve a (modified version of) the Laplace equation on the mesh
  SUBROUTINE solve_modified_Laplace_equation( mesh)
    ! Test the discretisation by solving (a modified version of) the Laplace equation:
    !
    ! d/dx ( N * df/dx) + d/dx ( N * df/dy) = 0
    !
    ! (Modified to include a spatially variable stiffness, making it much more similar to
    ! the SSA/DIVA, and necessitating the use of a staggered grid/mesh to solve it)
    
    USE mesh_operators_module
    USE utilities_module
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, ti
    REAL(dp), DIMENSION(:    ), POINTER                ::  f_a,  N_b,  b_a
    INTEGER                                            :: wf_a, wN_b, wb_a
    TYPE(type_sparse_matrix_CSR)                       :: M_a, BC_a, mN_b, mNddx_b
    REAL(dp)                                           :: x, y, xp, yp
    REAL(dp), PARAMETER                                :: amp = 0._dp
    INTEGER                                            :: nnz_BC
        
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Solving the modified Laplace equation...'
    CALL sync
    
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
    
    ! Convert stiffness from vector to diagonal matrix
    CALL convert_vector_to_diag_CSR( N_b, mN_b)
    
    ! Create the operator matrix M
    CALL multiply_matrix_matrix_CSR( mN_b          , mesh%M_ddx_a_b, mNddx_b)
    CALL multiply_matrix_matrix_CSR( mesh%M_ddx_b_a, mNddx_b       , M_a    )
    
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
        
    IF (par%master) WRITE(0,*) ' Finished solving the modified Laplace equation!'
    CALL sync
    
  END SUBROUTINE solve_modified_Laplace_equation
  SUBROUTINE solve_modified_Laplace_equation_combi( mesh)
    ! Test the discretisation by solving (a modified version of) the Laplace equation:
    !
    ! d/dx ( N * df/dx) + d/dx ( N * df/dy) = 0
    !
    ! (Modified to include a spatially variable stiffness, making it much more similar to
    ! the SSA/DIVA, and necessitating the use of a staggered grid/mesh to solve it)
    
    USE mesh_operators_module
    USE utilities_module
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:
    INTEGER                                            :: avi, ati, vi, aci, edge_index
    REAL(dp), DIMENSION(:    ), POINTER                ::  f_aca,  N_acb,  b_aca
    INTEGER                                            :: wf_aca, wN_acb, wb_aca
    TYPE(type_sparse_matrix_CSR)                       :: M_aca, BC_aca, mN_acb, mNddx_acb
    REAL(dp)                                           :: x, y, xp, yp
    REAL(dp), PARAMETER                                :: amp = 0._dp
    INTEGER                                            :: nnz_BC
        
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Solving the modified Laplace equation on the combined AC-mesh...'
    CALL sync
    
  ! Solve the equation on the a-grid (vertex)
  ! =========================================
    
    ! Allocate shared memory
    IF (par%master) WRITE(0,*) '    Allocate shared memory for f, b, and N...'
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   f_aca, wf_aca)
    CALL allocate_shared_dp_1D( mesh%nVAaAc,   b_aca, wb_aca)
    CALL allocate_shared_dp_1D( mesh%nTriAaAc, N_acb, wN_acb )
    
    ! Set up the spatially variable stiffness N on the b-grid (triangle)
    IF (par%master) WRITE(0,*) '    Set up the spatially variable stiffness on the b-grid...'
    DO ati = mesh%ati1, mesh%ati2
      x = mesh%TriGCAaAc( ati,1)
      y = mesh%TriGCAaAc( ati,2)
      xp = (x - mesh%xmin) / (mesh%xmax - mesh%xmin)
      yp = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
      N_acb( ati) = 1._dp + EXP(amp * SIN( 2._dp * pi * xp) * SIN(2._dp * pi * yp))
    END DO
    CALL sync
    
    ! Convert stiffness from vector to diagonal matrix
    IF (par%master) WRITE(0,*) '    Convert stiffness from vector to diagonal matrix...'
    CALL convert_vector_to_diag_CSR( N_acb, mN_acb)
    
    ! Create the operator matrix M
    IF (par%master) WRITE(0,*) '    Multiply N with ddx...'
    CALL multiply_matrix_matrix_CSR( mN_acb            , mesh%M_ddx_aca_acb, mNddx_acb)
    IF (par%master) WRITE(0,*) '    Multiply with ddx again to find the operator matrix M...'
    CALL multiply_matrix_matrix_CSR( mesh%M_ddx_acb_aca, mNddx_acb          , M_aca   )
    CALL write_CSR_matrix_to_NetCDF( M_aca, 'M_aca.nc')
    
    
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
    CALL allocate_matrix_CSR_dist( BC_aca, mesh%nVAaAc, mesh%nVAaAc, nnz_BC)
    
    ! Fill in values
    IF (par%master) WRITE(0,*) '    Filling in BC matrix...'
    BC_aca%nnz = 0
    BC_aca%ptr = 1
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
        BC_aca%nnz  = BC_aca%nnz+1
        BC_aca%index( BC_aca%nnz) = avi
        BC_aca%val(   BC_aca%nnz) = 1._dp
        
        ! Right-hand side vector
        y = mesh%VAaAc( avi,2)
        y = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
        b_aca( avi) = 0.5_dp * (1._dp - COS( 2._dp * pi * y))
        
      ELSE
        ! Otherwise: zero
        
        ! Matrix
        BC_aca%nnz  = BC_aca%nnz+1
        BC_aca%index( BC_aca%nnz) = avi
        BC_aca%val(   BC_aca%nnz) = 1._dp
        
        ! Right-hand side vector
        b_aca( avi) = 0._dp
        
      END IF
      
      ! Finalise matrix row
      BC_aca%ptr( avi+1 : mesh%nVAaAc) = BC_aca%nnz + 1
      
    END DO ! DO avi = mesh%avi1, mesh%avi2
    CALL sync
    
    ! Finalise matrix
    CALL finalise_matrix_CSR_dist( BC_aca, mesh%avi1, mesh%avi2)
    CALL write_CSR_matrix_to_NetCDF( BC_aca, 'BC_aca.nc')
    
    ! Add boundary conditions to the operator matrix
    IF (par%master) WRITE(0,*) '    Add boundary conditions to the operator matrix...'
    CALL overwrite_rows_CSR( M_aca, BC_aca)
    CALL write_CSR_matrix_to_NetCDF( M_aca, 'M_aca_final.nc')
    
    ! Solve the equation
    IF (par%master) WRITE(0,*) '    Solve matrix equation...'
    CALL solve_matrix_equation_CSR( M_aca, b_aca, f_aca, C%DIVA_choice_matrix_solver, &
      SOR_nit = 5000, SOR_tol = 0.0001_dp, SOR_omega = 1.7_dp, &
      PETSc_rtol = C%DIVA_PETSc_rtol, PETSc_abstol = C%DIVA_PETSc_abstol)
      
    ! Write result to debug file
    IF (par%master) WRITE(0,*) '    Map solution to debug field...'
    IF (par%master) THEN
      debug%dp_2D_ac_03 = f_aca
      CALL write_to_debug_file
    END IF
    CALL sync
        
    IF (par%master) WRITE(0,*) ' Finished solving the modified Laplace equation on the combined AC-mesh!!'
    CALL sync
    
  END SUBROUTINE solve_modified_Laplace_equation_combi

END MODULE tests_and_checks_module
