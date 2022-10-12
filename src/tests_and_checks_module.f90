MODULE tests_and_checks_module

  ! A bunch of tests and checks of the CSR matrix operators

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
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
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

  ! Import specific functionality
  USE data_types_module
  USE mesh_operators_module
  USE petsc_module
  USE petscksp
  USE sparse_matrix_module
  USE netcdf_debug_module, ONLY: debug, write_to_debug_file
  USE netcdf_input_module, ONLY: read_field_from_file_2D, read_field_from_file_2D_monthly, read_field_from_file_3D

  IMPLICIT NONE

CONTAINS

! ===== Test matrix stuff =====
! =============================

  SUBROUTINE run_all_matrix_tests( mesh)
    ! Test all the matrix operators

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_all_matrix_tests'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) '========================================'
    IF (par%master) WRITE(0,*) '=== Testing all the matrix operators ==='
    IF (par%master) WRITE(0,*) '========================================'
    IF (par%master) WRITE(0,*) ''

    CALL test_petsc_to_csr

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_matrix_tests

! == Test petsc-to-CSR conversion

  SUBROUTINE test_petsc_to_csr
    ! Test petsc-to-CSR conversion

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_petsc_to_csr'
    TYPE(type_sparse_matrix_CSR_dp)                    :: A_CSR
    TYPE(tMat)                                         :: A
    INTEGER                                            :: nrows, ncols, nnz
    REAL(dp), DIMENSION(:    ), POINTER                ::  xx,  yy
    INTEGER                                            :: wxx, wyy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up CSR matrix
    nrows = 5
    ncols = 5
    nnz   = 10

    CALL allocate_matrix_CSR_shared( A_CSR, nrows, ncols, nnz)

    IF (par%master) THEN
      A_CSR%nnz   = 10
      A_CSR%ptr   = [1,3,5,7,10,11]
      A_CSR%index = [1,3,2,4,1,5,1,2,4,5]
      A_CSR%val   = [1._dp,2._dp,3._dp,4._dp,5._dp,6._dp,7._dp,8._dp,9._dp,10._dp]
    END IF
    CALL sync

    ! Convert to PETSc format
    CALL mat_CSR2petsc( A_CSR, A)

    ! Text multiplication
    CALL allocate_shared_dp_1D( 5, xx, wxx)
    CALL allocate_shared_dp_1D( 5, yy, wyy)
    IF (par%master) THEN
      xx = [1._dp, 2._dp, 3._dp, 4._dp, 5._dp]
    END IF
    CALL sync

    CALL multiply_PETSc_matrix_with_vector_1D( A, xx, yy)

    ! Check the result
    IF (yy(1) /= 7._dp .OR. yy(2) /= 22._dp .OR. yy(3) /= 35._dp .OR. yy(4) /= 59._dp .OR. yy(5) /= 50._dp) THEN
      WRITE(0,*) 'test_petsc_to_csr - ERROR: did not find the correct answer!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_petsc_to_csr

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_matrix_operators_abc_2D'
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

    ! Add routine to path
    CALL init_routine( routine_name)

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
    CALL sync

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
!      !! 2nd-order
!      !debug%dp_2D_a_01 = d2dx2_a_to_a
!      !debug%dp_2D_a_02 = d2dxdy_a_to_a
!      !debug%dp_2D_a_03 = d2dy2_a_to_a
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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

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
    ! d2dx2_a_to_a_3D
    ! d2dxdy_a_to_a_3D
    ! d2dy2_a_to_a_3D

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_matrix_operators_abc_3D'
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

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_matrix_operators_abc_3D

  ! 2nd-order accurate matrix operators on the b-grid
  SUBROUTINE test_matrix_operators_2nd_order_b_to_b_2D( mesh)
    ! Test all the 2nd-order accurate matrix operators on the b-grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'test_matrix_operators_2nd_order_b_to_b_2D'
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

    ! Add routine to path
    CALL init_routine( routine_name)

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

    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M2_ddx_b_b   , d_b_ex, ddx_b   )
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M2_ddy_b_b   , d_b_ex, ddy_b   )
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M2_d2dx2_b_b , d_b_ex, d2dx2_b )
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M2_d2dxdy_b_b, d_b_ex, d2dxdy_b)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M2_d2dy2_b_b , d_b_ex, d2dy2_b )

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

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
    !   d/dx ( N * df/dx) + d/dx ( N * df/dy) = 0
    !
    ! (Modified to include a spatially variable stiffness, making it much more similar to
    ! the SSA/DIVA, and necessitating the use of a staggered grid/mesh to solve it)
    !
    ! Using the product rule, this expression can be expanded to read:
    !
    !   N * d2f/dx2 + dN/dx * df/dx + N * d2f/dy2 + dN/dy * df/dy = 0

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_modified_Laplace_equation_b'
    INTEGER                                            :: vi, ti
    REAL(dp)                                           :: x, y, xp, yp
    REAL(dp), PARAMETER                                :: amp = 7._dp
    REAL(dp), DIMENSION(:    ), POINTER                ::  N_a,  N_b,  dN_dx_b,  dN_dy_b
    INTEGER                                            :: wN_a, wN_b, wdN_dx_b, wdN_dy_b
    TYPE(tVec)                                         :: vN_b, vdN_dx_b, vdN_dy_b
    REAL(dp), DIMENSION(:    ), POINTER                ::  x_b,  b_b
    INTEGER                                            :: wx_b, wb_b
    TYPE(tMat)                                         :: M1, M2, M3, M4
    TYPE(tMat)                                         :: A_b
    INTEGER                                            :: nTri_BC
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            ::  Tri_BC_loc
    INTEGER,  DIMENSION(:    ), POINTER                ::  Tri_BC
    INTEGER                                            :: wTri_BC
    REAL(dp)                                           :: rtol, abstol

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Convert them to PETSc vectors
    CALL vec_double2petsc( N_b    , vN_b    )
    CALL vec_double2petsc( dN_dx_b, vdN_dx_b)
    CALL vec_double2petsc( dN_dy_b, vdN_dy_b)

  ! Create the stiffness matrix A
  ! =============================

    ! Initialise with known non-zero structure
    CALL MatDuplicate( mesh%M2_ddx_b_b, MAT_SHARE_NONZERO_PATTERN, A_b, perr)

    ! Add the different terms

    ! M1 = N * d2/dx2
    CALL MatDuplicate( mesh%M2_d2dx2_b_b, MAT_COPY_VALUES, M1, perr)
    CALL MatDiagonalScale( M1, vN_b    , PETSC_NULL_VEC, perr)

    ! M2 = N * d2/dy2
    CALL MatDuplicate( mesh%M2_d2dy2_b_b, MAT_COPY_VALUES, M2, perr)
    CALL MatDiagonalScale( M2, vN_b    , PETSC_NULL_VEC, perr)

    ! M3 = dN/dx * d/dx
    CALL MatDuplicate( mesh%M2_ddx_b_b  , MAT_COPY_VALUES, M3, perr)
    CALL MatDiagonalScale( M3, vdN_dx_b, PETSC_NULL_VEC, perr)

    ! M4 = dN/dy * d/dy
    CALL MatDuplicate( mesh%M2_ddy_b_b  , MAT_COPY_VALUES, M4, perr)
    CALL MatDiagonalScale( M4, vdN_dy_b, PETSC_NULL_VEC, perr)

    ! Add them to A
    CALL MatAXPY( A_b, 1._dp, M1, SAME_NONZERO_PATTERN, perr)
    CALL MatAXPY( A_b, 1._dp, M2, SAME_NONZERO_PATTERN, perr)
    CALL MatAXPY( A_b, 1._dp, M3, SAME_NONZERO_PATTERN, perr)
    CALL MatAXPY( A_b, 1._dp, M4, SAME_NONZERO_PATTERN, perr)

  ! Add boundary conditions
  ! =======================

    ! List all triangles where boundary conditions should be applied
    ALLOCATE( Tri_BC_loc( mesh%nTri))
    IF (par%master) THEN
      nTri_BC = 0
      DO ti = 1, mesh%nTri
        IF (mesh%Tri_edge_index( ti) > 0) THEN
          nTri_BC = nTri_BC + 1
          Tri_BC_loc( nTri_BC) = ti
        END IF
      END DO
    END IF
    CALL MPI_BCAST( nTri_BC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL allocate_shared_int_1D( nTri_BC, Tri_BC, wTri_BC)
    IF (par%master) THEN
      Tri_BC = Tri_BC_loc( 1:nTri_BC)
    END IF
    DEALLOCATE( Tri_BC_loc)
    CALL sync

    ! For those triangles, set matrix row to zeros and diagonal element to one
    IF (par%master) Tri_BC = Tri_BC - 1 ! Because PETSc indexes from 0...
    CALL sync
    CALL MatZeroRows( A_b, nTri_BC, Tri_BC, 1._dp, PETSC_NULL_VEC, PETSC_NULL_VEC, perr)
    CALL deallocate_shared( wTri_BC)

    ! Fill in right-hand side and initial guess
    DO ti = mesh%ti1, mesh%ti2

      IF (mesh%Tri_edge_index( ti) == 0) THEN
        ! Free triangle

        b_b( ti) = 0._dp
        x_b( ti) = 0._dp

      ELSE ! IF (mesh%Tri_edge_index( ti) == 0) THEN
        ! Border triangle

        IF (mesh%Tri_edge_index( ti) == 6 .OR. mesh%Tri_edge_index( ti) == 7 .OR. mesh%Tri_edge_index( ti) == 8) THEN
          ! Western border: bump
          y = mesh%TriGC( ti,2)
          y = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
          b_b( ti) = 0.5_dp * (1._dp - COS( 2._dp * pi * y))
          x_b( ti) = b_b( ti)
        ELSE
          ! Other borders: zero
          b_b( ti) = 0._dp
          x_b( ti) = 0._dp
        END IF

      END IF ! IF (mesh%Tri_edge_index( ti) == 0) THEN

    END DO
    CALL sync

    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.

    CALL MatAssemblyBegin( A_b, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   A_b, MAT_FINAL_ASSEMBLY, perr)

    ! Solve the equation
    rtol   = 0.001_dp
    abstol = 0.001_dp
    CALL solve_matrix_equation_PETSc( A_b, b_b, x_b, rtol, abstol)

    ! Write result to debug file
    IF (par%master) THEN
      debug%dp_2D_b_01 = x_b
      CALL write_to_debug_file
    END IF
    CALL sync

    ! Clean up after yourself
    CALL VecDestroy( vN_b    , perr)
    CALL VecDestroy( vdN_dx_b, perr)
    CALL VecDestroy( vdN_dy_b, perr)
    CALL MatDestroy( A_b     , perr)
    CALL deallocate_shared( wN_a)
    CALL deallocate_shared( wN_b)
    CALL deallocate_shared( wdN_dx_b)
    CALL deallocate_shared( wdN_dy_b)
    CALL deallocate_shared( wx_b)
    CALL deallocate_shared( wb_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_modified_Laplace_equation_b


  ! ===== Tests =====
  ! =================

  SUBROUTINE NetCDF_input_test( region)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_input_test'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Text flexible data reading
    CALL NetCDF_test_xy_2D(             region%mesh, region%name)
    CALL NetCDF_test_xy_2D_monthly(     region%mesh, region%name)
    CALL NetCDF_test_xy_3D(             region%mesh, region%name)
    CALL NetCDF_test_lonlat_2D(         region%mesh, region%name)
    CALL NetCDF_test_lonlat_2D_monthly( region%mesh, region%name)
    CALL NetCDF_test_lonlat_3D(         region%mesh, region%name)
    CALL NetCDF_test_mesh_2D(           region%mesh, region%name)
    CALL NetCDF_test_mesh_2D_monthly(   region%mesh, region%name)
    CALL NetCDF_test_mesh_3D(           region%mesh, region%name)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_input_test

  SUBROUTINE NetCDF_test_xy_2D( mesh, region_name)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_test_xy_2D'
    CHARACTER(LEN=256)                                 :: filename, field_name_options
    REAL(dp), DIMENSION(:    ), POINTER                ::  d
    INTEGER                                            :: wd

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, d, wd)

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hi_xy_2D.nc'
    field_name_options = 'default_options_Hi'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_01( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hi_xy_2D_xflip.nc'
    field_name_options = 'default_options_Hi'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_02( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hi_xy_2D_yflip.nc'
    field_name_options = 'default_options_Hi'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_03( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hi_xy_2D_xyflip.nc'
    field_name_options = 'default_options_Hi'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_04( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hi_xy_2D_yx.nc'
    field_name_options = 'default_options_Hi'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_05( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hi_xy_2D_yx_xflip.nc'
    field_name_options = 'default_options_Hi'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_06( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hi_xy_2D_yx_yflip.nc'
    field_name_options = 'default_options_Hi'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_07( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hi_xy_2D_yx_xyflip.nc'
    field_name_options = 'default_options_Hi'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_08( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    ! Clean up after yourself
    CALL deallocate_shared( wd)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_test_xy_2D

  SUBROUTINE NetCDF_test_xy_2D_monthly( mesh, region_name)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_test_xy_2D_monthly'
    CHARACTER(LEN=256)                                 :: filename, field_name_options
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d
    INTEGER                                            :: wd

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, 12, d, wd)

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_xy_2D_monthly.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_01( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_xy_2D_monthly_xflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_02( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_xy_2D_monthly_yflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_03( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_xy_2D_monthly_xyflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_04( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_xy_2D_monthly_yx.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_05( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_xy_2D_monthly_yx_xflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_06( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_xy_2D_monthly_yx_yflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_07( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_xy_2D_monthly_yx_xyflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_08( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    ! Clean up after yourself
    CALL deallocate_shared( wd)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_test_xy_2D_monthly

  SUBROUTINE NetCDF_test_xy_3D( mesh, region_name)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_test_xy_3D'
    CHARACTER(LEN=256)                                 :: filename, field_name_options
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d
    INTEGER                                            :: wd

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, d, wd)

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_xy_3D_zeta1.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_01( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_xy_3D_zeta2.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_02( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_xy_3D_zeta3.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_03( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_xy_3D_zeta4.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_04( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    ! Clean up after yourself
    CALL deallocate_shared( wd)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_test_xy_3D

  SUBROUTINE NetCDF_test_lonlat_2D( mesh, region_name)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_test_lonlat_2D'
    CHARACTER(LEN=256)                                 :: filename, field_name_options
    REAL(dp), DIMENSION(:    ), POINTER                ::  d
    INTEGER                                            :: wd

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, d, wd)

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_lonlat_2D.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_01( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_lonlat_2D_lonflip.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_02( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_lonlat_2D_latflip.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_03( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_lonlat_2D_lonlatflip.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_04( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_lonlat_2D_latlon.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_05( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_lonlat_2D_latlon_lonflip.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_06( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_lonlat_2D_latlon_latflip.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_07( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_lonlat_2D_latlon_lonlatflip.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_08( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_lonlat_2D_lonshift.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_09( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_lonlat_2D_lon180.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_10( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    ! Clean up after yourself
    CALL deallocate_shared( wd)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_test_lonlat_2D

  SUBROUTINE NetCDF_test_lonlat_2D_monthly( mesh, region_name)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_test_lonlat_2D_monthly'
    CHARACTER(LEN=256)                                 :: filename, field_name_options
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d
    INTEGER                                            :: wd

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, 12, d, wd)

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_lonlat_2D_monthly.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_01( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_lonlat_2D_monthly_lonflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_02( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_lonlat_2D_monthly_latflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_03( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_lonlat_2D_monthly_lonlatflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_04( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_lonlat_2D_monthly_latlon.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_05( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_lonlat_2D_monthly_latlon_lonflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_06( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_lonlat_2D_monthly_latlon_latflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_07( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_lonlat_2D_monthly_latlon_lonlatflip.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_08( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_lonlat_2D_monthly_lonshift.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_09( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_lonlat_2D_monthly_lon180.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_10( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    ! Clean up after yourself
    CALL deallocate_shared( wd)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_test_lonlat_2D_monthly

  SUBROUTINE NetCDF_test_lonlat_3D( mesh, region_name)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_test_lonlat_3D'
    CHARACTER(LEN=256)                                 :: filename, field_name_options
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d
    INTEGER                                            :: wd

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, d, wd)

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_xy_3D_zeta1.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_01( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_xy_3D_zeta2.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_02( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_xy_3D_zeta3.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_03( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_xy_3D_zeta4.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_04( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    ! Clean up after yourself
    CALL deallocate_shared( wd)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_test_lonlat_3D

  SUBROUTINE NetCDF_test_mesh_2D( mesh, region_name)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_test_mesh_2D'
    CHARACTER(LEN=256)                                 :: filename, field_name_options
    REAL(dp), DIMENSION(:    ), POINTER                ::  d
    INTEGER                                            :: wd

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, d, wd)

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Hs_mesh_2D.nc'
    field_name_options = 'default_options_Hs'
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_a_01( mesh%vi1:mesh%vi2) = d( mesh%vi1:mesh%vi2)
    CALL write_to_debug_file

    ! Clean up after yourself
    CALL deallocate_shared( wd)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_test_mesh_2D

  SUBROUTINE NetCDF_test_mesh_2D_monthly( mesh, region_name)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_test_mesh_2D_monthly'
    CHARACTER(LEN=256)                                 :: filename, field_name_options
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d
    INTEGER                                            :: wd

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, 12, d, wd)

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/T2m_mesh_2D_monthly.nc'
    field_name_options = 'T2m'
    CALL read_field_from_file_2D_monthly( filename, field_name_options, mesh, d, region_name)
    debug%dp_2D_monthly_a_01( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    ! Clean up after yourself
    CALL deallocate_shared( wd)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_test_mesh_2D_monthly

  SUBROUTINE NetCDF_test_mesh_3D( mesh, region_name)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_test_mesh_3D'
    CHARACTER(LEN=256)                                 :: filename, field_name_options
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d
    INTEGER                                            :: wd

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, d, wd)

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_mesh_3D_zeta1.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_01( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_mesh_3D_zeta2.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_02( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_mesh_3D_zeta3.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_03( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/Ti_mesh_3D_zeta4.nc'
    field_name_options = 'Ti'
    CALL read_field_from_file_3D( filename, field_name_options, mesh, d, region_name)
    debug%dp_3D_a_04( mesh%vi1:mesh%vi2,:) = d( mesh%vi1:mesh%vi2,:)
    CALL write_to_debug_file

    ! Clean up after yourself
    CALL deallocate_shared( wd)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_test_mesh_3D

END MODULE tests_and_checks_module
