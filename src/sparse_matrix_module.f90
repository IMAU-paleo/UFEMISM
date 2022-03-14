MODULE sparse_matrix_module

  ! Operations on sparse matrices in CSR format

  ! Import basic functionality
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
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
  USE parallel_module,                 ONLY: allocate_shared_dist_int_0D, allocate_shared_dist_dp_0D, &
                                             allocate_shared_dist_int_1D, allocate_shared_dist_dp_1D, &
                                             allocate_shared_dist_int_2D, allocate_shared_dist_dp_2D, &
                                             allocate_shared_dist_int_3D, allocate_shared_dist_dp_3D, &
                                             allocate_shared_dist_bool_1D, &
                                             adapt_shared_dist_int_1D,    adapt_shared_dist_dp_1D, &
                                             adapt_shared_dist_int_2D,    adapt_shared_dist_dp_2D, &
                                             adapt_shared_dist_int_3D,    adapt_shared_dist_dp_3D, &
                                             adapt_shared_dist_bool_1D
  USE data_types_module,               ONLY: type_mesh, type_grid, type_sparse_matrix_CSR_dp
  USE petsc_module,                    ONLY: solve_matrix_equation_CSR_PETSc

CONTAINS
  
! == Solve a matrix equation Ax = b, where A is provided as a sparse matrix in CSR format
  SUBROUTINE solve_matrix_equation_CSR( AA, b, x, choice_matrix_solver, SOR_nit, SOR_tol, SOR_omega, PETSc_rtol, PETSc_abstol, colour_v1, colour_v2, colour_vi)
    ! Solve the matrix equation Ax = b
    ! The matrix A is provided in Compressed Sparse Row (CSR) format
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: x
    CHARACTER(LEN=256),                  INTENT(IN)    :: choice_matrix_solver
    INTEGER,                             INTENT(IN)    :: SOR_nit
    REAL(dp),                            INTENT(IN)    :: SOR_tol
    REAL(dp),                            INTENT(IN)    :: SOR_omega
    REAL(dp),                            INTENT(IN)    :: PETSc_rtol
    REAL(dp),                            INTENT(IN)    :: PETSc_abstol
    INTEGER,  DIMENSION(5)    , OPTIONAL, INTENT(IN)   :: colour_v1, colour_v2
    INTEGER,  DIMENSION(:,:  ), OPTIONAL, INTENT(IN)   :: colour_vi
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_matrix_equation_CSR'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Safety
    IF (AA%m /= AA%n .OR. AA%m /= SIZE( b,1) .OR. AA%m /= SIZE( x,1)) THEN
      CALL crash('matrix sizes dont match!')
    END IF
    
    CALL check_for_NaN_dp_1D(  AA%val, 'AA%val')
    CALL check_for_NaN_dp_1D(  b,      'b'      )
    CALL check_for_NaN_dp_1D(  x,      'x'      )
    
    IF (choice_matrix_solver == 'SOR') THEN
      ! Use the old simple SOR solver
      
      CALL solve_matrix_equation_CSR_SOR( AA, b, x, SOR_nit, SOR_tol, SOR_omega, colour_v1, colour_v2, colour_vi)
      
    ELSEIF (choice_matrix_solver == 'PETSc') THEN
      ! Use the PETSc solver (much preferred, this is way faster and more stable!)
    
      CALL solve_matrix_equation_CSR_PETSc( AA, b, x, PETSc_rtol, PETSc_abstol)
    
    ELSE
      CALL crash('unknown choice_matrix_solver "' // TRIM( choice_matrix_solver) // '"!')
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE solve_matrix_equation_CSR
  SUBROUTINE solve_matrix_equation_CSR_SOR( AA, b, x, nit, tol, omega, colour_v1, colour_v2, colour_vi)
    ! Solve the matrix equation Ax = b using successive over-relaxation (SOR)
    ! The matrix A is provided in Compressed Sparse Row (CSR) format
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: x
    INTEGER,                             INTENT(IN)    :: nit
    REAL(dp),                            INTENT(IN)    :: tol
    REAL(dp),                            INTENT(IN)    :: omega
    INTEGER,  DIMENSION(5)    , OPTIONAL, INTENT(IN)   :: colour_v1, colour_v2
    INTEGER,  DIMENSION(:,:  ), OPTIONAL, INTENT(IN)   :: colour_vi
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_matrix_equation_CSR_SOR'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (PRESENT( colour_v1)) THEN
      ! Safety
      IF ((.NOT. PRESENT( colour_v2)) .OR. (.NOT. PRESENT( colour_vi))) THEN
        CALL crash('needs all three colour arguments!')
      END IF
      CALL solve_matrix_equation_CSR_SOR_coloured( AA, b, x, nit, tol, omega, colour_v1, colour_v2, colour_vi)
    ELSE
      CALL solve_matrix_equation_CSR_SOR_regular(  AA, b, x, nit, tol, omega)
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE solve_matrix_equation_CSR_SOR
  SUBROUTINE solve_matrix_equation_CSR_SOR_regular( AA, b, x, nit, tol, omega)
    ! Solve the matrix equation Ax = b using successive over-relaxation (SOR)
    ! The matrix A is provided in Compressed Sparse Row (CSR) format
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: x
    INTEGER,                             INTENT(IN)    :: nit
    REAL(dp),                            INTENT(IN)    :: tol
    REAL(dp),                            INTENT(IN)    :: omega
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_matrix_equation_CSR_SOR_regular'
    INTEGER                                            :: i,j,k,it,i1,i2
    REAL(dp)                                           :: lhs, res, res_max, omega_dyn
    REAL(dp), DIMENSION(:    ), POINTER                ::  diagA
    INTEGER                                            :: wdiagA
    LOGICAL                                            :: found_it
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( AA%m, diagA, wdiagA)
    
    ! Partition equations over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    ! Store the central diagonal separately for faster access
    DO i = i1, i2
      found_it = .FALSE.
      DO k = AA%ptr( i), AA%ptr( i+1)-1
        j = AA%index( k)
        IF (j == i) THEN
          diagA( i) = AA%val( k)
          found_it = .TRUE.
        END IF
      END DO
      IF (.NOT. found_it) THEN
        CALL crash('matrix is missing a diagonal element in row {int_01} !', int_01 = i)
      END IF
    END DO
    CALL sync
    
    ! Perform the successive over-relaxation
    omega_dyn = omega
    
    res_max = tol * 2._dp
    it = 0
    SOR_iterate: DO WHILE (res_max > tol .AND. it < nit)
      it = it+1
      
      res_max = 0._dp

      DO i = i1, i2
      
        lhs = 0._dp
        DO k = AA%ptr( i), AA%ptr( i+1)-1
          j = AA%index( k)
          lhs = lhs + AA%val( k) * x( j)
        END DO
        
        res = (lhs - b( i)) / diagA( i)
        res_max = MAX( res_max, ABS(res))
        
        x( i) = x( i) - omega_dyn * res
        
      END DO ! DO i = i1, i2
      CALL sync
      
      ! Check if we've reached a stable solution
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, res_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      
      !IF (par%master) WRITE(0,*) '      SOR iteration ', it, ': res_max = ', res_max
      
      IF (it > 100 .AND. res_max > 1E3_dp ) THEN
        
        ! Divergence detected - decrease omega, reset solution to zero, restart SOR.
        IF (par%master) CALL warning('divergence detected; decrease omega, reset solution to zero, restart SOR')
        omega_dyn = omega_dyn - 0.1_dp
        it = 0
        x( i1:i2) = 0._dp
        CALL sync
        
        IF (omega_dyn <= 0.1_dp) THEN
        END IF
      END IF
      
    END DO SOR_iterate
    
    ! Clean up after yourself
    CALL deallocate_shared( wdiagA)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE solve_matrix_equation_CSR_SOR_regular
  SUBROUTINE solve_matrix_equation_CSR_SOR_coloured( AA, b, x, nit, tol, omega, colour_v1, colour_v2, colour_vi)
    ! Solve the matrix equation Ax = b using successive over-relaxation (SOR)
    ! The matrix A is provided in Compressed Sparse Row (CSR) format
    ! 
    ! Use a five-colouring map for perfect parallelisation
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: x
    INTEGER,                             INTENT(IN)    :: nit
    REAL(dp),                            INTENT(IN)    :: tol
    REAL(dp),                            INTENT(IN)    :: omega
    INTEGER,  DIMENSION(5)    ,          INTENT(IN)    :: colour_v1, colour_v2
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: colour_vi
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_matrix_equation_CSR_SOR_coloured'
    INTEGER                                            :: fci,fcvi,i,j,k,it,i1,i2
    REAL(dp)                                           :: lhs, res, res_max, omega_dyn
    REAL(dp), DIMENSION(:    ), POINTER                ::  diagA
    INTEGER                                            :: wdiagA
    LOGICAL                                            :: found_it
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( AA%m, diagA, wdiagA)
    
    ! Partition equations over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    ! Store the central diagonal separately for faster access
    DO i = i1, i2
      found_it = .FALSE.
      DO k = AA%ptr( i), AA%ptr( i+1)-1
        j = AA%index( k)
        IF (j == i) THEN
          diagA( i) = AA%val( k)
          found_it = .TRUE.
        END IF
      END DO
      IF (.NOT. found_it) THEN
        CALL crash('matrix is missing a diagonal element!')
      END IF
    END DO
    CALL sync
    
    ! Perform the successive over-relaxation
    omega_dyn = omega
    
    res_max = tol * 2._dp
    it = 0
    SOR_iterate: DO WHILE (res_max > tol .AND. it < nit)
      it = it+1
      
      res_max = 0._dp

      DO fci = 1, 5
        DO fcvi = colour_v1( fci), colour_v2( fci)
        
          i = colour_vi( fcvi, fci)
      
          lhs = 0._dp
          DO k = AA%ptr( i), AA%ptr( i+1)-1
            j = AA%index( k)
            lhs = lhs + AA%val( k) * x( j)
          END DO
          
          res = (lhs - b( i)) / diagA( i)
          res_max = MAX( res_max, ABS(res))
          
          x( i) = x( i) - omega_dyn * res
        
        END DO ! DO fcvi = fcvi1, fcvi2
        CALL sync
      END DO ! DO fci = 1, 5
      
      ! Check if we've reached a stable solution
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, res_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      
      !IF (par%master) WRITE(0,*) '      SOR iteration ', it, ': res_max = ', res_max
      
      IF (it > 100 .AND. res_max > 1E3_dp ) THEN
        
        ! Divergence detected - decrease omega, reset solution to zero, restart SOR.
        IF (par%master) CALL warning('divergence detected; decrease omega, reset solution to zero, restart SOR')
        omega_dyn = omega_dyn - 0.1_dp
        it = 0
        x( i1:i2) = 0._dp
        CALL sync
        
        IF (omega_dyn <= 0.1_dp) THEN
          CALL crash('divergence detected even with extremely low relaxation parameter!')
        END IF
      END IF
      
    END DO SOR_iterate
    
    ! Clean up after yourself
    CALL deallocate_shared( wdiagA)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE solve_matrix_equation_CSR_SOR_coloured
  
! == Mathematical operations on sparse matrices in CSR format

  ! Multiplication: C = A*B
  SUBROUTINE multiply_matrix_matrix_CSR( AA, BB, CC, nz_template)
    ! Perform the matrix multiplication C = A*B, where all three
    ! matrices are provided in CSR format (parallelised)
    !
    ! NOTE: probably not very efficient, better to use PETSc in the future!
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),           INTENT(IN)    :: AA, BB
    TYPE(type_sparse_matrix_CSR_dp), OPTIONAL, INTENT(IN)    :: nz_template
    TYPE(type_sparse_matrix_CSR_dp),           INTENT(INOUT) :: CC
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'multiply_matrix_matrix_CSR'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Safety
    IF (AA%n /= BB%m) THEN
      CALL crash('sizes of A and B dont match!')
    END IF
    
    IF (PRESENT( nz_template)) THEN
      CALL multiply_matrix_matrix_CSR_template( AA, BB, CC, nz_template)
    ELSE
      CALL multiply_matrix_matrix_CSR_notemplate( AA, BB, CC)
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE multiply_matrix_matrix_CSR
  SUBROUTINE multiply_matrix_matrix_CSR_template( AA, BB, CC, nz_template)
    ! Perform the matrix multiplication C = A*B, where all three
    ! matrices are provided in CSR format (parallelised)
    ! 
    ! Only calculate entries in the provided non-zero-structure
    
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),           INTENT(IN)    :: AA, BB
    TYPE(type_sparse_matrix_CSR_dp), OPTIONAL, INTENT(IN)    :: nz_template
    TYPE(type_sparse_matrix_CSR_dp),           INTENT(INOUT) :: CC
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'multiply_matrix_matrix_CSR_template'
    TYPE(type_sparse_matrix_CSR_dp)                    :: BBT
    INTEGER                                            :: i1,i2,k1,k2,ic,kc1,kc2,kc,jc
    INTEGER                                            :: ka1 , ka2,  nnz_row_a,  ja1,  ja2
    INTEGER                                            :: kbt1, kbt2, nnz_row_bt, jbt1, jbt2
    INTEGER                                            :: ka, kbt, ja, jbt
    LOGICAL                                            :: finished
    REAL(dp)                                           :: Cij
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Safety
    IF (AA%n /= BB%m .OR. nz_template%m /= AA%m .OR. nz_template%n /= BB%n) THEN
      CALL crash('sizes of A and B dont match!')
    END IF
    
    ! Allocate shared distributed memory for C
    CALL allocate_matrix_CSR_shared( CC, AA%m, BB%n, nz_template%nnz)
    
    ! Copy non-zero structure from template
    CALL partition_list( CC%m           , par%i, par%n, i1, i2)
    CALL partition_list( nz_template%nnz, par%i, par%n, k1, k2)
    CC%ptr(   i1:i2) = nz_template%ptr(   i1:i2)
    CC%index( k1:k2) = nz_template%index( k1:k2)
    IF (par%master) THEN
      CC%ptr( CC%m+1) = nz_template%ptr( CC%m+1)
      CC%nnz = nz_template%nnz
    END IF
    
    ! Calculate the transpose BT of B
    CALL transpose_matrix_CSR( BB, BBT)
    
    ! Calculate entries in C
    DO ic = i1, i2
    
      kc1 = CC%ptr( ic)
      kc2 = CC%ptr( ic+1) - 1
    
      DO kc = kc1, kc2
        
        jc = CC%index( kc)
      
        Cij = 0._dp
        
        ka1 = AA%ptr( ic)
        ka2 = AA%ptr( ic+1) - 1
        nnz_row_a = ka2 + 1 - ka1
        
        kbt1 = BBT%ptr( jc)
        kbt2 = BBT%ptr( jc+1) - 1
        nnz_row_bt = kbt2 + 1 - kbt1
        
        IF (nnz_row_a == 0 .OR. nnz_row_bt == 0) THEN
          ! Either A nor BT has no entries in this row, so the dot product is zero
          CYCLE
        END IF
        
        ja1  = AA%index(  ka1)
        ja2  = AA%index(  ka2)
        jbt1 = BBT%index( kbt1)
        jbt2 = BBT%index( kbt2)
        
        IF (ja1 > jbt2 .OR. jbt1 > ja2) THEN
          ! No overlap in column ranges of A and BT, so the dot product is zero
          CYCLE
        END IF 
        
        ka        = ka1
        kbt       = kbt1
        finished  = .FALSE.
        
        DO WHILE (.NOT. finished)
          
          ja  = AA%index(  ka )
          jbt = BBT%index( kbt)
          
          IF (ja < jbt) THEN
            ka = ka + 1
          ELSEIF (jbt < ja) THEN
            kbt = kbt + 1
          ELSEIF (jbt == ja) THEN
            Cij = Cij + (AA%val( ka) * BBT%val( kbt))
            ka  = ka + 1
            kbt = kbt + 1
          END IF
          
          IF (ka > ka2 .OR. kbt > kbt2 .OR. ja > jbt2 .OR. jbt > ja2) finished = .TRUE.
          
        END DO ! DO WHILE ((.NOT. finished_a) .OR. (.NOT. finished_b))
        
        ! Write the result to C
        CC%val( kc) = Cij
        
      END DO ! DO kc = 1, CC%n
      
    END DO ! DO ic = i1, i2
    CALL sync
    
    ! Deallocate BT
    CALL deallocate_matrix_CSR( BBT)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE multiply_matrix_matrix_CSR_template
  SUBROUTINE multiply_matrix_matrix_CSR_notemplate( AA, BB, CC)
    ! Perform the matrix multiplication C = A*B, where all three
    ! matrices are provided in CSR format (parallelised)
    ! 
    ! When no non-zero-structure is supplied, perform a general (= slow) multiplication
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),           INTENT(IN)    :: AA, BB
    TYPE(type_sparse_matrix_CSR_dp),           INTENT(INOUT) :: CC
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'multiply_matrix_matrix_CSR_notemplate'
    INTEGER                                            :: nnz_max_loc
    TYPE(type_sparse_matrix_CSR_dp)                    :: BBT
    INTEGER                                            :: i1, i2, ic, jc
    INTEGER                                            :: ka1 , ka2 , nnz_row_a,  ja1,  ja2
    INTEGER                                            :: kbt1, kbt2, nnz_row_bt, jbt1, jbt2
    INTEGER                                            :: ka, kbt, ja, jbt
    LOGICAL                                            :: finished_a, finished_bt
    LOGICAL                                            :: is_nonzero
    REAL(dp)                                           :: Cij
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Safety
    IF (AA%n /= BB%m) THEN
      CALL crash('sizes of A and B dont match!')
    END IF
    
    ! Allocate shared distributed memory for C
    ! Assume number of non-zeros in A, B, and C will be approximately equal
    nnz_max_loc = CEILING( REAL(MAX( AA%nnz, BB%nnz),dp) / REAL( par%n,dp))
    CALL allocate_matrix_CSR_dist( CC, AA%m, BB%n, nnz_max_loc)
    
    ! Calculate the transpose BT of B
    CALL transpose_matrix_CSR( BB, BBT)
    
    ! Partition rows of over the processors
    CALL partition_list( CC%m, par%i, par%n, i1, i2)
    
    ! Initialise
    CC%ptr = 1
    CC%nnz = 0
    
    DO ic = i1, i2
    
      DO jc = 1, CC%n
      
        is_nonzero = .FALSE.
        Cij        = 0._dp
        
        ka1 = AA%ptr( ic)
        ka2 = AA%ptr( ic+1) - 1
        nnz_row_a = ka2 + 1 - ka1
        
        kbt1 = BBT%ptr( jc)
        kbt2 = BBT%ptr( jc+1) - 1
        nnz_row_bt = kbt2 + 1 - kbt1
        
        IF (nnz_row_a == 0 .OR. nnz_row_bt == 0) THEN
          ! Either A nor BT has no entries in this row, so the dot product is zero
          CYCLE
        END IF
        
        ja1  = AA%index(  ka1)
        ja2  = AA%index(  ka2)
        jbt1 = BBT%index( kbt1)
        jbt2 = BBT%index( kbt2)
        
        IF (ja1 > jbt2 .OR. jbt1 > ja2) THEN
          ! No overlap in column ranges of A and BT, so the dot product is zero
          CYCLE
        END IF 
        
        ka          = ka1
        kbt         = kbt1
        finished_a  = .FALSE.
        finished_bt = .FALSE.
        
        DO WHILE ((.NOT. finished_a) .AND. (.NOT. finished_bt))
          
          ja  = AA%index(  ka)
          jbt = BBT%index( kbt)
          
          IF (ja < jbt) THEN
            ka = ka + 1
          ELSEIF (jbt < ja) THEN
            kbt = kbt + 1
          ELSEIF (jbt == ja) THEN
            is_nonzero = .TRUE.
            Cij = Cij + (AA%val( ka) * BBT%val( kbt))
            ka  = ka + 1
            kbt = kbt + 1
          END IF
          
          IF (ka  > ka2  .OR. ja  > jbt2) finished_a  = .TRUE.
          IF (kbt > kbt2 .OR. jbt > ja2 ) finished_bt = .TRUE.
          
        END DO ! DO WHILE ((.NOT. finished_a) .OR. (.NOT. finished_b))
        
        ! If the dot product is non-zero, add the result to C
        IF (is_nonzero) CALL add_entry_CSR_dist( CC, ic, jc, Cij)
        
      END DO ! DO jc = 1, CC%n
      
    END DO ! DO ic = i1, i2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( CC, i1, i2)
    
    ! Deallocate BT
    CALL deallocate_matrix_CSR( BBT)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE multiply_matrix_matrix_CSR_notemplate
  
  ! Addition: C = A+B
  SUBROUTINE add_matrix_matrix_CSR( AA, BB, CC, alpha, beta, same_structure)
    ! Perform the matrix addition C = (alpha*A) + (beta*B), where all three
    ! matrices are provided in CSR format (parallelised)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA, BB
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: CC
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: alpha, beta
    LOGICAL,  OPTIONAL,                  INTENT(IN)    :: same_structure
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_matrix_matrix_CSR'
    LOGICAL                                            :: same_structurep
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Safety
    IF (AA%m /= BB%m .OR. AA%n /= BB%n) THEN
      CALL crash('sizes of A and B dont match!')
    END IF
    
    IF (PRESENT( same_structure)) THEN
      same_structurep = same_structure
    ELSE
      same_structurep = .FALSE.
    END IF
    
    IF (same_structurep) THEN
      CALL add_matrix_matrix_CSR_template(   AA, BB, CC, alpha, beta)
    ELSE
      CALL add_matrix_matrix_CSR_notemplate( AA, BB, CC, alpha, beta)
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE add_matrix_matrix_CSR
  SUBROUTINE add_matrix_matrix_CSR_template( AA, BB, CC, alpha, beta)
    ! Perform the matrix addition C = (alpha*A) + (beta*B), where all three
    ! matrices are provided in CSR format (parallelised)
    ! 
    ! NOTE: assumes A and B (and therefore C) have the same non-zero structure
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA, BB
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: CC
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: alpha, beta
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_matrix_matrix_CSR_template'
    REAL(dp)                                           :: alphap, betap
    INTEGER                                            :: k1,k2
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Safety
    IF (AA%m /= BB%m .OR. AA%n /= BB%n) THEN
      CALL crash('sizes of A and B dont match!')
    END IF
    IF (AA%nnz /= BB%nnz) THEN
      CALL crash('A and B dont have the same non-zero structure!')
    END IF
    
    ! If no coefficients are provided, assume unity
    IF (PRESENT( alpha)) THEN
      alphap = alpha
    ELSE
      alphap = 1._dp
    END IF
    IF (PRESENT( beta )) THEN
      betap  = beta
    ELSE
      betap  = 1._dp
    END IF
    
    ! Allocate shared memory for C
    CALL initialise_matrix_CSR_from_nz_template( CC, AA)
    
    ! Add data from A and B
    CALL partition_list( CC%nnz, par%i, par%n, k1, k2)
    CC%val( k1:k2) = alphap * AA%val( k1:k2) + betap * BB%val( k1:k2)
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE add_matrix_matrix_CSR_template
  SUBROUTINE add_matrix_matrix_CSR_notemplate( AA, BB, CC, alpha, beta)
    ! Perform the matrix addition C = (alpha*A) + (beta*B), where all three
    ! matrices are provided in CSR format (parallelised)
    ! 
    ! NOTE: assumes column entries and A and B are sorted ascending!
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA, BB
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: CC
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: alpha, beta
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_matrix_matrix_CSR_notemplate'
    REAL(dp)                                           :: alphap, betap
    INTEGER                                            :: i1, i2, ic
    INTEGER                                            :: ka1, ka2, nnz_row_a
    INTEGER                                            :: kb1, kb2, nnz_row_b
    INTEGER                                            :: kk, ka, kb, ja, jb, jc
    REAL(dp)                                           :: vc
    LOGICAL                                            :: finished_a, finished_b
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Safety
    IF (AA%m /= BB%m .OR. AA%n /= BB%n) THEN
      CALL crash('sizes of A and B dont match!')
    END IF
    
    ! If no coefficients are provided, assume unity
    IF (PRESENT( alpha)) THEN
      alphap = alpha
    ELSE
      alphap = 1._dp
    END IF
    IF (PRESENT( beta )) THEN
      betap  = beta
    ELSE
      betap  = 1._dp
    END IF
    
    ! Partition rows over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    ! Allocate distributed shared memory for C
    CALL allocate_matrix_CSR_dist( CC, AA%m, AA%n, AA%nnz + BB%nnz)
    
    ! Initialise
    CC%ptr = 1
    CC%nnz = 0
    
    DO ic = i1, i2
      
      ka1 = AA%ptr( ic)
      ka2 = AA%ptr( ic+1) - 1
      nnz_row_a = ka2 + 1 - ka1
      
      kb1 = BB%ptr( ic)
      kb2 = BB%ptr( ic+1) - 1
      nnz_row_b = kb2 + 1 - kb1
      
      IF (nnz_row_a == 0 .AND. nnz_row_b == 0) THEN
        ! Neither A nor B has entries in this row
        CYCLE
      ELSEIF (nnz_row_a == 0 .AND. nnz_row_b > 0) THEN
        ! A has no entries in this row, but B does; copy data from B
        DO kk = 1, nnz_row_b
          kb = kb1 + kk - 1
          jc = BB%index( kb)
          vc = BB%val(   kb) * betap
          CALL add_entry_CSR_dist( CC, ic, jc, vc)
        END DO
      ELSEIF (nnz_row_a > 0 .AND. nnz_row_b == 0) THEN
        ! B has no entries in this row, but A does; copy data from A
        DO kk = 1, nnz_row_a
          ka = ka1 + kk - 1
          jc = AA%index( ka)
          vc = AA%val(   ka) * alphap
          CALL add_entry_CSR_dist( CC, ic, jc, vc)
        END DO
      ELSE
        ! Both A and B have entries in this row
        
        ka = ka1
        kb = kb1
        finished_a = .FALSE.
        finished_b = .FALSE.
        
        DO WHILE ((.NOT. finished_a) .OR. (.NOT. finished_b))
          
          IF ((.NOT. finished_a) .AND. (.NOT. finished_b)) THEN
          
            ja = AA%index( ka)
            jb = BB%index( kb)
            
            IF (ja < jb) THEN
              jc = AA%index( ka)
              vc = AA%val(   ka) * alphap
              CALL add_entry_CSR_dist( CC, ic, jc, vc)
              ka = ka + 1
            ELSEIF (jb < ja) THEN
              jc = BB%index( kb)
              vc = BB%val(   kb) * betap
              CALL add_entry_CSR_dist( CC, ic, jc, vc)
              kb = kb + 1
            ELSEIF (jb == ja) THEN
              jc = AA%index( ka)
              vc = AA%val(   ka) * alphap + BB%val( kb) * betap
              CALL add_entry_CSR_dist( CC, ic, jc, vc)
              ka = ka + 1
              kb = kb + 1
            END IF
          
          ELSEIF (.NOT. finished_a) THEN
          
            jc = AA%index( ka)
            vc = AA%val(   ka) * alphap
            CALL add_entry_CSR_dist( CC, ic, jc, vc)
            ka = ka + 1
            
          ELSEIF (.NOT. finished_b) THEN
          
            jc = BB%index( kb)
            vc = BB%val(   kb) * betap
            CALL add_entry_CSR_dist( CC, ic, jc, vc)
            kb = kb + 1
            
          END IF ! IF ((.NOT. finished_a) .AND. (.NOT. finished_b)) THEN
          
          IF (ka == ka2 + 1) finished_a = .TRUE.
          IF (kb == kb2 + 1) finished_b = .TRUE.
          
        END DO ! DO WHILE ((.NOT. finished_a) .OR. (.NOT. finished_b))
        
      END IF ! IF (nnz_row_a == 0 .AND. nnz_row_b == 0) THEN
      
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( CC, i1, i2)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE add_matrix_matrix_CSR_notemplate
  
  ! Vector multiplication c = Ab
  SUBROUTINE multiply_matrix_vector_CSR( AA, BB, CC)
    ! Perform the matrix multiplication C = A*B, where 
    ! A is a CSR-format matrix, and B and C are regular vectors.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: BB
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: CC
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'multiply_matrix_vector_CSR'
    INTEGER                                            :: mb, mc, i1, i2, i, k, j
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    mb = SIZE( BB,1)
    mc = SIZE( CC,1)
    
    ! Safety
    IF (mb /= AA%n .OR. mc /= AA%m) THEN
      CALL crash('sizes of A, B, and C dont match!')
    END IF
    
    ! Partition rows over the processors
    CALL partition_list( mc, par%i, par%n, i1, i2)
    
    DO i = i1, i2
      CC( i) = 0._dp
      DO k = AA%ptr( i), AA%ptr( i+1)-1
        j = AA%index( k)
        CC( i) = CC( i) + AA%val( k) * BB( j)
      END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE multiply_matrix_vector_CSR
  SUBROUTINE multiply_matrix_vector_2D_CSR( AA, BB, CC)
    ! Perform the matrix multiplication C = A*B, where 
    ! A is a CSR-format matrix, and B and C are regular vectors.
    ! 
    ! NOTE: A = [m-by-n], B,C = [n-by-nz]
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: BB
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: CC
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'multiply_matrix_vector_2D_CSR'
    INTEGER                                            :: mb, mc, i1, i2, i, k, j
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    mb = SIZE( BB,1)
    mc = SIZE( CC,1)
    
    ! Safety
    IF (mb /= AA%n .OR. mc /= AA%m) THEN
      CALL crash('sizes of A, B, and C dont match!')
    END IF
    
    ! Partition rows over the processors
    CALL partition_list( mc, par%i, par%n, i1, i2)
    
    DO i = i1, i2
      CC( i,:) = 0._dp
      DO k = AA%ptr( i), AA%ptr( i+1)-1
        j = AA%index( k)
        CC( i,:) = CC( i,:) + AA%val( k) * BB( j,:)
      END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE multiply_matrix_vector_2D_CSR
  
  ! Transposition: AT = A^T
  SUBROUTINE transpose_matrix_CSR( AA, AAT)
    ! Calculate the transpose AT of the matrix A, where both
    ! matrices are provided in CSR format (parallelised)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: AAT
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'transpose_matrix_CSR'
    INTEGER                                            :: nnz_max_loc, i1, i2, i, ii, kk1, kk2, kk, jj
    
    ! Add routine to path
    CALL init_routine( routine_name)
        
    ! Allocate shared distributed memory for AAT
    nnz_max_loc = CEILING( REAL( AA%nnz,dp) / REAL( par%n,dp))
    CALL allocate_matrix_CSR_dist( AAT, AA%n, AA%m, nnz_max_loc)
    
    ! Partition rows of AT (so columns of A) over the processors
    CALL partition_list( AAT%m, par%i, par%n, i1, i2)
    
    DO i = i1, i2
      
      ! Find all entries in column i of A
      DO ii = 1, AA%m ! ii loops over all the rows of A
        kk1 = AA%ptr( ii)
        kk2 = AA%ptr( ii+1) - 1
        DO kk = kk1, kk2 ! kk loops over all the entries of A(ii,:)
          jj = AA%index( kk)
          IF (jj == i) THEN
            ! Found an entry in the i-th column, ii-th row of A; add it in the i-th row, ii-th column of AT
            CALL add_entry_CSR_dist( AAT, i, ii, AA%val( kk))
          END IF
        END DO
      END DO
      
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( AAT, i1, i2)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE transpose_matrix_CSR
  
  ! Sort the column entries in a CSR matrix in ascending order
  SUBROUTINE sort_columns_in_CSR( AA)
    ! Sort the columns in each row of CSR-formatted matrix A in ascending order
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: AA
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'sort_columns_in_CSR'
    INTEGER                                            :: i1,i2,i,k1,k2,k,kk,jk,jkk
    REAL(dp)                                           :: vk, vkk
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Partition rows over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    DO i = i1, i2
      
      k1 = AA%ptr( i)
      k2 = AA%ptr( i+1) - 1
      
      ! Sort the data
      DO k = k1, k2-1
        jk = AA%index( k)
        vk = AA%val(   k)
        DO kk = k+1, k2
          jkk = AA%index( kk)
          vkk = AA%val(   kk)
          IF (jkk < jk) THEN
            ! Switch columns
            AA%index( k ) = jkk
            AA%index( kk) = jk
            AA%val(   k ) = vkk
            AA%val(   kk) = vk
            jk = jkk
            vk = vkk
          END IF
        END DO
      END DO ! DO k = 1, nnz_row
      
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE sort_columns_in_CSR
  
! == Basic memory operations on sparse matrices in CSR format
  SUBROUTINE allocate_matrix_CSR_shared( AA, m, n, nnz_max)
    ! Allocate shared memory for a CSR-format sparse m-by-n matrix A
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: AA
    INTEGER,                             INTENT(IN)    :: m, n, nnz_max
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_matrix_CSR_shared'
    INTEGER                                            :: i1,i2
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    CALL allocate_shared_int_0D( AA%m,       AA%wm      )
    CALL allocate_shared_int_0D( AA%n,       AA%wn      )
    CALL allocate_shared_int_0D( AA%nnz_max, AA%wnnz_max)
    CALL allocate_shared_int_0D( AA%nnz,     AA%wnnz    )
    
    IF (par%master) THEN
      AA%m       = m
      AA%n       = n
      AA%nnz_max = nnz_max
      AA%nnz     = 0
    END IF
    CALL sync
    
    CALL allocate_shared_int_1D( AA%m+1,     AA%ptr,   AA%wptr  )
    CALL allocate_shared_int_1D( AA%nnz_max, AA%index, AA%windex)
    CALL allocate_shared_dp_1D(  AA%nnz_max, AA%val,   AA%wval  )
    
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    AA%ptr(i1:i2) = 1
    IF (par%master) AA%ptr( AA%m+1) = 1
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 7)
    
  END SUBROUTINE allocate_matrix_CSR_shared
  SUBROUTINE initialise_matrix_CSR_from_nz_template( AA, nz_template)
    ! Allocate shared memory and initialise row/column indices from a non-zero-structure template
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: AA
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: nz_template
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_matrix_CSR_from_nz_template'
    INTEGER                                            :: i1,i2,k1,k2
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate shared memory
    CALL allocate_matrix_CSR_shared( AA, nz_template%m, nz_template%n, nz_template%nnz)
    
    ! Copy row/column indices
    CALL partition_list( nz_template%m  , par%i, par%n, i1, i2)
    CALL partition_list( nz_template%nnz, par%i, par%n, k1, k2)
    
    AA%ptr(   i1:i2) = nz_template%ptr(   i1:i2)
    AA%index( k1:k2) = nz_template%index( k1:k2)
    AA%val(   k1:k2) = 0._dp
    
    IF (par%master) THEN
      AA%nnz = nz_template%nnz
      AA%ptr( AA%m+1) = AA%nnz+1
    END IF
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE initialise_matrix_CSR_from_nz_template
  SUBROUTINE deallocate_matrix_CSR( AA)
    ! Deallocate the memory used by the CSR-format sparse m-by-n matrix A
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: AA
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_matrix_CSR'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    CALL deallocate_shared( AA%wm)
    CALL deallocate_shared( AA%wn)
    CALL deallocate_shared( AA%wnnz_max)
    CALL deallocate_shared( AA%wnnz)
    CALL deallocate_shared( AA%wptr)
    CALL deallocate_shared( AA%windex)
    CALL deallocate_shared( AA%wval)
    
    NULLIFY( AA%m      )
    NULLIFY( AA%n      )
    NULLIFY( AA%nnz_max)
    NULLIFY( AA%nnz    )
    NULLIFY( AA%ptr    )
    NULLIFY( AA%index  )
    NULLIFY( AA%val    )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE deallocate_matrix_CSR
  
  SUBROUTINE allocate_matrix_CSR_dist( AA, m, n, nnz_max_proc)
    ! Allocate shared memory for a CSR-format sparse m-by-n matrix A
    ! 
    ! NOTE: uses local allocatable memory, so that the processes can create their own
    !       CSR lists in parallel; use "finalise_matrix_CSR" after this is done!
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: AA
    INTEGER,                             INTENT(IN)    :: m, n, nnz_max_proc
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_matrix_CSR_dist'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Matrix dimensions should be the same everywhere, so use the shared version for that.
    CALL allocate_shared_int_0D( AA%m, AA%wm)
    CALL allocate_shared_int_0D( AA%n, AA%wn)
    
    IF (par%master) THEN
      AA%m       = m
      AA%n       = n
    END IF
    CALL sync

    ! Allocate local memory
    ALLOCATE( AA%nnz_max)
    ALLOCATE( AA%nnz    )
    
    AA%nnz_max = nnz_max_proc
    AA%nnz     = 0
    
    ALLOCATE( AA%ptr(   AA%m+1    ))
    ALLOCATE( AA%index( AA%nnz_max))
    ALLOCATE( AA%val(   AA%nnz_max))
    
    AA%ptr   = 1
    AA%index = 0
    AA%val   = 0._dp
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 2)
    
  END SUBROUTINE allocate_matrix_CSR_dist
  SUBROUTINE add_entry_CSR_dist( AA, i, j, v)
    ! Add value v to row i, column j of CSR-formatted matrix A
    ! 
    ! NOTE: assumes all rows before i are finished and nothing exists yet for rows after i!
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: AA
    INTEGER,                             INTENT(IN)    :: i,j
    REAL(dp),                            INTENT(IN)    :: v
    
    ! Increase number of non-zeros
    AA%nnz = AA%nnz + 1
    
    ! List entry
    AA%index( AA%nnz) = j
    AA%val(   AA%nnz) = v
    
    ! Update pointer list
    AA%ptr( i+1 : AA%m+1) = AA%nnz+1
    
    ! Extend memory if necessary
    IF (AA%nnz > AA%nnz_max - 1000) CALL extend_matrix_CSR_dist( AA, 1000)
    
  END SUBROUTINE add_entry_CSR_dist
  SUBROUTINE extend_matrix_CSR_dist( AA, nnz_extra)
    ! Extend shared memory for a CSR-format sparse m-by-n matrix A
    ! 
    ! NOTE: uses local allocatable memory, so that the processes can create their own
    !       CSR lists in parallel; use "finalise_matrix_CSR" after this is done!
    !
    ! NOTE: also means this doesn't need to be called by all processes at once.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: AA
    INTEGER,                             INTENT(IN)    :: nnz_extra
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extend_matrix_CSR_dist'
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: index_temp
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: val_temp
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate temporary memory
    ALLOCATE( index_temp( AA%nnz))
    ALLOCATE( val_temp(   AA%nnz))
    
    ! Copy data to temporary memory
    index_temp = AA%index( 1:AA%nnz)
    val_temp   = AA%val(   1:AA%nnz)
    
    ! Deallocate memory
    DEALLOCATE( AA%index)
    DEALLOCATE( AA%val  )
    
    ! Allocate new, extended memory
    AA%nnz_max = AA%nnz_max + nnz_extra
    ALLOCATE( AA%index( AA%nnz_max))
    ALLOCATE( AA%val(   AA%nnz_max))
    AA%index = 0
    AA%val   = 0._dp
    
    ! Copy data back from temporary memory
    AA%index( 1:AA%nnz) = index_temp
    AA%val(   1:AA%nnz) = val_temp
    
    ! Deallocate temporary memory
    DEALLOCATE( index_temp)
    DEALLOCATE( val_temp  )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE extend_matrix_CSR_dist
  SUBROUTINE finalise_matrix_CSR_dist( AA, i1, i2)
    ! Finalise a CSR-format sparse m-by-n matrix A from distributed to shared memory
    !
    ! NOTE: each process has data for rows i1-i2
    
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(INOUT) :: AA
    INTEGER,                             INTENT(IN)    :: i1, i2
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'finalise_matrix_CSR_dist'
    INTEGER                                            :: status(MPI_STATUS_SIZE)
    INTEGER                                            :: nnz_tot
    INTEGER,  DIMENSION(:    ), POINTER                :: ptr, index
    REAL(dp), DIMENSION(:    ), POINTER                :: val
    INTEGER                                            :: wptr, windex, wval
    INTEGER                                            :: p, k1, k2, k_loc, k_sum
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Determine total number of non-zero elements in A
    CALL MPI_ALLREDUCE( AA%nnz, nnz_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    ! Allocate shared memory
    CALL allocate_shared_int_1D( AA%m+1,  ptr,   wptr  )
    CALL allocate_shared_int_1D( nnz_tot, index, windex)
    CALL allocate_shared_dp_1D(  nnz_tot, val,   wval  )
    
    ! Determine range of indices for each process
    k1 = 0
    k2 = 0
    IF (par%master) THEN
      k1    = 1
      k2    = AA%nnz
      k_sum = AA%nnz
    END IF
    DO p = 1, par%n-1
      IF (par%master) THEN
        CALL MPI_SEND( k_sum, 1, MPI_INTEGER, p, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_RECV( k_loc, 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        k_sum = k_sum + k_loc
      ELSEIF (p == par%i) THEN
        CALL MPI_RECV( k_sum, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        CALL MPI_SEND( AA%nnz, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        k1 = k_sum + 1
        k2 = k_sum + AA%nnz
        AA%ptr = AA%ptr + k_sum
      END IF
      CALL sync
    END DO
    
    ! Copy data from process-local memory to shared memory
    ptr(   i1:i2) = AA%ptr(   i1:i2  )
    IF (par%master) ptr( AA%m+1) = nnz_tot+1
    index( k1:k2) = AA%index( 1:AA%nnz)
    val(   k1:k2) = AA%val(   1:AA%nnz)
    CALL sync
    
    ! Deallocate the process-local memory
    DEALLOCATE( AA%nnz_max)
    DEALLOCATE( AA%nnz    )
    DEALLOCATE( AA%ptr    )
    DEALLOCATE( AA%index  )
    DEALLOCATE( AA%val    )
    
    ! Allocate shared memory for the sparse matrix size
    CALL allocate_shared_int_0D( AA%nnz_max, AA%wnnz_max)
    CALL allocate_shared_int_0D( AA%nnz    , AA%wnnz    )
    
    IF (par%master) THEN
      AA%nnz_max = nnz_tot
      AA%nnz     = nnz_tot
    END IF
    CALL sync
    
    ! Associate the pointers in A with the newly allocated shared memory
    
    ! ptr
    AA%wptr = wptr
    CALL MPI_WIN_SHARED_QUERY( AA%wptr,   0, windowsize, disp_unit, baseptr, ierr)
    CALL C_F_POINTER( baseptr, AA%ptr,   [AA%m+1])
    ! index
    AA%windex = windex
    CALL MPI_WIN_SHARED_QUERY( AA%windex, 0, windowsize, disp_unit, baseptr, ierr)
    CALL C_F_POINTER( baseptr, AA%index, [AA%nnz])
    ! val
    AA%wval = wval
    CALL MPI_WIN_SHARED_QUERY( AA%wval,   0, windowsize, disp_unit, baseptr, ierr)
    CALL C_F_POINTER( baseptr, AA%val,   [AA%nnz])
    
    ! Sort the column entries in the assembled matrix
    CALL sort_columns_in_CSR( AA)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 5)
  
  END SUBROUTINE finalise_matrix_CSR_dist
  
! == Safety checks
  SUBROUTINE check_CSR( AA, matname)
    ! Check a CSR matrix for consistency
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA
    CHARACTER(LEN=*),                    INTENT(IN)    :: matname
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_CSR'
    INTEGER                                            :: i1,i2,i,k1,k2,k,j1,j2
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Dimensions
    IF (SIZE( AA%ptr,1)-1 /= AA%m) THEN
      CALL crash('matname = ' // TRIM( matname) // '; m = {int_01}, but SIZE( ptr,1) = {int_02}', int_01 = AA%m, int_02 = SIZE( AA%ptr,1))
    END IF
    IF (AA%ptr( AA%m+1) /= AA%nnz+1) THEN
      CALL crash('matname = ' // TRIM( matname) // '; nnz = {int_01}, but ptr( end) = {int_02}', int_01 = AA%nnz, int_02 = AA%ptr( AA%m+1))
    END IF
    IF (SIZE( AA%index,1) /= AA%nnz) THEN
      CALL crash('matname = ' // TRIM( matname) // '; nnz = {int_01}, but SIZE( index,1) = {int_02}', int_01 = AA%nnz, int_02 = SIZE( AA%index,1))
    END IF
    IF (SIZE( AA%val,1) /= AA%nnz) THEN
      CALL crash('matname = ' // TRIM( matname) // '; nnz = {int_01}, but SIZE( val,1) = {int_02}', int_01 = AA%nnz, int_02 = SIZE( AA%val,1))
    END IF
    
    ! Check if the ptr values are ascending
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    DO i = i1, i2
      IF (AA%ptr( i+1) < AA%ptr( i)) THEN
        CALL crash('matname = ' // TRIM( matname) // '; ptr not ascending! (row {int_01}; k1 = {int_02}, k2 = {int_03})', int_01 = i, int_02 = k1, int_03 = k2)
      END IF
    END DO
    
    ! Check if the columns of each row are listed in ascending order, and if they're not too large or too small
    DO i = i1, i2
      
      k1 = AA%ptr( i)
      k2 = AA%ptr( i+1)-1
      
      DO k = k1, k2-1
      
        j1 = AA%index( k)
        j2 = AA%index( k+1)
        
        IF (j2 <= j1) THEN
          CALL crash('matname = ' // TRIM( matname) // '; columns not sorted in ascending order!')
        END IF
        
        IF (j1 > AA%n .OR. j2 > AA%n) THEN
          CALL crash('matname = ' // TRIM( matname) // '; A%index({int_01}) = {int_02} > A%n = {int_03}!', int_01 = k, int_02 = j1, int_03 = AA%n)
        END IF
        
      END DO
      
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE check_CSR
  SUBROUTINE check_CSR_dist( AA, matname, i1, i2)
    ! Check a (distributed, local memory) CSR matrix for consistency
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA
    CHARACTER(LEN=*),                    INTENT(IN)    :: matname
    INTEGER,                             INTENT(IN)    :: i1,i2
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_CSR_dist'
    INTEGER                                            :: i,k1,k2,k,j
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Dimensions
    IF (SIZE( AA%ptr,1)-1 /= AA%m) THEN
      CALL crash('matname = ' // TRIM( matname) // '; m = {int_01}, but SIZE( ptr,1) = {int_02}', int_01 = AA%m, int_02 = SIZE( AA%ptr,1))
    END IF
    IF (AA%ptr( AA%m+1) /= AA%nnz+1) THEN
      CALL crash('matname = ' // TRIM( matname) // '; nnz = {int_01}, but ptr( end)) = {int_02}', int_01 = AA%nnz, int_02 = AA%ptr( AA%m+1))
    END IF
    
    ! Check if the ptr values are ascending
    DO i = i1, i2
      IF (AA%ptr( i+1) < AA%ptr( i)) THEN
        CALL crash('matname = ' // TRIM( matname) // '; ptr not ascending! (row {int_01}; k1 = {int_02}, k2 = {int_03})', int_01 = i, int_02 = k1, int_03 = k2)
      END IF
    END DO
    
    ! Check column entries
    DO i = i1, i2
      
      k1 = AA%ptr( i)
      k2 = AA%ptr( i+1)-1
      
      DO k = k1, k2
      
        j = AA%index( k)
        
        IF (j > AA%n) THEN
          CALL crash('matname = ' // TRIM( matname) // '; A%index({int_01}) = {int_02} > A%n = {int_03}!', int_01 = k, int_02 = j, int_03 = AA%n)
        END IF
        
        IF (j <= 0) THEN
          CALL crash('matname = ' // TRIM( matname) // '; A%index({int_01}) = {int_02}!', int_01 = k, int_02 = j)
        END IF
        
      END DO
      
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE check_CSR_dist
  SUBROUTINE are_identical_matrices_CSR( AA, BB, isso)
    ! Check if CSR-formatted matrices A and B are identical
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: AA, BB
    LOGICAL,                             INTENT(OUT)   :: isso
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'are_identical_matrices_CSR'
    INTEGER                                            :: i1, i2, k1, k2, i, k
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    isso = .TRUE.
    
    ! Simple dimension check
    IF (AA%m /= BB%m .OR. AA%n /= BB%n .OR. AA%nnz /= BB%nnz) THEN
      isso = .FALSE.
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! ptr
    CALL partition_list( AA%m  , par%i, par%n, i1, i2)
    DO i = i1, i2
      IF (AA%ptr( i) /= BB%ptr( i)) THEN
        isso = .FALSE.
        EXIT
      END IF
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, isso, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. isso) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! index, val
    CALL partition_list( AA%nnz, par%i, par%n, k1, k2)
    DO k = k1, k2
      IF (AA%index( k) /= BB%index( k) .OR. AA%val( k) /= BB%val( k)) THEN
        isso = .FALSE.
        EXIT
      END IF
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, isso, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. isso) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE are_identical_matrices_CSR

END MODULE sparse_matrix_module
