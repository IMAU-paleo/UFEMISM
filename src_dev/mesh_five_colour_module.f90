MODULE mesh_five_colour_module

  USE mpi
  USE configuration_module,          ONLY: dp, C
  USE parallel_module,               ONLY: par, sync, allocate_shared_int_0D, allocate_shared_dp_0D, &
                                                      allocate_shared_int_1D, allocate_shared_dp_1D, &
                                                      allocate_shared_int_2D, allocate_shared_dp_2D, &
                                                      allocate_shared_int_3D, allocate_shared_dp_3D, &
                                                      allocate_shared_bool_1D, deallocate_shared
  USE data_types_module,             ONLY: type_mesh
  USE mesh_help_functions_module,    ONLY: partition_list, write_mesh_to_text_file

  IMPLICIT NONE
  
  CONTAINS
  
  SUBROUTINE calculate_five_colouring_AaAc( mesh)
    ! Calculate a five-colouring of the combined Aa-Ac mesh (to be used in parallelised SOR)
    ! Based on Williams, M.H.: "A Linear Algorithm for Colouring
    ! Planar Graphs With Five Colours", The Computer Journal 28, 78-81, 1985.
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: deg
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE       :: L
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: Q4
    INTEGER                                       :: Q4n
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: Q5
    INTEGER                                       :: Q5n
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: LDP
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE       :: mark
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: S_vi
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE       :: S_L
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: S_u
    INTEGER                                       :: Sn
    INTEGER                                       :: noofvert
    
    INTEGER                                       :: vi, ci
    
    ! Allocate shared memory
    CALL allocate_shared_int_1D( mesh%nVAaAc,    mesh%colour,    mesh%wcolour   )
    CALL allocate_shared_int_2D( mesh%nVAaAc, 5, mesh%colour_vi, mesh%wcolour_vi)
    CALL allocate_shared_int_1D(              5, mesh%colour_nV, mesh%wcolour_nV)
    
    ! Not parallelised (but not needed, very fast)
    IF (.NOT.par%master) THEN
    
      ALLOCATE( deg(  1  ))
      ALLOCATE( L(    1,1))
      ALLOCATE( Q4(   1  ))
      ALLOCATE( Q5(   1  ))
      ALLOCATE( LDP(  1  ))
      ALLOCATE( mark( 1  ))
      ALLOCATE( S_vi( 1  ))
      ALLOCATE( S_L(  1,1))
      ALLOCATE( S_u(  1  ))
    
    ELSE
    
      ! DENK DROM
      CALL write_mesh_to_text_file( mesh, 'mesh_01.txt')
      
      ! Initialise
      ALLOCATE( deg(  mesh%nVAaAc             ))
      ALLOCATE( L(    mesh%nVAaAc, mesh%nC_mem))
      ALLOCATE( Q4(   mesh%nVAaAc             ))
      ALLOCATE( Q5(   mesh%nVAaAc             ))
      ALLOCATE( LDP(  mesh%nVAaAc             ))
      ALLOCATE( mark( mesh%nVAaAc             ))
      ALLOCATE( S_vi( mesh%nVAaAc             ))
      ALLOCATE( S_L(  mesh%nVAaAc, mesh%nC_mem))
      ALLOCATE( S_u(  mesh%nVAaAc             ))
      
      deg      = mesh%nCAaAc
      L        = mesh%CAaAc
      Q4       = 0
      Q4n      = 0
      Q5       = 0
      Q5n      = 0
      LDP      = 0
      mark     = .FALSE.
      S_vi     = 0
      S_L      = 0
      S_u      = 0
      Sn       = 0
      noofvert = mesh%nVAaAc
      
      ! Initialise Q4,Q5
      DO vi = 1, mesh%nVAaAc
        IF     (deg( vi) <= 4) THEN
          CALL add_to_Q4( deg, LDP, Q4, Q4n, Q5, Q5n, vi)
        ELSEIF (deg( vi) == 5) THEN
          CALL add_to_Q5( deg, LDP, Q4, Q4n, Q5, Q5n, vi)
        END IF
      END DO
      
      ! Reduce the mesh
      CALL reduce( mesh, deg, L, LDP, mark, Q4, Q4n, Q5, Q5n, S_vi, S_L, S_u, Sn, noofvert)
      
      ! Colour the mesh
      CALL colour( mesh, Q4, S_vi, S_L, S_u, Sn)
      
      ! Check if the solution is valid
      CALL check_solution( mesh)
      
      ! Create same-coloured parallelisation domains
      mesh%colour_vi = 0
      mesh%colour_nV = 0
      
      DO vi = 1, mesh%nVAaAc
        ci = mesh%colour( vi)
        mesh%colour_nV( ci) = mesh%colour_nV( ci) + 1
        mesh%colour_vi( mesh%colour_nV( ci), ci) = vi
      END DO
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Assign processor domains to the five colours
    DO ci = 1, 5
      CALL partition_list( mesh%colour_nV( ci), par%i, par%n, mesh%colour_v1( ci), mesh%colour_v2( ci))
    END DO
    
    ! Clean up after yourself
    DEALLOCATE( deg )
    DEALLOCATE( L   )
    DEALLOCATE( Q4  )
    DEALLOCATE( Q5  )
    DEALLOCATE( LDP )
    DEALLOCATE( mark)
    DEALLOCATE( S_vi)
    DEALLOCATE( S_L )
    DEALLOCATE( S_u )
    
  END SUBROUTINE calculate_five_colouring_AaAc
  
  SUBROUTINE reduce( mesh, deg, L, LDP, mark, Q4, Q4n, Q5, Q5n, S_vi, S_L, S_u, Sn, noofvert)
    ! Iteratively reduce the mesh until 5 vertices remain
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: deg
    INTEGER,  DIMENSION(:,:  ), INTENT(INOUT)     :: L
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: LDP
    LOGICAL,  DIMENSION(:    ), INTENT(INOUT)     :: mark
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q4
    INTEGER,                    INTENT(INOUT)     :: Q4n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q5
    INTEGER,                    INTENT(INOUT)     :: Q5n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: S_vi
    INTEGER,  DIMENSION(:,:  ), INTENT(INOUT)     :: S_L
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: S_u
    INTEGER,                    INTENT(INOUT)     :: Sn
    INTEGER,                    INTENT(INOUT)     :: noofvert
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: vi, k, ci1, xi, ci2, yi, ci3, qii, ai, it
    LOGICAL                                       :: found_pair, xy_are_adjacent
    
    it = 0
    DO WHILE (noofvert > 5)
      
      it = it + 1
      IF (it > mesh%nVAaAc*2) THEN
        WRITE(0,*) '  ERROR - reduce got stuck!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      IF (Q4n > 0) THEN
        ! DELETE( top entry from Q4)
        
        vi = Q4( Q4n)
        CALL delete( deg, L, LDP, Q4, Q4n, Q5, Q5n, S_vi, S_L, S_u, Sn, noofvert, vi)
        
      ELSE
        ! take top entry v from Q5
        vi = Q5( Q5n)
        
        ! Look DO pair of non-adjacent vertices x,y from N(v)
        found_pair = .FALSE.
        k = C%nconmax
        DO ci1 = 1, deg( vi)
          xi = L( vi,ci1)
          DO ci2 = ci1+1, deg( vi)
            yi = L( vi,ci2)
            
            ! Check IF they are adjacent
            xy_are_adjacent = .FALSE.
            DO ci3 = 1, deg( yi)
              IF (L( yi,ci3) == xi) THEN
                xy_are_adjacent = .TRUE.
                EXIT
              END IF
            END DO
            IF (xy_are_adjacent) CYCLE
            
            IF (deg( xi) < k .AND. deg( yi) < k) THEN
              found_pair = .TRUE.
              EXIT
            END IF
          END DO
          IF (found_pair) EXIT
        END DO ! DO ci1 = 1: deg( vi)
        
        ! IF (two non-adjacent vertices x,y from N(v) can be found such
        ! that DEG(x) < k and DEG(y) < k)
        !   DELETE(v)
        !   IDENTIFY(x,y)
        ! ELSE
        !   return v to rear of Q5
        ! END IF
        
        IF (found_pair) THEN
          CALL delete(   deg, L, LDP,       Q4, Q4n, Q5, Q5n, S_vi, S_L, S_u, Sn, noofvert, vi)
          CALL identify( deg, L, LDP, mark, Q4, Q4n, Q5, Q5n, S_vi, S_L, S_u, Sn, noofvert, xi, yi)
        ELSE
          ! Return v to the rear of Q5
          Q5( 2:Q5n) = Q5( 1:Q5n-1)
          Q5( 1) = vi
          ! Update LDP
          DO qii = 1, Q5n
            ai = Q5( qii)
            LDP( ai) = qii
          END DO
        END IF
        
      END IF ! IF (Q4n > 0)
      
      !WRITE(0,*) 'noofvert = ', noofvert, ', Q4n = ', Q4n, ', Q5n = ', Q5n, ', Sn = ', Sn

    END DO ! while (noofvert > 5)
    
  END SUBROUTINE reduce
  SUBROUTINE colour( mesh, Q4, S_vi, S_L, S_u, Sn)
    ! Colour the mesh by reversing the reduction
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: Q4
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: S_vi
    INTEGER,  DIMENSION(:,:  ), INTENT(IN)        :: S_L
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: S_u
    INTEGER,                    INTENT(INOUT)     :: Sn
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: vi, ui, col, ci, it
    INTEGER, DIMENSION(5)                         :: L_colours
    INTEGER, DIMENSION(:     ), ALLOCATABLE       :: Lv
    
    ALLOCATE( Lv( mesh%nC_mem))
    
    ! Colour vertices remaining in Q4
    mesh%colour( Q4(1)) = 1
    mesh%colour( Q4(2)) = 2
    mesh%colour( Q4(3)) = 3
    mesh%colour( Q4(4)) = 4
    mesh%colour( Q4(5)) = 5
    
    ! Colour vertices in stack
    it = 0
    DO WHILE (Sn > 0)
      
      it = it + 1
      IF (it > mesh%nVAaAc*2) THEN
        WRITE(0,*) '  ERROR - colour got stuck!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
      vi = S_vi( Sn  )
      Lv = S_L(  Sn,:)
      ui = S_u(  Sn  )
      
      IF (S_u( Sn) > 0) THEN
        ! Paint u and v with the same colour
        
        mesh%colour( ui) = mesh%colour( vi)
        
      ELSE ! if (S_u( Sn) > 0)
        
        ! Find what colours are already used in neighbourhood
        L_colours = 0
        DO ci = 1, mesh%nC_mem
          ui = Lv( ci)
          IF (ui == 0) EXIT
          col = mesh%colour( ui)
          IF (col == 0) CYCLE
          L_colours( col) = 1
        END DO
        IF (SUM(L_colours) > 4) THEN
          WRITE(0,*) '    colour - ERROR: SUM(L_colours) > 4'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
        ! Paint vi with unused colour
        DO ci = 1, 5
          IF (L_colours( ci) == 0) THEN
            mesh%colour( vi) = ci
            EXIT
          END IF
        END DO
        
      END IF ! if (S_u( Sn) > 0)
      
      Sn = Sn - 1
      
    END DO ! while (Sn > 0)
    
    DEALLOCATE( Lv)
    
  END SUBROUTINE colour
  
  SUBROUTINE check_solution( mesh)
    ! Check if the solution is a valid five-colouring
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: vi, ci, vj
    
    DO vi = 1, mesh%nVAaAc
      IF (mesh%colour( vi) == 0) THEN
        WRITE(0,*) '  ERROR - the resulting five-colouring is invalid!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      DO ci = 1, mesh%nCAaAc( vi)
        vj = mesh%CAaAc( vi,ci)
        IF (mesh%colour( vi) == mesh%colour(vj)) THEN
          WRITE(0,*) '  ERROR - the resulting five-colouring is invalid!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END DO
    END DO
    
  END SUBROUTINE check_solution
  
  SUBROUTINE delete( deg, L, LDP, Q4, Q4n, Q5, Q5n, S_vi, S_L, S_u, Sn, noofvert, vi)
    ! Delete vertex vi from the reduced mesh
    
    IMPLICIT NONE
  
    ! Input variables
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: deg
    INTEGER,  DIMENSION(:,:  ), INTENT(INOUT)     :: L
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: LDP
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q4
    INTEGER,                    INTENT(INOUT)     :: Q4n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q5
    INTEGER,                    INTENT(INOUT)     :: Q5n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: S_vi
    INTEGER,  DIMENSION(:,:  ), INTENT(INOUT)     :: S_L
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: S_u
    INTEGER,                    INTENT(INOUT)     :: Sn
    INTEGER,                    INTENT(INOUT)     :: noofvert
    INTEGER,                    INTENT(IN)        :: vi
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: ci, wi, ci2
    
    ci = S_u( 1)
    
    ! for all w in L(v): remove v from L(w)
    DO ci = 1, deg( vi)
      wi = L( vi,ci)
      DO ci2 = 1, deg( wi)
        IF (L( wi,ci2) == vi) THEN
          L( wi,ci2:deg( wi)-1) = L( wi,ci2+1:deg( wi))
          L( wi, deg( wi)) = 0
          deg( wi) = deg( wi) - 1
          CALL check( deg, LDP, Q4, Q4n, Q5, Q5n, wi)
          EXIT
        END IF
      END DO
    END DO
    
    ! push v and pointer to L(v) onto S
    Sn = Sn + 1
    S_vi( Sn  ) = vi
    S_L(  Sn,:) = L( vi,:)
    
    ! remove v from either Q4 or Q5
    IF     (deg( vi) <= 4) THEN
      ! v is in Q4 remove it
      CALL remove_from_Q4( LDP, Q4, Q4n, Q5, Q5n, vi)
    ELSEIF (deg( vi) == 5) THEN
      ! v is in Q5 remove it
      CALL remove_from_Q5( LDP, Q4, Q4n, Q5, Q5n, vi)
    ELSE
      ! We're not allowed to delete vertices with higher degrees than this!
      WRITE(0,*) '    delete - ERROR: deg > 5'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    noofvert = noofvert - 1
    
  END SUBROUTINE delete
  SUBROUTINE identify( deg, L, LDP, mark, Q4, Q4n, Q5, Q5n, S_vi, S_L, S_u, Sn, noofvert, ui, vi)
    ! "identify" (merge) vertices ui and vi in the reduced mesh
    
    IMPLICIT NONE
  
    ! Input variables
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: deg
    INTEGER,  DIMENSION(:,:  ), INTENT(INOUT)     :: L
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: LDP
    LOGICAL,  DIMENSION(:    ), INTENT(INOUT)     :: mark
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q4
    INTEGER,                    INTENT(INOUT)     :: Q4n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q5
    INTEGER,                    INTENT(INOUT)     :: Q5n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: S_vi
    INTEGER,  DIMENSION(:,:  ), INTENT(INOUT)     :: S_L
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: S_u
    INTEGER,                    INTENT(INOUT)     :: Sn
    INTEGER,                    INTENT(INOUT)     :: noofvert
    INTEGER,                    INTENT(IN)        :: ui
    INTEGER,                    INTENT(IN)        :: vi
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: ci, wi, ciu, cii
    
    WRITE(0,*) 'beep'
    
    ci = S_L(1,1)
    
    ! mark neighbourhoog of v
    DO ci = 1, deg( vi)
      wi = L( vi,ci)
      mark( wi) = .TRUE.
    END DO
    
    ! disconnect neighbours of u from u, connect them to v instead
    DO ci = 1, deg( ui)
      wi = L( ui,ci)
      
      ! delete u from L(w)
      ciu = 0
      DO cii = 1, deg( wi)
        IF (L( wi,cii) == ui) THEN
          ciu = cii
          EXIT
        END IF
      END DO
      IF (ciu == 0) WRITE(0,*) '    identify - ERROR: ciu = 0'; CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      L( wi,ciu:deg( wi)-1) = L( wi,ciu+1:deg( wi))
      L( wi,deg( wi)) = 0
      deg( wi) = deg( wi) - 1
      CALL check( deg, LDP, Q4, Q4n, Q5, Q5n, wi)
      
      IF (mark( wi)) THEN
        ! w is already connected to v
      ELSE
        ! w is not yet connected to v
        
        ! add w to L(v)
        deg( vi) = deg( vi) + 1
        L( vi, deg( vi)) = wi
        CALL check( deg, LDP, Q4, Q4n, Q5, Q5n, vi)
        
        ! add v to L(w)
        deg( wi) = deg( wi) + 1
        L( wi, deg( wi)) = vi
        CALL check( deg, LDP, Q4, Q4n, Q5, Q5n, wi)
        
      END IF ! IF (~mark( wi))
    
    END DO ! FOR (w in L(u)) DO
    
    ! FOR (w in L(v)) MARK(w) = false
    DO ci = 1, deg( vi)
      wi = L( vi,ci)
      mark( wi) = .FALSE.
    END DO
    
    ! delete u from Qi
    IF     (deg( ui) <= 4) THEN
      ! remove u from Q4
      CALL remove_from_Q4( LDP, Q4, Q4n, Q5, Q5n, ui)
    ELSEIF (deg( ui) == 5) THEN
      ! remove u from Q5
      CALL remove_from_Q4( LDP, Q4, Q4n, Q5, Q5n, ui)
    ELSE
      WRITE(0,*) '    identify - ERROR: deg(ui) > 5'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! push u,v onto S
    Sn = Sn + 1
    S_vi( Sn) = ui
    S_u(  Sn) = vi
    
    noofvert = noofvert - 1
    
  END SUBROUTINE identify
  SUBROUTINE check( deg, LDP, Q4, Q4n, Q5, Q5n, vi)
    ! Update status of vi in Q4,Q5
    
    IMPLICIT NONE
  
    ! Input variables
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: deg
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: LDP
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q4
    INTEGER,                    INTENT(INOUT)     :: Q4n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q5
    INTEGER,                    INTENT(INOUT)     :: Q5n
    INTEGER,                    INTENT(IN)        :: vi
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    
    IF     (deg( vi) <= 4) THEN
      ! vi should be in Q4
      IF (      is_in_Q5( LDP, Q4, Q5, vi)) CALL remove_from_Q5(      LDP, Q4, Q4n, Q5, Q5n, vi)
      IF (.NOT. is_in_Q4( LDP, Q4, Q5, vi)) CALL add_to_Q4(      deg, LDP, Q4, Q4n, Q5, Q5n, vi)
    ELSEIF (deg( vi) == 5) THEN
      ! vi should be in Q5
      IF (      is_in_Q4( LDP, Q4, Q5, vi)) CALL remove_from_Q4(      LDP, Q4, Q4n, Q5, Q5n, vi)
      IF (.NOT. is_in_Q5( LDP, Q4, Q5, vi)) CALL add_to_Q5(      deg, LDP, Q4, Q4n, Q5, Q5n, vi)
    ELSE
      ! vi should be in neither
      IF (      is_in_Q4( LDP, Q4, Q5, vi)) THEN
        WRITE(0,*) '    check - ERROR: deg > 5 and vi is in Q4'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      IF (      is_in_Q5( LDP, Q4, Q5, vi)) CALL remove_from_Q5(      LDP, Q4, Q4n, Q5, Q5n, vi)
    END IF
    
  END SUBROUTINE check
  
  SUBROUTINE add_to_Q4( deg, LDP, Q4, Q4n, Q5, Q5n, vi)
    ! Add vertex vi to Q4
    
    IMPLICIT NONE
  
    ! Input variables
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: deg
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: LDP
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q4
    INTEGER,                    INTENT(INOUT)     :: Q4n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q5
    INTEGER,                    INTENT(INOUT)     :: Q5n
    INTEGER,                    INTENT(IN)        :: vi
    
    ! Local variables
    INTEGER                                       :: cerr, ierr
    
    cerr = Q5( 1)
    cerr = Q5n
    
    IF (deg( vi) > 4) THEN
      WRITE(0,*) '    add_to_Q4 - ERROR: deg > 4'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    IF (LDP( vi) > 0) THEN
      WRITE(0,*) '    add_to_Q4 - ERROR: LDP > 0'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    Q4n = Q4n + 1
    Q4( Q4n) = vi
    LDP( vi) = Q4n
    
  END SUBROUTINE add_to_Q4
  SUBROUTINE add_to_Q5( deg, LDP, Q4, Q4n, Q5, Q5n, vi)
    ! Add vertex vi to Q5
    
    IMPLICIT NONE
  
    ! Input variables
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: deg
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: LDP
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q4
    INTEGER,                    INTENT(INOUT)     :: Q4n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q5
    INTEGER,                    INTENT(INOUT)     :: Q5n
    INTEGER,                    INTENT(IN)        :: vi
    
    ! Local variables
    INTEGER                                       :: cerr, ierr
    
    cerr = Q4(1)
    cerr = Q4n
    
    IF (deg( vi) /= 5) THEN
      WRITE(0,*) '    add_to_Q5 - ERROR: deg /= 5'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    IF (LDP( vi) > 0) THEN
      WRITE(0,*) '    add_to_Q5 - ERROR: LDP > 0'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    Q5n = Q5n + 1
    Q5( Q5n) = vi
    LDP( vi) = Q5n
    
  END SUBROUTINE add_to_Q5
  SUBROUTINE remove_from_Q4( LDP, Q4, Q4n, Q5, Q5n, vi)
    ! Remove vertex vi from Q4
    
    IMPLICIT NONE
  
    ! Input variables
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: LDP
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q4
    INTEGER,                    INTENT(INOUT)     :: Q4n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q5
    INTEGER,                    INTENT(INOUT)     :: Q5n
    INTEGER,                    INTENT(IN)        :: vi
    
    ! Local variables
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: qi, qii, ai
    
    ai = Q5n
    
    qi = LDP( vi)
    
    IF (qi == 0) THEN
      WRITE(0,*) '    remove_from_Q4 - ERROR: qi == 0'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    IF (Q5( qi) == vi) THEN
      WRITE(0,*) '    remove_from_Q4 - ERROR: vi in Q5'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    IF (Q4( qi) /= vi) THEN
      WRITE(0,*) '    remove_from_Q4 - ERROR: vi not in Q4'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    Q4( qi:Q4n-1) = Q4( qi+1:Q4n)
    Q4( Q4n) = 0
    Q4n = Q4n - 1
    LDP( vi) = 0
    
    DO qii = qi, Q4n
      ai = Q4( qii)
      LDP( ai) = qii
    END DO
    
  END SUBROUTINE remove_from_Q4
  SUBROUTINE remove_from_Q5( LDP, Q4, Q4n, Q5, Q5n, vi)
    ! Remove vertex vi from Q5
    
    IMPLICIT NONE
  
    ! Input variables
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: LDP
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q4
    INTEGER,                    INTENT(INOUT)     :: Q4n
    INTEGER,  DIMENSION(:    ), INTENT(INOUT)     :: Q5
    INTEGER,                    INTENT(INOUT)     :: Q5n
    INTEGER,                    INTENT(IN)        :: vi
    
    ! Local variables
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: qi, qii, ai
    
    ai = Q4n
    
    qi = LDP( vi)
    
    IF (qi == 0) THEN
      WRITE(0,*) '    remove_from_Q5 - ERROR: qi == 0'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    IF (Q4( qi) == vi) THEN
      WRITE(0,*) '    remove_from_Q5 - ERROR: vi in Q4'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    IF (Q5( qi) /= vi) THEN
      WRITE(0,*) '    remove_from_Q5 - ERROR: vi not in Q5'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    Q5( qi:Q5n-1) = Q5( qi+1:Q5n)
    Q5( Q5n) = 0
    Q5n = Q5n - 1
    LDP( vi) = 0
    
    DO qii = qi, Q5n
      ai = Q5( qii)
      LDP( ai) = qii
    END DO
    
  END SUBROUTINE remove_from_Q5
  FUNCTION is_in_Q4( LDP, Q4, Q5, vi) RESULT( isso)
    ! Check if vertex vi is in Q4
    
    IMPLICIT NONE
  
    ! Input variables
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: LDP
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: Q4
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: Q5
    INTEGER,                    INTENT(IN)        :: vi
    LOGICAL                                       :: isso
    
    ! Local variables
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: qi
    
    isso = .FALSE.
    qi = LDP( vi)
    IF (qi==0) RETURN
    IF (Q4( qi) == vi) isso = .TRUE.
    IF (Q5( qi) == vi .AND. isso) THEN
      WRITE(0,*) '    is_in_Q4 - ERROR: vi is in both Q4 and Q5'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END FUNCTION is_in_Q4
  FUNCTION is_in_Q5( LDP, Q4, Q5, vi) RESULT( isso)
    ! Check if vertex vi is in Q5
    
    IMPLICIT NONE
  
    ! Input variables
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: LDP
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: Q4
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: Q5
    INTEGER,                    INTENT(IN)        :: vi
    LOGICAL                                       :: isso
    
    ! Local variables
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: qi
    
    isso = .FALSE.
    qi = LDP( vi)
    IF (qi==0) RETURN
    IF (Q5( qi) == vi) isso = .TRUE.
    IF (Q4( qi) == vi .AND. isso) THEN
      WRITE(0,*) '    is_in_Q5 - ERROR: vi is in both Q4 and Q5'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END FUNCTION is_in_Q5

END MODULE mesh_five_colour_module
