MODULE mesh_ArakawaC_module

  ! Routines used in creating the Arakawa C staggered mesh.

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
  USE data_types_module,               ONLY: type_mesh
  USE mesh_help_functions_module,      ONLY: is_boundary_segment

  IMPLICIT NONE

  CONTAINS
    
  SUBROUTINE make_Ac_mesh( mesh)
    ! Determine the Arakawa C mesh vertices and their neighbour functions.
    ! For a detailed explanation, read the documentation (not kidding, its in there).    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    INTEGER                                       :: vi, ci, vj, cj, iti, ti, n1, n2, n3, vl, vr
    INTEGER,  DIMENSION(:,:), ALLOCATABLE         :: Aci_temp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: VAc_temp
                
    ! Allocate memory for the mapping arrays. Aci is stored in temporary memory on the master process first,
    ! once we know how much is needed, we will allocate shared memory and transfer the data there.
    CALL allocate_shared_int_2D( mesh%nV,     mesh%nC_mem, mesh%iAci,            mesh%wiAci          )
    CALL allocate_shared_int_0D(                           mesh%nAC,             mesh%wnAC           )
    
    ! Go through all vertex connections, determine the location of the Ac vertex (on the connection midpoint),
    ! the indices of the left and right Aa vertices, and the neighbour functions on the Ac vertex
    ! Only done by master
    ! ====================================
    
    vr = 0
    vl = 0
    
    IF (.NOT. par%master) THEN
    
      ! Allocate a little memory to prevent compiler warnings
      ALLOCATE(VAc_temp( 1,1))
      ALLOCATE(Aci_temp( 1,1))
    
    ELSE
    
      ! Allocate temporary memory for VAc, Aci and neighbour functions
      ALLOCATE(VAc_temp( mesh%nTri * 3,2))
      ALLOCATE(Aci_temp( mesh%nTri * 3,4))
      
      mesh%nAc   = 0
      mesh%iAci  = 0
      Aci_temp   = 0
      VAc_temp   = 0._dp
      
      ! Go over all Aa vertices
      DO vi = 1, mesh%nV
      
        ! Go over all the vertex connections
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          
          ! Skip connections that were already considered in the opposite direction
          IF (mesh%iAci( vi,ci)>0) CYCLE
          
          ! List Ac vertex number in reference array
          mesh%nAc = mesh%nAc + 1
          mesh%iAci( vi,ci) = mesh%nAc
          VAc_temp(mesh%nAc,:) = (mesh%V( vi,:) + mesh%V(vj,:)) / 2._dp
            
          ! Find reverse connection
          DO cj = 1, mesh%nC( vj)
            IF (mesh%C( vj,cj) == vi) THEN
              mesh%iAci( vj,cj) = mesh%nAC
              EXIT
            END IF
          END DO
          
          IF (.NOT. is_boundary_segment( mesh, vi, vj)) THEN
              
            ! Find indices of vl and vr
            DO iti = 1, mesh%niTri( vi)
              ti = mesh%iTri( vi,iti)
              DO n1 = 1, 3
                n2 = n1+1
                IF (n2==4) n2=1
                n3 = n2+1
                IF (n3==4) n3=1
                
                IF (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj) THEN
                  vl = mesh%Tri( ti,n3)
                ELSE IF (mesh%Tri( ti,n1) == vj .AND. mesh%Tri( ti,n2) == vi) THEN
                  vr = mesh%Tri( ti,n3)
                END IF
              END DO
            END DO
            
            ! List relevant Aa vertex indices vor Ac vertex in reference array
            Aci_temp( mesh%nAc,:) = [vi, vj, vl, vr]
            
          ELSE
            ! Surface slope on an edge segment Ac vertex equals that of the (single) triangle
              
            ! Find index of vl
            DO iti = 1, mesh%niTri( vi)
              ti = mesh%iTri( vi,iti)
              DO n1 = 1, 3
                n2 = n1+1
                IF (n2==4) n2=1
                n3 = n2+1
                IF (n3==4) n3=1
                
                IF ((mesh%Tri( ti,n1)==vi .AND. mesh%Tri( ti,n2)==vj) .OR. (mesh%Tri( ti,n1)==vj .AND. mesh%Tri( ti,n2)==vi)) THEN
                  vl = mesh%Tri( ti,n3)
                END IF
              END DO
            END DO
            
            ! List relevant Aa vertex indices vor Ac vertex in reference array
            ! NOTE: fourth entry is set to 1 because it cant be zero, but the fourth
            !       neighbour function is set to zero, so it doesn't matter...
            Aci_temp( mesh%nAc,:) = [vi, vj, vl, 1]          
          
          END IF
            
        END DO ! DO ci = 1, mesh%nC(vi)
      END DO ! DO vi = 1, mesh%nV
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate shared memory, move data there
    ! =======================================
    
    CALL allocate_shared_dp_2D(  mesh%nAc, 2, mesh%VAc, mesh%wVAc)
    CALL allocate_shared_int_2D( mesh%nAc, 4, mesh%Aci, mesh%wAci)
    
    IF (par%master) THEN
      mesh%VAc   = VAc_temp(   1:mesh%nAc,:)
      mesh%Aci   = Aci_temp(   1:mesh%nAc,:)
    END IF
    CALl sync
    
    ! Deallocate temporary memory
    ! ===========================
    
    DEALLOCATE( VAc_temp)
    DEALLOCATE( Aci_temp)
    
    ! Determine Ac vertex domains
    CALL partition_list( mesh%nAc, par%i, par%n, mesh%ac1, mesh%ac2)
    
    ! Find Ac edge indices
    CALL find_Ac_edge_indices( mesh)
    
    ! Create the combined Aa/Ac mesh
    CALL make_combined_AaAc_mesh( mesh)   
    
  END SUBROUTINE make_Ac_mesh
  SUBROUTINE find_Ac_edge_indices( mesh)
    ! Find the edge indices of the Ac vertices
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    INTEGER                                       :: aci, vi, vj
    
    ! Allocate shared memory
    CALl allocate_shared_int_1D( mesh%nAc, mesh%edge_index_Ac, mesh%wedge_index_Ac)
    
    DO aci = mesh%ac1, mesh%ac2
  
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      
      IF     ((mesh%edge_index( vi)==8 .OR. mesh%edge_index( vi) == 1 .OR. mesh%edge_index( vi) == 2) .AND. &
              (mesh%edge_index( vj)==8 .OR. mesh%edge_index( vj) == 1 .OR. mesh%edge_index( vj) == 2)) THEN
        ! North
        mesh%edge_index_Ac( aci) = 1
      ELSEIF ((mesh%edge_index( vi)==2 .OR. mesh%edge_index( vi) == 3 .OR. mesh%edge_index( vi) == 4) .AND. &
              (mesh%edge_index( vj)==2 .OR. mesh%edge_index( vj) == 3 .OR. mesh%edge_index( vj) == 4)) THEN
        ! East
        mesh%edge_index_Ac( aci) = 3
      ELSEIF ((mesh%edge_index( vi)==4 .OR. mesh%edge_index( vi) == 5 .OR. mesh%edge_index( vi) == 6) .AND. &
              (mesh%edge_index( vj)==4 .OR. mesh%edge_index( vj) == 5 .OR. mesh%edge_index( vj) == 6)) THEN
        ! East
        mesh%edge_index_Ac( aci) = 5
      ELSEIF ((mesh%edge_index( vi)==6 .OR. mesh%edge_index( vi) == 7 .OR. mesh%edge_index( vi) == 8) .AND. &
              (mesh%edge_index( vj)==6 .OR. mesh%edge_index( vj) == 7 .OR. mesh%edge_index( vj) == 8)) THEN
        ! East
        mesh%edge_index_Ac( aci) = 7
      ELSE
        ! Not an edge Ac vertex
        mesh%edge_index_Ac( aci) = 0
      END IF
      
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync
    
  END SUBROUTINE find_Ac_edge_indices
  SUBROUTINE make_combined_AaAc_mesh( mesh)   
    ! Create connectivity list and neigbour functions for the combined Aa-Ac mesh
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    INTEGER                                       :: vi, vj, ci, aci, ai, vk, vl, vr, aci1, aci2, aci3, aci4
    LOGICAL                                       :: switchthem
    INTEGER                                       :: ti, vip, viq, vir, aci_pq, aci_qr, aci_rp, acj, ack, iti
    INTEGER                                       :: n
    
    ! Allocate shared memory
    CALL allocate_shared_int_0D( mesh%nVAaAc,   mesh%wnVAaAc  )
    CALL allocate_shared_int_0D( mesh%nTriAaAc, mesh%wnTriAaAc)
    
    IF (par%master) THEN
      mesh%nVAaAc   = mesh%nV + mesh%nAc
      mesh%nTriAaAc = mesh%nTri + SUM(mesh%niTri)
    END IF
    CALL sync
    
    CALL allocate_shared_dp_2D(  mesh%nVAaAc,   2,             mesh%VAaAc,    mesh%wVAaAc   )
    CALL allocate_shared_int_1D( mesh%nVAaAc,                  mesh%nCAaAc,   mesh%wnCAaAc  )
    CALL allocate_shared_int_2D( mesh%nVAaAc,   mesh%nC_mem,   mesh%CAaAc,    mesh%wCAaAc   )
    CALL allocate_shared_int_2D( mesh%nTriAaAc, 3,             mesh%TriAaAc,  mesh%wTriAaAc )
    
    ! List coordinates and connectivity for the regular (Aa) vertices
    DO vi = mesh%v1, mesh%v2
      ai = vi
      mesh%VAaAc(  ai,:) = mesh%V(  vi,:)
      mesh%nCAaAc( ai  ) = mesh%nC( vi)
      DO ci = 1, mesh%nC( vi)
        aci = mesh%iAci( vi,ci)
        mesh%CAaAc( ai,ci) = aci + mesh%nV
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! List coordinates and connectivity for the staggered (Ac) vertices
    DO aci = mesh%ac1, mesh%ac2
    
      ai = aci + mesh%nV
      
      mesh%VAaAc( ai,:) = mesh%VAc( aci,:)
      
      IF (mesh%edge_index_Ac( aci) > 0) THEN
        ! Edge Ac vertex has only 1 triangle, so 4 neighbours: [vi, aci1, aci2, vj]
        
        ! The two regular vertices spanning this Ac vertex
        vi = mesh%Aci( aci,1)
        vj = mesh%Aci( aci,2)
        
        ! Sort them anticlockwise
        switchthem = .FALSE.
        IF     (mesh%edge_index_Ac( aci) == 1) THEN
          IF (mesh%V( vi,1) > mesh%V( vj,1)) switchthem = .TRUE.
        ELSEIF (mesh%edge_index_Ac( aci) == 3) THEN
          IF (mesh%V( vi,2) < mesh%V( vj,2)) switchthem = .TRUE.
        ELSEIF (mesh%edge_index_Ac( aci) == 5) THEN
          IF (mesh%V( vi,1) < mesh%V( vj,1)) switchthem = .TRUE.
        ELSEIF (mesh%edge_index_Ac( aci) == 7) THEN
          IF (mesh%V( vi,2) > mesh%V( vj,2)) switchthem = .TRUE.
        ELSE
          WRITE(0,*) 'make_combined_AaAc_mesh - ERROR: unknown Ac edge index'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        IF (switchthem) THEN
          vi = vi + vj
          vj = vi - vj
          vi = vi - vj
        END IF
        
        ! The third vertex is already listed in the Aci table
        vk = mesh%Aci( aci,3)
        
        ! Find the two Ac vertices
        aci1 = 0
        aci2 = 0
        DO ci = 1, mesh%nC( vk)
          IF     (mesh%C( vk,ci) == vi) THEN
            aci1 = mesh%iAci( vk,ci)
          ELSEIF (mesh%C( vk,ci) == vj) THEN
            aci2 = mesh%iAci( vk,ci)
          END IF
        END DO
        IF (aci1 == 0 .OR. aci2 == 0) THEN
          WRITE(0,*) 'make_combined_AaAc_mesh - ERROR: couldnt find Ac vertex'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
        ! The Ac vertex aci is surrounded CCW by [vi, aci1, aci2, vj]
        mesh%nCAaAc( ai     ) = 4
        mesh%CAaAc(  ai, 1:4) = [vi, aci1+mesh%nV, aci2+mesh%nV, vj]
        
      ELSE ! IF (mesh%edge_index_Ac( aci) > 0) THEN
        ! Non-edge Ac vertex has 2 triangles, so 6 neighbours: [vi, aci1, aci2, vj, aci3, aci4]
        
        ! The four regular vertices
        vi = mesh%Aci( aci,1)
        vj = mesh%Aci( aci,2)
        vl = mesh%Aci( aci,3)
        vr = mesh%Aci( aci,4)
        
        ! Find the four surrounding aci vertices
        aci1 = 0
        aci2 = 0
        aci3 = 0
        aci4 = 0
        DO ci = 1, mesh%nC( vr)
          IF     (mesh%C( vr,ci) == vi) THEN
            aci1 = mesh%iAci( vr,ci)
          ELSEIF (mesh%C( vr,ci) == vj) THEN
            aci2 = mesh%iAci( vr,ci)
          END IF
        END DO
        DO ci = 1, mesh%nC( vl)
          IF     (mesh%C( vl,ci) == vj) THEN
            aci3 = mesh%iAci( vl,ci)
          ELSEIF (mesh%C( vl,ci) == vi) THEN
            aci4 = mesh%iAci( vl,ci)
          END IF
        END DO
        IF (aci1 == 0 .OR. aci2 == 0 .OR. aci3 == 0 .OR. aci4 == 0) THEN
          WRITE(0,*) 'make_combined_AaAc_mesh - ERROR: couldnt find Ac vertex'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
        ! The Ac vertex aci is surrounded CCW by vi, aci1, aci2, vj, aci3, aci4]
        mesh%nCAaAc( ai     ) = 6
        mesh%CAaAc(  ai, 1:6) = [vi, aci1+mesh%nV, aci2+mesh%nV, vj, aci3+mesh%nV, aci4+mesh%nV]
        
      END IF ! IF (mesh%edge_index_Ac( aci) > 0) THEN
      
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync

    ! Triangles
    DO ti = mesh%t1, mesh%t2
    
      vip = mesh%Tri( ti,1)
      viq = mesh%Tri( ti,2)
      vir = mesh%Tri( ti,3)
      
      aci_pq = 0
      aci_qr = 0
      aci_rp = 0
      DO ci = 1, mesh%nC( vip)
        IF (mesh%C( vip,ci) == viq) THEN
          aci_pq = mesh%iAci( vip,ci)
          EXIT
        END IF
      END DO
      DO ci = 1, mesh%nC( viq)
        IF (mesh%C( viq,ci) == vir) THEN
          aci_qr = mesh%iAci( viq,ci)
          EXIT
        END IF
      END DO
      DO ci = 1, mesh%nC( vir)
        IF (mesh%C( vir,ci) == vip) THEN
          aci_rp = mesh%iAci( vir,ci)
          EXIT
        END IF
      END DO
      
      mesh%TriAaAc( ti,:) = [aci_pq, aci_qr, aci_rp] + mesh%nV
      
    END DO ! DO ti = mesh%t1, mesh%t2
    
    IF (par%master) THEN
    
      ti = mesh%nTri
      DO vi = 1, mesh%nV
        DO iti = 1, mesh%niTri( vi)
          
          vj = 0
          vk = 0
          DO n = 1, 3
            IF (mesh%Tri( mesh%iTri( vi,iti),n) /= vi) THEN
              IF (vj == 0) THEN
                vj = mesh%Tri( mesh%iTri( vi,iti),n)
              ELSE
                vk = mesh%Tri( mesh%iTri( vi,iti),n)
              END IF
            END IF
          END DO
          
          acj = 0
          ack = 0
          DO ci = 1, mesh%nC( vi)
            IF (mesh%C( vi,ci) == vj) THEN
              acj = mesh%iAci( vi,ci)
            ELSEIF (mesh%C( vi,ci) == vk) THEN
              ack = mesh%iAci( vi,ci)
            END IF
          END DO
          
          ti = ti + 1
          mesh%TriAaAc( ti,:) = [vi, acj+mesh%nV, ack+mesh%nV]
          
        END DO ! DO iti = 1, mesh%niTri( vi)
      END DO ! DO vi = 1, mesh%nV

    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Determine process domains
    CALL partition_list( mesh%nVAaAc, par%i, par%n, mesh%a1, mesh%a2)
    
  END SUBROUTINE make_combined_AaAc_mesh

END MODULE mesh_ArakawaC_module