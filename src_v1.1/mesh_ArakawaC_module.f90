MODULE mesh_ArakawaC_module
  ! Routines used in creating the Arakawa C staggered mesh.

  USE mpi
  USE configuration_module,          ONLY: dp, C
  USE parallel_module,               ONLY: par, sync, allocate_shared_int_0D, allocate_shared_dp_0D, &
                                                      allocate_shared_int_1D, allocate_shared_dp_1D, &
                                                      allocate_shared_int_2D, allocate_shared_dp_2D, &
                                                      allocate_shared_int_3D, allocate_shared_dp_3D, &
                                                      allocate_shared_bool_1D, deallocate_shared, &
                                                      adapt_shared_int_1D,    adapt_shared_dp_1D, &
                                                      adapt_shared_int_2D,    adapt_shared_dp_2D, &
                                                      adapt_shared_int_3D,    adapt_shared_dp_3D, &
                                                      adapt_shared_bool_1D               
  USE data_types_module,             ONLY: type_mesh
  USE mesh_help_functions_module,    ONLY: is_boundary_segment, partition_list
  USE mesh_derivatives_module,       ONLY: get_neighbour_functions_vertex_gr

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
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Nx_Ac_temp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Ny_Ac_temp
    REAL(dp), DIMENSION(:  ), ALLOCATABLE         :: Np_Ac_temp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: No_Ac_temp
    REAL(dp), DIMENSION(4)                        :: Nxl, Nxr, Nyl, Nyr
    REAL(dp)                                      :: Nzl, Nzr
    REAL(dp), DIMENSION(4)                        :: Nx, Ny
    REAL(dp)                                      :: Ux, Uy, U, Np
    REAL(dp), DIMENSION(4)                        :: No
                
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
      ALLOCATE(VAc_temp(   1,1))
      ALLOCATE(Aci_temp(   1,1))
      ALLOCATE(Nx_Ac_temp( 1,1))
      ALLOCATE(Ny_Ac_temp( 1,1))
      ALLOCATE(Np_Ac_temp( 1  ))
      ALLOCATE(No_Ac_temp( 1,1))
    
    ELSE
    
      ! Allocate temporary memory for VAc, Aci and neighbour functions
      ALLOCATE(VAc_temp(   mesh%nTri * 3,2))
      ALLOCATE(Aci_temp(   mesh%nTri * 3,4))
      ALLOCATE(Nx_Ac_temp( mesh%nTri * 3,4))
      ALLOCATE(Ny_Ac_temp( mesh%nTri * 3,4))
      ALLOCATE(Np_Ac_temp( mesh%nTri * 3  ))
      ALLOCATE(No_Ac_temp( mesh%nTri * 3,4))
      
      mesh%nAc   = 0
      mesh%iAci  = 0
      Aci_temp   = 0
      VAc_temp   = 0._dp
      Nx_Ac_temp = 0._dp
      Ny_Ac_temp = 0._dp
      Np_Ac_temp = 0._dp
      No_Ac_temp = 0._dp
      
      ! Go over all Aa vertices
      DO vi = 1, mesh%nV
      
        ! Go over all the vertex connections
        DO ci = 1, mesh%nC(vi)
          vj = mesh%C(vi,ci)
          
          ! Skip connections that were already considered in the opposite direction
          IF (mesh%iAci(vi,ci)>0) CYCLE
          
          ! List Ac vertex number in reference array
          mesh%nAc = mesh%nAc + 1
          mesh%iAci(vi,ci) = mesh%nAc
          VAc_temp(mesh%nAc,:) = (mesh%V(vi,:) + mesh%V(vj,:)) / 2._dp
            
          ! Find reverse connection
          DO cj = 1, mesh%nC(vj)
            IF (mesh%C(vj,cj)==vi) THEN
              mesh%iAci(vj,cj) = mesh%nAC
              EXIT
            END IF
          END DO
          
          IF (.NOT. is_boundary_segment( mesh, vi, vj)) THEN
              
            ! Find indices of vl and vr
            DO iti = 1, mesh%niTri(vi)
              ti = mesh%iTri(vi,iti)
              DO n1 = 1, 3
                n2 = n1+1
                IF (n2==4) n2=1
                n3 = n2+1
                IF (n3==4) n3=1
                
                IF (mesh%Tri(ti,n1)==vi .AND. mesh%Tri(ti,n2)==vj) THEN
                  vl = mesh%Tri(ti,n3)
                ELSE IF (mesh%Tri(ti,n1)==vj .AND. mesh%Tri(ti,n2)==vi) THEN
                  vr = mesh%Tri(ti,n3)
                END IF
              END DO
            END DO
            
            ! List relevant Aa vertex indices vor Ac vertex in reference array
            Aci_temp(mesh%nAc,:) = [vi, vj, vl, vr]
            
            ! Calculate neighbour functions - first on adjacent triangles, then on Ac vertex
            ! x an y just as average of the left and right triangles
            Nxl = [mesh%V(vl,2)-mesh%V(vj,2), mesh%V(vi,2)-mesh%V(vl,2), mesh%V(vj,2)-mesh%V(vi,2), 0._dp]
            Nyl = [mesh%V(vj,1)-mesh%V(vl,1), mesh%V(vl,1)-mesh%V(vi,1), mesh%V(vi,1)-mesh%V(vj,1), 0._dp]
            Nxr = [mesh%V(vj,2)-mesh%V(vr,2), mesh%V(vr,2)-mesh%V(vi,2), 0._dp, mesh%V(vi,2)-mesh%V(vj,2)]
            Nyr = [mesh%V(vr,1)-mesh%V(vj,1), mesh%V(vi,1)-mesh%V(vr,1), 0._dp, mesh%V(vj,1)-mesh%V(vi,1)]
            Nzl = ( (mesh%V(vj,1)-mesh%V(vi,1)) * (mesh%V(vl,2)-mesh%V(vi,2)) ) - &
                  ( (mesh%V(vj,2)-mesh%V(vi,2)) * (mesh%V(vl,1)-mesh%V(vi,1)) )
            Nzr = ( (mesh%V(vr,1)-mesh%V(vi,1)) * (mesh%V(vj,2)-mesh%V(vi,2)) ) - &
                  ( (mesh%V(vr,2)-mesh%V(vi,2)) * (mesh%V(vj,1)-mesh%V(vi,1)) )
                  
            Nx = -((Nxl / Nzl) + (Nxr / Nzr))/2._dp
            Ny = -((Nyl / Nzl) + (Nyr / Nzr))/2._dp
            
          ELSE
            ! Surface slope on an edge segment Ac vertex equals that of the (single) triangle
              
            ! Find index of vl
            DO iti = 1, mesh%niTri(vi)
              ti = mesh%iTri(vi,iti)
              DO n1 = 1, 3
                n2 = n1+1
                IF (n2==4) n2=1
                n3 = n2+1
                IF (n3==4) n3=1
                
                IF ((mesh%Tri(ti,n1)==vi .AND. mesh%Tri(ti,n2)==vj) .OR. (mesh%Tri(ti,n1)==vj .AND. mesh%Tri(ti,n2)==vi)) THEN
                  vl = mesh%Tri(ti,n3)
                END IF
              END DO
            END DO
            
            ! List relevant Aa vertex indices vor Ac vertex in reference array
            ! NOTE: fourth entry is set to 1 because it cant be zero, but the fourth
            !       neighbour function is set to zero, so it doesn't matter...
            Aci_temp(mesh%nAc,:) = [vi, vj, vl, 1]
            
            ! Calculate neighbour functions on single adjacent triangle (equal to Ac vertex values)
            ! x an y just as average of the left and right triangles
            Nxl = [mesh%V(vl,2)-mesh%V(vj,2), mesh%V(vi,2)-mesh%V(vl,2), mesh%V(vj,2)-mesh%V(vi,2), 0._dp]
            Nyl = [mesh%V(vj,1)-mesh%V(vl,1), mesh%V(vl,1)-mesh%V(vi,1), mesh%V(vi,1)-mesh%V(vj,1), 0._dp]
            Nzl = ( (mesh%V(vj,1)-mesh%V(vi,1)) * (mesh%V(vl,2)-mesh%V(vi,2)) ) - &
                  ( (mesh%V(vj,2)-mesh%V(vi,2)) * (mesh%V(vl,1)-mesh%V(vi,1)) )
                  
            Nx = -Nxl / Nzl
            Ny = -Nyl / Nzl           
          
          END IF
            
          ! Then p and o by rotating x and y
          Ux = mesh%V( vj,1) - mesh%V( vi,1)
          Uy = mesh%V( vj,2) - mesh%V( vi,2)
          U  = SQRT(Ux**2 + Uy**2)
          
          Np = 1._dp / U          ! Slope along line from j to i is just f(j)-f(i)/dist(ji), only need to save that single number
                                  ! Checked that this gives the same result as just rotating Nx, Ny
          No = (Ny*Ux - Nx*Uy)/U  ! This one still depends on all four vertices, no symmetries...
          
          ! Save everything in the temporary memory
          Nx_Ac_temp( mesh%nAc,:) = Nx
          Ny_Ac_temp( mesh%nAc,:) = Ny
          Np_Ac_temp( mesh%nAc  ) = Np
          No_Ac_temp( mesh%nAc,:) = No 
            
        END DO ! DO ci = 1, mesh%nC(vi)
      END DO ! DO vi = 1, mesh%nV
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate shared memory, move data there
    ! =======================================
    
    CALL allocate_shared_dp_2D(  mesh%nAc,       2,            mesh%VAc,             mesh%wVAc           )
    CALL allocate_shared_int_2D( mesh%nAc,       4,            mesh%Aci,             mesh%wAci           )
    CALL allocate_shared_dp_2D(  mesh%nAc,       4,            mesh%Nx_Ac,           mesh%wNx_Ac         )
    CALL allocate_shared_dp_2D(  mesh%nAc,       4,            mesh%Ny_Ac,           mesh%wNy_Ac         )
    CALL allocate_shared_dp_1D(  mesh%nAc,                     mesh%Np_Ac,           mesh%wNp_Ac         )
    CALL allocate_shared_dp_2D(  mesh%nAc,       4,            mesh%No_Ac,           mesh%wNo_Ac         )
    
    IF (par%master) THEN
      mesh%VAc   = VAc_temp(   1:mesh%nAc,:)
      mesh%Aci   = Aci_temp(   1:mesh%nAc,:)
      mesh%Nx_Ac = Nx_Ac_temp( 1:mesh%nAc,:)
      mesh%Ny_Ac = Ny_Ac_temp( 1:mesh%nAc,:)
      mesh%Np_Ac = Np_Ac_temp( 1:mesh%nAc  )
      mesh%No_Ac = No_Ac_temp( 1:mesh%nAc,:)
    END IF
    CALl sync
    
    ! Deallocate temporary memory
    ! ===========================
    
    DEALLOCATE( VAc_temp)
    DEALLOCATE( Aci_temp)
    DEALLOCATE( Nx_Ac_temp)
    DEALLOCATE( Ny_Ac_temp)
    DEALLOCATE( Np_Ac_temp)
    DEALLOCATE( No_Ac_temp)
    
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
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: vi, vj, ci, aci, ai, vk, vl, vr, aci1, aci2, aci3, aci4
    LOGICAL                                       :: switchthem
    INTEGER                                       :: ti, vip, viq, vir, aci_pq, aci_qr, aci_rp, acj, ack, iti
    REAL(dp), DIMENSION(2)                        :: V_vi
    INTEGER                                       :: n, aj
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: V_vc
    LOGICAL                                       :: is_edge
    
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
    
    CALL allocate_shared_dp_2D(  mesh%nVAaAc,   mesh%nC_mem+1, mesh%Nx_AaAc,  mesh%wNx_AaAc )
    CALL allocate_shared_dp_2D(  mesh%nVAaAc,   mesh%nC_mem+1, mesh%Ny_AaAc,  mesh%wNy_AaAc )
    CALL allocate_shared_dp_2D(  mesh%nVAaAc,   mesh%nC_mem+1, mesh%Nxx_AaAc, mesh%wNxx_AaAc)
    CALL allocate_shared_dp_2D(  mesh%nVAaAc,   mesh%nC_mem+1, mesh%Nxy_AaAc, mesh%wNxy_AaAc)
    CALL allocate_shared_dp_2D(  mesh%nVAaAc,   mesh%nC_mem+1, mesh%Nyy_AaAc, mesh%wNyy_AaAc)
    
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

    ! Triangles (not sure if these are needed but they're nice for plotting)
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

    ! Neighbour functions
    
    ALLOCATE( V_vc( mesh%nC_mem, 2))
    
    DO ai = mesh%a1, mesh%a2
    
      vi  = ai
      aci = ai - mesh%nV
      
      V_vi = mesh%VAaAc(  ai,:)
      n    = mesh%nCAaAc( ai  )
      DO ci = 1, n
        aj = mesh%CAaAc( ai,ci)
        V_vc( ci,:) = mesh%VAaAc( aj,:)
      END DO
      
      is_edge = .FALSE.
      IF (ai <= mesh%nV) THEN
        IF (mesh%edge_index(    vi ) > 0) is_edge = .TRUE.
      ELSE
        IF (mesh%edge_index_Ac( aci) > 0) is_edge = .TRUE.
      END IF
      
      CALL get_neighbour_functions_vertex_gr( V_vi, n, V_vc, is_edge, &
        mesh%Nx_AaAc(  ai,:), &
        mesh%Ny_AaAc(  ai,:), &
        mesh%Nxx_AaAc( ai,:), &
        mesh%Nxy_AaAc( ai,:), &
        mesh%Nyy_AaAc( ai,:))
      
    END DO ! DO ai = mesh%a1, mesh%%a2
    CALL sync
    
  END SUBROUTINE make_combined_AaAc_mesh
  
  SUBROUTINE get_mesh_derivatives_Ac( mesh, d_Aa, ddx_Ac, ddy_Ac, ddp_Ac, ddo_Ac)    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d_Aa
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: ddx_Ac, ddy_Ac, ddp_Ac, ddo_Ac

    INTEGER                                       :: aci

    DO aci = mesh%ac1, mesh%ac2
      CALL get_mesh_derivatives_vertex_Ac( mesh, d_Aa, ddx_Ac( aci), ddy_Ac( aci), ddp_Ac( aci), ddo_Ac( aci), aci)
    END DO ! DO ci = mesh%ac1, mesh%ac2
    CALL sync
    
  END SUBROUTINE get_mesh_derivatives_Ac
  SUBROUTINE get_mesh_derivatives_vertex_Ac( mesh, d_Aa, ddx_Ac, ddy_Ac, ddp_Ac, ddo_Ac, aci)
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d_Aa
    REAL(dp),                   INTENT(OUT)       :: ddx_Ac, ddy_Ac, ddp_Ac, ddo_Ac
    INTEGER,                    INTENT(IN)        :: aci

    INTEGER                                       :: n

    ddx_Ac  = 0._dp
    ddy_Ac  = 0._dp
    ddp_Ac  = 0._dp
    ddo_Ac  = 0._dp

    DO n = 1, 4
      ddx_Ac = ddx_Ac + mesh%Nx_Ac( aci,n) * d_Aa( mesh%Aci( aci,n))
      ddy_Ac = ddy_Ac + mesh%Ny_Ac( aci,n) * d_Aa( mesh%Aci( aci,n))
      ddo_Ac = ddo_Ac + mesh%No_Ac( aci,n) * d_Aa( mesh%Aci( aci,n))
    END DO
    ddp_Ac = mesh%Np_Ac( aci) * (d_Aa( mesh%Aci( aci,2)) - d_Aa( mesh%Aci( aci,1)))
    
  END SUBROUTINE get_mesh_derivatives_vertex_Ac
  
  SUBROUTINE get_mesh_derivatives_AaAc( mesh, d_AaAc, ddx_AaAc, ddy_AaAc)    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d_AaAc
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: ddx_AaAc, ddy_AaAc

    INTEGER                                       :: ai

    DO ai = mesh%a1, mesh%a2
      CALL get_mesh_derivatives_vertex_AaAc( mesh, d_AaAc, ddx_AaAc( ai), ddy_AaAc( ai), ai)
    END DO ! DO ai = mesh%a1, mesh%a2
    CALL sync
    
  END SUBROUTINE get_mesh_derivatives_AaAc
  SUBROUTINE get_mesh_derivatives_vertex_AaAc( mesh, d_AaAc, ddx_AaAc, ddy_AaAc, ai)
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d_AaAc
    REAL(dp),                   INTENT(OUT)       :: ddx_AaAc, ddy_AaAc
    INTEGER,                    INTENT(IN)        :: ai

    INTEGER                                       :: ci
    
    ddx_AaAc    = mesh%Nx_AaAc( ai,mesh%nCAaAc( ai)+1) * d_AaAc( ai)
    ddy_AaAc    = mesh%Ny_AaAc( ai,mesh%nCAaAc( ai)+1) * d_AaAc( ai)
    
    DO ci = 1, mesh%nCAaAc( ai)
      ddx_AaAc  = ddx_AaAc  + mesh%Nx_AaAc( ai,ci) * d_AaAc( mesh%CAaAc( ai,ci))
      ddy_AaAc  = ddy_AaAc  + mesh%Ny_AaAc( ai,ci) * d_AaAc( mesh%CAaAc( ai,ci))
    END DO
    
  END SUBROUTINE get_mesh_derivatives_vertex_AaAc
  
  SUBROUTINE get_mesh_curvatures_AaAc( mesh, d_AaAc, ddxx_AaAc, ddxy_AaAc, ddyy_AaAc)    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d_AaAc
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: ddxx_AaAc, ddxy_AaAc, ddyy_AaAc

    INTEGER                                       :: ai

    DO ai = mesh%a1, mesh%a2
      CALL get_mesh_curvatures_vertex_AaAc( mesh, d_AaAc, ddxx_AaAc( ai), ddxy_AaAc( ai), ddyy_AaAc( ai), ai)
    END DO ! DO ai = mesh%a1, mesh%a2
    CALL sync
    
  END SUBROUTINE get_mesh_curvatures_AaAc
  SUBROUTINE get_mesh_curvatures_vertex_AaAc( mesh, d_AaAc, ddxx_AaAc, ddxy_AaAc, ddyy_AaAc, ai)
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d_AaAc
    REAL(dp),                   INTENT(OUT)       :: ddxx_AaAc, ddxy_AaAc, ddyy_AaAc
    INTEGER,                    INTENT(IN)        :: ai

    INTEGER                                       :: ci, ac

    ddxx_AaAc = d_AaAc( ai) * mesh%Nxx_AaAc( ai, mesh%nCAaAc( ai)+1)
    ddxy_AaAc = d_AaAc( ai) * mesh%Nxy_AaAc( ai, mesh%nCAaAc( ai)+1)
    ddyy_AaAc = d_AaAc( ai) * mesh%Nyy_AaAc( ai, mesh%nCAaAc( ai)+1)
    
    DO ci = 1, mesh%nCAaAc( ai)
      ac = mesh%CAaAc( ai,ci)
      ddxx_AaAc = ddxx_AaAc + d_AaAc( ai) * mesh%Nxx_AaAc( ai,ci)
      ddxy_AaAc = ddxy_AaAc + d_AaAc( ai) * mesh%Nxy_AaAc( ai,ci)
      ddyy_AaAc = ddyy_AaAc + d_AaAc( ai) * mesh%Nyy_AaAc( ai,ci)
    END DO
    
  END SUBROUTINE get_mesh_curvatures_vertex_AaAc
  
  SUBROUTINE apply_Neumann_boundary_AaAc( mesh, d_AaAc)
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(INOUT)     :: d_AaAc
    
    INTEGER                                       :: ai, ci, ac
    REAL(dp), DIMENSION(mesh%nC_mem)              :: vals
    INTEGER                                       :: nvals
        
    ! Apply Neumann boundary conditions: make sure ddx and ddy are zero at domain boundary
    ! Achieved by setting d for edge vertices to average value of non-edge neighbours
    ! ==========================================================================================
            
    DO ai = MAX(5,mesh%a1), mesh%a2   ! Edge vertices but not the four corners
    
      IF (ai <= mesh%nV) THEN
        IF (mesh%edge_index(    ai          ) == 0) CYCLE
      ELSE
        IF (mesh%edge_index_Ac( ai - mesh%nV) == 0) CYCLE
      END IF

      vals = 0._dp
      nvals = 0
  
      DO ci = 1, mesh%nCAaAc( ai)
      
        ac = mesh%CAaAc( ai,ci)
        
        IF (ac <= mesh%nV) THEN
          IF (mesh%edge_index(    ac          ) > 0) CYCLE
        ELSE
          IF (mesh%edge_index_Ac( ac - mesh%nV) > 0) CYCLE
        END IF
        
        nvals = nvals+1
        vals( nvals) = d_AaAc( ac)
      END DO
  
      d_AaAc( ai) = SUM(vals) / nvals
      
    END DO ! DO ai = MAX(5,mesh%a1), mesh%a2   ! Edge vertices but not the four corners
    CALL sync

    ! And the four corner vertices separately
    IF (par%master) THEN
      DO ai = 1, 4
  
        vals = 0._dp
        nvals = 0
  
        DO ci = 1, mesh%nCAaAc( ai)
          ac = mesh%CAaAc( ai,ci)
          nvals = nvals+1
          vals( nvals) = d_AaAc( ac)
        END DO
  
        d_AaAc( ai) = SUM(vals) / nvals
           
      END DO ! DO ai = 1, 4
    END IF
    CALL sync
    
  END SUBROUTINE apply_Neumann_boundary_AaAc 
  
  SUBROUTINE map_Aa_to_Ac(    mesh, d_Aa, d_Ac)
    ! Map data from the regular Aa mesh to the staggered Ac mesh    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d_Aa
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: d_Ac

    INTEGER                                       :: aci, vi, vj
    
    DO aci = mesh%ac1, mesh%ac2
    
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      
      d_Ac( aci) = (d_Aa( vi) + d_Aa( vj)) / 2._dp
    
    END DO ! DO ci = mesh%ac1, mesh%ac2
    CALL sync
    
  END SUBROUTINE map_Aa_to_Ac  
  SUBROUTINE map_Aa_to_Ac_3D( mesh, d_Aa, d_Ac)
    ! Map data from the regular Aa mesh to the staggered Ac mesh    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: d_Aa
    REAL(dp), DIMENSION(:,:),   INTENT(OUT)       :: d_Ac

    INTEGER                                       :: aci, vi, vj
    
    DO aci = mesh%ac1, mesh%ac2
    
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      
      d_Ac( aci,:) = (d_Aa( vi,:) + d_Aa( vj,:)) / 2._dp
    
    END DO ! DO ci = mesh%ac1, mesh%ac2
    CALL sync
    
  END SUBROUTINE map_Aa_to_Ac_3D  
  SUBROUTINE map_Ac_to_Aa(    mesh, d_Ac, d_Aa)
    ! Map data from the regular Aa mesh to the staggered Ac mesh    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d_Ac
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: d_Aa

    INTEGER                                       :: vi, ci, aci
    
    d_Aa( mesh%v1:mesh%v2) = 0._dp
    
    DO vi = mesh%v1, mesh%v2
      DO ci = 1, mesh%nC( vi)
        aci = mesh%iAci( vi,ci)
        d_Aa( vi) = d_Aa( vi) + d_Ac( aci) / mesh%nC( vi)
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
  END SUBROUTINE map_Ac_to_Aa  
  SUBROUTINE map_Ac_to_Aa_3D( mesh, d_Ac, d_Aa)
    ! Map data from the regular Aa mesh to the staggered Ac mesh    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: d_Ac
    REAL(dp), DIMENSION(:,:),   INTENT(OUT)       :: d_Aa

    INTEGER                                       :: vi, ci, aci
    
    d_Aa( mesh%v1:mesh%v2,:) = 0._dp
    
    DO vi = mesh%v1, mesh%v2
      DO ci = 1, mesh%nC( vi)
        aci = mesh%iAci( vi,ci)
        d_Aa( vi,:) = d_Aa(vi,:) + d_Ac( aci,:) / mesh%nC( vi)
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
  END SUBROUTINE map_Ac_to_Aa_3D
  
  SUBROUTINE rotate_xy_to_po( mesh, d_Ac_x, d_Ac_y, d_Ac_p, d_Ac_o)
    ! Rotate a vector field on the Ac mesh from [x,y] to [p,o] components
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d_Ac_x, d_Ac_y
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: d_Ac_p, d_Ac_o
    
    ! Local variables:
    INTEGER                                       :: ci, vi, vj
    REAL(dp)                                      :: Dx, Dy, D
    
    DO ci = mesh%ac1, mesh%ac2
    
      vi = mesh%Aci( ci,1)
      vj = mesh%Aci( ci,2)
      
      Dx = mesh%V( vj,1) - mesh%V( vi,1)
      Dy = mesh%V( vj,2) - mesh%V( vi,2)
      D  = SQRT(Dx**2 + Dy**2)
      
      d_Ac_p( ci) = d_Ac_x( ci) * Dx/D + d_Ac_y( ci) * Dy/D
      d_Ac_o( ci) = d_Ac_y( ci) * Dx/D - d_Ac_x( ci) * Dy/D
      
    END DO ! DO ci = mesh%ac1, mesh%ac2
    CALL sync
    
  END SUBROUTINE rotate_xy_to_po

END MODULE mesh_ArakawaC_module