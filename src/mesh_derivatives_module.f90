MODULE mesh_derivatives_module
  ! Routines for calculating the neighbour functions on a mesh.

  USE mpi
  USE configuration_module,          ONLY: dp, C
  USE parallel_module,               ONLY: par, sync, ierr, cerr, write_to_memory_log
  USE data_types_module,             ONLY: type_mesh
  USE mesh_help_functions_module,    ONLY: is_in_triangle

  IMPLICIT NONE

  CONTAINS
  
  SUBROUTINE get_neighbour_functions( mesh)
    ! Calculate the first and second order neighbour functions    
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    INTEGER                                       :: ti, vi, ci, vc
    REAL(dp)                                      :: ax, bx, cx, ay, by, cy, D
    REAL(dp), DIMENSION(2)                        :: V_vi
    INTEGER                                       :: n
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: V_vc
    LOGICAL                                       :: is_edge

    ! == First order on the triangles (used for some interpolation routines)
    mesh%NxTri( mesh%t1:mesh%t2,:) = 0._dp
    mesh%NyTri( mesh%t1:mesh%t2,:) = 0._dp
    CALL sync

    DO ti = mesh%t1, mesh%t2
      ax = mesh%V(mesh%Tri(ti,1),1)
      ay = mesh%V(mesh%Tri(ti,1),2)
      bx = mesh%V(mesh%Tri(ti,2),1)
      by = mesh%V(mesh%Tri(ti,2),2)
      cx = mesh%V(mesh%Tri(ti,3),1)
      cy = mesh%V(mesh%Tri(ti,3),2)

      D = ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)

      mesh%NxTri(ti,:) = [(by - cy), (cy - ay), (ay - by)] / D
      mesh%NyTri(ti,:) = [(cx - bx), (ax - cx), (bx - ax)] / D
    END DO ! DO ti = mesh%t1, mesh%t2
    CALL sync
    
    ! == First and second order on the vertices
    ALLOCATE( V_vc( mesh%nC_mem, 2))
    
    DO vi = mesh%v1, mesh%v2
    
      ! The vertex itself
      V_vi = mesh%V( vi,:)
      
      ! The neighbours
      n = mesh%nC( vi)
      DO ci = 1, n
        vc = mesh%C( vi,ci)
        V_vc( ci,:) = mesh%V( vc,:)
      END DO
      
      is_edge = .FALSE.
      IF (mesh%edge_index( vi) > 0) is_edge = .TRUE.
    
      CALL get_neighbour_functions_vertex_gr( V_vi, n, V_vc, is_edge, mesh%Nx( vi,:), mesh%Ny( vi,:), mesh%Nxx( vi,:), mesh%Nxy( vi,:), mesh%Nyy( vi,:))
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    DEALLOCATE( V_vc)

  END SUBROUTINE get_neighbour_functions
  SUBROUTINE get_neighbour_functions_vertex_gr( V_vi, n, V_vc, is_edge, Nx, Ny, Nxx, Nxy, Nyy)
    ! Calculate the first and second order neighbour functions on the vertices using the average gradient approach.
    ! All Equation numbers refer to the appendix of the original UFEMISM publication
    ! (Berends et al., GMD, 2021)
    
    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: V_vi
    INTEGER,                    INTENT(IN)        :: n
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: V_vc
    LOGICAL,                    INTENT(IN)        :: is_edge
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: Nx, Ny, Nxx, Nxy, Nyy
    
    ! Local variables:
    INTEGER                                       :: nTri, nSub, ti, tip1s, si, sip1s, sip2s, ci, cim1s, sim1s, sim2s
    REAL(dp), DIMENSION( n,3)                     :: NxTri, NyTri, NxSub, NySub
    REAL(dp)                                      :: xi, yi, xt, yt, xtp1s, ytp1s, nzt
    REAL(dp)                                      :: xs, ys, xsp1s, ysp1s, xsp2s, ysp2s, nzs
    REAL(dp)                                      :: Axxis, Bxxis, Cxxis
    REAL(dp)                                      :: Axyis, Bxyis, Cxyis
    REAL(dp)                                      :: Ayyis, Byyis, Cyyis
    
    Nx  = 0._dp
    Ny  = 0._dp
    Nxx = 0._dp
    Nxy = 0._dp
    Nyy = 0._dp

    IF (is_edge) THEN
      nTri = n-1
      nSub = n-2
    ELSE
      nTri = n
      nSub = n
    END IF

  ! Calculate neighbour functions on the triangles and subtriangles
  ! ===============================================================

    DO ti = 1, nTri

      ! Star notation
      tip1s = ti+1
      IF (tip1s > n) tip1s = tip1s - n

      ! The coordinates of the three vertices spanning this triangle
      xi    = V_vi( 1)
      yi    = V_vi( 2)
      xt    = V_vc( ti,1)
      yt    = V_vc( ti,2)
      xtp1s = V_vc( tip1s,1)
      ytp1s = V_vc( tip1s,2)

      ! The z-component of the normal vector
      nzt = ((xt - xi) * (ytp1s - yi)) - ((yt - yi) * (xtp1s - xi))  ! Eq. A1

      ! The neighbour functions
      NxTri( ti,:) = [(yt - ytp1s), (ytp1s - yi), (yi - yt)] / nzt   ! Eq. A3a
      NyTri( ti,:) = [(xtp1s - xt), (xi - xtp1s), (xt - xi)] / nzt   ! Eq. A3b

    END DO ! DO ti = 1, nTri

    DO si = 1, nSub

      ! Star notation
      sip1s = si+1
      IF (sip1s > n) sip1s = sip1s - n
      sip2s = sip1s+1
      IF (sip2s > n) sip2s = sip2s - n

      ! The coordinates of the four vertices involved in this sub-triangle
      xi    = V_vi( 1)
      yi    = V_vi( 2)
      xs    = V_vc( si,1)
      ys    = V_vc( si,2)
      xsp1s = V_vc( sip1s,1)
      ysp1s = V_vc( sip1s,2)
      xsp2s = V_vc( sip2s,1)
      ysp2s = V_vc( sip2s,2)

      ! The z-component of the normal vector
      nzs = (1._dp/3._dp * (xs + xsp1s - 2*xi) * (ysp1s + ysp2s - 2*yi)) - &                   ! Eq. A11
            (1._dp/3._dp * (ys + ysp1s - 2*yi) * (xsp1s + xsp2s - 2*xi))

      ! The neighbour functions
      NxSub( si,:) = [(ys - ysp2s), (ysp1s + ysp2s - 2._dp*yi), (2._dp*yi - ys - ysp1s)] / nzs ! Eq. A13a
      NySub( si,:) = [(xsp2s - xs), (2._dp*xi - xsp1s - xsp2s), (xs + xsp1s - 2._dp*xi)] / nzs ! Eq. A13b

    END DO ! DO si = 1, nSub

  ! Calculate 1st order neighbour functions on the vertex
  ! =====================================================

    IF (.NOT. is_edge) THEN
      ! Free vertex

      ! Home vertex
      Nx( n+1) = 1._dp/REAL(n,dp) * SUM( NxTri( 1:n,1))   ! Eq. A7a
      Ny( n+1) = 1._dp/REAL(n,dp) * SUM( NyTri( 1:n,1))

      ! Neighbours
      DO ci = 1, n

        ! Star notation
        cim1s = ci-1
        IF (cim1s == 0) cim1s = cim1s + n

        Nx( ci) = 1._dp/n * (NxTri( ci,2) + NxTri( cim1s,3))     ! Eq. A7b
        Ny( ci) = 1._dp/n * (NyTri( ci,2) + NyTri( cim1s,3))
      END DO ! DO ci = 1, n

    ELSE ! IF (.NOT. is_edge) THEN
      ! Boundary vertex

      ! Home vertex
      Nx( n+1) = 1._dp/REAL(n-1,dp) * SUM( NxTri( 1:n-1,1))   ! Eq. A9a
      Ny( n+1) = 1._dp/REAL(n-1,dp) * SUM( NyTri( 1:n-1,1))

      ! Neighbours
      DO ci = 1, n
        IF     (ci == 1) THEN
          Nx( ci) = 1._dp/REAL(n-1,dp) * NxTri( 1,2)                      ! Eq. A9b
          Ny( ci) = 1._dp/REAL(n-1,dp) * NyTri( 1,2)                      ! Eq. A9b
        ELSEIF (ci == n) THEN
          Nx( ci) = 1._dp/REAL(n-1,dp) * NxTri( n-1,3)                    ! Eq. A9b
          Ny( ci) = 1._dp/REAL(n-1,dp) * NyTri( n-1,3)                    ! Eq. A9b
        ELSE
          Nx( ci) = 1._dp/REAL(n-1,dp) * (NxTri( ci,2) + NxTri( ci-1,3))  ! Eq. A9b
          Ny( ci) = 1._dp/REAL(n-1,dp) * (NyTri( ci,2) + NyTri( ci-1,3))  ! Eq. A9b
        END IF
      END DO ! DO ci = 1, n

    END IF ! IF (.NOT. is_edge) THEN

  ! Calculate 2nd order neighbour functions on the vertex
  ! =====================================================

    IF (.NOT. is_edge) THEN
      ! Free vertex

      ! Home vertex
      DO si = 1, n

        ! Star notation
        sip1s = si + 1
        IF (sip1s > n) sip1s = sip1s - n

        Nxx( n+1) = Nxx( n+1) + (1._dp/REAL(n,dp) * ((NxSub( si,1) * Nx( n+1)) + &         ! Eq. A23a
                                                     (NxSub( si,2) * NxTri( si,   1)) + &
                                                     (NxSub( si,3) * NxTri( sip1s,1))))
        Nxy( n+1) = Nxy( n+1) + (1._dp/REAL(n,dp) * ((NySub( si,1) * Nx( n+1)) + &         ! Eq. A23b
                                                     (NySub( si,2) * NxTri( si,   1)) + &
                                                     (NySub( si,3) * NxTri( sip1s,1))))
        Nyy( n+1) = Nyy( n+1) + (1._dp/REAL(n,dp) * ((NySub( si,1) * Ny( n+1)) + &         ! Eq. A23c
                                                     (NySub( si,2) * NyTri( si,   1)) + &
                                                     (NySub( si,3) * NyTri( sip1s,1))))
      END DO ! DO si = 1, n

      ! Neighbours
      DO si = 1, n

        ! Star notation
        sim1s = si - 1
        IF (sim1s < 1) sim1s = sim1s + n
        sim2s = sim1s - 1
        IF (sim2s < 1) sim2s = sim2s + n

        Nxx( si) = 1._dp/REAL(n,dp) * ((NxTri( si,   2) * (NxSub( si,   2) + NxSub( sim1s,3)) + &  ! Eq. A23d
                                       (NxTri( sim1s,3) * (NxSub( sim1s,2) + NxSub( sim2s,3)) + &
                                       (Nx( si) * SUM( NxSub( 1:n,1))))))
        Nxy( si) = 1._dp/REAL(n,dp) * ((NxTri( si,   2) * (NySub( si,   2) + NySub( sim1s,3)) + &  ! Eq. A23e
                                       (NxTri( sim1s,3) * (NySub( sim1s,2) + NySub( sim2s,3)) + &
                                       (Nx( si) * SUM( NySub( 1:n,1))))))
        Nyy( si) = 1._dp/REAL(n,dp) * ((NyTri( si,   2) * (NySub( si,   2) + NySub( sim1s,3)) + &  ! Eq. A23f
                                       (NyTri( sim1s,3) * (NySub( sim1s,2) + NySub( sim2s,3)) + &
                                       (Ny( si) * SUM( NySub( 1:n,1))))))

      END DO ! DO si = 1, n

    ELSE ! IF (~is_edge)
      ! Boundary vertex

      ! Home vertex
      DO si = 1, n-2
        Nxx( n+1) = Nxx( n+1) + (1._dp/REAL(n-2,dp) * ((NxSub( si,1) * Nx( n+1)) + &         ! Eq. A24a
                                                       (NxSub( si,2) * NxTri( si,  1)) + &
                                                       (NxSub( si,3) * NxTri( si+1,1))))
        Nxy( n+1) = Nxy( n+1) + (1._dp/REAL(n-2,dp) * ((NySub( si,1) * Nx( n+1)) + &         ! Eq. A24b
                                                       (NySub( si,2) * NxTri( si,  1)) + &
                                                       (NySub( si,3) * NxTri( si+1,1))))
        Nyy( n+1) = Nyy( n+1) + (1._dp/REAL(n-2,dp) * ((NySub( si,1) * Ny( n+1)) + &         ! Eq. A24c
                                                       (NySub( si,2) * NyTri( si,  1)) + &
                                                       (NySub( si,3) * NyTri( si+1,1))))
      END DO ! DO si = 1, n-2

      ! Neighbours
      DO si = 1, n

        IF (si < n-1) THEN
          Axxis = NxSub( si,2) * NxTri( si,2)  ! Eq. A25a
          Axyis = NySub( si,2) * NxTri( si,2)  ! Eq. A25b
          Ayyis = NySub( si,2) * NyTri( si,2)  ! Eq. A25c
        ELSE
          Axxis = 0._dp
          Axyis = 0._dp
          Ayyis = 0._dp
        END IF

        IF (si > 1 .AND. si < n) THEN
          Bxxis = (NxSub( si-1,2) * NxTri( si-1,3)) + (NxSub( si-1,3) * NxTri( si,2))  ! Eq. A25d
          Bxyis = (NySub( si-1,2) * NxTri( si-1,3)) + (NySub( si-1,3) * NxTri( si,2))  ! Eq. A25e
          Byyis = (NySub( si-1,2) * NyTri( si-1,3)) + (NySub( si-1,3) * NyTri( si,2))  ! Eq. A25f
        ELSE
          Bxxis = 0._dp
          Bxyis = 0._dp
          Byyis = 0._dp
        END IF

        IF (si > 2) THEN
          Cxxis = NxSub( si-2,3) * NxTri( si-1,3)  ! Eq. A25g
          Cxyis = NySub( si-2,3) * NxTri( si-1,3)  ! Eq. A25h
          Cyyis = NySub( si-2,3) * NyTri( si-1,3)  ! Eq. A25i
        ELSE
          Cxxis = 0._dp
          Cxyis = 0._dp
          Cyyis = 0._dp
        END IF

        Nxx( si) = 1._dp/REAL(n-2,dp) * (Axxis + Bxxis + Cxxis + (Nx( si) * SUM( NxSub( 1:n-2,1)))) ! Eq. A24d
        Nxy( si) = 1._dp/REAL(n-2,dp) * (Axyis + Bxyis + Cxyis + (Nx( si) * SUM( NySub( 1:n-2,1)))) ! Eq. A24e
        Nyy( si) = 1._dp/REAL(n-2,dp) * (Ayyis + Byyis + Cyyis + (Ny( si) * SUM( NySub( 1:n-2,1)))) ! Eq. A24f

      END DO ! DO si = 1, n

    END IF ! IF (.NOT. is_edge) THEN
      
  END SUBROUTINE get_neighbour_functions_vertex_gr
  
  SUBROUTINE get_mesh_derivatives(                  mesh, d, ddx, ddy)
    ! Get first order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: ddx, ddy

    INTEGER                                       :: vi

    DO vi = mesh%v1, mesh%v2
      CALL get_mesh_derivatives_vertex( mesh, d, ddx( vi), ddy( vi), vi)
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE get_mesh_derivatives  
  SUBROUTINE get_mesh_derivatives_3D(               mesh, d, ddx, ddy, nz)
    ! Get first order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: d
    REAL(dp), DIMENSION(:,:),   INTENT(OUT)       :: ddx, ddy
    INTEGER,                    INTENT(IN)        :: nz

    INTEGER                                       :: vi, k

    DO vi = mesh%v1, mesh%v2
      DO k = 1, nZ
        CALL get_mesh_derivatives_vertex_3D( mesh, d, ddx( vi,k), ddy( vi,k), vi, k)
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE get_mesh_derivatives_3D
  SUBROUTINE get_mesh_derivatives_vertex(           mesh, d, ddx, ddy, vi)
    ! Get first order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d
    REAL(dp),                   INTENT(OUT)       :: ddx, ddy
    INTEGER                                       :: vi

    INTEGER                                       :: ci

    ddx    = mesh%Nx( vi,mesh%nC( vi)+1) * d( vi)
    ddy    = mesh%Ny( vi,mesh%nC( vi)+1) * d( vi)
    DO ci = 1, mesh%nC( vi)
      ddx  = ddx  + mesh%Nx( vi,ci) * d( mesh%C( vi,ci))
      ddy  = ddy  + mesh%Ny( vi,ci) * d( mesh%C( vi,ci))
    END DO

  END SUBROUTINE get_mesh_derivatives_vertex  
  SUBROUTINE get_mesh_derivatives_vertex_3D(        mesh, d, ddx, ddy, vi, k)
    ! Get first order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: d
    REAL(dp),                   INTENT(OUT)       :: ddx, ddy
    INTEGER                                       :: vi, k

    INTEGER                                       :: ci

    ddx    = mesh%Nx( vi,mesh%nC( vi)+1) * d( vi,k)
    ddy    = mesh%Ny( vi,mesh%nC( vi)+1) * d( vi,k)

    DO ci = 1, mesh%nC( vi)
      ddx  = ddx  + mesh%Nx( vi,ci) * d( mesh%C( vi,ci),k)
      ddy  = ddy  + mesh%Ny( vi,ci) * d( mesh%C( vi,ci),k)
    END DO

  END SUBROUTINE get_mesh_derivatives_vertex_3D
  SUBROUTINE get_mesh_curvatures(                   mesh, d, ddxx, ddxy, ddyy)
    ! Get second order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: ddxx, ddxy,ddyy

    INTEGER                                       :: vi

    DO vi = mesh%v1, mesh%v2
      CALL get_mesh_curvatures_vertex( mesh, d, ddxx( vi), ddxy( vi), ddyy( vi), vi)
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE get_mesh_curvatures 
  SUBROUTINE get_mesh_curvatures_vertex(            mesh, d, ddxx, ddxy, ddyy, vi)
    ! Get second order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d
    INTEGER,                    INTENT(IN)        :: vi
    REAL(dp),                   INTENT(OUT)       :: ddxx, ddxy,ddyy

    INTEGER                                       :: ci

    ddxx = mesh%Nxx( vi,mesh%nC( vi)+1) * d( vi)
    ddxy = mesh%Nxy( vi,mesh%nC( vi)+1) * d( vi)
    ddyy = mesh%Nyy( vi,mesh%nC( vi)+1) * d( vi)
      
    DO ci = 1, mesh%nC( vi)
      ddxx = ddxx + mesh%Nxx( vi,ci) * d( mesh%C( vi,ci))
      ddxy = ddxy + mesh%Nxy( vi,ci) * d( mesh%C( vi,ci))
      ddyy = ddyy + mesh%Nyy( vi,ci) * d( mesh%C( vi,ci))
    END DO

  END SUBROUTINE get_mesh_curvatures_vertex
  
  SUBROUTINE get_upwind_derivative_vertex_3D( mesh, U, V, d, vi, k, ddx, ddy)
    ! Get the "upwind" derivatives of d along U,V (used for thermodynamics)
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: U, V
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: d
    INTEGER,                    INTENT(IN)        :: vi, k
    REAL(dp),                   INTENT(OUT)       :: ddx, ddy
    
    INTEGER                                       :: iti, ti, tup
    REAL(dp), DIMENSION(2)                        :: W, p, pa, pb, pc
      
    ddx = 0._dp
    ddy = 0._dp
    
    ! If there's no "wind", there's no upwind direction.
    IF (ABS(U(vi,k)) < 1E-10_dp .AND. ABS(V(vi,k)) < 1E-10_dp) THEN
      CALL get_mesh_derivatives_vertex_3D( mesh, d, ddx, ddy, vi, k)
      RETURN
    END IF
    
    ! Find the upwind triangle
    W = [U(vi,k), V(vi,k)] * mesh%R(vi) / (4._dp * SQRT(U(vi,k)**2 + V(vi,k)**2))
    p = mesh%V(vi,:) - W
        
    tup = 0
    DO iti = 1, mesh%niTri(vi)
      ti = mesh%iTri(vi,iti)
      pa = mesh%V(mesh%Tri(ti,1),:)
      pb = mesh%V(mesh%Tri(ti,2),:)
      pc = mesh%V(mesh%Tri(ti,3),:)
      IF (is_in_triangle(pa, pb, pc, p)) THEN
        tup = ti
        EXIT
      END IF
    END DO
    
    IF (tup==0) THEN
      WRITE(0,*) 'ERROR: couldnt find upwind triangle for vi = ', vi, ' at [', mesh%V(vi,1), ',', mesh%V(vi,2), ']'
      WRITE(0,*) '       [U,V] = [', U(vi,k), ',', V(vi,k), ']'
      WRITE(0,*) '           W = [', W(1), ',', W(2), ']'
    END IF
    
    ! Get derivatives on this triangle
    ddx = mesh%NxTri(tup,1) * d(mesh%Tri(tup,1),k) + mesh%NxTri(tup,2) * d(mesh%Tri(tup,2),k) + mesh%NxTri(tup,3) * d(mesh%Tri(tup,3),k)
    ddy = mesh%NyTri(tup,1) * d(mesh%Tri(tup,1),k) + mesh%NyTri(tup,2) * d(mesh%Tri(tup,2),k) + mesh%NyTri(tup,3) * d(mesh%Tri(tup,3),k)
    
  END SUBROUTINE get_upwind_derivative_vertex_3D
  
  SUBROUTINE apply_Neumann_boundary(    mesh, d)
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(INOUT)     :: d
    
    INTEGER                                       :: vi, ci, vc
    REAL(dp), DIMENSION(mesh%nC_mem)              :: vals
    INTEGER                                       :: nvals
        
    ! Apply Neumann boundary conditions: make sure ddx and ddy are zero at domain boundary
    ! Achieved by setting d for edge vertices to average value of non-edge neighbours
    ! ==========================================================================================
            
    DO vi = MAX(5,mesh%v1), mesh%v2   ! Edge vertices but not the four corners
      IF (mesh%edge_index(vi) == 0) CYCLE

      vals = 0._dp
      nvals = 0
  
      DO ci = 1, mesh%nC(vi)
        vc = mesh%C(vi,ci)
        IF (mesh%edge_index(vc)>0) CYCLE
        nvals = nvals+1
        vals(nvals) = d(vc)
      END DO
  
      d(vi) = SUM(vals) / nvals
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

    ! And the four corner vertices separately
    IF (par%master) THEN
      DO vi = 1, 4
  
        vals = 0._dp
        nvals = 0
  
        DO ci = 1, mesh%nC(vi)
          vc = mesh%C(vi,ci)
          nvals = nvals+1
          vals(nvals) = d(vc)
        END DO
  
        d(vi) = SUM(vals) / nvals
           
      END DO ! DO vi = 1, 4
    END IF
    CALL sync
    
  END SUBROUTINE apply_Neumann_boundary  
  SUBROUTINE apply_Neumann_boundary_3D( mesh, d, nz)
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:,:),   INTENT(INOUT)     :: d
    INTEGER,                    INTENT(IN)        :: nz
    
    INTEGER                                       :: vi, ci, vc, k
    REAL(dp), DIMENSION(mesh%nC_mem)             :: vals
    INTEGER                                       :: nvals
        
    ! Apply Neumann boundary conditions: make sure ddx and ddy are zero at domain boundary
    ! Achieved by setting d for edge vertices to average value of non-edge neighbours
    ! ==========================================================================================
            
    DO vi = MAX(5,mesh%v1), mesh%v2   ! Edge vertices but not the four corners
      IF (mesh%edge_index(vi) == 0) CYCLE
      
      DO k = 1, nz

        vals = 0._dp
        nvals = 0
  
        DO ci = 1, mesh%nC(vi)
          vc = mesh%C(vi,ci)
          IF (mesh%edge_index(vc)>0) CYCLE
          nvals = nvals+1
          vals(nvals) = d(vc,k)
        END DO
  
        d(vi,k) = SUM(vals) / nvals
      
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

    ! And the four corner vertices separately
    IF (par%master) THEN
      DO vi = 1, 4
      
        DO k = 1, nz
  
          vals = 0._dp
          nvals = 0
  
          DO ci = 1, mesh%nC(vi)
            vc = mesh%C(vi,ci)
            nvals = nvals+1
            vals(nvals) = d(vc,k)
          END DO
  
          d(vi,k) = SUM(vals) / nvals
        
        END DO     
      END DO ! DO vi = 1, 4
    END IF
    CALL sync
    
  END SUBROUTINE apply_Neumann_boundary_3D

END MODULE mesh_derivatives_module
