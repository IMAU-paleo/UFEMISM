MODULE mesh_derivatives_module
  ! Routines for calculating the neighbour functions on a mesh.

  ! Given a function f on the mesh, its first order derivatives on
  ! vertex i are called fx(i) and fy(i). These can be expressed as a linear
  ! combination of the values of f on i and its neighbours using the
  ! neighbour functions Nx and Ny, which depend only on mesh geometry:
  !
  ! fx(i) = Nx(i) * f(i) + sum_over_neighbours( Nx(c) * f(c))
  ! fy(i) = Ny(i) * f(i) + sum_over_neighbours( Ny(c) * f(c))
  !
  ! Similarly, the second order derivatives are found using the second
  ! order neighbour functions Nxx, Nxy and Nyy:
  !
  ! fxx(i) = Nxx(i) * f(i) + sum_over_neighbours( Nxx(c) * f(c))
  ! fxy(i) = Nxy(i) * f(i) + sum_over_neighbours( Nxy(c) * f(c))
  ! fyy(i) = Nyy(i) * f(i) + sum_over_neighbours( Nyy(c) * f(c))
  !
  ! To find out how neighbour functions are derived, see the original
  ! UFEMISM paper (Berends et al., GMD, 2020)

  USE mpi
  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync                
  USE data_types_module,           ONLY: type_mesh
  USE mesh_help_functions_module,  ONLY: IsInTriangle

  IMPLICIT NONE

  CONTAINS
  
  SUBROUTINE GetNeighbourFunctions( mesh)
    ! Calculate the first and second order neighbour functions    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER                                       :: vi, ti, iti, vn, n, cn, t1, t2, it1, it2
    REAL(dp)                                      :: ax, bx, cx, ay, by, cy, D
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: NxxTri, NxyTri, NyyTri, NxTri_ordered, NyTri_ordered
    REAL(dp), DIMENSION(:),   ALLOCATABLE         :: NxTri_ordered_t1, NxTri_ordered_t2, NyTri_ordered_t1, NyTri_ordered_t2
    REAL(dp), DIMENSION(3)                        :: NxTrisub, NyTrisub

    ! == First order on the triangles
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

    ! == First order on the vertices
    mesh%Nx(mesh%v1:mesh%v2,:) = 0._dp
    mesh%Ny(mesh%v1:mesh%v2,:) = 0._dp
    CALL sync

    DO vi = mesh%v1, mesh%v2

      ! Go through all triangles this vertex is a part of
      DO iti = 1, mesh%niTri(vi)

        ti = mesh%iTri(vi,iti)

        ! Go through that triangle's vertices
        DO n = 1, 3
          vn = mesh%Tri(ti,n)

          IF (vn==vi) THEN
            ! It is the home vertex
            mesh%Nx(vi,mesh%nC(vi)+1) = mesh%Nx(vi,mesh%nC(vi)+1) + mesh%NxTri(ti,n) / REAL(mesh%niTri(vi),dp)
            mesh%Ny(vi,mesh%nC(vi)+1) = mesh%Ny(vi,mesh%nC(vi)+1) + mesh%NyTri(ti,n) / REAL(mesh%niTri(vi),dp)
          ELSE
            ! Find this vertex's position in the home vertex's connectivity list
            DO cn = 1, mesh%nC(vi)
              IF (mesh%C(vi,cn) == vn) THEN
                mesh%Nx(vi,cn) = mesh%Nx(vi,cn) + mesh%NxTri(ti,n) / REAL(mesh%niTri(vi),dp)
                mesh%Ny(vi,cn) = mesh%Ny(vi,cn) + mesh%NyTri(ti,n) / REAL(mesh%niTri(vi),dp)
              END IF
            END DO
          END IF ! IF (vn==vi) THEN

        END DO ! DO n = 1, 3
      END DO ! DO iti = 1, mesh%niTri(vi)
      
    END DO ! DO vi = 1, mesh%nV
    CALL sync

    ! == Second order
    mesh%Nxx(mesh%v1:mesh%v2,:) = 0._dp
    mesh%Nxy(mesh%v1:mesh%v2,:) = 0._dp
    mesh%Nyy(mesh%v1:mesh%v2,:) = 0._dp
    CALL sync

    ALLOCATE(NxxTri(           mesh%nC_mem, mesh%nC_mem+1))
    ALLOCATE(NxyTri(           mesh%nC_mem, mesh%nC_mem+1))
    ALLOCATE(NyyTri(           mesh%nC_mem, mesh%nC_mem+1))
    ALLOCATE(NxTri_ordered(    mesh%nC_mem, mesh%nC_mem+1))
    ALLOCATE(NyTri_ordered(    mesh%nC_mem, mesh%nC_mem+1))
    ALLOCATE(NxTri_ordered_t1(               mesh%nC_mem+1))
    ALLOCATE(NxTri_ordered_t2(               mesh%nC_mem+1))
    ALLOCATE(NyTri_ordered_t1(               mesh%nC_mem+1))
    ALLOCATE(NyTri_ordered_t2(               mesh%nC_mem+1))

    DO vi = mesh%v1, mesh%v2
      ! Can't do this for edge vertices; skip them.
      IF (mesh%edge_index(vi)>0) CYCLE

      ! First for the subtriangles (see Documentation)
      NxxTri           = 0._dp
      NxyTri           = 0._dp
      NyyTri           = 0._dp
      NxTri_ordered    = 0._dp
      NyTri_ordered    = 0._dp
      NxTri_ordered_t1 = 0._dp
      NxTri_ordered_t2 = 0._dp
      NyTri_ordered_t1 = 0._dp
      NyTri_ordered_t2 = 0._dp

      ! Order NxTri and NyTri in terms of neighbour vertices like NxV and NyV!
      DO iti = 1, mesh%niTri(vi)
        ti = mesh%iTri(vi,iti)

        DO cn = 1, mesh%nC(vi)
          DO n = 1, 3
            IF (mesh%Tri(ti,n) == mesh%C(vi,cn)) THEN
              NxTri_ordered(iti,cn) = mesh%NxTri(ti,n)
              NyTri_ordered(iti,cn) = mesh%NyTri(ti,n)
            ELSE IF (mesh%Tri(ti,n) == vi) THEN
              NxTri_ordered(iti,mesh%nC(vi)+1) = mesh%NxTri(ti,n)
              NyTri_ordered(iti,mesh%nC(vi)+1) = mesh%NyTri(ti,n)
            END IF
          END DO ! DO n = 1, 3
        END DO ! DO cn = 1, mesh%nC(vi)
      END DO ! DO iti = 1, mesh%niTri(vi)

      ! NxTri_ordered: for every iTri, how fx on that triangle depends on all
      ! neighbours of vi

      DO iti = 1, mesh%niTri(vi)
        ! Get two adjacent triangles (t2 ahead of t1, CCW)
        it2 = iti
        it1 = iti-1
        IF (it1==0) it1 = mesh%niTri(vi)

        t1 = mesh%iTri(vi,it1)
        t2 = mesh%iTri(vi,it2)

        ! How fx and fy on t1 and t2 depend on f on all neighbours of vi and f
        ! on vi itself
        NxTri_ordered_t1 = NxTri_ordered(it1,:)
        NyTri_ordered_t1 = NyTri_ordered(it1,:)
        NxTri_ordered_t2 = NxTri_ordered(it2,:)
        NyTri_ordered_t2 = NyTri_ordered(it2,:)

        ! A: vertex vi, z determined by Nx(vi)
        ! B: midpoint of t1, z determined by Nxtri(t1)
        ! C: midpoint of t2, z determined by Nxtri(t2)
        ax = mesh%V(vi,1)
        ay = mesh%V(vi,2)
        bx = (mesh%V(mesh%Tri(t1,1),1) + mesh%V(mesh%Tri(t1,2),1) + mesh%V(mesh%Tri(t1,3),1))/3
        by = (mesh%V(mesh%Tri(t1,1),2) + mesh%V(mesh%Tri(t1,2),2) + mesh%V(mesh%Tri(t1,3),2))/3
        cx = (mesh%V(mesh%Tri(t2,1),1) + mesh%V(mesh%Tri(t2,2),1) + mesh%V(mesh%Tri(t2,3),1))/3
        cy = (mesh%V(mesh%Tri(t2,1),2) + mesh%V(mesh%Tri(t2,2),2) + mesh%V(mesh%Tri(t2,3),2))/3

        D = ax * (by - cy) + bx * (cy - ay) + cx * (ay - by);

        ! How fx and fy on subtriangle [V(vi), Tricc(t1), Tricc(t2)] depend
        ! on f on those three points
        NxTrisub = [(by - cy), (cy - ay), (ay - by)] / D
        NyTrisub = [(cx - bx), (ax - cx), (bx - ax)] / D

        ! How fxx, fxy and fyy on subtriangle [V(vi), Tricc(t1), Tricc(t2)]
        ! depend on f on all neighbours of vi and f on vi itself
        NxxTri(iti,:) =  NxTrisub(1) * mesh%Nx(vi,:) + &
                         NxTrisub(2) * NxTri_ordered_t1 + &
                         NxTrisub(3) * NxTri_ordered_t2
        NxyTri(iti,:) =  NxTrisub(1) * mesh%Ny(vi,:) + &
                         NxTrisub(2) * NyTri_ordered_t1 + &
                         NxTrisub(3) * NyTri_ordered_t2
        NyyTri(iti,:) =  NyTrisub(1) * mesh%Ny(vi,:) + &
                         NyTrisub(2) * NyTri_ordered_t1 + &
                         NyTrisub(3) * NyTri_ordered_t2
      END DO ! DO iti = 1, mesh%niTri(vi)

      ! How fxx, fxy and fyy on fi depend on f on all neighbours of vi and f
      ! on vi itself
      mesh%Nxx(vi,:) = SUM(NxxTri,1) / REAL(mesh%niTri(vi),dp)
      mesh%Nxy(vi,:) = SUM(NxyTri,1) / REAL(mesh%niTri(vi),dp)
      mesh%Nyy(vi,:) = SUM(NyyTri,1) / REAL(mesh%niTri(vi),dp)

    END DO ! DO vi = 1, mesh%nV
    CALL sync

    DEALLOCATE(NxxTri)
    DEALLOCATE(NxyTri)
    DEALLOCATE(NyyTri)
    DEALLOCATE(NxTri_ordered)
    DEALLOCATE(NyTri_ordered)
    DEALLOCATE(NxTri_ordered_t1)
    DEALLOCATE(NxTri_ordered_t2)
    DEALLOCATE(NyTri_ordered_t1)
    DEALLOCATE(NyTri_ordered_t2)

  END SUBROUTINE GetNeighbourFunctions  
  SUBROUTINE GetNeighbourFunctions_masked( mesh, mask)
    ! Get "masked" neighbour functions. Essentially the mesh equivalent of a one-sided
    ! differencing scheme; only average over those triangles that do not contain unmasked elements.
    ! Should make sure gradients at the ice margin are properly steep

    ! Assumes NxTri and NyTri have already been calculated.
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,  DIMENSION(:),     INTENT(IN)        :: mask

    INTEGER                                       :: vi, ti, iti, vn, n, cn, t1, t2, it1, it2
    REAL(dp)                                      :: ax, bx, cx, ay, by, cy, D
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: NxxTri, NxyTri, NyyTri, NxTri_ordered, NyTri_ordered
    REAL(dp), DIMENSION(:),   ALLOCATABLE         :: NxTri_ordered_t1, NxTri_ordered_t2, NyTri_ordered_t1, NyTri_ordered_t2
    REAL(dp), DIMENSION(3)                        :: NxTrisub, NyTrisub
    INTEGER,  DIMENSION(:),   ALLOCATABLE         :: niTri_masked
    INTEGER,  DIMENSION(:,:), ALLOCATABLE         :: iTri_masked
    LOGICAL                                       :: contains_unmasked_vertex

    ! Allocate memory
    ALLOCATE(niTri_masked(mesh%nV))
    ALLOCATE(iTri_masked( mesh%nV, mesh%nC_mem))

    niTri_masked = 0
    iTri_masked  = 0

    ! Determine which triangles should be used for each vertex
    DO vi = mesh%v1, mesh%v2
      IF (mask(vi)==1) THEN
        ! The vertex is masked; only use fully masked triangles (one-sided differencing)

        DO iti = 1, mesh%niTri(vi)
          ti = mesh%iTri(vi,iti)
          contains_unmasked_vertex = .FALSE.
          DO n = 1, 3
            vn = mesh%Tri(ti,n)
            IF (mask(vn)==0) THEN
              contains_unmasked_vertex = .TRUE.
            END IF
          END DO

          IF (.NOT. contains_unmasked_vertex) THEN
            niTri_masked(vi) = niTri_masked(vi)+1
            iTri_masked(vi,niTri_masked(vi)) = ti
          END IF

        END DO ! DO iti = 1, mesh%niTri(vi)

      ELSE ! IF (mask(vi)==1)
        ! The vertex is unmasked; use all triangles (two-sided differencing)

        niTri_masked(vi)  = mesh%niTri(vi)
        iTri_masked(vi,:) = mesh%iTri(vi,:)

      END IF
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

    ! == First order on the vertices
    mesh%Nxm(mesh%v1:mesh%v2,:) = 0._dp
    mesh%Nxm(mesh%v1:mesh%V2,:) = 0._dp
    CALL sync

    DO vi = mesh%v1, mesh%v2

      ! Go through all triangles this vertex is a part of
      DO iti = 1, niTri_masked(vi)

        ti = iTri_masked(vi,iti)

        ! Go through that triangle's vertices
        DO n = 1, 3
          vn = mesh%Tri(ti,n)

          IF (vn==vi) THEN
            ! It is the home vertex
            mesh%Nxm(vi,mesh%nC(vi)+1) = mesh%Nxm(vi,mesh%nC(vi)+1) + mesh%NxTri(ti,n) / REAL(mesh%niTri(vi),dp)
            mesh%Nym(vi,mesh%nC(vi)+1) = mesh%Nym(vi,mesh%nC(vi)+1) + mesh%NyTri(ti,n) / REAL(mesh%niTri(vi),dp)
          ELSE
            ! Find this vertex's position in the home vertex's connectivity list
            DO cn = 1, mesh%nC(vi)
              IF (mesh%C(vi,cn) == vn) THEN
                mesh%Nxm(vi,cn) = mesh%Nxm(vi,cn) + mesh%NxTri(ti,n) / REAL(mesh%niTri(vi),dp)
                mesh%Nym(vi,cn) = mesh%Nym(vi,cn) + mesh%NyTri(ti,n) / REAL(mesh%niTri(vi),dp)
              END IF
            END DO
          END IF ! IF (vn==vi) THEN

        END DO ! DO n = 1, 3

      END DO ! DO iti = 1, mesh%niTri(vi)

    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

    ! == Second order
    mesh%Nxxm(mesh%v1:mesh%v1,:) = 0._dp
    mesh%Nxym(mesh%v1:mesh%v1,:) = 0._dp
    mesh%Nyym(mesh%v1:mesh%v1,:) = 0._dp
    CALL sync

    ALLOCATE(NxxTri(       mesh%nC_mem, mesh%nC_mem+1))
    ALLOCATE(NxyTri(       mesh%nC_mem, mesh%nC_mem+1))
    ALLOCATE(NyyTri(       mesh%nC_mem, mesh%nC_mem+1))
    ALLOCATE(NxTri_ordered(mesh%nC_mem, mesh%nC_mem+1))
    ALLOCATE(NyTri_ordered(mesh%nC_mem, mesh%nC_mem+1))
    ALLOCATE(NxTri_ordered_t1(           mesh%nC_mem+1))
    ALLOCATE(NxTri_ordered_t2(           mesh%nC_mem+1))
    ALLOCATE(NyTri_ordered_t1(           mesh%nC_mem+1))
    ALLOCATE(NyTri_ordered_t2(           mesh%nC_mem+1))

    DO vi = mesh%v1, mesh%v2
      ! Can't do this for edge vertices; skip them.
      IF (mesh%edge_index(vi)>0) CYCLE

      ! First for the subtriangles
      NxxTri           = 0._dp
      NxyTri           = 0._dp
      NyyTri           = 0._dp
      NxTri_ordered    = 0._dp
      NyTri_ordered    = 0._dp
      NxTri_ordered_t1 = 0._dp
      NxTri_ordered_t2 = 0._dp
      NyTri_ordered_t1 = 0._dp
      NyTri_ordered_t2 = 0._dp

      ! Order NxTri and NyTri in terms of neighbour vertices like NxV and NyV!
      DO iti = 1, mesh%niTri(vi)
        ti = mesh%iTri(vi,iti)

        DO cn = 1, mesh%nC(vi)
          DO n = 1, 3
            IF (mesh%Tri(ti,n) == mesh%C(vi,cn)) THEN
              NxTri_ordered(iti,cn) = mesh%NxTri(ti,n)
              NyTri_ordered(iti,cn) = mesh%NyTri(ti,n)
            ELSE IF (mesh%Tri(ti,n) == vi) THEN
              NxTri_ordered(iti,mesh%nC(vi)+1) = mesh%NxTri(ti,n)
              NyTri_ordered(iti,mesh%nC(vi)+1) = mesh%NyTri(ti,n)
            END IF
          END DO ! DO n = 1, 3
        END DO ! DO cn = 1, mesh%nC(vi)
      END DO ! DO iti = 1, mesh%niTri(vi)

      ! NxTri_ordered: for every iTri, how fx on that triangle depends on all
      ! neighbours of vi

      IF (mask(vi)==0) THEN
        n = 1
      ELSE
        n = 2
      END IF

      DO iti = n, niTri_masked(vi)
        ! Get two adjacent triangles (t2 ahead of t1, CCW)
        it2 = iti
        it1 = iti-1
        IF (it1==0) it1 = niTri_masked(vi)

        t1 = iTri_masked(vi,it1)
        t2 = iTri_masked(vi,it2)

        ! How fx and fy on t1 and t2 depend on f on all neighbours of vi and f
        ! on vi itself
        NxTri_ordered_t1 = NxTri_ordered(it1,:)
        NyTri_ordered_t1 = NyTri_ordered(it1,:)
        NxTri_ordered_t2 = NxTri_ordered(it2,:)
        NyTri_ordered_t2 = NyTri_ordered(it2,:)

        ! A: vertex vi, z determined by Nx(vi)
        ! B: midpoint of t1, z determined by Nxtri(t1)
        ! C: midpoint of t2, z determined by Nxtri(t2)
        ax = mesh%V(vi,1)
        ay = mesh%V(vi,2)
        bx = (mesh%V(mesh%Tri(t1,1),1) + mesh%V(mesh%Tri(t1,2),1) + mesh%V(mesh%Tri(t1,3),1))/3
        by = (mesh%V(mesh%Tri(t1,1),2) + mesh%V(mesh%Tri(t1,2),2) + mesh%V(mesh%Tri(t1,3),2))/3
        cx = (mesh%V(mesh%Tri(t2,1),1) + mesh%V(mesh%Tri(t2,2),1) + mesh%V(mesh%Tri(t2,3),1))/3
        cy = (mesh%V(mesh%Tri(t2,1),2) + mesh%V(mesh%Tri(t2,2),2) + mesh%V(mesh%Tri(t2,3),2))/3

        D = ax * (by - cy) + bx * (cy - ay) + cx * (ay - by);

        ! How fx and fy on subtriangle [V(vi), Tricc(t1), Tricc(t2)] depend
        ! on f on those three points
        NxTrisub = [(by - cy), (cy - ay), (ay - by)] / D
        NyTrisub = [(cx - bx), (ax - cx), (bx - ax)] / D

        ! How fxx, fxy and fyy on subtriangle [V(vi), Tricc(t1), Tricc(t2)]
        ! depend on f on all neighbours of vi and f on vi itself
        NxxTri(iti,:) = NxTrisub(1) * mesh%Nx(vi,:) + &
                        NxTrisub(2) * NxTri_ordered_t1 + &
                        NxTrisub(3) * NxTri_ordered_t2
        NxyTri(iti,:) = NxTrisub(1) * mesh%Ny(vi,:) + &
                        NxTrisub(2) * NyTri_ordered_t1 + &
                        NxTrisub(3) * NyTri_ordered_t2
        NyyTri(iti,:) = NyTrisub(1) * mesh%Ny(vi,:) + &
                        NyTrisub(2) * NyTri_ordered_t1 + &
                        NyTrisub(3) * NyTri_ordered_t2
      END DO ! DO iti = 1, mesh%niTri(vi)

      ! How fxx, fxy and fyy on fi depend on f on all neighbours of vi and f
      ! on vi itself
      IF (niTri_masked(vi)>0) THEN
          mesh%Nxxm(vi,:) = SUM(NxxTri,1) / REAL(niTri_masked(vi),dp)
          mesh%Nxym(vi,:) = SUM(NxyTri,1) / REAL(niTri_masked(vi),dp)
          mesh%Nyym(vi,:) = SUM(NyyTri,1) / REAL(niTri_masked(vi),dp)
      END IF

    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

    DEALLOCATE(NxxTri)
    DEALLOCATE(NxyTri)
    DEALLOCATE(NyyTri)
    DEALLOCATE(NxTri_ordered)
    DEALLOCATE(NyTri_ordered)
    DEALLOCATE(NxTri_ordered_t1)
    DEALLOCATE(NxTri_ordered_t2)
    DEALLOCATE(NyTri_ordered_t1)
    DEALLOCATE(NyTri_ordered_t2)
    DEALLOCATE(niTri_masked)
    DEALLOCATE(iTri_masked)

  END SUBROUTINE GetNeighbourFunctions_masked
  
  SUBROUTINE GetMeshDerivatives(                  mesh, d, ddx, ddy)
    ! Get first order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: ddx, ddy

    INTEGER                                       :: vi, ci

    DO vi = mesh%v1, mesh%v2
      ddx( vi)   =            mesh%Nx( vi,mesh%nC( vi)+1) * d( vi)
      ddy( vi)   =            mesh%Ny( vi,mesh%nC( vi)+1) * d( vi)
      DO ci = 1, mesh%nC( vi)
        ddx( vi) = ddx( vi) + mesh%Nx( vi,ci) * d( mesh%C( vi,ci))
        ddy( vi) = ddy( vi) + mesh%Ny( vi,ci) * d( mesh%C( vi,ci))
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE GetMeshDerivatives  
  SUBROUTINE GetMeshDerivatives_3D(               mesh, d, ddx, ddy)
    ! Get first order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: d
    REAL(dp), DIMENSION(:,:),   INTENT(OUT)       :: ddx, ddy

    INTEGER                                       :: vi, ci

    DO vi = mesh%v1, mesh%v2
      ddx( vi,:)   =              mesh%Nx( vi,mesh%nC( vi)+1) * d( vi,:)
      ddy( vi,:)   =              mesh%Ny( vi,mesh%nC( vi)+1) * d( vi,:)
      DO ci = 1, mesh%nC( vi)
        ddx( vi,:) = ddx( vi,:) + mesh%Nx( vi,ci) * d( mesh%C( vi,ci),:)
        ddy( vi,:) = ddy( vi,:) + mesh%Ny( vi,ci) * d( mesh%C( vi,ci),:)
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE GetMeshDerivatives_3D
  SUBROUTINE GetMeshDerivativesVertex(            mesh, d, ddx, ddy, vi)
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

  END SUBROUTINE GetMeshDerivativesVertex  
  SUBROUTINE GetMeshDerivativesVertex_3D(         mesh, d, ddx, ddy, vi, k)
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

  END SUBROUTINE GetMeshDerivativesVertex_3D
  SUBROUTINE GetMeshCurvatures(                   mesh, d, ddxx, ddxy, ddyy)
    ! Get second order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: ddxx, ddxy,ddyy

    INTEGER                                       :: vi, ci

    DO vi = mesh%v1, mesh%v2
      ddxx( vi)    = mesh%Nxx( vi,mesh%nC( vi)+1) * d( vi)
      ddxy( vi)    = mesh%Nxy( vi,mesh%nC( vi)+1) * d( vi)
      ddyy( vi)    = mesh%Nyy( vi,mesh%nC( vi)+1) * d( vi)
      DO ci = 1, mesh%nC(vi)
        ddxx( vi)  = ddxx( vi) + mesh%Nxx( vi,ci) * d( mesh%C( vi,ci))
        ddxy( vi)  = ddxy( vi) + mesh%Nxy( vi,ci) * d( mesh%C( vi,ci))
        ddyy( vi)  = ddyy( vi) + mesh%Nyy( vi,ci) * d( mesh%C( vi,ci))
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE GetMeshCurvatures 
  SUBROUTINE GetMeshCurvaturesVertex(             mesh, d, ddxx, ddxy, ddyy, vi)
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

  END SUBROUTINE GetMeshCurvaturesVertex  
  SUBROUTINE GetMeshDerivatives_masked(           mesh, d, ddx, ddy)
    ! Get second order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: ddx, ddy

    INTEGER                                       :: vi, ci

    DO vi = mesh%v1, mesh%v2
      ddx( vi)    = mesh%Nxm( vi,mesh%nC( vi)+1) * d( vi)
      ddy( vi)    = mesh%Nym( vi,mesh%nC( vi)+1) * d( vi)
      DO ci = 1, mesh%nC( vi)
        ddx( vi)  = ddx( vi)  + mesh%Nxm( vi,ci) * d(mesh%C( vi,ci))
        ddy( vi)  = ddy( vi)  + mesh%Nym( vi,ci) * d(mesh%C( vi,ci))
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE GetMeshDerivatives_masked  
  SUBROUTINE GetMeshDerivatives_masked_vertex_3D( mesh, d, ddx, ddy, vi, k)
    ! Get second order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: d
    REAL(dp),                   INTENT(OUT)       :: ddx, ddy
    INTEGER,                    INTENT(IN)        :: vi,k

    INTEGER                                       :: ci

    ddx    = mesh%Nxm( vi,mesh%nC( vi)+1) * d( vi,k)
    ddy    = mesh%Nym( vi,mesh%nC( vi)+1) * d( vi,k)

    DO ci = 1, mesh%nC( vi)
      ddx  = ddx  + mesh%Nxm( vi,ci) * d( mesh%C( vi,ci),k)
      ddy  = ddy  + mesh%Nym( vi,ci) * d( mesh%C( vi,ci),k)
    END DO

  END SUBROUTINE GetMeshDerivatives_masked_vertex_3D 
  SUBROUTINE GetMeshCurvatures_masked(            mesh, d, ddxx, ddxy, ddyy)
    ! Get second order spatial derivatives of a function on the mesh
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: ddxx, ddxy, ddyy

    INTEGER                                       :: vi, ci

    DO vi = mesh%V1, mesh%v2
      ddxx( vi) = mesh%Nxxm( vi,mesh%nC( vi)+1) * d( vi)
      ddxy( vi) = mesh%Nxym( vi,mesh%nC( vi)+1) * d( vi)
      ddyy( vi) = mesh%Nyym( vi,mesh%nC( vi)+1) * d( vi)
      DO ci = 1, mesh%nC( vi)
        ddxx( vi) = ddxx( vi) + mesh%Nxxm( vi,ci) * d( mesh%C( vi,ci))
        ddxy( vi) = ddxy( vi) + mesh%Nxym( vi,ci) * d( mesh%C( vi,ci))
        ddyy( vi) = ddyy( vi) + mesh%Nyym( vi,ci) * d( mesh%C( vi,ci))
      END DO
    END DO ! DO vi = mesh%V1, mesh%v2
    CALL sync

  END SUBROUTINE GetMeshCurvatures_masked  
  
  SUBROUTINE GetUpwindDerivativeVertex_3D( mesh, U, V, d, vi, k, ddx, ddy)
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
      CALL GetMeshDerivativesVertex_3D( mesh, d, ddx, ddy, vi, k)
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
      IF (IsInTriangle(pa, pb, pc, p)) THEN
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
    
  END SUBROUTINE GetUpwindDerivativeVertex_3D
  
  SUBROUTINE ApplyNeumannBoundary( mesh, d)
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(INOUT)     :: d
    
    INTEGER                                       :: vi, ci, vc
    REAL(dp), DIMENSION(mesh%nC_mem)             :: vals
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
    
  END SUBROUTINE ApplyNeumannBoundary  
  SUBROUTINE ApplyNeumannBoundary_3D( mesh, d, nz)
    
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
    
  END SUBROUTINE ApplyNeumannBoundary_3D

END MODULE mesh_derivatives_module
