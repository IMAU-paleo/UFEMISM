MODULE mesh_help_functions_module
  ! General help functions used in mesh creation and updating.

  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync, allocate_shared_memory_int_0D, allocate_shared_memory_dp_0D, &
                                                    allocate_shared_memory_int_1D, allocate_shared_memory_dp_1D, &
                                                    allocate_shared_memory_int_2D, allocate_shared_memory_dp_2D, &
                                                    allocate_shared_memory_int_3D, allocate_shared_memory_dp_3D, &
                                                    allocate_shared_memory_bool_1D, deallocate_shared_memory, &
                                                    adapt_shared_memory_int_1D,    adapt_shared_memory_dp_1D, &
                                                    adapt_shared_memory_int_2D,    adapt_shared_memory_dp_2D, &
                                                    adapt_shared_memory_int_3D,    adapt_shared_memory_dp_3D, &
                                                    adapt_shared_memory_bool_1D                   
  USE data_types_module,           ONLY: type_mesh, type_model_region

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'

  CONTAINS
    
! == Calculating some extra mesh data: Voronoi cell areas, connection widths,
!    triangle areas, resolution, lat/lon-coordinates
  SUBROUTINE FindVoronoiCellAreas( mesh)
    ! Find the areas of the Voronoi cells of all the vertices

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: vi, nVor, n
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Vor
    REAL(dp)                                      :: Aerr
    
    ALLOCATE(Vor(mesh%nconmax+2,2))

    mesh%A(mesh%v1:mesh%v2) = 0._dp
    DO vi = mesh%v1, mesh%v2

      CALL FindVoronoiCellVertices(mesh, vi, Vor, nVor)

      DO n = 2, nVor
        mesh%A(vi) = mesh%A(vi) + norm2(Cross( &
        [Vor(n  ,1)-mesh%V(vi,1), Vor(n  ,2)-mesh%V(vi,2), 0._dp], &
        [Vor(n-1,1)-mesh%V(vi,1), Vor(n-1,2)-mesh%V(vi,2), 0._dp] )) / 2._dp
      END DO
    END DO
    CALL sync

    DEALLOCATE(Vor)

    ! Check if everything went alright
    IF (par%master) THEN
      Aerr = ABS(1._dp - SUM(mesh%A ) / ((mesh%xmax-mesh%xmin)*(mesh%ymax-mesh%ymin)))
      IF (Aerr > 0.0001_dp) WRITE(0,*) 'FindVoronoiCellAreas - WARNING: sum of Voronoi cell areas doesnt match square area of mesh! (error of ', Aerr/100._dp ,' %)'
    END IF
    
  END SUBROUTINE FindVoronoiCellAreas 
  SUBROUTINE FindVoronoiCellGeometricCentres( mesh)
    ! Find the geometric centres of the Voronoi cells of all the vertices

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables
    INTEGER                                       :: vi, nvi
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Vor
    INTEGER                                       :: nVor
    REAL(dp), DIMENSION(2)                        :: p, q
    REAL(dp)                                      :: LI_mxydx, LI_xydy
        
    mesh%VorGC( mesh%v1:mesh%V2,:) = 0._dp
    
    ALLOCATE(Vor(mesh%nconmax+2,2))
    
    DO vi = mesh%v1, mesh%v2
    
      CALL FindVoronoiCellVertices(mesh, vi, Vor, nVor)
      
      LI_mxydx = 0._dp
      LI_xydy  = 0._dp
      
      DO nvi = 2, nVor
        p = Vor( nvi-1,:)
        q = Vor( nvi,:)
        LI_mxydx = LI_mxydx + LineIntegral_mxydx( p, q, mesh%tol_dist)
        LI_xydy  = LI_xydy  + LineIntegral_xydy(  p, q, mesh%tol_dist)    
      END DO 
      
      IF (mesh%edge_index(vi)>0) THEN
        p = Vor( nVor,:)
        q = mesh%V(vi,:)
        LI_mxydx = LI_mxydx + LineIntegral_mxydx( p, q, mesh%tol_dist)
        LI_xydy  = LI_xydy  + LineIntegral_xydy(  p, q, mesh%tol_dist)
        p = mesh%V(vi,:)
        q = Vor(1,:)
        LI_mxydx = LI_mxydx + LineIntegral_mxydx( p, q, mesh%tol_dist)
        LI_xydy  = LI_xydy  + LineIntegral_xydy(  p, q, mesh%tol_dist)
      END IF
       
      mesh%VorGC(vi,:) = [LI_mxydx/mesh%A(vi), LI_xydy/mesh%A(vi)]
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    DEALLOCATE(Vor)
    
  END SUBROUTINE FindVoronoiCellGeometricCentres
  SUBROUTINE FindConnectionWidths(mesh)
    ! Find the width of the line separating two connected vertices (equal to the distance
    ! between the circumcenters of their two shared triangles)

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: v1, nv2, v2, t1, t2, iti, ti, n
    LOGICAL                                       :: hasv2
    
    mesh%Cw(mesh%v2:mesh%v2,:) = 0._dp
    CALL sync

    ! The way to go: for each vertex pair, find the two triangles that both
    ! are a part of. If there's only one, there's one Voronoi vertex and its
    ! projection on the map edge.
    DO v1 = mesh%v1, mesh%v2

      DO nv2 = 1, mesh%nC(v1)
        v2 = mesh%C(v1,nv2)

        t1 = 0
        t2 = 0
        DO iti = 1, mesh%niTri(v1)
          ti = mesh%iTri(v1,iti)
          hasv2 = .FALSE.
          DO n = 1, 3
            IF (mesh%Tri(ti,n)==v2) hasv2 = .TRUE.
          END DO
          IF (hasv2) THEN
            IF (t1==0) THEN
              t1 = ti
            ELSE
              t2 = ti
            END IF
          END IF
        END DO ! DO iti = 1, mesh%niTri(v1)

        ! We should now have at least a single shared triangle.
        IF (t1==0) THEN
          WRITE(0,*) 'FindConnectionWidths: ERROR - Couldnt find a single shared triangle!'
          STOP
        END IF

        IF (t2>0) THEN
          ! Two shared triangles; Cw equals the distance between the two
          ! triangles' circumcenters
          mesh%Cw(v1,nv2) = norm2(mesh%Tricc(t1,:) - mesh%Tricc(t2,:))
        ELSE
          ! One shared triangle; Cw equals the distance between that triangle's
          ! circumcenter and the domain boundary
          IF     (mesh%Tri_edge_index(t1)==1) THEN
            mesh%Cw(v1,nv2) = MAX(0._dp, mesh%ymax - mesh%Tricc(t1,2))
          ELSEIF (mesh%Tri_edge_index(t1)==3) THEN
            mesh%Cw(v1,nv2) = MAX(0._dp, mesh%xmax - mesh%Tricc(t1,1))
          ELSEIF (mesh%Tri_edge_index(t1)==5) THEN
            mesh%Cw(v1,nv2) = MAX(0._dp, mesh%Tricc(t1,2) - mesh%ymin)
          ELSEIF (mesh%Tri_edge_index(t1)==7) THEN
            mesh%Cw(v1,nv2) = MAX(0._dp, mesh%Tricc(t1,1) - mesh%xmin)
          ELSE
            WRITE(0,*) 'FindConnectionWidths: ERROR - The only shared triangle isnt a Boundary triangle - this cannot be!'
            STOP
          END IF
        END IF ! IF (t2>0) THEN

      END DO ! DO nv2 = 1, mesh%nC(v1)
    END DO ! DO v1 = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE FindConnectionWidths   
  SUBROUTINE FindTriangleAreas(mesh)
    ! Find the areas of all the triangles

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: ti
    REAL(dp), DIMENSION(2)                        :: pa, pb, pc

    DO ti = mesh%t1, mesh%t2
      pa = mesh%V(mesh%Tri(ti,1),:)
      pb = mesh%V(mesh%Tri(ti,2),:)
      pc = mesh%V(mesh%Tri(ti,3),:)
      mesh%TriA(ti) = FindTriangleArea(pa,pb,pc)
    END DO ! DO ti = mesh%t1, mesh%t2
    CALL sync

  END SUBROUTINE FindTriangleAreas
  SUBROUTINE DetermineMeshResolution(mesh)
    ! Find the areas of all the triangles

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: vi, vj, ci
    
    mesh%R(mesh%v1:mesh%v2) = mesh%xmax - mesh%xmin
    
    DO vi = mesh%v1, mesh%v2
      DO ci = 1, mesh%nC(vi)
        vj = mesh%C(vi,ci)
        mesh%R(vi) = MIN(mesh%R(vi), SQRT((mesh%V(vj,1)-mesh%V(vi,1))**2 + (mesh%V(vj,2)-mesh%V(vi,2))**2))
      END DO
    END DO
    CALL sync
      
    mesh%resolution_min = MINVAL(mesh%R)
    mesh%resolution_max = MAXVAL(mesh%R)

  END SUBROUTINE DetermineMeshResolution  
  SUBROUTINE CalculateDynamicOmega(mesh)
    ! Calculate the relaxation parameter omega, used in SOR solvers like the SSA, as a function of mesh geometry
    
    USE parameters_module, ONLY : n_flow, q_plastic, u_threshold, ice_density, grav, epsilon_sq_0, delta_v
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:    
    REAL(dp), PARAMETER                                :: alpha_min_min = 0.3_dp
    REAL(dp), PARAMETER                                :: alpha_min_max = 1.1_dp
    REAL(dp), PARAMETER                                :: omega_min = 1.0_dp
    REAL(dp), PARAMETER                                :: omega_max = 1.4_dp
    REAL(dp), DIMENSION(:), POINTER                    :: Tri_alpha_min
    INTEGER                                            :: wTri_alpha_min, ti, vp, vq, vr
    REAL(dp), DIMENSION(2)                             :: p, q, r, pq, qr, rp
    REAL(dp)                                           :: ap, aq, ar, omega_dyn
    
    mesh%omega(mesh%v1:mesh%v2) = 100.0_dp
  
    CALL allocate_shared_memory_dp_1D( mesh%nTri, Tri_alpha_min, wTri_alpha_min)
    Tri_alpha_min(mesh%t1:mesh%t1) = 0._dp
    
    DO ti = mesh%t1, mesh%t2
    
      ! Triangle vertex indices
      vp = mesh%Tri(ti,1)
      vq = mesh%Tri(ti,2)
      vr = mesh%Tri(ti,3)
    
      ! Triangle vertex coordinates
      p = mesh%V(vp,:)
      q = mesh%V(vq,:)
      r = mesh%V(vr,:)
          
      ! Triangle legs
      pq = p-q
      qr = q-r
      rp = r-p
  
      ! Internal angles
      ap = ACOS(-(rp(1)*pq(1) + rp(2)*pq(2))/(NORM2(rp)*NORM2(pq)))
      aq = ACOS(-(pq(1)*qr(1) + pq(2)*qr(2))/(NORM2(pq)*NORM2(qr)))
      ar = ACOS(-(rp(1)*qr(1) + rp(2)*qr(2))/(NORM2(rp)*NORM2(qr)))
      
      ! Smallest internal angle
      Tri_alpha_min(ti) = MINVAL([ap, aq, ar])
      
      ! Dynamic omega
      omega_dyn = MIN(omega_max, MAX(omega_min, omega_min + (Tri_alpha_min(ti) - alpha_min_min) * (omega_max - omega_min) / (alpha_min_max - alpha_min_min)))
      
      mesh%omega(vp) = MIN(mesh%omega(vp), omega_dyn)
      mesh%omega(vq) = MIN(mesh%omega(vq), omega_dyn)
      mesh%omega(vr) = MIN(mesh%omega(vr), omega_dyn)
            
    END DO ! DO ti = mesh%t1, mesh%t2
    CALL sync
    
    CALL deallocate_shared_memory( wTri_alpha_min)
    NULLIFY( Tri_alpha_min)
      
  END SUBROUTINE CalculateDynamicOmega
  
! == Finding the vertices of a vertex' Voronoi cell
  SUBROUTINE FindVoronoiCellVertices(mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a vertex's Voronoi cell

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: vi
    REAL(dp), DIMENSION(:,:), INTENT(INOUT)       :: Vor
    INTEGER,                  INTENT(OUT)         :: nVor
    
    ! Local variables
    INTEGER                                       :: vvi

    IF (mesh%edge_index(vi)==0) THEN
      CALL FindVoronoiCellVertices_free(mesh, vi, Vor, nVor)
    ELSEIF (mesh%edge_index(vi)==2 .OR. mesh%edge_index(vi)==4 .OR. mesh%edge_index(vi)==6 .OR. mesh%edge_index(vi)==8) THEN
      CALL FindVoronoiCellVertices_corner(mesh, vi, Vor, nVor)
    ELSE
      CALL FindVoronoiCellVertices_edge(mesh, vi, Vor, nVor)
    END IF
    
    DO vvi = 1, nVor
      IF (Vor(vvi,1) < mesh%xmin - mesh%tol_dist .OR. &
          Vor(vvi,1) > mesh%xmax + mesh%tol_dist .OR. &
          Vor(vvi,2) < mesh%ymin - mesh%tol_dist .OR. &
          Vor(vvi,2) > mesh%ymax + mesh%tol_dist) THEN
        WRITE(0,*) '   FindVoronoiCellVertices - ERROR: found Voronoi cell vertex outside of mesh domain!'
        STOP
      END IF
      Vor(vvi,1) = MAX( MIN( Vor(vvi,1), mesh%xmax), mesh%xmin)
      Vor(vvi,2) = MAX( MIN( Vor(vvi,2), mesh%ymax), mesh%ymin)
    END DO

  END SUBROUTINE FindVoronoiCellVertices 
  SUBROUTINE FindVoronoiCellVertices_free(mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a free vertex's Voronoi cell

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: vi
    REAL(dp), DIMENSION(:,:), INTENT(INOUT)       :: Vor
    INTEGER,                  INTENT(OUT)         :: nVor

    ! Local variables
    INTEGER                                       :: iti, iti_clock, iti_anti, ti, ti_clock, ti_anti
    REAL(dp), DIMENSION(2)                        :: cc, cc_clock, cc_anti

    Vor  = 0._dp
    nVor = 0

    DO iti = 1, mesh%niTri(vi)

      ! Find indices of current, clockwise neighbouring and anticlockwise
      ! neighbouring triangles.
      iti_clock = iti - 1
      IF (iti_clock==0) iti_clock = mesh%niTri(vi)
      iti_anti  = iti + 1
      IF (iti_anti>mesh%niTri(vi)) iti_anti = 1

      ti       = mesh%iTri(vi,iti)
      ti_clock = mesh%iTri(vi,iti_clock)
      ti_anti  = mesh%iTri(vi,iti_anti)

      ! If necessary, crop (split) the circumcenter of the current triangle.
      cc       = mesh%Tricc(ti,:)
      CALL CropCircumcenter(mesh, ti, ti_clock, cc_clock)
      CALL CropCircumcenter(mesh, ti, ti_anti,  cc_anti)

      ! Add the resulting Voronoi vertex/vertices
      IF (cc_clock(1) /= cc_anti(1) .OR. cc_clock(2) /= cc_anti(2)) THEN
        nVor = nVor+1
        Vor(nVor,:) = cc_clock
        nVor = nVor+1
        Vor(nVor,:) = cc_anti
      ELSE
        nVor = nVor+1
        Vor(nVor,:) = cc
      END IF

    END DO ! DO t = 1, mesh%niTri(vi)

    ! Repeat the first Voronoi vertex
    nVor = nVor+1
    Vor(nVor,:) = Vor(1,:)

  END SUBROUTINE FindVoronoiCellVertices_free
  SUBROUTINE FindVoronoiCellVertices_edge(mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up an edge vertex's Voronoi cell

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: vi
    REAL(dp), DIMENSION(:,:), INTENT(INOUT)       :: Vor
    INTEGER,                  INTENT(OUT)         :: nVor

    ! Local variables
    INTEGER                                       :: iti
    REAL(dp), DIMENSION(2)                        :: cc, cc_clock, cc_anti, cc_cropped

    Vor  = 0._dp
    nVor = 0

    ! == Boundary cell ==
    ! If the first or last circumcenter lies outside of the grid, crop it.
    ! If not, add the point on the edge closest to that circumcenter as an additional Voronoi cell vertex.

    DO iti = 1, mesh%niTri(vi)

      cc = mesh%Tricc(mesh%iTri(vi,iti),:)

      IF (iti == 1) THEN
        ! Start by possibly adding the boundary projection of the vertex
        IF     ((mesh%edge_index(vi)==1 .OR. mesh%edge_index(vi)==2) .AND. cc(2)<mesh%ymax) THEN
          nVor = nVor+1
          Vor(nVor,:) = [cc(1), mesh%ymax]
        ELSEIF ((mesh%edge_index(vi)==3 .OR. mesh%edge_index(vi)==4) .AND. cc(1)<mesh%xmax) THEN
          nVor = nVor+1
          Vor(nVor,:) = [mesh%xmax, cc(2)]
        ELSEIF ((mesh%edge_index(vi)==5 .OR. mesh%edge_index(vi)==6) .AND. cc(2)>mesh%ymin) THEN
          nVor = nVor+1
          Vor(nVor,:) = [cc(1), mesh%ymin]
        ELSEIF ((mesh%edge_index(vi)==7 .OR. mesh%edge_index(vi)==8) .AND. cc(1)>mesh%xmin) THEN
          nVor = nVor+1
          Vor(nVor,:) = [mesh%xmin, cc(2)]
        END IF

        ! Then add the (possibly cropped) vertex
        CALL CropCircumcenter(mesh, mesh%iTri(vi,1), mesh%iTri(vi,2), cc_cropped)
        nVor = nVor+1
        Vor(nVor,:) = cc_cropped
      END IF ! IF (iti == 1) THEN

      IF (iti > 1 .AND. iti < mesh%niTri(vi)) THEN
        ! Split the circumcenter
        CALL CropCircumcenter(mesh, mesh%iTri(vi,iti), mesh%iTri(vi,iti+1), cc_anti)
        CALL CropCircumcenter(mesh, mesh%iTri(vi,iti), mesh%iTri(vi,iti-1), cc_clock)
        IF (cc_anti(1)/=cc_clock(1) .OR. cc_anti(2)/=cc_clock(2)) THEN
          nVor = nVor+1
          Vor(nVor,:) = cc_clock
          nVor = nVor+1
          Vor(nVor,:) = cc_anti
        ELSE
          nVor = nVor+1
          Vor(nVor,:) = cc
        END IF
      END IF ! IF (iti > 1 .AND. iti < mesh%niTri(vi)) THEN

      IF (iti == mesh%niTri(vi)) THEN
        ! First add the (possibly cropped) vertex
        CALL CropCircumcenter(mesh, mesh%iTri(vi,iti), mesh%iTri(vi,iti-1), cc_cropped)
        nVor = nVor+1
        Vor(nVor,:) = cc_cropped

        ! Then possibly add the boundary projection of the vertex
        IF     ((mesh%edge_index(vi)==1 .OR. mesh%edge_index(vi)==8) .AND. cc(2)<mesh%ymax) THEN
          nVor = nVor+1
          Vor(nVor,:) = [cc(1), mesh%ymax]
        ELSEIF ((mesh%edge_index(vi)==3 .OR. mesh%edge_index(vi)==2) .AND. cc(1)<mesh%xmax) THEN
          nVor = nVor+1
          Vor(nVor,:) = [mesh%xmax, cc(2)]
        ELSEIF ((mesh%edge_index(vi)==5 .OR. mesh%edge_index(vi)==4) .AND. cc(2)>mesh%ymin) THEN
          nVor = nVor+1
          Vor(nVor,:) = [cc(1), mesh%ymin]
        ELSEIF ((mesh%edge_index(vi)==7 .OR. mesh%edge_index(vi)==6) .AND. cc(1)>mesh%xmin) THEN
          nVor = nVor+1
          Vor(nVor,:) = [mesh%xmin, cc(2)]
        END IF

      END IF ! IF (iti == mesh%niTri(vi)) THEN

    END DO ! DO n = 1, mesh%niTri(vi)

  END SUBROUTINE FindVoronoiCellVertices_edge
  SUBROUTINE FindVoronoiCellVertices_corner(mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a corner vertex's Voronoi cell

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: vi
    REAL(dp), DIMENSION(:,:), INTENT(INOUT)       :: Vor
    INTEGER,                  INTENT(OUT)         :: nVor

    ! Local variables
    REAL(dp), DIMENSION(2)                        :: cc

    Vor  = 0._dp
    nVor = 0
    
    IF (mesh%niTri(vi) > 1) THEN
      ! This corner vertex has more than one triangle, can be handled by Edge version  
          
      CALL FindVoronoiCellVertices_edge(mesh, vi, Vor, nVor)
      
    ELSE
      ! This corner vertex has only a single triangle, best handled manually

      cc = mesh%Tricc(mesh%iTri(vi,1),:)

      IF     (mesh%edge_index(vi)==2) THEN
        ! Northeast corner
        nVor = 3
        Vor(1,:) = [cc(1), mesh%ymax]
        Vor(2,:) = cc
        Vor(3,:) = [mesh%xmax, cc(2)]
      ELSEIF (mesh%edge_index(vi)==4) THEN
        ! Southeast corner
        nVor = 3
        Vor(1,:) = [mesh%xmax, cc(2)]
        Vor(2,:) = cc
        Vor(3,:) = [cc(1), mesh%ymin]
      ELSEIF (mesh%edge_index(vi)==6) THEN
        ! Southwest corner
        nVor = 3
        Vor(1,:) = [cc(1), mesh%ymin]
        Vor(2,:) = cc
        Vor(3,:) = [mesh%xmin, cc(2)]
      ELSEIF (mesh%edge_index(vi)==8) THEN
        ! Northwest corner
        nVor = 3
        Vor(1,:) = [mesh%xmin, cc(2)]
        Vor(2,:) = cc
        Vor(3,:) = [cc(1), mesh%ymax]
      ELSE
        WRITE(0,*) 'A non-corner vertex has only one triangle? This cannot be!'
        STOP
      END IF ! IF (mesh%edge_index(vi)==2) THEN
      
    END IF

  END SUBROUTINE FindVoronoiCellVertices_corner
  SUBROUTINE CropCircumcenter(mesh,t1,t2,ccc)
    ! Crop the circumcenter of triangle t1 in the direction of t2

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: t1, t2
    REAL(dp), DIMENSION(2),   INTENT(OUT)         :: ccc

    REAL(dp)                                      :: la, lb, lc, le, lf, lg
    REAL(dp), DIMENSION(2)                        :: p, q

    ccc  = mesh%Tricc(t1,:)

    IF     (mesh%Tri_edge_index(t1)==1 .AND. mesh%Tricc(t1,2)>mesh%ymax) THEN
      ! North boundary triangle
      CALL LineFromPoints([mesh%xmin, mesh%ymax], [mesh%xmax, mesh%ymax], la, lb, lc)
    ELSEIF (mesh%Tri_edge_index(t1)==3 .AND. mesh%Tricc(t1,1)>mesh%xmax) THEN
      ! East boundary triangle
      CALL LineFromPoints([mesh%xmax, mesh%ymax], [mesh%xmax, mesh%ymin], la, lb, lc)
    ELSEIF (mesh%Tri_edge_index(t1)==5 .AND. mesh%Tricc(t1,2)<mesh%ymin) THEN
      ! South boundary triangle
      CALL LineFromPoints([mesh%xmin, mesh%ymin], [mesh%xmax, mesh%ymin], la, lb, lc)
    ELSEIF (mesh%Tri_edge_index(t1)==7 .AND. mesh%Tricc(t1,1)<mesh%xmin) THEN
      ! West boundary triangle
      CALL LineFromPoints([mesh%xmin, mesh%ymax], [mesh%xmin, mesh%ymin], la, lb, lc)
    ELSE
      RETURN
    END IF

   p = mesh%Tricc(t1,:)
   q = mesh%Tricc(t2,:)
   CALL LineFromPoints(p, q, le, lf, lg)
   CALL LineLineIntersection(la, lb, lc, le, lf, lg, ccc);

  END SUBROUTINE CropCircumcenter 
  
! == The oblique stereographic projection
  SUBROUTINE GetLatLonCoordinates(mesh)
   ! Use the inverse stereographic projection for the mesh model region to calculate
   ! lat/lon coordinates for all the vertices

    TYPE(type_mesh),          INTENT(INOUT)       :: mesh
    
    INTEGER                                       :: vi
    
    ! Calculate lat and lon directly from X and Y using inverse projection
    IF     (mesh%region_name == 'NAM') THEN
      mesh%lambda_M               = C%lambda_M_NAM
      mesh%phi_M                  = C%phi_M_NAM
      mesh%alpha_stereographic    = C%alpha_NAM
    ELSEIF (mesh%region_name == 'EAS') THEN
      mesh%lambda_M               = C%lambda_M_EAS
      mesh%phi_M                  = C%phi_M_EAS
      mesh%alpha_stereographic    = C%alpha_EAS
    ELSEIF (mesh%region_name == 'GRL') THEN
      mesh%lambda_M               = C%lambda_M_GRL
      mesh%phi_M                  = C%phi_M_GRL
      mesh%alpha_stereographic    = C%alpha_GRL
    ELSEIF (mesh%region_name == 'ANT') THEN
      mesh%lambda_M               = C%lambda_M_ANT
      mesh%phi_M                  = C%phi_M_ANT
      mesh%alpha_stereographic    = C%alpha_ANT
    END IF    
  
    DO vi = mesh%v1, mesh%v2
      CALL inverse_oblique_sg_projection(mesh%V(vi,1), mesh%V(vi,2), mesh%lambda_M, mesh%phi_M, mesh%alpha_stereographic, mesh%lon(vi), mesh%lat(vi))
    END DO
    CALL sync
    
  END SUBROUTINE GetLatLonCoordinates  
  SUBROUTINE oblique_sg_projection(lambda, phi, lambda_M_deg, phi_M_deg, alpha_deg, x_IM_P_prime, y_IM_P_prime)
    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates which coincide with the GCM grid points to the rectangular IM coordinate
    ! system, with coordinates (x,y).
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)            :: lambda        ! lon in degrees
    REAL(dp), INTENT(IN)            :: phi           ! lat in degrees
    
    ! Polar stereographic projection parameters
    REAL(dp), INTENT(IN)  :: lambda_M_deg  ! in degrees
    REAL(dp), INTENT(IN)  :: phi_M_deg     ! in degrees
    REAL(dp), INTENT(IN)  :: alpha_deg     ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT)           :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT)           :: y_IM_P_prime  ! in meter

    ! Local variables:
    REAL(dp)                        :: phi_P         ! in radians
    REAL(dp)                        :: lambda_P      ! in radians
    REAL(dp)                        :: t_P_prime
    
    REAL(dp)                        :: lambda_M, phi_M, alpha
    
    lambda_M = C%deg2rad * lambda_M_deg
    phi_M    = C%deg2rad * phi_M_deg
    alpha    = C%deg2rad * alpha_deg

    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = C%deg2rad * phi
    lambda_P = C%deg2rad * lambda

    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = ((1._dp + COS(alpha)) / (1._dp + COS(phi_P) * COS(phi_M) * COS(lambda_P - lambda_M) + SIN(phi_P) * SIN(phi_M))) / C%deg2rad

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x_IM_P_prime =  C%earth_radius * (COS(phi_P) * SIN(lambda_P - lambda_M)) * t_P_prime
    y_IM_P_prime =  C%earth_radius * (SIN(phi_P) * COS(phi_M) - (COS(phi_P) * SIN(phi_M)) * COS(lambda_P - lambda_M)) * t_P_prime   
    
    
  END SUBROUTINE oblique_sg_projection
  SUBROUTINE inverse_oblique_sg_projection(x_IM_P_prime, y_IM_P_prime, lambda_M_deg, phi_M_deg, alpha_deg, lambda_P, phi_P)
    ! This subroutine projects with an inverse oblique stereographic projection the
    ! (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, alpha, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter
    
    ! Polar stereographic projection parameters
    REAL(dp), INTENT(IN)  :: lambda_M_deg  ! in degrees
    REAL(dp), INTENT(IN)  :: phi_M_deg     ! in degrees
    REAL(dp), INTENT(IN)  :: alpha_deg     ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: x_3D_P_prime  ! in meter
    REAL(dp)              :: y_3D_P_prime  ! in meter
    REAL(dp)              :: z_3D_P_prime  ! in meter
    REAL(dp)              :: a
    REAL(dp)              :: t_P
    REAL(dp)              :: x_3D_P        ! in meter
    REAL(dp)              :: y_3D_P        ! in meter
    REAL(dp)              :: z_3D_P        ! in meter
    
    REAL(dp)              :: lambda_M, phi_M, alpha
    
    lambda_M = C%deg2rad * lambda_M_deg
    phi_M    = C%deg2rad * phi_M_deg
    alpha    = C%deg2rad * alpha_deg

    ! See equations (2.14-2.16) or equations (B.21-B.23) in Reerink et al. (2010):
    x_3D_P_prime = C%earth_radius * COS(alpha) * COS(lambda_M) * COS(phi_M) - SIN(lambda_M) * (x_IM_P_prime*C%deg2rad) - COS(lambda_M) * SIN(phi_M) * (y_IM_P_prime*C%deg2rad)
    y_3D_P_prime = C%earth_radius * COS(alpha) * SIN(lambda_M) * COS(phi_M) + COS(lambda_M) * (x_IM_P_prime*C%deg2rad) - SIN(lambda_M) * SIN(phi_M) * (y_IM_P_prime*C%deg2rad)
    z_3D_P_prime = C%earth_radius * COS(alpha) *                 SIN(phi_M)                                          +                   COS(phi_M) * (y_IM_P_prime*C%deg2rad)

    ! See equation (2.13) or equation (B.20) in Reerink et al. (2010):
    a = COS(lambda_M) * COS(phi_M) * x_3D_P_prime  +  SIN(lambda_M) * COS(phi_M) * y_3D_P_prime  +  SIN(phi_M) * z_3D_P_prime

    ! See equation (2.12) or equation (B.19) in Reerink et al. (2010):
    t_P = (2._dp * C%earth_radius**2 + 2._dp * C%earth_radius * a) / (C%earth_radius**2 + 2._dp * C%earth_radius * a + x_3D_P_prime**2 + y_3D_P_prime**2 + z_3D_P_prime**2)

    ! See equations (2.9-2.11) or equations (B.16-B.18) in Reerink et al. (2010):
    x_3D_P =  C%earth_radius * COS(lambda_M) * COS(phi_M) * (t_P - 1._dp) + x_3D_P_prime * t_P
    y_3D_P =  C%earth_radius * SIN(lambda_M) * COS(phi_M) * (t_P - 1._dp) + y_3D_P_prime * t_P
    z_3D_P =  C%earth_radius *                   SIN(phi_M) * (t_P - 1._dp) + z_3D_P_prime * t_P

    ! See equation (2.7) or equation (B.24) in Reerink et al. (2010):
    IF(x_3D_P <  0._dp                      ) THEN
      lambda_P = 180._dp + C%rad2deg * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P >= 0._dp) THEN
      lambda_P =           C%rad2deg * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P <  0._dp) THEN
      lambda_P = 360._dp + C%rad2deg * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P >  0._dp) THEN
      lambda_P =  90._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P <  0._dp) THEN
      lambda_P = 270._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P == 0._dp) THEN
      lambda_P =   0._dp
    END IF

    ! See equation (2.8) or equation (B.25) in Reerink et al. (2010):
    IF(x_3D_P /= 0._dp .OR. y_3D_P /= 0._dp) THEN
      phi_P = C%rad2deg * ATAN(z_3D_P / sqrt(x_3D_P**2 + y_3D_P**2))
    ELSE IF(z_3D_P >  0._dp) THEN
      phi_P =   90._dp
    ELSE IF(z_3D_P <  0._dp) THEN
      phi_P =  -90._dp
    END IF
  END SUBROUTINE inverse_oblique_sg_projection
     
! == Some basic geometrical operations
  FUNCTION   IsInTriangle( pa, pb, pc, p) RESULT(isso)
    ! Check if the point p lies inside the triangle abc, or within distance tol_dist of its edges

    REAL(dp), DIMENSION(2),   INTENT(IN)          :: pa, pb, pc, p
    LOGICAL                                       :: isso
    REAL(dp)                                      :: as_x, as_y, s1, s2, s3  
    REAL(dp), PARAMETER                           :: tol = 1E-8_dp

    ! Check if p lies in the interior of the triangle abc
    as_x = p(1)-pa(1)
    as_y = p(2)-pa(2)

    s1 = ((pb(1)-pa(1))*as_y-(pb(2)-pa(2))*as_x)
    s2 = ((pc(1)-pa(1))*as_y-(pc(2)-pa(2))*as_x)
    s3 = ((pc(1)-pb(1))*(p(2)-pb(2))-(pc(2)-pb(2))*(p(1)-pb(1)))
   
    isso = .FALSE.   

    IF (s1 > -tol .AND. s2 < tol .AND. s3 > -tol) THEN
      isso = .TRUE.
      RETURN
    END IF
   
  END FUNCTION IsInTriangle  
  FUNCTION   IsBoundarySegment(mesh,v1,v2) RESULT(isso)
   ! Determine whether or not the line between two vertices is an Edge segment

   TYPE(type_mesh),          INTENT(IN)        :: mesh
   INTEGER,                  INTENT(IN)          :: v1, v2
   LOGICAL                                       :: isso

   IF (mesh%edge_index(v1)==0 .OR. mesh%edge_index(v2)==0) THEN
    isso = .FALSE.
    RETURN
   END IF

   isso = .FALSE.

   IF (mesh%edge_index(v1)==1) THEN
    IF (mesh%edge_index(v2)==8 .OR. &
        mesh%edge_index(v2)==1 .OR. &
        mesh%edge_index(v2)==2) THEN
      isso = .TRUE.
    END IF
   ELSEIF (mesh%edge_index(v1)==2) THEN
    IF (mesh%edge_index(v2)==8 .OR. &
       mesh%edge_index(v2)==1 .OR. &
       mesh%edge_index(v2)==2 .OR. &
       mesh%edge_index(v2)==3 .OR. &
       mesh%edge_index(v2)==4) THEN
      isso = .TRUE.
    END IF
   ELSEIF (mesh%edge_index(v1)==3) THEN
    IF (mesh%edge_index(v2)==2 .OR. &
       mesh%edge_index(v2)==3 .OR. &
       mesh%edge_index(v2)==4) THEN
      isso = .TRUE.
    END IF
   ELSEIF (mesh%edge_index(v1)==4) THEN
    IF (mesh%edge_index(v2)==2 .OR. &
       mesh%edge_index(v2)==3 .OR. &
       mesh%edge_index(v2)==4 .OR. &
       mesh%edge_index(v2)==5 .OR. &
       mesh%edge_index(v2)==6) THEN
      isso = .TRUE.
    END IF
   ELSEIF (mesh%edge_index(v1)==5) THEN
    IF (mesh%edge_index(v2)==4 .OR. &
       mesh%edge_index(v2)==5 .OR. &
       mesh%edge_index(v2)==6) THEN
      isso = .TRUE.
    END IF
   ELSEIF (mesh%edge_index(v1)==6) THEN
    IF (mesh%edge_index(v2)==4 .OR. &
       mesh%edge_index(v2)==5 .OR. &
       mesh%edge_index(v2)==6 .OR. &
       mesh%edge_index(v2)==7 .OR. &
       mesh%edge_index(v2)==8) THEN
      isso = .TRUE.
    END IF
   ELSEIF (mesh%edge_index(v1)==7) THEN
    IF (mesh%edge_index(v2)==6 .OR. &
       mesh%edge_index(v2)==7 .OR. &
       mesh%edge_index(v2)==8) THEN
      isso = .TRUE.
    END IF
   ELSEIF (mesh%edge_index(v1)==8) THEN
    IF (mesh%edge_index(v2)==6 .OR. &
       mesh%edge_index(v2)==7 .OR. &
       mesh%edge_index(v2)==8 .OR. &
       mesh%edge_index(v2)==1 .OR. &
       mesh%edge_index(v2)==2) THEN
      isso = .TRUE.
    END IF
   END IF

  END FUNCTION IsBoundarySegment
  SUBROUTINE IsEncroachedUpon(mesh,v1,v2,isso)
    ! TRUE if the line between v1 and v2 is encroached upon by any other
    ! vertex

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                    INTENT(IN)        :: v1,v2
    LOGICAL,                    INTENT(OUT)       :: isso
    INTEGER                                       :: v

    isso = .FALSE.
    DO v = 1, mesh%nV
      IF (v==v1 .OR. v==v2) CYCLE
      IF (NORM2(mesh%V(v,:) - (mesh%V(v1,:) + mesh%V(v2,:))/2._dp) < NORM2(mesh%V(v1,:) - mesh%V(v2,:))/2._dp) THEN
        isso = .TRUE.
        RETURN
      END IF
    END DO
  END SUBROUTINE IsEncroachedUpon
  SUBROUTINE EncroachesUpon(mesh,x,y,isso,va,vb)
    ! TRUE if the point [x,y] encroaches upon any boundary segment [va,vb]

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp),                   INTENT(IN)        :: x, y
    LOGICAL,                    INTENT(OUT)       :: isso
    INTEGER,                    INTENT(OUT)       :: va,vb

    INTEGER                                       :: ti, p1, p2, p3

    isso = .FALSE.
    va   = 0
    vb   = 0

    ! Go through all existing triangles - maybe only the triangle containing
    ! vi needs to be checked?
    DO ti = 1, mesh%nTri
      ! Only check boundary triangles
      IF (mesh%Tri_edge_index(ti)==0) CYCLE

      p1 = mesh%Tri(ti,1)
      p2 = mesh%Tri(ti,2)
      p3 = mesh%Tri(ti,3)

      ! Only check boundary segments
      IF (IsBoundarySegment(mesh,p1,p2)) THEN
        IF (NORM2([x,y] - (mesh%V(p1,:) + mesh%V(p2,:))/2._dp) < NORM2(mesh%V(p1,:) - mesh%V(p2,:))/2._dp) THEN
          va = p1
          vb = p2
          isso = .TRUE.
          RETURN
        END IF
      END IF

      IF (IsBoundarySegment(mesh,p2,p3)) THEN
        IF (NORM2([x,y] - (mesh%V(p2,:) + mesh%V(p3,:))/2._dp) < NORM2(mesh%V(p2,:) - mesh%V(p3,:))/2._dp) THEN
          va = p2
          vb = p3
          isso = .TRUE.
          RETURN
        END IF
      END IF

      IF (IsBoundarySegment(mesh,p3,p1)) THEN
        IF (NORM2([x,y] - (mesh%V(p3,:) + mesh%V(p1,:))/2._dp) < NORM2(mesh%V(p3,:) - mesh%V(p1,:))/2._dp) THEN
          va = p3
          vb = p1
          isso = .TRUE.
          RETURN
        END IF
      END IF
    END DO ! DO ti = 1, mesh%nTri
  END SUBROUTINE EncroachesUpon
  SUBROUTINE UpdateTriangleCircumcenter(mesh,ti)
    ! Calculate the circumcenter of mesh triangle ti
    
    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti
    
    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: v1, v2, v3, cc
    
    v1 = mesh%V( mesh%Tri( ti,1),:)
    v2 = mesh%V( mesh%Tri( ti,2),:)
    v3 = mesh%V( mesh%Tri( ti,3),:)     
    CALL FindCircumcenter(v1, v2, v3, cc)
    
    ! If FindCircumcenter yields infinity, it's because p and q have the
    ! same y-coordinate. Rearrange vertices in triangle matrix (maintaining
    ! counter-clockwise orientation) and try again.
    
    IF (cc(1) > (mesh%xmax-mesh%xmin)*1000._dp .OR. cc(2) > (mesh%ymax-mesh%ymin)*1000._dp) THEN
      mesh%Tri(ti,:) = [mesh%Tri(ti,2), mesh%Tri(ti,3), mesh%Tri(ti,1)]
      v1 = mesh%V( mesh%Tri( ti,1),:)
      v2 = mesh%V( mesh%Tri( ti,2),:)
      v3 = mesh%V( mesh%Tri( ti,3),:) 
      CALL FindCircumcenter(v1, v2, v3, cc)
    END IF
    
    IF (cc(1) > (mesh%xmax-mesh%xmin)*1000._dp .OR. cc(2) > (mesh%ymax-mesh%ymin)*1000._dp) THEN
      mesh%Tri(ti,:) = [mesh%Tri(ti,2), mesh%Tri(ti,3), mesh%Tri(ti,1)]
      v1 = mesh%V( mesh%Tri( ti,1),:)
      v2 = mesh%V( mesh%Tri( ti,2),:)
      v3 = mesh%V( mesh%Tri( ti,3),:) 
      CALL FindCircumcenter(v1, v2, v3, cc)
    END IF
    
    IF (cc(1) > (mesh%xmax-mesh%xmin)*1000._dp .OR. cc(2) > (mesh%ymax-mesh%ymin)*1000._dp) THEN
      WRITE(0,*) '  UpdateTriangleCircumcenter - ERROR: triangle  doesn''t yield a valid circumcenter!'
    END IF
    
    mesh%Tricc(ti,:) = cc
    
  END SUBROUTINE UpdateTriangleCircumcenter
  
! == Routines for merging meshes created by parallel processes
  SUBROUTINE MergeVertices( mesh, nVl, nVr, nTril, nTrir, T, nT, vil, vir, orientation)
    ! Merge overlapping vertices vil and vir
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER, DIMENSION(:,:),    INTENT(INOUT)     :: T
    INTEGER,                    INTENT(IN)        :: nVl, nVr, nTril, nTrir, nT
    INTEGER,                    INTENT(IN)        :: vil, vir
    INTEGER,                    INTENT(IN)        :: orientation
    
    ! Local variables
    LOGICAL                                       :: IsNorth, IsSouth
    INTEGER                                       :: tieu, tiel, tiwl, tiwu, vilu, vill, viru, virl
    INTEGER                                       :: n, nnext, nprev
    INTEGER                                       :: nCl, nCr, niTril, niTrir, vi, ci, vc, ci2, vcr, tii, ti, vi2

    IF ( ABS( mesh%V( vil,1) - mesh%V( vir,1)) > mesh%tol_dist .OR. ABS( mesh%V( vil,2) - mesh%V( vir,2)) > mesh%tol_dist) THEN
      WRITE(0,'(A,I5,A,I5)') 'Cannot merge vertices ', vil, ' and ', vir, ' - no overlap!'
      STOP
    END IF

    IsNorth = .FALSE.
    IsSouth = .FALSE.

    IF (orientation == 0) THEN
    
      IF   ( mesh%edge_index( vil) == 2) THEN
        IF ( mesh%edge_index( vir) == 8) THEN
          IsNorth = .TRUE.
        ELSE
          WRITE(0,*) 'Cannot merge north-east corner of mesh_left with non-north-west-corner of mesh_right!'
          STOP
        END IF
      END IF
      IF   ( mesh%edge_index( vil) == 4) THEN
        IF ( mesh%edge_index( vir) == 6) THEN
          IsSouth = .TRUE.
        ELSE
          WRITE(0,*) 'Cannot merge south-east corner of mesh_left with non-south-west-corner of mesh_right!'
          STOP
        END IF
      END IF
    
    ELSEIF (orientation == 1) THEN
    
      IF   ( mesh%edge_index( vil) == 8) THEN
        IF ( mesh%edge_index( vir) == 6) THEN
          IsNorth = .TRUE.
        ELSE
          WRITE(0,*) 'Cannot merge north-east corner of mesh_left with non-north-west-corner of mesh_right!'
          STOP
        END IF
      END IF
      IF   ( mesh%edge_index( vil) == 2) THEN
        IF ( mesh%edge_index( vir) == 4) THEN
          IsSouth = .TRUE.
        ELSE
          WRITE(0,*) 'Cannot merge south-east corner of mesh_left with non-south-west-corner of mesh_right!'
          STOP
        END IF
      END IF
      
    END IF

    ! == Triangle connectivity 
    ! make sure newly adjacent triangles are connected
    tieu = mesh%iTri( vil, 1)
    tiel = mesh%iTri( vil, mesh%niTri( vil))
    tiwl = mesh%iTri( vir, 1)
    tiwu = mesh%iTri( vir, mesh%niTri( vir))
    vilu = mesh%C(    vil, 1)
    vill = mesh%C(    vil, mesh%nC( vil))
    virl = mesh%C(    vir, 1)
    viru = mesh%C(    vir, mesh%nC( vir))
    IF     ( IsNorth) THEN
      DO n = 1, 3
        nnext = n+1
        nprev = n-1
        IF (nnext == 4) nnext = 1
        IF (nprev == 0) nprev = 3

        IF ( mesh%Tri( tiel,n) == vil .AND. mesh%Tri( tiel,nprev) == vill) mesh%TriC( tiel,nnext) = tiwl
        IF ( mesh%Tri( tiwl,n) == vir .AND. mesh%Tri( tiwl,nnext) == virl) mesh%TriC( tiwl,nprev) = tiel
      END DO
    ELSEIF ( IsSouth) THEN
      DO n = 1, 3
        nnext = n+1
        nprev = n-1
        IF (nnext == 4) nnext = 1
        IF (nprev == 0) nprev = 3

        IF ( mesh%Tri( tieu,n) == vil .AND. mesh%Tri( tieu,nnext) == vilu) mesh%TriC( tieu,nprev) = tiwu
        IF ( mesh%Tri( tiwu,n) == vir .AND. mesh%Tri( tiwu,nprev) == viru) mesh%TriC( tiwu,nnext) = tieu
      END DO
    ELSE
      DO n = 1, 3
        nnext = n+1
        nprev = n-1
        IF ( nnext == 4) nnext = 1
        IF ( nprev == 0) nprev = 3

        IF ( mesh%Tri( tieu,n) == vil .AND. mesh%Tri( tieu,nnext) == vilu) mesh%TriC( tieu,nprev) = tiwu
        IF ( mesh%Tri( tiel,n) == vil .AND. mesh%Tri( tiel,nprev) == vill) mesh%TriC( tiel,nnext) = tiwl
        IF ( mesh%Tri( tiwl,n) == vir .AND. mesh%Tri( tiwl,nnext) == virl) mesh%TriC( tiwl,nprev) = tiel
        IF ( mesh%Tri( tiwu,n) == vir .AND. mesh%Tri( tiwu,nprev) == viru) mesh%TriC( tiwu,nnext) = tieu
      END DO
    END IF

    ! == Vertex connectivity - add connections from vir to vil
    nCl = mesh%nC( vil)
    nCr = mesh%nC( vir)
    IF      ( IsNorth) THEN
      mesh%nC( vil) = nCl + nCr - 1
      mesh%C(  vil, nCl+1:nCl+nCr-1) = mesh%C( vir, 2:nCr)
    ELSEIF  (IsSouth) THEN
      mesh%nC( vil) = nCl + nCr - 1
      mesh%C(  vil, nCr:nCl+nCr-1) = mesh%C( vil, 1:nCl)
      mesh%C(  vil,   1:    nCr-1) = mesh%C( vir, 1:nCr-1)
    ELSE
      mesh%nC( vil) = nCl + nCr - 2
      mesh%C(  vil, nCl+1:nCl+nCr-2) = mesh%C( vir,2:nCr-1)
    END IF

    ! == Reverse vertex connectivity 
    ! Make sure neighbours of vir now connect to vil
    DO ci = 1, mesh%nC( vir)
      vc = mesh%C( vir,ci)
      DO ci2 = 1, mesh%nC( vc)
        vcr = mesh%C( vc,ci2)
        IF (vcr == vir) THEN
          mesh%C( vc,ci2) = vil
          EXIT
        END IF
      END DO
    END DO

    ! == Inverse triangles
    ! Add iTri from vir to vil
    niTril = mesh%niTri( vil)
    niTrir = mesh%niTri( vir)
    mesh%niTri( vil) = niTril + niTrir
    IF     (IsNorth) THEN
      mesh%iTri( vil, niTril+1:niTril+niTrir) = mesh%iTri( vir, 1:niTrir)
    ELSEIF (IsSouth) THEN
      mesh%iTri( vil, niTrir+1:niTril+niTrir) = mesh%iTri( vil, 1:niTril)
      mesh%iTri( vil,        1:       niTrir) = mesh%iTri( vir, 1:niTrir)
    ELSE
      mesh%iTri( vil, niTril+1:niTril+niTrir) = mesh%iTri( vir, 1:niTrir)
    END IF

    ! == Triangles 
    ! Make sure triangles containing vir now contain vil
    DO tii = 1, mesh%niTri( vil)
      ti = mesh%iTri( vil,tii)
      DO n = 1, 3
        IF (mesh%Tri( ti,n) == vir) THEN
          mesh%Tri( ti,n) = vil
          EXIT
        END IF
      END DO
    END DO

    ! == Edge index
    IF (orientation == 0) THEN
    
      ! 1 for north, 5 for south, 0 elsewhere
      IF     (IsNorth) THEN
        mesh%edge_index( vil) = 1
      ELSEIF (IsSouth) THEN
        mesh%edge_index( vil) = 5
      ELSE
        mesh%edge_index( vil) = 0
      END IF
    
    ELSEIF (orientation == 1) THEN
    
      ! 7 for north, 3 for south, 0 elsewhere
      IF     (IsNorth) THEN
        mesh%edge_index( vil) = 7
      ELSEIF (IsSouth) THEN
        mesh%edge_index( vil) = 3
      ELSE
        mesh%edge_index( vil) = 0
      END IF
    
    END IF
    
    ! == Remove vir from the mesh
    mesh%V(          vir:mesh%nV-1 ,:) = mesh%V(          vir+1:mesh%nV ,:)
    mesh%nC(         vir:mesh%nV-1   ) = mesh%nC(         vir+1:mesh%nV   )
    mesh%C(          vir:mesh%nV-1 ,:) = mesh%C(          vir+1:mesh%nV ,:)
    mesh%niTri(      vir:mesh%nV-1   ) = mesh%niTri(      vir+1:mesh%nV   )
    mesh%iTri(       vir:mesh%nV-1 ,:) = mesh%iTri(       vir+1:mesh%nV ,:)
    mesh%edge_index( vir:mesh%nV-1   ) = mesh%edge_index( vir+1:mesh%nV   )

    mesh%V(          mesh%nV,:) = [0, 0]
    mesh%nC(         mesh%nV  ) = 0
    mesh%niTri(      mesh%nV  ) = 0
    mesh%edge_index( mesh%nV  ) = 0

    mesh%nV = mesh%nV - 1

    ! == Update all vertex connections and triangles pointing to vertices higher than vir
    ! Since those vertex indices have just decreased by one
    DO vi2 = nVl+1, nVl+nVr
      DO ci = 1, mesh%nC( vi2)
        vc = mesh%C( vi2,ci)
        IF (vc > vir)  mesh%C( vi2,ci) = vc-1
      END DO
    END DO

    DO ti = 1, nT
      vi = T(ti,1)
      DO ci = 1, mesh%nC( vi)
        vc = mesh%C( vi,ci)
        IF (vc > vir) mesh%C( vi,ci) = vc-1    
      END DO
    END DO

    DO ti = nTril, nTril+nTrir
      DO n = 1, 3
        vi2 = mesh%Tri( ti,n)
        IF (vi2 > vir) mesh%Tri( ti,n) = vi2-1
      END DO
    END DO

    !  Update the boundary translation table
    DO ti = 1, nT
      IF ( T(ti,2) > vir) T(ti,2) = T(ti,2) - 1
    END DO
    
  END SUBROUTINE MergeVertices
  SUBROUTINE SwitchVertices( mesh, vi1, vi2)
    ! Switch vertices vi1 and vi2
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: vi1, vi2
    INTEGER, DIMENSION(:,:), ALLOCATABLE          :: ctovi1, ctovi2
    INTEGER                                       :: ci, vc, ti, iti, n, ci2
    
    ALLOCATE( ctovi1( mesh%nconmax, 2))
    ALLOCATE( ctovi2( mesh%nconmax, 2))
    
    ! == Vertex coordinates
    mesh%V(vi1,:) = mesh%V(vi1,:) + mesh%V(vi2,:)
    mesh%V(vi2,:) = mesh%V(vi1,:) - mesh%V(vi2,:)
    mesh%V(vi1,:) = mesh%V(vi1,:) - mesh%V(vi2,:)

    ! == Vertex connectivity
    ! Find connections that point to vi1 and vi2
    ctovi1 = 0
    ctovi2 = 0
    DO ci = 1, mesh%nC(vi1)
      vc = mesh%C(vi1,ci)
      IF (vc==vi2) CYCLE
      DO ci2 = 1, mesh%nC(vc)
        IF (mesh%C(vc,ci2)==vi1) THEN
          ctovi1(ci,:) = [vc,ci2]
          EXIT
        END IF
      END DO
    END DO
    DO ci = 1, mesh%nC(vi2)
      vc = mesh%C(vi2,ci)
      IF (vc==vi1) CYCLE
      DO ci2 = 1, mesh%nC(vc)
        IF (mesh%C(vc,ci2)==vi2) THEN
          ctovi2(ci,:) = [vc,ci2]
          EXIT
        END IF
      END DO
    END DO

    ! Switch their own connectivity lists
    mesh%C(  vi1,:) = mesh%C(  vi1,:) + mesh%C(  vi2,:)
    mesh%C(  vi2,:) = mesh%C(  vi1,:) - mesh%C(  vi2,:)
    mesh%C(  vi1,:) = mesh%C(  vi1,:) - mesh%C(  vi2,:)
    mesh%nC( vi1  ) = mesh%nC( vi1  ) + mesh%nC( vi2  )
    mesh%nC( vi2  ) = mesh%nC( vi1  ) - mesh%nC( vi2  )
    mesh%nC( vi1  ) = mesh%nC( vi1  ) - mesh%nC( vi2  )

    ! If they are interconnected, change those too
    DO ci = 1, mesh%nC(vi1)
      IF (mesh%C(vi1,ci)==vi1) mesh%C(vi1,ci) = vi2
    END DO
    DO ci = 1, mesh%nC(vi2)
      IF (mesh%C(vi2,ci)==vi2) mesh%C(vi2,ci) = vi1
    END DO

    ! Switch the other connections
    DO ci = 1, mesh%nC(vi2)
      IF (ctovi1(ci,1)==0) CYCLE
      mesh%C( ctovi1(ci,1), ctovi1(ci,2)) = vi2
    END DO
    DO ci = 1, mesh%nC(vi1)
      IF (ctovi2(ci,1)==0) CYCLE
      mesh%C( ctovi2(ci,1), ctovi2(ci,2)) = vi1
    END DO

    ! == Triangles
    ! First update the surrounding triangles
    DO iti = 1, mesh%niTri(vi1)
      ti = mesh%iTri(vi1,iti)
      DO n = 1, 3
        IF (mesh%Tri(ti,n)==vi1) mesh%Tri(ti,n) = -1
      END DO
    END DO
    DO iti = 1, mesh%niTri(vi2)
      ti = mesh%iTri(vi2,iti)
      DO n = 1, 3
        IF (mesh%Tri(ti,n)==vi2) mesh%Tri(ti,n) = -2
      END DO
    END DO
    DO iti = 1, mesh%niTri(vi1)
      ti = mesh%iTri(vi1,iti)
      DO n = 1, 3
        IF (mesh%Tri(ti,n)== -1) mesh%Tri(ti,n) = vi2
      END DO
    END DO
    DO iti = 1, mesh%niTri(vi2)
      ti = mesh%iTri(vi2,iti)
      DO n = 1, 3
        IF (mesh%Tri(ti,n)== -2) mesh%Tri(ti,n) = vi1
      END DO
    END DO

    ! Then switch their own triangles
    mesh%iTri(  vi1,:) = mesh%iTri(  vi1,:) + mesh%iTri(  vi2,:)
    mesh%iTri(  vi2,:) = mesh%iTri(  vi1,:) - mesh%iTri(  vi2,:)
    mesh%iTri(  vi1,:) = mesh%iTri(  vi1,:) - mesh%iTri(  vi2,:)
    mesh%niTri( vi1  ) = mesh%niTri( vi1  ) + mesh%niTri( vi2  )
    mesh%niTri( vi2  ) = mesh%niTri( vi1  ) - mesh%niTri( vi2  )
    mesh%niTri( vi1  ) = mesh%niTri( vi1  ) - mesh%niTri( vi2  )

    ! == Edge indices
    mesh%edge_index( vi1  ) = mesh%edge_index( vi1  ) + mesh%edge_index( vi2  )
    mesh%edge_index( vi2  ) = mesh%edge_index( vi1  ) - mesh%edge_index( vi2  )
    mesh%edge_index( vi1  ) = mesh%edge_index( vi1  ) - mesh%edge_index( vi2  )
    
    DEALLOCATE( ctovi1)
    DEALLOCATE( ctovi2)
    
  END SUBROUTINE SwitchVertices
  SUBROUTINE RedoTriEdgeIndices( mesh)
    ! Called by MergeSubmeshes, ran by a single processor for that processor's submesh
    ! Be careful not to include any sync statements here!
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: vi, vi_prev, ti, vi_NW, vi_NE, vi_SE, vi_SW
        
    mesh%Tri_edge_index = 0
    
    ! West (SN)
    ! =========
    vi_prev = 1
    vi = mesh%C( vi_prev, mesh%nC( vi_prev))
    DO WHILE (ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist)
      vi = mesh%C( vi_prev, mesh%nC( vi_prev))
      ti = mesh%iTri( vi, 1)
      mesh%Tri_edge_index(ti) = 7
      ti = mesh%iTri( vi, mesh%niTri(vi))
      mesh%Tri_edge_index(ti) = 7
      vi_prev = vi
    END DO
    vi_NW = vi
    
    ! North (WE)
    ! ==========
    DO WHILE (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist)
      vi = mesh%C( vi_prev, mesh%nC( vi_prev))
      ti = mesh%iTri( vi, 1)
      mesh%Tri_edge_index(ti) = 1
      ti = mesh%iTri( vi, mesh%niTri(vi))
      mesh%Tri_edge_index(ti) = 1
      vi_prev = vi
    END DO
    vi_NE = vi
    
    ! East (NS)
    ! =========
    DO WHILE (ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist)
      vi = mesh%C( vi_prev, mesh%nC( vi_prev))
      ti = mesh%iTri( vi, 1)
      mesh%Tri_edge_index(ti) = 3
      ti = mesh%iTri( vi, mesh%niTri(vi))
      mesh%Tri_edge_index(ti) = 3
      vi_prev = vi
    END DO
    vi_SE = vi
    
    ! South (EW)
    ! ==========
    DO WHILE (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist)
      vi = mesh%C( vi_prev, mesh%nC( vi_prev))
      ti = mesh%iTri( vi, 1)
      mesh%Tri_edge_index(ti) = 5
      ti = mesh%iTri( vi, mesh%niTri(vi))
      mesh%Tri_edge_index(ti) = 5
      vi_prev = vi
    END DO
    vi_SW = vi
    
    ! Correct the last ones on each side
    ti = mesh%iTri( vi_SW, mesh%niTri(vi_SW))
    mesh%Tri_edge_index(ti) = 7
    ti = mesh%iTri( vi_SE, mesh%niTri(vi_SE))
    mesh%Tri_edge_index(ti) = 5
    ti = mesh%iTri( vi_NE, mesh%niTri(vi_NE))
    mesh%Tri_edge_index(ti) = 3
    ti = mesh%iTri( vi_NW, mesh%niTri(vi_NW))
    mesh%Tri_edge_index(ti) = 1
        
  END SUBROUTINE RedoTriEdgeIndices
  SUBROUTINE DetermineProcessXYDomain_regular(region, i, orientation, xmin, xmax, ymin, ymax)
    ! Generate symmetrical domains (used for generating the first mesh, when no previous
    ! mesh exists for determining workload-balanced domains)
  
    ! In/output variables
    TYPE(type_model_region),    INTENT(IN)        :: region
    INTEGER,                    INTENT(IN)        :: i
    INTEGER,                    INTENT(IN)        :: orientation
    REAL(dp),                   INTENT(OUT)       :: xmin, xmax, ymin, ymax
    
    ! Local variables:
    REAL(dp)                                      :: xminmin, yminmin, xmaxmax, ymaxmax, xrange, yrange
        
    xminmin = MINVAL(region%init%x)
    xmaxmax = MAXVAL(region%init%x)
    yminmin = MINVAL(region%init%y)
    ymaxmax = MAXVAL(region%init%y)
    
!    ! Exception for Antarctica. Since the domain in the initial file is perfectly square,
!    ! this gives problems with mesh generation (especially remapping, since it leads to a
!    ! lof of cocircular vertex groups, which give ill-defined Voronoi cells)
!    IF (region%name == 'ANT') THEN
!      yminmin = yminmin + 10000._dp
!      ymaxmax = ymaxmax - 10000._dp      
!    END IF
    
    xrange  = xmaxmax - xminmin
    yrange  = ymaxmax - yminmin
    
    IF (orientation == 0) THEN
      ! North-south domain borders, used everywhere except for GRL
      xmin = xminmin + REAL(i  ) * xrange / REAL(par%n)
      xmax = xminmin + REAL(i+1) * xrange / REAL(par%n)
      ymin = yminmin
      ymax = ymaxmax
    ELSE
      ! East-west domain borders, used only for GRL
      xmin = xminmin
      xmax = xmaxmax
      ymin = yminmin + REAL(i  ) * yrange / REAL(par%n)
      ymax = yminmin + REAL(i+1) * yrange / REAL(par%n)
    END IF
        
    
  END SUBROUTINE DetermineProcessXYDomain_regular
  SUBROUTINE DetermineProcessXYDomain_balanced(mesh, i, orientation, xmin, xmax, ymin, ymax)
    ! Generate workload-balanced domains, based on vertices of the previous mesh
  
    ! In/output variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                    INTENT(IN)        :: i
    INTEGER,                    INTENT(IN)        :: orientation
    REAL(dp),                   INTENT(OUT)       :: xmin, xmax, ymin, ymax
    
    ! Local variables:
    REAL(dp)                                      :: xminmin, yminmin, xmaxmax, ymaxmax, xrange, yrange
        
    xminmin = mesh%xmin
    xmaxmax = mesh%xmax
    yminmin = mesh%ymin
    ymaxmax = mesh%ymax
    
    xrange  = xmaxmax - xminmin
    yrange  = ymaxmax - yminmin
    
    IF (orientation == 0) THEN
      ! North-south domain borders, used everywhere except for GRL
      ymin = yminmin
      ymax = ymaxmax
      CALL DivideDomain_x(mesh, ymin, ymax, par%n, i, xmin, xmax)
    ELSE
      ! East-west domain borders, used only for GRL
      xmin = xminmin
      xmax = xmaxmax
      CALL DivideDomain_y(mesh, xmin, xmax, par%n, i, ymin, ymax)
    END IF
    
  END SUBROUTINE DetermineProcessXYDomain_balanced
  SUBROUTINE DivideDomain_x(mesh, ymin, ymax, n, i, xmin, xmax)
    ! Given a mesh, divide the region between ymin and ymax into n parts with
    ! equal numbers of vertices. Return the [xmin,xmax] limits of the i-th part.
  
    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp),                   INTENT(IN)        :: ymin, ymax
    INTEGER,                    INTENT(IN)        :: n, i
    REAL(dp),                   INTENT(OUT)       :: xmin, xmax
    
    ! Local variables:
    INTEGER,  PARAMETER                           :: nbins = 1000
    INTEGER                                       :: nV_sub
    INTEGER,  DIMENSION(:  ), ALLOCATABLE         :: vbins, intvbins
    REAL(dp), DIMENSION(:  ), ALLOCATABLE         :: xvals, xvals_opt
    INTEGER                                       :: bi, vi, p
    REAL(dp)                                      :: xrange, sumv
      
    ! Determine x ranges of bins
    xrange = mesh%xmax - mesh%xmin
    ALLOCATE( xvals( nbins))
    DO bi = 1, nbins
      xvals(bi) = mesh%xmin + (xrange * REAL(bi-1)/REAL(nbins-1))
    END DO
    
    ! Fill histogram
    ALLOCATE( vbins( nbins))
    vbins = 0
    vbins = CEILING( REAL(mesh%nV) * 0.01_dp/ REAL(par%n))
    DO vi = 1, mesh%nV
      IF (mesh%V(vi,2) >= ymin .AND. mesh%V(vi,2) <= ymax) THEN
        bi = 1 + NINT( REAL(nbins-1) * (mesh%V(vi,1) - mesh%xmin) / xrange)
        vbins(bi) = vbins(bi) + 1
      END IF
    END DO
    
    nV_sub = SUM(vbins)
    
    ! Integrate histogram
    ALLOCATE( intvbins( nbins))
    intvbins = 0
    intvbins(1) = vbins(1)
    DO bi = 2, nbins
      intvbins(bi) = intvbins(bi-1) + vbins(bi)
    END DO
    
    ! Determine domain range
    ALLOCATE( xvals_opt( n+1))
    xvals_opt(1)   = mesh%xmin
    xvals_opt(n+1) = mesh%xmax
    
    bi = 1
    DO p = 1, n-1
      sumv = REAL(p * nV_sub) / REAL(n)
      DO WHILE ( REAL(intvbins(bi)) < sumv .AND. bi < nbins)
        bi = bi + 1
      END DO
      xvals_opt(p+1) = xvals(bi)
    END DO
    
    xmin = xvals_opt(i+1)
    xmax = xvals_opt(i+2)
    
    DEALLOCATE( vbins)
    DEALLOCATE( xvals)
    DEALLOCATE( intvbins)
    DEALLOCATE( xvals_opt)
    
  END SUBROUTINE DivideDomain_x
  SUBROUTINE DivideDomain_y(mesh, xmin, xmax, n, i, ymin, ymax)
    ! Given a mesh, divide the region between xmin and xmax into n parts with
    ! equal numbers of vertices. Return the [ymin,ymax] limits of the i-th part.
  
    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp),                   INTENT(IN)        :: xmin, xmax
    INTEGER,                    INTENT(IN)        :: n, i
    REAL(dp),                   INTENT(OUT)       :: ymin, ymax
    
    ! Local variables:
    INTEGER,  PARAMETER                           :: nbins = 1000
    INTEGER                                       :: nV_sub
    INTEGER,  DIMENSION(:  ), ALLOCATABLE         :: vbins, intvbins
    REAL(dp), DIMENSION(:  ), ALLOCATABLE         :: yvals, yvals_opt
    INTEGER                                       :: bi, vi, p
    REAL(dp)                                      :: yrange, sumv
      
    ! Determine y ranges of bins
    yrange = mesh%ymax - mesh%ymin
    ALLOCATE( yvals( nbins))
    DO bi = 1, nbins
      yvals(bi) = mesh%ymin + (yrange * REAL(bi-1)/REAL(nbins-1))
    END DO
    
    ! Fill histogram
    ALLOCATE( vbins( nbins))
    vbins = 0
    vbins = CEILING( REAL(mesh%nV) * 0.01_dp/ REAL(par%n))
    DO vi = 1, mesh%nV
      IF (mesh%V(vi,1) >= xmin .AND. mesh%V(vi,1) <= xmax) THEN
        bi = 1 + NINT( REAL(nbins-1) * (mesh%V(vi,2) - mesh%ymin) / yrange)
        vbins(bi) = vbins(bi) + 1
      END IF
    END DO
    
    nV_sub = SUM(vbins)
    
    ! Integrate histogram
    ALLOCATE( intvbins( nbins))
    intvbins = 0
    intvbins(1) = vbins(1)
    DO bi = 2, nbins
      intvbins(bi) = intvbins(bi-1) + vbins(bi)
    END DO
    
    ! Determine domain range
    ALLOCATE( yvals_opt( n+1))
    yvals_opt(1)   = mesh%ymin
    yvals_opt(n+1) = mesh%ymax
    
    bi = 1
    DO p = 1, n-1
      sumv = REAL(p * nV_sub) / REAL(n)
      DO WHILE ( REAL(intvbins(bi)) < sumv .AND. bi < nbins)
        bi = bi + 1
      END DO
      yvals_opt(p+1) = yvals(bi)
    END DO
    
    ymin = yvals_opt(i+1)
    ymax = yvals_opt(i+2)
    
    DEALLOCATE( vbins)
    DEALLOCATE( yvals)
    DEALLOCATE( intvbins)
    DEALLOCATE( yvals_opt)
    
  END SUBROUTINE DivideDomain_y
  
! == Subroutines for creating Points of Interest (POI)
  SUBROUTINE FindPOIXYCoordinates( mesh)
    ! Use the stereographic projection parameters for the relevant model region
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: POIi
    REAL(dp)                                      :: lambda_M, phi_M, alpha_stereographic, lat, lon, x, y
    
    IF (mesh%region_name == 'NAM') THEN  
      
      lambda_M               = C%lambda_M_NAM
      phi_M                  = C%phi_M_NAM
      alpha_stereographic    = C%alpha_NAM
      
      DO POIi = 1, mesh%nPOI
        mesh%POI_coordinates(POIi,1) = C%mesh_POI_NAM_coordinates((POIi*2)-1)
        mesh%POI_coordinates(POIi,2) = C%mesh_POI_NAM_coordinates( POIi*2   )
        mesh%POI_resolutions(POIi  ) = C%mesh_POI_NAM_resolutions( POIi     )
      END DO
      
    ELSEIF (mesh%region_name == 'EAS') THEN
    
      lambda_M               = C%lambda_M_EAS
      phi_M                  = C%phi_M_EAS
      alpha_stereographic    = C%alpha_EAS
      
      DO POIi = 1, mesh%nPOI
        mesh%POI_coordinates(POIi,1) = C%mesh_POI_EAS_coordinates((POIi*2)-1)
        mesh%POI_coordinates(POIi,2) = C%mesh_POI_EAS_coordinates( POIi*2   )
        mesh%POI_resolutions(POIi  ) = C%mesh_POI_EAS_resolutions( POIi     )
      END DO
      
    ELSEIF (mesh%region_name == 'GRL') THEN
    
      lambda_M               = C%lambda_M_GRL
      phi_M                  = C%phi_M_GRL
      alpha_stereographic    = C%alpha_GRL
      
      DO POIi = 1, mesh%nPOI
        mesh%POI_coordinates(POIi,1) = C%mesh_POI_GRL_coordinates((POIi*2)-1)
        mesh%POI_coordinates(POIi,2) = C%mesh_POI_GRL_coordinates( POIi*2   )
        mesh%POI_resolutions(POIi  ) = C%mesh_POI_GRL_resolutions( POIi     )
      END DO
      
    ELSEIF (mesh%region_name == 'ANT') THEN
    
      lambda_M               = C%lambda_M_ANT
      phi_M                  = C%phi_M_ANT
      alpha_stereographic    = C%alpha_ANT
      
      DO POIi = 1, mesh%nPOI
        mesh%POI_coordinates(POIi,1) = C%mesh_POI_ANT_coordinates((POIi*2)-1)
        mesh%POI_coordinates(POIi,2) = C%mesh_POI_ANT_coordinates( POIi*2   )
        mesh%POI_resolutions(POIi  ) = C%mesh_POI_ANT_resolutions( POIi     )
      END DO
      
    END IF
    
    DO POIi = 1, mesh%nPOI
    
      lat = mesh%POI_coordinates(POIi,1)
      lon = mesh%POI_coordinates(POIi,2)
      
      CALL oblique_sg_projection(lon, lat, lambda_M, phi_M, alpha_stereographic, x, y)
    
      mesh%POI_XY_coordinates(POIi,1) = x
      mesh%POI_XY_coordinates(POIi,2) = y
    
    END DO  
  
  END SUBROUTINE FIndPOIXYCoordinates
  SUBROUTINE FindPOIVerticesAndWeights( mesh)
  ! For each POI in this region, find the indices of the three vertices spanning the triangle
  ! containing it, and find their respective weights for a trilinear interpolation.
  ! Called by all, executed by the master.
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: POIi, ti, vi_a, vi_b, vi_c
    REAL(dp), DIMENSION(2)                        :: p, pa, pb, pc
    REAL(dp)                                      :: Atot, Aa, Ab, Ac
  
    IF (par%master) THEN
    
    ti = 1
    DO POIi = 1, mesh%nPOI
    
      ! The location of the POI
      p = mesh%POI_XY_coordinates(POIi,:)
      
      ! Find the triangle containing it
      CALL FindContainingTriangle( mesh, p, ti)
      
      ! The three vertices spanning this triangle
      vi_a = mesh%Tri(ti,1)
      vi_b = mesh%Tri(ti,2)
      vi_c = mesh%Tri(ti,3)
      pa   = mesh%V(vi_a,:)
      pb   = mesh%V(vi_b,:)
      pc   = mesh%V(vi_c,:)
      
      ! Calculate their interpolation weights
      Atot = NORM2( CROSS( [pb(1)-pa(1), pb(2)-pa(2), 0._dp], [pc(1)-pa(1), pc(2)-pa(2), 0._dp] )) / 2._dp
      Aa   = NORM2( CROSS( [pb(1)-p( 1), pb(2)-p( 2), 0._dp], [pc(1)-p( 1), pc(2)-p( 2), 0._dp] )) / 2._dp
      Ab   = NORM2( CROSS( [pc(1)-p( 1), pc(2)-p( 2), 0._dp], [pa(1)-p( 1), pa(2)-p( 2), 0._dp] )) / 2._dp
      Ac   = NORM2( CROSS( [pa(1)-p( 1), pa(2)-p( 2), 0._dp], [pb(1)-p( 1), pb(2)-p( 2), 0._dp] )) / 2._dp
      
      ! Save the indices and weights
      mesh%POI_vi( POIi,:) = [vi_a,    vi_b,    vi_c   ]
      MESH%POI_w(  POIi,:) = [Aa/Atot, Ab/Atot, Ac/Atot]
      
    END DO ! DO POIi = 1, mesh%nPOI
    
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE FindPOIVerticesAndWeights
  
! == Some basic search operations on a mesh
  SUBROUTINE FindContainingTriangle( mesh, p, t)
   ! Start at initial guess t, search outward from there using a
   ! flood-fill algorithm.

    TYPE(type_mesh),          INTENT(INOUT)       :: mesh
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: t

    REAL(dp), DIMENSION(2)                        :: q, r, s
    INTEGER                                       :: ncycle, t_prev
    REAL(dp), DIMENSION(2)                        :: gcti, gctc
    REAL(dp)                                      :: d, dc, dcmin
    INTEGER                                       :: tc, tcmin
    LOGICAL                                       :: FoundIt
    INTEGER                                       :: n, ti, n2, tin
    
    ! If p lies outside the mesh domain, throw an error
    IF (p(1) < mesh%xmin .OR. p(1) > mesh%xmax .OR. p(2) < mesh%ymin .OR. p(2) > mesh%ymax) THEN
      WRITE(0,*) 'FindContainingTriangle - ERROR: point lies outside mesh domain!'
      STOP
    END IF

    ! See if the initial guess is correct.
    q = mesh%V(mesh%Tri(t,1),:)
    r = mesh%V(mesh%Tri(t,2),:)
    s = mesh%V(mesh%Tri(t,3),:)
    IF (IsInTriangle(q, r, s, p)) RETURN

    ! If not, start with a linear search.
    ncycle = 0
    t_prev = t
    DO WHILE (ncycle < mesh%nTri)

      gcti = (mesh%V(mesh%Tri(t,1),:) + mesh%V(mesh%Tri(t,2),:) + mesh%V(mesh%Tri(t,3),:)) / 3._dp
      d = NORM2(gcti - p)

      dcmin = d + 10._dp
      tcmin = 0
      DO n = 1, 3
        tc   = mesh%TriC(t,n)
        IF (tc==0)      CYCLE ! This triangle neighbour doesn't exist
        IF (tc==t_prev) CYCLE ! This triangle neighbour is the one we just came from
        gctc = (mesh%V(mesh%Tri(tc,1),:) + mesh%V(mesh%Tri(tc,2),:) + mesh%V(mesh%Tri(tc,3),:)) / 3._dp
        dc = NORM2( gctc - p)
        IF (dc < dcmin) THEN
          dcmin = dc
          tcmin = tc
        END IF
      END DO

      IF (dcmin < d) THEN
        t_prev = t
        t = tcmin
      ELSE
        EXIT
      END IF

    END DO ! DO WHILE (ncycle < mesh%nTri)

    ! Check if the result from the linear search is correct.
    q = mesh%V(mesh%Tri(t,1),:)
    r = mesh%V(mesh%Tri(t,2),:)
    s = mesh%V(mesh%Tri(t,3),:)
    IF (IsInTriangle(q, r, s, p)) RETURN
    IF (LiesOnLineSegment( q, r, p, mesh%tol_dist)) RETURN
    IF (LiesOnLineSegment( r, s, p, mesh%tol_dist)) RETURN
    IF (LiesOnLineSegment( s, q, p, mesh%tol_dist)) RETURN

    ! It's not. Perform a flood-fill style outward search.
    
    ! Initialise map and stack.
    mesh%TriMap     = 0
    mesh%TriStack1  = 0
    mesh%TriStack2  = 0
    mesh%TriStackN1 = 0
    mesh%TriStackN2 = 0
    mesh%TriMap(t)  = 1 ! We checked that one.

    ! Add t' neighbours to the stack.
    DO n = 1, 3
      IF (mesh%TriC(t,n) > 0) THEN
        mesh%TriStackN1 = mesh%TriStackN1+1
        mesh%TriStack1( mesh%TriStackN1) = mesh%TriC(t,n)
      END IF
    END DO

    FoundIt = .FALSE.
    DO WHILE (.NOT. FoundIt)
      ! Check all triangles in the stack. If they're not it, add their
      ! non-checked neighbours to the new stack.

      mesh%TriStack2  = 0
      mesh%TriStackN2 = 0

      DO n = 1, mesh%TriStackN1
        ti = mesh%TriStack1(n)
        q = mesh%V(mesh%Tri(ti,1),:)
        r = mesh%V(mesh%Tri(ti,2),:)
        s = mesh%V(mesh%Tri(ti,3),:)
        IF (IsInTriangle(q, r, s, p)) THEN

          ! Found it!
          t = ti

          RETURN

        ELSE ! if (IsInTriangle(ti,p))
          ! Did not find it. And add this triangle's non-checked neighbours to the new stack.

          DO n2 = 1, 3
            tin = mesh%TriC(ti,n2)
            IF (tin==0)              CYCLE ! This neighbour doesn't exist.
            IF (mesh%TriMap(tin)==1) CYCLE ! This neighbour has already been checked or is already in the stack.
            mesh%TriStackN2 = mesh%TriStackN2 + 1
            mesh%TriStack2( mesh%TriStackN2) = tin
            mesh%TriMap(tin) = 1
          END DO

        END IF ! IF (IsInTriangle(q, r, s, p, tol)) THEN
      END DO ! DO n = 1, mesh%triStackN1

      ! Cycle stacks.
      mesh%TriStack1  = mesh%TriStack2
      mesh%TriStackN1 = mesh%TriStackN2

      ! If no more non-checked neighbours could be found, terminate and throw an error.
      IF (mesh%TriStackN2==0) THEN
        WRITE(0,*) 'FindContainingTriangle - ERROR: Couldnt find triangle containing this point!'
        STOP
      END IF

    END DO ! DO WHILE (.NOT. FoundIt)


  END SUBROUTINE FindContainingTriangle
  SUBROUTINE FindContainingTriangle_old( mesh, p, t)
   ! Start at initial guess t, search outward from there using a
   ! flood-fill algorithm.

    TYPE(type_mesh),          INTENT(INOUT)       :: mesh
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: t

    REAL(dp), DIMENSION(2)                        :: q, r, s
    INTEGER                                       :: ncycle, t_prev
    REAL(dp), DIMENSION(2)                        :: gcti, gctc
    REAL(dp)                                      :: d, dc, dcmin
    INTEGER                                       :: tc, tcmin
    INTEGER                                       :: n, n2, ti, tin
    LOGICAL                                       :: FoundIt
    
    ! If p lies outside the mesh domain, throw and error
    IF (p(1) < mesh%xmin .OR. p(1) > mesh%xmax .OR. p(2) < mesh%ymin .OR. p(2) > mesh%ymax) THEN
      WRITE(0,*) 'FindContainingTriangle - ERROR: point lies outside mesh domain!'
      STOP
    END IF

    ! See if the initial guess is correct.
    q = mesh%V(mesh%Tri(t,1),:)
    r = mesh%V(mesh%Tri(t,2),:)
    s = mesh%V(mesh%Tri(t,3),:)
    IF (IsInTriangle(q, r, s, p)) RETURN

    ! If not, start with a linear search.
    ncycle = 0
    t_prev = t
    DO WHILE (ncycle < mesh%nTri)

      gcti = (mesh%V(mesh%Tri(t,1),:) + mesh%V(mesh%Tri(t,2),:) + mesh%V(mesh%Tri(t,3),:)) / 3._dp
      d = NORM2(gcti - p)

      dcmin = d + 10._dp
      tcmin = 0
      DO n = 1, 3
        tc   = mesh%TriC(t,n)
        IF (tc==0)      CYCLE ! This triangle neighbour doesn't exist
        IF (tc==t_prev) CYCLE ! This triangle neighbour is the one we just came from
        gctc = (mesh%V(mesh%Tri(tc,1),:) + mesh%V(mesh%Tri(tc,2),:) + mesh%V(mesh%Tri(tc,3),:)) / 3._dp
        dc = NORM2( gctc - p)
        IF (dc < dcmin) THEN
          dcmin = dc
          tcmin = tc
        END IF
      END DO

      IF (dcmin < d) THEN
        t_prev = t
        t = tcmin
      ELSE
        EXIT
      END IF

    END DO ! DO WHILE (ncycle < mesh%nTri)

    ! Check if the result from the linear search is correct (for triangles,
    ! this is not necessarily true, though it will be very close).
    q = mesh%V(mesh%Tri(t,1),:)
    r = mesh%V(mesh%Tri(t,2),:)
    s = mesh%V(mesh%Tri(t,3),:)
    IF (IsInTriangle(q, r, s, p)) RETURN

    ! It's not. Perform a flood-fill style outward search.
    
    ! Initialise map and stack.
    mesh%TriMap     = 0
    mesh%TriStack1  = 0
    mesh%TriStack2  = 0
    mesh%TriStackN1 = 0
    mesh%TriStackN2 = 0
    mesh%TriMap(t)  = 1 ! We checked that one.

    ! Add tguess' neighbours to the stack.
    DO n = 1, 3
      IF (mesh%TriC(t,n) > 0) THEN
        mesh%TriStackN1 = mesh%TriStackN1+1
        mesh%TriStack1( mesh%TriStackN1) = mesh%TriC(t,n)
      END IF
    END DO

    FoundIt = .FALSE.
    DO WHILE (.NOT. FoundIt)
      ! Check all triangles in the stack. If they're not it, add their
      ! non-checked neighbours to the new stack.

      mesh%TriStack2  = 0
      mesh%TriStackN2 = 0

      DO n = 1, mesh%TriStackN1
        ti = mesh%TriStack1(n)
        q = mesh%V(mesh%Tri(ti,1),:)
        r = mesh%V(mesh%Tri(ti,2),:)
        s = mesh%V(mesh%Tri(ti,3),:)
        IF (IsInTriangle(q, r, s, p)) THEN

          ! Found it!
          t = ti

          RETURN

        ELSE ! if (IsInTriangle(ti,p))
          ! Did not find it. And add this triangle's non-checked neighbours to the new stack.

          DO n2 = 1, 3
            tin = mesh%TriC(ti,n2)
            IF (tin==0)              CYCLE ! This neighbour doesn't exist.
            IF (mesh%TriMap(tin)==1) CYCLE ! This neighbour has already been checked or is already in the stack.
            mesh%TriStackN2 = mesh%TriStackN2 + 1
            mesh%TriStack2( mesh%TriStackN2) = tin
            mesh%TriMap(tin) = 1
          END DO

        END IF ! IF (IsInTriangle(q, r, s, p, tol)) THEN
      END DO ! DO n = 1, mesh%triStackN1

      ! Cycle stacks.
      mesh%TriStack1  = mesh%TriStack2
      mesh%TriStackN1 = mesh%TriStackN2

      ! If no more non-checked neighbours could be found, terminate and throw an error.
      IF (mesh%TriStackN2==0) THEN
        WRITE(0,*) 'FindContainingTriangle - ERROR: Couldnt find triangle containing this point!'
        STOP
      END IF

    END DO ! DO WHILE (.NOT. FoundIt)

  END SUBROUTINE FindContainingTriangle_old
  SUBROUTINE FindContainingVertex(mesh, p, vi)
    ! Find the vertex whose Voronoi cell contains the point p, using a "linear search"
    ! Start at initial guess vi. Check all neighbours of vi, find the one
    ! closest to p, select that one as the new vi. Repeat until all neighbours
    ! of vi are further away from p than vi itself.
    
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(  2), INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: vi

    INTEGER                                       :: ncycle, vi_prev, ci, vc, vcmin
    REAL(dp)                                      :: d, dc, dcmin

    
    ncycle = 0
    vi_prev = vi
    DO WHILE (ncycle < mesh%nV)
      
      d = NORM2( mesh%V(vi,:) - p)
      
      dcmin = d + 10._dp
      vcmin = 0
      DO ci = 1, mesh%nC(vi)
        vc = mesh%C(vi,ci)
        IF (vc==vi_prev) CYCLE ! This is the neighbour we just came from
        dc = NORM2( mesh%V(vc,:) - p)
        IF (dc < dcmin) THEN
          dcmin = dc
          vcmin = vc
        END IF
      END DO
      
      IF (dcmin < d) THEN
        vi_prev = vi
        vi = vcmin
      ELSE
        RETURN
      END IF
      
    END DO ! DO WHILE (ncycle < mesh%nV)
    
    ! If we reach this point, we didnt find the containing vertex - should not be possible, so throw an error
    WRITE(0,*) ' FindContainingVertex - ERROR: couldnt find closest vertex!'
    STOP

  END SUBROUTINE FindContainingVertex
  FUNCTION   IsInVoronoiCell(mesh, p, vi) RESULT(isso)
    ! Checks whether of not the point p lies in the Voronoi cell of vertex vi

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(  2), INTENT(IN)          :: p
    INTEGER,                  INTENT(IN)          :: vi
    
    ! Local variables:
    LOGICAL                                       :: isso
    INTEGER                                       :: nVor, n, nprev
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Vor
    REAL(dp), DIMENSION(2)                        :: ta, tb, tc

    isso = .FALSE.

    ALLOCATE(Vor(mesh%nconmax+2,2))

    CALL FindVoronoiCellVertices(mesh, vi, Vor, nVor)

    DO n = 1, nVor
      nprev = n-1
      IF (nprev==0) nprev = nVor

      ta = Vor(n,:)
      tb = Vor(nprev,:)
      tc = mesh%V(vi,:)
      IF (IsInTriangle(ta, tb, tc, p)) THEN
        isso = .TRUE.
        RETURN
      END IF
    END DO

    DEALLOCATE(Vor)

  END FUNCTION IsInVoronoiCell
  FUNCTION LiesOnLineSegment(pa, pb, pc, tol_dist) RESULT(isso)
    ! Test if the point pc lies on the line pb-pc, or within tol_dist of it
    
    ! In/output variables:
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: pa, pb, pc
    REAL(dp),                 INTENT(IN)          :: tol_dist
    LOGICAL                                       :: isso
    
    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: d, e, d_norm, e_par, e_ort
    
    d = pb - pa
    e = pc - pa

    d_norm = d / NORM2(d)

    e_par = (e(1)*d_norm(1) + e(2)*d_norm(2)) * d_norm
    e_ort = e - e_par

    isso = .TRUE.
    IF (NORM2(e_ort) > tol_dist) THEN
      isso = .FALSE.
      RETURN
    END IF

    IF ((e(1)*d(1) + e(2)*d(2)) > 0._dp) THEN
      IF (NORM2(e_par) > (NORM2(d))) THEN
        isso = .FALSE.
        RETURN
      END IF
    ELSE
      IF (NORM2(e_par) > 0._dp) THEN
        isso = .FALSE.
        RETURN
      END IF
    END IF

  END FUNCTION LiesOnLineSegment
  SUBROUTINE SegmentIntersection( p, q, r, s, llis, do_cross, tol_dist)
    ! Find out if the line segments [pq] and [rs] intersect. If so, return
    ! the coordinates of the point of intersection
    
    ! In/output variables:
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p, q, r, s
    REAL(dp), DIMENSION(2),   INTENT(OUT)         :: llis
    LOGICAL,                  INTENT(OUT)         :: do_cross
    REAL(dp),                 INTENT(IN)          :: tol_dist
    
    ! Local variables:
    REAL(dp), DIMENSION(2,2)                      :: A
    REAL(dp), DIMENSION(2)                        :: x, b
    INTEGER,  DIMENSION(2)                        :: IPIV
    INTEGER                                       :: info

    ! External subroutines:      
    EXTERNAL DGESV ! Lapack routine that solves matrix equation Ax=b for x (in double precision)
    
    ! If pq and rs are colinear, define them as not intersecting
    IF ( NORM2( CROSS( [q(1)-p(1),q(2)-p(2),0._dp], [s(1)-r(1),s(2)-r(2),0._dp] )) < tol_dist) THEN
      llis = [0._dp, 0._dp]
      do_cross = .FALSE.
      RETURN
    END IF
    
    A(1,:) = [(p(1)-q(1)), (r(1)-s(1))]
    A(2,:) = [(p(2)-q(2)), (r(2)-s(2))]
    b = [(r(1)-q(1)), (r(2)-q(2))]
    
    ! The LAPACK solver will overwrite the right-hand side b with the solution x. Therefore we 
    ! first copy the rhs in the solution vector x:
    x = b
    
    ! Solve Ax = b using LAPACK
    CALL DGESV( 2, 1, A, 2, IPIV, x, 2, info)
    
    llis = [q(1) + x(1) * (p(1)-q(1)), q(2) + x(1) * (p(2)-q(2))]
    
    IF (x(1)>0._dp .AND. x(1)<1._dp .AND. x(2)>0._dp .AND. x(2)<1._dp) THEN
      do_cross = .TRUE.
    ELSE
      do_cross = .FALSE.
    END IF
    
  END SUBROUTINE SegmentIntersection
    
! == Some very basic geometry, mainly used for determining triangle circumcenters
  SUBROUTINE LineFromPoints(p,q,la,lb,lc)

    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p, q
    REAL(dp),                   INTENT(OUT)       :: la, lb, lc

    ! Line PQ is represented as la*x + lb*y = lc
    la = q(2) - p(2)
    lb = p(1) - q(1)
    lc = la*(p(1))+ lb*(p(2));

  END SUBROUTINE LineFromPoints  
  SUBROUTINE PerpendicularBisectorFromLine(p,q,la1,lb1,la2,lb2,lc2)

    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p, q
    REAL(dp),                   INTENT(IN)        :: la1, lb1
    REAL(dp),                   INTENT(OUT)       :: la2, lb2, lc2
    REAL(dp)                                      :: temp
    REAL(dp), DIMENSION(2)                        :: m

    m = (p+q)/2
    lc2 = -lb1*m(1) + la1*m(2)

    temp = la1
    la2 = -lb1
    lb2 = temp

  END SUBROUTINE PerpendicularBisectorFromLine  
  SUBROUTINE LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis)
    ! Find the intersection llis of the lines la1*x+lb1*y=lc1 and la2*x+lb2*y=lc2

    REAL(dp),                   INTENT(IN)        :: la1, lb1, lc1, la2, lb2, lc2
    REAL(dp), DIMENSION(2),     INTENT(OUT)       :: llis
    REAL(dp)                                      :: d

    d = la1*lb2 - la2*lb1
    IF (d == 0) THEN
      ! The lines are parallel.
      llis = [1E30, 1E30]
    ELSE
      llis = [(lb2*lc1 - lb1*lc2), (la1*lc2 - la2*lc1)]/d
    END IF
  END SUBROUTINE LineLineIntersection  
  SUBROUTINE FindCircumcenter(p,q,r,cc)
    ! Find the circumcenter cc of the triangle pqr

    ! Some basic vector operations
    ! Find the circumcenter cc of the triangle pqr
    ! If pqr are colinear, returns [1e30,1e30]

    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p, q, r
    REAL(dp), DIMENSION(2),     INTENT(OUT)       :: cc
    REAL(dp)                                      :: la1,lb1,lc1,le1,lf1,lg1
    REAL(dp)                                      :: la2,lb2,lc2,le2,lf2,lg2

    cc = [0._dp, 0._dp]

    ! Line PQ is represented as ax + by = c, Line QR is represented as ex + fy = g
    CALL LineFromPoints(p,q,la1,lb1,lc1)
    CALL LineFromPoints(q,r,le1,lf1,lg1);

    ! Converting lines PQ and QR to perpendicular
    ! bisectors. After this, L = ax + by = c
    ! M = ex + fy = g
    CALL PerpendicularBisectorFromLine(p,q,la1,lb1,la2,lb2,lc2);
    CALL PerpendicularBisectorFromLine(q,r,le1,lf1,le2,lf2,lg2);

    ! The point of intersection of L and M gives
    ! the circumcenter
    CALL LineLineIntersection(la2,lb2,lc2,le2,lf2,lg2,cc);

  END SUBROUTINE FindCircumcenter
  FUNCTION Cross(va,vb) RESULT(vc)
    ! Vector product vc between vectors va and vb
    REAL(dp), DIMENSION(3),     INTENT(IN)        :: va, vb
    REAL(dp), DIMENSION(3)                        :: vc

    vc(1) = (va(2)*vb(3)) - (va(3)*vb(2))
    vc(2) = (va(3)*vb(1)) - (va(1)*vb(3))
    vc(3) = (va(1)*vb(2)) - (va(2)*vb(1))

  END FUNCTION Cross
  FUNCTION Cross2(a,b) RESULT(z)
    ! Vector product z between 2-dimensional vectors a and b
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: a, b
    REAL(dp)                                      :: z

    z = (a(1)*b(2)) - (a(2)*b(1))

  END FUNCTION Cross2
  FUNCTION FindTriangleArea(pq,pr,ps) RESULT(TriA)
    ! Find the area of the triangle [pq,pr,ps]
    
    REAL(dp), DIMENSION(2), INTENT(IN)  :: pq, pr, ps
    REAL(dp)                            :: TriA
    
    TriA = NORM2( Cross( [pr(1)-pq(1), pr(2)-pq(2), 0._dp], [ps(1)-pq(1), ps(2)-pq(2), 0._dp] )) / 2._dp
    
  END FUNCTION FindTriangleArea

! == Three basic line integrals used for conservative remapping
  FUNCTION LineIntegral_xdy( p, q, tol_dist) RESULT(I_pq)
    ! Calculate the line integral x dy from p to q
    
    ! In/output variables:
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
    REAL(dp),                                INTENT(IN)    :: tol_dist
    REAL(dp)                                               :: I_pq
    
    ! Local variables:
    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy
        
    xp = p(1)
    yp = p(2)
    xq = q(1)
    yq = q(2)
    
    IF (ABS(yp-yq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF
    
    dx = q(1)-p(1)
    dy = q(2)-p(2)
    
    I_pq = xp*dy - yp*dx + (dx / (2._dp*dy)) * (yq**2 - yp**2)
    
  END FUNCTION LineIntegral_xdy
  FUNCTION LineIntegral_mxydx( p, q, tol_dist) RESULT(I_pq)
    ! Calculate the line integral -xy dx from p to q
    
    ! In/output variables:
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
    REAL(dp),                                INTENT(IN)    :: tol_dist
    REAL(dp)                                               :: I_pq
    
    ! Local variables:
    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy
        
    xp = p(1)
    yp = p(2)
    xq = q(1)
    yq = q(2)
    
    IF (ABS(xp-xq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF
    
    dx = q(1)-p(1)
    dy = q(2)-p(2)
    
    I_pq = (1._dp/2._dp * (xp*dy/dx - yp) * (xq**2-xp**2)) - (1._dp/3._dp * dy/dx * (xq**3-xp**3))
    
  END FUNCTION LineIntegral_mxydx
  FUNCTION LineIntegral_xydy( p, q, tol_dist) RESULT(I_pq)
    ! Calculate the line integral xy dy from p to q
    
    ! In/output variables:
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
    REAL(dp),                                INTENT(IN)    :: tol_dist
    REAL(dp)                                               :: I_pq
    
    ! Local variables:
    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy
        
    xp = p(1)
    yp = p(2)
    xq = q(1)
    yq = q(2)
    
    IF (ABS(yp-yq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF
    
    dx = q(1)-p(1)
    dy = q(2)-p(2)
    
    I_pq = (1._dp/2._dp * (xp - yp*dx/dy) * (yq**2-yp**2)) + (1._dp/3._dp * dx/dy * (yq**3-yp**3))
    
  END FUNCTION LineIntegral_xydy
      
! == Help routines for use in mesh updating
  SUBROUTINE MeanCartOverVoronoiCell_dp( mesh, d, x, y, nx, ny, vi, v)
    ! Find the mean value of the data field d, given on the cartesian
    ! grid x,y, over the Voronoi cell of vertex vi

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(:,:), INTENT(IN)          :: d          ! Data field
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: x          ! Data x grid
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: y          ! Data y grid
    INTEGER,                  INTENT(IN)          :: nx         ! Number of x elements
    INTEGER,                  INTENT(IN)          :: ny         ! Number of y elements
    INTEGER,                  INTENT(IN)          :: vi

    REAL(dp),                 INTENT(OUT)         :: v

    REAL(dp)                                      :: sumel, trisumel
    INTEGER                                       :: nel, trinel
    INTEGER                                       :: nVor, n
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Vor

    sumel = 0._dp
    nel   = 0

    ALLOCATE(Vor(mesh%nconmax+2,2))

    CALL FindVoronoiCellVertices(mesh, vi, Vor, nVor)

    DO n = 2, nVor
      CALL SumCartOverTriangle_dp( mesh%V(vi,:), Vor(n-1,:), Vor(n,:), d, x, y, nx, ny, trisumel, trinel)
      sumel = sumel + trisumel
      nel   = nel   + trinel
    END DO

    IF (nel>4) THEN
      v = sumel / nel
    ELSE
      ! Too few elements for a proper mean - just do bicubic interpolation to
      ! the vertex location
      CALL Cart_Bilinear_dp( d, x, y, nx, ny, mesh%V(vi,:), v)
    END IF

    DEALLOCATE(Vor)

  END SUBROUTINE MeanCartOverVoronoiCell_dp
  SUBROUTINE MeanCartOverVoronoiCell_int( mesh, d, x, y, nx, ny, vi, v)
    ! Find the mean value of the data field d, given on the cartesian
    ! grid x,y, over the Voronoi cell of vertex vi

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER, DIMENSION(:,:),  INTENT(IN)          :: d          ! Data field
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: x          ! Data x grid
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: y          ! Data y grid
    INTEGER,                  INTENT(IN)          :: nx         ! Number of x elements
    INTEGER,                  INTENT(IN)          :: ny         ! Number of y elements
    INTEGER,                  INTENT(IN)          :: vi

    REAL(dp),                 INTENT(OUT)         :: v
    
    INTEGER                                       :: trisumel
    REAL(dp)                                      :: sumel
    INTEGER                                       :: nel, trinel
    INTEGER                                       :: nVor, n
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Vor

    sumel = 0._dp
    nel   = 0

    ALLOCATE(Vor(mesh%nconmax+2,2))

    CALL FindVoronoiCellVertices(mesh, vi, Vor, nVor)

    DO n = 2, nVor
      CALL SumCartOverTriangle_int( mesh%V(vi,:), Vor(n-1,:), Vor(n,:), d, x, y, nx, ny, trisumel, trinel)
      sumel = sumel + REAL(trisumel)
      nel   = nel   + trinel
    END DO

    IF (nel>4) THEN
      v = sumel / REAL(nel)
    ELSE
      ! Too few elements for a proper mean - just do bicubic interpolation to
      ! the vertex location
      CALL Cart_Bilinear_int( d, x, y, nx, ny, mesh%V(vi,:), v)
    END IF

    DEALLOCATE(Vor)

  END SUBROUTINE MeanCartOverVoronoiCell_int  
  SUBROUTINE SumCartOverTriangle_dp( p1, p2, p3, d, x, y, nx, ny, trisumel, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear

    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p1,p2,p3   ! The points spanning the triangle
    REAL(dp), DIMENSION(:,:), INTENT(IN)          :: d          ! Data field
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: x          ! Data x grid
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: y          ! Data y grid
    INTEGER,                  INTENT(IN)          :: nx         ! Number of x elements
    INTEGER,                  INTENT(IN)          :: ny         ! Number of y elements

    REAL(dp),                 INTENT(OUT)         :: trisumel
    INTEGER,                  INTENT(OUT)         :: trinel

    INTEGER                                       :: il, iu, jl, ju, i, j

    trisumel = 0._dp
    trinel   = 0

    il = MAX(1,  1 +FLOOR((MINVAL([p1(1), p2(1), p3(1)])-MINVAL(x)) / (x(2)-x(1))))
    iu = MIN(nx, nx-FLOOR((MAXVAL(x)-MAXVAL([p1(1), p2(1), p3(1)])) / (x(2)-x(1))))
    jl = MAX(1,  1 +FLOOR((MINVAL([p1(2), p2(2), p3(2)])-MINVAL(y)) / (y(2)-y(1))))
    ju = MIN(ny, ny-FLOOR((MAXVAL(y)-MAXVAL([p1(2), p2(2), p3(1)])) / (y(2)-y(1))))

    DO i = il, iu
    DO j = jl, ju
     IF (IsInTriangle(p1,p2,p3,[x(i), y(j)])) THEN
       trisumel = trisumel + d(i,j)
       trinel   = trinel   + 1
     END IF
    END DO
    END DO

  END SUBROUTINE SumCartOverTriangle_dp
  SUBROUTINE SumCartOverTriangle_int( p1, p2, p3, d, x, y, nx, ny, trisumel, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear

    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p1,p2,p3   ! The points spanning the triangle
    INTEGER,  DIMENSION(:,:), INTENT(IN)          :: d          ! Data field
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: x          ! Data x grid
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: y          ! Data y grid
    INTEGER,                  INTENT(IN)          :: nx         ! Number of x elements
    INTEGER,                  INTENT(IN)          :: ny         ! Number of y elements

    INTEGER,                  INTENT(OUT)         :: trisumel
    INTEGER,                  INTENT(OUT)         :: trinel

    INTEGER                                       :: il, iu, jl, ju, i, j

    trisumel = 0
    trinel   = 0

    il = MAX(1,  1 +FLOOR((MINVAL([p1(1), p2(1), p3(1)])-MINVAL(x)) / (x(2)-x(1))))
    iu = MIN(nx, nx-FLOOR((MAXVAL(x)-MAXVAL([p1(1), p2(1), p3(1)])) / (x(2)-x(1))))
    jl = MAX(1,  1 +FLOOR((MINVAL([p1(2), p2(2), p3(2)])-MINVAL(y)) / (y(2)-y(1))))
    ju = MIN(ny, ny-FLOOR((MAXVAL(y)-MAXVAL([p1(2), p2(2), p3(1)])) / (y(2)-y(1))))

    DO i = il, iu
    DO j = jl, ju
     IF (IsInTriangle(p1,p2,p3,[x(i), y(j)])) THEN
       trisumel = trisumel + d(i,j)
       trinel   = trinel   + 1
     END IF
    END DO
    END DO

  END SUBROUTINE SumCartOverTriangle_int
  SUBROUTINE MaxCartOverTriangle_dp( p1, p2, p3, d, x, y, nx, ny, vmax, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear

    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p1,p2,p3   ! The points spanning the triangle
    REAL(dp), DIMENSION(:,:), INTENT(IN)          :: d          ! Data field
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: x          ! Data x grid
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: y          ! Data y grid
    INTEGER,                  INTENT(IN)          :: nx         ! Number of x elements
    INTEGER,                  INTENT(IN)          :: ny         ! Number of y elements

    REAL(dp),                 INTENT(OUT)         :: vmax
    INTEGER,                  INTENT(OUT)         :: trinel

    INTEGER                                       :: il, iu, jl, ju, i, j

    vmax   = -1e30
    trinel = 0

    il = MAX(1,  1 +FLOOR((MINVAL([p1(1), p2(1), p3(1)])-MINVAL(x)) / (x(2)-x(1))))
    iu = MIN(nx, nx-FLOOR((MAXVAL(x)-MAXVAL([p1(1), p2(1), p3(1)])) / (x(2)-x(1))))
    jl = MAX(1,  1 +FLOOR((MINVAL([p1(2), p2(2), p3(2)])-MINVAL(y)) / (y(2)-y(1))))
    ju = MIN(ny, ny-FLOOR((MAXVAL(y)-MAXVAL([p1(2), p2(2), p3(1)])) / (y(2)-y(1))))

    DO i = il, iu
    DO j = jl, ju
     IF (IsInTriangle(p1,p2,p3,[x(i), y(j)])) THEN
      IF (d(i,j) > vmax) THEN
       vmax   = d(i,j)
       trinel = 1
      END IF
     END IF
    END DO
    END DO

  END SUBROUTINE MaxCartOverTriangle_dp
  SUBROUTINE MaxCartOverTriangle_int( p1, p2, p3, d, x, y, nx, ny, vmax, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear

    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p1,p2,p3   ! The points spanning the triangle
    INTEGER,  DIMENSION(:,:), INTENT(IN)          :: d          ! Data field
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: x          ! Data x grid
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: y          ! Data y grid
    INTEGER,                  INTENT(IN)          :: nx         ! Number of x elements
    INTEGER,                  INTENT(IN)          :: ny         ! Number of y elements

    INTEGER,                  INTENT(OUT)         :: vmax
    INTEGER,                  INTENT(OUT)         :: trinel

    INTEGER                                       :: il, iu, jl, ju, i, j

    vmax   = -100000000
    trinel = 0

    il = MAX(1,  1 +FLOOR((MINVAL([p1(1), p2(1), p3(1)])-MINVAL(x)) / (x(2)-x(1))))
    iu = MIN(nx, nx-FLOOR((MAXVAL(x)-MAXVAL([p1(1), p2(1), p3(1)])) / (x(2)-x(1))))
    jl = MAX(1,  1 +FLOOR((MINVAL([p1(2), p2(2), p3(2)])-MINVAL(y)) / (y(2)-y(1))))
    ju = MIN(ny, ny-FLOOR((MAXVAL(y)-MAXVAL([p1(2), p2(2), p3(1)])) / (y(2)-y(1))))

    DO i = il, iu
    DO j = jl, ju
     IF (IsInTriangle(p1,p2,p3,[x(i), y(j)])) THEN
      IF (d(i,j) > vmax) THEN
       vmax   = d(i,j)
       trinel = 1
      END IF
     END IF
    END DO
    END DO

  END SUBROUTINE MaxCartOverTriangle_int
  SUBROUTINE MinCartOverTriangle_dp( p1, p2, p3, d, x, y, nx, ny, vmin, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear

    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p1,p2,p3   ! The points spanning the triangle
    REAL(dp), DIMENSION(:,:), INTENT(IN)          :: d          ! Data field
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: x          ! Data x grid
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: y          ! Data y grid
    INTEGER,                  INTENT(IN)          :: nx         ! Number of x elements
    INTEGER,                  INTENT(IN)          :: ny         ! Number of y elements

    REAL(dp),                 INTENT(OUT)         :: vmin
    INTEGER,                  INTENT(OUT)         :: trinel

    INTEGER                                       :: il, iu, jl, ju, i, j

    vmin   = 1e30
    trinel = 0

    il = MAX(1,  1 +FLOOR((MINVAL([p1(1), p2(1), p3(1)])-MINVAL(x)) / (x(2)-x(1))))
    iu = MIN(nx, nx-FLOOR((MAXVAL(x)-MAXVAL([p1(1), p2(1), p3(1)])) / (x(2)-x(1))))
    jl = MAX(1,  1 +FLOOR((MINVAL([p1(2), p2(2), p3(2)])-MINVAL(y)) / (y(2)-y(1))))
    ju = MIN(ny, ny-FLOOR((MAXVAL(y)-MAXVAL([p1(2), p2(2), p3(1)])) / (y(2)-y(1))))

    DO i = il, iu
    DO j = jl, ju
     IF (IsInTriangle(p1,p2,p3,[x(i), y(j)])) THEN
      IF (d(i,j) < vmin) THEN
       vmin   = d(i,j)
       trinel = 1
      END IF
     END IF
    END DO
    END DO

  END SUBROUTINE MinCartOverTriangle_dp
  SUBROUTINE MinCartOverTriangle_int( p1, p2, p3, d, x, y, nx, ny, vmin, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear

    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p1,p2,p3   ! The points spanning the triangle
    INTEGER,  DIMENSION(:,:), INTENT(IN)          :: d          ! Data field
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: x          ! Data x grid
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: y          ! Data y grid
    INTEGER,                  INTENT(IN)          :: nx         ! Number of x elements
    INTEGER,                  INTENT(IN)          :: ny         ! Number of y elements

    INTEGER,                  INTENT(OUT)         :: vmin
    INTEGER,                  INTENT(OUT)         :: trinel

    INTEGER                                       :: il, iu, jl, ju, i, j

    vmin   = 100000000
    trinel = 0

    il = MAX(1,  1 +FLOOR((MINVAL([p1(1), p2(1), p3(1)])-MINVAL(x)) / (x(2)-x(1))))
    iu = MIN(nx, nx-FLOOR((MAXVAL(x)-MAXVAL([p1(1), p2(1), p3(1)])) / (x(2)-x(1))))
    jl = MAX(1,  1 +FLOOR((MINVAL([p1(2), p2(2), p3(2)])-MINVAL(y)) / (y(2)-y(1))))
    ju = MIN(ny, ny-FLOOR((MAXVAL(y)-MAXVAL([p1(2), p2(2), p3(1)])) / (y(2)-y(1))))

    DO i = il, iu
    DO j = jl, ju
     IF (IsInTriangle(p1,p2,p3,[x(i), y(j)])) THEN
      IF (d(i,j) < vmin) THEN
       vmin   = d(i,j)
       trinel = 1
      END IF
     END IF
    END DO
    END DO

  END SUBROUTINE MinCartOverTriangle_int
  SUBROUTINE Cart_Bilinear_dp( d, x, y, nx, ny, p, v)
    ! Bicubic interpolation of Cartesian data
    ! Interpolates the data field d of size nx,ny, given on grid x,y,
    ! at the point p, resulting in value v

    REAL(dp), DIMENSION(:,:), INTENT(IN)          :: d   ! Data to interpolate
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: x   ! Data x grid
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: y   ! Data y grid
    INTEGER,                  INTENT(IN)          :: nx  ! Number of x elements
    INTEGER,                  INTENT(IN)          :: ny  ! Number of y elements
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p   ! Location where we want to know interpolated value

    REAL(dp),                 INTENT(OUT)         :: v   ! Interpolated value

    INTEGER                                       :: il, iu, jl, ju
    REAL(dp)                                      :: wil, wiu, wjl, wju

    ! x and y indices of data grid surrounding p
    il = MAX(1,MIN(nx-1, 1 + FLOOR((p(1)-MINVAL(x)) / (x(2)-x(1)))))
    iu = il+1
    jl = MAX(1,MIN(ny-1, 1 + FLOOR((p(2)-MINVAL(y)) / (y(2)-y(1)))))
    ju = jl+1

    ! Relative contributions from these four points
    wil = (x(iu) - p(1))/(x(2)-x(1))
    wiu = 1-wil
    wjl = (y(ju) - p(2))/(y(2)-y(1))
    wju = 1-wjl

    ! Bicubic interpolation
    v =  (wil * wjl * d(il,jl)) + &
         (wil * wju * d(il,ju)) + &
         (wiu * wjl * d(iu,jl)) + &
         (wiu * wju * d(iu,ju))

  END SUBROUTINE Cart_Bilinear_dp
  SUBROUTINE Cart_Bilinear_int( d, x, y, nx, ny, p, v)
    ! Bicubic interpolation of Cartesian data
    ! Interpolates the data field d of size nx,ny, given on grid x,y,
    ! at the point p, resulting in value v

    INTEGER, DIMENSION(:,:), INTENT(IN)          :: d   ! Data to interpolate
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: x   ! Data x grid
    REAL(dp), DIMENSION(:  ), INTENT(IN)          :: y   ! Data y grid
    INTEGER,                  INTENT(IN)          :: nx  ! Number of x elements
    INTEGER,                  INTENT(IN)          :: ny  ! Number of y elements
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p   ! Location where we want to know interpolated value

    REAL(dp),                 INTENT(OUT)         :: v   ! Interpolated value

    INTEGER                                       :: il, iu, jl, ju
    REAL(dp)                                      :: wil, wiu, wjl, wju

    ! x and y indices of data grid surrounding p
    il = MAX(1,MIN(nx-1, 1 + FLOOR((p(1)-MINVAL(x)) / (x(2)-x(1)))))
    iu = il+1
    jl = MAX(1,MIN(ny-1, 1 + FLOOR((p(2)-MINVAL(y)) / (y(2)-y(1)))))
    ju = jl+1

    ! Relative contributions from these four points
    wil = (x(iu) - p(1))/(x(2)-x(1))
    wiu = 1-wil
    wjl = (y(ju) - p(2))/(y(2)-y(1))
    wju = 1-wjl

    ! Bicubic interpolation
    v =  (wil * wjl * REAL(d(il,jl))) + &
         (wil * wju * REAL(d(il,ju))) + &
         (wiu * wjl * REAL(d(iu,jl))) + &
         (wiu * wju * REAL(d(iu,ju)))

  END SUBROUTINE Cart_Bilinear_int
  SUBROUTINE Mesh_Bilinear_dp(mesh, d, p, ti, v)
    ! Interpolate data field d, given on some mesh, to point p, resulting in v
    ! Guess that p is in triangle ti, as a start for FindContainingTriangle_oldmesh

    TYPE(type_mesh),          INTENT(INOUT)       :: mesh
    REAL(dp), DIMENSION(:),   INTENT(IN)          :: d
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: ti
    REAL(dp),                 INTENT(OUT)         :: v

    INTEGER                                       :: v1, v2, v3
    REAL(dp)                                      :: dx, dy, ddx, ddy

    ! Find the triangle containing p
    CALL FindContainingTriangle(mesh, p, ti)

    ! Find function values and derivatives of d on the points of ti
    v1 = mesh%Tri(ti,1)
    v2 = mesh%Tri(ti,2)
    v3 = mesh%Tri(ti,3)

    dx = p(1) - mesh%V(v1,1)
    dy = p(2) - mesh%V(v1,2)
     
    ddx = mesh%NxTri(ti,1) * d(v1) + mesh%NxTri(ti,2) * d(v2) + mesh%NxTri(ti,3) * d(v3)
    ddy = mesh%NyTri(ti,1) * d(v1) + mesh%NyTri(ti,2) * d(v2) + mesh%NyTri(ti,3) * d(v3)
     
    v = d(v1) + (dx * ddx) + (dy * ddy)

  END SUBROUTINE Mesh_Bilinear_dp
  SUBROUTINE Mesh_Bilinear_int(mesh, d, p, ti, v)
    ! Interpolate data field d, given on some mesh, to point p, resulting in v
    ! Guess that p is in triangle ti, as a start for FindContainingTriangle_oldmesh

    TYPE(type_mesh),          INTENT(INOUT)       :: mesh
    INTEGER,  DIMENSION(:),   INTENT(IN)          :: d
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: ti
    REAL(dp),                 INTENT(OUT)         :: v

    INTEGER                                       :: v1, v2, v3
    REAL(dp)                                      :: dx, dy, ddx, ddy

    ! Find the triangle containing p
    CALL FindContainingTriangle(mesh, p, ti)

    ! Find function values and derivatives of d on the points of ti
    v1 = mesh%Tri(ti,1)
    v2 = mesh%Tri(ti,2)
    v3 = mesh%Tri(ti,3)

    dx = p(1) - mesh%V(v1,1)
    dy = p(2) - mesh%V(v1,2)
     
    ddx = mesh%NxTri(ti,1) * REAL(d(v1)) + mesh%NxTri(ti,2) * REAL(d(v2)) + mesh%NxTri(ti,3) * REAL(d(v3))
    ddy = mesh%NyTri(ti,1) * REAL(d(v1)) + mesh%NyTri(ti,2) * REAL(d(v2)) + mesh%NyTri(ti,3) * REAL(d(v3))
     
    v = d(v1) + (dx * ddx) + (dy * ddy)

  END SUBROUTINE Mesh_Bilinear_int
  SUBROUTINE NewTriangleContainsOldMask( mesh, pa, pb, pc, mask, vi_closest_to_gc, isso)
    ! Used in mesh updating. Given a (new mesh) triangle [pa,pb,pc], check if that triangle
    ! contains any (old) mesh vertices where the (old mesh) mask has value 1

    ! In/output variables:
    TYPE(type_mesh),          INTENT(INOUT)       :: mesh
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: pa, pb, pc
    INTEGER,  DIMENSION(:),   INTENT(IN)          :: mask
    INTEGER,                  INTENT(INOUT)       :: vi_closest_to_gc
    LOGICAL,                  INTENT(OUT)         :: isso
    
    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: gc, p
    INTEGER                                       :: ti_containing_gc, vi1, vi2, vi3
    REAL(dp)                                      :: dx, dy, mdp1, mdp2, mdp3, ddx, ddy, mdpgc
    INTEGER                                       :: nvi, vi, ci, vc
    
    isso = .FALSE.
    
    ! Use a linear search to find the (old) mesh vertex closest to the geometric center of the triangle
    gc = (pa+pb+pc) / 3._dp
    CALL FindContainingVertex( mesh, gc, vi_closest_to_gc)
    
    ! If that vertex doesn't lie inside the triangle, then probably none of them do - use trilinear interpolation instead.
    p = mesh%V(vi_closest_to_gc,:)
    IF (.NOT. IsInTriangle( pa, pb, pc, p)) THEN
    
      ! Find the (old) mesh triangle containing the (new mesh) triangle's geometric center
      ti_containing_gc = mesh%iTri(vi_closest_to_gc,1)
      CALL FindContainingTriangle( mesh, gc, ti_containing_gc)

      ! The three vertices spanning this (old) mesh triangle
      vi1 = mesh%Tri(ti_containing_gc,1)
      vi2 = mesh%Tri(ti_containing_gc,2)
      vi3 = mesh%Tri(ti_containing_gc,3)

      ! Find the mask value at gc by interpolating between these three vertices.
      dx = gc(1) - mesh%V(vi1,1)
      dy = gc(2) - mesh%V(vi1,2)
      
      mdp1 = REAL(mask(vi1))
      mdp2 = REAL(mask(vi2))
      mdp3 = REAL(mask(vi3))
     
      ddx = mesh%NxTri(ti_containing_gc,1) * mdp1 + mesh%NxTri(ti_containing_gc,2) * mdp2 + mesh%NxTri(ti_containing_gc,3) * mdp3
      ddy = mesh%NyTri(ti_containing_gc,1) * mdp1 + mesh%NyTri(ti_containing_gc,2) * mdp2 + mesh%NyTri(ti_containing_gc,3) * mdp3
     
      mdpgc = mdp1 + (dx * ddx) + (dy * ddy)
      
      IF (mdpgc > 0.1_dp) THEN
        isso = .TRUE.
      ELSE
        isso = .FALSE.
      END IF
      
      RETURN
    
    END IF
    
    ! That (old) mesh vertex lies inside the (new mesh) triangle. Use a FloodFill search to find all
    ! the (old) mesh vertices that lie inside the (new mesh) triangle. If any of them have a mask value of 1, return immediately.
    
    IF (mask(vi_closest_to_gc)==1) THEN
      isso = .TRUE.
      RETURN
    END IF
    
    mesh%VMap    = 0
    mesh%VStack1 = 0
    mesh%VStack2 = 0
    
    mesh%VMap( vi_closest_to_gc) = 1
    mesh%VStack1( 1) = vi_closest_to_gc
    mesh%VStackN1 = 1
    
    DO WHILE (mesh%VStackN1 > 0)
    
      ! Clean up the new stack
      mesh%VStack2  = 0
      mesh%VStackN2 = 0
    
      ! Go over the entire old stack. Add any non-checked neighbours of these stack elements
      ! to the new stack. When finished, cycle the stacks.
      DO nvi = 1, mesh%VStackN1
        vi = mesh%VStack1(nvi)
        
        DO ci = 1, mesh%nC(vi)
        
          vc = mesh%C(vi,ci)
          p  = mesh%V(vc,:)
          
          IF (mesh%VMap(vc)==1) CYCLE ! This neighbour has already been checked.
          
          ! Mark this neighbour as checked
          mesh%VMap(vc) = 1
          
          ! If it lies inside the triangle, add it to the new stack.
          IF (IsInTriangle(pa, pb, pc, p)) THEN
          
            ! If it has what we're looking for, stop the search.
            IF (mask(vc) == 1) THEN
              isso = .TRUE.
              RETURN
            END IF
            mesh%VStackN2 = mesh%VStackN2 + 1
            mesh%VStack2( mesh%VStackN2) = vc
          END IF
          
        END DO ! DO ci = 1, mesh%nC(vi)
        
      END DO ! DO nvi = 1, mesh%VStackN1
      
      ! Cycle the stacks
      mesh%VStack1  = mesh%VStack2
      mesh%VStackN1 = mesh%VStackN2
    
    END DO ! DO WHILE (mesh%VStackN1) > 0    
    
  END SUBROUTINE NewTriangleContainsOldMask
  
! == Diagnostic tools: write a (small) mesh to the screen, and check if mesh data is self-consistent
  SUBROUTINE WriteMeshToScreen( mesh)

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER                                       :: vi, ti
    
    WRITE(0,*) '============================================================================'
    WRITE(0,*) ' xmin = ', mesh%xmin
    WRITE(0,*) ' xmax = ', mesh%xmax
    WRITE(0,*) ' ymin = ', mesh%ymin
    WRITE(0,*) ' ymax = ', mesh%ymax
    WRITE(0,*) ' vi    nC             C            niTri         iTri         edge_index         x              y'
    DO vi = 1, mesh%nV
      WRITE(0,'(A,I3,A,I3,A,6I3,A,I3,A,6I3,A,I3,A,F12.1,A,F12.1)') &
      ' ', vi, '   ', mesh%nC(vi), '    ', mesh%C(vi,1:6), '    ', mesh%niTri(vi), '    ', mesh%iTri(vi,1:6), '    ', mesh%edge_index(vi), &
      '    ', mesh%V(vi,1), '    ', mesh%V(vi,2)
    END DO
    
    WRITE(0,*) ' ti       Tri         TriC     Tri_edge_index'
    DO ti = 1, mesh%nTri
      WRITE(0,'(A,I3,A,3I3,A,3I3,A,I3)') &
      ' ', ti, '   ', mesh%Tri(ti,:), '    ', mesh%TriC(ti,:), '    ', mesh%Tri_edge_index(ti)
    END DO
    WRITE(0,*) '============================================================================'
  
  END SUBROUTINE WriteMeshToScreen
  SUBROUTINE WriteMeshToTextFile( mesh, filename)

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    CHARACTER(LEN=*),         INTENT(IN)          :: filename
    INTEGER                                       :: vi, ci, ti
       
    ! Create a new text file
    OPEN(UNIT  = 1337, FILE = (TRIM(C%output_dir) // filename), STATUS = 'REPLACE')
    
    ! Header
    WRITE(UNIT = 1337, FMT = '(A)')       ' Mesh data'
    WRITE(UNIT = 1337, FMT = '(A,F14.4)') '  xmin    = ', mesh%xmin
    WRITE(UNIT = 1337, FMT = '(A,F14.4)') '  xmax    = ', mesh%xmax
    WRITE(UNIT = 1337, FMT = '(A,F14.4)') '  ymin    = ', mesh%ymin
    WRITE(UNIT = 1337, FMT = '(A,F14.4)') '  ymax    = ', mesh%ymax
    WRITE(UNIT = 1337, FMT = '(A,I2)')    '  nconmax = ', mesh%nconmax
    WRITE(UNIT = 1337, FMT = '(A,I6)')    '  nV      = ', mesh%nV
    WRITE(UNIT = 1337, FMT = '(A,I6)')    '  nTri    = ', mesh%nTri
    WRITE(UNIT = 1337, FMT = '(A)')       ''
    WRITE(UNIT = 1337, FMT = '(A,I6,A,I3,A)')       'Vertex data: ', mesh%nV, ' rows, ', 2 + 1 + mesh%nconmax + 1 + mesh%nconmax + 1, ' columns'
    WRITE(UNIT = 1337, FMT = '(A)')       ''
    
    ! Vertex data
    WRITE(UNIT = 1337, FMT = '(A)')       'V  nC  C  niTri  iTri  edge_index'    
    DO vi = 1, mesh%nV
      WRITE(UNIT = 1337, FMT = '(2F24.14,I3)', ADVANCE = 'NO') mesh%V(vi,1), mesh%V(vi,2), mesh%nC(vi)
      DO ci = 1, mesh%nconmax
        WRITE(UNIT = 1337, FMT = '(I6)', ADVANCE = 'NO') mesh%C(vi,ci)
      END DO
      WRITE(UNIT = 1337, FMT = '(I3)', ADVANCE = 'NO') mesh%niTri(vi)
      DO ci = 1, mesh%nconmax
        WRITE(UNIT = 1337, FMT = '(I6)', ADVANCE = 'NO') mesh%iTri(vi,ci)
      END DO
      WRITE(UNIT = 1337, FMT = '(I3)', ADVANCE = 'NO') mesh%edge_index(vi)
      WRITE(UNIT = 1337, FMT = '(A)') ''
    END DO
    WRITE(UNIT = 1337, FMT = '(A)')       ''
    
    ! Triangle data
    WRITE(UNIT = 1337, FMT = '(A)')       'Tri  TriC  Tri_edge_index'    
    DO ti = 1, mesh%nTri
      WRITE(UNIT = 1337, FMT = '(6I6,I3)') mesh%Tri(ti,1), mesh%Tri(ti,2), mesh%Tri(ti,3), mesh%TriC(ti,1), mesh%TriC(ti,2), mesh%TriC(ti,3), mesh%Tri_edge_index(ti)
    END DO
    
    ! Close the text file
    CLOSE(UNIT = 1337)
    
  END SUBROUTINE WriteMeshToTextFile
  SUBROUTINE CheckMesh( mesh)
    ! Check if the mesh data is self-consistent

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER                                       :: vi, ci, vc, ci2, vc2, iti, iti2, ti, n, v1, v2, v3, ti2, n2
    LOGICAL                                       :: FoundIt
    
    IF (.NOT. par%master) RETURN    
    
    ! == Array sizes
    ! =============================================================
    
!    IF (SIZE(mesh%V)              /= mesh%nV   * 2            ) WRITE(0,*) ' CheckMesh - ERROR: array V              has incorrect size!'
!    IF (SIZE(mesh%nC)             /= mesh%nV   * 1            ) WRITE(0,*) ' CheckMesh - ERROR: array nC             has incorrect size!'
!    IF (SIZE(mesh%C)              /= mesh%nV   * mesh%nconmax ) WRITE(0,*) ' CheckMesh - ERROR: array C              has incorrect size!'
!    IF (SIZE(mesh%niTri)          /= mesh%nV   * 1            ) WRITE(0,*) ' CheckMesh - ERROR: array niTri          has incorrect size!'
!    IF (SIZE(mesh%iTri)           /= mesh%nV   * mesh%nconmax ) WRITE(0,*) ' CheckMesh - ERROR: array iTri           has incorrect size!'
!    IF (SIZE(mesh%edge_index)     /= mesh%nV   * 1            ) WRITE(0,*) ' CheckMesh - ERROR: array edge_index     has incorrect size!'
!    
!    IF (SIZE(mesh%Tri)            /= mesh%nTri * 3            ) WRITE(0,*) ' CheckMesh - ERROR: array Tri            has incorrect size!'
!    IF (SIZE(mesh%Tricc)          /= mesh%nTri * 2            ) WRITE(0,*) ' CheckMesh - ERROR: array Tricc          has incorrect size!'
!    IF (SIZE(mesh%TriC)           /= mesh%nTri * 3            ) WRITE(0,*) ' CheckMesh - ERROR: array TriC           has incorrect size!'
!    IF (SIZE(mesh%Tri_edge_index) /= mesh%nTri * 1            ) WRITE(0,*) ' CheckMesh - ERROR: array Tri_edge_index has incorrect size!'
        
    ! == V
    ! =============================================================
    DO vi = 1, mesh%nV
      IF (mesh%V(vi,1) < mesh%xmin .OR. mesh%V(vi,1) > mesh%xmax .OR. mesh%V(vi,2) < mesh%ymin .OR. mesh%V(vi,2) > mesh%ymax) THEN
        WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' outside mesh domain!'
      END IF
    END DO
    
    ! == nC
    ! =============================================================
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC(vi)
        IF (mesh%C(vi,ci) == 0) WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has fewer connections than nC says!'
      END DO
      DO ci = mesh%nC(vi)+1, mesh%nconmax
        IF (mesh%C(vi,ci) > 0)  WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has more connections than nC says!'
      END DO
    END DO
    
    ! == C
    ! =============================================================
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC(vi)
        vc = mesh%C(vi,ci)
        FoundIt = .FALSE.
        DO ci2 = 1, mesh%nC(vc)
          vc2 = mesh%C(vc,ci2)
          IF (vc2==vi) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' is connected to ', vc, ', but not the other way round!'
      END DO
    END DO
    
    ! == niTri
    ! =============================================================
    DO vi = 1, mesh%nV
      DO iti = 1, mesh%niTri(vi)
        IF (mesh%iTri(vi,iti) == 0) WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has fewer iTriangles than niTri says!'
      END DO
      DO iti = mesh%niTri(vi)+1, mesh%nconmax
        IF (mesh%iTri(vi,iti) > 0)  WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has more iTriangles than nC says!'
      END DO
    END DO
    
    ! == iTri
    ! =============================================================
    DO vi = 1, mesh%nV
      DO iti = 1, mesh%niTri(vi)
        ti = mesh%iTri(vi,iti)
        FoundIt = .FALSE.
        DO n = 1, 3
          IF (mesh%Tri(ti,n)==vi) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' lists triangle ', ti, ' in iTri, but that triangle ', ti, ' doesnt contain vertex ', vi, '!'
      END DO
      
      IF (mesh%edge_index(vi) == 0) THEN
      
        DO ci = 1, mesh%nC(vi)
          vc = mesh%C(vi,ci)          
          n = 0
          DO iti = 1, mesh%niTri(vi)
            ti = mesh%iTri(vi,iti)
            DO iti2 = 1, mesh%niTri(vc)
              ti2 = mesh%iTri(vc,iti2)
              IF (ti==ti2) THEN
                n = n+1
                EXIT
              END IF
            END DO
          END DO
          IF (.NOT. (n==2)) WRITE(0,*) ' CheckMesh - ERROR: non-edge vertices ', vi, ' and ', vc, ' share ', n, ' triangles'        
        END DO
        
      ELSE
      
        DO ci = 1, mesh%nC(vi)
          vc = mesh%C(vi,ci)
          IF (mesh%edge_index(vc)==0) THEN
                
            n = 0
            DO iti = 1, mesh%niTri(vi)
              ti = mesh%iTri(vi,iti)
              DO iti2 = 1, mesh%niTri(vc)
                ti2 = mesh%iTri(vc,iti2)
                IF (ti==ti2) THEN
                  n = n+1
                  EXIT
                END IF
              END DO
            END DO
            IF (.NOT. (n==2)) WRITE(0,*) ' CheckMesh - ERROR: non-edge vertices ', vi, ' and ', vc, ' share ', n, ' triangles'
          
          ELSE
                
            n = 0
            DO iti = 1, mesh%niTri(vi)
              ti = mesh%iTri(vi,iti)
              DO iti2 = 1, mesh%niTri(vc)
                ti2 = mesh%iTri(vc,iti2)
                IF (ti==ti2) THEN
                  n = n+1
                  EXIT
                END IF
              END DO
            END DO
            IF ((mesh%edge_index(vi)==1 .AND. mesh%edge_index(vc)==7) .OR. &
                (mesh%edge_index(vi)==1 .AND. mesh%edge_index(vc)==3) .OR. &
                (mesh%edge_index(vi)==3 .AND. mesh%edge_index(vc)==1) .OR. &
                (mesh%edge_index(vi)==3 .AND. mesh%edge_index(vc)==5) .OR. &
                (mesh%edge_index(vi)==5 .AND. mesh%edge_index(vc)==3) .OR. &
                (mesh%edge_index(vi)==5 .AND. mesh%edge_index(vc)==7) .OR. &
                (mesh%edge_index(vi)==7 .AND. mesh%edge_index(vc)==5) .OR. &
                (mesh%edge_index(vi)==7 .AND. mesh%edge_index(vc)==1)) CYCLE
            IF (.NOT. (n==1)) WRITE(0,*) ' CheckMesh - ERROR: edge vertices ', vi, ' and ', vc, ' share ', n, ' triangles'
            
          END IF
        END DO
        
      END IF
    END DO
    
    ! == edge_index
    ! =============================================================
    DO vi = 1, mesh%nV
    
      IF (mesh%edge_index(vi) == 0) THEN
      
        IF (mesh%V(vi,1) <= mesh%xmin .OR. mesh%V(vi,1) >= mesh%xmax .OR. mesh%V(vi,2) <= mesh%ymin .OR. mesh%V(vi,2) >= mesh%ymax) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 0 but lies on or beyond the mesh domain boundary!'
        END IF
        
        ! First and last neighbours must be connected
        vc  = mesh%C(vi,1)
        vc2 = mesh%C(vi,mesh%nC(vi))
        
        FoundIt = .FALSE.
        DO ci = 1, mesh%nC(vc)
          IF (mesh%C(vc,ci)==vc2) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 0, but its first and last neighbours are not connected!'
        
      ELSEIF (mesh%edge_index(vi) == 1) THEN
      
        IF (ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 1 but does not lie on the N boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==8 .OR. mesh%edge_index(vc)==1)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 1 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==1 .OR. mesh%edge_index(vc)==2)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 1 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 1 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 1 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 2) THEN
      
        IF (.NOT. vi==3) WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' is listed as NE corner!'
      
        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 2 but does not lie on the NE corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==8 .OR. mesh%edge_index(vc)==1)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 2 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==3 .OR. mesh%edge_index(vc)==4)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 2 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 2 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 2 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 3) THEN
      
        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 3 but does not lie on the E boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==2 .OR. mesh%edge_index(vc)==3)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 3 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==3 .OR. mesh%edge_index(vc)==4)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 3 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 3 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 3 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 4) THEN
      
        IF (.NOT. vi==2) WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' is listed as SE corner!'
      
        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 4 but does not lie on the SE corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==2 .OR. mesh%edge_index(vc)==3)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 4 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==5 .OR. mesh%edge_index(vc)==6)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 4 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 4 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 4 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 5) THEN
      
        IF (ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 5 but does not lie on the S boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==4 .OR. mesh%edge_index(vc)==5)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 5 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==5 .OR. mesh%edge_index(vc)==6)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 5 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 5 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 5 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 6) THEN
      
        IF (.NOT. vi==1) WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' is listed as SW corner!'
      
        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 6 but does not lie on the SW corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==4 .OR. mesh%edge_index(vc)==5)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 6 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==7 .OR. mesh%edge_index(vc)==8)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 6 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 6 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 6 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 7) THEN
      
        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 7 but does not lie on the W boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==6 .OR. mesh%edge_index(vc)==7)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 7 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==7 .OR. mesh%edge_index(vc)==8)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 7 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 7 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 7 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 8) THEN
      
        IF (.NOT. vi==4) WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' is listed as NW corner!'
      
        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 8 but does not lie on the NW corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==6 .OR. mesh%edge_index(vc)==7)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 8 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==1 .OR. mesh%edge_index(vc)==2)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 8 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 8 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
          WRITE(0,*) ' CheckMesh - ERROR: vertex ', vi, ' has edge_index 8 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      END IF
      
    END DO
    
    ! == Tri
    ! =============================================================
    
    DO ti = 1, mesh%nTri
      DO n = 1, 3
        vi = mesh%Tri(ti,n)
        FoundIt = .FALSE.
        DO iti = 1, mesh%niTri(vi)
          IF (mesh%iTri(vi,iti) == ti) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' CheckMesh - ERROR: triangle ', ti, ' contains vertex ', vi, ', but that vertex doesnt list ti as an iTri!'
      END DO
      
      v1 = mesh%Tri(ti,1)
      v2 = mesh%Tri(ti,2)
      v3 = mesh%Tri(ti,3)
      
      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v1)
        vc = mesh%C(v1,ci)
        IF (vc==v2) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) WRITE(0,*) ' CheckMesh - ERROR: triangle ', ti, ' contains unconnected vertices ', v1, ' and ', v2, '!'
      
      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v1)
        vc = mesh%C(v1,ci)
        IF (vc==v3) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) WRITE(0,*) ' CheckMesh - ERROR: triangle ', ti, ' contains unconnected vertices ', v1, ' and ', v3, '!'
      
      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v2)
        vc = mesh%C(v2,ci)
        IF (vc==v3) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) WRITE(0,*) ' CheckMesh - ERROR: triangle ', ti, ' contains unconnected vertices ', v2, ' and ', v3, '!'
    END DO
    
    ! == TriC
    ! =============================================================
    
    DO ti = 1, mesh%nTri
      DO n = 1, 3
        ti2 = mesh%TriC(ti,n)
        IF (ti2 == 0) THEN
          IF (mesh%Tri_edge_index(ti) == 0) WRITE(0,*) ' CheckMesh - ERROR: non-edge triangle ', ti, ' misses a neighbour!'
          CYCLE
        END IF
        FoundIt = .FALSE.
        DO n2 = 1, 3
          IF (mesh%TriC(ti2,n2) == ti) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' CheckMesh - ERROR: triangle ', ti, ' is connected to ', ti2, ', but not the other way round!'
      END DO
    END DO
        
    
  END SUBROUTINE CheckMesh

END MODULE mesh_help_functions_module
