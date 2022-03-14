MODULE mesh_help_functions_module

  ! General help functions used in mesh creation and updating.

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
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
  USE data_types_module,               ONLY: type_mesh, type_model_region
  USE utilities_module,                ONLY: line_integral_mxydx, line_integral_xydy

  IMPLICIT NONE

  ! Interfaces to LAPACK, which are otherwise implicitly generated (taken from
  ! LAPACK source)
  !  *
  !  *  -- LAPACK routine (version 3.1) --
  !  *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  !  *     November 2006
  interface 
    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
    END SUBROUTINE
  end interface

  CONTAINS
    
! == Calculating some extra mesh data: Voronoi cell areas, connection widths,
!    triangle areas, resolution, lat/lon-coordinates
  SUBROUTINE find_Voronoi_cell_areas( mesh)
    ! Find the areas of the Voronoi cells of all the vertices    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_Voronoi_cell_areas'
    INTEGER                                       :: vi, nVor, n
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Vor
    REAL(dp)                                      :: Aerr
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ALLOCATE(Vor(mesh%nC_mem+2,2))

    mesh%A(mesh%vi1:mesh%vi2) = 0._dp
    DO vi = mesh%vi1, mesh%vi2

      CALL find_Voronoi_cell_vertices(mesh, vi, Vor, nVor)

      DO n = 2, nVor
        mesh%A(vi) = mesh%A(vi) + ABS(cross2( &
        [Vor(n  ,1)-mesh%V(vi,1), Vor(n  ,2)-mesh%V(vi,2)], &
        [Vor(n-1,1)-mesh%V(vi,1), Vor(n-1,2)-mesh%V(vi,2)] )) / 2._dp
      END DO
    END DO
    CALL sync

    DEALLOCATE(Vor)

    ! Check if everything went alright
    IF (par%master) THEN
      Aerr = ABS(1._dp - SUM(mesh%A ) / ((mesh%xmax-mesh%xmin)*(mesh%ymax-mesh%ymin))) / 100._dp
      IF (Aerr > 0.0001_dp) CALL warning('sum of Voronoi cell areas doesnt match square area of mesh! (error of {dp_01} %)', dp_01 = Aerr)
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE find_Voronoi_cell_areas 
  SUBROUTINE find_Voronoi_cell_geometric_centres( mesh)
    ! Find the geometric centres of the Voronoi cells of all the vertices    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_Voronoi_cell_geometric_centres'
    INTEGER                                       :: vi, nvi
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Vor
    INTEGER                                       :: nVor
    REAL(dp), DIMENSION(2)                        :: p, q
    REAL(dp)                                      :: LI_mxydx, LI_xydy
    REAL(dp)                                      :: LI_mxydx_seg, LI_xydy_seg
    
    ! Add routine to path
    CALL init_routine( routine_name)
        
    mesh%VorGC( mesh%vi1:mesh%vi2,:) = 0._dp
    
    ALLOCATE(Vor(mesh%nC_mem+2,2))
    
    DO vi = mesh%vi1, mesh%vi2
    
      CALL find_Voronoi_cell_vertices(mesh, vi, Vor, nVor)
      
      LI_mxydx = 0._dp
      LI_xydy  = 0._dp
      
      DO nvi = 2, nVor
        p = Vor( nvi-1,:)
        q = Vor( nvi,:)
        CALL line_integral_mxydx( p, q, mesh%tol_dist, LI_mxydx_seg)
        CALL line_integral_xydy(  p, q, mesh%tol_dist, LI_xydy_seg )
        LI_mxydx = LI_mxydx + LI_mxydx_seg
        LI_xydy  = LI_xydy  + LI_xydy_seg
      END DO 
      
      IF (mesh%edge_index( vi) > 0) THEN
      
        p = Vor( nVor,:)
        q = mesh%V( vi,:)
        CALL line_integral_mxydx( p, q, mesh%tol_dist, LI_mxydx_seg)
        CALL line_integral_xydy(  p, q, mesh%tol_dist, LI_xydy_seg)
        LI_mxydx = LI_mxydx + LI_mxydx_seg
        LI_xydy  = LI_xydy  + LI_xydy_seg
        
        p = mesh%V( vi,:)
        q = Vor( 1,:)
        CALL line_integral_mxydx( p, q, mesh%tol_dist, LI_mxydx_seg)
        CALL line_integral_xydy(  p, q, mesh%tol_dist, LI_xydy_seg)
        LI_mxydx = LI_mxydx + LI_mxydx_seg
        LI_xydy  = LI_xydy  + LI_xydy_seg
        
      END IF
       
      mesh%VorGC( vi,:) = [LI_mxydx / mesh%A( vi), LI_xydy / mesh%A( vi)]
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    DEALLOCATE(Vor)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE find_Voronoi_cell_geometric_centres
  SUBROUTINE find_connection_widths( mesh)
    ! Find the width of the line separating two connected vertices (equal to the distance
    ! between the circumcenters of their two shared triangles)    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_connection_widths'
    INTEGER                                       :: v1, nv2, v2, t1, t2, iti, ti, n
    LOGICAL                                       :: hasv2
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    mesh%Cw(mesh%vi2:mesh%vi2,:) = 0._dp
    CALL sync

    ! The way to go: for each vertex pair, find the two triangles that both
    ! are a part of. If there's only one, there's one Voronoi vertex and its
    ! projection on the map edge.
    DO v1 = mesh%vi1, mesh%vi2

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
          CALL crash('couldnt find a single shared triangle!')
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
            CALL crash('the only shared triangle isnt a Boundary triangle - this cannot be!')
          END IF
        END IF ! IF (t2>0) THEN

      END DO ! DO nv2 = 1, mesh%nC(v1)
    END DO ! DO v1 = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE find_connection_widths   
  SUBROUTINE find_triangle_areas( mesh)
    ! Find the areas of all the triangles    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_triangle_areas'
    INTEGER                                       :: ti
    REAL(dp), DIMENSION(2)                        :: pa, pb, pc
    
    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2
      pa = mesh%V( mesh%Tri( ti,1),:)
      pb = mesh%V( mesh%Tri( ti,2),:)
      pc = mesh%V( mesh%Tri( ti,3),:)
      CALL find_triangle_area( pa, pb, pc, mesh%TriA(ti))
    END DO ! DO ti = mesh%ti1, mesh%ti2
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE find_triangle_areas
  SUBROUTINE determine_mesh_resolution( mesh)
    ! Find the areas of all the triangles    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_mesh_resolution'
    INTEGER                                       :: vi, vj, ci
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    mesh%R(mesh%vi1:mesh%vi2) = mesh%xmax - mesh%xmin
    
    DO vi = mesh%vi1, mesh%vi2
      DO ci = 1, mesh%nC(vi)
        vj = mesh%C(vi,ci)
        mesh%R(vi) = MIN(mesh%R(vi), SQRT((mesh%V(vj,1)-mesh%V(vi,1))**2 + (mesh%V(vj,2)-mesh%V(vi,2))**2))
      END DO
    END DO
    CALL sync
      
    mesh%resolution_min = MINVAL(mesh%R)
    mesh%resolution_max = MAXVAL(mesh%R)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_mesh_resolution
  SUBROUTINE calc_triangle_geometric_centres( mesh)
    ! Find the geometric centres of all the triangles    
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_triangle_geometric_centres'
    INTEGER                                       :: ti
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    DO ti = mesh%ti1, mesh%ti2
      CALL update_triangle_geometric_center( mesh, ti)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE calc_triangle_geometric_centres
  
! == Finding the vertices of a vertex' Voronoi cell
  SUBROUTINE find_Voronoi_cell_vertices(        mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a vertex's Voronoi cell    
    
    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: vi
    REAL(dp), DIMENSION(:,:), INTENT(INOUT)       :: Vor
    INTEGER,                  INTENT(OUT)         :: nVor
    
    ! Local variables
    INTEGER                                       :: vvi

    IF (mesh%edge_index(vi)==0) THEN
      CALL find_Voronoi_cell_vertices_free(mesh, vi, Vor, nVor)
    ELSEIF (mesh%edge_index(vi)==2 .OR. mesh%edge_index(vi)==4 .OR. mesh%edge_index(vi)==6 .OR. mesh%edge_index(vi)==8) THEN
      CALL find_Voronoi_cell_vertices_corner(mesh, vi, Vor, nVor)
    ELSE
      CALL find_Voronoi_cell_vertices_edge(mesh, vi, Vor, nVor)
    END IF
    
    DO vvi = 1, nVor
      IF (Vor(vvi,1) < mesh%xmin - mesh%tol_dist .OR. &
          Vor(vvi,1) > mesh%xmax + mesh%tol_dist .OR. &
          Vor(vvi,2) < mesh%ymin - mesh%tol_dist .OR. &
          Vor(vvi,2) > mesh%ymax + mesh%tol_dist) THEN
        WRITE(0,*) '   find_Voronoi_cell_vertices - ERROR: found Voronoi cell vertex outside of mesh domain!'
        STOP
      END IF
      Vor(vvi,1) = MAX( MIN( Vor(vvi,1), mesh%xmax), mesh%xmin)
      Vor(vvi,2) = MAX( MIN( Vor(vvi,2), mesh%ymax), mesh%ymin)
    END DO

  END SUBROUTINE find_Voronoi_cell_vertices 
  SUBROUTINE find_Voronoi_cell_vertices_free(   mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a free vertex's Voronoi cell    
    
    IMPLICIT NONE

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
      CALL crop_circumcenter(mesh, ti, ti_clock, cc_clock)
      CALL crop_circumcenter(mesh, ti, ti_anti,  cc_anti)

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

  END SUBROUTINE find_Voronoi_cell_vertices_free
  SUBROUTINE find_Voronoi_cell_vertices_edge(   mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up an edge vertex's Voronoi cell    
    
    IMPLICIT NONE

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
        CALL crop_circumcenter(mesh, mesh%iTri(vi,1), mesh%iTri(vi,2), cc_cropped)
        nVor = nVor+1
        Vor(nVor,:) = cc_cropped
      END IF ! IF (iti == 1) THEN

      IF (iti > 1 .AND. iti < mesh%niTri(vi)) THEN
        ! Split the circumcenter
        CALL crop_circumcenter(mesh, mesh%iTri(vi,iti), mesh%iTri(vi,iti+1), cc_anti)
        CALL crop_circumcenter(mesh, mesh%iTri(vi,iti), mesh%iTri(vi,iti-1), cc_clock)
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
        CALL crop_circumcenter(mesh, mesh%iTri(vi,iti), mesh%iTri(vi,iti-1), cc_cropped)
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

  END SUBROUTINE find_Voronoi_cell_vertices_edge
  SUBROUTINE find_Voronoi_cell_vertices_corner( mesh, vi, Vor, nVor)
    ! Find the coordinates of the points making up a corner vertex's Voronoi cell    
    
    IMPLICIT NONE

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
          
      CALL find_Voronoi_cell_vertices_edge(mesh, vi, Vor, nVor)
      
      IF     (mesh%edge_index( vi) == 2) THEN
        ! Northeast corner
        nVor = nVor + 1
        Vor( nVor,:) = [mesh%xmax, mesh%ymax]
      ELSEIF (mesh%edge_index( vi) == 4) THEN
        ! Southeast corner
        nVor = nVor + 1
        Vor( nVor,:) = [mesh%xmax, mesh%ymin]
      ELSEIF (mesh%edge_index( vi) == 6) THEN
        ! Southwest corner
        nVor = nVor + 1
        Vor( nVor,:) = [mesh%xmin, mesh%ymin]
      ELSEIF (mesh%edge_index( vi) == 8) THEN
        ! Northwest corner
        nVor = nVor + 1
        Vor( nVor,:) = [mesh%xmin, mesh%ymax]
      END IF
      
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

  END SUBROUTINE find_Voronoi_cell_vertices_corner
  SUBROUTINE crop_circumcenter( mesh, t1, t2, ccc)
    ! Crop the circumcenter of triangle t1 in the direction of t2    
    
    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: t1, t2
    REAL(dp), DIMENSION(2),   INTENT(OUT)         :: ccc

    REAL(dp)                                      :: la, lb, lc, le, lf, lg
    REAL(dp), DIMENSION(2)                        :: p, q

    ccc  = mesh%Tricc(t1,:)

    IF     (mesh%Tri_edge_index(t1)==1 .AND. mesh%Tricc(t1,2)>mesh%ymax) THEN
      ! North boundary triangle
      CALL line_from_points([mesh%xmin, mesh%ymax], [mesh%xmax, mesh%ymax], la, lb, lc)
    ELSEIF (mesh%Tri_edge_index(t1)==3 .AND. mesh%Tricc(t1,1)>mesh%xmax) THEN
      ! East boundary triangle
      CALL line_from_points([mesh%xmax, mesh%ymax], [mesh%xmax, mesh%ymin], la, lb, lc)
    ELSEIF (mesh%Tri_edge_index(t1)==5 .AND. mesh%Tricc(t1,2)<mesh%ymin) THEN
      ! South boundary triangle
      CALL line_from_points([mesh%xmin, mesh%ymin], [mesh%xmax, mesh%ymin], la, lb, lc)
    ELSEIF (mesh%Tri_edge_index(t1)==7 .AND. mesh%Tricc(t1,1)<mesh%xmin) THEN
      ! West boundary triangle
      CALL line_from_points([mesh%xmin, mesh%ymax], [mesh%xmin, mesh%ymin], la, lb, lc)
    ELSE
      RETURN
    END IF

   p = mesh%Tricc(t1,:)
   q = mesh%Tricc(t2,:)
   CALL line_from_points(p, q, le, lf, lg)
   CALL line_line_intersection(la, lb, lc, le, lf, lg, ccc);

  END SUBROUTINE crop_circumcenter
  SUBROUTINE find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
    ! Return the endpoints of the shared Voronoi cell boundary represented by edge aci
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: aci
    REAL(dp), DIMENSION(2),   INTENT(OUT)         :: cc1, cc2
    
    ! Local variables
    INTEGER                                       :: til,tir
    
    til = mesh%Aci( aci,5)
    tir = mesh%Aci( aci,6)
  
    IF (mesh%edge_index_Ac( aci) > 0) THEN
      ! Boundary segments have only one adjacent triangle
      
      IF (til > 0) THEN
        cc1 = mesh%Tricc( til,:)
      ELSE
        cc1 = mesh%Tricc( tir,:)
      END IF
      IF     (mesh%edge_index_Ac( aci) == 1) THEN
        ! North
        cc2 = [cc1(1), mesh%ymax]
      ELSEIF (mesh%edge_index_Ac( aci) == 3) THEN
        ! East
        cc2 = [mesh%xmax, cc1(2)]
      ELSEIF (mesh%edge_index_Ac( aci) == 5) THEN
        ! South
        cc2 = [cc1(1), mesh%ymin]
      ELSEIF (mesh%edge_index_Ac( aci) == 7) THEN
        ! West
        cc2 = [mesh%xmin, cc1(2)]
      END IF
      
    ELSE ! IF (mesh%edge_index_Ac( aci) > 0) THEN
      
      cc1 = mesh%Tricc( til,:)
      cc2 = mesh%Tricc( tir,:)
      
    END IF ! IF (mesh%edge_index_Ac( aci) > 0) THEN
    
  END SUBROUTINE find_shared_Voronoi_boundary
  
! == The oblique stereographic projection
  SUBROUTINE get_lat_lon_coordinates( mesh)
   ! Use the inverse stereographic projection for the mesh model region to calculate
   ! lat/lon coordinates for all the vertices    
    
    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(INOUT)       :: mesh
    
    INTEGER                                       :: vi
    
    ! Calculate lat and lon directly from X and Y using inverse projection 
    DO vi = mesh%vi1, mesh%vi2
      CALL inverse_oblique_sg_projection(mesh%V(vi,1), mesh%V(vi,2), mesh%lambda_M, mesh%phi_M, mesh%alpha_stereo, mesh%lon(vi), mesh%lat(vi))
    END DO
    CALL sync
    
  END SUBROUTINE get_lat_lon_coordinates  
  SUBROUTINE oblique_sg_projection( lambda, phi, lambda_M_deg, phi_M_deg, alpha_deg, x_IM_P_prime, y_IM_P_prime)

    ! lambda = lon, phi = lat

    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates which coincide with the GCM grid points to the rectangular IM coordinate
    ! system, with coordinates (x,y).
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    
    USE parameters_module, ONLY: pi, earth_radius
    
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
    
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = (pi / 180._dp) * phi
    lambda_P = (pi / 180._dp) * lambda

    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = ((1._dp + COS(alpha)) / (1._dp + COS(phi_P) * COS(phi_M) * COS(lambda_P - lambda_M) + SIN(phi_P) * SIN(phi_M))) / (pi / 180._dp)

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x_IM_P_prime =  earth_radius * (COS(phi_P) * SIN(lambda_P - lambda_M)) * t_P_prime
    y_IM_P_prime =  earth_radius * (SIN(phi_P) * COS(phi_M) - (COS(phi_P) * SIN(phi_M)) * COS(lambda_P - lambda_M)) * t_P_prime   
    
    
  END SUBROUTINE oblique_sg_projection
  SUBROUTINE inverse_oblique_sg_projection( x_IM_P_prime, y_IM_P_prime, lambda_M_deg, phi_M_deg, alpha_deg, lambda_P, phi_P)
    ! This subroutine projects with an inverse oblique stereographic projection the
    ! (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, alpha, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    
    USE parameters_module, ONLY: pi, earth_radius
    
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
    
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! See equations (2.14-2.16) or equations (B.21-B.23) in Reerink et al. (2010):
    x_3D_P_prime = earth_radius * COS(alpha) * COS(lambda_M) * COS(phi_M) - SIN(lambda_M) * (x_IM_P_prime*(pi / 180._dp)) - COS(lambda_M) * SIN(phi_M) * (y_IM_P_prime*(pi / 180._dp))
    y_3D_P_prime = earth_radius * COS(alpha) * SIN(lambda_M) * COS(phi_M) + COS(lambda_M) * (x_IM_P_prime*(pi / 180._dp)) - SIN(lambda_M) * SIN(phi_M) * (y_IM_P_prime*(pi / 180._dp))
    z_3D_P_prime = earth_radius * COS(alpha) *                 SIN(phi_M)                                          +                   COS(phi_M) * (y_IM_P_prime*(pi / 180._dp))

    ! See equation (2.13) or equation (B.20) in Reerink et al. (2010):
    a = COS(lambda_M) * COS(phi_M) * x_3D_P_prime  +  SIN(lambda_M) * COS(phi_M) * y_3D_P_prime  +  SIN(phi_M) * z_3D_P_prime

    ! See equation (2.12) or equation (B.19) in Reerink et al. (2010):
    t_P = (2._dp * earth_radius**2 + 2._dp * earth_radius * a) / (earth_radius**2 + 2._dp * earth_radius * a + x_3D_P_prime**2 + y_3D_P_prime**2 + z_3D_P_prime**2)

    ! See equations (2.9-2.11) or equations (B.16-B.18) in Reerink et al. (2010):
    x_3D_P =  earth_radius * COS(lambda_M) * COS(phi_M) * (t_P - 1._dp) + x_3D_P_prime * t_P
    y_3D_P =  earth_radius * SIN(lambda_M) * COS(phi_M) * (t_P - 1._dp) + y_3D_P_prime * t_P
    z_3D_P =  earth_radius *                   SIN(phi_M) * (t_P - 1._dp) + z_3D_P_prime * t_P

    ! See equation (2.7) or equation (B.24) in Reerink et al. (2010):
    IF(x_3D_P <  0._dp                      ) THEN
      lambda_P = 180._dp + (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P >= 0._dp) THEN
      lambda_P =           (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P <  0._dp) THEN
      lambda_P = 360._dp + (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P >  0._dp) THEN
      lambda_P =  90._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P <  0._dp) THEN
      lambda_P = 270._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P == 0._dp) THEN
      lambda_P =   0._dp
    END IF

    ! See equation (2.8) or equation (B.25) in Reerink et al. (2010):
    IF(x_3D_P /= 0._dp .OR. y_3D_P /= 0._dp) THEN
      phi_P = (180._dp / pi) * ATAN(z_3D_P / sqrt(x_3D_P**2 + y_3D_P**2))
    ELSE IF(z_3D_P >  0._dp) THEN
      phi_P =   90._dp
    ELSE IF(z_3D_P <  0._dp) THEN
      phi_P =  -90._dp
    END IF
  END SUBROUTINE inverse_oblique_sg_projection
     
! == Some basic geometrical operations
  FUNCTION   is_in_triangle( pa, pb, pc, p) RESULT(isso)
    ! Check if the point p lies inside the triangle abc, or within distance tol_dist of its edges    
    
    IMPLICIT NONE

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
   
  END FUNCTION is_in_triangle  
  FUNCTION   is_boundary_segment( mesh, v1, v2) RESULT(isso)
   ! Determine whether or not the line between two vertices is an Edge segment    
    
    IMPLICIT NONE

   TYPE(type_mesh),          INTENT(IN)          :: mesh
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

  END FUNCTION is_boundary_segment
  SUBROUTINE is_encroached_upon( mesh, v1, v2, isso)
    ! TRUE if the line between v1 and v2 is encroached upon by any other
    ! vertex    
    
    IMPLICIT NONE

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
  END SUBROUTINE is_encroached_upon
  SUBROUTINE encroaches_upon( mesh, x, y, isso, va, vb)
    ! TRUE if the point [x,y] encroaches upon any boundary segment [va,vb]    
    
    IMPLICIT NONE

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
      IF (is_boundary_segment(mesh,p1,p2)) THEN
        IF (NORM2([x,y] - (mesh%V(p1,:) + mesh%V(p2,:))/2._dp) < NORM2(mesh%V(p1,:) - mesh%V(p2,:))/2._dp) THEN
          va = p1
          vb = p2
          isso = .TRUE.
          RETURN
        END IF
      END IF

      IF (is_boundary_segment(mesh,p2,p3)) THEN
        IF (NORM2([x,y] - (mesh%V(p2,:) + mesh%V(p3,:))/2._dp) < NORM2(mesh%V(p2,:) - mesh%V(p3,:))/2._dp) THEN
          va = p2
          vb = p3
          isso = .TRUE.
          RETURN
        END IF
      END IF

      IF (is_boundary_segment(mesh,p3,p1)) THEN
        IF (NORM2([x,y] - (mesh%V(p3,:) + mesh%V(p1,:))/2._dp) < NORM2(mesh%V(p3,:) - mesh%V(p1,:))/2._dp) THEN
          va = p3
          vb = p1
          isso = .TRUE.
          RETURN
        END IF
      END IF
    END DO ! DO ti = 1, mesh%nTri
  END SUBROUTINE encroaches_upon
  FUNCTION is_walltowall( mesh, ti) RESULT(isso)
   ! Determine whether or not a triangle is "wall to wall"
   ! (i.e. contains vertices lying on opposite domain boundaries)
    
    IMPLICIT NONE

   TYPE(type_mesh),          INTENT(IN)          :: mesh
   INTEGER,                  INTENT(IN)          :: ti
   LOGICAL                                       :: isso
   
   LOGICAL                                       :: has_north, has_south, has_east, has_west
   INTEGER                                       :: n, vi
   
   has_north = .FALSE.
   has_south = .FALSE.
   has_east  = .FALSE.
   has_west  = .FALSE.
   
   DO n = 1, 3
     vi = mesh%Tri( ti,n)
     IF (mesh%edge_index( vi) == 1) has_north = .TRUE.
     IF (mesh%edge_index( vi) == 3) has_east  = .TRUE.
     IF (mesh%edge_index( vi) == 5) has_south = .TRUE.
     IF (mesh%edge_index( vi) == 7) has_west  = .TRUE.
   END DO
   
   isso = .FALSE.
   IF (has_north .AND. has_south) isso = .TRUE.
   IF (has_west  .AND. has_east ) isso = .TRUE.
   
  END FUNCTION is_walltowall
  SUBROUTINE update_triangle_circumcenter( mesh, ti)
    ! Calculate the circumcenter of mesh triangle ti    
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti
    
    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: v1, v2, v3, cc
    
    v1 = mesh%V( mesh%Tri( ti,1),:)
    v2 = mesh%V( mesh%Tri( ti,2),:)
    v3 = mesh%V( mesh%Tri( ti,3),:)     
    CALL find_circumcenter( v1, v2, v3, cc)
    
    ! If find_circumcenter yields infinity, it's because p and q have the
    ! same y-coordinate. Rearrange vertices in triangle matrix (maintaining
    ! counter-clockwise orientation) and try again.
    
    IF (cc(1) > (mesh%xmax-mesh%xmin)*100000._dp .OR. cc(2) > (mesh%ymax-mesh%ymin)*10000._dp) THEN
      mesh%Tri(ti,:) = [mesh%Tri(ti,2), mesh%Tri(ti,3), mesh%Tri(ti,1)]
      v1 = mesh%V( mesh%Tri( ti,1),:)
      v2 = mesh%V( mesh%Tri( ti,2),:)
      v3 = mesh%V( mesh%Tri( ti,3),:) 
      CALL find_circumcenter( v1, v2, v3, cc)
    END IF
    
    IF (cc(1) > (mesh%xmax-mesh%xmin)*100000._dp .OR. cc(2) > (mesh%ymax-mesh%ymin)*100000._dp) THEN
      mesh%Tri(ti,:) = [mesh%Tri(ti,2), mesh%Tri(ti,3), mesh%Tri(ti,1)]
      v1 = mesh%V( mesh%Tri( ti,1),:)
      v2 = mesh%V( mesh%Tri( ti,2),:)
      v3 = mesh%V( mesh%Tri( ti,3),:) 
      CALL find_circumcenter( v1, v2, v3, cc)
    END IF
    
    IF (cc(1) > (mesh%xmax-mesh%xmin)*100000._dp .OR. cc(2) > (mesh%ymax-mesh%ymin)*100000._dp) THEN
      WRITE(0,*) '  update_triangle_circumcenter - ERROR: triangle  doesn''t yield a valid circumcenter!'
    END IF
    
    mesh%Tricc(ti,:) = cc
    
  END SUBROUTINE update_triangle_circumcenter
  SUBROUTINE update_triangle_geometric_center( mesh, ti)
    ! Calculate the geometric centre of mesh triangle ti    
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti
    
    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: v1, v2, v3
    
    v1 = mesh%V( mesh%Tri( ti,1),:)
    v2 = mesh%V( mesh%Tri( ti,2),:)
    v3 = mesh%V( mesh%Tri( ti,3),:) 
    
    mesh%TriGC( ti,:) = (v1 + v2 + v3) / 3._dp
    
  END SUBROUTINE update_triangle_geometric_center
  
! == Routines for merging meshes created by parallel processes
  SUBROUTINE merge_vertices( mesh, nVl, nVr, nTril, nTrir, T, nT, vil, vir, orientation)
    ! Merge overlapping vertices vil and vir    
    
    IMPLICIT NONE
  
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
    
  END SUBROUTINE merge_vertices
  SUBROUTINE switch_vertices( mesh, vi1, vi2)
    ! Switch vertices vi1 and vi2    
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: vi1, vi2
    INTEGER, DIMENSION(:,:), ALLOCATABLE          :: ctovi1, ctovi2
    INTEGER                                       :: ci, vc, ti, iti, n, ci2
    
    ALLOCATE( ctovi1( mesh%nC_mem, 2))
    ALLOCATE( ctovi2( mesh%nC_mem, 2))
    
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
    
  END SUBROUTINE switch_vertices
  SUBROUTINE redo_Tri_edge_indices( mesh)
    ! Called by MergeSubmeshes, ran by a single processor for that processor's submesh
    ! Be careful not to include any sync statements here!    
    
    IMPLICIT NONE
  
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
        
  END SUBROUTINE redo_Tri_edge_indices
  
! == Routines for partitioning a domain or list, either regularly or load-balanced
  SUBROUTINE partition_domain_regular( xmin, xmax, i, n, xl, xr)
    ! Given the domain [xmin,xmax], partition it into n equally-sized parts, and return the range [xl,xr] for the i-th part
    
    ! In/output variables:
    REAL(dp),                   INTENT(IN)        :: xmin, xmax
    INTEGER,                    INTENT(IN)        :: i, n
    REAL(dp),                   INTENT(OUT)       :: xl, xr
    
    xl = xmin + REAL(i  ) * (xmax - xmin) / REAL(n)
    xr = xmin + REAL(i+1) * (xmax - xmin) / REAL(n)
    
  END SUBROUTINE partition_domain_regular
  SUBROUTINE partition_domain_x_balanced( mesh, ymin, ymax, i, n, xmin, xmax)
    ! Given a mesh, divide the region between ymin and ymax into n parts with
    ! (approximately) equal numbers of vertices. Return the [xmin,xmax] limits of the i-th part.
  
    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp),                   INTENT(IN)        :: ymin, ymax
    INTEGER,                    INTENT(IN)        :: i, n
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
    vbins = CEILING( REAL(mesh%nV) * 0.01_dp/ REAL(n))
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
    
  END SUBROUTINE partition_domain_x_balanced
  SUBROUTINE partition_domain_y_balanced( mesh, xmin, xmax, i, n, ymin, ymax)
    ! Given a mesh, divide the region between xmin and xmax into n parts with
    ! (approximately) equal numbers of vertices. Return the [ymin,ymax] limits of the i-th part.
  
    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp),                   INTENT(IN)        :: xmin, xmax
    INTEGER,                    INTENT(IN)        :: i, n
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
    vbins = CEILING( REAL(mesh%nV) * 0.01_dp/ REAL(n))
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
    
  END SUBROUTINE partition_domain_y_balanced
  
! == Subroutines for creating Points of Interest (POI)
  SUBROUTINE find_POI_xy_coordinates( mesh)
    ! Use the stereographic projection parameters for the relevant model region    
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: POIi
    REAL(dp)                                      :: lat, lon, x, y
    
    IF (mesh%region_name == 'NAM') THEN
      DO POIi = 1, mesh%nPOI
        mesh%POI_coordinates(POIi,1) = C%POI_NAM_coordinates((POIi*2)-1)
        mesh%POI_coordinates(POIi,2) = C%POI_NAM_coordinates( POIi*2   )
        mesh%POI_resolutions(POIi  ) = C%POI_NAM_resolutions( POIi     )
      END DO
    ELSEIF (mesh%region_name == 'EAS') THEN
      DO POIi = 1, mesh%nPOI
        mesh%POI_coordinates(POIi,1) = C%POI_EAS_coordinates((POIi*2)-1)
        mesh%POI_coordinates(POIi,2) = C%POI_EAS_coordinates( POIi*2   )
        mesh%POI_resolutions(POIi  ) = C%POI_EAS_resolutions( POIi     )
      END DO
    ELSEIF (mesh%region_name == 'GRL') THEN
      DO POIi = 1, mesh%nPOI
        mesh%POI_coordinates(POIi,1) = C%POI_GRL_coordinates((POIi*2)-1)
        mesh%POI_coordinates(POIi,2) = C%POI_GRL_coordinates( POIi*2   )
        mesh%POI_resolutions(POIi  ) = C%POI_GRL_resolutions( POIi     )
      END DO
    ELSEIF (mesh%region_name == 'ANT') THEN
      DO POIi = 1, mesh%nPOI
        mesh%POI_coordinates(POIi,1) = C%POI_ANT_coordinates((POIi*2)-1)
        mesh%POI_coordinates(POIi,2) = C%POI_ANT_coordinates( POIi*2   )
        mesh%POI_resolutions(POIi  ) = C%POI_ANT_resolutions( POIi     )
      END DO
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR - region name unknown!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    DO POIi = 1, mesh%nPOI
    
      lat = mesh%POI_coordinates(POIi,1)
      lon = mesh%POI_coordinates(POIi,2)
      
      CALL oblique_sg_projection( lon, lat, mesh%lambda_M, mesh%phi_M, mesh%alpha_stereo, x, y)
    
      mesh%POI_XY_coordinates(POIi,1) = x
      mesh%POI_XY_coordinates(POIi,2) = y
    
    END DO  
  
  END SUBROUTINE find_POI_xy_coordinates
  SUBROUTINE find_POI_vertices_and_weights( mesh)
  ! For each POI in this region, find the indices of the three vertices spanning the triangle
  ! containing it, and find their respective weights for a trilinear interpolation.
  ! Called by all, executed by the master.    
    
    IMPLICIT NONE
  
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
        CALL find_containing_triangle( mesh, p, ti)
        
        ! The three vertices spanning this triangle
        vi_a = mesh%Tri( ti,1)
        vi_b = mesh%Tri( ti,2)
        vi_c = mesh%Tri( ti,3)
        
        pa   = mesh%V( vi_a,:)
        pb   = mesh%V( vi_b,:)
        pc   = mesh%V( vi_c,:)
        
        ! Calculate their interpolation weights
        CALL find_triangle_area( pa, pb, pc, Atot)
        CALL find_triangle_area( pb, pc, p , Aa  )
        CALL find_triangle_area( pc, pa, p , Ab  )
        CALL find_triangle_area( pa, pb, p , Ac  )
        
        ! Save the indices and weights
        mesh%POI_vi( POIi,:) = [vi_a,    vi_b,    vi_c   ]
        MESH%POI_w(  POIi,:) = [Aa/Atot, Ab/Atot, Ac/Atot]
        
      END DO ! DO POIi = 1, mesh%nPOI
    
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE find_POI_vertices_and_weights
  
! == Some basic search operations on a mesh
  SUBROUTINE find_containing_triangle( mesh, p, t)
   ! Start at initial guess t, search outward from there using a
   ! flood-fill algorithm.    
    
    IMPLICIT NONE

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
      WRITE(0,*) 'find_containing_triangle - ERROR: point lies outside mesh domain!'
      STOP
    END IF

    ! See if the initial guess is correct.
    q = mesh%V(mesh%Tri(t,1),:)
    r = mesh%V(mesh%Tri(t,2),:)
    s = mesh%V(mesh%Tri(t,3),:)
    IF (is_in_triangle(q, r, s, p)) RETURN

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
    IF (is_in_triangle(q, r, s, p)) RETURN
    IF (lies_on_line_segment( q, r, p, mesh%tol_dist)) RETURN
    IF (lies_on_line_segment( r, s, p, mesh%tol_dist)) RETURN
    IF (lies_on_line_segment( s, q, p, mesh%tol_dist)) RETURN

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
        IF (is_in_triangle(q, r, s, p)) THEN

          ! Found it!
          t = ti

          RETURN

        ELSE ! if (is_in_triangle(ti,p))
          ! Did not find it. And add this triangle's non-checked neighbours to the new stack.

          DO n2 = 1, 3
            tin = mesh%TriC(ti,n2)
            IF (tin==0)              CYCLE ! This neighbour doesn't exist.
            IF (mesh%TriMap(tin)==1) CYCLE ! This neighbour has already been checked or is already in the stack.
            mesh%TriStackN2 = mesh%TriStackN2 + 1
            mesh%TriStack2( mesh%TriStackN2) = tin
            mesh%TriMap(tin) = 1
          END DO

        END IF ! IF (is_in_triangle(q, r, s, p, tol)) THEN
      END DO ! DO n = 1, mesh%triStackN1

      ! Cycle stacks.
      mesh%TriStack1  = mesh%TriStack2
      mesh%TriStackN1 = mesh%TriStackN2

      ! If no more non-checked neighbours could be found, terminate and throw an error.
      IF (mesh%TriStackN2==0) THEN
        WRITE(0,*) 'find_containing_triangle - ERROR: couldnt find triangle containing this point!'
        STOP
      END IF

    END DO ! DO WHILE (.NOT. FoundIt)


  END SUBROUTINE find_containing_triangle
  SUBROUTINE find_containing_vertex( mesh, p, vi)
    ! Find the vertex whose Voronoi cell contains the point p, using a "linear search"
    ! Start at initial guess vi. Check all neighbours of vi, find the one
    ! closest to p, select that one as the new vi. Repeat until all neighbours
    ! of vi are further away from p than vi itself.    
    
    IMPLICIT NONE
    
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
    WRITE(0,*) ' find_containing_vertex - ERROR: couldnt find closest vertex!'
    STOP

  END SUBROUTINE find_containing_vertex
  FUNCTION is_in_Voronoi_cell( mesh, p, vi) RESULT(isso)
    ! If the point p lies closer to vertex vi than to any other vertex, then
    ! by definition it lies inside the Voronoi cell of vertex vi. 
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(  2), INTENT(IN)          :: p
    INTEGER,                  INTENT(IN)          :: vi
    
    ! Local variables:
    LOGICAL                                       :: isso
    REAL(dp)                                      :: dist_vi
    INTEGER                                       :: vvi, vj
    
    isso = .TRUE.
    
    dist_vi = NORM2( mesh%V( vi,:) - p)
    
    DO vvi = 1, mesh%nC( vi)
      vj = mesh%C( vi,vvi)
      IF (NORM2( mesh%V( vj,:) - p) < dist_vi) THEN
        isso = .FALSE.
        RETURN
      END IF
    END DO

  END FUNCTION is_in_Voronoi_cell
  FUNCTION lies_on_line_segment( pa, pb, pc, tol_dist) RESULT(isso)
    ! Test if the point pc lies on the line pb-pc, or within tol_dist of it    
    
    IMPLICIT NONE
    
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

  END FUNCTION lies_on_line_segment
  SUBROUTINE segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
    ! Find out if the line segments [pq] and [rs] intersect. If so, return
    ! the coordinates of the point of intersection    
    
    IMPLICIT NONE
    
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

    ! If pq and rs are colinear, define them as not intersecting
    IF ( ABS( cross2( [q(1)-p(1), q(2)-p(2)], [s(1)-r(1), s(2)-r(2)] )) < tol_dist) THEN
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
    
  END SUBROUTINE segment_intersection
    
! == Some very basic geometry, mainly used for determining triangle circumcenters
  SUBROUTINE line_from_points( p, q, la, lb, lc)    
    
    IMPLICIT NONE

    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p, q
    REAL(dp),                   INTENT(OUT)       :: la, lb, lc

    ! Line PQ is represented as la*x + lb*y = lc
    la = q(2) - p(2)
    lb = p(1) - q(1)
    lc = la*(p(1))+ lb*(p(2));

  END SUBROUTINE line_from_points  
  SUBROUTINE perpendicular_bisector_from_line( p, q, la1, lb1, la2, lb2, lc2)    
    
    IMPLICIT NONE

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

  END SUBROUTINE perpendicular_bisector_from_line  
  SUBROUTINE line_line_intersection( la1, lb1, lc1, la2, lb2, lc2, llis)
    ! Find the intersection llis of the lines la1*x+lb1*y=lc1 and la2*x+lb2*y=lc2    
    
    IMPLICIT NONE

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
    
  END SUBROUTINE line_line_intersection  
  SUBROUTINE find_circumcenter( p, q, r, cc)
    ! Find the circumcenter cc of the triangle pqr    
    
    IMPLICIT NONE

    ! Some basic vector operations
    ! Find the circumcenter cc of the triangle pqr
    ! If pqr are colinear, returns [1e30,1e30]

    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p, q, r
    REAL(dp), DIMENSION(2),     INTENT(OUT)       :: cc
    REAL(dp)                                      :: la1,lb1,lc1,le1,lf1,lg1
    REAL(dp)                                      :: la2,lb2,lc2,le2,lf2,lg2

    cc = [0._dp, 0._dp]

    ! Line PQ is represented as ax + by = c, Line QR is represented as ex + fy = g
    CALL line_from_points( p, q, la1, lb1, lc1)
    CALL line_from_points( q, r, le1, lf1, lg1)

    ! Converting lines PQ and QR to perpendicular
    ! bisectors. After this, L = ax + by = c
    ! M = ex + fy = g
    CALL perpendicular_bisector_from_line( p, q, la1, lb1, la2, lb2, lc2)
    CALL perpendicular_bisector_from_line( q, r, le1, lf1, le2, lf2, lg2)

    ! The point of intersection of L and M gives
    ! the circumcenter
    CALL line_line_intersection( la2, lb2, lc2, le2, lf2, lg2, cc)

  END SUBROUTINE find_circumcenter
  FUNCTION cross2( a,b) RESULT(z)
    ! Vector product z between 2-dimensional vectors a and b    
    
    IMPLICIT NONE
    
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: a, b
    REAL(dp)                                      :: z

    z = (a(1)*b(2)) - (a(2)*b(1))

  END FUNCTION cross2
  SUBROUTINE find_triangle_area( pq, pr, ps, TriA)
    ! Find the area of the triangle [pq,pr,ps]    
    
    IMPLICIT NONE
    
    REAL(dp), DIMENSION(2), INTENT(IN)  :: pq, pr, ps
    REAL(dp),               INTENT(OUT) :: TriA
    
    TriA = ABS( cross2( [pr(1)-pq(1), pr(2)-pq(2)], [ps(1)-pq(1), ps(2)-pq(2)] )) / 2._dp
    
  END SUBROUTINE find_triangle_area
      
! == Help routines for use in mesh updating
  SUBROUTINE mean_cart_over_Voronoi_cell_dp(  mesh, d, x, y, nx, ny, vi, v)
    ! Find the mean value of the data field d, given on the cartesian
    ! grid x,y, over the Voronoi cell of vertex vi    
    
    IMPLICIT NONE

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

    ALLOCATE(Vor(mesh%nC_mem+2,2))

    CALL find_Voronoi_cell_vertices(mesh, vi, Vor, nVor)

    DO n = 2, nVor
      CALL sum_cart_over_triangle_dp( mesh%V(vi,:), Vor(n-1,:), Vor(n,:), d, x, y, nx, ny, trisumel, trinel)
      sumel = sumel + trisumel
      nel   = nel   + trinel
    END DO

    IF (nel>4) THEN
      v = sumel / nel
    ELSE
      ! Too few elements for a proper mean - just do bicubic interpolation to
      ! the vertex location
      CALL cart_bilinear_dp( d, x, y, nx, ny, mesh%V(vi,:), v)
    END IF

    DEALLOCATE(Vor)

  END SUBROUTINE mean_cart_over_Voronoi_cell_dp
  SUBROUTINE mean_cart_over_Voronoi_cell_int( mesh, d, x, y, nx, ny, vi, v)
    ! Find the mean value of the data field d, given on the cartesian
    ! grid x,y, over the Voronoi cell of vertex vi    
    
    IMPLICIT NONE

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

    ALLOCATE(Vor(mesh%nC_mem+2,2))

    CALL find_Voronoi_cell_vertices(mesh, vi, Vor, nVor)

    DO n = 2, nVor
      CALL sum_cart_over_triangle_int( mesh%V(vi,:), Vor(n-1,:), Vor(n,:), d, x, y, nx, ny, trisumel, trinel)
      sumel = sumel + REAL(trisumel)
      nel   = nel   + trinel
    END DO

    IF (nel>4) THEN
      v = sumel / REAL(nel)
    ELSE
      ! Too few elements for a proper mean - just do bicubic interpolation to
      ! the vertex location
      CALL cart_bilinear_int( d, x, y, nx, ny, mesh%V(vi,:), v)
    END IF

    DEALLOCATE(Vor)

  END SUBROUTINE mean_cart_over_Voronoi_cell_int  
  SUBROUTINE sum_cart_over_triangle_dp(  p1, p2, p3, d, x, y, nx, ny, trisumel, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear    
    
    IMPLICIT NONE

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
     IF (is_in_triangle(p1,p2,p3,[x(i), y(j)])) THEN
       trisumel = trisumel + d(i,j)
       trinel   = trinel   + 1
     END IF
    END DO
    END DO

  END SUBROUTINE sum_cart_over_triangle_dp
  SUBROUTINE sum_cart_over_triangle_int( p1, p2, p3, d, x, y, nx, ny, trisumel, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear    
    
    IMPLICIT NONE

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
     IF (is_in_triangle(p1,p2,p3,[x(i), y(j)])) THEN
       trisumel = trisumel + d(i,j)
       trinel   = trinel   + 1
     END IF
    END DO
    END DO

  END SUBROUTINE sum_cart_over_triangle_int
  SUBROUTINE max_cart_over_triangle_dp(  p1, p2, p3, d, x, y, nx, ny, vmax, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear    
    
    IMPLICIT NONE

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
     IF (is_in_triangle(p1,p2,p3,[x(i), y(j)])) THEN
      IF (d(i,j) > vmax) THEN
       vmax   = d(i,j)
       trinel = 1
      END IF
     END IF
    END DO
    END DO

  END SUBROUTINE max_cart_over_triangle_dp
  SUBROUTINE max_cart_over_triangle_int( p1, p2, p3, d, x, y, nx, ny, vmax, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear    
    
    IMPLICIT NONE

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
     IF (is_in_triangle(p1,p2,p3,[x(i), y(j)])) THEN
      IF (d(i,j) > vmax) THEN
       vmax   = d(i,j)
       trinel = 1
      END IF
     END IF
    END DO
    END DO

  END SUBROUTINE max_cart_over_triangle_int
  SUBROUTINE min_cart_over_triangle_dp(  p1, p2, p3, d, x, y, nx, ny, vmin, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear    
    
    IMPLICIT NONE

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
     IF (is_in_triangle(p1,p2,p3,[x(i), y(j)])) THEN
      IF (d(i,j) < vmin) THEN
       vmin   = d(i,j)
       trinel = 1
      END IF
     END IF
    END DO
    END DO

  END SUBROUTINE min_cart_over_triangle_dp
  SUBROUTINE min_cart_over_triangle_int( p1, p2, p3, d, x, y, nx, ny, vmin, trinel)
    ! Find the lowest value v that a data field d on
    ! cartesian grid x,y assumes within the triangle [p1,p2,p3]
    ! If trinel==0, there is no cartesian gridpoint inside the triangle,
    ! probably best to use Cart_Bilinear    
    
    IMPLICIT NONE

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
     IF (is_in_triangle(p1,p2,p3,[x(i), y(j)])) THEN
      IF (d(i,j) < vmin) THEN
       vmin   = d(i,j)
       trinel = 1
      END IF
     END IF
    END DO
    END DO

  END SUBROUTINE min_cart_over_triangle_int
  SUBROUTINE cart_bilinear_dp(  d, x, y, nx, ny, p, v)
    ! Bicubic interpolation of Cartesian data
    ! Interpolates the data field d of size nx,ny, given on grid x,y,
    ! at the point p, resulting in value v    
    
    IMPLICIT NONE

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

  END SUBROUTINE cart_bilinear_dp
  SUBROUTINE cart_bilinear_int( d, x, y, nx, ny, p, v)
    ! Bicubic interpolation of Cartesian data
    ! Interpolates the data field d of size nx,ny, given on grid x,y,
    ! at the point p, resulting in value v    
    
    IMPLICIT NONE

    INTEGER,  DIMENSION(:,:), INTENT(IN)          :: d   ! Data to interpolate
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

  END SUBROUTINE cart_bilinear_int
  SUBROUTINE mesh_bilinear_dp(  mesh, d, p, ti, v)
    ! Interpolate data field d, given on some mesh, to point p, resulting in v
    ! Guess that p is in triangle ti, as a start for find_containing_triangle_oldmesh    
    
    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(INOUT)       :: mesh
    REAL(dp), DIMENSION(:),   INTENT(IN)          :: d
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: ti
    REAL(dp),                 INTENT(OUT)         :: v

    INTEGER                                       :: via, vib, vic
    REAL(dp), DIMENSION(2)                        :: pa, pb, pc
    REAL(dp)                                      :: Atri_abp, Atri_bcp, Atri_cap, Atri_abc

    ! Find the triangle containing p
    CALL find_containing_triangle( mesh, p, ti)

    ! Trilinearly interpolate d to p
    via = mesh%Tri( ti,1)
    vib = mesh%Tri( ti,2)
    vic = mesh%Tri( ti,3)
    
    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)
    
    CALL find_triangle_area( pa, pb, p, Atri_abp)
    CALL find_triangle_area( pb, pc, p, Atri_bcp)
    CALL find_triangle_area( pc, pa, p, Atri_cap)
    Atri_abc = Atri_abp + Atri_bcp + Atri_cap
    
    v = (d( via) * Atri_bcp + d( vib) * Atri_cap + d( vic) * Atri_abp) / Atri_abc

  END SUBROUTINE mesh_bilinear_dp
  SUBROUTINE mesh_bilinear_int( mesh, d, p, ti, v)
    ! Interpolate data field d, given on some mesh, to point p, resulting in v
    ! Guess that p is in triangle ti, as a start for find_containing_triangle_oldmesh    
    
    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(INOUT)       :: mesh
    INTEGER,  DIMENSION(:),   INTENT(IN)          :: d
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: ti
    REAL(dp),                 INTENT(OUT)         :: v

    INTEGER                                       :: via, vib, vic
    REAL(dp), DIMENSION(2)                        :: pa, pb, pc
    REAL(dp)                                      :: Atri_abp, Atri_bcp, Atri_cap, Atri_abc

    ! Find the triangle containing p
    CALL find_containing_triangle( mesh, p, ti)

    ! Trilinearly interpolate d to p
    via = mesh%Tri( ti,1)
    vib = mesh%Tri( ti,2)
    vic = mesh%Tri( ti,3)
    
    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)
    
    CALL find_triangle_area( pa, pb, p, Atri_abp)
    CALL find_triangle_area( pb, pc, p, Atri_bcp)
    CALL find_triangle_area( pc, pa, p, Atri_cap)
    Atri_abc = Atri_abp + Atri_bcp + Atri_cap
    
    v = (REAL(d( via),dp) * Atri_bcp + REAL(d( vib),dp) * Atri_cap + REAL(d( vic),dp) * Atri_abp) / Atri_abc

  END SUBROUTINE mesh_bilinear_int
  SUBROUTINE new_triangle_contains_old_mask( mesh, pa, pb, pc, mask, vi_closest_to_gc, isso)
    ! Used in mesh updating. Given a (new mesh) triangle [pa,pb,pc], check if that triangle
    ! contains any (old) mesh vertices where the (old mesh) mask has value 1    
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(INOUT)       :: mesh
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: pa, pb, pc
    INTEGER,  DIMENSION(:),   INTENT(IN)          :: mask
    INTEGER,                  INTENT(INOUT)       :: vi_closest_to_gc
    LOGICAL,                  INTENT(OUT)         :: isso
    
    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: gc, p
    INTEGER                                       :: ti_containing_gc, via, vib, vic
    REAL(dp), DIMENSION(2)                        :: pia, pib, pic
    REAL(dp)                                      :: Atri_abp, Atri_bcp, Atri_cap, Atri_abc, mdpgc
    INTEGER                                       :: nvi, vi, ci, vc
    
    isso = .FALSE.
    
    ! Use a linear search to find the (old) mesh vertex closest to the geometric center of the triangle
    gc = (pa+pb+pc) / 3._dp
    CALL find_containing_vertex( mesh, gc, vi_closest_to_gc)
    
    ! If that vertex doesn't lie inside the triangle, then probably none of them do - use trilinear interpolation instead.
    p = mesh%V( vi_closest_to_gc,:)
    IF (.NOT. is_in_triangle( pa, pb, pc, p)) THEN
    
      ! Find the (old) mesh triangle containing the (new mesh) triangle's geometric center
      ti_containing_gc = mesh%iTri( vi_closest_to_gc,1)
      CALL find_containing_triangle( mesh, gc, ti_containing_gc)

      ! The three vertices spanning this (old) mesh triangle
      via = mesh%Tri( ti_containing_gc,1)
      vib = mesh%Tri( ti_containing_gc,2)
      vic = mesh%Tri( ti_containing_gc,3)
    
      pia = mesh%V( via,:)
      pib = mesh%V( vib,:)
      pic = mesh%V( vic,:)
      
      CALL find_triangle_area( pia, pib, gc, Atri_abp)
      CALL find_triangle_area( pib, pic, gc, Atri_bcp)
      CALL find_triangle_area( pic, pia, gc, Atri_cap)
      Atri_abc = Atri_abp + Atri_bcp + Atri_cap
    
      mdpgc = (REAL(mask( via),dp) * Atri_bcp + REAL(mask( vib),dp) * Atri_cap + REAL(mask( vic),dp) * Atri_abp) / Atri_abc
      
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
        vi = mesh%VStack1( nvi)
        
        DO ci = 1, mesh%nC( vi)
        
          vc = mesh%C( vi,ci)
          p  = mesh%V( vc,:)
          
          IF (mesh%VMap( vc)==1) CYCLE ! This neighbour has already been checked.
          
          ! Mark this neighbour as checked
          mesh%VMap( vc) = 1
          
          ! If it lies inside the triangle, add it to the new stack.
          IF (is_in_triangle( pa, pb, pc, p)) THEN
          
            ! If it has what we're looking for, stop the search.
            IF (mask( vc) == 1) THEN
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
    
  END SUBROUTINE new_triangle_contains_old_mask
  
! == Rotate a vector field [u,v] to local [p,o] components on the staggered c (edge) grid
  SUBROUTINE rotate_xy_to_po_stag( mesh, u_c, v_c, p_c, o_c)
    ! Rotate a vector field [u,v] to local [p,o] components on the staggered c (edge) grid
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:    ), INTENT(IN)        :: u_c, v_c
    REAL(dp), DIMENSION(:    ), INTENT(OUT)       :: p_c, o_c
    
    ! Local variables:
    INTEGER                                       :: ci, vi, vj
    REAL(dp)                                      :: Dx, Dy, D
    
    DO ci = mesh%ci1, mesh%ci2
    
      vi = mesh%Aci( ci,1)
      vj = mesh%Aci( ci,2)
      
      Dx = mesh%V( vj,1) - mesh%V( vi,1)
      Dy = mesh%V( vj,2) - mesh%V( vi,2)
      D  = SQRT(Dx**2 + Dy**2)
      
      p_c( ci) = u_c( ci) * Dx/D + v_c( ci) * Dy/D
      o_c( ci) = v_c( ci) * Dx/D - u_c( ci) * Dy/D
      
    END DO ! DO ci = mesh%ci1, mesh%ci2
    CALL sync
    
  END SUBROUTINE rotate_xy_to_po_stag
  
! == Diagnostic tools: write a (small) mesh to the screen, and check if mesh data is self-consistent
  SUBROUTINE write_mesh_to_screen( mesh)    
    
    IMPLICIT NONE

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
  
  END SUBROUTINE write_mesh_to_screen
  SUBROUTINE write_mesh_to_text_file( mesh, filename)    
    
    IMPLICIT NONE

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
    WRITE(UNIT = 1337, FMT = '(A,I2)')    '  nC_mem = ', mesh%nC_mem
    WRITE(UNIT = 1337, FMT = '(A,I6)')    '  nV      = ', mesh%nV
    WRITE(UNIT = 1337, FMT = '(A,I6)')    '  nTri    = ', mesh%nTri
    WRITE(UNIT = 1337, FMT = '(A)')       ''
    WRITE(UNIT = 1337, FMT = '(A,I6,A,I3,A)')       'Vertex data: ', mesh%nV, ' rows, ', 2 + 1 + mesh%nC_mem + 1 + mesh%nC_mem + 1, ' columns'
    WRITE(UNIT = 1337, FMT = '(A)')       ''
    
    ! Vertex data
    WRITE(UNIT = 1337, FMT = '(A)')       'V  nC  C  niTri  iTri  edge_index'    
    DO vi = 1, mesh%nV
      WRITE(UNIT = 1337, FMT = '(2F24.14,I3)', ADVANCE = 'NO') mesh%V(vi,1), mesh%V(vi,2), mesh%nC(vi)
      DO ci = 1, mesh%nC_mem
        WRITE(UNIT = 1337, FMT = '(I6)', ADVANCE = 'NO') mesh%C(vi,ci)
      END DO
      WRITE(UNIT = 1337, FMT = '(I3)', ADVANCE = 'NO') mesh%niTri(vi)
      DO ci = 1, mesh%nC_mem
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
    
  END SUBROUTINE write_mesh_to_text_file
  SUBROUTINE check_mesh( mesh)
    ! Check if the mesh data is self-consistent    
    
    IMPLICIT NONE

    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER                                       :: vi, ci, vc, ci2, vc2, iti, iti2, ti, n, v1, v2, v3, ti2, n2
    LOGICAL                                       :: FoundIt
    
   ! IF (.NOT. par%master) RETURN
        
    ! == V
    ! =============================================================
    DO vi = 1, mesh%nV
      IF (mesh%V(vi,1) < mesh%xmin - mesh%tol_dist .OR. mesh%V(vi,1) > mesh%xmax + mesh%tol_dist .OR. &
          mesh%V(vi,2) < mesh%ymin - mesh%tol_dist .OR. mesh%V(vi,2) > mesh%ymax + mesh%tol_dist) THEN
        WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' outside mesh domain! (x = [', &
          mesh%xmin, ', ', mesh%V(vi,1), ',', mesh%xmax, '], y = [', mesh%ymin, ', ', mesh%V(vi,2), ',', mesh%ymax, ']'
      END IF
    END DO
    
    ! == nC
    ! =============================================================
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC(vi)
        IF (mesh%C(vi,ci) == 0) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has fewer connections than nC says!'
      END DO
      DO ci = mesh%nC(vi)+1, mesh%nC_mem
        IF (mesh%C(vi,ci) > 0)  WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has more connections than nC says!'
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
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is connected to ', vc, ', but not the other way round!'
      END DO
    END DO
    
    ! == niTri
    ! =============================================================
    DO vi = 1, mesh%nV
      DO iti = 1, mesh%niTri(vi)
        IF (mesh%iTri(vi,iti) == 0) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has fewer iTriangles than niTri says!'
      END DO
      DO iti = mesh%niTri(vi)+1, mesh%nC_mem
        IF (mesh%iTri(vi,iti) > 0)  WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has more iTriangles than nC says!'
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
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' lists triangle ', ti, ' in iTri, but that triangle ', ti, ' doesnt contain vertex ', vi, '!'
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
          IF (.NOT. (n==2)) WRITE(0,*) ' check_mesh - ERROR: non-edge vertices ', vi, ' and ', vc, ' share ', n, ' triangles'        
        END DO
        
      ELSE ! IF (mesh%edge_index(vi) == 0) THEN
      
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
            IF (.NOT. (n==2)) WRITE(0,*) ' check_mesh - ERROR: edge vertex ', vi, ' and non-edge vertex ', vc, ' share ', n, ' triangles'
          
          ELSE ! IF (mesh%edge_index(vc)==0) THEN
                
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
            IF (.NOT. is_boundary_segment( mesh, vi, vc)) CYCLE
            IF (.NOT. (n==1)) WRITE(0,*) ' check_mesh - ERROR: edge vertices ', vi, ' and ', vc, ' share ', n, ' triangles'
            
          END IF
        END DO
        
      END IF
    END DO
    
    ! == edge_index
    ! =============================================================
    DO vi = 1, mesh%nV
    
      IF (mesh%edge_index(vi) == 0) THEN
      
        IF (mesh%V(vi,1) <= mesh%xmin .OR. mesh%V(vi,1) >= mesh%xmax .OR. mesh%V(vi,2) <= mesh%ymin .OR. mesh%V(vi,2) >= mesh%ymax) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 0 but lies on or beyond the mesh domain boundary!'
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
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 0, but its first and last neighbours are not connected!'
        
      ELSEIF (mesh%edge_index(vi) == 1) THEN
      
        IF (ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 1 but does not lie on the N boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==8 .OR. mesh%edge_index(vc)==1)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 1 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==1 .OR. mesh%edge_index(vc)==2)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 1 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 1 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 1 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 2) THEN
      
        IF (.NOT. vi==3) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as NE corner!'
      
        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 2 but does not lie on the NE corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==8 .OR. mesh%edge_index(vc)==1)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 2 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==3 .OR. mesh%edge_index(vc)==4)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 2 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 2 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 2 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 3) THEN
      
        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 3 but does not lie on the E boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==2 .OR. mesh%edge_index(vc)==3)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 3 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==3 .OR. mesh%edge_index(vc)==4)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 3 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==1 .OR. mesh%Tri_edge_index(ti)==2 .OR. mesh%Tri_edge_index(ti)==3)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 3 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 3 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 4) THEN
      
        IF (.NOT. vi==2) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as SE corner!'
      
        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 4 but does not lie on the SE corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==2 .OR. mesh%edge_index(vc)==3)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 4 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==5 .OR. mesh%edge_index(vc)==6)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 4 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 4 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 4 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 5) THEN
      
        IF (ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 5 but does not lie on the S boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==4 .OR. mesh%edge_index(vc)==5)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 5 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==5 .OR. mesh%edge_index(vc)==6)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 5 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==3 .OR. mesh%Tri_edge_index(ti)==4 .OR. mesh%Tri_edge_index(ti)==5)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 5 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 5 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 6) THEN
      
        IF (.NOT. vi==1) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as SW corner!'
      
        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 6 but does not lie on the SW corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==4 .OR. mesh%edge_index(vc)==5)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 6 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==7 .OR. mesh%edge_index(vc)==8)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 6 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 6 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 6 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 7) THEN
      
        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 7 but does not lie on the W boundary!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==6 .OR. mesh%edge_index(vc)==7)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 7 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==7 .OR. mesh%edge_index(vc)==8)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 7 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==5 .OR. mesh%Tri_edge_index(ti)==6 .OR. mesh%Tri_edge_index(ti)==7)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 7 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 7 but its last iTri doesnt have a matching Tri_edge_index!'
        END IF
        
      ELSEIF (mesh%edge_index(vi) == 8) THEN
      
        IF (.NOT. vi==4) WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' is listed as NW corner!'
      
        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 8 but does not lie on the NW corner!'
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%edge_index(vc)==6 .OR. mesh%edge_index(vc)==7)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 8 but its first connection doesnt have a matching edge_index!'
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%edge_index(vc)==1 .OR. mesh%edge_index(vc)==2)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 8 but its last connection doesnt have a matching edge_index!'
        END IF
        ti = mesh%iTri(vi,1)
        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 8 but its first iTri doesnt have a matching Tri_edge_index!'
        END IF
        ti = mesh%iTri(vi,mesh%niTri(vi))
        IF (.NOT. (mesh%Tri_edge_index(ti)==7 .OR. mesh%Tri_edge_index(ti)==8 .OR. mesh%Tri_edge_index(ti)==1)) THEN
          WRITE(0,*) ' check_mesh - ERROR: vertex ', vi, ' has edge_index 8 but its last iTri doesnt have a matching Tri_edge_index!'
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
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: triangle ', ti, ' contains vertex ', vi, ', but that vertex doesnt list ti as an iTri!'
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
      IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: triangle ', ti, ' contains unconnected vertices ', v1, ' and ', v2, '!'
      
      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v1)
        vc = mesh%C(v1,ci)
        IF (vc==v3) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: triangle ', ti, ' contains unconnected vertices ', v1, ' and ', v3, '!'
      
      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v2)
        vc = mesh%C(v2,ci)
        IF (vc==v3) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: triangle ', ti, ' contains unconnected vertices ', v2, ' and ', v3, '!'
    END DO
    
    ! == TriC
    ! =============================================================
    
    DO ti = 1, mesh%nTri
      DO n = 1, 3
        ti2 = mesh%TriC(ti,n)
        IF (ti2 == 0) THEN
          IF (mesh%Tri_edge_index(ti) == 0) WRITE(0,*) ' check_mesh - ERROR: non-edge triangle ', ti, ' misses a neighbour!'
          CYCLE
        END IF
        FoundIt = .FALSE.
        DO n2 = 1, 3
          IF (mesh%TriC(ti2,n2) == ti) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) WRITE(0,*) ' check_mesh - ERROR: triangle ', ti, ' is connected to ', ti2, ', but not the other way round!'
      END DO
    END DO
        
    
  END SUBROUTINE check_mesh

END MODULE mesh_help_functions_module
