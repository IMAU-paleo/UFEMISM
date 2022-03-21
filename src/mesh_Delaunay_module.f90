MODULE mesh_Delaunay_module

  ! Routines used for updating a Delaunay triangulation by splitting a triangle, a line or a segment.

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
  USE data_types_module,               ONLY: type_mesh
  USE mesh_help_functions_module,      ONLY: is_in_triangle, is_boundary_segment, find_circumcenter, find_containing_triangle, &
                                             line_from_points, line_line_intersection, perpendicular_bisector_from_line, &
                                             encroaches_upon, update_triangle_circumcenter

  IMPLICIT NONE

  CONTAINS
  
  SUBROUTINE split_triangle( mesh, ti, p_new)
     ! Split mesh triangle ti into three new ones, adding a vertex at p_new    
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,                             INTENT(IN)    :: ti
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p_new

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'split_triangle'
    INTEGER                                            :: v1, v2, t1, n, nc, nnext, t_old
    INTEGER                                            :: vi, vj, p1, p2, p3, tn1, tn2, tn3, t_new1, t_new2, t_new3
    REAL(dp)                                           :: la, lb, lc, le, lf, lg, la2, lb2, lc2, le2, lf2, lg2
    REAL(dp), DIMENSION(2)                             :: p, q, cc2
    INTEGER                                            :: nf, e
    INTEGER,  DIMENSION(5)                             :: edgevals, edgevals_left, edgevals_right
    LOGICAL                                            :: isencroached
    INTEGER                                            :: va, vb
    LOGICAL                                            :: did_flip
    
    ! Add routine to path
    CALL init_routine( routine_name)
     
    ! Check fi this triangle's circumcenter encroaches upon a boundary segment.
    ! If so, split that segment.
    CALL encroaches_upon( mesh, p_new(1), p_new(2), isencroached, va, vb)
    IF (isencroached) THEN
      p = (mesh%V(va,:) + mesh%V(vb,:)) / 2._dp
      CALL split_segment( mesh, va, vb, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == If the circumcenter lies outside of the grid, split a segment instead.
    IF (p_new(1)<mesh%xmin) THEN
      IF (p_new(2)<mesh%ymin .OR. p_new(2)>mesh%ymax) THEN
        CALL crash('circumcenter lies way outside of grid!')
      END IF
      ! Find the two vertices of the segment.
      v1 = 0
      v2 = 0
      DO t1 = 1, mesh%nTri
        IF (.NOT. (mesh%Tri_edge_index(t1)==7 .OR. mesh%Tri_edge_index(t1)==0)) CYCLE
        DO n = 1, 3
          nnext = n+1
          IF (nnext==4) nnext=1
          IF (abs(mesh%V(mesh%Tri(t1,n    ),1) - mesh%xmin)<0.0001_dp .AND. mesh%V(mesh%Tri(t1,n),2    )>p_new(2) .AND. &
              abs(mesh%V(mesh%Tri(t1,nnext),1) - mesh%xmin)<0.0001_dp .AND. mesh%V(mesh%Tri(t1,nnext),2)<p_new(2)) THEN
            v1 = mesh%Tri(t1,n)
            v2 = mesh%Tri(t1,nnext)
          END IF
        END DO
      END DO
      IF (v1==0) THEN
        CALL crash('circumcenter lies outside of grid; couldnt find segment to split!')
      END IF
      p = (mesh%V(v1,:) + mesh%V(v2,:)) / 2._dp
      CALL split_segment( mesh, v1, v2, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    IF (p_new(1)>mesh%xmax) THEN
      IF (p_new(2)<mesh%ymin .OR. p_new(2)>mesh%ymax) THEN
        CALL crash('circumcenter lies way outside of grid!')
      END IF
      ! Find the two vertices of the segment.
      v1 = 0;
      v2 = 0;
      DO t1 = 1, mesh%nTri
        IF (.NOT. (mesh%Tri_edge_index(t1)==3 .OR. mesh%Tri_edge_index(t1)==0)) CYCLE
        DO n = 1, 3
          nnext = n+1
          IF (nnext==4) nnext=1
          IF (abs(mesh%V(mesh%Tri(t1,n    ),1) - mesh%xmax)<0.0001_dp .AND. mesh%V(mesh%Tri(t1,n    ),2)<p_new(2) .AND. &
              abs(mesh%V(mesh%Tri(t1,nnext),1) - mesh%xmax)<0.0001_dp .AND. mesh%V(mesh%Tri(t1,nnext),2)>p_new(2)) THEN
            v1 = mesh%Tri(t1,n)
            v2 = mesh%Tri(t1,nnext)
          END IF
        END DO
      END DO
      IF (v1==0) THEN
        CALL crash('circumcenter lies outside of grid; couldnt find segment to split!')
      END IF
      p = (mesh%V(v1,:) + mesh%V(v2,:)) / 2._dp
      CALL split_segment( mesh, v1, v2, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    IF (p_new(2)<mesh%ymin) THEN
      IF (p_new(1)<mesh%xmin .OR. p_new(1)>mesh%xmax) THEN
        CALL crash('circumcenter lies way outside of grid!')
      END IF
      ! Find the two vertices of the segment.
      v1 = 0
      v2 = 0
      DO t1 = 1, mesh%nTri
        IF (.NOT. (mesh%Tri_edge_index(t1)==5 .OR. mesh%Tri_edge_index(t1)==0)) CYCLE
        DO n = 1, 3
          nnext = n+1
          IF (nnext==4) nnext=1
          IF (abs(mesh%V(mesh%Tri(t1,n    ),2) - mesh%ymin)<0.0001_dp .AND. mesh%V(mesh%Tri(t1,n),1    )<p_new(1) .AND. &
              abs(mesh%V(mesh%Tri(t1,nnext),2) - mesh%ymin)<0.0001_dp .AND. mesh%V(mesh%Tri(t1,nnext),1)>p_new(1)) THEN
            v1 = mesh%Tri(t1,n)
            v2 = mesh%Tri(t1,nnext)
          END IF
        END DO
      END DO
      IF (v1==0) THEN
        CALL crash('circumcenter lies outside of grid; couldnt find segment to split!')
      END IF
      p = (mesh%V(v1,:) + mesh%V(v2,:)) / 2._dp
      CALL split_segment( mesh, v1, v2, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    IF (p_new(2)>mesh%ymax) THEN
      IF (p_new(1)<mesh%xmin .OR. p_new(1)>mesh%xmax) THEN
        CALL crash('circumcenter lies way outside of grid!')
      END IF
      ! Find the two vertices of the segment.
      v1 = 0
      v2 = 0
      DO t1 = 1, mesh%nTri
        IF (.NOT. (mesh%Tri_edge_index(t1)==1 .OR. mesh%Tri_edge_index(t1)==0)) CYCLE
        DO n = 1, 3
          nnext = n+1
          IF (nnext==4) nnext=1
          IF (abs(mesh%V(mesh%Tri(t1,n    ),2) - mesh%ymax)<0.0001_dp .AND. mesh%V(mesh%Tri(t1,n),1    )>p_new(1) .AND. &
              abs(mesh%V(mesh%Tri(t1,nnext),2) - mesh%ymax)<0.0001_dp .AND. mesh%V(mesh%Tri(t1,nnext),1)<p_new(1)) THEN
            v1 = mesh%Tri(t1,n)
            v2 = mesh%Tri(t1,nnext)
          END IF
        END DO
      END DO
      IF (v1==0) THEN
        CALL crash('circumcenter lies outside of grid; couldnt find segment to split!')
      END IF
      p = (mesh%V(v1,:) + mesh%V(v2,:)) / 2._dp
      CALL split_segment( mesh, v1, v2, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == Find the triangle containing this circumcenter
    ! Use the triangle whose circumcenter it is as a first guess, search
    ! outward from there using a flood-fill algorithm.
    t_old = ti
    CALL find_containing_triangle( mesh, p_new, t_old)

    ! == If the new vertex is (almost) colinear with two vertices of
    ! the containing triangle, split the line between them instead.
    DO n = 1, 3
      nnext = n+1
      IF (nnext==4) nnext = 1
      vi = mesh%Tri(t_old,n)
      vj = mesh%Tri(t_old,nnext)

      p = mesh%V(vi,:)
      q = mesh%V(vj,:)
      CALL line_from_points(p, p_new, la, lb, lc)
      CALL line_from_points(p_new, q, le, lf, lg)

      CALL perpendicular_bisector_from_line(p,p_new,la,lb,la2,lb2,lc2)
      CALL perpendicular_bisector_from_line(p_new,q,le,lf,le2,lf2,lg2);

      CALL line_line_intersection(la2,lb2,lc2,le2,lf2,lg2,cc2)

      IF (abs(cc2(1))>(mesh%xmax-mesh%xmin)*1000._dp .OR. abs(cc2(2))>(mesh%ymax-mesh%ymin)*1000._dp) THEN
        CALL split_line( mesh, vi, vj, p_new)
        CALL finalise_routine( routine_name)
        RETURN
      END IF
    END DO
    
    !WRITE(0,*) 'Split triangle ', t_old

    ! == Find that triangle's vertices and neighbours
    p1  = mesh%Tri(t_old,1)
    p2  = mesh%Tri(t_old,2)
    p3  = mesh%Tri(t_old,3)
    tn1 = mesh%TriC(t_old,1)
    tn2 = mesh%TriC(t_old,2)
    tn3 = mesh%TriC(t_old,3)

    ! == Add the new vertex to the mesh
    mesh%nV = mesh%nV+1
    mesh%V(mesh%nV,:) = p_new
    mesh%edge_index(mesh%nV) = 0

    ! == Replace ti_old by three new triangles
    t_new1 = t_old
    t_new2 = mesh%nTri+1
    t_new3 = mesh%nTri+2
    mesh%Tri(t_new1,:) = [p1, p2, mesh%nV]
    mesh%Tri(t_new2,:) = [p2, p3, mesh%nV]
    mesh%Tri(t_new3,:) = [p3, p1, mesh%nV]
    mesh%nTri = mesh%nTri+2

    ! Add these to the mesh%RefStack (if they're not in there already)
    IF (mesh%RefMap(t_new1)==0) THEN
      mesh%RefMap(t_new1) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = t_new1
    END IF
    IF (mesh%RefMap(t_new2)==0) THEN
      mesh%RefMap(t_new2) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = t_new2
    END IF
    IF (mesh%RefMap(t_new3)==0) THEN
      mesh%RefMap(t_new3) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = t_new3
    END IF

    ! == Update vertex connectivity matrix
    ! p1
    DO n = 1, mesh%nC(p1)
      IF (mesh%C(p1,n)==p2) THEN
        mesh%C(p1,:) = [mesh%C(p1,1:n), mesh%nV, mesh%C(p1,n+1:mesh%nC_mem-1)]
        mesh%nC(p1) = mesh%nC(p1)+1
        EXIT
      END IF
    END DO
    ! p2
    DO n = 1, mesh%nC(p2)
      IF (mesh%C(p2,n)==p3) THEN
        mesh%C(p2,:) = [mesh%C(p2,1:n), mesh%nV, mesh%C(p2,n+1:mesh%nC_mem-1)]
        mesh%nC(p2) = mesh%nC(p2)+1
        EXIT
      END IF
    END DO
    ! p3
    DO n = 1, mesh%nC(p3)
      IF (mesh%C(p3,n)==p1) THEN
        mesh%C(p3,:) = [mesh%C(p3,1:n), mesh%nV, mesh%C(p3,n+1:mesh%nC_mem-1)]
        mesh%nC(p3) = mesh%nC(p3)+1
        EXIT
      END IF
    END DO
    ! new vertex
    mesh%C(mesh%nV,1:3) = [p1, p2, p3]
    mesh%nC(mesh%nV)    = 3

    ! == Update triangle connectivity matrix
    ! (existing) neighbours
    DO n = 1, 3
      IF (tn1>0) THEN
        IF (mesh%TriC(tn1,n)==t_old) mesh%TriC(tn1,n) = t_new2
      END IF
      IF (tn2>0) THEN
        IF (mesh%TriC(tn2,n)==t_old) mesh%TriC(tn2,n) = t_new3
      END IF
      IF (tn3>0) THEN
        IF (mesh%TriC(tn3,n)==t_old) mesh%TriC(tn3,n) = t_new1
      END IF
    END DO
    ! new triangles
    mesh%TriC(t_new1,:) = [t_new2, t_new3, tn3]
    mesh%TriC(t_new2,:) = [t_new3, t_new1, tn1]
    mesh%TriC(t_new3,:) = [t_new1, t_new2, tn2]

    ! == Update triangle edge indices
    edgevals       = [0, 1, 3, 5, 7]
    edgevals_left  = [0, 8, 2, 4, 6]
    edgevals_right = [0, 2, 4, 6, 8]
    DO e = 1, 5
      nc = 0;
      DO n = 1, 3
        IF (mesh%edge_index(mesh%Tri(t_new1,n)) == edgevals(e) .OR. &
            mesh%edge_index(mesh%Tri(t_new1,n)) == edgevals_left(e) .OR. &
            mesh%edge_index(mesh%Tri(t_new1,n)) == edgevals_right(e)) THEN
         nc = nc+1
        END IF
      END DO
      IF (nc==2) mesh%Tri_edge_index(t_new1) = edgevals(e)

      nc = 0;
      DO n = 1, 3
        IF (mesh%edge_index(mesh%Tri(t_new2,n)) == edgevals(e) .OR. &
            mesh%edge_index(mesh%Tri(t_new2,n)) == edgevals_left(e) .OR. &
            mesh%edge_index(mesh%Tri(t_new2,n)) == edgevals_right(e)) THEN
         nc = nc+1
        END IF
      END DO
      IF (nc==2) mesh%Tri_edge_index(t_new2) = edgevals(e)

      nc = 0;
      DO n = 1, 3
        IF (mesh%edge_index(mesh%Tri(t_new3,n)) == edgevals(e) .OR. &
            mesh%edge_index(mesh%Tri(t_new3,n)) == edgevals_left(e) .OR. &
            mesh%edge_index(mesh%Tri(t_new3,n)) == edgevals_right(e)) THEN
         nc = nc+1
        END IF
      END DO
      IF (nc==2) mesh%Tri_edge_index(t_new3) = edgevals(e)
    END DO

    ! == Update inverse triangle matrix
    ! p1
    DO n = 1, mesh%niTri(p1)
      IF (mesh%iTri(p1,n)==t_old) THEN
        mesh%iTri(p1,:) = [mesh%iTri(p1,1:n-1), t_new1, t_new3, mesh%iTri(p1,n+1:mesh%nC_mem-1)]
        mesh%niTri(p1) = mesh%niTri(p1)+1
        EXIT
      END IF
    END DO
    ! p2
    DO n = 1, mesh%niTri(p2)
      IF (mesh%iTri(p2,n)==t_old) THEN
        mesh%iTri(p2,:) = [mesh%iTri(p2,1:n-1), t_new2, t_new1, mesh%iTri(p2,n+1:mesh%nC_mem-1)]
        mesh%niTri(p2) = mesh%niTri(p2)+1
        EXIT
      END IF
    END DO
    ! p3
    DO n = 1, mesh%niTri(p3)
      IF (mesh%iTri(p3,n)==t_old) THEN
        mesh%iTri(p3,:) = [mesh%iTri(p3,1:n-1), t_new3, t_new2, mesh%iTri(p3,n+1:mesh%nC_mem-1)]
        mesh%niTri(p3) = mesh%niTri(p3)+1
        EXIT
      END IF
    END DO
    ! new vertex
    mesh%iTri(mesh%nV,1:3) = [t_new1, t_new2, t_new3]
    mesh%niTri(mesh%nV) = 3

    ! == Update triangle circumcenters
    CALL update_triangle_circumcenter( mesh, t_new1)
    CALL update_triangle_circumcenter( mesh, t_new2)
    CALL update_triangle_circumcenter( mesh, t_new3)

    ! == Propagate flip operations outward
    ! Start with the newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    mesh%Triflip = 0
    nf           = 0

    IF (tn1>0) THEN
      nf=nf+1
      mesh%Triflip(nf,:) = [t_new2, tn1]
    END IF
    IF (tn2>0) THEN
      nf=nf+1
      mesh%Triflip(nf,:) = [t_new3, tn2]
    END IF
    IF (tn3>0) THEN
      nf=nf+1
      mesh%Triflip(nf,:) = [t_new1, tn3]
    END IF

    DO WHILE (nf>0)
      CALL flip_triangle_pairs( mesh, nf, did_flip)
    END DO
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE split_triangle
  SUBROUTINE split_line(     mesh, v1a, v2a, p_new)
     ! Split the line between vertices v1a and v2a at point p_new    
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,                             INTENT(IN)    :: v1a, v2a
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p_new

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'split_line'
    INTEGER                                            :: v1, v2, ti, t1, t2, n, nc, nnext, vo1, vo2
    INTEGER                                            :: t1new1, t1new2, t2new1, t2new2, t1nv1, t1nv2, t2nv1, t2nv2
    LOGICAL                                            :: AreConnected, SwitchThem
    REAL(dp), DIMENSION(2)                             :: p
    INTEGER                                            :: nf, e
    INTEGER,  DIMENSION(5)                             :: edgevals, edgevals_left, edgevals_right
    LOGICAL                                            :: did_flip
    
    ! Add routine to path
    CALL init_routine( routine_name)

    v1 = v1a
    v2 = v2a

    ! Check if they are even connected.
    AreConnected = .FALSE.
    DO n = 1, mesh%nC(v1)
      IF (mesh%C(v1,n)==v2) AreConnected = .TRUE.
    END DO
    IF (.NOT. AreConnected) THEN
      CALL crash('trying to split a non-existing line!')
    END IF

    ! If both of those are boundary vertices, the line is a Segment - refer
    ! to that function instead.
    IF (is_boundary_segment(mesh,v1,v2)) THEN
      p = (mesh%V(v1,:) + mesh%V(v2,:)) / 2._dp
      CALL split_segment( mesh, v1, v2, p)
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Find the triangles t1 and t2 that contain v1 and v2
    t1 = 0
    t2 = 0
    DO ti = 1, mesh%nTri
      nc = 0
      DO n = 1, 3
        IF (mesh%Tri(ti,n)==v1 .OR. mesh%Tri(ti,n)==v2) nc = nc+1
      END DO
      IF (nc==2) THEN
        IF (t1==0) THEN
         t1=ti
        ELSE
         t2=ti
        END IF
      END IF
    END DO

    IF (t1==0 .OR. t2==0) THEN
      CALL crash('couldnt find two triangles containing both vertices!')
    END IF

    ! Order v1 and v2 anticlockwise in triangle t1
    SwitchThem = .FALSE.
    DO n = 1, 3
      nnext=n+1
      IF (nnext==4) nnext=1
      IF (mesh%Tri(t1,n)==v1) THEN
        IF (mesh%Tri(t1,nnext)/=v2) SwitchThem = .TRUE.
      END IF
    END DO
    IF (SwitchThem) THEN
      v1 = v1+v2
      v2 = v1-v2
      v1 = v1-v2
    END IF
    
    !WRITE(0,*) 'Split line [', v1, ' - ', v2, ']'

    ! == Find the other relevant indices - non-shared vertices and neighbour triangles
    ! t1nv1: triangle next to t1, across from v1
    ! t1nv2: triangle next to t1, across from v2
    ! t2nv1: triangle next to t2, across from v1
    ! t2nv2: triangle next to t2, across from v2
    t1nv1 = 0
    t1nv2 = 0
    t2nv1 = 0
    t2nv2 = 0
    vo1   = 0
    vo2   = 0
    DO n = 1, 3
      IF (mesh%Tri(t1,n)==v1) THEN
        t1nv1 = mesh%TriC(t1,n)
      ELSEIF (mesh%Tri(t1,n)==v2) THEN
        t1nv2 = mesh%TriC(t1,n)
      ELSE
        vo1 = mesh%Tri(t1,n);
      END IF
      IF (mesh%Tri(t2,n)==v1) THEN
        t2nv1 = mesh%TriC(t2,n)
      ELSEIF (mesh%Tri(t2,n)==v2) THEN
        t2nv2 = mesh%TriC(t2,n)
      ELSE
        vo2 = mesh%Tri(t2,n);
      END IF
    END DO

    ! == Add new vertex to mesh
    mesh%nV = mesh%nV+1
    mesh%V(mesh%nV,:) = p_new
    mesh%edge_index(mesh%nV) = 0

    ! == Create four new triangles
    t1new1 = t1
    t1new2 = mesh%nTri+1
    t2new1 = t2
    t2new2 = mesh%nTri+2
    mesh%Tri(t1new1,:) = [vo1, mesh%nV, v2]
    mesh%Tri(t1new2,:) = [vo1, v1, mesh%nV]
    mesh%Tri(t2new1,:) = [vo2, mesh%nV, v1]
    mesh%Tri(t2new2,:) = [vo2, v2, mesh%nV]
    mesh%nTri = mesh%nTri+2

    ! Add these to the mesh%RefStack (if they're not in there already)
    IF (mesh%RefMap(t1new1)==0) THEN
      mesh%RefMap(t1new1) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = t1new1
    END IF
    IF (mesh%RefMap(t1new2)==0) THEN
      mesh%RefMap(t1new2) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = t1new2
    END IF
    IF (mesh%RefMap(t2new1)==0) THEN
      mesh%RefMap(t2new1) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = t2new1
    END IF
    IF (mesh%RefMap(t2new2)==0) THEN
      mesh%RefMap(t2new2) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = t2new2
    END IF

    ! == Find circumcenters
    CALL update_triangle_circumcenter( mesh, t1new1)
    CALL update_triangle_circumcenter( mesh, t1new2)
    CALL update_triangle_circumcenter( mesh, t2new1)
    CALL update_triangle_circumcenter( mesh, t2new2)

    ! == Update inverse triangle matrix
    ! vo1
    DO n = 1, mesh%niTri(vo1)
      IF (mesh%iTri(vo1,n)==t1) THEN
        mesh%iTri(vo1,:) = [mesh%iTri(vo1,1:n-1), t1new2, t1new1, mesh%iTri(vo1,n+1:mesh%nC_mem-1)]
        mesh%niTri(vo1) = mesh%niTri(vo1)+1
        EXIT
      END IF
    END DO
    ! v1
    DO n = 1, mesh%niTri(v1)
      IF (mesh%iTri(v1,n)==t1) mesh%iTri(v1,n) = t1new2
      IF (mesh%iTri(v1,n)==t2) mesh%iTri(v1,n) = t2new1
    END DO
    ! vo2
    DO n = 1, mesh%niTri(vo2)
      IF (mesh%iTri(vo2,n)==t2) THEN
        mesh%iTri(vo2,:) = [mesh%iTri(vo2,1:n-1), t2new2, t2new1, mesh%iTri(vo2,n+1:mesh%nC_mem-1)]
        mesh%niTri(vo2) = mesh%niTri(vo2)+1
        EXIT
      END IF
    END DO
    ! v2
    DO n = 1, mesh%niTri(v2)
      IF (mesh%iTri(v2,n)==t1) mesh%iTri(v2,n) = t1new1
      IF (mesh%iTri(v2,n)==t2) mesh%iTri(v2,n) = t2new2
    END DO
    ! new vertex
    mesh%iTri(mesh%nV,1:4) = [t1new2, t2new1, t2new2, t1new1]
    mesh%niTri(mesh%nV) = 4

    ! == Update triangle connectivity matrix
    ! t1nv2
    IF (t1nv2>0) THEN
      DO n = 1, 3
        IF (mesh%TriC(t1nv2,n) == t1) mesh%TriC(t1nv2,n) = t1new2
      END DO
    END IF
    ! t1nv1
    IF (t1nv1>0) THEN
      DO n = 1, 3
        IF (mesh%TriC(t1nv1,n) == t1) mesh%TriC(t1nv1,n) = t1new1
      END DO
    END IF
    ! t2nv2
    IF (t2nv2>0) THEN
      DO n = 1, 3
        IF (mesh%TriC(t2nv2,n) == t2) mesh%TriC(t2nv2,n) = t2new1
      END DO
    END IF
    ! t2nv1
    IF (t2nv1>0) THEN
      DO n = 1, 3
        IF (mesh%TriC(t2nv1,n) == t2) mesh%TriC(t2nv1,n) = t2new2
      END DO
    END IF

    ! The four new triangles
    mesh%TriC(t1new1,:) = [t2new2, t1nv1, t1new2]
    mesh%TriC(t1new2,:) = [t2new1, t1new1, t1nv2]
    mesh%TriC(t2new1,:) = [t1new2, t2nv2, t2new2]
    mesh%TriC(t2new2,:) = [t1new1, t2new1, t2nv1]

    ! == Update triangle edge indices
    edgevals       = [0, 1, 3, 5, 7]
    edgevals_left  = [0, 8, 2, 4, 6]
    edgevals_right = [0, 2, 4, 6, 8]
    !for ti = [t1new1, t1new2, t2new1, t2new2]
    DO e = 1, 5
      nc = 0;
      DO n = 1, 3
        IF (mesh%edge_index(mesh%Tri(t1new1,n)) == edgevals(e) .OR. &
            mesh%edge_index(mesh%Tri(t1new1,n)) == edgevals_left(e) .OR. &
            mesh%edge_index(mesh%Tri(t1new1,n)) == edgevals_right(e)) THEN
         nc = nc+1
        END IF
      END DO
      IF (nc==2) mesh%Tri_edge_index(t1new1) = edgevals(e)

      nc = 0;
      DO n = 1, 3
        IF (mesh%edge_index(mesh%Tri(t1new2,n)) == edgevals(e) .OR. &
            mesh%edge_index(mesh%Tri(t1new2,n)) == edgevals_left(e) .OR. &
            mesh%edge_index(mesh%Tri(t1new2,n)) == edgevals_right(e)) THEN
         nc = nc+1
        END IF
      END DO
      IF (nc==2) mesh%Tri_edge_index(t1new2) = edgevals(e)

      nc = 0;
      DO n = 1, 3
        IF (mesh%edge_index(mesh%Tri(t2new1,n)) == edgevals(e) .OR. &
            mesh%edge_index(mesh%Tri(t2new1,n)) == edgevals_left(e) .OR. &
            mesh%edge_index(mesh%Tri(t2new1,n)) == edgevals_right(e)) THEN
         nc = nc+1
        END IF
      END DO
      IF (nc==2) mesh%Tri_edge_index(t2new1) = edgevals(e)

      nc = 0;
      DO n = 1, 3
        IF (mesh%edge_index(mesh%Tri(t2new2,n)) == edgevals(e) .OR. &
            mesh%edge_index(mesh%Tri(t2new2,n)) == edgevals_left(e) .OR. &
            mesh%edge_index(mesh%Tri(t2new2,n)) == edgevals_right(e)) THEN
         nc = nc+1
        END IF
      END DO
      IF (nc==2) mesh%Tri_edge_index(t2new2) = edgevals(e)
    END DO

    ! == Update vertex connectivity matrix
    ! po1
    DO n = 1, mesh%nC(vo1)
      IF (mesh%C(vo1,n)==v1) THEN
        mesh%C(vo1,:) = [mesh%C(vo1,1:n), mesh%nV, mesh%C(vo1,n+1:mesh%nC_mem-1)]
        mesh%nC(vo1) = mesh%nC(vo1)+1
        EXIT
      END IF
    END DO
    ! v1
    DO n = 1, mesh%nC(v1)
      IF (mesh%C(v1,n)==v2) THEN
        mesh%C(v1,n) = mesh%nV
        EXIT
      END IF
    END DO
    ! vo2
    DO n = 1, mesh%nC(vo2)
      IF (mesh%C(vo2,n)==v2) THEN
        mesh%C(vo2,:) = [mesh%C(vo2,1:n), mesh%nV, mesh%C(vo2,n+1:mesh%nC_mem-1)]
        mesh%nC(vo2) = mesh%nC(vo2)+1
        EXIT
      END IF
    END DO
    ! v2
    DO n = 1, mesh%nC(v2)
      IF (mesh%C(v2,n)==v1) THEN
        mesh%C(v2,n) = mesh%nV
        EXIT
      END IF
    END DO
    ! new vertex
    mesh%C(mesh%nV,1:4) = [v1, vo2, v2, vo1]
    mesh%nC(mesh%nV) = 4

    ! == Propagate flip operations outward
    ! Start with the newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    mesh%Triflip = 0
    nf           = 0

    IF (t1nv1>0) THEN
      nf=nf+1
      mesh%Triflip(nf,:) = [t1nv1, t1new1]
    END IF
    IF (t1nv2>0) THEN
      nf=nf+1
      mesh%Triflip(nf,:) = [t1nv2, t1new2]
    END IF
    IF (t2nv1>0) THEN
      nf=nf+1
      mesh%Triflip(nf,:) = [t2nv1, t2new2]
    END IF
    IF (t2nv2>0) THEN
      nf=nf+1
      mesh%Triflip(nf,:) = [t2nv2, t2new1]
    END IF

    DO WHILE (nf>0)
      CALL flip_triangle_pairs( mesh, nf, did_flip)
    END DO
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE split_line
  SUBROUTINE split_segment(  mesh, v1a, v2a, p_new)
    ! Split an Edge segment in two, adding a vertex halfway.    
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,                             INTENT(IN)    :: v1a, v2a
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p_new

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'split_segment'
    INTEGER                                            :: v1, v2, ti, t1, n, nc, nnext, tnv1, tnv2, vo
    INTEGER                                            :: tnew1, tnew2
    LOGICAL                                            :: SwitchThem
    INTEGER                                            :: nf
    LOGICAL                                            :: did_flip
    
    ! Add routine to path
    CALL init_routine( routine_name)

    v1 = v1a
    v2 = v2a
     
    !IF (par%master) WRITE(0,*) 'Split segment [', v1, ' - ', v2, ']'

    ! Find the triangle t1 that contains v1 and v2
    t1 = 0
    DO ti = 1, mesh%nTri
      nc = 0
      DO n = 1, 3
        IF (mesh%Tri(ti,n)==v1 .OR. mesh%Tri(ti,n)==v2) nc = nc+1
      END DO
      IF (nc==2) t1=ti
    END DO

    IF (t1==0) THEN
      CALL crash('couldnt find triangle containing both vertices!')
    END IF
    
    IF (mesh%edge_index(v1)==0 .OR. mesh%edge_index(v2)==0) THEN
      CALL crash('segment isnt made up of boundary vertices!')
    END IF

    ! Order v1 and v2 anticlockwise in triangle t1
    SwitchThem = .FALSE.
    DO n = 1, 3
      nnext=n+1
      IF (nnext==4) nnext=1
      IF (mesh%Tri(t1,n)==v1) THEN
        IF (mesh%Tri(t1,nnext)/=v2) SwitchThem = .TRUE.
      END IF
    END DO
    IF (SwitchThem) THEN
      v1 = v1+v2
      v2 = v1-v2
      v1 = v1-v2
    END IF

    ! == Find the other relevant indices - non-shared vertex and neighbour triangles
    ! tnv1: triangle next to t1, across from v1
    ! tnv2: triangle next to t1, across from v2
    tnv1 = 0
    tnv2 = 0
    vo   = 0
    DO n = 1, 3
      IF (mesh%Tri(t1,n)==v1) THEN
        tnv1 = mesh%TriC(t1,n)
      ELSEIF (mesh%Tri(t1,n)==v2) THEN
        tnv2 = mesh%TriC(t1,n)
      ELSE
        vo = mesh%Tri(t1,n)
      END IF
    END DO

    ! == Add new vertex to mesh
    mesh%nV = mesh%nV+1
    mesh%V(mesh%nV,:) = p_new

    ! == Determine edge index
    IF (       mesh%edge_index(v1)==1) THEN
      mesh%edge_index(mesh%nV) = 1
    ELSEIF (   mesh%edge_index(v1)==2) THEN
      IF (     mesh%edge_index(v2)==1 .OR. mesh%edge_index(v2)==8) THEN
        mesh%edge_index(mesh%nV) = 1
      ELSEIF ( mesh%edge_index(v2)==3 .OR. mesh%edge_index(v2)==4) THEN
        mesh%edge_index(mesh%nV) = 3
      ELSE
        CALL crash('edge indices of v1 and v2 dont make sense!')
      END IF
    ELSEIF (   mesh%edge_index(v1)==3) THEN
      mesh%edge_index(mesh%nV) = 3
    ELSEIF (   mesh%edge_index(v1)==4) THEN
      IF (     mesh%edge_index(v2)==3 .OR. mesh%edge_index(v2)==2) THEN
        mesh%edge_index(mesh%nV) = 3
      ELSEIF ( mesh%edge_index(v2)==5 .OR. mesh%edge_index(v2)==6) THEN
        mesh%edge_index(mesh%nV) = 5
      ELSE
        CALL crash('edge indices of v1 and v2 dont make sense!')
      END IF
    ELSEIF (   mesh%edge_index(v1)==5) THEN
      mesh%edge_index(mesh%nV) = 5
    ELSEIF (   mesh%edge_index(v1)==6) THEN
      IF (     mesh%edge_index(v2)==5 .OR. mesh%edge_index(v2)==4) THEN
        mesh%edge_index(mesh%nV) = 5
      ELSEIF ( mesh%edge_index(v2)==7 .OR. mesh%edge_index(v2)==8) THEN
        mesh%edge_index(mesh%nV) = 7
      ELSE
        CALL crash('edge indices of v1 and v2 dont make sense!')
      END IF
    ELSEIF (   mesh%edge_index(v1)==7) THEN
      mesh%edge_index(mesh%nV) = 7
    ELSEIF (   mesh%edge_index(v1)==8) THEN
      IF (     mesh%edge_index(v2)==7 .OR. mesh%edge_index(v2)==6) THEN
        mesh%edge_index(mesh%nV) = 7
      ELSEIF ( mesh%edge_index(v2)==1 .OR. mesh%edge_index(v2)==2) THEN
        mesh%edge_index(mesh%nV) = 1
      ELSE
        CALL crash('edge indices of v1 and v2 dont make sense!')
      END IF
    ELSE
      CALL crash('edge indices of v1 and v2 dont make sense!')
    END IF

    ! == Create two new triangles
    tnew1 = t1
    tnew2 = mesh%nTri+1
    mesh%Tri(tnew1,:) = [v1, mesh%nV, vo]
    mesh%Tri(tnew2,:) = [v2, vo, mesh%nV]
    mesh%nTri = mesh%nTri+1

    ! Add these to the mesh%RefStack (if they're not in there already)
    IF (mesh%RefMap(tnew1)==0) THEN
      mesh%RefMap(tnew1) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = tnew1
    END IF
    IF (mesh%RefMap(tnew2)==0) THEN
      mesh%RefMap(tnew2) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = tnew2
    END IF

    ! == Tri_edge_index
    mesh%Tri_edge_index(tnew1) = mesh%Tri_edge_index(t1)
    mesh%Tri_edge_index(tnew2) = mesh%Tri_edge_index(t1)

    ! == Update triangle Circumcenters
    CALL update_triangle_circumcenter( mesh, tnew1)
    CALL update_triangle_circumcenter( mesh, tnew2)

    ! == Update inverse triangle matrix
    ! vo
    DO n = 1, mesh%niTri(vo)
      IF (mesh%iTri(vo,n)==t1) THEN
        mesh%iTri(vo,:) = [mesh%iTri(vo,1:n-1), tnew1, tnew2, mesh%iTri(vo,n+1:mesh%nC_mem-1)];
        mesh%niTri(vo) = mesh%niTri(vo)+1
        EXIT
      END IF
    END DO
    ! v1 - nothing changes here
    ! v2
    DO n = 1, mesh%niTri(v2)
      IF (mesh%iTri(v2,n)==t1) mesh%iTri(v2,n) = tnew2;
    END DO
    ! new vertex
    mesh%iTri(mesh%nV,1:2) = [tnew2, tnew1]
    mesh%niTri(mesh%nV) = 2;

    ! == Update triangle connectivity matrix
    ! tnv1
    IF (tnv1>0) THEN
      DO n = 1, 3
        IF (mesh%TriC(tnv1,n) == t1) mesh%TriC(tnv1,n) = tnew2
      END DO
    END IF
    ! tnv2
    IF (tnv2>0) THEN
      DO n = 1, 3
        IF (mesh%TriC(tnv2,n) == t1) mesh%TriC(tnv2,n) = tnew1
      END DO
    END IF

    ! The two new triangles
    mesh%TriC(tnew1,:) = [tnew2, tnv2, 0]
    mesh%TriC(tnew2,:) = [tnew1, 0, tnv1]

    ! == Update vertex connectivity matrix
    ! vo
    DO n = 1, mesh%nC(vo)
      IF (mesh%C(vo,n)==v1) THEN
        mesh%C(vo,:) = [mesh%C(vo,1:n), mesh%nV, mesh%C(vo,n+1:mesh%nC_mem-1)]
        mesh%nC(vo) = mesh%nC(vo)+1
        EXIT
      END IF
    END DO
    ! v1
    DO n = 1, mesh%nC(v1)
      IF (mesh%C(v1,n)==v2) THEN
        mesh%C(v1,n) = mesh%nV
        EXIT
      END IF
    END DO
    ! v2
    DO n = 1, mesh%nC(v2)
      IF (mesh%C(v2,n)==v1) THEN
        mesh%C(v2,n) = mesh%nV
        EXIT
      END IF
    END DO
    ! new vertex
    mesh%C(mesh%nV,1:3) = [v2, vo, v1]
    mesh%nC(mesh%nV) = 3

    ! == Propagate flip operations outward
    ! Start with the newly created triangles and their neighbours. Any
    ! new possible flip pairs generated by a flip pair are added to the
    ! list by flip_triangle_pairs, making sure the loop runs until all flip
    ! operations have been done.

    mesh%Triflip = 0
    nf           = 0

    IF (tnv2>0) THEN
      nf = nf+1
      mesh%Triflip(nf,:) = [tnew1, tnv2]
    END IF
    IF (tnv1>0) THEN
      nf = nf+1
      mesh%Triflip(nf,:) = [tnew2, tnv1]
    END IF

    DO WHILE (nf>0)
      CALL flip_triangle_pairs( mesh, nf, did_flip)
    END DO
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE split_segment
  SUBROUTINE flip_triangle_pairs( mesh, nf, did_flip)
    ! Flip adjacent triangles, if possible and neccesary. Add new triangle
    ! pairs to the list.    
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,                             INTENT(INOUT) :: nf
    LOGICAL,                             INTENT(OUT)   :: did_flip

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'flip_triangle_pairs'
    INTEGER                                            :: t1, t2, n, vo1, vo2, v1, v2, t1nv1, t1nv2, t2nv1, t2nv2
    INTEGER                                            :: e, contains_edge
    LOGICAL                                            :: n1to2, n2to1, FlipThem
    INTEGER,  DIMENSION(5)                             :: edgevals, edgevals_left, edgevals_right
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    did_flip = .FALSE.
    
    t1 = mesh%Triflip(1,1)
    t2 = mesh%Triflip(1,2)
    
    IF (t1 == 0 .OR. t2 == 0) THEN
      CALL crash('received t=0!')
    END IF
    
    ! == First, check if the two are really adjacent. If not, that's because of an earlier flip operation.
    ! The current one is now redundant, so remove it from the list.
    n1to2 = .FALSE.
    n2to1 = .FALSE.
    DO n = 1, 3
      IF (mesh%TriC(t1,n)==t2) n1to2 = .TRUE.
      IF (mesh%TriC(t2,n)==t1) n2to1 = .TRUE.
    END DO

    IF ((n1to2 .AND. .NOT. n2to1) .OR. (n2to1 .AND. .NOT. n1to2)) THEN
      CALL crash('somethings really wrong with the triangle connectivity matrix!')
    END IF
    IF (.NOT. n1to2 .AND. .NOT. n2to1) THEN
      ! The two triangles are no longer connected; remove them from the list and return.
      mesh%Triflip(1:mesh%nTri-1,:) = mesh%Triflip(2:mesh%nTri,:)
      nf = nf-1
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == Check if a flip is necessary
    ! If not, remove the pair from the flip list.
    CALL need_flipping( mesh, t1, t2, FlipThem, vo1, vo2, v1, v2, t1nv1, t1nv2, t2nv1, t2nv2)

    IF (.NOT. FlipThem) THEN
      ! The two triangles do not need to be flipped; remove them from the list and return.
      mesh%Triflip(1:mesh%nTri-1,:) = mesh%Triflip(2:mesh%nTri,:)
      nf = nf-1
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! == Flip them
    
    did_flip = .TRUE.
    
    ! WRITE(0,'(A,I6,A,I6)') '     Flipping triangles ', t1, ' and ', t2

    ! == Update the triangle matrix
    mesh%Tri(t1,:) = [vo1, v1, vo2]
    mesh%Tri(t2,:) = [vo2, v2, vo1]

    ! Add these to the mesh%RefStack (if they're not in there already)
    IF (mesh%RefMap(t1)==0) THEN
      mesh%RefMap(t1) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = t1
    END IF
    IF (mesh%RefMap(t2)==0) THEN
      mesh%RefMap(t2) = 1
      mesh%RefStackN = mesh%RefStackN+1
      mesh%RefStack(mesh%RefStackN) = t2
    END IF

    ! == Update the triangle connectivity matrix
    ! t1nv1
    IF (t1nv1>0) THEN
      DO n = 1, 3
        IF (mesh%TriC(t1nv1,n)==t1) mesh%TriC(t1nv1,n) = t2
      END DO
    END IF
    ! t1nv2: nothing changes
    ! t2nv1: nothing changes
    ! t2nv2
    IF (t2nv2>0) THEN
      DO n = 1, 3
        IF (mesh%TriC(t2nv2,n)==t2) mesh%TriC(t2nv2,n) = t1
      END DO
    END IF
    ! The two new triangles
    mesh%TriC(t1,:) = [t2nv2, t2, t1nv2]
    mesh%TriC(t2,:) = [t1nv1, t1, t2nv1]

    ! == Update inverse triangle matrix
    ! v1
    DO n = 1, mesh%niTri(v1)
      IF (mesh%iTri(v1,n)==t2) THEN
        mesh%iTri(v1,:) = [mesh%iTri(v1,1:n-1), mesh%iTri(v1,n+1:mesh%nC_mem), 0]
        mesh%niTri(v1) = mesh%niTri(v1)-1
        EXIT
      END IF
    END DO
    ! v2
    DO n = 1, mesh%niTri(v2)
      IF (mesh%iTri(v2,n)==t1) THEN
        mesh%iTri(v2,:) = [mesh%iTri(v2,1:n-1), mesh%iTri(v2,n+1:mesh%nC_mem), 0]
        mesh%niTri(v2) = mesh%niTri(v2)-1
        EXIT
      END IF
    END DO
    ! vo1
    DO n = 1, mesh%niTri(vo1)
      IF (mesh%iTri(vo1,n)==t1) THEN
        mesh%iTri(vo1,:) = [mesh%iTri(vo1,1:n), t2, mesh%iTri(vo1,n+1:mesh%nC_mem-1)]
        mesh%niTri(vo1) = mesh%niTri(vo1)+1
        EXIT
      END IF
    END DO
    ! vo2
    DO n = 1, mesh%niTri(vo2)
      IF (mesh%iTri(vo2,n)==t2) THEN
        mesh%iTri(vo2,:) = [mesh%iTri(vo2,1:n), t1, mesh%iTri(vo2,n+1:mesh%nC_mem-1)]
        mesh%niTri(vo2) = mesh%niTri(vo2)+1
        EXIT
      END IF
    END DO

    ! == Update vertex connectivity matrix
    ! v1
    DO n = 1, mesh%nC(v1)
      IF (mesh%C(v1,n)==v2) THEN
        mesh%C(v1,:) = [mesh%C(v1,1:n-1), mesh%C(v1,n+1:mesh%nC_mem), 0]
        mesh%nC(v1) = mesh%nC(v1)-1
        EXIT
      END IF
    END DO
    ! v2
    DO n = 1, mesh%nC(v2)
      IF (mesh%C(v2,n)==v1) THEN
        mesh%C(v2,:) = [mesh%C(v2,1:n-1), mesh%C(v2,n+1:mesh%nC_mem), 0]
        mesh%nC(v2) = mesh%nC(v2)-1
        EXIT
      END IF
    END DO
    ! vo1
    DO n = 1, mesh%nC(vo1)
      IF (mesh%C(vo1,n)==v1) THEN
        mesh%C(vo1,:) = [mesh%C(vo1,1:n), vo2, mesh%C(vo1,n+1:mesh%nC_mem-1)]
        mesh%nC(vo1) = mesh%nC(vo1)+1
        EXIT
      END IF
    END DO
    ! vo2
    DO n = 1, mesh%nC(vo2)
      IF (mesh%C(vo2,n)==v2) THEN
        mesh%C(vo2,:) = [mesh%C(vo2,1:n), vo1, mesh%C(vo2,n+1:mesh%nC_mem-1)]
        mesh%nC(vo2) = mesh%nC(vo2)+1
        EXIT
      END IF
    END DO

    ! == Update triangle circumcenters
   CALL update_triangle_circumcenter( mesh, t1)
   CALL update_triangle_circumcenter( mesh, t2)

    ! == Update triangle edge indices
    mesh%Tri_edge_index(t1) = 0
    mesh%Tri_edge_index(t2) = 0
    edgevals       = [0, 1, 3, 5, 7]
    edgevals_left  = [0, 8, 2, 4, 6]
    edgevals_right = [0, 2, 4, 6, 8]
    DO e = 1, 5
     contains_edge = 0
     DO n = 1, 3
      IF (mesh%edge_index(mesh%Tri(t1,n)) == edgevals(e) .OR. &
          mesh%edge_index(mesh%Tri(t1,n)) == edgevals_left(e) .OR. &
          mesh%edge_index(mesh%Tri(t1,n)) == edgevals_right(e)) THEN
       contains_edge = contains_edge+1
      END IF
     END DO
     IF (contains_edge==2) mesh%Tri_edge_index(t1) = edgevals(e);

     contains_edge = 0
     DO n = 1, 3
      IF (mesh%edge_index(mesh%Tri(t2,n)) == edgevals(e) .OR. &
          mesh%edge_index(mesh%Tri(t2,n)) == edgevals_left(e) .OR. &
          mesh%edge_index(mesh%Tri(t2,n)) == edgevals_right(e)) THEN
       contains_edge = contains_edge+1
      END IF
     END DO
     IF (contains_edge==2) mesh%Tri_edge_index(t2) = edgevals(e);
    END DO ! Do e = 1, 5

    ! == Remove current triangle pair from flip list, add 4 new ones
    mesh%Triflip(1:mesh%nTri-1,:) = mesh%Triflip(2:mesh%nTri,:)
    nf = nf-1

    IF (t1nv1>0) THEN
      mesh%Triflip(2:mesh%nTri,:) = mesh%Triflip(1:mesh%nTri-1,:)
      mesh%Triflip(1,:) = [t1nv1, t2]
      nf = nf+1
    END IF
    IF (t1nv2>0) THEN
      mesh%Triflip(2:mesh%nTri,:) = mesh%Triflip(1:mesh%nTri-1,:)
      mesh%Triflip(1,:) = [t1nv2, t1]
      nf = nf+1
    END IF
    IF (t2nv1>0) THEN
      mesh%Triflip(2:mesh%nTri,:) = mesh%Triflip(1:mesh%nTri-1,:)
      mesh%Triflip(1,:) = [t2nv1, t2]
      nf = nf+1
    END IF
    IF (t2nv2>0) THEN
      mesh%Triflip(2:mesh%nTri,:) = mesh%Triflip(1:mesh%nTri-1,:)
      mesh%Triflip(1,:) = [t2nv2, t1]
      nf = nf+1
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE flip_triangle_pairs
  SUBROUTINE need_flipping( mesh, t1, t2, isso, vo1, vo2, v1, v2, t1nv1, t1nv2, t2nv1, t2nv2)
    ! == First, determine general info
    ! Shared vertices v1 and v2 (sorted clockwise in t1), non-shared
    ! vertices vo1 and vo2, neigbours to t1 t1nv1 (across from v1) and
    ! t1nv2 (across from v2) and neighbours to t2 t2nv1 and t2nv2 (idem)    
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: t1, t2
    LOGICAL,                             INTENT(OUT)   :: isso
    INTEGER,                             INTENT(OUT)   :: vo1, vo2, v1, v2, t1nv1, t1nv2, t2nv1, t2nv2

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'need_flipping'
    REAL(dp), DIMENSION(2)                             :: p, q, r, s
    INTEGER                                            :: n, n1, n2, nnext
    LOGICAL                                            :: isint2, SwitchThem
    
    ! Add routine to path
    CALL init_routine( routine_name)

    v1  = 0
    v2  = 0
    vo1 = 0
    DO n1 = 1, 3
      isint2 = .FALSE.
      DO n2 = 1, 3
        IF (mesh%Tri(t1,n1)==mesh%Tri(t2,n2)) THEN
          isint2 = .TRUE.
          IF (v1==0) THEN
            v1 = mesh%Tri(t1,n1)
          ELSE
            v2 = mesh%Tri(t1,n1)
          END IF
        END IF
      END DO ! DO n2 = 1, 3
      IF (.NOT. isint2) vo1 = mesh%Tri(t1,n1)
    END DO ! DO n1 = 1, 3
    DO n = 1, 3
      IF (mesh%Tri(t2,n)/=v1 .AND. mesh%Tri(t2,n)/=v2) vo2 = mesh%Tri(t2,n)
    END DO ! DO n = 1, 3

    ! Order v1 and v2 anticlockwise in triangle t1
    SwitchThem = .FALSE.
    DO n = 1, 3
      nnext=n+1
      IF (nnext==4) nnext=1
      IF (mesh%Tri(t1,n)==v1) THEN
        IF (mesh%Tri(t1,nnext)/=v2) SwitchThem = .TRUE.
      END IF
    END DO ! DO n = 1, 3
    IF (SwitchThem) THEN
      v1 = v1+v2
      v2 = v1-v2
      v1 = v1-v2
    END IF

    ! == Find neighbour triangles
    DO n = 1, 3
      IF (mesh%Tri(t1,n)==v1) t1nv1 = mesh%TriC(t1,n)
      IF (mesh%Tri(t1,n)==v2) t1nv2 = mesh%TriC(t1,n)
      IF (mesh%Tri(t2,n)==v1) t2nv1 = mesh%TriC(t2,n)
      IF (mesh%Tri(t2,n)==v2) t2nv2 = mesh%TriC(t2,n)
    END DO

    ! == Determine if triangle pair t1,t2 requires flipping
    isso = .FALSE.
!    IF (norm2(mesh%V(vo2,:) - mesh%Tricc(t1,:)) < norm2(mesh%V(mesh%Tri(t1,1),:) - mesh%Tricc(t1,:)) - mesh%tol_dist .OR. &
!        norm2(mesh%V(vo1,:) - mesh%Tricc(t2,:)) < norm2(mesh%V(mesh%Tri(t2,1),:) - mesh%Tricc(t2,:)) - mesh%tol_dist) THEN
    IF (norm2(mesh%V(vo2,:) - mesh%Tricc(t1,:)) < norm2(mesh%V(mesh%Tri(t1,1),:) - mesh%Tricc(t1,:)) .OR. &
        norm2(mesh%V(vo1,:) - mesh%Tricc(t2,:)) < norm2(mesh%V(mesh%Tri(t2,1),:) - mesh%Tricc(t2,:))) THEN
      isso = .TRUE.
    END IF

    ! If the outer angle at v1 or v2 is concave, don't flip.
    ! Check this by checking if v1 lies inside the triangle
    ! [vo2,vo2,v2], or the other way round.
    p = mesh%V(vo1,:)
    q = mesh%V(vo2,:)
    r = mesh%V(v1,:)
    s = mesh%V(v2,:)
    IF  (is_in_triangle(p, q, s, r) .OR. &
         is_in_triangle(p, q, r, s)) THEN
      isso = .FALSE.
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE need_flipping
  
  SUBROUTINE move_vertex( mesh, vi, p)
    ! Move vertex vi of the mesh to point p
    
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    INTEGER,                         INTENT(IN)        :: vi
    REAL(dp), DIMENSION(2),          INTENT(IN)        :: p
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'move_vertex'
    INTEGER                                            :: iti, ti, t1, t2, n, nf
    LOGICAL                                            :: did_flip, did_flip_pair
    
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Move the vertex
    mesh%V( vi,:) = p
    
    ! Update surrounding triangle circumcentres
    DO iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      CALL update_triangle_circumcenter( mesh, ti)
    END DO
    
    ! Update triangulation
    did_flip = .TRUE.
    DO WHILE (did_flip)
    
      did_flip = .FALSE.
    
      mesh%Triflip = 0
      nf           = 0
      
      DO iti = 1, mesh%niTri( vi)
        t1 = mesh%iTri( vi,iti)
        DO n = 1, 3
          t2 = mesh%TriC( t1,n)
          IF (t2 > 0) THEN
            nf = nf + 1
            mesh%Triflip( nf,:) = [t1,t2]
          END IF
        END DO
      END DO
      
      DO WHILE (nf > 0)
        CALL flip_triangle_pairs( mesh, nf, did_flip_pair)
        IF (did_flip_pair) did_flip = .TRUE.
      END DO
      
    END DO ! DO WHILE (did_flip)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE move_vertex

END MODULE mesh_Delaunay_module
