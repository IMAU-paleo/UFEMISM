MODULE mesh_smooth_module
  ! Routines for applying a 2D Gaussian smoothing filter to a data field on a mesh.

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
  USE data_types_module,           ONLY: type_mesh
  USE mesh_help_functions_module,  ONLY: LineFromPoints, LineLineIntersection

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'

  CONTAINS
    
  ! == 2D Gaussian smoothing
  SUBROUTINE smooth_Gaussian(mesh, r_smooth, d, d_smooth)
    ! Smooth the data d on the mesh using a Gaussian filter with standard deviation r_smooth

    ! In/output variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh       ! The mesh
    REAL(dp),                   INTENT(IN)        :: r_smooth   ! The smoothing radius (m)
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d          ! The input data to be smoothed
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: d_smooth   ! The smoothed output data
    
    ! Local variables
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: v_int
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: w_int
    INTEGER                                       :: ni, vi, i_mid
    REAL(dp), DIMENSION(:), POINTER               :: d_smooth_x
    INTEGER                                       :: wd_smooth_x
    
    ! Initialise with unsmoothed values
    d_smooth(mesh%v1:mesh%v2) = d(mesh%v1:mesh%v2)
    CALL sync
    
    ! Allocate shared memory for the x-smoothed data
    CALL allocate_shared_memory_dp_1D( mesh%nV, d_smooth_x, wd_smooth_x)
    d_smooth_x(mesh%v1:mesh%v2) = d(mesh%v1:mesh%v2)
    CALL sync
    
    ! Allocate (non-shared) memory for the transects v_int and w_int
    ni    =               CEILING(20._dp*r_smooth / mesh%resolution_min)
    i_mid =  CEILING(REAL(CEILING(20._dp*r_smooth / mesh%resolution_min))/2._dp)
    ALLOCATE(v_int(ni))
    ALLOCATE(w_int(ni))
    v_int = 0._dp
    w_int = 0._dp
    
    ! Actual smoothing
    ! ================
    ! Smooth along the X dimension
    DO vi = mesh%v1, mesh%v2
      IF (mesh%edge_index(vi)>0) CYCLE
      CALL subsmooth_x(mesh, d, r_smooth, vi, i_mid, v_int, w_int, d_smooth_x(vi))
    END DO
    CALL sync
    ! Smooth along the Y dimension
    DO vi = mesh%v1, mesh%v2
      IF (mesh%edge_index(vi)>0) CYCLE
      CALL subsmooth_y(mesh, d_smooth_x, r_smooth, vi, i_mid, v_int, w_int, d_smooth(vi))
    END DO
    CALL sync
    
    ! Deallocate memory for the transects v_int and w_int
    DEALLOCATE(v_int)
    DEALLOCATE(w_int)
    
    ! Deallocate shared memory for the x-smoothed data
    CALL deallocate_shared_memory( wd_smooth_x)
    NULLIFY( d_smooth_x)
    
  END SUBROUTINE smooth_Gaussian
  SUBROUTINE subsmooth_x(mesh, d, r_smooth, vi, i_mid, v_int, w_int, d_smooth)
    ! Gaussian smoothing along the horizontal dimension.
    ! Determines profile along an X transect in both directions, then integrates
    ! over that to determine the smoothed value using a normal distribution G(x)dx
    ! d_smooth = int(d_profile * G(x)dx) / int(G(x)dx)

    ! In/output variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh       ! The mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d          ! The data to be smoothed
    REAL(dp),                   INTENT(IN)        :: r_smooth   ! The smoothing radius (m)
    INTEGER,                    INTENT(IN)        :: vi         ! The inde of the vertex on which we're calculating the smoothed value
    INTEGER,                    INTENT(IN)        :: i_mid      ! Midpoint of v_int and w_int (since those are preallocated in the master routine)
    REAL(dp), DIMENSION(:),     INTENT(INOUT)     :: v_int      ! List of data values along the smoothing dimension
    REAL(dp), DIMENSION(:),     INTENT(INOUT)     :: w_int      ! List of x    values along the smoothing dimension
    REAL(dp),                   INTENT(OUT)       :: d_smooth   ! Final smoothed value of d on vi
    
    ! Local variables
    INTEGER                                       :: i_min, i_max, i_int, n, n1
    INTEGER                                       :: ti, iti, ti_prev, ti_next
    INTEGER                                       :: v1_prev, v2_prev, v1_next, v2_next, v1_next_a, v1_next_b, v2_next_a, v2_next_b
    INTEGER                                       :: ci1, ci2
    REAL(dp)                                      :: xn1, xn2, yn1, yn2
    INTEGER                                       :: quad1, quad2, ncontains
    REAL(dp)                                      :: la1, lb1, lc1, la2, lb2, lc2
    REAL(dp), DIMENSION(2)                        :: p, q, llis_prev, llis_next, llis_next_a, llis_next_b
    REAL(dp)                                      :: w_v1, w_v2, v_llis, v_llis_next, w_v1_next, w_v2_next
    REAL(dp)                                      :: x, dx, dw, sum_vw, sum_w

    ! Initialise output at zero
    d_smooth     = 0._dp
    v_int        = 0._dp
    w_int        = 0._dp
    v_int(i_mid) = d(vi)   ! Central value at the vertex we're interested in
    
    i_min = i_mid-1
    i_max = i_mid+1
    
    v1_prev = 0
    v2_prev = 0
    
    ! =====  Eastward transect

    ! =====  Find out which triangle to start in
    ti_prev = 0
    DO ci1 = 1, mesh%nC(vi)
      ci2 = ci1+1
      IF (ci2==mesh%nC(vi)+1) ci2 = 1

      v1_prev = mesh%C(vi,ci1)
      v2_prev = mesh%C(vi,ci2)

      ! Find x and y coordinates of both neighbour vertices
      xn1 = mesh%V(v1_prev,1) - mesh%V(vi,1)
      yn1 = mesh%V(v1_prev,2) - mesh%V(vi,2)
      xn2 = mesh%V(v2_prev,1) - mesh%V(vi,1)
      yn2 = mesh%V(v2_prev,2) - mesh%V(vi,2)

      ! Find quadrants of both neighbour vertices
      IF (xn1>0) THEN
        IF (yn1>0) THEN
          quad1 = 1
        ELSE
          quad1 = 4
        END IF
      ELSE
        IF (yn1>0) THEN
          quad1 = 2
        ELSE
          quad1 = 3
        END IF
      END IF
      IF (xn2>0) THEN
        IF (yn2>0) THEN
          quad2 = 1
        ELSE
          quad2 = 4
        END IF
      ELSE
        IF (yn2>0) THEN
          quad2 = 2
        ELSE
          quad2 = 3
        END IF
      END IF

      ! Find which triangle contaings both these neighbours
      IF ((quad1==4 .AND. quad2==1) .OR. (quad1==3 .AND. quad2==1) .OR. (quad1==4 .AND. quad2==2)) THEN
        DO iti = 1,mesh%niTri(vi)
          ti = mesh%iTri(vi,iti)
          ncontains = 0
          DO n=1,3
            IF (mesh%Tri(ti,n)==v1_prev .OR. mesh%Tri(ti,n)==v2_prev) THEN
              ncontains = ncontains+1
            END IF
          END DO
          IF (ncontains==2) THEN
            ti_prev = ti
            EXIT
          END IF
        END DO
      END IF

      IF (ti_prev > 0) EXIT

    END DO ! DO ci1 = 1, mesh%nC(vi)

    ! =====  Move eastward through the mesh

    ! The eastward line
    la1 = 0._dp
    lb1 = 1._dp
    lc1 = mesh%V(vi,2)

    ! Start at first line intersection
    p = mesh%V(v1_prev,:)
    q = mesh%V(v2_prev,:)
    CALL LineFromPoints(p,q,la2,lb2,lc2)
    CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_prev);

    ! Add value at that intersection to the integral
    w_v1   = (mesh%V(v2_prev,2) - llis_prev(2)) / (mesh%V(v2_prev,2) - mesh%V(v1_prev,2))
    w_v2   = 1._dp - w_v1
    v_llis = w_v1 * d(v1_prev) + w_v2 * d(v2_prev)

    v_int(i_mid+1) = v_llis
    w_int(i_mid+1) = llis_prev(1) - mesh%V(vi,1)

    ti_next   = ti_prev
    v1_next   = v1_prev
    v2_next   = v2_prev
    llis_next = llis_prev

    i_int = i_mid+2

    DO WHILE (llis_next(1) < min(mesh%V(vi,1) + 3*r_smooth, mesh%xmax) .AND. (mesh%Tri_edge_index(ti_prev)==0))

      ! Find the next triangle
      DO iti = 1, mesh%niTri(v1_prev)
        ti_next = mesh%iTri(v1_prev,iti)
        IF (ti_next==ti_prev) CYCLE
        ncontains = 0
        DO n=1,3
          IF (mesh%Tri(ti_next,n)==v1_prev .OR. mesh%Tri(ti_next,n)==v2_prev) THEN
            ncontains = ncontains+1
          END IF
        END DO
        IF (ncontains==2) EXIT
      END DO

      ! Find the two other line intersections
      v1_next   = 0
      v1_next_a = 0
      v1_next_b = 0
      v2_next   = 0
      v2_next_a = 0
      v2_next_b = 0
      
      DO n=1,3
        n1 = n+1
        IF (n1==4) n1=1
        v1_next = mesh%Tri(ti_next,n)
        v2_next = mesh%Tri(ti_next,n1)
        p = mesh%V(v1_next,:)
        q = mesh%V(v2_next,:)
        CALL LineFromPoints(p,q,la2,lb2,lc2)
        IF (v1_next==v1_prev) THEN
          CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_next_a)
          v1_next_a = v1_next
          v2_next_a = v2_next
        ELSEIF (v2_next==v2_prev) THEN
          CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_next_b)
          v1_next_b = v1_next
          v2_next_b = v2_next
        END IF
      END DO

      ! Determine which if these is correct
      IF     (llis_next_a(1) < llis_prev(1)) THEN
        v1_next =     v1_next_b
        v2_next =     v2_next_b
        llis_next = llis_next_b
      ELSEIF (llis_next_b(1) < llis_prev(1)) THEN
        v1_next =     v1_next_a;
        v2_next =     v2_next_a;
        llis_next = llis_next_a;
      ELSEIF (llis_next_b(1) > llis_next_a(1)) THEN
        v1_next =     v1_next_a
        v2_next =     v2_next_a
        llis_next = llis_next_a
      ELSE
        v1_next =     v1_next_b
        v2_next =     v2_next_b
        llis_next = llis_next_b
      END IF

      ! Calculate value and weight at new intersection
      w_v1_next    = (mesh%V(v2_next,2) - llis_next(2)) / (mesh%V(v2_next,2) - mesh%V(v1_next,2))
      w_v2_next    = 1._dp - w_v1_next
      v_llis_next  = w_v1_next * d(v1_next) + w_v2_next * d(v2_next)

      v_int(i_int) = v_llis_next
      w_int(i_int) = llis_next(1)-mesh%V(vi,1)
      i_int = i_int+1
      i_max = i_max+1

      ! Cycle triangle, vertices and llis
      ti_prev   = ti_next
      v1_prev   = v1_next
      v2_prev   = v2_next
      llis_prev = llis_next

    END DO ! DO WHILE (llis_next(1) < min(mesh%V(vi,1) + 3*r_smooth, mesh%xmax) .AND. (mesh%Tri_edge_index(ti_prev)==0))

    ! =====  Westward transect

    ! =====  Find out which triangle to start in
    ti_prev = 0
    DO ci1 = 1,mesh%nC(vi)
      ci2 = ci1+1
      IF (ci2==mesh%nC(vi)+1) ci2 = 1

      v1_prev = mesh%C(vi,ci1)
      v2_prev = mesh%C(vi,ci2)

      ! Find x and y coordinates of both neighbour vertices
      xn1 = mesh%V(v1_prev,1) - mesh%V(vi,1)
      yn1 = mesh%V(v1_prev,2) - mesh%V(vi,2)
      xn2 = mesh%V(v2_prev,1) - mesh%V(vi,1)
      yn2 = mesh%V(v2_prev,2) - mesh%V(vi,2)

      ! Find quadrants of both neighbour vertices
      IF (xn1>0) THEN
        IF (yn1>0) THEN
          quad1 = 1
        ELSE
          quad1 = 4
        END IF
      ELSE
        IF (yn1>0) THEN
          quad1 = 2
        ELSE
          quad1 = 3
        END IF
      END IF
      IF (xn2>0) THEN
        IF (yn2>0) THEN
          quad2 = 1
        ELSE
          quad2 = 4
        END IF
      ELSE
        IF (yn2>0) THEN
          quad2 = 2
        ELSE
          quad2 = 3
        END IF
      END IF

      ! Find which triangle contaings both these neighbours
      IF ((quad1==2 .AND. quad2==3) .OR. (quad1==1 .AND. quad2==3) .OR. (quad1==2 .AND. quad2==4)) THEN
        DO iti = 1,mesh%niTri(vi)
          ti = mesh%iTri(vi,iti)
          ncontains = 0
          DO n=1,3
            IF (mesh%Tri(ti,n)==v1_prev .OR. mesh%Tri(ti,n)==v2_prev) THEN
              ncontains = ncontains+1
            END IF
          END DO
          IF (ncontains==2) THEN
            ti_prev = ti
            EXIT
          END IF
        END DO
      END IF

      IF (ti_prev > 0) EXIT

    END DO ! DO ci1 = 1,mesh%nC(vi)

    ! =====  Move westward through the mesh% Skip the value at vi since we already have that one

    ! The eastward line
    la1 = 0._dp
    lb1 = 1._dp
    lc1 = mesh%V(vi,2)

    ! Start at first line intersection
    p = mesh%V(v1_prev,:)
    q = mesh%V(v2_prev,:)
    CALL LineFromPoints(p,q,la2,lb2,lc2)
    CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_prev)

    ! Add value at that intersection to the integral
    w_v1   = (mesh%V(v2_prev,2) - llis_prev(2)) / (mesh%V(v2_prev,2) - mesh%V(v1_prev,2))
    w_v2   = 1._dp - w_v1 
    v_llis = w_v1 * d(v1_prev) + w_v2 * d(v2_prev)

    v_int(i_mid-1) = v_llis
    w_int(i_mid-1) = llis_prev(1) - mesh%V(vi,1)

    ti_next   = ti_prev
    v1_next   = v1_prev
    v2_next   = v2_prev
    llis_next = llis_prev

    i_int = i_mid-2

    DO WHILE (llis_next(1) > max(mesh%V(vi,1) - 3*r_smooth, mesh%xmin) .AND. (mesh%Tri_edge_index(ti_prev)==0))

      ! Find the next triangle
      DO iti = 1,mesh%niTri(v1_prev)
        ti_next = mesh%iTri(v1_prev,iti)
        IF (ti_next==ti_prev) CYCLE
        ncontains = 0
        DO n=1,3
          IF (mesh%Tri(ti_next,n)==v1_prev .OR. mesh%Tri(ti_next,n)==v2_prev) THEN
            ncontains = ncontains+1
          END IF
        END DO
        IF (ncontains==2) EXIT
      END DO

      ! Find the two other line intersections
      v1_next   = 0
      v1_next_a = 0
      v1_next_b = 0
      v2_next   = 0
      v2_next_a = 0
      v2_next_b = 0
      DO n=1,3
        n1 = n+1
        IF (n1==4) n1=1
        v1_next = mesh%Tri(ti_next,n)
        v2_next = mesh%Tri(ti_next,n1)
        p = mesh%V(v1_next,:)
        q = mesh%V(v2_next,:)
        CALL LineFromPoints(p,q,la2,lb2,lc2)
        IF (v1_next==v1_prev) THEN
          CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_next_a)
          v1_next_a = v1_next
          v2_next_a = v2_next
        ELSEIF (v2_next==v2_prev) THEN
          CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_next_b)
          v1_next_b = v1_next
          v2_next_b = v2_next
        END IF
      END DO

      ! Determine which if these is correct
      IF     (llis_next_a(1) > llis_prev(1)) THEN
        v1_next =     v1_next_b
        v2_next =     v2_next_b
        llis_next = llis_next_b
      ELSEIF (llis_next_b(1) > llis_prev(1)) THEN
        v1_next =     v1_next_a
        v2_next =     v2_next_a
        llis_next = llis_next_a
      ELSEIF (llis_next_b(1) < llis_next_a(1)) THEN
        v1_next =     v1_next_a
        v2_next =     v2_next_a
        llis_next = llis_next_a
      ELSE
        v1_next =     v1_next_b
        v2_next =     v2_next_b
        llis_next = llis_next_b
      END IF

      ! Calculate value and weight at new intersection
      w_v1_next    = (mesh%V(v2_next,2) - llis_next(2)) / (mesh%V(v2_next,2) - mesh%V(v1_next,2))
      w_v2_next    = 1._dp - w_v1_next
      v_llis_next  = w_v1_next * d(v1_next) + w_v2_next * d(v2_next)

      v_int(i_int) = v_llis_next
      w_int(i_int) = llis_next(1)-mesh%V(vi,1)
      i_int = i_int-1
      i_min = i_min-1

      ! Cycle triangle, vertices and llis
      ti_prev   = ti_next
      v1_prev   = v1_next
      v2_prev   = v2_next
      llis_prev = llis_next

    END DO !  DO WHILE (llis_next(1) > max(mesh%V(vi,1) - 3*r_smooth, mesh%xmin) .AND. (mesh%Tri_edge_index(ti_prev)==0))

    ! =====  Integrate over transect
    sum_vw = 0._dp
    sum_w  = 0._dp

    DO i_int = i_min,i_max
      x  = w_int(i_int)

      IF     (i_int==i_min) THEN
        dx = w_int(i_int+1)-w_int(i_int)
      ELSEIF (i_int==i_max) THEN
        dx = w_int(i_int)-w_int(i_int-1)
      ELSE
        dx = (w_int(i_int+1)-w_int(i_int-1))/2._dp
      END IF

      dw = dx * EXP(-(x**2)/(2*r_smooth**2));

      sum_vw = sum_vw + (v_int(i_int) * dw)
      sum_w  = sum_w  + dw
    END DO

    ! Calculate final smoothed value
    d_smooth = sum_vw / sum_w

  END SUBROUTINE subsmooth_x
  SUBROUTINE subsmooth_y(mesh, d, r_smooth, vi, i_mid, v_int, w_int, d_smooth)
    ! Gaussian smoothing along the vertical dimension.
    ! Determines profile along an Y transect in both directions, then integrates
    ! over that to determine the smoothed value using a normal distribution G(y)dy
    ! d_smooth = int(d_profile * G(y)dy) / int(G(y)dy)

    ! In/output variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh       ! The mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: d          ! The data to be smoothed
    REAL(dp),                   INTENT(IN)        :: r_smooth   ! The smoothing radius (m)
    INTEGER,                    INTENT(IN)        :: vi         ! The inde of the vertex on which we're calculating the smoothed value
    INTEGER,                    INTENT(IN)        :: i_mid      ! Midpoint of v_int and w_int (since those are preallocated in the master routine)
    REAL(dp), DIMENSION(:),     INTENT(INOUT)     :: v_int      ! List of data values along the smoothing dimension
    REAL(dp), DIMENSION(:),     INTENT(INOUT)     :: w_int      ! List of x    values along the smoothing dimension
    REAL(dp),                   INTENT(OUT)       :: d_smooth   ! Final smoothed value of d on vi
    
    ! Local variables
    INTEGER                                       :: i_min, i_max, i_int, n, n1
    INTEGER                                       :: ti, iti, ti_prev, ti_next
    INTEGER                                       :: v1_prev, v2_prev, v1_next, v2_next, v1_next_a, v1_next_b, v2_next_a, v2_next_b
    INTEGER                                       :: ci1, ci2
    REAL(dp)                                      :: xn1, xn2, yn1, yn2
    INTEGER                                       :: quad1, quad2, ncontains
    REAL(dp)                                      :: la1, lb1, lc1, la2, lb2, lc2
    REAL(dp), DIMENSION(2)                        :: p, q, llis_prev, llis_next, llis_next_a, llis_next_b
    REAL(dp)                                      :: w_v1, w_v2, v_llis, v_llis_next, w_v1_next, w_v2_next
    REAL(dp)                                      :: x, dx, dw, sum_vw, sum_w

    ! Initialise output at zero
    d_smooth     = 0._dp
    v_int        = 0._dp
    w_int        = 0._dp
    v_int(i_mid) = d(vi)   ! Central value at the vertex we're interested in
    
    i_min = i_mid-1
    i_max = i_mid+1
    
    v1_prev = 0
    v2_prev = 0

    ! =====  Northward transect

    ! =====  Find out which triangle to start in
    ti_prev = 0
    DO ci1 = 1, mesh%nC(vi)
      ci2 = ci1+1
      IF (ci2==mesh%nC(vi)+1) ci2 = 1

      v1_prev = mesh%C(vi,ci1)
      v2_prev = mesh%C(vi,ci2)

      ! Find x and y coordinates of both neighbour vertices
      xn1 = mesh%V(v1_prev,1) - mesh%V(vi,1)
      yn1 = mesh%V(v1_prev,2) - mesh%V(vi,2)
      xn2 = mesh%V(v2_prev,1) - mesh%V(vi,1)
      yn2 = mesh%V(v2_prev,2) - mesh%V(vi,2)

      ! Find quadrants of both neighbour vertices
      IF (xn1>0) THEN
        IF (yn1>0) THEN
          quad1 = 1
        ELSE
          quad1 = 4
        END IF
      ELSE
        IF (yn1>0) THEN
          quad1 = 2
        ELSE
          quad1 = 3
        END IF
      END IF
      IF (xn2>0) THEN
        IF (yn2>0) THEN
          quad2 = 1
        ELSE
          quad2 = 4
        END IF
      ELSE
        IF (yn2>0) THEN
          quad2 = 2
        ELSE
          quad2 = 3
        END IF
      END IF

      ! Find which triangle contaings both these neighbours
      IF ((quad1==1 .AND. quad2==2) .OR. (quad1==4 .AND. quad2==2) .OR. (quad1==1 .AND. quad2==3)) THEN
        DO iti = 1,mesh%niTri(vi)
          ti = mesh%iTri(vi,iti)
          ncontains = 0
          DO n=1,3
            IF (mesh%Tri(ti,n)==v1_prev .OR. mesh%Tri(ti,n)==v2_prev) THEN
              ncontains = ncontains+1
            END IF
          END DO
          IF (ncontains==2) THEN
            ti_prev = ti
            EXIT
          END IF
        END DO
      END IF

      IF (ti_prev > 0) EXIT

    END DO ! DO ci1 = 1, mesh%nC(vi)

    ! =====  Move eastward through the mesh

    ! The northward line
    la1 = 1._dp
    lb1 = 0._dp
    lc1 = mesh%V(vi,1)

    ! Start at first line intersection
    p = mesh%V(v1_prev,:)
    q = mesh%V(v2_prev,:)
    CALL LineFromPoints(p,q,la2,lb2,lc2)
    CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_prev);

    ! Add value at that intersection to the integral
    w_v1   = (mesh%V(v2_prev,1) - llis_prev(1)) / (mesh%V(v2_prev,1) - mesh%V(v1_prev,1))
    w_v2   = 1._dp - w_v1
    v_llis = w_v1 * d(v1_prev) + w_v2 * d(v2_prev)

    v_int(i_mid+1) = v_llis
    w_int(i_mid+1) = llis_prev(2) - mesh%V(vi,2)

    ti_next   = ti_prev
    v1_next   = v1_prev
    v2_next   = v2_prev
    llis_next = llis_prev

    i_int = i_mid+2

    DO WHILE (llis_next(2) < min(mesh%V(vi,2) + 3*r_smooth, mesh%ymax) .AND. (mesh%Tri_edge_index(ti_prev)==0))

      ! Find the next triangle
      DO iti = 1, mesh%niTri(v1_prev)
        ti_next = mesh%iTri(v1_prev,iti)
        IF (ti_next==ti_prev) CYCLE
        ncontains = 0
        DO n=1,3
          IF (mesh%Tri(ti_next,n)==v1_prev .OR. mesh%Tri(ti_next,n)==v2_prev) THEN
            ncontains = ncontains+1
          END IF
        END DO
        IF (ncontains==2) EXIT
      END DO

      ! Find the two other line intersections
      v1_next   = 0
      v1_next_a = 0
      v1_next_b = 0
      v2_next   = 0
      v2_next_a = 0
      v2_next_b = 0
      
      DO n=1,3
        n1 = n+1
        IF (n1==4) n1=1
        v1_next = mesh%Tri(ti_next,n)
        v2_next = mesh%Tri(ti_next,n1)
        p = mesh%V(v1_next,:)
        q = mesh%V(v2_next,:)
        CALL LineFromPoints(p,q,la2,lb2,lc2)
        IF (v1_next==v1_prev) THEN
          CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_next_a)
          v1_next_a = v1_next
          v2_next_a = v2_next
        ELSEIF (v2_next==v2_prev) THEN
          CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_next_b)
          v1_next_b = v1_next
          v2_next_b = v2_next
        END IF
      END DO

      ! Determine which if these is correct
      IF     (llis_next_a(2) < llis_prev(2)) THEN
        v1_next =     v1_next_b
        v2_next =     v2_next_b
        llis_next = llis_next_b
      ELSEIF (llis_next_b(2) < llis_prev(2)) THEN
        v1_next =     v1_next_a;
        v2_next =     v2_next_a;
        llis_next = llis_next_a;
      ELSEIF (llis_next_b(2) > llis_next_a(2)) THEN
        v1_next =     v1_next_a
        v2_next =     v2_next_a
        llis_next = llis_next_a
      ELSE
        v1_next =     v1_next_b
        v2_next =     v2_next_b
        llis_next = llis_next_b
      END IF

      ! Calculate value and weight at new intersection
      w_v1_next    = (mesh%V(v2_next,1) - llis_next(1)) / (mesh%V(v2_next,1) - mesh%V(v1_next,1))
      w_v2_next    = 1._dp - w_v1_next
      v_llis_next  = w_v1_next * d(v1_next) + w_v2_next * d(v2_next)

      v_int(i_int) = v_llis_next
      w_int(i_int) = llis_next(2)-mesh%V(vi,2)
      i_int = i_int+1
      i_max = i_max+1

      ! Cycle triangle, vertices and llis
      ti_prev   = ti_next
      v1_prev   = v1_next
      v2_prev   = v2_next
      llis_prev = llis_next

    END DO ! DO WHILE (llis_next(1) < min(mesh%V(vi,1) + 3*r_smooth, mesh%xmax) .AND. (mesh%Tri_edge_index(ti_prev)==0))

    ! =====  Southward transect

    ! =====  Find out which triangle to start in
    ti_prev = 0
    DO ci1 = 1,mesh%nC(vi)
      ci2 = ci1+1
      IF (ci2==mesh%nC(vi)+1) ci2 = 1

      v1_prev = mesh%C(vi,ci1)
      v2_prev = mesh%C(vi,ci2)

      ! Find x and y coordinates of both neighbour vertices
      xn1 = mesh%V(v1_prev,1) - mesh%V(vi,1)
      yn1 = mesh%V(v1_prev,2) - mesh%V(vi,2)
      xn2 = mesh%V(v2_prev,1) - mesh%V(vi,1)
      yn2 = mesh%V(v2_prev,2) - mesh%V(vi,2)

      ! Find quadrants of both neighbour vertices
      IF (xn1>0) THEN
        IF (yn1>0) THEN
          quad1 = 1
        ELSE
          quad1 = 4
        END IF
      ELSE
        IF (yn1>0) THEN
          quad1 = 2
        ELSE
          quad1 = 3
        END IF
      END IF
      IF (xn2>0) THEN
        IF (yn2>0) THEN
          quad2 = 1
        ELSE
          quad2 = 4
        END IF
      ELSE
        IF (yn2>0) THEN
          quad2 = 2
        ELSE
          quad2 = 3
        END IF
      END IF

      ! Find which triangle contaings both these neighbours
      IF ((quad1==3 .AND. quad2==4) .OR. (quad1==2 .AND. quad2==4) .OR. (quad1==3 .AND. quad2==1)) THEN
        DO iti = 1,mesh%niTri(vi)
          ti = mesh%iTri(vi,iti)
          ncontains = 0
          DO n=1,3
            IF (mesh%Tri(ti,n)==v1_prev .OR. mesh%Tri(ti,n)==v2_prev) THEN
              ncontains = ncontains+1
            END IF
          END DO
          IF (ncontains==2) THEN
            ti_prev = ti
            EXIT
          END IF
        END DO
      END IF

      IF (ti_prev > 0) EXIT

    END DO ! DO ci1 = 1,mesh%nC(vi)

    ! =====  Move westward through the mesh% Skip the value at vi since we already have that one

    ! The southward line
    la1 = 1._dp
    lb1 = 0._dp
    lc1 = mesh%V(vi,1)

    ! Start at first line intersection
    p = mesh%V(v1_prev,:)
    q = mesh%V(v2_prev,:)
    CALL LineFromPoints(p,q,la2,lb2,lc2)
    CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_prev)

    ! Add value at that intersection to the integral
    w_v1   = (mesh%V(v2_prev,1) - llis_prev(1)) / (mesh%V(v2_prev,1) - mesh%V(v1_prev,1))
    w_v2   = 1._dp - w_v1 
    v_llis = w_v1 * d(v1_prev) + w_v2 * d(v2_prev)

    v_int(i_mid-1) = v_llis
    w_int(i_mid-1) = llis_prev(2) - mesh%V(vi,2)

    ti_next   = ti_prev
    v1_next   = v1_prev
    v2_next   = v2_prev
    llis_next = llis_prev

    i_int = i_mid-2

    DO WHILE (llis_next(2) > max(mesh%V(vi,2) - 3*r_smooth, mesh%ymin) .AND. (mesh%Tri_edge_index(ti_prev)==0))

      ! Find the next triangle
      DO iti = 1,mesh%niTri(v1_prev)
        ti_next = mesh%iTri(v1_prev,iti)
        IF (ti_next==ti_prev) CYCLE
        ncontains = 0
        DO n=1,3
          IF (mesh%Tri(ti_next,n)==v1_prev .OR. mesh%Tri(ti_next,n)==v2_prev) THEN
            ncontains = ncontains+1
          END IF
        END DO
        IF (ncontains==2) EXIT
      END DO

      ! Find the two other line intersections
      v1_next   = 0
      v1_next_a = 0
      v1_next_b = 0
      v2_next   = 0
      v2_next_a = 0
      v2_next_b = 0
      
      DO n=1,3
        n1 = n+1
        IF (n1==4) n1=1
        v1_next = mesh%Tri(ti_next,n)
        v2_next = mesh%Tri(ti_next,n1)
        p = mesh%V(v1_next,:)
        q = mesh%V(v2_next,:)
        CALL LineFromPoints(p,q,la2,lb2,lc2)
        IF (v1_next==v1_prev) THEN
          CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_next_a)
          v1_next_a = v1_next
          v2_next_a = v2_next
        ELSEIF (v2_next==v2_prev) THEN
          CALL LineLineIntersection(la1,lb1,lc1,la2,lb2,lc2,llis_next_b)
          v1_next_b = v1_next
          v2_next_b = v2_next
        END IF
      END DO

      ! Determine which if these is correct
      IF     (llis_next_a(2) > llis_prev(2)) THEN
        v1_next =     v1_next_b
        v2_next =     v2_next_b
        llis_next = llis_next_b
      ELSEIF (llis_next_b(2) > llis_prev(2)) THEN
        v1_next =     v1_next_a
        v2_next =     v2_next_a
        llis_next = llis_next_a
      ELSEIF (llis_next_b(2) < llis_next_a(2)) THEN
        v1_next =     v1_next_a
        v2_next =     v2_next_a
        llis_next = llis_next_a
      ELSE
        v1_next =     v1_next_b
        v2_next =     v2_next_b
        llis_next = llis_next_b
      END IF

      ! Calculate value and weight at new intersection
      w_v1_next    = (mesh%V(v2_next,1) - llis_next(1)) / (mesh%V(v2_next,1) - mesh%V(v1_next,1))
      w_v2_next    = 1._dp - w_v1_next
      v_llis_next  = w_v1_next * d(v1_next) + w_v2_next * d(v2_next)

      v_int(i_int) = v_llis_next
      w_int(i_int) = llis_next(2)-mesh%V(vi,2)
      i_int = i_int-1
      i_min = i_min-1

      ! Cycle triangle, vertices and llis
      ti_prev   = ti_next
      v1_prev   = v1_next
      v2_prev   = v2_next
      llis_prev = llis_next

    END DO !  DO WHILE (llis_next(1) > max(mesh%V(vi,1) - 3*r_smooth, mesh%xmin) .AND. (mesh%Tri_edge_index(ti_prev)==0))

    ! =====  Integrate over transect
    sum_vw = 0._dp
    sum_w  = 0._dp

    DO i_int = i_min,i_max
      x  = w_int(i_int)

      IF     (i_int==i_min) THEN
        dx = w_int(i_int+1)-w_int(i_int)
      ELSEIF (i_int==i_max) THEN
        dx = w_int(i_int)-w_int(i_int-1)
      ELSE
        dx = (w_int(i_int+1)-w_int(i_int-1))/2._dp
      END IF

      dw = dx * EXP(-(x**2)/(2*r_smooth**2));

      sum_vw = sum_vw + (v_int(i_int) * dw)
      sum_w  = sum_w  + dw
    END DO

    ! Calculate final smoothed value
    d_smooth = sum_vw / sum_w

  END SUBROUTINE subsmooth_y

END MODULE mesh_smooth_module
