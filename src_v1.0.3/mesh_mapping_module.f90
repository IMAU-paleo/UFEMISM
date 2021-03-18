MODULE mesh_mapping_module
  ! Routines for creating mapping arrays and mapping data between meshes and grids

  USE mpi
  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync, allocate_shared_int_0D, allocate_shared_dp_0D, &
                                                    allocate_shared_int_1D, allocate_shared_dp_1D, &
                                                    allocate_shared_int_2D, allocate_shared_dp_2D, &
                                                    allocate_shared_int_3D, allocate_shared_dp_3D, &
                                                    allocate_shared_bool_1D, deallocate_shared, &
                                                    adapt_shared_int_1D,    adapt_shared_dp_1D, &
                                                    adapt_shared_int_2D,    adapt_shared_dp_2D, &
                                                    adapt_shared_int_3D,    adapt_shared_dp_3D, &
                                                    adapt_shared_bool_1D                   
  USE data_types_module,           ONLY: type_mesh, type_remapping, type_remapping_trilin, type_remapping_nearest_neighbour, &
                                         type_remapping_conservative
  USE mesh_help_functions_module,  ONLY: is_in_triangle, find_Voronoi_cell_vertices, find_containing_triangle, find_triangle_area, find_containing_vertex, &
                                         line_integral_xdy, line_integral_mxydx, line_integral_xydy, is_boundary_segment, cross2, &
                                         lies_on_line_segment, segment_intersection, partition_domain_x_balanced, write_mesh_to_text_file, &
                                         line_from_points, line_line_intersection
  USE mesh_derivatives_module,     ONLY: get_mesh_derivatives, get_mesh_derivatives_3D

  IMPLICIT NONE

  CONTAINS

  ! == Subroutines for mapping data between a cartesian grid and the mesh ==
  SUBROUTINE map_cart_to_mesh( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, ny, i1, i2, cartdata, meshdata)
    ! Map a single cartesian data field to the mesh, using the pre-calculated mapping arrays
    
    IMPLICIT NONE

    TYPE(type_mesh),                         INTENT(IN)  :: mesh
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)  :: map_mean_i
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)  :: map_mean_w
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)  :: map_bilin_i
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)  :: map_bilin_w
    INTEGER,                                 INTENT(IN)  :: ny
    INTEGER,                                 INTENT(IN)  :: i1, i2
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)  :: cartdata

    REAL(dp), DIMENSION(:    ),              INTENT(OUT) :: meshdata

    INTEGER                                              :: vi, i, j, il, iu, jl, ju
    REAL(dp)                                             :: wi, wiljl, wilju, wiujl, wiuju

    meshdata( mesh%v1:mesh%v2) = 0._dp
    CALL sync

    ! First the means
    DO i = i1, i2
      DO j = 1, ny
        IF (map_mean_i(i,j)==0) CYCLE
  
        vi = map_mean_i(i,j)
        wi = map_mean_w(i,j)
        meshdata(vi) = meshdata(vi) + cartdata(i,j) * wi
  
      END DO
    END DO ! DO i = i1, i2
    CALL sync

    ! Then the bilinear interpolation
    DO vi = mesh%v1, mesh%v2
      IF (map_bilin_i(vi,1)==0) CYCLE

      il = map_bilin_i( vi,1)
      iu = map_bilin_i( vi,2)
      jl = map_bilin_i( vi,3)
      ju = map_bilin_i( vi,4)

      wiljl = map_bilin_w( vi,1)
      wilju = map_bilin_w( vi,2)
      wiujl = map_bilin_w( vi,3)
      wiuju = map_bilin_w( vi,4)

      meshdata(vi) = cartdata( il,jl) * wiljl + &
                     cartdata( il,ju) * wilju + &
                     cartdata( iu,jl) * wiujl + &
                     cartdata( iu,ju) * wiuju

    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE map_cart_to_mesh
  SUBROUTINE map_cart_to_mesh_3D( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, ny, nz, i1, i2, cartdata, meshdata)
    ! Map a single cartesian data field to the mesh, using the pre-calculated mapping arrays
    
    IMPLICIT NONE

    TYPE(type_mesh),                         INTENT(IN)  :: mesh
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)  :: map_mean_i
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)  :: map_mean_w
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)  :: map_bilin_i
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)  :: map_bilin_w
    INTEGER,                                 INTENT(IN)  :: ny, nz
    INTEGER,                                 INTENT(IN)  :: i1, i2
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)  :: cartdata

    REAL(dp), DIMENSION(:,:  ),              INTENT(OUT) :: meshdata

    INTEGER                                              :: vi, i, j, k, il, iu, jl, ju
    REAL(dp)                                             :: wi, wiljl, wilju, wiujl, wiuju

    meshdata( mesh%v1:mesh%v2,:) = 0._dp
    CALL sync

    ! First the means
    DO i = i1, i2
      DO j = 1, ny
        IF (map_mean_i(i,j)==0) CYCLE
  
        vi = map_mean_i(i,j)
        wi = map_mean_w(i,j)
        
        DO k = 1, nz
          meshdata(vi,k) = meshdata(vi,k) + cartdata(i,j,k) * wi
        END DO
  
      END DO
    END DO ! DO i = i1, i2
    CALL sync

    ! Then the bilinear interpolation
    DO vi = mesh%v1, mesh%v2
      IF (map_bilin_i(vi,1)==0) CYCLE

      il = map_bilin_i(vi,1)
      iu = map_bilin_i(vi,2)
      jl = map_bilin_i(vi,3)
      ju = map_bilin_i(vi,4)

      wiljl = map_bilin_w(vi,1)
      wilju = map_bilin_w(vi,2)
      wiujl = map_bilin_w(vi,3)
      wiuju = map_bilin_w(vi,4)

      DO k = 1, nz
        meshdata(vi,k) = cartdata(il,jl,k) * wiljl + &
                         cartdata(il,ju,k) * wilju + &
                         cartdata(iu,jl,k) * wiujl + &
                         cartdata(iu,ju,k) * wiuju
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

  END SUBROUTINE map_cart_to_mesh_3D
  SUBROUTINE get_cart_to_mesh_map( mesh, x, y, nx, ny, map_mean_i,  map_mean_w,  map_bilin_i,  map_bilin_w, &
                                            wmap_mean_i, wmap_mean_w, wmap_bilin_i, wmap_bilin_w)
    ! For those triangles that contain enough cartesian gridcells for a mean, the index of that triangle and
    ! the relative contribution of each contained gridcell (= 1 / number of contained gridcells) are stored
    ! on the cartesian grid (since each gridcell will only contribute to the mean of one triangle).

    ! For those triangles that are too small, bilinear interpolation between four cartesian gridcells must be
    ! done. The i and j indices of those gridcells, and their relative contributions. are stored on the mesh.
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:  ),   INTENT(IN)        :: x          ! Data x grid
    REAL(dp), DIMENSION(:  ),   INTENT(IN)        :: y          ! Data y grid
    INTEGER,                    INTENT(IN)        :: nx         ! Number of x elements
    INTEGER,                    INTENT(IN)        :: ny         ! Number of y elements

    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(OUT)       :: map_mean_i  ! On cart grid: vertex indices
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)       :: map_mean_w  ! On cart grid: weights
    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(OUT)       :: map_bilin_i ! On mesh     : cart indices (4x i and j))
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)       :: map_bilin_w ! On mesh     : weights      (4x)    
    INTEGER,                    INTENT(OUT)       :: wmap_mean_i, wmap_mean_w, wmap_bilin_i, wmap_bilin_w

    REAL(dp), DIMENSION(2)                        :: p,q,r,s
    INTEGER                                       :: vi, nVor, n, nmax
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Vor
    REAL(dp)                                      :: cxmin,cxmax,cymin,cymax
    INTEGER                                       :: il,iu,jl,ju,i,j,il2,iu2,jl2,ju2
    REAL(dp)                                      :: wil,wiu,wjl,wju

    INTEGER,  DIMENSION(nx,ny)                    :: cart_i, cart_i_temp
    REAL(dp), DIMENSION(nx,ny)                    :: cart_w, cart_w_temp
    INTEGER,  DIMENSION(mesh%nV, 4)               :: mesh_i
    REAL(dp), DIMENSION(mesh%nV, 4)               :: mesh_w
    
    REAL(dp)                                      :: x_range_mesh, y_range_mesh
    
    ! Allocate shared memory for the mapping arrays
    CALL allocate_shared_int_2D(nx,      ny,     map_mean_i,   wmap_mean_i)
    CALL allocate_shared_dp_2D( nx,      ny,     map_mean_w,   wmap_mean_w)
    CALL allocate_shared_int_2D(mesh%nV, 4,      map_bilin_i,  wmap_bilin_i)
    CALL allocate_shared_dp_2D( mesh%nV, 4,      map_bilin_w,  wmap_bilin_w)

    IF (par%master) THEN
      map_mean_i  = 0
      map_mean_w  = 0._dp
      map_bilin_i = 0
      map_bilin_w = 0._dp
    END IF
    CALL sync
    
    cart_i = 0
    cart_w = 0._dp
    mesh_i = 0
    mesh_w = 0._dp

    cart_i_temp = 0
    cart_w_temp = 0._dp
    
    il2 = 1
    iu2 = 1
    jl2 = 1
    ju2 = 1
    
    x_range_mesh = mesh%xmax - mesh%xmin
    y_range_mesh = mesh%ymax - mesh%ymin

    ALLOCATE(Vor(mesh%nC_mem+2,2))

    DO vi = mesh%v1, mesh%v2

      cart_i_temp(il2:iu2,jl2:ju2) = 0
      cart_w_temp(il2:iu2,jl2:ju2) = 0._dp
      nmax        = 0

      ! Find Voronoi cell vertices      
      CALL find_Voronoi_cell_vertices(mesh, vi, Vor, nVor)
      
      il2 = MAX(1,  1 +FLOOR((MINVAL(Vor(:,1))-MINVAL(x)) / (x(2)-x(1))))
      iu2 = MIN(nx, nx-FLOOR((MAXVAL(x)-MAXVAL(Vor(:,1))) / (x(2)-x(1))))
      jl2 = MAX(1,  1 +FLOOR((MINVAL(Vor(:,2))-MINVAL(y)) / (y(2)-y(1))))
      ju2 = MIN(ny, ny-FLOOR((MAXVAL(y)-MAXVAL(Vor(:,2))) / (y(2)-y(1))))

      DO n = 2, nVor

        ! For all the subtriangles of the Voronoi cell, find the Cartesian
        ! gridcells that lie within it.

        ! Find x and y bounds of Voronoi subtriangle
        cxmin = MINVAL([Vor(n,1), Vor(n-1,1), mesh%V(vi,1)])
        cxmax = MAXVAL([Vor(n,1), Vor(n-1,1), mesh%V(vi,1)])
        cymin = MINVAL([Vor(n,2), Vor(n-1,2), mesh%V(vi,2)])
        cymax = MAXVAL([Vor(n,2), Vor(n-1,2), mesh%V(vi,2)])

        ! Find corresponding i and j bounds on Cartesian grid
        il = MAX(1,  1 +FLOOR((cxmin-MINVAL(x)) / (x(2)-x(1))))
        iu = MIN(nx, nx-FLOOR((MAXVAL(x)-cxmax) / (x(2)-x(1))))
        jl = MAX(1,  1 +FLOOR((cymin-MINVAL(y)) / (y(2)-y(1))))
        ju = MIN(ny, ny-FLOOR((MAXVAL(y)-cymax) / (y(2)-y(1))))

        ! Determine which Cartesian grid cells lie within the subtriangle
        p = Vor(n,:)
        q = Vor(n-1,:)
        r = mesh%V(vi,:)
        DO i = il, iu
        DO j = jl, ju
          s = [MIN(mesh%xmax - x_range_mesh/1E6, MAX(mesh%xmin + x_range_mesh/1E6, x(i))), &
               MIN(mesh%ymax - y_range_mesh/1E6, MAX(mesh%ymin + y_range_mesh/1E6, y(j)))]
          IF (is_in_triangle( p, q, r, s)) THEN
            nmax = nmax + 1
            cart_i_temp(i,j) = vi
            cart_w_temp(i,j) = 1._dp
          END IF
        END DO
        END DO

      END DO

      ! Determine if the triangle contains enough elements for a Mean. If not,
      ! use bilinear interpolation instead

      IF (nmax>4) THEN
        ! Mean over Cartesian gridcells
        DO i = il2, iu2
        DO j = jl2, ju2
          IF (cart_i_temp(i,j)==0) CYCLE
          cart_i(i,j) = cart_i_temp(i,j)
          cart_w(i,j) = cart_w_temp(i,j) / REAL(nmax)
        END DO
        END DO
        
      ELSE
        ! Do a bilinear interpolation instead
        il = MAX(1,MIN(nx-1, 1 + FLOOR((mesh%V(vi,1)-MINVAL(x)) / (x(2)-x(1)))))
        iu = il+1
        jl = MAX(1,MIN(ny-1, 1 + FLOOR((mesh%V(vi,2)-MINVAL(y)) / (y(2)-y(1)))))
        ju = jl+1

        wil = (x(iu) - mesh%V(vi,1))/(x(2)-x(1))
        wiu = 1-wil
        wjl = (y(ju) - mesh%V(vi,2))/(y(2)-y(1))
        wju = 1-wjl

        mesh_i(vi,:) = [il, iu, jl, ju]
        mesh_w(vi,:) = [wil*wjl, wil*wju, wiu*wjl, wiu*wju]
      END IF

    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

    DEALLOCATE(Vor)
    
    ! Write to the actual shared array - one by one
    DO n = 0, par%n-1
      IF (n==par%i) THEN
        DO i = 1, nx
        DO j = 1, NY
          IF (cart_i(i,j)==0) CYCLE
          map_mean_i(i,j) = cart_i(i,j)
          map_mean_w(i,j) = cart_w(i,j)
        END DO
        END DO
        DO vi = 1, mesh%nV
          IF (mesh_i(vi,1)==0) CYCLE
          map_bilin_i(vi,:) = mesh_i(vi,:)
          map_bilin_w(vi,:) = mesh_w(vi,:)
        END DO
      END IF
      CALL sync
    END DO

  END SUBROUTINE get_cart_to_mesh_map
  SUBROUTINE get_glob_to_mesh_map( mesh, lat, lon, nlat, nlon, map_mean_i,  map_mean_w,  map_bilin_i,  map_bilin_w, &
                                                              wmap_mean_i, wmap_mean_w, wmap_bilin_i, wmap_bilin_w)
    ! Global GCM data fields have such a low resolution that we only ever do bilinear interpolation.
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:  ),   INTENT(IN)        :: lat         ! Data lat grid
    REAL(dp), DIMENSION(:  ),   INTENT(IN)        :: lon         ! Data lon grid
    INTEGER,                    INTENT(IN)        :: nlat        ! Number of lat elements
    INTEGER,                    INTENT(IN)        :: nlon        ! Number of lon elements

    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(OUT)       :: map_mean_i  ! On cart grid: vertex indices
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)       :: map_mean_w  ! On cart grid: weights
    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(OUT)       :: map_bilin_i ! On mesh     : cart indices (4x i and j))
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)       :: map_bilin_w ! On mesh     : weights      (4x)
    INTEGER,                             INTENT(OUT)       :: wmap_mean_i, wmap_mean_w, wmap_bilin_i, wmap_bilin_w
    
    INTEGER                                       :: vi, n
    INTEGER                                       :: il,iu,jl,ju
    REAL(dp)                                      :: wil,wiu,wjl,wju

    INTEGER,  DIMENSION(mesh%nV, 4)               :: mesh_i
    REAL(dp), DIMENSION(mesh%nV, 4)               :: mesh_w
    
    REAL(dp)                                      :: x_range_mesh, y_range_mesh
    
    ! Allocate shared memory for the mapping arrays
    CALL allocate_shared_int_2D(nlon,      nlat,     map_mean_i,   wmap_mean_i)
    CALL allocate_shared_dp_2D( nlon,      nlat,     map_mean_w,   wmap_mean_w)
    CALL allocate_shared_int_2D(mesh%nV,   4,        map_bilin_i,  wmap_bilin_i)
    CALL allocate_shared_dp_2D( mesh%nV,   4,        map_bilin_w,  wmap_bilin_w)
    
    IF (par%master) THEN
      map_mean_i  = 0
      map_mean_w  = 0._dp
      map_bilin_i = 0
      map_bilin_w = 0._dp
    END IF
    CALL sync
    
    mesh_i = 0
    mesh_w = 0._dp
    
    x_range_mesh = mesh%xmax - mesh%xmin
    y_range_mesh = mesh%ymax - mesh%ymin

    DO vi = mesh%v1, mesh%v2
      
      ! Find enveloping lat-lon indices
      il  = MAX(1,MIN(nlon-1, 1 + FLOOR((mesh%lon(vi)-MINVAL(lon)) / (lon(2)-lon(1)))))
      iu  = il+1        
      wil = (lon(iu) - mesh%lon(vi))/(lon(2)-lon(1))
      wiu = 1-wil
        
      jl  = MAX(1,MIN(nlat-1, 1 + FLOOR((mesh%lat(vi)-MINVAL(lat)) / (lat(2)-lat(1)))))
      ju  = jl+1        
      wjl = (lat(ju) - mesh%lat(vi))/(lat(2)-lat(1))
      wju = 1-wjl

      mesh_i(vi,:) = [il, iu, jl, ju]
      mesh_w(vi,:) = [wil*wjl, wil*wju, wiu*wjl, wiu*wju]

    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! Write to the actual shared array - one by one
    DO n = 0, par%n-1
      IF (n==par%i) THEN
        DO vi = 1, mesh%nV
          IF (mesh_i(vi,1)==0) CYCLE
          map_bilin_i(vi,:) = mesh_i(vi,:)
          map_bilin_w(vi,:) = mesh_w(vi,:)
        END DO
      END IF
      CALL sync
    END DO

    
  END SUBROUTINE get_glob_to_mesh_map

  ! == Subroutines for remapping data between an old mesh and a new mesh ==
  SUBROUTINE create_remapping_arrays( mesh_src, mesh_dst, map)
    ! Create remapping arrays for remapping data from mesh_src to mesh_dst, for all remapping methods
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_src
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_dst
    TYPE(type_remapping),                    INTENT(INOUT) :: map
    
    REAL(dp), DIMENSION(:), POINTER                        ::  A_src
    INTEGER                                                :: wA_src, vi, nV_max
        
    ! Create all the remapping arrays
    CALL create_remapping_arrays_trilin(            mesh_src, mesh_dst, map%trilin)
    CALL create_remapping_arrays_nearest_neighbour( mesh_src, mesh_dst, map%nearest_neighbour)
    
    ! Trick - in order to save on memory in conservative remapping, compare local resolutions of
    ! the two meshes to determine how many source mesh vertices can contribute to
    ! a destination mesh vertex.
    
    CALL allocate_shared_dp_1D( mesh_dst%nV, A_src, wA_src)    
    CALL remap_nearest_neighbour_2D( mesh_dst, map%nearest_neighbour, mesh_src%A, A_src)    
    nV_max = 0
    DO vi = 1, mesh_dst%nV
      nV_max = MAX(MAX( nV_max, CEILING( mesh_dst%A(vi) / A_src(vi))), CEILING( A_src(vi) / mesh_dst%A(vi)))
    END DO
    nV_max = nV_max * 3    
    CALL deallocate_shared( wA_src)
    
    CALL create_remapping_arrays_conservative(      mesh_src, mesh_dst, nV_max, map%conservative)
    
  END SUBROUTINE create_remapping_arrays
  SUBROUTINE deallocate_remapping_arrays( map)
    
    IMPLICIT NONE
  
    TYPE(type_remapping),                    INTENT(INOUT) :: map
    
    CALL deallocate_shared( map%trilin%wvi)
    CALL deallocate_shared( map%trilin%ww)
    
    NULLIFY( map%trilin%vi)
    NULLIFY( map%trilin%w)    
    
    CALL deallocate_shared( map%nearest_neighbour%wvi)
    
    NULLIFY( map%nearest_neighbour%vi)
    
    CALL deallocate_shared( map%conservative%wnV)
    CALL deallocate_shared( map%conservative%wvi)
    CALL deallocate_shared( map%conservative%ww0)
    CALL deallocate_shared( map%conservative%ww1x)
    CALL deallocate_shared( map%conservative%ww1y)
        
    NULLIFY( map%conservative%nV)
    NULLIFY( map%conservative%vi)
    NULLIFY( map%conservative%w0) 
    NULLIFY( map%conservative%w1x) 
    NULLIFY( map%conservative%w1y) 
    
  END SUBROUTINE deallocate_remapping_arrays
  
  SUBROUTINE create_remapping_arrays_trilin( mesh_src, mesh_dst, map)
    ! Create remapping arrays for remapping data from mesh_src to mesh_dst using trilinear interpolation
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_src
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_dst
    TYPE(type_remapping_trilin),             INTENT(INOUT) :: map
    
    ! Local variables
    INTEGER                                                :: vi_dst, ti_src, via_src, vib_src, vic_src
    REAL(dp), DIMENSION(2)                                 :: p, pa, pb, pc
    REAL(dp)                                               :: Atot, Aa, Ab, Ac
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( mesh_dst%nV, 3, map%vi, map%wvi)
    CALL allocate_shared_dp_2D(  mesh_dst%nV, 3, map%w,  map%ww)
    
    ! For all mesh_dst vertices: find the mesh_src triangle containing that vertex,
    ! and interpolate between the three vertices spanning that triangle
    
    ti_src = 1
    
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
    
      ! The vertex coordinates
      p = mesh_dst%V(vi_dst,:)
    
      ! Find the mesh_src triangle containing this vertex
      CALL find_containing_triangle( mesh_src, p, ti_src)
      
      ! Find the three vertices spanning this triangle
      via_src = mesh_src%Tri(ti_src,1)
      vib_src = mesh_src%Tri(ti_src,2)
      vic_src = mesh_src%Tri(ti_src,3)
      
      pa  = mesh_src%V(via_src,:)
      pb  = mesh_src%V(vib_src,:)
      pc  = mesh_src%V(vic_src,:)
      
      ! Calculate their interpolation weights
      Atot = find_triangle_area( pa, pb, pc)
      Aa   = find_triangle_area( pb, pc, p )
      Ab   = find_triangle_area( pc, pa, p )
      Ac   = find_triangle_area( pa, pb, p )
      
      map%vi( vi_dst,:) = [via_src, vib_src, vic_src]
      map%w(  vi_dst,:) = [Aa/Atot, Ab/Atot, Ac/Atot]
       
    END DO ! DO vi_dst = mesh_dst%v1, mesh_dst%v2
    CALL sync
    
  END SUBROUTINE create_remapping_arrays_trilin
  SUBROUTINE create_remapping_arrays_nearest_neighbour( mesh_src, mesh_dst, map)
    ! Create remapping arrays for remapping data from mesh_src to mesh_dst using nearest neighbour interpolation
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst
    TYPE(type_remapping_nearest_neighbour),  INTENT(INOUT) :: map
    
    ! Local variables
    INTEGER                                                :: vi_dst, vi_src
    REAL(dp), DIMENSION(2)                                 :: p
    
    ! Allocate shared memory
    CALL allocate_shared_int_1D( mesh_dst%nV, map%vi, map%wvi)
    
    ! For all mesh_dst vertices: find the mesh_src vertex whose Voronoi cell
    ! contains it.
    
    vi_src = 1
    
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
    
      ! The vertex coordinates
      p = mesh_dst%V(vi_dst,:)
    
      ! Find the mesh_src triangle containing this vertex
      CALL find_containing_vertex(mesh_src, p, vi_src)
      
      map%vi( vi_dst) = vi_src
      
    END DO ! DO vi_dst = mesh_dst%v1, mesh_dst%v2
    CALL sync
    
  END SUBROUTINE create_remapping_arrays_nearest_neighbour
  
  ! == Subroutines for conservative remapping
  SUBROUTINE create_remapping_arrays_conservative( mesh_src, mesh_dst, nV_max, map)
    ! Create remapping arrays for remapping data from mesh_src to mesh_dst using 1st and 2nd order conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_src ! INOUT instead of IN because the search maps and stacks are used
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_dst
    INTEGER,                                 INTENT(IN)    :: nV_max    ! ! Maximum number of src vertices contributing to a dst vertex, based on local resolution difference
    TYPE(type_remapping_conservative),       INTENT(INOUT) :: map                                           
        
    REAL(dp), DIMENSION(:,:,:), POINTER                    :: Vor_lines_src, Vor_lines_dst
    INTEGER,  DIMENSION(:,:), POINTER                      :: Vor_vi_ti_src, Vor_vi_ti_dst
    INTEGER                                                :: wVor_lines_src, wVor_lines_dst, wVor_vi_ti_src, wVor_vi_ti_dst
    
    INTEGER                                                :: vi_src_start_proc, vi_dst_start_proc
    INTEGER,  DIMENSION(:  ), ALLOCATABLE                  :: proc_domain_Ac_src
    INTEGER,  DIMENSION(:  ), ALLOCATABLE                  :: proc_domain_Ac_dst
    
    ! Results from integrating over Voronoi boundary lines for both meshes
    INTEGER,  DIMENSION(:  ), POINTER                      :: r_Ac_src_nV
    INTEGER,  DIMENSION(:,:), POINTER                      :: r_Ac_src_vi_dst_left
    INTEGER,  DIMENSION(:,:), POINTER                      :: r_Ac_src_vi_dst_right
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_Ac_src_LI_xdy
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_Ac_src_LI_mxydx
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_Ac_src_LI_xydy
    INTEGER                                                :: wr_Ac_src_nV, wr_Ac_src_vi_dst_left, wr_Ac_src_vi_dst_right, wr_Ac_src_LI_xdy, wr_Ac_src_LI_mxydx, wr_Ac_src_LI_xydy
    
    INTEGER,  DIMENSION(:  ), POINTER                      :: r_Ac_dst_nV
    INTEGER,  DIMENSION(:,:), POINTER                      :: r_Ac_dst_vi_src_left
    INTEGER,  DIMENSION(:,:), POINTER                      :: r_Ac_dst_vi_src_right
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_Ac_dst_LI_xdy
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_Ac_dst_LI_mxydx
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_Ac_dst_LI_xydy
    INTEGER                                                :: wr_Ac_dst_nV, wr_Ac_dst_vi_src_left, wr_Ac_dst_vi_src_right, wr_Ac_dst_LI_xdy, wr_Ac_dst_LI_mxydx, wr_Ac_dst_LI_xydy
    
    LOGICAL                                                :: CountCoincidences
    
    ! Integration results rearranged to actual vertices
    INTEGER,  DIMENSION(:  ), POINTER                      :: r_src_nV
    INTEGER,  DIMENSION(:,:), POINTER                      :: r_src_vi_dst
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_src_LI_xdy
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_src_LI_mxydx
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_src_LI_xydy
    INTEGER                                                :: wr_src_nV, wr_src_vi_dst, wr_src_LI_xdy, wr_src_LI_mxydx, wr_src_LI_xydy
    
    INTEGER,  DIMENSION(:  ), POINTER                      :: r_dst_nV
    INTEGER,  DIMENSION(:,:), POINTER                      :: r_dst_vi_src
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_dst_LI_xdy
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_dst_LI_mxydx
    REAL(dp), DIMENSION(:,:), POINTER                      :: r_dst_LI_xydy
    INTEGER                                                :: wr_dst_nV, wr_dst_vi_src, wr_dst_LI_xdy, wr_dst_LI_mxydx, wr_dst_LI_xydy
    
    ! ================================================================================================================================
    
    ! Allocate memory for all the temporary data used in conservative remapping
    CALL create_remapping_arrays_conservative_allocate_memory( mesh_src, mesh_dst, nV_max, &
      proc_domain_Ac_src, proc_domain_Ac_dst, &
       Vor_lines_src,  Vor_lines_dst,  Vor_vi_ti_src,  Vor_vi_ti_dst, &
      wVor_lines_src, wVor_lines_dst, wVor_vi_ti_src, wVor_vi_ti_dst, &
       r_Ac_src_nV,  r_Ac_src_vi_dst_left,  r_Ac_src_vi_dst_right,  r_Ac_src_LI_xdy,  r_Ac_src_LI_mxydx,  r_Ac_src_LI_xydy, &
      wr_Ac_src_nV, wr_Ac_src_vi_dst_left, wr_Ac_src_vi_dst_right, wr_Ac_src_LI_xdy, wr_Ac_src_LI_mxydx, wr_Ac_src_LI_xydy, &
       r_Ac_dst_nV,  r_Ac_dst_vi_src_left,  r_Ac_dst_vi_src_right,  r_Ac_dst_LI_xdy,  r_Ac_dst_LI_mxydx,  r_Ac_dst_LI_xydy, &
      wr_Ac_dst_nV, wr_Ac_dst_vi_src_left, wr_Ac_dst_vi_src_right, wr_Ac_dst_LI_xdy, wr_Ac_dst_LI_mxydx, wr_Ac_dst_LI_xydy, &
       r_src_nV,  r_src_vi_dst,  r_src_LI_xdy,  r_src_LI_mxydx,  r_src_LI_xydy, &
      wr_src_nV, wr_src_vi_dst, wr_src_LI_xdy, wr_src_LI_mxydx, wr_src_LI_xydy, &
       r_dst_nV,  r_dst_vi_src,  r_dst_LI_xdy,  r_dst_LI_mxydx,  r_dst_LI_xydy, &
      wr_dst_nV, wr_dst_vi_src, wr_dst_LI_xdy, wr_dst_LI_mxydx, wr_dst_LI_xydy, map)
    
    ! Find the coordinates and relevant indices of the Voronoi boundary lines.
    ! Since each line describes the boundary between the Voronoi cells of two
    ! connected vertices, each Voronoi line can be uniquely described by an Aci index.
    
    CALL find_Voronoi_boundary_lines( mesh_src, Vor_lines_src, Vor_vi_ti_src)
    CALL find_Voronoi_boundary_lines( mesh_dst, Vor_lines_dst, Vor_vi_ti_dst)
    
    ! Determine process vertex domains for both meshes
  
    CALL determine_process_domains( mesh_src, proc_domain_Ac_src, vi_src_start_proc)
    CALL determine_process_domains( mesh_dst, proc_domain_Ac_dst, vi_dst_start_proc)

    ! Integrate over all Voronoi cell boundary lines for both meshes
  
    CountCoincidences = .FALSE.  
    CALL integrate_over_Voronoi_boundaries( mesh_src, mesh_dst, proc_domain_Ac_src, &
      vi_src_start_proc, nV_max, Vor_lines_src, Vor_lines_dst, Vor_vi_ti_dst, &
      r_Ac_src_nV, r_Ac_src_vi_dst_left, r_Ac_src_vi_dst_right, r_Ac_src_LI_xdy, r_Ac_src_LI_mxydx, r_Ac_src_LI_xydy, CountCoincidences)
  
    CountCoincidences = .TRUE.  
    CALL integrate_over_Voronoi_boundaries( mesh_dst, mesh_src, proc_domain_Ac_dst, &
      vi_dst_start_proc, nV_max, Vor_lines_dst, Vor_lines_src, Vor_vi_ti_src, &
      r_Ac_dst_nV, r_Ac_dst_vi_src_left, r_Ac_dst_vi_src_right, r_Ac_dst_LI_xdy, r_Ac_dst_LI_mxydx, r_Ac_dst_LI_xydy, CountCoincidences)

    ! Rearrange integral contributions from Aci to vertices
  
    CALL rearrange_contributions_from_lines_to_vertices( mesh_src, r_Ac_src_nV, r_Ac_src_vi_dst_left, r_Ac_src_vi_dst_right, &
      r_Ac_src_LI_xdy, r_Ac_src_LI_mxydx, r_Ac_src_LI_xydy, r_src_nV, r_src_vi_dst, r_src_LI_xdy, r_src_LI_mxydx, r_src_LI_xydy)
  
    CALL rearrange_contributions_from_lines_to_vertices( mesh_dst, r_Ac_dst_nV, r_Ac_dst_vi_src_left, r_Ac_dst_vi_src_right, &
      r_Ac_dst_LI_xdy, r_Ac_dst_LI_mxydx, r_Ac_dst_LI_xydy, r_dst_nV, r_dst_vi_src, r_dst_LI_xdy, r_dst_LI_mxydx, r_dst_LI_xydy)

    ! Integrate around domain boundary
    
    IF (par%master) CALL integrate_around_domain_boundary( mesh_dst, mesh_src, Vor_lines_dst, Vor_lines_src, &
      r_dst_nV, r_dst_vi_src, r_dst_LI_xdy, r_dst_LI_mxydx, r_dst_LI_xydy)
    CALL sync

    ! Add contributions from mesh_src to mesh_dst
    
    CALL add_contributions_from_opposite_mesh( mesh_dst, r_dst_nV, r_dst_vi_src, r_dst_LI_xdy, r_dst_LI_mxydx, r_dst_LI_xydy, &
      r_src_nV, r_src_vi_dst, r_src_LI_xdy, r_src_LI_mxydx, r_src_LI_xydy)

    ! Finish incomplete mesh_dst vertices
  
    CALL finish_incomplete_vertices( mesh_dst, mesh_src, nV_max, Vor_lines_dst, r_dst_nV, r_dst_vi_src, r_dst_LI_xdy, r_dst_LI_mxydx, r_dst_LI_xydy, &
      r_src_nV, r_src_vi_dst, r_src_LI_xdy, r_src_LI_mxydx, r_src_LI_xydy)

    ! Convert line integrals to remapping weights
    
    CALL calculate_remapping_weights_from_line_integrals( mesh_dst, mesh_src, r_dst_nV, r_dst_vi_src, r_dst_LI_xdy, r_dst_LI_mxydx, r_dst_LI_xydy, map)
      
    ! Check if everything worked
    
    CALL check_if_remapping_is_conservative( mesh_src, mesh_dst, map)
    
    ! Clean up after yourself 
       
    CALL create_remapping_arrays_conservative_deallocate_memory( &
      proc_domain_Ac_src, proc_domain_Ac_dst, &
      wVor_lines_src, wVor_lines_dst, wVor_vi_ti_src, wVor_vi_ti_dst, &
      wr_Ac_src_nV, wr_Ac_src_vi_dst_left, wr_Ac_src_vi_dst_right, wr_Ac_src_LI_xdy, wr_Ac_src_LI_mxydx, wr_Ac_src_LI_xydy, &
      wr_Ac_dst_nV, wr_Ac_dst_vi_src_left, wr_Ac_dst_vi_src_right, wr_Ac_dst_LI_xdy, wr_Ac_dst_LI_mxydx, wr_Ac_dst_LI_xydy, &
      wr_src_nV, wr_src_vi_dst, wr_src_LI_xdy, wr_src_LI_mxydx, wr_src_LI_xydy, &
      wr_dst_nV, wr_dst_vi_src, wr_dst_LI_xdy, wr_dst_LI_mxydx, wr_dst_LI_xydy)
    
  END SUBROUTINE create_remapping_arrays_conservative  
  SUBROUTINE create_remapping_arrays_conservative_allocate_memory( mesh_src, mesh_dst, nV_max, &
    proc_domain_Ac_src, proc_domain_Ac_dst, &
     Vor_lines_src,  Vor_lines_dst,  Vor_vi_ti_src,  Vor_vi_ti_dst, &
    wVor_lines_src, wVor_lines_dst, wVor_vi_ti_src, wVor_vi_ti_dst, &
     r_Ac_src_nV,  r_Ac_src_vi_dst_left,  r_Ac_src_vi_dst_right,  r_Ac_src_LI_xdy,  r_Ac_src_LI_mxydx,  r_Ac_src_LI_xydy, &
    wr_Ac_src_nV, wr_Ac_src_vi_dst_left, wr_Ac_src_vi_dst_right, wr_Ac_src_LI_xdy, wr_Ac_src_LI_mxydx, wr_Ac_src_LI_xydy, &
     r_Ac_dst_nV,  r_Ac_dst_vi_src_left,  r_Ac_dst_vi_src_right,  r_Ac_dst_LI_xdy,  r_Ac_dst_LI_mxydx,  r_Ac_dst_LI_xydy, &
    wr_Ac_dst_nV, wr_Ac_dst_vi_src_left, wr_Ac_dst_vi_src_right, wr_Ac_dst_LI_xdy, wr_Ac_dst_LI_mxydx, wr_Ac_dst_LI_xydy, &
     r_src_nV,  r_src_vi_dst,  r_src_LI_xdy,  r_src_LI_mxydx,  r_src_LI_xydy, &
    wr_src_nV, wr_src_vi_dst, wr_src_LI_xdy, wr_src_LI_mxydx, wr_src_LI_xydy, &
     r_dst_nV,  r_dst_vi_src,  r_dst_LI_xdy,  r_dst_LI_mxydx,  r_dst_LI_xydy, &
    wr_dst_nV, wr_dst_vi_src, wr_dst_LI_xdy, wr_dst_LI_mxydx, wr_dst_LI_xydy, map)
    
    ! Allocate memory for all the temporary arrays used in conservative remapping.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst
    INTEGER,                                 INTENT(IN)    :: nV_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE, INTENT(INOUT) :: proc_domain_Ac_src, proc_domain_Ac_dst
    REAL(dp), DIMENSION(:,:,:), POINTER,     INTENT(INOUT) ::  Vor_lines_src,  Vor_lines_dst
    INTEGER,                                 INTENT(INOUT) :: wVor_lines_src, wVor_lines_dst
    INTEGER,  DIMENSION(:,:  ), POINTER,     INTENT(INOUT) ::  Vor_vi_ti_src,  Vor_vi_ti_dst
    INTEGER,                                 INTENT(INOUT) :: wVor_vi_ti_src, wVor_vi_ti_dst
    INTEGER,  DIMENSION(:    ), POINTER,     INTENT(INOUT) ::  r_Ac_src_nV,  r_Ac_dst_nV
    INTEGER,                                 INTENT(INOUT) :: wr_Ac_src_nV, wr_Ac_dst_nV
    INTEGER,  DIMENSION(:,:  ), POINTER,     INTENT(INOUT) ::  r_Ac_src_vi_dst_left,  r_Ac_src_vi_dst_right,  r_Ac_dst_vi_src_left,  r_Ac_dst_vi_src_right
    INTEGER,                                 INTENT(INOUT) :: wr_Ac_src_vi_dst_left, wr_Ac_src_vi_dst_right, wr_Ac_dst_vi_src_left, wr_Ac_dst_vi_src_right
    REAL(dp), DIMENSION(:,:  ), POINTER,     INTENT(INOUT) ::  r_Ac_src_LI_xdy,  r_Ac_src_LI_mxydx,  r_Ac_src_LI_xydy,  r_Ac_dst_LI_xdy,  r_Ac_dst_LI_mxydx,  r_Ac_dst_LI_xydy
    INTEGER,                                 INTENT(INOUT) :: wr_Ac_src_LI_xdy, wr_Ac_src_LI_mxydx, wr_Ac_src_LI_xydy, wr_Ac_dst_LI_xdy, wr_Ac_dst_LI_mxydx, wr_Ac_dst_LI_xydy
    INTEGER,  DIMENSION(:    ), POINTER,     INTENT(INOUT) ::  r_src_nV,  r_dst_nV
    INTEGER,                                 INTENT(INOUT) :: wr_src_nV, wr_dst_nV
    INTEGER,  DIMENSION(:,:  ), POINTER,     INTENT(INOUT) ::  r_src_vi_dst,  r_dst_vi_src
    INTEGER,                                 INTENT(INOUT) :: wr_src_vi_dst, wr_dst_vi_src
    REAL(dp), DIMENSION(:,:  ), POINTER,     INTENT(INOUT) ::  r_src_LI_xdy,  r_src_LI_mxydx,  r_src_LI_xydy,  r_dst_LI_xdy,  r_dst_LI_mxydx,  r_dst_LI_xydy
    INTEGER,                                 INTENT(INOUT) :: wr_src_LI_xdy, wr_src_LI_mxydx, wr_src_LI_xydy, wr_dst_LI_xdy, wr_dst_LI_mxydx, wr_dst_LI_xydy
    TYPE(type_remapping_conservative),       INTENT(INOUT) :: map
    
    ! Process domain maps
    ALLOCATE( proc_domain_Ac_src( mesh_src%nAc))
    ALLOCATE( proc_domain_Ac_dst( mesh_dst%nAc))
    
    ! Coordinates and different indices of the Voronoi boundary lines (one for each Ac vertex)
    CALL allocate_shared_dp_3D(  mesh_src%nAc, 2, 2, Vor_lines_src, wVor_lines_src)
    CALL allocate_shared_dp_3D(  mesh_dst%nAc, 2, 2, Vor_lines_dst, wVor_lines_dst)
    CALL allocate_shared_int_2D( mesh_src%nAc, 6,    Vor_vi_ti_src, wVor_vi_ti_src)
    CALL allocate_shared_int_2D( mesh_dst%nAc, 6,    Vor_vi_ti_dst, wVor_vi_ti_dst)
    
    ! Results from integrating over all Voronoi cell boundary lines - Aci data
    
    CALL allocate_shared_int_1D( mesh_src%nAc,          r_Ac_src_nV,            wr_Ac_src_nV           )
    CALL allocate_shared_int_2D( mesh_src%nAc, nV_max,  r_Ac_src_vi_dst_left,   wr_Ac_src_vi_dst_left  )
    CALL allocate_shared_int_2D( mesh_src%nAc, nV_max,  r_Ac_src_vi_dst_right,  wr_Ac_src_vi_dst_right )
    CALL allocate_shared_dp_2D(  mesh_src%nAc, nV_max,  r_Ac_src_LI_xdy,        wr_Ac_src_LI_xdy       )
    CALL allocate_shared_dp_2D(  mesh_src%nAc, nV_max,  r_Ac_src_LI_mxydx,      wr_Ac_src_LI_mxydx     )
    CALL allocate_shared_dp_2D(  mesh_src%nAc, nV_max,  r_Ac_src_LI_xydy,       wr_Ac_src_LI_xydy      )
    
    CALL allocate_shared_int_1D( mesh_dst%nAc,          r_Ac_dst_nV,            wr_Ac_dst_nV           )
    CALL allocate_shared_int_2D( mesh_dst%nAc, nV_max,  r_Ac_dst_vi_src_left,   wr_Ac_dst_vi_src_left  )
    CALL allocate_shared_int_2D( mesh_dst%nAc, nV_max,  r_Ac_dst_vi_src_right,  wr_Ac_dst_vi_src_right )
    CALL allocate_shared_dp_2D(  mesh_dst%nAc, nV_max,  r_Ac_dst_LI_xdy,        wr_Ac_dst_LI_xdy       )
    CALL allocate_shared_dp_2D(  mesh_dst%nAc, nV_max,  r_Ac_dst_LI_mxydx,      wr_Ac_dst_LI_mxydx     )
    CALL allocate_shared_dp_2D(  mesh_dst%nAc, nV_max,  r_Ac_dst_LI_xydy,       wr_Ac_dst_LI_xydy      )

    ! Results from integrating over all Voronoi cell boundary lines - vertex data
    
    CALL allocate_shared_int_1D( mesh_src%nV,           r_src_nV,               wr_src_nV              )
    CALL allocate_shared_int_2D( mesh_src%nV,  nV_max,  r_src_vi_dst,           wr_src_vi_dst          )
    CALL allocate_shared_dp_2D(  mesh_src%nV,  nV_max,  r_src_LI_xdy,           wr_src_LI_xdy          )
    CALL allocate_shared_dp_2D(  mesh_src%nV,  nV_max,  r_src_LI_mxydx,         wr_src_LI_mxydx        )
    CALL allocate_shared_dp_2D(  mesh_src%nV,  nV_max,  r_src_LI_xydy,          wr_src_LI_xydy         )
     
    CALL allocate_shared_int_1D( mesh_dst%nV,           r_dst_nV,               wr_dst_nV              )
    CALL allocate_shared_int_2D( mesh_dst%nV,  nV_max,  r_dst_vi_src,           wr_dst_vi_src          )
    CALL allocate_shared_dp_2D(  mesh_dst%nV,  nV_max,  r_dst_LI_xdy,           wr_dst_LI_xdy          )
    CALL allocate_shared_dp_2D(  mesh_dst%nV,  nV_max,  r_dst_LI_mxydx,         wr_dst_LI_mxydx        )
    CALL allocate_shared_dp_2D(  mesh_dst%nV,  nV_max,  r_dst_LI_xydy,          wr_dst_LI_xydy         )
    
    ! Actual mapping weights
     
    CALL allocate_shared_int_1D( mesh_dst%nV,           map%nV,                 map%wnV                )
    CALL allocate_shared_int_2D( mesh_dst%nV,  nV_max,  map%vi,                 map%wvi                )
    CALL allocate_shared_dp_2D(  mesh_dst%nV,  nV_max,  map%w0,                 map%ww0                )
    CALL allocate_shared_dp_2D(  mesh_dst%nV,  nV_max,  map%w1x,                map%ww1x               )
    CALL allocate_shared_dp_2D(  mesh_dst%nV,  nV_max,  map%w1y,                map%ww1y               )
    
  END SUBROUTINE create_remapping_arrays_conservative_allocate_memory
  SUBROUTINE find_Voronoi_boundary_lines( mesh, Vor_lines, Vor_vi_ti)
    ! Create list of Voronoi cell boundary lines, stored in Ac data
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:,:),              INTENT(OUT)   :: Vor_lines
    INTEGER,  DIMENSION(:,:),                INTENT(OUT)   :: Vor_vi_ti
    
    ! Local variables
    INTEGER                                                :: aci
    INTEGER                                                :: via, vib, vil, vir
    INTEGER                                                :: til, tir
    INTEGER                                                :: iti, iti1, iti2, ti, ti1, ti2
    INTEGER                                                :: n1, n2, n3
    REAL(dp), DIMENSION(2)                                 :: ccl, ccr, p1, p2
    REAL(dp)                                               :: l1a, l1b, l1c, l2a, l2b, l2c

    DO aci = mesh%ac1, mesh%ac2
    
      via = mesh%Aci( aci, 1)
      vib = mesh%Aci( aci, 2)
      vil = mesh%Aci( aci, 3)
      vir = mesh%Aci( aci, 4)

      IF (.NOT. is_boundary_segment( mesh, via, vib)) THEN
        ! Both of the two vertices are free vertices, both the left and right
        ! vertex and triangle should exist

        ! Find the left and right triangles
        til = 0
        tir = 0
        DO iti1 = 1, mesh%niTri( via)
          iti2 = iti1 + 1
          IF (iti2 == mesh%niTri( via)+1) iti2 = 1
          ti1 = mesh%iTri( via, iti1)
          ti2 = mesh%iTri( via, iti2)
          DO n1 = 1, 3
            n2 = n1 + 1
            IF (n2==4) n2 = 1
            n3 = n2 + 1
            IF (n3==4) n3 = 1
            IF (mesh%Tri( ti1, n1) == via .AND. mesh%Tri( ti1, n2) == vir .AND. mesh%Tri( ti1, n3) == vib) THEN
              tir = ti1
              til = ti2
            END IF
          END DO
        END DO

        Vor_vi_ti( aci, :) = [via, vib, vir, vil, tir, til]

        ! Find the circumcenters of these two triangles. Crop them IF necessary.
        ccl = mesh%Tricc( til, :)
        ccr = mesh%Tricc( tir, :)

        IF     (ccl(1) < mesh%xmin) THEN
          p1 = [mesh%xmin, mesh%ymin]
          p2 = [mesh%xmin, mesh%ymax]
          CALL line_from_points( ccr, ccl, l1a, l1b, l1c)
          CALL line_from_points( p1,  p2 , l2a, l2b, l2c)
          CALL line_line_intersection( l1a, l1b, l1c, l2a, l2b, l2c, ccl)
        ELSEIF (ccl(1) > mesh%xmax) THEN
          p1 = [mesh%xmax, mesh%ymin]
          p2 = [mesh%xmax, mesh%ymax]
          CALL line_from_points( ccr, ccl, l1a, l1b, l1c)
          CALL line_from_points( p1,  p2 , l2a, l2b, l2c)
          CALL line_line_intersection( l1a, l1b, l1c, l2a, l2b, l2c, ccl)
        ELSEIF (ccl(2) < mesh%ymin) THEN
          p1 = [mesh%xmin, mesh%ymin]
          p2 = [mesh%xmax, mesh%ymin]
          CALL line_from_points( ccr, ccl, l1a, l1b, l1c)
          CALL line_from_points( p1,  p2 , l2a, l2b, l2c)
          CALL line_line_intersection( l1a, l1b, l1c, l2a, l2b, l2c, ccl)
        ELSEIF (ccl(2) > mesh%ymax) THEN
          p1 = [mesh%xmin, mesh%ymax]
          p2 = [mesh%xmax, mesh%ymax]
          CALL line_from_points( ccr, ccl, l1a, l1b, l1c)
          CALL line_from_points( p1,  p2 , l2a, l2b, l2c)
          CALL line_line_intersection( l1a, l1b, l1c, l2a, l2b, l2c, ccl)
        END IF

        IF     (ccr(1) < mesh%xmin) THEN
          p1 = [mesh%xmin, mesh%ymin]
          p2 = [mesh%xmin, mesh%ymax]
          CALL line_from_points( ccr, ccl, l1a, l1b, l1c)
          CALL line_from_points( p1,  p2 , l2a, l2b, l2c)
          CALL line_line_intersection( l1a, l1b, l1c, l2a, l2b, l2c, ccr)
        ELSEIF (ccr(1) > mesh%xmax) THEN
          p1 = [mesh%xmax, mesh%ymin]
          p2 = [mesh%xmax, mesh%ymax]
          CALL line_from_points( ccr, ccl, l1a, l1b, l1c)
          CALL line_from_points( p1,  p2 , l2a, l2b, l2c)
          CALL line_line_intersection( l1a, l1b, l1c, l2a, l2b, l2c, ccr)
        ELSEIF (ccr(2) < mesh%ymin) THEN
          p1 = [mesh%xmin, mesh%ymin]
          p2 = [mesh%xmax, mesh%ymin]
          CALL line_from_points( ccr, ccl, l1a, l1b, l1c)
          CALL line_from_points( p1,  p2 , l2a, l2b, l2c)
          CALL line_line_intersection( l1a, l1b, l1c, l2a, l2b, l2c, ccr)
        ELSEIF (ccr(2) > mesh%ymax) THEN
          p1 = [mesh%xmin, mesh%ymax]
          p2 = [mesh%xmax, mesh%ymax]
          CALL line_from_points( ccr, ccl, l1a, l1b, l1c)
          CALL line_from_points( p1,  p2 , l2a, l2b, l2c)
          CALL line_line_intersection( l1a, l1b, l1c, l2a, l2b, l2c, ccr)
        END IF

        Vor_lines( aci, 1, 1) = ccr(1)
        Vor_lines( aci, 1, 2) = ccr(2)
        Vor_lines( aci, 2, 1) = ccl(1)
        Vor_lines( aci, 2, 2) = ccl(2)

      ELSE ! IF (.NOT. is_boundary_segment( mesh, via, vib)) THEN
        ! Both vertices are boundary vertices, so only one triangle exists.

        ! Find the triangle (can be either left or right, dIFference matters!)
        til = 0
        tir = 0
        vil = 0
        vir = 0
        DO iti = 1, mesh%niTri( via)
          ti = mesh%iTri( via, iti)
          DO n1 = 1, 3
            n2 = n1 + 1
            IF (n2==4) n2 = 1
            n3 = n2 + 1
            IF (n3==4) n3 = 1
            IF     (mesh%Tri( ti, n1) == via .AND. mesh%Tri( ti, n2) == vib) THEN
              til = ti
              vil = mesh%Tri( ti, n3)
            ELSEIF (mesh%Tri( ti, n1) == vib .AND. mesh%Tri( ti, n2) == via) THEN
              tir = ti
              vir = mesh%Tri( ti, n3)
            END IF
          END DO
        END DO

        Vor_vi_ti( aci, :) = [via, vib, vir, vil, tir, til]

        ! Find the endpoints of this Voronoi boundary line
        IF (til > 0) THEN
          ! Triangle lies to the left of the Aci line, Voronoi line runs from
          ! domain boundary to circumcenter

          p2 = mesh%Tricc(til,:)

          IF     ((mesh%edge_index(via)==8 .OR. mesh%edge_index(via)==1 .OR. mesh%edge_index(via)==2) .AND. &
                  (mesh%edge_index(vib)==8 .OR. mesh%edge_index(vib)==1 .OR. mesh%edge_index(vib)==2)) THEN
            ! North
            p1 = [p2(1), mesh%ymax]
          ELSEIF ((mesh%edge_index(via)==2 .OR. mesh%edge_index(via)==3 .OR. mesh%edge_index(via)==4) .AND. &
                  (mesh%edge_index(vib)==2 .OR. mesh%edge_index(vib)==3 .OR. mesh%edge_index(vib)==4)) THEN
            ! East
            p1 = [mesh%xmax, p2(2)]
          ELSEIF ((mesh%edge_index(via)==4 .OR. mesh%edge_index(via)==5 .OR. mesh%edge_index(via)==6) .AND. &
                  (mesh%edge_index(vib)==4 .OR. mesh%edge_index(vib)==5 .OR. mesh%edge_index(vib)==6)) THEN
            ! South
            p1 = [p2(1), mesh%ymin]
          ELSEIF ((mesh%edge_index(via)==6 .OR. mesh%edge_index(via)==7 .OR. mesh%edge_index(via)==8) .AND. &
                  (mesh%edge_index(vib)==6 .OR. mesh%edge_index(vib)==7 .OR. mesh%edge_index(vib)==8)) THEN
            ! West
            p1 = [mesh%xmin, p2(2)]
          END IF

        ELSEIF (tir > 0) THEN
          ! Triangle lies to the right of the Aci line, Voronoi line runs from
          ! circumcenter tot domain boundary

          p1 = mesh%Tricc(tir,:)

          IF     ((mesh%edge_index(via)==8 .OR. mesh%edge_index(via)==1 .OR. mesh%edge_index(via)==2) .AND. &
                  (mesh%edge_index(vib)==8 .OR. mesh%edge_index(vib)==1 .OR. mesh%edge_index(vib)==2)) THEN
            ! North
            p2 = [p1(1), mesh%ymax]
          ELSEIF ((mesh%edge_index(via)==2 .OR. mesh%edge_index(via)==3 .OR. mesh%edge_index(via)==4) .AND. &
                  (mesh%edge_index(vib)==2 .OR. mesh%edge_index(vib)==3 .OR. mesh%edge_index(vib)==4)) THEN
            ! East
            p2 = [mesh%xmax, p1(2)]
          ELSEIF ((mesh%edge_index(via)==4 .OR. mesh%edge_index(via)==5 .OR. mesh%edge_index(via)==6) .AND. &
                  (mesh%edge_index(vib)==4 .OR. mesh%edge_index(vib)==5 .OR. mesh%edge_index(vib)==6)) THEN
            ! South
            p2 = [p1(1), mesh%ymin]
          ELSEIF ((mesh%edge_index(via)==6 .OR. mesh%edge_index(via)==7 .OR. mesh%edge_index(via)==8) .AND. &
                  (mesh%edge_index(vib)==6 .OR. mesh%edge_index(vib)==7 .OR. mesh%edge_index(vib)==8)) THEN
            ! West
            p2 = [mesh%xmin, p1(2)]
          END IF

        END IF
    
        ! The circumcenter of an edge triangle might lie outside the domain.
        p1 = [MAX( MIN( p1(1), mesh%xmax), mesh%xmin), MAX( MIN( p1(2), mesh%ymax), mesh%ymin)]
        p2 = [MAX( MIN( p2(1), mesh%xmax), mesh%xmin), MAX( MIN( p2(2), mesh%ymax), mesh%ymin)]

        Vor_lines( aci, 1, 1) = p1(1)
        Vor_lines( aci, 1, 2) = p1(2)
        Vor_lines( aci, 2, 1) = p2(1)
        Vor_lines( aci, 2, 2) = p2(2)

      END IF ! IF (.NOT. is_boundary_segment( mesh, via, vib)) THEN
      
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync

  END SUBROUTINE find_Voronoi_boundary_lines
  SUBROUTINE determine_process_domains( mesh, proc_domain_Ac, vi_start_proc)
    ! Determine which vertices and Voronoi boundary lines (Ac vertices) may be treated by which process,
    ! for optimal load-balanced parallelisation
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),              INTENT(OUT)   :: proc_domain_Ac
    INTEGER,                                 INTENT(OUT)   :: vi_start_proc
    
    ! Local variables:
    INTEGER                                                :: i, vi, ci, aci
    REAL(dp)                                               :: xmin, xmax
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                :: proc_domain_V
    
    ALLOCATE( proc_domain_V( mesh%nV))
    
    proc_domain_V  = 0
    proc_domain_Ac = 0
    vi_start_proc  = 0
    
    ! Distribute vertices over processes using load-balanced domain sizes
    DO i = 0, par%n-1
  
      CALL partition_domain_x_balanced( mesh, mesh%ymin, mesh%ymax, i, par%n, xmin, xmax)
      
      ! Make domains very slightly wider to prevent vertices that lie exactly on the boundary from being excluded
      xmin = xmin - ((mesh%xmax - mesh%xmin) * 0.001_dp)
      xmax = xmax + ((mesh%xmax - mesh%xmin) * 0.001_dp)
    
      DO vi = 1, mesh%nV
        IF (mesh%V(vi,1)>=xmin .AND. mesh%V(vi,1) <= xmax) THEN
        
          ! This vertex lies in this process' domain; mark it
          proc_domain_V(vi) = i
    
          ! Mark all its Ac connections
          DO ci = 1, mesh%nC( vi)
            aci = mesh%iAci( vi, ci)
            proc_domain_Ac( aci) = i
          END DO
       
        END IF
      END DO ! DO vi = 1, mesh%nV
      
    END DO ! DO i = 0, par%n-1
    
    ! Find starting vertex
    vi = 0
    DO WHILE (vi_start_proc == 0)
      vi = vi + 1
      IF (proc_domain_V(vi)==par%i) vi_start_proc = vi
    END DO
    
    DEALLOCATE( proc_domain_V)
    
    CALL sync
    
  END SUBROUTINE determine_process_domains
  SUBROUTINE integrate_over_Voronoi_boundaries( mesh_top, mesh_bot,  proc_domain_Ac_top, &
    vi_top_start_proc, nV_max, Vor_lines_top, Vor_lines_bot, Vor_vi_ti_bot, &
    r_Ac_nV, r_Ac_vi_bot_left, r_Ac_vi_bot_right, r_Ac_LI_xdy, r_Ac_LI_mxydx, r_Ac_LI_xydy, CountCoincidences)
    ! Integrate over all the Voronoi boundarie lines of mesh_top, following them through mesh_bot.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_top
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_bot
    INTEGER,  DIMENSION(:    ),              INTENT(IN)    :: proc_domain_Ac_top
    INTEGER,                                 INTENT(IN)    :: vi_top_start_proc
    INTEGER,                                 INTENT(IN)    :: nV_max
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_top, Vor_lines_bot
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: Vor_vi_ti_bot
    INTEGER,  DIMENSION(:    ),              INTENT(INOUT) :: r_Ac_nV
    INTEGER,  DIMENSION(:,:  ),              INTENT(INOUT) :: r_Ac_vi_bot_left, r_Ac_vi_bot_right
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: r_Ac_LI_xdy, r_Ac_LI_mxydx, r_Ac_LI_xydy
    LOGICAL,                                 INTENT(IN)    :: CountCoincidences
    
    ! Local variables:
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                :: Aci_map_top
    INTEGER                                                :: vi_top, vi_bot, ci, aci, vc, vii
    INTEGER                                                :: r_sng_nV                 ! How many Voronoi cells of the opposite mesh does this line pass through
    INTEGER,  DIMENSION(:  ), ALLOCATABLE                  :: r_sng_vi_left            ! Which vertices of the opposite mesh lie to the left  of this line
    INTEGER,  DIMENSION(:  ), ALLOCATABLE                  :: r_sng_vi_right           ! Which vertices of the opposite mesh lie to the right of this line
    REAL(dp), DIMENSION(:  ), ALLOCATABLE                  :: r_sng_LI_xdy             ! The three line integrals
    REAL(dp), DIMENSION(:  ), ALLOCATABLE                  :: r_sng_LI_mxydx
    REAL(dp), DIMENSION(:  ), ALLOCATABLE                  :: r_sng_LI_xydy

    ! Results from integrating over a single Voronoi cell boundary line
    ! (process-local memory, since this step is parallelised)    
    ALLOCATE( r_sng_vi_left(  nV_max))
    ALLOCATE( r_sng_vi_right( nV_max))
    ALLOCATE( r_sng_LI_xdy(   nV_max))
    ALLOCATE( r_sng_LI_mxydx( nV_max))
    ALLOCATE( r_sng_LI_xydy(  nV_max))

    ! Move flood-fill style through the mesh vertices, to make searching
    ! for containing opposite mesh vertices more efficient.
    
    mesh_top%VMap     = 0
    mesh_top%VStack1  = 0
    mesh_top%VStackN1 = 0

    mesh_top%VStack1(1) = vi_top_start_proc
    mesh_top%VStackN1   = 1
  
    ! Keep track of which Voronoi boundary lines (Ac vertices) we've already treated
    ALLOCATE( Aci_map_top( mesh_top%nAc))
    Aci_map_top = 0

    vi_bot            = 1

    DO WHILE (mesh_top%VStackN1 > 0)

      ! Flood-fill
      ! Take the last vertex from the stack, mark it as checked
      vi_top = mesh_top%VStack1( mesh_top%VStackN1)
      mesh_top%VStackN1 = mesh_top%VStackN1 - 1
      mesh_top%VMap( vi_top) = 1
      ! Add all its unchecked neighbours to the stack
      DO ci = 1, mesh_top%nC( vi_top)
        vc = mesh_top%C( vi_top, ci)
        IF (mesh_top%VMap( vc) == 1) CYCLE ! This neighbour has already been passed
   !     IF (proc_domain_V_top( vi_top) /= par%i) CYCLE ! This neighbour is outside of this process domain
        mesh_top%VStackN1 = mesh_top%VStackN1 + 1
        mesh_top%VStack1( mesh_top%VStackN1) = vc
        mesh_top%VMap( vc) = 1
      END DO

      ! Integrate over all aci vertices that haven't been treated yet
      DO ci = 1, mesh_top%nC( vi_top)
        
        aci = mesh_top%iAci( vi_top, ci)
        
        ! Skip lines that we've already treated
        IF (Aci_map_top( aci) == 1) CYCLE
        Aci_map_top( aci) = 1
        
        ! Skip lines that are outside the process domain
        IF (proc_domain_Ac_top( aci) /= par%i) CYCLE

        ! Integrate
        CALL calculate_line_integral_contributions( mesh_bot, Vor_lines_top, Vor_lines_bot, Vor_vi_ti_bot, aci, &
          r_sng_nV, r_sng_vi_left, r_sng_vi_right, r_sng_LI_xdy, r_sng_LI_mxydx, r_sng_LI_xydy, vi_bot, CountCoincidences, nV_max)

        ! Store data
        r_Ac_nV( aci) = r_sng_nV
        DO vii = 1, r_sng_nV
          r_Ac_vi_bot_left(  aci, vii) = r_sng_vi_left(  vii)
          r_Ac_vi_bot_right( aci, vii) = r_sng_vi_right( vii)
          r_Ac_LI_xdy(       aci, vii) = r_sng_LI_xdy(   vii)
          r_Ac_LI_mxydx(     aci, vii) = r_sng_LI_mxydx( vii)
          r_Ac_LI_xydy(      aci, vii) = r_sng_LI_xydy(  vii)
        END DO
      END DO ! DO ci = 1, mesh_top%nC( vi_top)
      
    END DO ! DO WHILE (mesh_top%VStackN > 0)
    CALL sync   
    
    DEALLOCATE( r_sng_vi_left)
    DEALLOCATE( r_sng_vi_right)
    DEALLOCATE( r_sng_LI_xdy)
    DEALLOCATE( r_sng_LI_mxydx)
    DEALLOCATE( r_sng_LI_xydy) 
    DEALLOCATE( Aci_map_top)   
    
  END SUBROUTINE integrate_over_Voronoi_boundaries
  SUBROUTINE rearrange_contributions_from_lines_to_vertices( mesh, r_Ac_nV, r_Ac_vi_bot_left, r_Ac_vi_bot_right, &
    r_Ac_LI_xdy, r_Ac_LI_mxydx, r_Ac_LI_xydy, r_nV, r_vi_bot, r_LI_xdy, r_LI_mxydx, r_LI_xydy)
    ! Rearrange the line integral contributions from Voronoi boundary lines (Ac vertices) to vertices.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),              INTENT(IN)    :: r_Ac_nV
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: r_Ac_vi_bot_left
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: r_Ac_vi_bot_right
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_Ac_LI_xdy
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_Ac_LI_mxydx
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_Ac_LI_xydy
    INTEGER,  DIMENSION(:    ),              INTENT(INOUT) :: r_nV
    INTEGER,  DIMENSION(:,:  ),              INTENT(INOUT) :: r_vi_bot
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: r_LI_xdy
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: r_LI_mxydx
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: r_LI_xydy
    
    ! Local variables:
    INTEGER                                                :: vi, ci, aci, dir, vii, vii2, vii3, vi_bot
    REAL(dp)                                               :: LI_xdy, LI_mxydx, LI_xydy
    LOGICAl                                                :: IsListed

    DO vi = mesh%v1, mesh%v2

      ! For each vertex, go over all its Ac connections, and add their
      ! contributions to this vertex

      DO ci = 1, mesh%nC( vi)
        aci = mesh%iAci( vi, ci)

        dir = 0
        IF     (mesh%Aci( aci,1) == vi) THEN
          dir = 0
        ELSEIF (mesh%Aci( aci,2) == vi) THEN
          dir = 1
        END IF

        ! Go over all opposite mesh vertices contributing to this Voronoi
        ! boundary line, and check IF they're already listed. If not, add
        ! them.

        DO vii = 1, r_Ac_nV( aci)

          ! Make sure we get the correct direction of integration
          IF (dir == 0) THEN
            vi_bot   =  r_Ac_vi_bot_left(  aci, vii)
            LI_xdy   =  r_Ac_LI_xdy(       aci, vii)
            LI_mxydx =  r_Ac_LI_mxydx(     aci, vii)
            LI_xydy  =  r_Ac_LI_xydy(      aci, vii)
          ELSE
            vi_bot   =  r_Ac_vi_bot_right( aci, vii)
            LI_xdy   = -r_Ac_LI_xdy(       aci, vii)
            LI_mxydx = -r_Ac_LI_mxydx(     aci, vii)
            LI_xydy  = -r_Ac_LI_xydy(      aci, vii)
          END IF

          IsListed = .FALSE.
          vii3     = 0
          DO vii2 = 1, r_nV( vi)+1
            IF (r_vi_bot( vi, vii2) == vi_bot) THEN
              ! It is listed here.
              IsListed = .TRUE.
              vii3 = vii2
              EXIT
            END IF
          END DO
          ! If it was not listed, increase number of contributing vertices by one.
          IF (.NOT. IsListed) THEN
            r_nV( vi) = r_nV( vi)+1
            vii3 = r_nV( vi)
          END IF

          ! Add contributions
          r_vi_bot(   vi, vii3) = vi_bot
          r_LI_xdy(   vi, vii3) = r_LI_xdy(   vi, vii3) + LI_xdy
          r_LI_mxydx( vi, vii3) = r_LI_mxydx( vi, vii3) + LI_mxydx
          r_LI_xydy(  vi, vii3) = r_LI_xydy(  vi, vii3) + LI_xydy
        END DO ! DO vii = 1, r_Ac%nV( aci)
      END DO ! DO ci = 1, mesh%nC( vi)

    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
  END SUBROUTINE rearrange_contributions_from_lines_to_vertices
  SUBROUTINE integrate_around_domain_boundary( mesh_top, mesh_bot, Vor_lines_top, Vor_lines_bot, r_dst_nV, r_dst_vi_src, r_dst_LI_xdy, r_dst_LI_mxydx, r_dst_LI_xydy)
    ! Integrate around the domain boundary - much easier to move this to a separate routine
    ! Only called by the Master; be careful not to include any sync statements!
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_top
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_bot
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_top
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_bot
    INTEGER,  DIMENSION(:  ),                INTENT(INOUT) :: r_dst_nV
    INTEGER,  DIMENSION(:,:),                INTENT(INOUT) :: r_dst_vi_src
    REAL(dp), DIMENSION(:,:),                INTENT(INOUT) :: r_dst_LI_xdy
    REAL(dp), DIMENSION(:,:),                INTENT(INOUT) :: r_dst_LI_mxydx
    REAL(dp), DIMENSION(:,:),                INTENT(INOUT) :: r_dst_LI_xydy
    
    ! Local variables:
    INTEGER                                                :: vi_top, vi_bot, vi_next_top, vi_next_bot, vi_top_int, vi_bot_int
    INTEGER                                                :: ci, vc, aci, vii, vii2, vi_src2
    REAL(dp), DIMENSION(2)                                 :: p, pc, ccr, ccl, p1, p_next_top, p_next_bot
    REAL(dp)                                               :: LI_xdy, LI_mxydx, LI_xydy
    LOGICAL                                                :: Finished, IsListed
    
  ! Southern boundary
  ! ==================

    vi_bot = 1
    vi_top = 1
    p = [mesh_bot%xmin, mesh_bot%ymin]

    Finished = .FALSE.
    DO WHILE (.NOT.Finished)

      ! Find the next vertex into whose Voronoi cell the domain boundary moves
      ! (this is not always the next adjacent edge vertex!)

      ! top
      vi_next_top = 0
      DO ci = 1, mesh_top%nC( vi_top)

        ! The Ac connection
        aci = mesh_top%iAci( vi_top, ci)

        ! The neighbour vertex
        IF (mesh_top%Aci( aci,1) == vi_top) THEN
          vc = mesh_top%Aci( aci,2)
        ELSE
          vc = mesh_top%Aci( aci,1)
        END IF
        pc  = mesh_top%V( vc,:)

        ! The endpoints of the shared Voronoi boundary
        ccr = [Vor_lines_top( aci,1,1), Vor_lines_top( aci,1,2)]
        ccl = [Vor_lines_top( aci,2,1), Vor_lines_top( aci,2,2)]    
        IF (NORM2(ccr-ccl) < mesh_top%tol_dist) CYCLE

        IF     ( ABS( ccr(2) - mesh_top%ymin) < mesh_top%tol_dist .AND. pc(1) > p(1)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further east - this is the one
          vi_next_top = vc
          p_next_top  = ccr
        ELSEIF ( ABS( ccl(2) - mesh_top%ymin) < mesh_top%tol_dist .AND. pc(1) > p(1)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further east - this is the one
          vi_next_top = vc
          p_next_top  = ccl
        END IF

      END DO ! DO ci = 1, mesh_top%nC( vi_bot)

      IF (vi_next_top == 0) THEN
        IF (vi_top == 2) THEN
          vi_next_top = mesh_top%C( vi_top, 1)
          p_next_top  = [mesh_top%xmax, mesh_top%ymin]
        ELSE
          WRITE(0,*) '    Remapping - integrate around domain boundary ERROR!'
          STOP
        END IF
      END IF

      ! bot
      vi_next_bot = 0
      DO ci = 1, mesh_bot%nC( vi_bot)

        ! The Ac connection
        aci = mesh_bot%iAci( vi_bot, ci)

        ! The neighbour vertex
        IF (mesh_bot%Aci( aci,1) == vi_bot) THEN
          vc = mesh_bot%Aci( aci,2)
        ELSE
          vc = mesh_bot%Aci( aci,1)
        END IF
        pc  = mesh_bot%V( vc,:)

        ! The endpoints of the shared Voronoi boundary
        ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
        ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]    
        IF (NORM2(ccr-ccl) < mesh_bot%tol_dist) CYCLE

        IF     ( ABS( ccr(2) - mesh_bot%ymin) < mesh_bot%tol_dist .AND. pc(1) > p(1)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further east - this is the one
          vi_next_bot = vc
          p_next_bot  = ccr
        ELSEIF ( ABS( ccl(2) - mesh_bot%ymin) < mesh_bot%tol_dist .AND. pc(1) > p(1)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further east - this is the one
          vi_next_bot = vc
          p_next_bot  = ccl
        END IF

      END DO ! DO ci = 1, mesh_bot%nC( vi_bot)

      IF (vi_next_bot == 0) THEN
        IF (vi_bot == 2) THEN
          vi_next_bot = mesh_bot%C( vi_bot, 1)
          p_next_bot  = [mesh_bot%xmax, mesh_bot%ymin]
        ELSE
          WRITE(0,*) '    Remapping - integrate around domain boundary ERROR!'
          STOP
        END IF
      END IF

      vi_bot_int = vi_bot
      vi_top_int = vi_top  

      IF (vi_bot==2 .AND. vi_top==2) THEN
        p_next_bot = [mesh_bot%xmax, mesh_bot%ymin]
        p_next_top = [mesh_bot%xmax, mesh_bot%ymin]
        Finished = .TRUE.
      END IF

      p1 = p

      IF (p_next_bot(1) < p_next_top(1) - mesh_top%tol_dist) THEN
        ! vi_bot comes first advance vi_bot
        p      = p_next_bot
        vi_bot = vi_next_bot
      ELSEIF (p_next_top(1) < p_next_bot(1) - mesh_top%tol_dist) THEN
        ! vi_top comes first advance vi_top
        p      = p_next_top
        vi_top = vi_next_top
      ELSE
        ! vi_bot and vi_top coincide advance both
        p      = p_next_top
        vi_bot = vi_next_bot
        vi_top = vi_next_top
      END IF

      ! Calculate the three %line integrals
      LI_xdy   = line_integral_xdy(   p1, p, mesh_top%tol_dist)
      LI_mxydx = line_integral_mxydx( p1, p, mesh_top%tol_dist)
      LI_xydy  = line_integral_xydy(  p1, p, mesh_top%tol_dist)

      ! Check if this contribution is already listed
      IsListed = .FALSE.
      vii2     = 0
      DO vii = 1, r_dst_nV( vi_top_int)+1
        vi_src2 = r_dst_vi_src( vi_top_int, vii)
        IF (vi_src2 == vi_bot_int) THEN
          ! It is listed here
          IsListed = .TRUE.
          vii2     = vii
          EXIT
        END IF
      END DO
      ! If it was not listed, increase number of contributing vertices by one.
      IF (.NOT. IsListed) THEN
        r_dst_nV( vi_top_int) = r_dst_nV( vi_top_int) + 1
        vii2 = r_dst_nV( vi_top_int)
      END IF

      ! Add contributions
      r_dst_vi_src(   vi_top_int, vii2) = vi_bot_int
      r_dst_LI_xdy(   vi_top_int, vii2) = r_dst_LI_xdy(   vi_top_int, vii2) + LI_xdy
      r_dst_LI_mxydx( vi_top_int, vii2) = r_dst_LI_mxydx( vi_top_int, vii2) + LI_mxydx
      r_dst_LI_xydy(  vi_top_int, vii2) = r_dst_LI_xydy(  vi_top_int, vii2) + LI_xydy

    END DO ! DO WHILE (.NOT. Finished)

  ! Eastern boundary
  ! =================

    vi_bot = 2
    vi_top = 2
    p = [mesh_bot%xmax, mesh_bot%ymin]

    Finished = .FALSE.
    DO WHILE (.NOT.Finished)

      ! Find the next vertex into whose Voronoi cell the domain boundary moves
      ! (this is not always the next adjacent edge vertex!)

      ! top
      vi_next_top = 0
      DO ci = 1, mesh_top%nC( vi_top)

        ! The Ac connection
        aci = mesh_top%iAci( vi_top, ci)

        ! The neighbour vertex
        IF (mesh_top%Aci( aci,1) == vi_top) THEN
          vc = mesh_top%Aci( aci,2)
        ELSE
          vc = mesh_top%Aci( aci,1)
        END IF
        pc  = mesh_top%V( vc,:)

        ! The endpoints of the shared Voronoi boundary
        ccr = [Vor_lines_top( aci,1,1), Vor_lines_top( aci,1,2)]
        ccl = [Vor_lines_top( aci,2,1), Vor_lines_top( aci,2,2)]    
        IF (NORM2(ccr-ccl) < mesh_top%tol_dist) CYCLE

        IF     ( ABS( ccr(1) - mesh_top%xmax) < mesh_top%tol_dist .AND. pc(2) > p(2)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further north - this is the one
          vi_next_top = vc
          p_next_top  = ccr
        ELSEIF ( ABS( ccl(1) - mesh_top%xmax) < mesh_top%tol_dist .AND. pc(2) > p(2)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further north - this is the one
          vi_next_top = vc
          p_next_top  = ccl
        END IF      

      END DO ! DO ci = 1, mesh_top%nC( vi_bot)

      IF (vi_next_top == 0) THEN
        IF (vi_top == 3) THEN
          vi_next_top = mesh_top%C( vi_top, 1)
          p_next_top  = [mesh_top%xmax, mesh_top%ymax]
        ELSE
          WRITE(0,*) '    Remapping - integrate around domain boundary ERROR!'
          STOP
        END IF
      END IF

      ! bot
      vi_next_bot = 0
      DO ci = 1, mesh_bot%nC( vi_bot)

        ! The Ac connection
        aci = mesh_bot%iAci( vi_bot, ci)

        ! The neighbour vertex
        IF (mesh_bot%Aci( aci,1) == vi_bot) THEN
          vc = mesh_bot%Aci( aci,2)
        ELSE
          vc = mesh_bot%Aci( aci,1)
        END IF
        pc  = mesh_bot%V( vc,:)

        ! The endpoints of the shared Voronoi boundary
        ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
        ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]    
        IF (NORM2(ccr-ccl) < mesh_bot%tol_dist) CYCLE

        IF     ( ABS( ccr(1) - mesh_bot%xmax) < mesh_bot%tol_dist .AND. pc(2) > p(2)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further north - this is the one
          vi_next_bot = vc
          p_next_bot  = ccr
        ELSEIF ( ABS( ccl(1) - mesh_bot%xmax) < mesh_bot%tol_dist .AND. pc(2) > p(2)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further north - this is the one
          vi_next_bot = vc
          p_next_bot  = ccl
        END IF

      END DO ! DO ci = 1, mesh_bot%nC( vi_bot)

      IF (vi_next_bot == 0) THEN
        IF (vi_bot == 3) THEN
          vi_next_bot = mesh_bot%C( vi_bot, 1)
          p_next_bot  = [mesh_bot%xmax, mesh_bot%ymax]
        ELSE
          WRITE(0,*) '    Remapping - integrate around domain boundary ERROR!'
          STOP
        END IF
      END IF

      vi_bot_int = vi_bot
      vi_top_int = vi_top  

      IF (vi_bot==3 .AND. vi_top==3) THEN
        p_next_bot = [mesh_bot%xmax, mesh_bot%ymax]
        p_next_top = [mesh_bot%xmax, mesh_bot%ymax]
        Finished = .TRUE.
      END IF

      p1 = p

      IF (p_next_bot(2) < p_next_top(2) - mesh_top%tol_dist) THEN
        ! vi_bot comes first advance vi_bot
        p      = p_next_bot
        vi_bot = vi_next_bot
      ELSEIF (p_next_top(2) < p_next_bot(2) - mesh_top%tol_dist) THEN
        ! vi_top comes first advance vi_top
        p      = p_next_top
        vi_top = vi_next_top
      ELSE
        ! vi_bot and vi_top coincide advance both
        p      = p_next_top
        vi_bot = vi_next_bot
        vi_top = vi_next_top
      END IF

      ! Calculate the three %line integrals
      LI_xdy   = line_integral_xdy(   p1, p, mesh_top%tol_dist)
      LI_mxydx = line_integral_mxydx( p1, p, mesh_top%tol_dist)
      LI_xydy  = line_integral_xydy(  p1, p, mesh_top%tol_dist)

      ! Check if this contribution is already listed
      IsListed = .FALSE.
      vii2     = 0
      DO vii = 1, r_dst_nV( vi_top_int)+1
        vi_src2 = r_dst_vi_src( vi_top_int, vii)
        IF (vi_src2 == vi_bot_int) THEN
          ! It is listed here
          IsListed = .TRUE.
          vii2     = vii
          EXIT
        END IF
      END DO
      ! If it was not listed, increase number of contributing vertices by one.
      IF (.NOT. IsListed) THEN
        r_dst_nV( vi_top_int) = r_dst_nV( vi_top_int) + 1
        vii2 = r_dst_nV( vi_top_int)
      END IF

      ! Add contributions
      r_dst_vi_src(   vi_top_int, vii2) = vi_bot_int
      r_dst_LI_xdy(   vi_top_int, vii2) = r_dst_LI_xdy(   vi_top_int, vii2) + LI_xdy
      r_dst_LI_mxydx( vi_top_int, vii2) = r_dst_LI_mxydx( vi_top_int, vii2) + LI_mxydx
      r_dst_LI_xydy(  vi_top_int, vii2) = r_dst_LI_xydy(  vi_top_int, vii2) + LI_xydy

    END DO ! DO WHILE (.NOT. Finished)

  ! Northern boundary
  ! ==================

    vi_bot = 3
    vi_top = 3
    p = [mesh_bot%xmax, mesh_bot%ymax]

    Finished = .FALSE.
    DO WHILE (.NOT.Finished)

      ! Find the next vertex into whose Voronoi cell the domain boundary moves
      ! (this is not always the next adjacent edge vertex!)

      ! top
      vi_next_top = 0
      DO ci = 1, mesh_top%nC( vi_top)

        ! The Ac connection
        aci = mesh_top%iAci( vi_top, ci)

        ! The neighbour vertex
        IF (mesh_top%Aci( aci,1) == vi_top) THEN
          vc = mesh_top%Aci( aci,2)
        ELSE
          vc = mesh_top%Aci( aci,1)
        END IF
        pc  = mesh_top%V( vc,:)

        ! The endpoints of the shared Voronoi boundary
        ccr = [Vor_lines_top( aci,1,1), Vor_lines_top( aci,1,2)]
        ccl = [Vor_lines_top( aci,2,1), Vor_lines_top( aci,2,2)]    
        IF (NORM2(ccr-ccl) < mesh_top%tol_dist) CYCLE

        IF     ( ABS( ccr(2) - mesh_top%ymax) < mesh_top%tol_dist .AND. pc(1) < p(1)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further west - this is the one
          vi_next_top = vc
          p_next_top  = ccr
        ELSEIF ( ABS( ccl(2) - mesh_top%ymax) < mesh_top%tol_dist .AND. pc(1) < p(1)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further west - this is the one
          vi_next_top = vc
          p_next_top  = ccl
        END IF      

      END DO ! DO ci = 1, mesh_top%nC( vi_bot)

      IF (vi_next_top == 0) THEN
        IF (vi_top == 4) THEN
          vi_next_top = mesh_top%C( vi_top, 1)
          p_next_top  = [mesh_top%xmin, mesh_top%ymax]
        ELSE
          WRITE(0,*) '    Remapping - integrate around domain boundary ERROR!'
          STOP
        END IF
      END IF

      ! bot
      vi_next_bot = 0
      DO ci = 1, mesh_bot%nC( vi_bot)

        ! The Ac connection
        aci = mesh_bot%iAci( vi_bot, ci)

        ! The neighbour vertex
        IF (mesh_bot%Aci( aci,1) == vi_bot) THEN
          vc = mesh_bot%Aci( aci,2)
        ELSE
          vc = mesh_bot%Aci( aci,1)
        END IF
        pc  = mesh_bot%V( vc,:)

        ! The endpoints of the shared Voronoi boundary
        ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
        ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]    
        IF (NORM2(ccr-ccl) < mesh_bot%tol_dist) CYCLE

        IF     ( ABS( ccr(2) - mesh_bot%ymax) < mesh_bot%tol_dist .AND. pc(1) < p(1)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further west - this is the one
          vi_next_bot = vc
          p_next_bot  = ccr
        ELSEIF ( ABS( ccl(2) - mesh_bot%ymax) < mesh_bot%tol_dist .AND. pc(1) < p(1)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further west - this is the one
          vi_next_bot = vc
          p_next_bot  = ccl
        END IF

      END DO ! DO ci = 1, mesh_bot%nC( vi_bot)

      IF (vi_next_bot == 0) THEN
        IF (vi_bot == 4) THEN
          vi_next_bot = mesh_bot%C( vi_bot, 1)
          p_next_bot  = [mesh_bot%xmin, mesh_bot%ymax]
        ELSE
          WRITE(0,*) '    Remapping - integrate around domain boundary ERROR!'
          STOP
        END IF
      END IF

      vi_bot_int = vi_bot
      vi_top_int = vi_top  

      IF (vi_bot==4 .AND. vi_top==4) THEN
        p_next_bot = [mesh_bot%xmin, mesh_bot%ymax]
        p_next_top = [mesh_bot%xmin, mesh_bot%ymax]
        Finished = .TRUE.
      END IF

      p1 = p

      IF (p_next_bot(1) > p_next_top(1) + mesh_top%tol_dist) THEN
        ! vi_bot comes first advance vi_bot
        p      = p_next_bot
        vi_bot = vi_next_bot
      ELSEIF (p_next_top(1) > p_next_bot(1) + mesh_top%tol_dist) THEN
        ! vi_top comes first advance vi_top
        p      = p_next_top
        vi_top = vi_next_top
      ELSE
        ! vi_bot and vi_top coincide advance both
        p      = p_next_top
        vi_bot = vi_next_bot
        vi_top = vi_next_top
      END IF

      ! Calculate the three %line integrals
      LI_xdy   = line_integral_xdy(   p1, p, mesh_top%tol_dist)
      LI_mxydx = line_integral_mxydx( p1, p, mesh_top%tol_dist)
      LI_xydy  = line_integral_xydy(  p1, p, mesh_top%tol_dist)

      ! Check if this contribution is already listed
      IsListed = .FALSE.
      vii2     = 0
      DO vii = 1, r_dst_nV( vi_top_int)+1
        vi_src2 = r_dst_vi_src( vi_top_int, vii)
        IF (vi_src2 == vi_bot_int) THEN
          ! It is listed here
          IsListed = .TRUE.
          vii2     = vii
          EXIT
        END IF
      END DO
      ! If it was not listed, increase number of contributing vertices by one.
      IF (.NOT. IsListed) THEN
        r_dst_nV( vi_top_int) = r_dst_nV( vi_top_int) + 1
        vii2 = r_dst_nV( vi_top_int)
      END IF

      ! Add contributions
      r_dst_vi_src(   vi_top_int, vii2) = vi_bot_int
      r_dst_LI_xdy(   vi_top_int, vii2) = r_dst_LI_xdy(   vi_top_int, vii2) + LI_xdy
      r_dst_LI_mxydx( vi_top_int, vii2) = r_dst_LI_mxydx( vi_top_int, vii2) + LI_mxydx
      r_dst_LI_xydy(  vi_top_int, vii2) = r_dst_LI_xydy(  vi_top_int, vii2) + LI_xydy

    END DO ! DO WHILE (.NOT. Finished)

  ! Western boundary
  ! =================

    vi_bot = 4
    vi_top = 4
    p = [mesh_bot%xmin, mesh_bot%ymax]

    Finished = .FALSE.
    DO WHILE (.NOT.Finished)

      ! Find the next vertex into whose Voronoi cell the domain boundary moves
      ! (this is not always the next adjacent edge vertex!)

      ! top
      vi_next_top = 0
      DO ci = 1, mesh_top%nC( vi_top)

        ! The Ac connection
        aci = mesh_top%iAci( vi_top, ci)

        ! The neighbour vertex
        IF (mesh_top%Aci( aci,1) == vi_top) THEN
          vc = mesh_top%Aci( aci,2)
        ELSE
          vc = mesh_top%Aci( aci,1)
        END IF
        pc  = mesh_top%V( vc,:)

        ! The endpoints of the shared Voronoi boundary
        ccr = [Vor_lines_top( aci,1,1), Vor_lines_top( aci,1,2)]
        ccl = [Vor_lines_top( aci,2,1), Vor_lines_top( aci,2,2)]    
        IF (NORM2(ccr-ccl) < mesh_top%tol_dist) CYCLE

        IF     ( ABS( ccr(1) - mesh_top%xmin) < mesh_top%tol_dist .AND. pc(2) < p(2)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further north - this is the one
          vi_next_top = vc
          p_next_top  = ccr
        ELSEIF ( ABS( ccl(1) - mesh_top%xmin) < mesh_top%tol_dist .AND. pc(2) < p(2)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further north - this is the one
          vi_next_top = vc
          p_next_top  = ccl
        END IF      

      END DO ! DO ci = 1, mesh_top%nC( vi_bot)

      IF (vi_next_top == 0) THEN
        IF (vi_top == 1) THEN
          vi_next_top = mesh_top%C( vi_top, 1)
          p_next_top  = [mesh_top%xmin, mesh_top%ymin]
        ELSE
          WRITE(0,*) '    Remapping - integrate around domain boundary ERROR!'
          STOP
        END IF
      END IF

      ! bot
      vi_next_bot = 0
      DO ci = 1, mesh_bot%nC( vi_bot)

        ! The Ac connection
        aci = mesh_bot%iAci( vi_bot, ci)

        ! The neighbour vertex
        IF (mesh_bot%Aci( aci,1) == vi_bot) THEN
          vc = mesh_bot%Aci( aci,2)
        ELSE
          vc = mesh_bot%Aci( aci,1)
        END IF
        pc  = mesh_bot%V( vc,:)

        ! The endpoints of the shared Voronoi boundary
        ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
        ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]    
        IF (NORM2(ccr-ccl) < mesh_bot%tol_dist) CYCLE

        IF     ( ABS( ccr(1) - mesh_bot%xmin) < mesh_bot%tol_dist .AND. pc(2) < p(2)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further north - this is the one
          vi_next_bot = vc
          p_next_bot  = ccr
        ELSEIF ( ABS( ccl(1) - mesh_bot%xmin) < mesh_bot%tol_dist .AND. pc(2) < p(2)) THEN
          ! This shared Voronoi boundary ends at the domain boundary, and the
          ! neighbour lies further north - this is the one
          vi_next_bot = vc
          p_next_bot  = ccl
        END IF

      END DO ! DO ci = 1, mesh_bot%nC( vi_bot)

      IF (vi_next_bot == 0) THEN
        IF (vi_bot == 1) THEN
          vi_next_bot = mesh_bot%C( vi_bot, 1)
          p_next_bot  = [mesh_bot%xmin, mesh_bot%ymin]
        ELSE
          WRITE(0,*) '    Remapping - integrate around domain boundary ERROR!'
          STOP
        END IF
      END IF

      vi_bot_int = vi_bot
      vi_top_int = vi_top  

      IF (vi_bot==1 .AND. vi_top==1) THEN
        p_next_bot = [mesh_bot%xmin, mesh_bot%ymin]
        p_next_top = [mesh_bot%xmin, mesh_bot%ymin]
        Finished = .TRUE.
      END IF

      p1 = p

      IF (p_next_bot(2) > p_next_top(2) + mesh_top%tol_dist) THEN
        ! vi_bot comes first advance vi_bot
        p      = p_next_bot
        vi_bot = vi_next_bot
      ELSEIF (p_next_top(2) > p_next_bot(2) + mesh_top%tol_dist) THEN
        ! vi_top comes first advance vi_top
        p      = p_next_top
        vi_top = vi_next_top
      ELSE
        ! vi_bot and vi_top coincide advance both
        p      = p_next_top
        vi_bot = vi_next_bot
        vi_top = vi_next_top
      END IF

      ! Calculate the three %line integrals
      LI_xdy   = line_integral_xdy(   p1, p, mesh_top%tol_dist)
      LI_mxydx = line_integral_mxydx( p1, p, mesh_top%tol_dist)
      LI_xydy  = line_integral_xydy(  p1, p, mesh_top%tol_dist)

      ! Check if this contribution is already listed
      IsListed = .FALSE.
      vii2     = 0
      DO vii = 1, r_dst_nV( vi_top_int)+1
        vi_src2 = r_dst_vi_src( vi_top_int, vii)
        IF (vi_src2 == vi_bot_int) THEN
          ! It is listed here
          IsListed = .TRUE.
          vii2     = vii
          EXIT
        END IF
      END DO
      ! If it was not listed, increase number of contributing vertices by one.
      IF (.NOT. IsListed) THEN
        r_dst_nV( vi_top_int) = r_dst_nV( vi_top_int) + 1
        vii2 = r_dst_nV( vi_top_int)
      END IF

      ! Add contributions
      r_dst_vi_src(   vi_top_int, vii2) = vi_bot_int
      r_dst_LI_xdy(   vi_top_int, vii2) = r_dst_LI_xdy(   vi_top_int, vii2) + LI_xdy
      r_dst_LI_mxydx( vi_top_int, vii2) = r_dst_LI_mxydx( vi_top_int, vii2) + LI_mxydx
      r_dst_LI_xydy(  vi_top_int, vii2) = r_dst_LI_xydy(  vi_top_int, vii2) + LI_xydy

    END DO ! DO WHILE (.NOT. Finished)
    
  END SUBROUTINE integrate_around_domain_boundary
  SUBROUTINE add_contributions_from_opposite_mesh( mesh_top, r_top_nV, r_top_vi_bot, r_top_LI_xdy, r_top_LI_mxydx, r_top_LI_xydy, &
    r_bot_nV, r_bot_vi_bot, r_bot_LI_xdy, r_bot_LI_mxydx, r_bot_LI_xydy)
    ! Add the line integral contributions from the botosite mesh to complete the loop integrals around the overlap regions
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_top
    INTEGER,  DIMENSION(:    ),              INTENT(INOUT) :: r_top_nV
    INTEGER,  DIMENSION(:,:  ),              INTENT(INOUT) :: r_top_vi_bot
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: r_top_LI_xdy
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: r_top_LI_mxydx
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: r_top_LI_xydy
    INTEGER,  DIMENSION(:    ),              INTENT(IN)    :: r_bot_nV
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: r_bot_vi_bot
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_bot_LI_xdy
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_bot_LI_mxydx
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_bot_LI_xydy
    
    ! Local variables:
    INTEGER                                                :: vi_top, vii_top, vi_bot, vii_bot, vii2_bot, vi_top_mirror
    LOGICAL                                                :: FoundIt

    ! Go over all listed contributing mesh_bot vertices, and add their
    ! contributions to this mesh_top vertex. This should complete the loop
    ! integrals for all mesh_top vertices that do not contain an entire
    ! mesh_bot Voronoi cell inside their own Voronoi cell.

    DO vi_top = mesh_top%v1, mesh_top%v2

      DO vii_top = 1, r_top_nV( vi_top)

        ! mesh_top vertex vi_top lists mesh_bot vertex vi_bot as a contribution in slot vii_top
        vi_bot = r_top_vi_bot( vi_top, vii_top)        

        ! Go over all mesh_top vertices listed as contributing to this
        ! mesh_bot vertex. If we find a match, add the contribution. If not, it
        ! might be that this mesh_ vertex' Voronoi cell is contained
        ! entirely inside the mesh_bot Voronoi cell. In that case, this
        ! mesh vertex should only list a single mesh_bot contribution.

        FoundIt  = .FALSE.
        vii2_bot = 0
        DO vii_bot = 1, r_bot_nV( vi_bot)
          ! mesh_bot vertex vi_bot lists mesh_top vertex vi_top_mirror as a contribution in slot vii_bot
          vi_top_mirror = r_bot_vi_bot( vi_bot, vii_bot)
          IF (vi_top_mirror == vi_top) THEN
            ! mesh_bot vertex vi_bot lists mesh_top vertex vi_top as a contribution in slot vii2_bot
            vii2_bot = vii_bot
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO

        ! Check
        IF (.NOT. FoundIt) THEN
          IF (r_top_nV( vi_top) > 1) THEN
            WRITE(0,*) '    Remapping - ERROR: couldnt find vi_top_mirror!'
            STOP
          END IF
        END IF

        IF (FoundIt) THEN
          ! Add the contribution
          r_top_LI_xdy(   vi_top, vii_top) = r_top_LI_xdy(   vi_top, vii_top) + r_bot_LI_xdy(   vi_bot, vii2_bot)
          r_top_LI_mxydx( vi_top, vii_top) = r_top_LI_mxydx( vi_top, vii_top) + r_bot_LI_mxydx( vi_bot, vii2_bot)
          r_top_LI_xydy(  vi_top, vii_top) = r_top_LI_xydy(  vi_top, vii_top) + r_bot_LI_xydy(  vi_bot, vii2_bot)
        END IF
        
      END DO ! DO vii_top = 1, r_top_nV( vi_top)

    END DO ! DO vi_top = mesh_top%v1, mesh_top%v2
    CALL sync
    
  END SUBROUTINE add_contributions_from_opposite_mesh
  SUBROUTINE finish_incomplete_vertices( mesh_top, mesh_bot, nV_max, Vor_lines_top, r_top_nV, r_top_vi_bot, r_top_LI_xdy, r_top_LI_mxydx, r_top_LI_xydy, &
    r_bot_nV, r_bot_vi_top, r_bot_LI_xdy, r_bot_LI_mxydx, r_bot_LI_xydy)
    ! When a mesh_bot Voronoi cell is enclosed inside a mesh_top Voronoi cell,
    ! that mesh_top vertex won't be listed as contributing to that mesh_top
    ! vertex in r_top. Although it will be listed in r_bot, searching through
    ! that array cannot be easily parallelised. A more efficient way is to find
    ! mesh_top vertices that are "incomplete" (i.e. the sum of LI_xdy
    ! contributions doesn't add up to the known Voronoi cell area). There
    ! should not be too many of them. For the ones that are there, use a
    ! flood-fill search to find all mesh_src vertices contained within them,
    ! find all of those that are complete enclosed, and add their contributions.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_top
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_bot
    INTEGER,                                 INTENT(IN)    :: nV_max
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_top
    INTEGER,  DIMENSION(:    ),              INTENT(INOUT) :: r_top_nV
    INTEGER,  DIMENSION(:,:  ),              INTENT(INOUT) :: r_top_vi_bot
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: r_top_LI_xdy
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: r_top_LI_mxydx
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: r_top_LI_xydy
    INTEGER,  DIMENSION(:    ),              INTENT(IN)    :: r_bot_nV
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: r_bot_vi_top
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_bot_LI_xdy
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_bot_LI_mxydx
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_bot_LI_xydy
    
    ! Local variables:
    INTEGER                                                :: vi_top, vi_bot, ci, vc_bot, vii, vii2, nV
    INTEGER                                                :: nV_bot_contained
    INTEGER,  DIMENSION(:  ), ALLOCATABLE                  :: vi_bot_contained
    REAL(dp), DIMENSION(2)                                 :: p, q
    REAL(dp)                                               :: Aint, Eint
    LOGICAL                                                :: IsListed

    nV_bot_contained = 0
    ALLOCATE( vi_bot_contained( nV_max))

    DO vi_top = mesh_top%v1, mesh_top%v2

      nV_bot_contained = 0
      vi_bot_contained = 0

      ! Check IF this mesh_top vertex is "complete"
      Aint = SUM( r_top_LI_xdy( vi_top, 1:r_top_nV( vi_top)))
      Eint = ABS( 1._dp - mesh_top%A( vi_top) / Aint)

      IF (Eint > 1e-6_dp) THEN
        ! This vertex is incomplete

        ! Find all contained mesh_bot vertices

        mesh_bot%VMap     = 0
        mesh_bot%VStack1  = 0
        mesh_bot%VStackN1 = 0

        ! Start with the one whose Voronoi cell contains vi_top, and its neighbours
        vi_bot = 5
        p = mesh_top%V( vi_top,:)
        CALL find_containing_vertex( mesh_bot, p, vi_bot)

        mesh_bot%VStackN1 = mesh_bot%VStackN1 + 1
        mesh_bot%VStack1( mesh_bot%VStackN1) = vi_bot
        mesh_bot%VMap( vi_bot) = 1

        DO ci = 1, mesh_bot%nC( vi_bot)
          vc_bot = mesh_bot%C( vi_bot, ci)
          mesh_bot%VStackN1 = mesh_bot%VStackN1 + 1
          mesh_bot%VStack1( mesh_bot%VStackN1) = vc_bot
          mesh_bot%VMap( vc_bot) = 1
        END DO

        ! Do a flood-fill search outward from these few seeds, add any vertices
        ! we find to be inside the mesh_top Voronoi cell to the list
        DO WHILE (mesh_bot%VStackN1 > 0)

          ! Take the last vertex from the stack
          vi_bot = mesh_bot%VStack1( mesh_bot%VStackN1)
          mesh_bot%VStackN1 = mesh_bot%VStackN1 - 1

          ! If it lies inside the mesh_top Voronoi cell, add its neighbours to the stack
          q = mesh_bot%V( vi_bot,:)
          IF (is_in_Voronoi_cell_remap( mesh_top, vi_top, Vor_lines_top, q)) THEN

            nV_bot_contained = nV_bot_contained + 1
            vi_bot_contained( nV_bot_contained) = vi_bot

            DO ci = 1, mesh_bot%nC( vi_bot)
              vc_bot = mesh_bot%C( vi_bot, ci)
              IF (mesh_bot%VMap( vc_bot)==0) THEN
                mesh_bot%VStackN1 = mesh_bot%VStackN1 + 1
                mesh_bot%VStack1( mesh_bot%VStackN1) = vc_bot
                mesh_bot%VMap( vc_bot) = 1
              END IF
            END DO
            
          END IF ! IF (is_in_Voronoi_cell_remap( mesh_top, vi_top, Vor_lines_top, q)) THEN
        END DO ! DO WHILE (mesh_bot%VStackN > 0)

        ! Go over all mesh_bot vertices that lie inside this mesh_top Voronoi cell.
        ! If we find one that only lists vi_top as a contribution, that means
        ! it is complete enclosed in this mesh_top Voronoi cell. In that case,
        ! add its contributions to vi_top%

        DO vii = 1, nV_bot_contained
          vi_bot = vi_bot_contained( vii)

          ! Exceptions for edge vertices enclosed inside edge vertices - those
          ! have already been picked up by the integration around the domain boundary
          IF (mesh_top%edge_index(vi_top)>0 .AND. mesh_bot%edge_index(vi_bot)>0) CYCLE

          IF (r_bot_nV( vi_bot) == 1 .AND. r_bot_vi_top( vi_bot, 1) == vi_top) THEN
            ! This vertex only lists vi_top as a contributor, so it is
            ! enclosed. However, if they share an edge, it might already be
            ! listed. Check for this.

            IsListed = .FALSE.
            DO vii2 = 1, r_top_nV( vi_top)
              IF (r_top_vi_bot( vi_top, vii2) == vi_bot) THEN
                IsListed = .TRUE.
                EXIT
              END IF
            END DO
            IF (IsListed) CYCLE

            ! Add the contribution
            r_top_nV( vi_top) = r_top_nV( vi_top) + 1
            nV = r_top_nV( vi_top)
            r_top_vi_bot(   vi_top, nV) = vi_bot
            r_top_LI_xdy(   vi_top, nV) = r_bot_LI_xdy(   vi_bot, 1)
            r_top_LI_mxydx( vi_top, nV) = r_bot_LI_mxydx( vi_bot, 1)
            r_top_LI_xydy(  vi_top, nV) = r_bot_LI_xydy(  vi_bot, 1)
          END IF
        END DO ! DO vii = 1, nV_bot_contained

      END IF ! IF (Eint > 1e-6_dp) THEN
      
    END DO ! DO vi_top = mesh_top%v1, mesh_top%v2
    CALL sync

    ! Check if all vertices really are complete now
    ! =============================================
  
    DO vi_top = mesh_top%v1, mesh_top%v2
    
      ! Ignore errors on the boundary, that's too tricky to get right,
      ! and not really important since there will be no ice there anyway.
      IF (mesh_top%edge_index(vi_top) > 0) CYCLE

      ! Check if this mesh_top vertex is "complete"
      Aint = SUM( r_top_LI_xdy( vi_top, 1:r_top_nV( vi_top)))
      Eint = ABS( 1._dp - mesh_top%A( vi_top) / Aint)

      IF (Eint > 1e-3_dp) THEN
        ! This vertex is incomplete    
        WRITE(0,*) '    Remapping - ERROR: vi = ', vi_top, ' is still incomplete!'
      END IF
      
    END DO
    CALL sync
    
  END SUBROUTINE finish_incomplete_vertices
  SUBROUTINE calculate_remapping_weights_from_line_integrals( mesh_top, mesh_bot, r_top_nV, r_top_vi, r_top_LI_xdy, r_top_LI_mxydx, r_top_LI_xydy, map)
    ! Calculate the actual remapping weights from the three line integrals
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_top
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_bot
    INTEGER,  DIMENSION(:    ),              INTENT(IN)    :: r_top_nV
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: r_top_vi
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_top_LI_xdy
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_top_LI_mxydx
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: r_top_LI_xydy
    TYPE(type_remapping_conservative),       INTENT(INOUT) :: map
    
    ! Local variables:
    INTEGER                                                :: vi, vvi, vi_opp

    map%w0(  mesh_top%v1:mesh_top%v2,:) = 0._dp
    map%w1x( mesh_top%v1:mesh_top%v2,:) = 0._dp
    map%w1y( mesh_top%v1:mesh_top%v2,:) = 0._dp

    map%nV(  mesh_top%v1:mesh_top%v2  ) = r_top_nV(  mesh_top%v1:mesh_top%v2)
    map%vi(  mesh_top%v1:mesh_top%v2,:) = r_top_vi(  mesh_top%v1:mesh_top%v2,:)

    DO vi = mesh_top%v1, mesh_top%v2
      DO vvi = 1, map%nV( vi)
        vi_opp = map%vi( vi, vvi)
        map%w0(  vi, vvi) = (r_top_LI_xdy(   vi, vvi) / mesh_top%A(vi))
        map%w1x( vi, vvi) = (r_top_LI_mxydx( vi, vvi) / mesh_top%A(vi)) - (mesh_bot%VorGC( vi_opp,1) * map%w0( vi, vvi))
        map%w1y( vi, vvi) = (r_top_LI_xydy(  vi, vvi) / mesh_top%A(vi)) - (mesh_bot%VorGC( vi_opp,2) * map%w0( vi, vvi))
      END DO
    END DO
    CALL sync
    
  END SUBROUTINE calculate_remapping_weights_from_line_integrals
  SUBROUTINE check_if_remapping_is_conservative( mesh_src, mesh_dst, map)
    ! Check if the remapping arrays really are conservative
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_src ! INOUT instead of IN because the search maps and stacks are used
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst
    TYPE(type_remapping_conservative),       INTENT(INOUT) :: map
    
    ! Local variables
    REAL(dp), DIMENSION(:), POINTER                        ::  d_src,  d_dst
    INTEGER                                                :: wd_src, wd_dst
    INTEGER                                                :: vi_dst, n_wrong, ti_src, via_src, vib_src, vic_src
    REAL(dp), DIMENSION(2)                                 :: p, pa, pb, pc
    REAL(dp)                                               :: Atot, Aa, Ab, Ac
       
    CALL allocate_shared_dp_1D( mesh_src%nV, d_src, wd_src)
    CALL allocate_shared_dp_1D( mesh_dst%nV, d_dst, wd_dst)
    
    d_src(mesh_src%v1:mesh_src%v2) = 1._dp
    CALL sync
    CALL remap_cons_2nd_order_2D( mesh_src, mesh_dst, map, d_src, d_dst)
    CALL sync
    
    ! Ignore edge vertices, as they're useless, and maybe 99% of remapping errors occur there
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      IF (mesh_dst%edge_index( vi_dst) > 0) d_dst( vi_dst) = 1._dp
    END DO
    CALL sync
    
    ! Check for vertices where conservative remapping went wrong. If there are not too many, replace them with linear interpolation.
    n_wrong = 0
    DO vi_dst = 1, mesh_dst%nV
      IF (ABS(1._dp - d_dst( vi_dst)) > 1E-1_dp) THEN
      
        n_wrong = n_wrong + 1
        
        map%nV(  vi_dst  ) = 3
        map%vi(  vi_dst,:) = 0
        map%w0(  vi_dst,:) = 0._dp
        map%w1x( vi_dst,:) = 0._dp
        map%w1y( vi_dst,:) = 0._dp
        
        ti_src = 1
        p      = mesh_dst%V( vi_dst, :)
        CALL find_containing_triangle( mesh_src, p, ti_src)
        
        ! Find the three vertices spanning this triangle
        via_src = mesh_src%Tri( ti_src,1)
        vib_src = mesh_src%Tri( ti_src,2)
        vic_src = mesh_src%Tri( ti_src,3)
        
        pa  = mesh_src%V( via_src,:)
        pb  = mesh_src%V( vib_src,:)
        pc  = mesh_src%V( vic_src,:)
        
        ! Calculate their interpolation weights
        Atot = find_triangle_area( pa, pb, pc)
        Aa   = find_triangle_area( pb, pc, p )
        Ab   = find_triangle_area( pc, pa, p )
        Ac   = find_triangle_area( pa, pb, p )
        
        map%vi( vi_dst,1:3) = [via_src, vib_src, vic_src]
        map%w0( vi_dst,1:3) = [Aa/Atot, Ab/Atot, Ac/Atot]        
                
      END IF      
    END DO ! DO vi_dst = 1, mesh_dst%nV   
    CALL sync
    
    IF (n_wrong > CEILING(REAL(mesh_dst%nV) / 100._dp)) THEN
      IF (par%master) WRITE(0,'(A,I6,A,I6,A)') '    ERROR: 2nd order conservative remapping goes wrong for ', n_wrong, ' out of ', mesh_dst%nV, ' vertices!'
      STOP
    END IF
    
!    ! Check for conservation and monotony
!    int_src = SUM(d_src * mesh_src%A)
!    int_dst = SUM(d_dst * mesh_dst%A)
!    Eint = ABS(1._dp - int_dst / int_src)
!    Emon = MAXVAL(ABS(1._dp - d_dst))
!    IF (Eint > 1E-3_dp) THEN  
!      IF (par%master) WRITE(0,'(A,E9.2)') '    ERROR: 2nd order conservative remapping isnt conservative! Integral error = ', Eint
!      STOP
!    ELSEIF (Emon > 1E-1_dp) THEN    
!      IF (par%master) WRITE(0,'(A,E9.2)') '    ERROR: 2nd order conservative remapping isnt monotonous! Maximum error = ', Emon
!      STOP
!    ELSE
!      IF (par%master) WRITE(0,'(A,E9.2,A,E9.2,A)') '    2nd order conservative remapping is conservative to ', Eint, ' (maximum local error: ', Emon, ')'
!    END IF
!    CALL sync     
    
    CALL deallocate_shared( wd_src)
    CALL deallocate_shared( wd_dst)    
    NULLIFY( d_src)
    NULLIFY( d_dst)
    
  END SUBROUTINE check_if_remapping_is_conservative
  SUBROUTINE create_remapping_arrays_conservative_deallocate_memory( &
    proc_domain_Ac_src, proc_domain_Ac_dst, &
    wVor_lines_src, wVor_lines_dst, wVor_vi_ti_src, wVor_vi_ti_dst, &
    wr_Ac_src_nV, wr_Ac_src_vi_dst_left, wr_Ac_src_vi_dst_right, wr_Ac_src_LI_xdy, wr_Ac_src_LI_mxydx, wr_Ac_src_LI_xydy, &
    wr_Ac_dst_nV, wr_Ac_dst_vi_src_left, wr_Ac_dst_vi_src_right, wr_Ac_dst_LI_xdy, wr_Ac_dst_LI_mxydx, wr_Ac_dst_LI_xydy, &
    wr_src_nV, wr_src_vi_dst, wr_src_LI_xdy, wr_src_LI_mxydx, wr_src_LI_xydy, &
    wr_dst_nV, wr_dst_vi_src, wr_dst_LI_xdy, wr_dst_LI_mxydx, wr_dst_LI_xydy)
    
    ! Clean up after yourself: deallocate memory for all the temporary arrays used in conservative remapping.
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,  DIMENSION(:    ), ALLOCATABLE, INTENT(INOUT) :: proc_domain_Ac_src, proc_domain_Ac_dst
    INTEGER,                                 INTENT(INOUT) :: wVor_lines_src, wVor_lines_dst
    INTEGER,                                 INTENT(INOUT) :: wVor_vi_ti_src, wVor_vi_ti_dst
    INTEGER,                                 INTENT(INOUT) :: wr_Ac_src_nV, wr_Ac_dst_nV
    INTEGER,                                 INTENT(INOUT) :: wr_Ac_src_vi_dst_left, wr_Ac_src_vi_dst_right, wr_Ac_dst_vi_src_left, wr_Ac_dst_vi_src_right
    INTEGER,                                 INTENT(INOUT) :: wr_Ac_src_LI_xdy, wr_Ac_src_LI_mxydx, wr_Ac_src_LI_xydy, wr_Ac_dst_LI_xdy, wr_Ac_dst_LI_mxydx, wr_Ac_dst_LI_xydy
    INTEGER,                                 INTENT(INOUT) :: wr_src_nV, wr_dst_nV
    INTEGER,                                 INTENT(INOUT) :: wr_src_vi_dst, wr_dst_vi_src
    INTEGER,                                 INTENT(INOUT) :: wr_src_LI_xdy, wr_src_LI_mxydx, wr_src_LI_xydy, wr_dst_LI_xdy, wr_dst_LI_mxydx, wr_dst_LI_xydy
    
    DEALLOCATE( proc_domain_Ac_src)
    DEALLOCATE( proc_domain_Ac_dst)   
            
    CALL deallocate_shared( wVor_lines_src)
    CALL deallocate_shared( wVor_lines_dst)
    CALL deallocate_shared( wVor_vi_ti_src)
    CALL deallocate_shared( wVor_vi_ti_dst) 
    
    CALL deallocate_shared( wr_Ac_src_nV          )
    CALL deallocate_shared( wr_Ac_src_vi_dst_left )
    CALL deallocate_shared( wr_Ac_src_vi_dst_right)
    CALL deallocate_shared( wr_Ac_src_LI_xdy      )
    CALL deallocate_shared( wr_Ac_src_LI_mxydx    )
    CALL deallocate_shared( wr_Ac_src_LI_xydy     )
    
    CALL deallocate_shared( wr_Ac_dst_nV          )
    CALL deallocate_shared( wr_Ac_dst_vi_src_left )
    CALL deallocate_shared( wr_Ac_dst_vi_src_right)
    CALL deallocate_shared( wr_Ac_dst_LI_xdy      )
    CALL deallocate_shared( wr_Ac_dst_LI_mxydx    )
    CALL deallocate_shared( wr_Ac_dst_LI_xydy     )
    
    CALL deallocate_shared( wr_src_nV      )
    CALL deallocate_shared( wr_src_vi_dst  )
    CALL deallocate_shared( wr_src_LI_xdy  )
    CALL deallocate_shared( wr_src_LI_mxydx)
    CALL deallocate_shared( wr_src_LI_xydy )
     
    CALL deallocate_shared( wr_dst_nV      )
    CALL deallocate_shared( wr_dst_vi_src  )
    CALL deallocate_shared( wr_dst_LI_xdy  )
    CALL deallocate_shared( wr_dst_LI_mxydx)
    CALL deallocate_shared( wr_dst_LI_xydy ) 
    
  END SUBROUTINE create_remapping_arrays_conservative_deallocate_memory
  
  SUBROUTINE calculate_line_integral_contributions( mesh_bot, Vor_lines_top, Vor_lines_bot, Vor_vi_ti_bot, aci, &
          r_sng_nV, r_sng_vi_left, r_sng_vi_right, r_sng_LI_xdy, r_sng_LI_mxydx, r_sng_LI_xydy, vi_bot, CountCoincidences, nV_max)
    ! Trace Voronoi boundary line aci through mesh_bot, calculate the three line integrals along all sections of this line.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_bot
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_top
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_bot
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: Vor_vi_ti_bot
    INTEGER,                                 INTENT(IN)    :: aci
    INTEGER,                                 INTENT(OUT)   :: r_sng_nV
    INTEGER,  DIMENSION(:    ),              INTENT(OUT)   :: r_sng_vi_left
    INTEGER,  DIMENSION(:    ),              INTENT(OUT)   :: r_sng_vi_right
    REAL(dp), DIMENSION(:    ),              INTENT(OUT)   :: r_sng_LI_xdy
    REAL(dp), DIMENSION(:    ),              INTENT(OUT)   :: r_sng_LI_mxydx
    REAL(dp), DIMENSION(:    ),              INTENT(OUT)   :: r_sng_LI_xydy
    INTEGER,                                 INTENT(INOUT) :: vi_bot
    LOGICAL,                                 INTENT(IN)    :: CountCoincidences
    INTEGER,                                 INTENT(IN)    :: nV_max
    
    ! Local variables:
    REAL(dp), DIMENSION(2)                                 :: p, q, p_next
    LOGICAL                                                :: Finished, DoAdd, Coincides
    INTEGER                                                :: aci_coincide, tri_coincide, vi_in, vi_bot_left, vi_bot_right
    INTEGER                                                :: ncycle

    ! Initialise results
    r_sng_nV       = 0
    r_sng_vi_left  = 0
    r_sng_vi_right = 0
    r_sng_LI_xdy   = 0._dp
    r_sng_LI_mxydx = 0._dp
    r_sng_LI_xydy  = 0._dp

    ! The start and end points of the line we must trace
    p = [Vor_lines_top( aci, 1, 1), Vor_lines_top( aci, 1, 2)]
    q = [Vor_lines_top( aci, 2, 1), Vor_lines_top( aci, 2, 2)]

    ! In rare cases, p and q coincide. In that case, don't bother.
    IF (NORM2(p-q) < mesh_bot%tol_dist) RETURN

  !  WRITE(0,*) '  Process ', par%i, ' - Integrating over the Voronoi boundary between vertices ', mesh_top%Aci(aci,1), ' and ', mesh_top%Aci(aci,2), ' (aci = ', aci, ')'

    ! The mesh vertex whose Voronoi cell contains the start point (might
    ! coincide with a line or vertex of that Voronoi cell's boundary!)
    CALL find_containing_vertex( mesh_bot, p, vi_bot)

    ! Find a proper description of the starting point
    CALL trace_line_through_mesh_start( mesh_bot, p, vi_bot, Vor_lines_bot, aci_coincide, tri_coincide, vi_in)

    ! Trace the line through the mesh
    Finished = .FALSE.
    ncycle   = 0
    DO WHILE (.NOT. Finished)

      ! Find the point p_next where pq crosses into a new Voronoi cell
      CALL trace_line_through_mesh( mesh_bot, p, q, Vor_lines_bot, Vor_vi_ti_bot, aci_coincide, tri_coincide, vi_in, p_next, vi_bot_left, vi_bot_right, Coincides, Finished)

      ! Integrate over the line section p-p_next
      DoAdd = .TRUE.
      IF (Coincides .AND. (.NOT. CountCoincidences)) DoAdd = .FALSE.
      IF (DoAdd) THEN
        CALL add_line_integral_contributions( p, p_next, vi_bot_left, vi_bot_right, r_sng_nV, r_sng_vi_left, r_sng_vi_right, r_sng_LI_xdy, r_sng_LI_mxydx, r_sng_LI_xydy, mesh_bot%tol_dist)
      END IF

      ! Cycle the pointer
      p = p_next

      ! Safety
      ncycle = ncycle + 1
      IF (ncycle > nV_max) THEN
        WRITE(0,*) '   Remapping - ERROR: line tracer got stuck!'
        STOP
      END IF

    END DO ! DO WHILE (.NOT. Finished)
          
  END SUBROUTINE calculate_line_integral_contributions  
  SUBROUTINE add_line_integral_contributions( p, q, vi_bot_left, vi_bot_right, r_sng_nV, r_sng_vi_left, r_sng_vi_right, r_sng_LI_xdy, r_sng_LI_mxydx, r_sng_LI_xydy, tol)
    ! Add the line integral contributions from a line segment to the list of contributions for the entire line
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: q
    INTEGER,                                 INTENT(IN)    :: vi_bot_left
    INTEGER,                                 INTENT(IN)    :: vi_bot_right
    INTEGER,                                 INTENT(INOUT) :: r_sng_nV
    INTEGER,  DIMENSION(:  ),                INTENT(INOUT) :: r_sng_vi_left
    INTEGER,  DIMENSION(:  ),                INTENT(INOUT) :: r_sng_vi_right
    REAL(dp), DIMENSION(:  ),                INTENT(INOUT) :: r_sng_LI_xdy
    REAL(dp), DIMENSION(:  ),                INTENT(INOUT) :: r_sng_LI_mxydx
    REAL(dp), DIMENSION(:  ),                INTENT(INOUT) :: r_sng_LI_xydy
    REAL(dp),                                INTENT(IN)    :: tol
    
    ! Local variables:
    INTEGER                                                :: nV
    
    r_sng_nV = r_sng_nV + 1
    nV = r_sng_nV
    
    r_sng_vi_left(  nV) = vi_bot_left
    r_sng_vi_right( nV) = vi_bot_right
    r_sng_LI_xdy(   nV) = line_integral_xdy(   p, q, tol)
    r_sng_LI_mxydx( nV) = line_integral_mxydx( p, q, tol)
    r_sng_LI_xydy(  nV) = line_integral_xydy(  p, q, tol)
    
  END SUBROUTINE add_line_integral_contributions
  
  SUBROUTINE trace_line_through_mesh_start(        mesh_bot, p, vi_bot, Vor_lines_bot, aci_coincide, tri_coincide, vi_in)
    ! Given a point p that lies in the interior, or on the boundary, of Voronoi cell vi_bot in mesh_bot, find out if p:
    ! - coincides with Voronoi boundary line aci_coincide
    ! -                triangle circumcenter (= Voronoi boundary vertex) tri_coincide
    ! - lies in the interior of Voronoi cell vi_in (= vi_bot if so, or = 0 if not)
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_bot
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p
    INTEGER,                                 INTENT(IN)    :: vi_bot
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_bot
    INTEGER,                                 INTENT(OUT)   :: aci_coincide
    INTEGER,                                 INTENT(OUT)   :: tri_coincide
    INTEGER,                                 INTENT(OUT)   :: vi_in
    
    ! Local variables:
    INTEGER                                                :: iti, ti, ci, aci
    REAL(dp), DIMENSION(2)                                 :: pa, pb
    
    aci_coincide = 0
    tri_coincide = 0
    vi_in        = 0

    ! Check if p lies on one of the Voronoi cell vertices (i.e. circumcenters of surrounding triangles)
    DO iti = 1, mesh_bot%niTri( vi_bot)
      ti = mesh_bot%iTri( vi_bot, iti)
      IF (NORM2( p - mesh_bot%Tricc( ti,:)) < mesh_bot%tol_dist) THEN
        ! p coincide with the circumcenter of triangle ti
        tri_coincide = ti
        RETURN
      END IF
    END DO

    ! Check if p lies on one of the Voronoi cell lines
    DO ci = 1, mesh_bot%nC( vi_bot)
      aci = mesh_bot%iAci( vi_bot, ci)
      pa  = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
      pb  = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]
      IF (NORM2(pa-pb) < mesh_bot%tol_dist) CYCLE
      IF (lies_on_line_segment( pa, pb, p, mesh_bot%tol_dist) .OR. NORM2(p-pa) < mesh_bot%tol_dist .OR. NORM2(p-pb) < mesh_bot%tol_dist) THEN
        ! p coincide with the cc-cc line specIFied by Ac vertex aci
        aci_coincide = aci
        RETURN
      END IF
    END DO

    ! If it doesn't lie on one of the lines or vertices of the Voronoi cell
    ! boundary, then it must lie in the interior
    vi_in = vi_bot
    
  END SUBROUTINE trace_line_through_mesh_start
  SUBROUTINE trace_line_through_mesh(              mesh_bot, p, q, Vor_lines_bot, Vor_vi_ti_bot, aci_coincide, tri_coincide, vi_in, p_next, vi_bot_left, vi_bot_right, Coincides, Finished)
    ! Given a point p that either:
    ! - coincides with mesh_bot Voronoi boundary line aci_coincide,
    ! - coincides with mesh_bot triangle circumcenter (= Voronoi boundary vertex) tri_coincide, or
    ! - lies in the interior of Voronoi cell vi_in,
    ! find the point p_next where the line pq crosses into a new Voronoi cell, and return indices
    ! of the mesh_bot vertices vi_bot_left and vi_bot_right whose Voronoi cells lie to the left and
    ! right of this line, respectively. If these are two different ones, pq Coincides with a Voronoi
    ! boundary line of mesh_bot. If q lies inside, or on the boundary of, the same Voronoi cell as p,
    ! we are Finished.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_bot
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: q
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_bot
    INTEGER,  DIMENSION(:,:),                INTENT(IN)    :: Vor_vi_ti_bot
    INTEGER,                                 INTENT(INOUT) :: aci_coincide
    INTEGER,                                 INTENT(INOUT) :: tri_coincide
    INTEGER,                                 INTENT(INOUT) :: vi_in
    REAL(dp), DIMENSION(2),                  INTENT(OUT)   :: p_next
    INTEGER,                                 INTENT(OUT)   :: vi_bot_left
    INTEGER,                                 INTENT(OUT)   :: vi_bot_right
    LOGICAL,                                 INTENT(OUT)   :: Coincides
    LOGICAL,                                 INTENT(OUT)   :: Finished
    
    IF (tri_coincide > 0) THEN
      ! p coincides with the circumcenter of mesh_bot triangle tri_coincide
      CALL trace_line_through_mesh_tri_coincide( mesh_bot, p, q, Vor_lines_bot, Vor_vi_ti_bot, aci_coincide, tri_coincide, vi_in, p_next, vi_bot_left, vi_bot_right, Coincides, Finished)
    ELSEIF (aci_coincide > 0) THEN
      ! p coincides with the mesh_bot triangle cc-cc line specified by the mesh_bot Ac vertex aci_coincide
      CALL trace_line_through_mesh_aci_coincide( mesh_bot, p, q, Vor_lines_bot, Vor_vi_ti_bot, aci_coincide, tri_coincide, vi_in, p_next, vi_bot_left, vi_bot_right, Coincides, Finished)
    ELSE
      ! p lies inside the interior of mesh_bot vertex vi_in
      CALL trace_line_through_mesh_vi_interior(  mesh_bot, p, q, Vor_lines_bot, Vor_vi_ti_bot, aci_coincide, tri_coincide, vi_in, p_next, vi_bot_left, vi_bot_right, Coincides, Finished)
    END IF
    
  END SUBROUTINE trace_line_through_mesh
  SUBROUTINE trace_line_through_mesh_aci_coincide( mesh_bot, p, q, Vor_lines_bot, Vor_vi_ti_bot, aci_coincide, tri_coincide, vi_in, p_next, vi_bot_left, vi_bot_right, Coincides, Finished)
    ! Given a point p that either:
    ! - coincides with mesh_bot Voronoi boundary line aci_coincide,
    ! - coincides with mesh_bot triangle circumcenter (= Voronoi boundary vertex) tri_coincide, or
    ! - lies in the interior of Voronoi cell vi_in,
    ! find the point p_next where the line pq crosses into a new Voronoi cell, and return indices
    ! of the mesh_bot vertices vi_bot_left and vi_bot_right whose Voronoi cells lie to the left and
    ! right of this line, respectively. If these are two different ones, pq Coincides with a Voronoi
    ! boundary line of mesh_bot. If q lies inside, or on the boundary of, the same Voronoi cell as p,
    ! we are Finished.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_bot
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: q
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_bot
    INTEGER,  DIMENSION(:,:),                INTENT(IN)    :: Vor_vi_ti_bot
    INTEGER,                                 INTENT(INOUT) :: aci_coincide
    INTEGER,                                 INTENT(INOUT) :: tri_coincide
    INTEGER,                                 INTENT(INOUT) :: vi_in
    REAL(dp), DIMENSION(2),                  INTENT(OUT)   :: p_next
    INTEGER,                                 INTENT(OUT)   :: vi_bot_left
    INTEGER,                                 INTENT(OUT)   :: vi_bot_right
    LOGICAL,                                 INTENT(OUT)   :: Coincides
    LOGICAL,                                 INTENT(OUT)   :: Finished
    
    ! Local variables:
    INTEGER                                                :: via, vib, vir, vil, tir, til, ci, aci, ti1, ti2, vi1, vi2, vic
    REAL(dp), DIMENSION(2)                                 :: ccr, ccl, cc1, cc2, llis
    LOGICAL                                                :: do_cross

    ! Safety
    IF (aci_coincide==0 .OR.tri_coincide>0 .OR.vi_in>0) THEN
      WRITE(0,*) '    Remapping - ERROR: coincidence indicators dont match!'
      STOP
    END IF

    ! Indices of the four vertices and two triangles describing aci
    via = Vor_vi_ti_bot( aci_coincide, 1)
    vib = Vor_vi_ti_bot( aci_coincide, 2)
    vir = Vor_vi_ti_bot( aci_coincide, 3)
    vil = Vor_vi_ti_bot( aci_coincide, 4)
    tir = Vor_vi_ti_bot( aci_coincide, 5)
    til = Vor_vi_ti_bot( aci_coincide, 6)

    ! Coordinates of the two triangle circumcenters (endpoints of the Voronoi boundary line)
    ccr = [Vor_lines_bot( aci_coincide,1,1), Vor_lines_bot( aci_coincide,1,2)]
    ccl = [Vor_lines_bot( aci_coincide,2,1), Vor_lines_bot( aci_coincide,2,2)]

    ! Safety
    IF (.NOT. lies_on_line_segment( ccr, ccl, p, mesh_bot%tol_dist)) THEN
      WRITE(0,*) '    Remapping - ERROR: coincidence indicators dont match!'
      STOP
    END IF

    ! Initialise output
    aci_coincide = 0
    tri_coincide = 0
    vi_in        = 0
    p_next       = p
    vi_bot_left  = 0
    vi_bot_right = 0
    Coincides    = .FALSE.
    Finished     = .FALSE.

  ! Check IF q lies on either of the two endpoints of the Voronoi boundary line
  ! ===========================================================================
    
    IF     (NORM2( q-ccr) < mesh_bot%tol_dist) THEN
      ! q lies on ccr
      p_next       = q
      vi_bot_left  = vib
      vi_bot_right = via
      Coincides    = .TRUE.
      Finished     = .TRUE.
      RETURN
    ELSEIF (NORM2( q-ccl) < mesh_bot%tol_dist) THEN
      ! q lies on ccl
      p_next       = q
      vi_bot_left  = via
      vi_bot_right = vib
      Coincides    = .TRUE.
      Finished     = .TRUE.
      RETURN
    END IF

  ! Check IF q lies on the line from p to ccl or ccr
  ! ================================================
  
    IF     (lies_on_line_segment( p, ccr, q, mesh_bot%tol_dist) .AND. NORM2(p-ccr) > mesh_bot%tol_dist) THEN
      ! q lies on the line from p to ccr
      p_next       = q
      vi_bot_left  = vib
      vi_bot_right = via
      Coincides    = .TRUE.
      Finished     = .TRUE.
      RETURN
    ELSEIF (lies_on_line_segment( p, ccl, q, mesh_bot%tol_dist) .AND. NORM2(p-ccl) > mesh_bot%tol_dist) THEN
      ! q lies on the line from p to ccl
      p_next       = q
      vi_bot_left  = via
      vi_bot_right = vib
      Coincides    = .TRUE.
      Finished     = .TRUE.
      RETURN
    END IF

  ! Check IF q lies on the boundary of the Voronoi cell of via
  ! ==========================================================
  
    DO ci = 1, mesh_bot%nC( via)
      aci = mesh_bot%iAci( via, ci)
      cc1 = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
      cc2 = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]
      IF (NORM2(cc1-cc2) < mesh_bot%tol_dist) CYCLE
      IF (NORM2(q-cc1) < mesh_bot%tol_dist .OR. NORM2(q-cc2) < mesh_bot%tol_dist .OR. lies_on_line_segment( cc1, cc2, q, mesh_bot%tol_dist)) THEN
        ! q lies on the Voronoi boundary of via
        p_next       = q
        vi_bot_left  = via
        vi_bot_right = via
        Finished     = .TRUE.
        RETURN    
      END IF
    END DO

  ! Check IF q lies on the boundary of the Voronoi cell of vib
  ! ==========================================================
  
    DO ci = 1, mesh_bot%nC( vib)
      aci = mesh_bot%iAci( vib, ci)
      cc1 = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
      cc2 = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]
      IF (NORM2(cc1-cc2) < mesh_bot%tol_dist) CYCLE
      IF (NORM2(q-cc1) < mesh_bot%tol_dist .OR. NORM2(q-cc2) < mesh_bot%tol_dist .OR. lies_on_line_segment( cc1, cc2, q, mesh_bot%tol_dist)) THEN
        ! q lies on the Voronoi boundary of via
        p_next       = q
        vi_bot_left  = vib
        vi_bot_right = vib
        Finished     = .TRUE.
        RETURN    
      END IF
    END DO

  ! Check IF q lies in the interior of either of the two adjacent Voronoi cells
  ! ===========================================================================
  
    IF     (is_in_Voronoi_cell_remap( mesh_bot, via, Vor_lines_bot, q)) THEN
      ! q lies in the interior of the Voronoi cell of via
      p_next       = q
      vi_bot_left  = via
      vi_bot_right = via
      Finished     = .TRUE.
      RETURN  
    ELSEIF (is_in_Voronoi_cell_remap( mesh_bot, vib, Vor_lines_bot, q)) THEN
      ! q lies in the interior of the Voronoi cell of vib
      p_next       = q
      vi_bot_left  = vib
      vi_bot_right = vib
      Finished     = .TRUE.
      RETURN  
    END IF

  ! Check IF q passes through either of the two endpoints of the cc-cc line
  ! =======================================================================
  
    IF     (lies_on_line_segment( p, q, ccr, mesh_bot%tol_dist) .AND. NORM2(p-ccr) > mesh_bot%tol_dist) THEN
      ! q passes through ccr
      p_next       = ccr
      vi_bot_left  = vib
      vi_bot_right = via
      tri_coincide = tir
      Coincides    = .TRUE.
      RETURN
    ELSEIF (lies_on_line_segment( p, q, ccl, mesh_bot%tol_dist) .AND. NORM2(p-ccl) > mesh_bot%tol_dist) THEN
      ! q passes through ccl
      p_next       = ccl
      vi_bot_left  = via
      vi_bot_right = vib
      tri_coincide = til
      Coincides    = .TRUE.
      RETURN
    END IF

  ! Check IF pq passes through the boundary of the Voronoi cell of via
  ! ==================================================================
  
    DO ci = 1, mesh_bot%nC( via)
      IF (mesh_bot%C( via, ci)==vib) CYCLE
      aci = mesh_bot%iAci( via, ci)
      vic = mesh_bot%C( via, ci)
      vi1 = Vor_vi_ti_bot( aci, 3)
      vi2 = Vor_vi_ti_bot( aci, 4)
      ti1 = Vor_vi_ti_bot( aci, 5)
      ti2 = Vor_vi_ti_bot( aci, 6)
      cc1 = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
      cc2 = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]

      IF     (lies_on_line_segment( p, q, cc1, mesh_bot%tol_dist)) THEN
        ! pq passes through this Voronoi vertex
        p_next       = cc1
        vi_bot_left  = via
        vi_bot_right = via
        tri_coincide = ti1
        RETURN
      ELSEIF (lies_on_line_segment( p, q, cc2, mesh_bot%tol_dist)) THEN
        ! pq passes through this Voronoi vertex
        p_next       = cc2
        vi_bot_left  = via
        vi_bot_right = via
        tri_coincide = ti2
        RETURN
      END IF

      CALL segment_intersection( cc1, cc2, p, q, llis, do_cross, mesh_bot%tol_dist)
      IF (do_cross) THEN
        ! pq cross this line
        p_next       = llis
        vi_bot_left  = via
        vi_bot_right = via
        aci_coincide = aci
        RETURN    
      END IF
    END DO

  ! Check IF pq passes through the boundary of the Voronoi cell of vib
  ! ==================================================================
  
    DO ci = 1, mesh_bot%nC( vib)
      IF (mesh_bot%C( vib, ci)==via) CYCLE
      aci = mesh_bot%iAci( vib, ci)
      vic = mesh_bot%C( vib, ci)
      vi1 = Vor_vi_ti_bot( aci, 3)
      vi2 = Vor_vi_ti_bot( aci, 4)
      ti1 = Vor_vi_ti_bot( aci, 5)
      ti2 = Vor_vi_ti_bot( aci, 6)
      cc1 = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
      cc2 = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]

      IF     (lies_on_line_segment( p, q, cc1, mesh_bot%tol_dist)) THEN
        ! pq passes through this Voronoi vertex
        p_next       = cc1
        vi_bot_left  = vib
        vi_bot_right = vib
        tri_coincide = ti1
        RETURN
      ELSEIF (lies_on_line_segment( p, q, cc2, mesh_bot%tol_dist)) THEN
        ! pq passes through this Voronoi vertex
        p_next       = cc2
        vi_bot_left  = vib
        vi_bot_right = vib
        tri_coincide = ti2
        RETURN
      END IF

      CALL segment_intersection( cc1, cc2, p, q, llis, do_cross, mesh_bot%tol_dist)
      IF (do_cross) THEN
        ! pq cross this line
        p_next       = llis
        vi_bot_left  = vib
        vi_bot_right = vib
        aci_coincide = aci
        RETURN    
      END IF
    END DO

    ! None of the available options are .TRUE. this shouldn't be possible
    WRITE(0,*) '     trace_line_through_mesh_aci_coincide - ERROR: reached the unreachable point!'
    STOP
    
  END SUBROUTINE trace_line_through_mesh_aci_coincide
  SUBROUTINE trace_line_through_mesh_tri_coincide( mesh_bot, p, q, Vor_lines_bot, Vor_vi_ti_bot, aci_coincide, tri_coincide, vi_in, p_next, vi_bot_left, vi_bot_right, Coincides, Finished)
    ! Given a point p that either:
    ! - coincides with mesh_bot Voronoi boundary line aci_coincide,
    ! - coincides with mesh_bot triangle circumcenter (= Voronoi boundary vertex) tri_coincide, or
    ! - lies in the interior of Voronoi cell vi_in,
    ! find the point p_next where the line pq crosses into a new Voronoi cell, and return indices
    ! of the mesh_bot vertices vi_bot_left and vi_bot_right whose Voronoi cells lie to the left and
    ! right of this line, respectively. If these are two different ones, pq Coincides with a Voronoi
    ! boundary line of mesh_bot. If q lies inside, or on the boundary of, the same Voronoi cell as p,
    ! we are Finished.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_bot
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: q
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_bot
    INTEGER,  DIMENSION(:,:),                INTENT(IN)    :: Vor_vi_ti_bot
    INTEGER,                                 INTENT(INOUT) :: aci_coincide
    INTEGER,                                 INTENT(INOUT) :: tri_coincide
    INTEGER,                                 INTENT(INOUT) :: vi_in
    REAL(dp), DIMENSION(2),                  INTENT(OUT)   :: p_next
    INTEGER,                                 INTENT(OUT)   :: vi_bot_left
    INTEGER,                                 INTENT(OUT)   :: vi_bot_right
    LOGICAL,                                 INTENT(OUT)   :: Coincides
    LOGICAL,                                 INTENT(OUT)   :: Finished
    
    ! Local variables:
    INTEGER                                                :: ti, vi, ci, aci, via, vib, n, til, tir, vil, vir, vvi
    REAL(dp), DIMENSION(2)                                 :: cc, ccr, ccl, llis
    LOGICAL                                                :: IsListed, do_cross

    ! Safety
    IF (aci_coincide>0 .OR.tri_coincide==0 .OR.vi_in>0) THEN
      WRITE(0,*) '    Remapping - ERROR: coincidence indicators dont match!'
      STOP
    END IF

    ! The circumcenters, and the three mesh vertices
    ti  = tri_coincide
    cc  = mesh_bot%Tricc( ti,:)

    ! Safety
    IF (NORM2(p-cc) > mesh_bot%tol_dist) THEN
      WRITE(0,*) '    Remapping - ERROR: coincidence indicators dont match!'
      STOP
    END IF

    ! Initialise output
    aci_coincide = 0
    tri_coincide = 0
    vi_in        = 0
    p_next       = p
    vi_bot_left  = 0
    vi_bot_right = 0
    Coincides    = .FALSE.
    Finished     = .FALSE.

    ! Triangle circumcenters can somrtimes overlap, making this a tricky one.
    ! First, list all vertices that have this point as a Voronoi vertex, and
    ! the indices of all Voronoi lines of which this point is an end point.
    ! Store them in in the mesh_bot vertex stacks DO memory efficiency

    mesh_bot%VStackN1   = 0
    mesh_bot%VStackN2   = 0

    DO n = 1, 3
      vi = mesh_bot%Tri( ti, n)

      DO ci = 1, mesh_bot%nC( vi)
        aci = mesh_bot%iAci( vi, ci)
        via = Vor_vi_ti_bot( aci,1)
        vib = Vor_vi_ti_bot( aci,2)
        ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
        ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]

        ! Skip zero-length lines
        IF (NORM2(ccr-ccl) < mesh_bot%tol_dist) CYCLE

        IF (NORM2(p-ccr) < mesh_bot%tol_dist .OR. NORM2(p-ccl) < mesh_bot%tol_dist) THEN
          ! p overlaps with an endpoint of this Voronoi boundary line
          ! Store the indices of the line and the vertices it separates in the stacks

          ! Store aci in stack 1
          IsListed = .FALSE.
          DO vvi = 1, mesh_bot%VStackN1
            IF (mesh_bot%VStack1(vvi) == aci) THEN
              IsListed = .TRUE.
              EXIT
            END IF
          END DO
          IF (.NOT. IsListed) THEN
            mesh_bot%VStackN1 = mesh_bot%VStackN1 + 1
            mesh_bot%VStack1( mesh_bot%VStackN1) = aci
          END IF

          ! Store via and vib in stack 2
          IsListed = .FALSE.
          DO vvi = 1, mesh_bot%VStackN2
            IF (mesh_bot%VStack2(vvi) == via) THEN
              IsListed = .TRUE.
              EXIT
            END IF
          END DO
          IF (.NOT. IsListed) THEN
            mesh_bot%VStackN2 = mesh_bot%VStackN2 + 1
            mesh_bot%VStack2( mesh_bot%VStackN2) = via
          END IF

          IsListed = .FALSE.
          DO vvi = 1, mesh_bot%VStackN2
            IF (mesh_bot%VStack2(vvi) == vib) THEN
              IsListed = .TRUE.
              EXIT
            END IF
          END DO
          IF (.NOT. IsListed) THEN
            mesh_bot%VStackN2 = mesh_bot%VStackN2 + 1
            mesh_bot%VStack2( mesh_bot%VStackN2) = vib
          END IF
        END IF ! IF (NORM2(p-ccr) < mesh_bot%tol_dist .OR. NORM2(p-ccl) < mesh_bot%tol_dist) THEN
      END DO ! DO ci = 1, mesh_bot%nC( vi)
    END DO ! DO n = 1, 3

    ! Check IF q lies on the endpoints of any of the Voronoi lines originating here
    DO ci = 1, mesh_bot%VStackN1

      aci = mesh_bot%VStack1( ci)
      via = Vor_vi_ti_bot( aci,1)
      vib = Vor_vi_ti_bot( aci,2)
      vir = Vor_vi_ti_bot( aci,3)
      vil = Vor_vi_ti_bot( aci,4)
      tir = Vor_vi_ti_bot( aci,5)
      til = Vor_vi_ti_bot( aci,6)
      ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
      ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]

      ! Check IF q lies on the endpoint of this Voronoi line
      IF     (NORM2(q-ccr) < mesh_bot%tol_dist) THEN
        ! q coincides with ccr
        p_next       = q
        vi_bot_left  = vib
        vi_bot_right = via
        Coincides    = .TRUE.
        Finished     = .TRUE.
        RETURN
      ELSEIF (NORM2(q-ccl) < mesh_bot%tol_dist) THEN
        ! q coincides with ccl
        p_next       = q
        vi_bot_left  = via
        vi_bot_right = vib
        Coincides    = .TRUE.
        Finished     = .TRUE.
        RETURN
      END IF

    END DO ! DO ci = 1, mesh_bot%VStackN1

    ! Check IF q lies on any of the Voronoi lines originating here
    DO ci = 1, mesh_bot%VStackN1

      aci = mesh_bot%VStack1( ci)
      via = Vor_vi_ti_bot( aci,1)
      vib = Vor_vi_ti_bot( aci,2)
      vir = Vor_vi_ti_bot( aci,3)
      vil = Vor_vi_ti_bot( aci,4)
      tir = Vor_vi_ti_bot( aci,5)
      til = Vor_vi_ti_bot( aci,6)
      ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
      ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]

      ! Check IF q lies on this Voronoi line
      IF     (lies_on_line_segment( p, ccr, q, mesh_bot%tol_dist) .AND. NORM2(p-ccr) > mesh_bot%tol_dist) THEN
        ! q coincides with ccr
        p_next       = q
        vi_bot_left  = vib
        vi_bot_right = via
        Coincides    = .TRUE.
        Finished     = .TRUE.
        RETURN
      ELSEIF (lies_on_line_segment( p, ccl, q, mesh_bot%tol_dist) .AND. NORM2(p-ccl) > mesh_bot%tol_dist) THEN
        ! q coincides with ccl
        p_next       = q
        vi_bot_left  = via
        vi_bot_right = vib
        Coincides    = .TRUE.
        Finished     = .TRUE.
        RETURN
      END IF

    END DO ! DO ci = 1, mesh_bot%VStackN1

    ! Check IF pq passes through the endpoints of any of the Voronoi lines originating here
    DO ci = 1, mesh_bot%VStackN1

      aci = mesh_bot%VStack1( ci)
      via = Vor_vi_ti_bot( aci,1)
      vib = Vor_vi_ti_bot( aci,2)
      vir = Vor_vi_ti_bot( aci,3)
      vil = Vor_vi_ti_bot( aci,4)
      tir = Vor_vi_ti_bot( aci,5)
      til = Vor_vi_ti_bot( aci,6)
      ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
      ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]

      IF (NORM2(ccl-ccr) < mesh_bot%tol_dist) CYCLE

      ! Check IF pq passes through the endpoint of this Voronoi line
      IF     (NORM2(p-ccl) < mesh_bot%tol_dist .AND. lies_on_line_segment( p, q, ccr, mesh_bot%tol_dist)) THEN
        ! q coincides with ccr
        p_next       = ccr
        vi_bot_left  = vib
        vi_bot_right = via
        Coincides    = .TRUE.
        tri_coincide = tir
        RETURN
      ELSEIF (NORM2(p-ccr) < mesh_bot%tol_dist .AND. lies_on_line_segment( p, q, ccl, mesh_bot%tol_dist)) THEN
        ! q coincides with ccl
        p_next       = ccl
        vi_bot_left  = via
        vi_bot_right = vib
        Coincides    = .TRUE.
        tri_coincide = til
        RETURN
      END IF

    END DO ! DO ci = 1, mesh_bot%VStackN1

    ! Check IF q lies on the boundary or in the interior of any of the surrounding Voronoi cells
    DO vvi = 1, mesh_bot%VStackN2

      vi = mesh_bot%VStack2( vvi)

      ! Check IF q lies on the boundary of this Voronoi cell
      DO ci = 1, mesh_bot%nC( vi)
        aci = mesh_bot%iAci( vi, ci)
        ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
        ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]

        IF (NORM2(ccl-ccr) < mesh_bot%tol_dist) CYCLE

        IF (NORM2(q-ccr) < mesh_bot%tol_dist .OR. NORM2(q-ccl) < mesh_bot%tol_dist .OR. lies_on_line_segment( ccl, ccr, q, mesh_bot%tol_dist)) THEN
          ! q lies on the boundary of this Voronoi cell
          p_next       = q
          vi_bot_left  = vi
          vi_bot_right = vi
          Finished     = .TRUE.
          RETURN      
        END IF
      END DO

      ! Check IF q lies in the interior of this Voronoi cell
      IF (is_in_Voronoi_cell_remap( mesh_bot, vi, Vor_lines_bot, q)) THEN
        ! q lies on the boundary of this Voronoi cell
        p_next       = q
        vi_bot_left  = vi
        vi_bot_right = vi
        Finished     = .TRUE.
        RETURN 
      END IF
    END DO ! DO vvi = 1, mesh_bot%VStackN2

    ! Check IF pq crosses through any vertices of the boundaries of the surrounding Voronoi cells
    DO vvi = 1, mesh_bot%VStackN2

      vi = mesh_bot%VStack2( vvi)

      ! Check IF pq crosses through a vertx of this Voronoi cell
      DO ci = 1, mesh_bot%nC( vi)
        aci = mesh_bot%iAci( vi, ci)
        tir = Vor_vi_ti_bot( aci,5)
        til = Vor_vi_ti_bot( aci,6)
        ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
        ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]
        
        IF (NORM2(ccr-ccl) < mesh_bot%tol_dist) CYCLE

        IF     (lies_on_line_segment( p, q, ccr, mesh_bot%tol_dist) .AND. NORM2(p-ccr) > mesh_bot%tol_dist) THEN
          ! pq passes through this Voronoi vertex
          p_next       = ccr
          vi_bot_left  = vi
          vi_bot_right = vi
          tri_coincide = tir
          RETURN
        ELSEIF (lies_on_line_segment( p, q, ccl, mesh_bot%tol_dist) .AND. NORM2(p-ccl) > mesh_bot%tol_dist) THEN
          ! pq passes through this Voronoi vertex
          p_next       = ccl
          vi_bot_left  = vi
          vi_bot_right = vi
          tri_coincide = til
          RETURN
        END IF
      END DO
    END DO ! DO vvi = 1, mesh_bot%VStackN2

    ! Check IF pq crosses the boundary of any of the surrounding Voronoi cells
    DO vvi = 1, mesh_bot%VStackN2

      vi = mesh_bot%VStack2( vvi)

      ! Check IF pq crosses the boundary of this Voronoi cell
      DO ci = 1, mesh_bot%nC( vi)
        aci = mesh_bot%iAci( vi, ci)
        ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
        ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]

        IF (NORM2(ccr-ccl) < mesh_bot%tol_dist) CYCLE

        CALL segment_intersection( p, q, ccl, ccr, llis, do_cross, mesh_bot%tol_dist)
        IF (do_cross .AND. NORM2(p-ccr) > mesh_bot%tol_dist .AND. NORM2(p-ccl) > mesh_bot%tol_dist) THEN
          ! pq crosses the boundary of this Voronoi cell
          p_next       = llis
          vi_bot_left  = vi
          vi_bot_right = vi
          aci_coincide = aci
          RETURN
        END IF
      END DO
    END DO ! DO vvi = 1, mesh_bot%VStackN2

    ! None of the available options are .TRUE. this shouldn't be possible
    WRITE(0,*) '     trace_line_through_mesh_tri_coincide - ERROR: reached the unreachable point!'
    STOP
    
  END SUBROUTINE trace_line_through_mesh_tri_coincide
  SUBROUTINE trace_line_through_mesh_vi_interior(  mesh_bot, p, q, Vor_lines_bot, Vor_vi_ti_bot, aci_coincide, tri_coincide, vi_in, p_next, vi_bot_left, vi_bot_right, Coincides, Finished)
    ! Given a point p that either:
    ! - coincides with mesh_bot Voronoi boundary line aci_coincide,
    ! - coincides with mesh_bot triangle circumcenter (= Voronoi boundary vertex) tri_coincide, or
    ! - lies in the interior of Voronoi cell vi_in,
    ! find the point p_next where the line pq crosses into a new Voronoi cell, and return indices
    ! of the mesh_bot vertices vi_bot_left and vi_bot_right whose Voronoi cells lie to the left and
    ! right of this line, respectively. If these are two different ones, pq Coincides with a Voronoi
    ! boundary line of mesh_bot. If q lies inside, or on the boundary of, the same Voronoi cell as p,
    ! we are Finished.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_bot
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: q
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_bot
    INTEGER,  DIMENSION(:,:),                INTENT(IN)    :: Vor_vi_ti_bot
    INTEGER,                                 INTENT(INOUT) :: aci_coincide
    INTEGER,                                 INTENT(INOUT) :: tri_coincide
    INTEGER,                                 INTENT(INOUT) :: vi_in
    REAL(dp), DIMENSION(2),                  INTENT(OUT)   :: p_next
    INTEGER,                                 INTENT(OUT)   :: vi_bot_left
    INTEGER,                                 INTENT(OUT)   :: vi_bot_right
    LOGICAL,                                 INTENT(OUT)   :: Coincides
    LOGICAL,                                 INTENT(OUT)   :: Finished
    
    ! Local variables:
    INTEGER                                                :: via, ci, aci, vib, vir, vil, tir, til
    REAL(dp), DIMENSION(2)                                 :: ccr, ccl, llis
    LOGICAL                                                :: do_cross

    ! Safety
    IF (aci_coincide>0 .OR.tri_coincide>0 .OR.vi_in==0) THEN
      WRITE(0,*) '    Remapping - ERROR: coincidence indicators dont match!'
      STOP
    END IF

    ! The mesh_bot vertex whose Voronoi cell contains p
    via = vi_in

    ! Safety
    via = vi_in
    CALL find_containing_vertex( mesh_bot, p, via)
    IF (vi_in /= via) THEN
      WRITE(0,*) '    Remapping - ERROR: coincidence indicators dont match!'
      STOP
    END IF

    ! Initialise output
    aci_coincide = 0
    tri_coincide = 0
    vi_in        = 0
    vi_bot_left  = via
    vi_bot_right = via
    Coincides    = .FALSE.
    Finished     = .FALSE.

    ! Check IF q lies on any of the vertices or lines of the boundary of the Voronoi cell
    DO ci = 1, mesh_bot%nC( via)
      aci = mesh_bot%iAci( via, ci)
      ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
      ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]
      IF (NORM2(ccr-ccl) < mesh_bot%tol_dist) CYCLE
      IF (NORM2(q-ccl) < mesh_bot%tol_dist .OR. NORM2(q-ccr) < mesh_bot%tol_dist .OR. lies_on_line_segment( ccr, ccl, q, mesh_bot%tol_dist)) THEN
        ! q lies on this circumcenter
        p_next       = q
        Finished     = .TRUE.
        RETURN
      END IF
    END DO

    ! Check IF q lies in the interior of the Voronoi cell
    IF (is_in_Voronoi_cell_remap( mesh_bot, via, Vor_lines_bot, q)) THEN
      ! q lies in the interior of the Voronoi cell of via
      p_next       = q
      Finished     = .TRUE.
      RETURN  
    END IF

    ! Check IF pq passes through any of the vertices, or crosses any of the lines, of the boundary of the Voronoi cell
    DO ci = 1, mesh_bot%nC( via)
      aci = mesh_bot%iAci( via, ci)
      ccr = [Vor_lines_bot( aci,1,1), Vor_lines_bot( aci,1,2)]
      ccl = [Vor_lines_bot( aci,2,1), Vor_lines_bot( aci,2,2)]
      vib = mesh_bot%C( via, ci)
      vir = Vor_vi_ti_bot( aci,3)
      vil = Vor_vi_ti_bot( aci,4)
      tir = Vor_vi_ti_bot( aci,5)
      til = Vor_vi_ti_bot( aci,6)

      IF (NORM2(ccr-ccl) < mesh_bot%tol_dist) CYCLE

      IF     (lies_on_line_segment(p, q, ccr, mesh_bot%tol_dist)) THEN
        ! q passes through this circumcenter
        p_next       = ccr
        tri_coincide = tir
        RETURN
      ELSEIF (lies_on_line_segment(p, q, ccl, mesh_bot%tol_dist)) THEN
        ! q passes through this circumcenter
        p_next       = ccl
        tri_coincide = til
        RETURN
      END IF

      CALL segment_intersection( ccr, ccl, p, q, llis, do_cross, mesh_bot%tol_dist)
      IF (do_cross) THEN
        ! pq cross this line
        p_next       = llis
        aci_coincide = aci
        RETURN
      END IF
    END DO

    ! None of the available options are .TRUE. this shouldn't be possible
    WRITE(0,*) '     trace_line_through_mesh_vi_interior - ERROR: reached the unreachable point!'
    STOP
    
  END SUBROUTINE trace_line_through_mesh_vi_interior
    
  FUNCTION is_in_Voronoi_cell_remap( mesh, vi, Vor_lines, q) RESULT(isso)
    ! Use a series of cross products with the Voronoi cell boundary lines
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    INTEGER,                                 INTENT(IN)    :: vi
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: q
    LOGICAL                                                :: isso
    
    ! Local variables:
    INTEGER                                                :: ci, aci
    REAL(dp), DIMENSION(2)                                 :: p1, p2, dl, dq
    
    isso = .TRUE.
    DO ci = 1, mesh%nC( vi)
      aci = mesh%iAci( vi, ci)
      IF (mesh%Aci( aci,1) == vi) THEN
        p1 = [Vor_lines( aci,1,1), Vor_lines( aci,1,2)]
        p2 = [Vor_lines( aci,2,1), Vor_lines( aci,2,2)]
      ELSE
        p2 = [Vor_lines( aci,1,1), Vor_lines( aci,1,2)]
        p1 = [Vor_lines( aci,2,1), Vor_lines( aci,2,2)]
      END IF
      IF (NORM2(p2-p1) < mesh%tol_dist) CYCLE
      dl = p2 - p1
      dq = q  - p1
      IF (cross2( dl, dq) < 0._dp) THEN
        isso = .FALSE.
        RETURN
      END IF
    END DO

  END FUNCTION is_in_Voronoi_cell_remap
  
! == Low-level remapping routines
  SUBROUTINE remap_trilin_2D(                           mesh_dst, map, d_src, d_dst    )
    ! Remap data from mesh_src to mesh_dst using trilinear interpolation
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst  ! Destination mesh
    TYPE(type_remapping_trilin),             INTENT(IN)    :: map       ! Remapping array
    REAL(dp), DIMENSION(:),                  INTENT(IN)    :: d_src     ! Data field on source      mesh
    REAL(dp), DIMENSION(:),                  INTENT(OUT)   :: d_dst     ! Data field on destination mesh
    
    ! Local variables
    INTEGER                                                :: vi_dst
        
    d_dst(mesh_dst%v1:mesh_dst%v2) = 0._dp
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      d_dst(vi_dst) = (d_src( map%vi( vi_dst,1)) * map%w( vi_dst,1)) + &
                      (d_src( map%vi( vi_dst,2)) * map%w( vi_dst,2)) + &
                      (d_src( map%vi( vi_dst,3)) * map%w( vi_dst,3))
    END DO
    CALL sync
    
  END SUBROUTINE remap_trilin_2D
  SUBROUTINE remap_trilin_3D(                           mesh_dst, map, d_src, d_dst, nz)
    ! Remap data from mesh_src to mesh_dst using trilinear interpolation
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst  ! Destination mesh
    TYPE(type_remapping_trilin),             INTENT(IN)    :: map       ! Remapping array
    REAL(dp), DIMENSION(:,:),                INTENT(IN)    :: d_src     ! Data field on source      mesh
    REAL(dp), DIMENSION(:,:),                INTENT(OUT)   :: d_dst     ! Data field on destination mesh
    INTEGER,                                 INTENT(IN)    :: nz
    
    ! Local variables
    INTEGER                                                :: vi_dst, k
        
    d_dst(mesh_dst%v1:mesh_dst%v2,:) = 0._dp
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      DO k = 1, nz
        d_dst(vi_dst,k) = (d_src( map%vi( vi_dst,1),k) * map%w( vi_dst,1)) + &
                          (d_src( map%vi( vi_dst,2),k) * map%w( vi_dst,2)) + &
                          (d_src( map%vi( vi_dst,3),k) * map%w( vi_dst,3))
      END DO
    END DO
    CALL sync
    
  END SUBROUTINE remap_trilin_3D
  SUBROUTINE remap_trilin_monthly(                      mesh_dst, map, d_src, d_dst    )
    ! Remap data from mesh_src to mesh_dst using trilinear interpolation
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst  ! Destination mesh
    TYPE(type_remapping_trilin),             INTENT(IN)    :: map       ! Remapping array
    REAL(dp), DIMENSION(:,:),                INTENT(IN)    :: d_src     ! Data field on source      mesh
    REAL(dp), DIMENSION(:,:),                INTENT(OUT)   :: d_dst     ! Data field on destination mesh
    
    ! Local variables
    INTEGER                                                :: vi_dst, m
        
    d_dst(mesh_dst%v1:mesh_dst%v2,:) = 0._dp
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      DO m = 1, 12
        d_dst(vi_dst,m) = (d_src( map%vi( vi_dst,1),m) * map%w( vi_dst,1)) + &
                          (d_src( map%vi( vi_dst,2),m) * map%w( vi_dst,2)) + &
                          (d_src( map%vi( vi_dst,3),m) * map%w( vi_dst,3))
      END DO
    END DO
    CALL sync
    
  END SUBROUTINE remap_trilin_monthly
  SUBROUTINE remap_nearest_neighbour_2D(                mesh_dst, map, d_src, d_dst    )
    ! Remap data from mesh_src to mesh_dst using nearest neighbour interpolation
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst  ! Destination mesh
    TYPE(type_remapping_nearest_neighbour),  INTENT(IN)    :: map       ! Remapping array
    REAL(dp), DIMENSION(:),                  INTENT(IN)    :: d_src     ! Data field on source      mesh
    REAL(dp), DIMENSION(:),                  INTENT(OUT)   :: d_dst     ! Data field on destination mesh
    
    ! Local variables
    INTEGER                                                :: vi_dst
        
    d_dst(mesh_dst%v1:mesh_dst%v2) = 0._dp
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      d_dst(vi_dst) = d_src( map%vi( vi_dst))
    END DO
    CALL sync
    
  END SUBROUTINE remap_nearest_neighbour_2D
  SUBROUTINE remap_nearest_neighbour_3D(                mesh_dst, map, d_src, d_dst, nz)
    ! Remap data from mesh_src to mesh_dst using nearest neighbour interpolation
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst  ! Destination mesh
    TYPE(type_remapping_nearest_neighbour),  INTENT(IN)    :: map       ! Remapping array
    REAL(dp), DIMENSION(:,:),                INTENT(IN)    :: d_src     ! Data field on source      mesh
    REAL(dp), DIMENSION(:,:),                INTENT(OUT)   :: d_dst     ! Data field on destination mesh
    INTEGER,                                 INTENT(IN)    :: nz
    
    ! Local variables
    INTEGER                                                :: vi_dst, k
        
    d_dst(mesh_dst%v1:mesh_dst%v2,:) = 0._dp
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      DO k = 1, nz
        d_dst(vi_dst,k) = d_src( map%vi( vi_dst),k)
      END DO
    END DO
    CALL sync
    
  END SUBROUTINE remap_nearest_neighbour_3D
  SUBROUTINE remap_cons_1st_order_2D(                   mesh_dst, map, d_src, d_dst    )
    ! Remap data from mesh_src to mesh_dst using 1st order conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst  ! Destination mesh
    TYPE(type_remapping_conservative),       INTENT(IN)    :: map       ! Remapping array
    REAL(dp), DIMENSION(:),                  INTENT(IN)    :: d_src     ! Data field on source      mesh
    REAL(dp), DIMENSION(:),                  INTENT(OUT)   :: d_dst     ! Data field on destination mesh
    
    ! Local variables
    INTEGER                                                :: vi_dst, nvi
    
    d_dst(mesh_dst%v1:mesh_dst%v2) = 0._dp        
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      DO nvi = 1, map%nV(vi_dst)
        d_dst(vi_dst) = d_dst(vi_dst) + (d_src(map%vi(vi_dst,nvi)) * map%w0( vi_dst,nvi))
      END DO
    END DO
    CALL sync
    
  END SUBROUTINE remap_cons_1st_order_2D
  SUBROUTINE remap_cons_1st_order_3D(                   mesh_dst, map, d_src, d_dst, nz)
    ! Remap data from mesh_src to mesh_dst using 1st order conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst  ! Destination mesh
    TYPE(type_remapping_conservative),       INTENT(IN)    :: map       ! Remapping array
    REAL(dp), DIMENSION(:,:),                INTENT(IN)    :: d_src     ! Data field on source      mesh
    REAL(dp), DIMENSION(:,:),                INTENT(OUT)   :: d_dst     ! Data field on destination mesh
    INTEGER,                                 INTENT(IN)    :: nz
    
    ! Local variables
    INTEGER                                                :: vi_dst, nvi, k
    
    d_dst(mesh_dst%v1:mesh_dst%v2,:) = 0._dp        
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      DO nvi = 1, map%nV(vi_dst)
        DO k = 1, nz
          d_dst(vi_dst,k) = d_dst(vi_dst,k) + (d_src(map%vi(vi_dst,nvi),k) * map%w0( vi_dst,nvi))
        END DO
      END DO
    END DO
    CALL sync
    
  END SUBROUTINE remap_cons_1st_order_3D
  SUBROUTINE remap_cons_2nd_order_2D(         mesh_src, mesh_dst, map, d_src, d_dst    )
    ! Remap data from mesh_src to mesh_dst using 2nd order conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_src  ! Source      mesh
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst  ! Destination mesh
    TYPE(type_remapping_conservative),       INTENT(IN)    :: map       ! Remapping array
    REAL(dp), DIMENSION(:),                  INTENT(IN)    :: d_src     ! Data field on source      mesh
    REAL(dp), DIMENSION(:),                  INTENT(OUT)   :: d_dst     ! Data field on destination mesh
    
    ! Local variables
    INTEGER                                                :: vi_dst, nvi
    REAL(dp), DIMENSION(:), POINTER                        :: ddx_src, ddy_src
    INTEGER                                                :: wddx_src, wddy_src
    
    CALL allocate_shared_dp_1D(  mesh_src%nV, ddx_src, wddx_src)
    CALL allocate_shared_dp_1D(  mesh_src%nV, ddy_src, wddy_src)
    
    d_dst(mesh_dst%v1:mesh_dst%v2) = 0._dp
    CALL sync
    
    CALL get_mesh_derivatives( mesh_src, d_src, ddx_src, ddy_src)
    CALL sync
        
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      DO nvi = 1, map%nV(vi_dst)
        d_dst(vi_dst) = d_dst(vi_dst) + (d_src(   map%vi( vi_dst,nvi)) * map%w0(  vi_dst,nvi)) + &
                                        (ddx_src( map%vi( vi_dst,nvi)) * map%w1x( vi_dst,nvi)) + &
                                        (ddy_src( map%vi( vi_dst,nvi)) * map%w1y( vi_dst,nvi))
      END DO
    END DO
    CALL sync
    
    CALL deallocate_shared( wddx_src)
    CALL deallocate_shared( wddy_src)
    
    NULLIFY( ddx_src)
    NULLIFY( ddy_src)
    
  END SUBROUTINE remap_cons_2nd_order_2D
  SUBROUTINE remap_cons_2nd_order_2D_monthly( mesh_src, mesh_dst, map, d_src, d_dst    )
    ! Remap data from mesh_src to mesh_dst using 2nd order conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_src  ! Source      mesh
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst  ! Destination mesh
    TYPE(type_remapping_conservative),       INTENT(IN)    :: map       ! Remapping array
    REAL(dp), DIMENSION(:,:),                INTENT(IN)    :: d_src     ! Data field on source      mesh
    REAL(dp), DIMENSION(:,:),                INTENT(OUT)   :: d_dst     ! Data field on destination mesh
    
    ! Local variables
    INTEGER                                                :: vi_dst, nvi, m
    REAL(dp), DIMENSION(:,:), POINTER                      :: ddx_src, ddy_src
    INTEGER                                                :: wddx_src, wddy_src
    
    CALL allocate_shared_dp_2D(  mesh_src%nV, 12, ddx_src, wddx_src)
    CALL allocate_shared_dp_2D(  mesh_src%nV, 12, ddy_src, wddy_src)
    
    d_dst(mesh_dst%v1:mesh_dst%v2,:) = 0._dp
    
    CALL get_mesh_derivatives_3D( mesh_src, d_src, ddx_src, ddy_src, 12)
        
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      DO nvi = 1, map%nV(vi_dst)
        DO m = 1, 12
          d_dst(vi_dst,m) = d_dst(vi_dst,m) + (d_src(   map%vi( vi_dst,nvi),m) * map%w0(  vi_dst,nvi)) + &
                                              (ddx_src( map%vi( vi_dst,nvi),m) * map%w1x( vi_dst,nvi)) + &
                                              (ddy_src( map%vi( vi_dst,nvi),m) * map%w1y( vi_dst,nvi))
        END DO
      END DO
    END DO
    CALL sync
    
    CALL deallocate_shared( wddx_src)
    CALL deallocate_shared( wddy_src)
    
    NULLIFY( ddx_src)
    NULLIFY( ddy_src)
    
  END SUBROUTINE remap_cons_2nd_order_2D_monthly
  SUBROUTINE remap_cons_2nd_order_3D(         mesh_src, mesh_dst, map, d_src, d_dst, nz)
    ! Remap data from mesh_src to mesh_dst using 2nd order conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_src  ! Source      mesh
    TYPE(type_mesh),                         INTENT(IN)    :: mesh_dst  ! Destination mesh
    TYPE(type_remapping_conservative),       INTENT(IN)    :: map       ! Remapping array
    REAL(dp), DIMENSION(:,:),                INTENT(IN)    :: d_src     ! Data field on source      mesh
    REAL(dp), DIMENSION(:,:),                INTENT(OUT)   :: d_dst     ! Data field on destination mesh
    INTEGER,                                 INTENT(IN)    :: nz
    
    ! Local variables
    INTEGER                                                :: vi_dst, nvi, k
    REAL(dp), DIMENSION(:,:), POINTER                      :: ddx_src, ddy_src
    INTEGER                                                :: wddx_src, wddy_src
    
    CALL allocate_shared_dp_2D(  mesh_src%nV, C%nz, ddx_src, wddx_src)
    CALL allocate_shared_dp_2D(  mesh_src%nV, C%nz, ddy_src, wddy_src)
    
    d_dst(mesh_dst%v1:mesh_dst%v2,:) = 0._dp
    
    CALL get_mesh_derivatives_3D( mesh_src, d_src, ddx_src, ddy_src, nz)
        
    DO vi_dst = mesh_dst%v1, mesh_dst%v2
      DO nvi = 1, map%nV(vi_dst)
        DO k = 1, nz
          d_dst(vi_dst,k) = d_dst(vi_dst,k) + (d_src(   map%vi( vi_dst,nvi),k) * map%w0(  vi_dst,nvi)) + &
                                              (ddx_src( map%vi( vi_dst,nvi),k) * map%w1x( vi_dst,nvi)) + &
                                              (ddy_src( map%vi( vi_dst,nvi),k) * map%w1y( vi_dst,nvi))
        END DO
      END DO
    END DO
    CALL sync
    
    CALL deallocate_shared( wddx_src)
    CALL deallocate_shared( wddy_src)
    
    NULLIFY( ddx_src)
    NULLIFY( ddy_src)
    
  END SUBROUTINE remap_cons_2nd_order_3D
  
! == Medium-level remapping routines, where the memory containing the data is reallocated
  SUBROUTINE remap_field_dp(         mesh_old, mesh_new, map, d, w, method)
    ! Map a single data fields from the old mesh to the new mesh. Includes memory reallocation.
    
    IMPLICIT NONE
        
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    REAL(dp), DIMENSION(:  ), POINTER,   INTENT(INOUT) :: d            ! Pointer to the data
    INTEGER,                             INTENT(INOUT) :: w            ! MPI window to the shared memory space containing that data
    CHARACTER(LEN=*),                    INTENT(IN)    :: method
    
    ! Local variables
    REAL(dp), DIMENSION(:  ), POINTER                  :: d_temp
    INTEGER                                            :: w_temp
    INTEGER                                            :: v1_old, v2_old
    
    ! Calculate old and new process vertex domains    
    v1_old = MAX(1,             FLOOR(REAL(mesh_old%nV *  par%i      / par%n)) + 1)
    v2_old = MIN(mesh_old%nV,   FLOOR(REAL(mesh_old%nV * (par%i + 1) / par%n)))
    
    ! Allocate shared memory to temporarily store the old data
    CALL allocate_shared_dp_1D( mesh_old%nV, d_temp, w_temp)
    
    ! Copy the old data there
    d_temp(v1_old:v2_old) = d(v1_old:v2_old)
    
    ! Deallocate and reallocate the shared memory space
    CALL deallocate_shared(w)
    NULLIFY(d)    
    CALL allocate_shared_dp_1D( mesh_new%nV, d, w)
    
    ! Map data field from old mesh to new mesh
    IF (method == 'trilin') THEN
      CALL remap_trilin_2D(                      mesh_new, map%trilin,            d_temp, d)
    ELSEIF (method == 'nearest_neighbour') THEN
      CALL remap_nearest_neighbour_2D(           mesh_new, map%nearest_neighbour, d_temp, d)
    ELSEIF (method == 'cons_1st_order') THEN
      CALL remap_cons_1st_order_2D(              mesh_new, map%conservative,      d_temp, d)
    ELSEIF (method == 'cons_2nd_order') THEN
      CALL remap_cons_2nd_order_2D(    mesh_old, mesh_new, map%conservative,      d_temp, d)
    ELSE
      WRITE(0,*) ' remap_field_dp - ERROR: "method" can only be "trilin", "nearest_neighbour", "cons_1st_order" or "cons_2nd_order"!'
      STOP
    END IF
    
    ! Deallocate temporary shared memory
    CALL deallocate_shared(w_temp)
    NULLIFY(d_temp)
    
    
  END SUBROUTINE remap_field_dp
  SUBROUTINE remap_field_dp_3D(      mesh_old, mesh_new, map, d, w, method)
    ! Map a single data fields from the old mesh to the new mesh. Includes memory reallocation.
    
    IMPLICIT NONE
        
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:), POINTER,   INTENT(INOUT) :: d            ! Pointer to the data
    INTEGER,                             INTENT(INOUT) :: w            ! MPI window to the shared memory space containing that data
    CHARACTER(LEN=*),                    INTENT(IN)    :: method
    
    ! Local variables
    REAL(dp), DIMENSION(:,:), POINTER                  :: d_temp
    INTEGER                                            :: w_temp
    INTEGER                                            :: v1_old, v2_old
    
    ! Calculate old and new process vertex domains    
    v1_old = MAX(1,             FLOOR(REAL(mesh_old%nV *  par%i      / par%n)) + 1)
    v2_old = MIN(mesh_old%nV,   FLOOR(REAL(mesh_old%nV * (par%i + 1) / par%n)))
    
    ! Allocate shared memory to temporarily store the old data
    CALL allocate_shared_dp_2D( mesh_old%nV, C%nz, d_temp, w_temp)
    
    ! Copy the old data there
    d_temp(v1_old:v2_old,:) = d(v1_old:v2_old,:)
    
    ! Deallocate and reallocate the shared memory space
    CALL deallocate_shared(w)
    NULLIFY(d)    
    CALL allocate_shared_dp_2D( mesh_new%nV, C%nz, d, w)
    
    ! Map data field from old mesh to new mesh
    IF (method == 'trilin') THEN
      CALL remap_trilin_3D(                      mesh_new, map%trilin,       d_temp, d, C%nZ)    
    ELSEIF (method == 'cons_1st_order') THEN
      CALL remap_cons_1st_order_3D(              mesh_new, map%conservative, d_temp, d, C%nZ)
    ELSEIF (method == 'cons_2nd_order') THEN
      CALL remap_cons_2nd_order_3D(    mesh_old, mesh_new, map%conservative, d_temp, d, C%nZ)
    ELSE
      WRITE(0,*) ' remap_field_dp_3D - ERROR: "method" can only be "trilin", "nearest_neighbour", "cons_1st_order" or "cons_2nd_order"!'
      STOP
    END IF
    
    ! Deallocate temporary shared memory
    CALL deallocate_shared(w_temp)
    NULLIFY(d_temp)
    
    
  END SUBROUTINE remap_field_dp_3D
  SUBROUTINE remap_field_dp_monthly( mesh_old, mesh_new, map, d, w, method)
    ! Map a single data fields from the old mesh to the new mesh. Includes memory reallocation.
    
    IMPLICIT NONE
        
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:), POINTER,   INTENT(INOUT) :: d            ! Pointer to the data
    INTEGER,                             INTENT(INOUT) :: w            ! MPI window to the shared memory space containing that data
    CHARACTER(LEN=*),                    INTENT(IN)    :: method
    
    ! Local variables
    REAL(dp), DIMENSION(:,:), POINTER                  :: d_temp
    INTEGER                                            :: w_temp
    INTEGER                                            :: v1_old, v2_old
    
    ! Calculate old and new process vertex domains    
    v1_old = MAX(1,             FLOOR(REAL(mesh_old%nV *  par%i      / par%n)) + 1)
    v2_old = MIN(mesh_old%nV,   FLOOR(REAL(mesh_old%nV * (par%i + 1) / par%n)))
    
    ! Allocate shared memory to temporarily store the old data
    CALL allocate_shared_dp_2D( mesh_old%nV, 12, d_temp, w_temp)
    
    ! Copy the old data there
    d_temp(v1_old:v2_old,:) = d(v1_old:v2_old,:)
    
    ! Deallocate and reallocate the shared memory space
    CALL deallocate_shared(w)
    NULLIFY(d)    
    CALL allocate_shared_dp_2D( mesh_new%nV, 12, d, w)
    
    ! Map data field from old mesh to new mesh
    IF (method == 'trilin') THEN
      CALL remap_trilin_monthly(                      mesh_new, map%trilin,            d_temp, d)
    ELSE
      WRITE(0,*) ' remap_field_dp_monthly - ERROR: "method" can only be "trilin", "nearest_neighbour", "cons_1st_order" or "cons_2nd_order"!'
      STOP
    END IF
    
    ! Deallocate temporary shared memory
    CALL deallocate_shared(w_temp)
    NULLIFY(d_temp)
    
    
  END SUBROUTINE remap_field_dp_monthly
  
! == Simple reallocation for when actual remapping is not required
  SUBROUTINE reallocate_field_dp(    nV, d, w)
    ! For data that doesn't need to be mapped because it will be recalculated anyway.
    
    IMPLICIT NONE
    
    ! In/output variables
    INTEGER,                             INTENT(IN)    :: nV
    REAL(dp), DIMENSION(:  ), POINTER,   INTENT(INOUT) :: d            ! Pointer to the data
    INTEGER,                             INTENT(INOUT) :: w            ! MPI window to the shared memory space containing that data
    
    ! Deallocate and reallocate the shared memory space
    CALL deallocate_shared(w)
    NULLIFY(d)    
    CALL allocate_shared_dp_1D( nV, d, w)
    
  END SUBROUTINE reallocate_field_dp
  SUBROUTINE reallocate_field_dp_3D( nV, d, w, nz)
    ! For data that doesn't need to be mapped because it will be recalculated anyway.
    
    IMPLICIT NONE
    
    ! In/output variables
    INTEGER,                             INTENT(IN)    :: nV
    REAL(dp), DIMENSION(:,:), POINTER,   INTENT(INOUT) :: d            ! Pointer to the data
    INTEGER,                             INTENT(INOUT) :: w            ! MPI window to the shared memory space containing that data
    INTEGER,                             INTENT(IN)    :: nz
    
    ! Deallocate and reallocate the shared memory space
    CALL deallocate_shared(w)
    NULLIFY(d)    
    CALL allocate_shared_dp_2D( nV, nz, d, w)
    
  END SUBROUTINE reallocate_field_dp_3D
  SUBROUTINE reallocate_field_int(   nV, d, w)
    ! For data that doesn't need to be mapped because it will be recalculated anyway.
    
    IMPLICIT NONE
    
    ! In/output variables
    INTEGER,                             INTENT(IN)    :: nV
    INTEGER,  DIMENSION(:  ), POINTER,   INTENT(INOUT) :: d            ! Pointer to the data
    INTEGER,                             INTENT(INOUT) :: w            ! MPI window to the shared memory space containing that data
    
    ! Deallocate and reallocate the shared memory space
    CALL deallocate_shared(w)
    NULLIFY(d)    
    CALL allocate_shared_int_1D( nV, d, w)
    
  END SUBROUTINE reallocate_field_int
  
END MODULE mesh_mapping_module
