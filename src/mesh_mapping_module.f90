MODULE mesh_mapping_module

  ! Routines for creating mapping arrays and mapping data between meshes and grids

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
  USE parallel_module,                 ONLY: adapt_shared_int_1D,    adapt_shared_dp_1D, &
                                             adapt_shared_int_2D,    adapt_shared_dp_2D, &
                                             adapt_shared_int_3D,    adapt_shared_dp_3D, &
                                             adapt_shared_bool_1D                   
  USE data_types_module,               ONLY: type_mesh, type_remapping, type_remapping_trilin, type_remapping_nearest_neighbour, &
                                             type_remapping_conservative, type_grid, type_latlongrid, type_remapping_latlon2mesh, &
                                             type_remapping_conservative_intermediate_Ac_local, type_remapping_conservative_intermediate_Ac_shared, &
                                             type_remapping_conservative_intermediate_shared, type_remapping_conservative_intermediate_local
  USE mesh_help_functions_module,      ONLY: is_in_triangle, find_Voronoi_cell_vertices, find_containing_triangle, find_triangle_area, find_containing_vertex, &
                                             line_integral_xdy, line_integral_mxydx, line_integral_xydy, is_boundary_segment, cross2, &
                                             lies_on_line_segment, segment_intersection, partition_domain_x_balanced, write_mesh_to_text_file, &
                                             line_from_points, line_line_intersection
  USE mesh_operators_module,           ONLY: ddx_a_to_a_2D, ddy_a_to_a_2D, ddx_a_to_a_3D, ddy_a_to_a_3D
  USE utilities_module,                ONLY: smooth_Gaussian_2D_grid, smooth_Gaussian_3D_grid

  IMPLICIT NONE

  CONTAINS

  ! == Subroutines for mapping data between a square grid and the mesh ==
  SUBROUTINE create_remapping_arrays_mesh_grid( mesh, grid)
    ! Create remapping arrays for remapping data between a square grid and the model mesh
    ! using pseudo-conservative remapping
    
    USE parallel_module, ONLY: allocate_shared_dist_int_1D, allocate_shared_dist_int_2D, &
                               allocate_shared_dist_dp_1D, allocate_shared_dist_dp_2D
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh
    TYPE(type_grid),                         INTENT(INOUT) :: grid
    
    ! Local variables
    CHARACTER(LEN=64), PARAMETER                           :: routine_name = 'create_remapping_arrays_mesh_grid'
    INTEGER                                                :: nm1, nm2
    INTEGER                                                :: nmax, nmax_list, n, vi
    INTEGER                                                :: n_cells_in_vertex, il, iu, jl, ju, i, j, ni
    REAL(dp)                                               :: V_dx, V_xl, V_xu, V_yl, V_yu
    REAL(dp)                                               :: xo_min, xo_max, yo_min, yo_max, A_overlap
    INTEGER,  DIMENSION(:    ), POINTER                    :: i_vi, j_vi
    REAL(dp), DIMENSION(:    ), POINTER                    :: A_vi
    INTEGER :: wi_vi, wj_vi, wA_vi
    INTEGER,  DIMENSION(:,:  ), POINTER                    :: ii
    REAL(dp), DIMENSION(:    ), POINTER                    :: A, w_m2g, w_g2m
    INTEGER :: wii, wA, ww_m2g, ww_g2m
    REAL(dp), DIMENSION(:    ), POINTER                    :: A_tot_mesh
    REAL(dp), DIMENSION(:,:  ), POINTER                    :: A_tot_grid
    INTEGER :: wA_tot_mesh, wA_tot_grid
    LOGICAL,  DIMENSION(:,:  ), ALLOCATABLE                :: grid_islisted
    LOGICAL                                                :: do_extend_memory
    REAL(dp), DIMENSION(2)                                 :: pv, pa, pb, pc
    INTEGER                                                :: ti, via, vib, vic
    REAL(dp)                                               :: Atri, Aa, Ab, Ac, wpa, wpb, wpc
    INTEGER                                                :: p, n1, n2, n_tot
    
    nm1 = par%mem%n
    
    ! Determine maximum number of contributions
    nmax = 4 * CEILING( MAX( MAXVAL(mesh%A) / grid%dx**2, grid%dx**2 / MINVAL(mesh%A)))
    
    CALL allocate_shared_dist_int_1D( nmax, i_vi, wi_vi)
    CALL allocate_shared_dist_int_1D( nmax, j_vi, wj_vi)
    CALL allocate_shared_dist_dp_1D(  nmax, A_vi, wA_vi)
    
    ! Allocate memory for process-local lists
    n = 0
    nmax_list = CEILING( REAL(mesh%nV + grid%nx * grid%ny, dp) / REAL(par%n,dp) )
    CALL allocate_shared_dist_int_2D( nmax_list, 3, ii,    wii   )
    CALL allocate_shared_dist_dp_1D(  nmax_list,    A,     wA    )
    CALL allocate_shared_dist_dp_1D(  nmax_list,    w_m2g, ww_m2g)
    CALL allocate_shared_dist_dp_1D(  nmax_list,    w_g2m, ww_g2m)
    
    ! Allocate memory for sum maps
    CALL allocate_shared_dist_dp_1D( mesh%nV,          A_tot_mesh, wA_tot_mesh)
    CALL allocate_shared_dist_dp_2D( grid%nx, grid%ny, A_tot_grid, wA_tot_grid)
    ALLOCATE( grid_islisted( grid%nx, grid%ny))
    
    A_tot_mesh    = 0._dp
    A_tot_grid    = 0._dp
    grid_islisted = .FALSE.
    
    ! Calculate overlaps for all vertices
    DO vi = 1, mesh%nV! mesh%vi1, mesh%vi2
    
      IF (vi >= mesh%vi1 .AND. vi <= mesh%vi2) THEN

      ! Reset temporary data
      n_cells_in_vertex = 0
      i_vi = 0
      j_vi = 0
      A_vi = 0._dp

      ! Draw a square of the same area as this vertex' Voronoi cell
      IF (mesh%edge_index( vi) == 0) THEN
        V_dx = SQRT(        mesh%A( vi))
      ELSE
        V_dx = SQRT(2._dp * mesh%A( vi)) ! To compensate for the fact that half of a boundary vertex' Voronoi cell is "missing"
      END IF
      
      V_xl = mesh%V( vi,1) - V_dx / 2._dp
      V_xu = mesh%V( vi,1) + V_dx / 2._dp
      V_yl = mesh%V( vi,2) - V_dx / 2._dp
      V_yu = mesh%V( vi,2) + V_dx / 2._dp

      ! Find all grid cells that overlap with this square
      il = MAX( 1,       CEILING(-1.5_dp + REAL(FLOOR(grid%nx / 2._dp),dp) + V_xl / grid%dx) )
      iu = MIN( grid%nx, CEILING( 1.5_dp + REAL(FLOOR(grid%nx / 2._dp),dp) + V_xu / grid%dx) )
      jl = MAX( 1,       CEILING(-1.5_dp + REAL(FLOOR(grid%ny / 2._dp),dp) + V_yl / grid%dx) )
      ju = MIN( grid%ny, CEILING( 1.5_dp + REAL(FLOOR(grid%ny / 2._dp),dp) + V_yu / grid%dx) )

      ! Find areas of overlap with all of these grid cells
      DO i = il, iu
      DO j = jl, ju

        xo_min = MAX( V_xl, grid%x( i) - grid%dx / 2._dp)
        xo_max = MIN( V_xu, grid%x( i) + grid%dx / 2._dp)
        yo_min = MAX( V_yl, grid%y( j) - grid%dx / 2._dp)
        yo_max = MIN( V_yu, grid%y( j) + grid%dx / 2._dp)

        ! Calculate area of overlap
        A_overlap = (xo_max - xo_min) * (yo_max - yo_min)

        ! If there's no overlap, skip this grid cell.
        IF (xo_max <= xo_min .OR. yo_max <= yo_min .OR. A_overlap == 0._dp) CYCLE
    
        ! Fill into temporary mapping array
        n_cells_in_vertex = n_cells_in_vertex + 1
        i_vi( n_cells_in_vertex) = i
        j_vi( n_cells_in_vertex) = j
        A_vi( n_cells_in_vertex) = A_overlap
        
        A_tot_mesh(    vi ) = A_tot_mesh( vi ) + A_overlap
        A_tot_grid(    i,j) = A_tot_grid( i,j) + A_overlap
        grid_islisted( i,j) = .TRUE.

      END DO ! DO j = jl, ju
      END DO ! DO i = il, iu

      ! Add data to list
      DO ni = 1, n_cells_in_vertex
        n = n + 1
        ii( n,:) = [vi, i_vi( ni), j_vi( ni)]
        A(  n  ) = A_vi( ni)
      END DO
      
      END IF ! IF (vi >= mesh%vi1 .AND. vi <= mesh%vi2) THEN

      ! Expend list memory if necessary
      do_extend_memory = .FALSE.
      IF (n > nmax_list - nmax) do_extend_memory = .TRUE.
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, do_extend_memory, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      IF (do_extend_memory) THEN
        nmax_list = MAX( n + 10*nmax, MAX( nmax_list, CEILING( REAL(n,dp) * 1.2_dp)))
        CALL create_remapping_arrays_mesh_grid_expand_list_memory( n, nmax_list, ii, wii, A, wA, w_m2g, ww_m2g, w_g2m, ww_g2m)
      END IF
      
    END DO
    CALL sync
    
    ! Calculate weights
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, A_tot_mesh, mesh%nV,         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, A_tot_grid, grid%nx*grid%ny, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    DO ni = 1, n
      vi = ii( ni,1)
      i  = ii( ni,2)
      j  = ii( ni,3)
      w_m2g( ni) = A( ni) / A_tot_grid( i,j)
      w_g2m( ni) = A( ni) / A_tot_mesh( vi )
    END DO
    CALL sync

    ! Find grid cells that are not yet listed; use trilinear interpolation for those.
    ! Only the master adds entries to the lists, but all processes must run the loops
    ! so the shared memory extension routines can be called.
    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, grid_islisted, grid%nx*grid%ny, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    
    ti = 1
    DO i = 1, grid%nx
    DO j = 1, grid%ny
      
      ! Only do this for grid cells that have not been linked to any mesh vertex.
      IF (grid_islisted( i,j)) CYCLE

      ! Find the mesh triangle containing this grid cell
      pv = [MIN( mesh%xmax, MAX( mesh%xmin, grid%x(i) )), &
            MIN( mesh%ymax, MAX( mesh%ymin, grid%y(j) ))]
      CALL find_containing_triangle( mesh, pv, ti)

      ! Trilinear interpolation
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      pa  = mesh%V( via,:)
      pb  = mesh%V( vib,:)
      pc  = mesh%V( vic,:)

      Atri = cross2( (pb-pa), (pc-pa) )
      Aa   = cross2( (pb-pv), (pc-pv) )
      Ab   = cross2( (pc-pv), (pa-pv) )
      Ac   = cross2( (pa-pv), (pb-pv) )

      wpa  = Aa / Atri
      wpb  = Ab / Atri
      wpc  = Ac / Atri

      ! List weights
      IF (par%master) THEN
        n = n + 1
        ii(    n,:) = [via,i,j]
        w_m2g( n  ) = wpa
        w_g2m( n  ) = 0._dp
        n = n + 1
        ii(    n,:) = [vib,i,j]
        w_m2g( n  ) = wpb
        w_g2m( n  ) = 0._dp
        n = n + 1
        ii(    n,:) = [vic,i,j]
        w_m2g( n  ) = wpc
        w_g2m( n  ) = 0._dp
      END IF
      CALL sync

      ! Expend list memory if necessary
      do_extend_memory = .FALSE.
      IF (n > nmax_list - nmax) do_extend_memory = .TRUE.
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, do_extend_memory, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      IF (do_extend_memory) THEN
        nmax_list = MAX(n + 10*nmax, MAX( nmax_list, CEILING( REAL(n,dp) * 1.2_dp)))
        CALL create_remapping_arrays_mesh_grid_expand_list_memory( n, nmax_list, ii, wii, A, wA, w_m2g, ww_m2g, w_g2m, ww_g2m)
      END IF

    END DO ! DO j = 1, grid%ny
    END DO ! DO i = 1, grid%nx
    
    ! Clean up after yourself
    CALL deallocate_shared( wi_vi)
    CALL deallocate_shared( wj_vi)
    CALL deallocate_shared( wA_vi)
    CALL deallocate_shared( wA)
    CALL deallocate_shared( wA_tot_mesh)
    CALL deallocate_shared( wA_tot_grid)
    DEALLOCATE( grid_islisted)
    
    ! Gather data from all the processes
    ! ==================================
    
    ! Allocate shared memory to hold the merged data
    n_tot = n
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    ! n
    n1 = 0
    n2 = 0
    CALL allocate_shared_int_0D(           grid%map%n,     grid%map%wn    )
    DO p = 0, par%n-1
      IF (p == par%i) THEN
        n1 = grid%map%n + 1
        n2 = grid%map%n + n
        grid%map%n = grid%map%n + n
      END IF
      CALL sync
    END DO
    
    ! ii
    CALL allocate_shared_int_2D( n_tot, 3, grid%map%ii,    grid%map%wii   )
    DO p = 0, par%n-1
      IF (p == par%i) THEN
        grid%map%ii(    n1:n2,:) = ii(    1:n,:)
      END IF
      CALL sync
    END DO
    CALL deallocate_shared( wii)
    
    ! w_m2g
    CALL allocate_shared_dp_1D(  n_tot,    grid%map%w_m2g, grid%map%ww_m2g)
    DO p = 0, par%n-1
      IF (p == par%i) THEN
        grid%map%w_m2g( n1:n2  ) = w_m2g( 1:n  )
      END IF
      CALL sync
    END DO
    CALL deallocate_shared( ww_m2g)
    
    ! w_g2m
    CALL allocate_shared_dp_1D(  n_tot,    grid%map%w_g2m, grid%map%ww_g2m)
    DO p = 0, par%n-1
      IF (p == par%i) THEN
        grid%map%w_g2m( n1:n2  ) = w_g2m( 1:n  )
      END IF
      CALL sync
    END DO
    CALL deallocate_shared( ww_g2m)
    
    ! Update memory tracker
    nm2 = par%mem%n
    CALL write_to_memory_log( routine_name, nm1, nm2)
    
  END SUBROUTINE create_remapping_arrays_mesh_grid
  SUBROUTINE create_remapping_arrays_mesh_grid_expand_list_memory( n, nmax_list, ii, wii, A, wA, w_m2g, ww_m2g, w_g2m, ww_g2m)
  
    USE parallel_module, ONLY: adapt_shared_dist_int_2D, adapt_shared_dist_dp_1D
    
    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                                 INTENT(IN   ) :: n
    INTEGER,                                 INTENT(IN   ) :: nmax_list
    INTEGER,  DIMENSION(:,:  ), POINTER,     INTENT(INOUT) :: ii
    REAL(dp), DIMENSION(:    ), POINTER,     INTENT(INOUT) :: A, w_m2g, w_g2m
    INTEGER,                                 INTENT(INOUT) :: wii, wA, ww_m2g, ww_g2m
    
    CALL adapt_shared_dist_int_2D(  n, nmax_list, 3, ii,    wii   )
    CALL adapt_shared_dist_dp_1D(   n, nmax_list,    A,     wA    )
    CALL adapt_shared_dist_dp_1D(   n, nmax_list,    w_m2g, ww_m2g)
    CALL adapt_shared_dist_dp_1D(   n, nmax_list,    w_g2m, ww_g2m)
  
  END SUBROUTINE create_remapping_arrays_mesh_grid_expand_list_memory
  SUBROUTINE map_grid2mesh_2D(     mesh, grid, d_grid, d_mesh)
    ! Remapping data from a square grid to the model mesh using pseudo-conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: d_grid
    REAL(dp), DIMENSION(:    ),              INTENT(OUT)   :: d_mesh
    
    ! Local variables:
    INTEGER                                                :: ii,vi,i,j
    
    ii = mesh%nV
    
    IF (par%master) THEN
      d_mesh = 0._dp
      DO ii = 1, grid%map%n
        vi = grid%map%ii( ii,1)
        i  = grid%map%ii( ii,2)
        j  = grid%map%ii( ii,3)
        d_mesh( vi) = d_mesh( vi) + grid%map%w_g2m( ii) * d_grid( i,j)
      END DO
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE map_grid2mesh_2D
  SUBROUTINE map_grid2mesh_3D(     mesh, grid, d_grid, d_mesh)
    ! Remapping data from a square grid to the model mesh using pseudo-conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: d_grid
    REAL(dp), DIMENSION(:,:  ),              INTENT(OUT)   :: d_mesh
    
    ! Local variables:
    INTEGER                                                :: ii,vi,i,j
    
    ii = mesh%nV
    
    IF (par%master) THEN
      d_mesh = 0._dp
      DO ii = 1, grid%map%n
        vi = grid%map%ii( ii,1)
        i  = grid%map%ii( ii,2)
        j  = grid%map%ii( ii,3)
        d_mesh( vi,:) = d_mesh( vi,:) + grid%map%w_g2m( ii) * d_grid( i,j,:)
      END DO
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE map_grid2mesh_3D
  SUBROUTINE map_mesh2grid_2D(     mesh, grid, d_mesh, d_grid)
    ! Remapping data from the model mesh to a square grid using pseudo-conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: d_mesh
    REAL(dp), DIMENSION(:,:  ),              INTENT(OUT)   :: d_grid
    
    ! Local variables:
    INTEGER                                                :: ii,vi,i,j
    
    ii = mesh%nV
    
    IF (par%master) THEN
      d_grid = 0._dp
      DO ii = 1, grid%map%n
        vi = grid%map%ii( ii,1)
        i  = grid%map%ii( ii,2)
        j  = grid%map%ii( ii,3)
        d_grid( i,j) = d_grid( i,j) + grid%map%w_m2g( ii) * d_mesh( vi)
      END DO
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE map_mesh2grid_2D
  SUBROUTINE map_mesh2grid_2D_min( mesh, grid, d_mesh, d_grid)
    ! Remapping data from the model mesh to a square grid using pseudo-conservative remapping
    ! Takes the minimum value over all contributing values, to be used for inspecting mesh resolution
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: d_mesh
    REAL(dp), DIMENSION(:,:  ),              INTENT(OUT)   :: d_grid
    
    ! Local variables:
    INTEGER                                                :: ii,vi,i,j
    
    ii = mesh%nV
    
    IF (par%master) THEN
      d_grid = MAXVAL( d_mesh)
      DO ii = 1, grid%map%n
        vi = grid%map%ii( ii,1)
        i  = grid%map%ii( ii,2)
        j  = grid%map%ii( ii,3)
        d_grid( i,j) = MIN(d_grid( i,j), d_mesh( vi))
      END DO
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE map_mesh2grid_2D_min
  SUBROUTINE map_mesh2grid_2D_max( mesh, grid, d_mesh, d_grid)
    ! Remapping data from the model mesh to a square grid using pseudo-conservative remapping
    ! Takes the minimum value over all contributing values, to be used for inspecting mesh resolution
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: d_mesh
    REAL(dp), DIMENSION(:,:  ),              INTENT(OUT)   :: d_grid
    
    ! Local variables:
    INTEGER                                                :: ii,vi,i,j
    
    ii = mesh%nV
    
    IF (par%master) THEN
      d_grid = MINVAL( d_mesh)
      DO ii = 1, grid%map%n
        vi = grid%map%ii( ii,1)
        i  = grid%map%ii( ii,2)
        j  = grid%map%ii( ii,3)
        d_grid( i,j) = MAX(d_grid( i,j), d_mesh( vi))
      END DO
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE map_mesh2grid_2D_max
  SUBROUTINE map_mesh2grid_3D(     mesh, grid, d_mesh, d_grid)
    ! Remapping data from the model mesh to a square grid using pseudo-conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: d_mesh
    REAL(dp), DIMENSION(:,:,:),              INTENT(OUT)   :: d_grid
    
    ! Local variables:
    INTEGER                                                :: ii,vi,i,j
    
    ii = mesh%nV
    
    IF (par%master) THEN
      d_grid = 0._dp
      DO ii = 1, grid%map%n
        vi = grid%map%ii( ii,1)
        i  = grid%map%ii( ii,2)
        j  = grid%map%ii( ii,3)
        d_grid( i,j,:) = d_grid( i,j,:) + grid%map%w_m2g( ii) * d_mesh( vi,:)
      END DO
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE map_mesh2grid_3D
  SUBROUTINE deallocate_remapping_arrays_mesh_grid( grid)
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                         INTENT(INOUT) :: grid
    
    CALL deallocate_shared( grid%map%wn)
    CALL deallocate_shared( grid%map%wii)
    CALL deallocate_shared( grid%map%ww_m2g)
    CALL deallocate_shared( grid%map%ww_g2m)
    
  END SUBROUTINE deallocate_remapping_arrays_mesh_grid
  
  ! == Subroutine for mapping data from a global lat-lon grid and the mesh ==
  SUBROUTINE create_remapping_arrays_glob_mesh( mesh, grid, map)
    ! Create remapping arrays for remapping data from a global lat-lon grid to the model mesh
    ! using bilinear interpolation
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_latlongrid),                   INTENT(IN)    :: grid
    TYPE(type_remapping_latlon2mesh),        INTENT(INOUT) :: map
    
    INTEGER                                                :: vi
    INTEGER                                                :: il,iu,jl,ju
    REAL(dp)                                               :: wil,wiu,wjl,wju
    
    ! Allocate shared memory for the mapping arrays
    CALL allocate_shared_int_1D( mesh%nV, map%ilat1, map%wilat1)
    CALL allocate_shared_int_1D( mesh%nV, map%ilat2, map%wilat2)
    CALL allocate_shared_int_1D( mesh%nV, map%ilon1, map%wilon1)
    CALL allocate_shared_int_1D( mesh%nV, map%ilon2, map%wilon2)
    CALL allocate_shared_dp_1D(  mesh%nV, map%wlat1, map%wwlat1)
    CALL allocate_shared_dp_1D(  mesh%nV, map%wlat2, map%wwlat2)
    CALL allocate_shared_dp_1D(  mesh%nV, map%wlon1, map%wwlon1)
    CALL allocate_shared_dp_1D(  mesh%nV, map%wlon2, map%wwlon2)

    DO vi = mesh%vi1, mesh%vi2
      
      ! Find enveloping lat-lon indices
      il  = MAX(1, MIN( grid%nlon-1, 1 + FLOOR((mesh%lon( vi) - MINVAL(grid%lon)) / (grid%lon(2) - grid%lon(1)))))
      iu  = il + 1        
      wil = (grid%lon(iu) - mesh%lon( vi)) / (grid%lon(2) - grid%lon(1))
      wiu = 1._dp - wil
      
      ! Exception for pixels near the zero meridian
      IF (mesh%lon( vi) < MINVAL(grid%lon)) THEN
        il  = grid%nlon
        iu  = 1      
        wil = (grid%lon( iu) - mesh%lon( vi)) / (grid%lon(2) - grid%lon(1))
        wiu = 1._dp - wil
      ELSEIF (mesh%lon( vi) > MAXVAL(grid%lon)) THEN
        il  = grid%nlon
        iu  = 1
        wiu = (mesh%lon( vi) - grid%lon( il)) / (grid%lon(2) - grid%lon(1))
        wil = 1._dp - wiu
      END IF
          
      jl  = MAX(1, MIN( grid%nlat-1, 1 + FLOOR((mesh%lat( vi) - MINVAL(grid%lat)) / (grid%lat(2) - grid%lat(1)))))
      ju  = jl + 1        
      wjl = (grid%lat( ju) - mesh%lat( vi)) / (grid%lat(2) - grid%lat(1))
      wju = 1 - wjl

      ! Write to mapping arrays
      map%ilon1( vi) = il
      map%ilon2( vi) = iu
      map%ilat1( vi) = jl
      map%ilat2( vi) = ju
      map%wlon1( vi) = wil
      map%wlon2( vi) = wiu
      map%wlat1( vi) = wjl
      map%wlat2( vi) = wju

    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
  END SUBROUTINE create_remapping_arrays_glob_mesh
  SUBROUTINE map_latlon2mesh_2D( mesh, map, d_grid, d_mesh)
    ! Map data from a global lat-lon grid to the model mesh using bilinear interpolation
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_remapping_latlon2mesh),        INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: d_grid
    REAL(dp), DIMENSION(:    ),              INTENT(OUT)   :: d_mesh
    
    ! Local variables:
    INTEGER                                                :: vi
    INTEGER                                                :: il,iu,jl,ju
    REAL(dp)                                               :: wil,wiu,wjl,wju
    
    DO vi = mesh%vi1, mesh%vi2
    
      il  = map%ilon1( vi)
      iu  = map%ilon2( vi)
      jl  = map%ilat1( vi)
      ju  = map%ilat2( vi)
      wil = map%wlon1( vi)
      wiu = map%wlon2( vi)
      wjl = map%wlat1( vi)
      wju = map%wlat2( vi)
      
      d_mesh( vi) = (wil * wjl * d_grid( il,jl)) + &
                    (wil * wju * d_grid( il,ju)) + &
                    (wiu * wjl * d_grid( iu,jl)) + &
                    (wiu * wju * d_grid( iu,ju))
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
  END SUBROUTINE map_latlon2mesh_2D
  SUBROUTINE map_latlon2mesh_3D( mesh, map, d_grid, d_mesh)
    ! Map data from a global lat-lon grid to the model mesh using bilinear interpolation
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_remapping_latlon2mesh),        INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: d_grid
    REAL(dp), DIMENSION(:,:  ),              INTENT(OUT)   :: d_mesh
    
    ! Local variables:
    INTEGER                                                :: vi
    INTEGER                                                :: il,iu,jl,ju
    REAL(dp)                                               :: wil,wiu,wjl,wju
    
    DO vi = mesh%vi1, mesh%vi2
    
      il  = map%ilon1( vi)
      iu  = map%ilon2( vi)
      jl  = map%ilat1( vi)
      ju  = map%ilat2( vi)
      wil = map%wlon1( vi)
      wiu = map%wlon2( vi)
      wjl = map%wlat1( vi)
      wju = map%wlat2( vi)
      
      d_mesh( vi,:) = (wil * wjl * d_grid( il,jl,:)) + &
                      (wil * wju * d_grid( il,ju,:)) + &
                      (wiu * wjl * d_grid( iu,jl,:)) + &
                      (wiu * wju * d_grid( iu,ju,:))
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
  END SUBROUTINE map_latlon2mesh_3D
  SUBROUTINE deallocate_remapping_arrays_glob_mesh( map)
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_remapping_latlon2mesh),        INTENT(INOUT) :: map
    
    CALL deallocate_shared( map%wilat1)
    CALL deallocate_shared( map%wilat2)
    CALL deallocate_shared( map%wilon1)
    CALL deallocate_shared( map%wilon2)
    CALL deallocate_shared( map%wwlat1)
    CALL deallocate_shared( map%wwlat2)
    CALL deallocate_shared( map%wwlon1)
    CALL deallocate_shared( map%wwlon2)
    
  END SUBROUTINE deallocate_remapping_arrays_glob_mesh

  ! == Subroutines creating different kinds of remapping arrays
  SUBROUTINE create_remapping_arrays( mesh_src, mesh_dst, map)
    ! Create remapping arrays for remapping data from mesh_src to mesh_dst, for all remapping methods
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_src
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_dst
    TYPE(type_remapping),                    INTENT(INOUT) :: map
    
    ! Local variables:
    CHARACTER(LEN=64), PARAMETER                           :: routine_name = 'create_remapping_arrays'
    INTEGER                                                :: n1, n2
    
    n1 = par%mem%n
        
    ! Create all remapping arrays
    CALL create_remapping_arrays_trilin(            mesh_src, mesh_dst, map%trilin)
    CALL create_remapping_arrays_nearest_neighbour( mesh_src, mesh_dst, map%nearest_neighbour)
    CALL create_remapping_arrays_conservative(      mesh_src, mesh_dst, map%conservative)
    
    n2 = par%mem%n
    CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE create_remapping_arrays
  SUBROUTINE deallocate_remapping_arrays( map)
    
    IMPLICIT NONE
  
    TYPE(type_remapping),                    INTENT(INOUT) :: map
    
    CALL deallocate_shared( map%trilin%wvi)
    CALL deallocate_shared( map%trilin%ww)   
    
    CALL deallocate_shared( map%nearest_neighbour%wvi)
    
    CALL deallocate_shared( map%conservative%wn_tot)
    CALL deallocate_shared( map%conservative%wvli1)
    CALL deallocate_shared( map%conservative%wvli2)
    CALL deallocate_shared( map%conservative%wvi)
    CALL deallocate_shared( map%conservative%ww0)
    CALL deallocate_shared( map%conservative%ww1x)
    CALL deallocate_shared( map%conservative%ww1y)
    
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
    
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
    
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
      CALL find_triangle_area( pa, pb, pc, Atot)
      CALL find_triangle_area( pb, pc, p , Aa  )
      CALL find_triangle_area( pc, pa, p , Ab  )
      CALL find_triangle_area( pa, pb, p , Ac  )
      
      map%vi( vi_dst,:) = [via_src, vib_src, vic_src]
      map%w(  vi_dst,:) = [Aa/Atot, Ab/Atot, Ac/Atot]
       
    END DO ! DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
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
    
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
    
      ! The vertex coordinates
      p = mesh_dst%V(vi_dst,:)
    
      ! Find the mesh_src triangle containing this vertex
      CALL find_containing_vertex(mesh_src, p, vi_src)
      
      map%vi( vi_dst) = vi_src
      
    END DO ! DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
    CALL sync
    
  END SUBROUTINE create_remapping_arrays_nearest_neighbour
  
! == Subroutines for conservative remapping
  SUBROUTINE create_remapping_arrays_conservative( mesh_src, mesh_dst, map)
    ! Create remapping arrays for remapping data from mesh_src to mesh_dst using 1st and 2nd order conservative remapping
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT)   :: mesh_src ! INOUT instead of IN because the search maps and stacks are used
    TYPE(type_mesh),                         INTENT(INOUT)   :: mesh_dst
    TYPE(type_remapping_conservative),       INTENT(INOUT)   :: map
    
        
    REAL(dp), DIMENSION(:,:,:), POINTER                      :: Vor_lines_src, Vor_lines_dst
    INTEGER,  DIMENSION(:,:  ), POINTER                      :: Vor_vi_ti_src, Vor_vi_ti_dst
    INTEGER                                                  :: wVor_lines_src, wVor_lines_dst, wVor_vi_ti_src, wVor_vi_ti_dst
    
    INTEGER                                                  :: vi_src_start_proc, vi_dst_start_proc
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                  :: proc_domain_Ac_src
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                  :: proc_domain_Ac_dst
    
    REAL(dp)                                                 :: A_min, A_max
    INTEGER                                                  :: nV_max, nS_max, n_max_src, n_max_dst
    
    ! Results from integrating over Voronoi boundary lines for both meshes
    TYPE(type_remapping_conservative_intermediate_Ac_local)  :: r_src_Ac_proc
    TYPE(type_remapping_conservative_intermediate_Ac_local)  :: r_dst_Ac_proc
    TYPE(type_remapping_conservative_intermediate_Ac_shared) :: r_src_Ac
    TYPE(type_remapping_conservative_intermediate_Ac_shared) :: r_dst_Ac
    
    LOGICAL                                                  :: CountCoincidences
    
    ! Integration results rearranged to actual vertices
    TYPE(type_remapping_conservative_intermediate_local)     :: r_src_proc
    TYPE(type_remapping_conservative_intermediate_local)     :: r_dst_proc
    TYPE(type_remapping_conservative_intermediate_shared)    :: r_src
    TYPE(type_remapping_conservative_intermediate_shared)    :: r_dst
    
    INTEGER                                                  :: n_tot_Ac_src, n_tot_Ac_dst, n_tot_src, n_tot_dst
    INTEGER(KIND=MPI_ADDRESS_KIND)                           :: windowsize
    
    ! ================================================================================================================================
    
    !CALL write_mesh_to_text_file( mesh_src, 'mesh_src.txt')
    !CALL write_mesh_to_text_file( mesh_dst, 'mesh_dst.txt')
    
    A_min  = MIN( MINVAL( mesh_src%A), MINVAL( mesh_dst%A))
    A_max  = MAX( MAXVAL( mesh_src%A), MAXVAL( mesh_dst%A))
    nV_max = CEILING( A_max / A_min)
    nS_max = CEILING( SQRT( REAL( nV_max, dp)))
    
    ! Find the coordinates and relevant indices of the Voronoi boundary lines.
    ! Since each line describes the boundary between the Voronoi cells of two
    ! connected vertices, each Voronoi line can be uniquely described by an Aci index.
    
    CALL allocate_shared_dp_3D(  mesh_src%nAc, 2, 2, Vor_lines_src, wVor_lines_src)
    CALL allocate_shared_dp_3D(  mesh_dst%nAc, 2, 2, Vor_lines_dst, wVor_lines_dst)
    CALL allocate_shared_int_2D( mesh_src%nAc, 6,    Vor_vi_ti_src, wVor_vi_ti_src)
    CALL allocate_shared_int_2D( mesh_dst%nAc, 6,    Vor_vi_ti_dst, wVor_vi_ti_dst)
    
    CALL find_Voronoi_boundary_lines( mesh_src, Vor_lines_src, Vor_vi_ti_src)
    CALL find_Voronoi_boundary_lines( mesh_dst, Vor_lines_dst, Vor_vi_ti_dst)
    
    ! Determine process vertex domains for both meshes
    
    ALLOCATE( proc_domain_Ac_src( mesh_src%nAc))
    ALLOCATE( proc_domain_Ac_dst( mesh_dst%nAc))
  
    CALL determine_process_domains( mesh_src, proc_domain_Ac_src, vi_src_start_proc)
    CALL determine_process_domains( mesh_dst, proc_domain_Ac_dst, vi_dst_start_proc)

    ! Integrate over all Voronoi cell boundary lines for both meshes
    
    n_max_src = CEILING( 0.7_dp * REAL(mesh_src%nAc,dp))
    n_max_dst = CEILING( 0.7_dp * REAL(mesh_dst%nAc,dp))
    
    CALL allocate_memory_Ac_local( r_src_Ac_proc, n_max_src, mesh_src%nAc)
    CALL allocate_memory_Ac_local( r_dst_Ac_proc, n_max_dst, mesh_dst%nAc)
  
    CountCoincidences = .FALSE.
    CALL integrate_over_Voronoi_boundaries( mesh_src, mesh_dst, nV_max, nS_max, proc_domain_Ac_src, &
      vi_src_start_proc, Vor_lines_src, Vor_lines_dst, Vor_vi_ti_dst, r_src_Ac_proc, CountCoincidences)
    ! Gather data from the parallel processes into shared memory (deallocates the process-local memory)
    CALL gather_memory_Ac( mesh_src, r_src_Ac_proc, r_src_Ac, proc_domain_Ac_src)
    n_tot_Ac_src = r_src_Ac%n_tot
  
    CountCoincidences = .TRUE.  
    CALL integrate_over_Voronoi_boundaries( mesh_dst, mesh_src, nV_max, nS_max, proc_domain_Ac_dst, &
      vi_dst_start_proc, Vor_lines_dst, Vor_lines_src, Vor_vi_ti_src, r_dst_Ac_proc, CountCoincidences)
    ! Gather data from the parallel processes into shared memory (deallocates the process-local memory)
    CALL gather_memory_Ac( mesh_dst, r_dst_Ac_proc, r_dst_Ac, proc_domain_Ac_dst)
    n_tot_Ac_dst = r_dst_Ac%n_tot
    
    ! Reorder src line integral lists from edge-based to vertex-based (deallocates edge-based data)
    CALL reorder_contributions_from_lines_to_vertices( mesh_src, mesh_dst, nV_max, r_src_Ac, r_src_proc)
    ! Gather data from the parallel processes into shared memory (deallocates the process-local memory)
    CALL gather_memory( mesh_src, r_src_proc, r_src)
    n_tot_src = r_src%n_tot
    
    ! Reorder dst line integral lists from edge-based to vertex-based (deallocates edge-based data)
    CALL reorder_contributions_from_lines_to_vertices( mesh_dst, mesh_src, nV_max, r_dst_Ac, r_dst_proc)
    ! Gather data from the parallel processes into shared memory (deallocates the process-local memory)
    CALL gather_memory( mesh_dst, r_dst_proc, r_dst)
    n_tot_dst = r_dst%n_tot
    
    DEALLOCATE( proc_domain_Ac_src)
    DEALLOCATE( proc_domain_Ac_dst)
    CALL deallocate_shared( wVor_lines_src)
    CALL deallocate_shared( wVor_lines_dst)
    CALL deallocate_shared( wVor_vi_ti_src)
    CALL deallocate_shared( wVor_vi_ti_dst)
    
    ! Add integral data from the mesh_src Voronoi boundaries to dst lists (deallocates r_src)
    CALL add_entries_from_opposite_mesh( mesh_dst, mesh_src, r_dst, r_src)
    
    ! Calculate the remapping weights from the line integrals (deallocates r_dst)
    CALL calculate_remapping_weights_from_line_integrals( mesh_dst, mesh_src, r_dst, map)
      
    ! Check if everything worked
    CALL check_if_remapping_is_conservative( mesh_src, mesh_dst, map)
    
    ! Update memory use tracker
    ! Since a lot of the memory used here it not MPI shared memory, the standard allocate_shared routines
    ! do not include this memory, and os memory use is underestimated. Adding separate statements next to
    ! all the allocate/deallocate calls is impractical, so instead just make a best estimate of how much is used.
    
    windowsize = (3 * mesh_src%nAc + 2 * n_tot_Ac_src) * 4_MPI_ADDRESS_KIND + &
                 (                   3 * n_tot_Ac_src) * 8_MPI_ADDRESS_KIND + &
                 (3 * mesh_dst%nAc + 2 * n_tot_Ac_dst) * 4_MPI_ADDRESS_KIND + &
                 (                   3 * n_tot_Ac_dst) * 8_MPI_ADDRESS_KIND + &
                 (3 * mesh_src%nV  + 2 * n_tot_src   ) * 4_MPI_ADDRESS_KIND + &
                 (                   3 * n_tot_src   ) * 8_MPI_ADDRESS_KIND + &
                 (3 * mesh_dst%nV  + 2 * n_tot_dst   ) * 4_MPI_ADDRESS_KIND + &
                 (                   3 * n_tot_dst   ) * 8_MPI_ADDRESS_KIND + &
                 (3 * mesh_dst%nV  + 2 * n_tot_dst   ) * 4_MPI_ADDRESS_KIND + &
                 (                   3 * n_tot_dst   ) * 8_MPI_ADDRESS_KIND
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, windowsize, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF (par%master) THEN
      par%mem%total = par%mem%total + windowsize
      par%mem%n = par%mem%n + 1
      par%mem%h( par%mem%n) = par%mem%total
      par%mem%total = par%mem%total - windowsize
      par%mem%n = par%mem%n + 1
      par%mem%h( par%mem%n) = par%mem%total
    END IF
    
  END SUBROUTINE create_remapping_arrays_conservative
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

    DO aci = mesh%ci1, mesh%ci2
    
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
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
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
  SUBROUTINE integrate_over_Voronoi_boundaries( mesh_top, mesh_bot, nV_max, nS_max, proc_domain_Ac_top, &
    vi_top_start_proc, Vor_lines_top, Vor_lines_bot, Vor_vi_ti_bot, r_Ac_proc, CountCoincidences)
    ! Integrate over all the Voronoi boundarie lines of mesh_top, following them through mesh_bot.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_top
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_bot
    INTEGER,                                 INTENT(IN)    :: nV_max
    INTEGER,                                 INTENT(IN)    :: nS_max
    INTEGER,  DIMENSION(:    ),              INTENT(IN)    :: proc_domain_Ac_top
    INTEGER,                                 INTENT(IN)    :: vi_top_start_proc
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_top
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_bot
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: Vor_vi_ti_bot
    TYPE(type_remapping_conservative_intermediate_Ac_local), INTENT(INOUT) :: r_Ac_proc
    LOGICAL,                                 INTENT(IN)    :: CountCoincidences
    
    ! Local variables:
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                :: Aci_map_top
    INTEGER                                                :: vi_top, vi_bot, ci, aci, vc
    INTEGER                                                :: r_sng_nS                 ! How many Voronoi cells of the opposite mesh does this line pass through
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                :: r_sng_vi_left            ! Which vertices of the opposite mesh lie to the left  of this line
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                :: r_sng_vi_right           ! Which vertices of the opposite mesh lie to the right of this line
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: r_sng_LI_xdy             ! The three line integrals
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: r_sng_LI_mxydx
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: r_sng_LI_xydy
    INTEGER                                                :: sli1, sli2, si, sli, n_max_new

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
          r_sng_nS, r_sng_vi_left, r_sng_vi_right, r_sng_LI_xdy, r_sng_LI_mxydx, r_sng_LI_xydy, vi_bot, CountCoincidences, nV_max)

        ! Store data
        sli1 = r_Ac_proc%n_tot + 1
        r_Ac_proc%n_tot = r_Ac_proc%n_tot + r_sng_nS
        sli2 = r_Ac_proc%n_tot
        r_Ac_proc%nS(   aci) = r_sng_nS
        r_Ac_proc%sli1( aci) = sli1
        r_Ac_proc%sli2( aci) = sli2
        DO si = 1, r_sng_nS
          sli = sli1 + si - 1
          r_Ac_proc%vi_opp_left(  sli) = r_sng_vi_left(  si)
          r_Ac_proc%vi_opp_right( sli) = r_sng_vi_right( si)
          r_Ac_proc%LI_xdy(       sli) = r_sng_LI_xdy(   si)
          r_Ac_proc%LI_mxydx(     sli) = r_sng_LI_mxydx( si)
          r_Ac_proc%LI_xydy(      sli) = r_sng_LI_xydy(  si)
        END DO
    
        ! Extend memory if necessary
        IF (r_Ac_proc%n_tot > r_Ac_proc%n_max - 20 * nS_max) THEN
          n_max_new = r_Ac_proc%n_tot + 100 * nS_max
          CALL extend_memory_Ac_local( r_Ac_proc, n_max_new)
        END IF
        
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
  SUBROUTINE reorder_contributions_from_lines_to_vertices( mesh_top, mesh_bot, nV_max, r_top_Ac, r_top_proc)
    ! Reorder the mesh_top integration data from an Ac-based list to a vertex-based list
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                                          INTENT(IN)    :: mesh_top
    TYPE(type_mesh),                                          INTENT(IN)    :: mesh_bot
    INTEGER,                                                  INTENT(IN)    :: nV_max
    TYPE(type_remapping_conservative_intermediate_Ac_shared), INTENT(INOUT) :: r_top_Ac
    TYPE(type_remapping_conservative_intermediate_local),     INTENT(INOUT) :: r_top_proc
    
    ! Local variables:
    INTEGER                                                :: n_max, n_max_new
    INTEGER                                                :: vi_top, ci_top, aci_top, sli_top, vi_bot
    REAL(dp)                                               :: LI_xdy, LI_mxydx, LI_xydy
    
    ! Allocate process-local memory to hold vertex-based lists
    n_max = mesh_top%nV + mesh_bot%nV
    CALL allocate_memory_local( r_top_proc, n_max, mesh_top%nV)
    
    ! Fill the local lists
    DO vi_top = mesh_top%vi1, mesh_top%vi2
      DO ci_top = 1, mesh_top%nC( vi_top)

        aci_top = mesh_top%iAci( vi_top, ci_top)

        DO sli_top = r_top_Ac%sli1( aci_top), r_top_Ac%sli2( aci_top)

          ! Make sure we're integrating in the right direction
          IF     (mesh_top%Aci( aci_top,1) == vi_top) THEN
            vi_bot   =  r_top_Ac%vi_opp_left(  sli_top)
            LI_xdy   =  r_top_Ac%LI_xdy(       sli_top)
            LI_mxydx =  r_top_Ac%LI_mxydx(     sli_top)
            LI_xydy  =  r_top_Ac%LI_xydy(      sli_top)
          ELSEIF (mesh_top%Aci( aci_top,2) == vi_top) THEN
            vi_bot   =  r_top_Ac%vi_opp_right( sli_top)
            LI_xdy   = -r_top_Ac%LI_xdy(       sli_top)
            LI_mxydx = -r_top_Ac%LI_mxydx(     sli_top)
            LI_xydy  = -r_top_Ac%LI_xydy(      sli_top)
          ELSE
            WRITE(0,*) '  reorder_contributions_from_lines_to_vertices ERROR - serious error in Aci indexing!'
            STOP
          END IF

          CALL add_entry_to_local_vertex_based_list( r_top_proc, vi_top, vi_bot, LI_xdy, LI_mxydx, LI_xydy)

        END DO ! for si = r_top_Ac%si1( aci_top): r_top_Ac%si2( aci_top)

      END DO ! for ci_top = 1: mesh_top%nC( vi_top)

      ! Extend memory if necessary
      IF (r_top_proc%n_tot > r_top_proc%n_max - 2 * nV_max) THEN
        n_max_new = r_top_proc%n_tot + 10 * nV_max
        CALL extend_memory_local( r_top_proc, n_max_new)
      END IF

    END DO ! for vi_top = 1: mesh_top%nV
    CALL sync
    
    ! Deallocate Ac-based data
    CALL deallocate_shared( r_top_Ac%wn_max)
    CALL deallocate_shared( r_top_Ac%wn_tot)
    CALL deallocate_shared( r_top_Ac%wnS)
    CALL deallocate_shared( r_top_Ac%wsli1)
    CALL deallocate_shared( r_top_Ac%wsli2)
    CALL deallocate_shared( r_top_Ac%wvi_opp_left)
    CALL deallocate_shared( r_top_Ac%wvi_opp_right)
    CALL deallocate_shared( r_top_Ac%wLI_xdy)
    CALL deallocate_shared( r_top_Ac%wLI_mxydx)
    CALL deallocate_shared( r_top_Ac%wLI_xydy)
    
  END SUBROUTINE reorder_contributions_from_lines_to_vertices
  SUBROUTINE add_entry_to_local_vertex_based_list( r_top, vi_top, vi_bot, LI_xdy, LI_mxydx, LI_xydy)
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_remapping_conservative_intermediate_local), INTENT(INOUT) :: r_top
    INTEGER,                                              INTENT(IN)    :: vi_top
    INTEGER,                                              INTENT(IN)    :: vi_bot
    REAL(dp),                                             INTENT(IN)    :: LI_xdy
    REAL(dp),                                             INTENT(IN)    :: LI_mxydx
    REAL(dp),                                             INTENT(IN)    :: LI_xydy
    
    ! Local variables:
    LOGICAL                                                             :: IsListed
    INTEGER                                                             :: vli_top, vli_top2

    ! Check if this opposite mesh vertex is already listed
    IsListed = .FALSE.
    vli_top      = 0
    IF (r_top%nV( vi_top) > 0) THEN
      DO vli_top2 = r_top%vli1( vi_top), r_top%vli2( vi_top)
        IF (r_top%vi_opp( vli_top2) == vi_bot) THEN
          IsListed = .TRUE.
          vli_top = vli_top2
          EXIT
        END IF
      END DO
    END IF

    IF (.NOT. IsListed) THEN
      r_top%n_tot = r_top%n_tot + 1
      IF (r_top%nV( vi_top) == 0) THEN
        r_top%nV(   vi_top) = 1
        r_top%vli1( vi_top) = r_top%n_tot
        r_top%vli2( vi_top) = r_top%n_tot
      ELSE
        r_top%nV(   vi_top) = r_top%nV(   vi_top) + 1
        r_top%vli2( vi_top) = r_top%vli2( vi_top) + 1
      END IF
      r_top%vi_opp(   r_top%n_tot) = vi_bot
      r_top%LI_xdy(   r_top%n_tot) = LI_xdy
      r_top%LI_mxydx( r_top%n_tot) = LI_mxydx
      r_top%LI_xydy(  r_top%n_tot) = LI_xydy
    ELSE
      r_top%LI_xdy(   vli_top) = r_top%LI_xdy(   vli_top) + LI_xdy
      r_top%LI_mxydx( vli_top) = r_top%LI_mxydx( vli_top) + LI_mxydx
      r_top%LI_xydy(  vli_top) = r_top%LI_xydy(  vli_top) + LI_xydy
    END IF

  END SUBROUTINE add_entry_to_local_vertex_based_list
  SUBROUTINE add_entries_from_opposite_mesh( mesh_top, mesh_bot, r_top, r_bot)
    ! Add data from r_bot to r_top
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                                          INTENT(IN)    :: mesh_top
    TYPE(type_mesh),                                          INTENT(IN)    :: mesh_bot
    TYPE(type_remapping_conservative_intermediate_shared),    INTENT(INOUT) :: r_top
    TYPE(type_remapping_conservative_intermediate_shared),    INTENT(INOUT) :: r_bot
    
    ! Local variables:
    TYPE(type_remapping_conservative_intermediate_local)                    :: r_temp
    INTEGER                                                                 :: vi_top, vli_top, vi_bot, vli_bot, vi_top_opp
    REAL(dp)                                                                :: LI_xdy, LI_mxydx, LI_xydy
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                                 :: nV_extra
    INTEGER                                                                 :: vli1_from, vli2_from, vli1_to, vli2_to
    
    IF (par%master) THEN
      ALLOCATE( nV_extra( mesh_top%nV))
    ELSE
      ALLOCATE( nV_extra(1))
    END IF
    nV_extra = 0
    
    IF (par%master) THEN
    
      ! First go over all entries in r_top. Find the contributing vi_top
      ! vertices in r_bot and add their line integrals to r_top.
      ! Remove those entries from r_bot.
      ! ===================================================================

      DO vi_top = 1, mesh_top%nV
        DO vli_top = r_top%vli1( vi_top), r_top%vli2( vi_top)

          vi_bot = r_top%vi_opp( vli_top)

          vli_bot = r_bot%vli1( vi_bot)
          DO WHILE (vli_bot <= r_bot%vli2( vi_bot))
            vi_top_opp = r_bot%vi_opp( vli_bot)
            IF (vi_top_opp == vi_top) THEN
            
              ! Add the line integrals to r_top
              r_top%LI_xdy(   vli_top) = r_top%LI_xdy(   vli_top) + r_bot%LI_xdy(   vli_bot)
              r_top%LI_mxydx( vli_top) = r_top%LI_mxydx( vli_top) + r_bot%LI_mxydx( vli_bot)
              r_top%LI_xydy(  vli_top) = r_top%LI_xydy(  vli_top) + r_bot%LI_xydy(  vli_bot)
              
              ! Remove this entry from r_bot
              r_bot%vi_opp(   r_bot%vli1( vi_bot): r_bot%vli2( vi_bot)) = [r_bot%vi_opp(   r_bot%vli1( vi_bot): vli_bot-1), r_bot%vi_opp(   vli_bot+1: r_bot%vli2( vi_bot)), 0    ]
              r_bot%LI_xdy(   r_bot%vli1( vi_bot): r_bot%vli2( vi_bot)) = [r_bot%LI_xdy(   r_bot%vli1( vi_bot): vli_bot-1), r_bot%LI_xdy(   vli_bot+1: r_bot%vli2( vi_bot)), 0._dp]
              r_bot%LI_mxydx( r_bot%vli1( vi_bot): r_bot%vli2( vi_bot)) = [r_bot%LI_mxydx( r_bot%vli1( vi_bot): vli_bot-1), r_bot%LI_mxydx( vli_bot+1: r_bot%vli2( vi_bot)), 0._dp]
              r_bot%LI_xydy ( r_bot%vli1( vi_bot): r_bot%vli2( vi_bot)) = [r_bot%LI_xydy(  r_bot%vli1( vi_bot): vli_bot-1), r_bot%LI_xydy(  vli_bot+1: r_bot%vli2( vi_bot)), 0._dp]
              r_bot%n_tot = r_bot%n_tot - 1
              r_bot%nV( vi_bot) = r_bot%nV( vi_bot) - 1
              r_bot%vli2( vi_bot) = r_bot%vli2( vi_bot) - 1
              IF (r_bot%vli2( vi_bot) == 0) r_bot%vli1( vi_bot) = 0
              
            ELSE
              vli_bot = vli_bot + 1
            END IF
          END DO

        END DO ! DO vli_top = r_top%vli1( vi_top), r_top%vli2( vi_top)
      END DO ! DO vi_top = 1, mesh_top%nV

      ! The remaining entries in r_bot have no corresponding entries in r_top.
      ! This can only happen when a mesh_bot Voronoi cell is completely enclosed
      ! within a mesh_top Voronoi cell. These line integrals will be added as
      ! separate entries to r_top by combining the two lists.
      ! ========================================================================

      ! First determine how many extra entries each vi_top will need
      DO vi_bot = 1, mesh_bot%nV
        IF (r_bot%nV( vi_bot) == 0) THEN
          CYCLE
        ELSEIF (r_bot%nV( vi_bot) > 1) THEN
          WRITE(0,*) 'add_entries_from_opposite_mesh - ERROR: reduced r_bot%nV should never be > 1!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSE
          vi_top = r_bot%vi_opp( r_bot%vli1( vi_bot))
          nV_extra( vi_top) = nV_extra( vi_top) + 1
        END IF
      END DO

      ! Allocate temporary memory for storing unmerged r_top
      CALL allocate_memory_local( r_temp, r_top%n_tot, mesh_top%nV)

      ! Copy data there
      r_temp%n_max    = r_top%n_tot
      r_temp%n_tot    = r_top%n_tot
      r_temp%nV       = r_top%nV
      r_temp%vli1     = r_top%vli1
      r_temp%vli2     = r_top%vli2
      r_temp%vi_opp   = r_top%vi_opp(   1:r_top%n_tot)
      r_temp%LI_xdy   = r_top%LI_xdy(   1:r_top%n_tot)
      r_temp%LI_mxydx = r_top%LI_mxydx( 1:r_top%n_tot)
      r_temp%LI_xydy  = r_top%LI_xydy(  1:r_top%n_tot)
      
    END IF ! IF (par%master) THEN
    CALL sync

    ! Deallocate r_top
    r_top%nV(   mesh_top%vi1:mesh_top%vi2) = 0
    r_top%vli1( mesh_top%vi1:mesh_top%vi2) = 0
    r_top%vli2( mesh_top%vi1:mesh_top%vi2) = 0
    CALL deallocate_shared( r_top%wvi_opp)
    CALL deallocate_shared( r_top%wLI_xdy)
    CALL deallocate_shared( r_top%wLI_mxydx)
    CALL deallocate_shared( r_top%wLI_xydy)

    ! Allocate new memory that can accomodate both lists
    IF (par%master) r_top%n_max = r_temp%n_tot + r_bot%n_tot
    IF (par%master) r_top%n_tot = 0
    CALL sync
    CALL allocate_shared_int_1D( r_top%n_max, r_top%vi_opp,   r_top%wvi_opp)
    CALL allocate_shared_dp_1D(  r_top%n_max, r_top%LI_xdy,   r_top%wLI_xdy)
    CALL allocate_shared_dp_1D(  r_top%n_max, r_top%LI_mxydx, r_top%wLI_mxydx)
    CALL allocate_shared_dp_1D(  r_top%n_max, r_top%LI_xydy,  r_top%wLI_xydy)
    
    IF (par%master) THEN

      ! Copy data back from temporary memory (adjust list indices accordingly)
      vli1_to = 1
      DO vi_top = 1, mesh_top%nV

        IF (r_temp%nV( vi_top) > 0) THEN

          ! List indices
          vli1_from = r_temp%vli1( vi_top)
          vli2_from = r_temp%vli2( vi_top)
          vli2_to   = vli1_to + r_temp%nV( vi_top) - 1

          r_top%n_tot = r_top%n_tot + r_temp%nV( vi_top)
          r_top%nV(   vi_top) = r_temp%nV( vi_top)
          r_top%vli1( vi_top) = vli1_to
          r_top%vli2( vi_top) = vli2_to

          ! Copy data
          r_top%vi_opp(   vli1_to: vli2_to) = r_temp%vi_opp(   vli1_from: vli2_from)
          r_top%LI_xdy(   vli1_to: vli2_to) = r_temp%LI_xdy(   vli1_from: vli2_from)
          r_top%LI_mxydx( vli1_to: vli2_to) = r_temp%LI_mxydx( vli1_from: vli2_from)
          r_top%LI_xydy(  vli1_to: vli2_to) = r_temp%LI_xydy(  vli1_from: vli2_from)

          vli1_to = vli1_to + r_top%nV( vi_top) + nV_extra( vi_top)
        END IF

      END DO ! DO vi_top = 1, mesh_top%nV

      ! Deallocate temporary memory
      r_temp%n_max     = 0
      r_temp%n_tot     = 0
      DEALLOCATE( r_temp%nV)
      DEALLOCATE( r_temp%vli1)
      DEALLOCATE( r_temp%vli2)
      DEALLOCATE( r_temp%vi_opp)
      DEALLOCATE( r_temp%LI_xdy)
      DEALLOCATE( r_temp%LI_mxydx)
      DEALLOCATE( r_temp%LI_xydy)

      ! Add entries from r_bot
      DO vi_bot = 1, mesh_bot%nV
        IF (r_bot%nV( vi_bot) == 1) THEN
          ! Find data
          vli_bot  = r_bot%vli1( vi_bot)
          vi_top   = r_bot%vi_opp( vli_bot)
          LI_xdy   = r_bot%LI_xdy( vli_bot)
          LI_mxydx = r_bot%LI_mxydx( vli_bot)
          LI_xydy  = r_bot%LI_xydy( vli_bot)
          ! Add to r_top (no need to check IF its listed, we already know its not)
          r_top%n_tot = r_top%n_tot + 1
          r_top%nV(   vi_top) = r_top%nV(   vi_top) + 1
          r_top%vli2( vi_top) = r_top%vli2( vi_top) + 1
          vli_top = r_top%vli2( vi_top)
          r_top%vi_opp(   vli_top) = vi_bot
          r_top%LI_xdy(   vli_top) = LI_xdy
          r_top%LI_mxydx( vli_top) = LI_mxydx
          r_top%LI_xydy(  vli_top) = LI_xydy
        END IF
      END DO ! DO vi_bot = 1, mesh_bot%nV
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Clean up after yourself
    DEALLOCATE( nV_extra)
    
    ! Deallocate r_src
    CALL deallocate_shared( r_bot%wn_max)
    CALL deallocate_shared( r_bot%wn_tot)
    CALL deallocate_shared( r_bot%wnV)
    CALL deallocate_shared( r_bot%wvli1)
    CALL deallocate_shared( r_bot%wvli2)
    CALL deallocate_shared( r_bot%wvi_opp)
    CALL deallocate_shared( r_bot%wLI_xdy)
    CALL deallocate_shared( r_bot%wLI_mxydx)
    CALL deallocate_shared( r_bot%wLI_xydy)

  END SUBROUTINE add_entries_from_opposite_mesh
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
    INTEGER                                                :: ci, vc, aci, vii, vii2, vi_src2, it
    REAL(dp), DIMENSION(2)                                 :: p, pc, ccr, ccl, p1, p_next_top, p_next_bot
    REAL(dp)                                               :: LI_xdy, LI_mxydx, LI_xydy
    LOGICAL                                                :: Finished, IsListed
    
  ! Southern boundary
  ! ==================

    vi_bot = 1
    vi_top = 1
    p = [mesh_bot%xmin, mesh_bot%ymin]

    Finished = .FALSE.
    it = 0
    DO WHILE (.NOT.Finished)
    
      it = it + 1
      IF (it > mesh_top%nV + mesh_bot%nV) THEN
        WRITE(0,*) '    Remapping - integrate around domain boundary got stuck (south)!'
        STOP
      END IF

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
      CALL line_integral_xdy(   p1, p, mesh_top%tol_dist, LI_xdy  )
      CALL line_integral_mxydx( p1, p, mesh_top%tol_dist, LI_mxydx)
      CALL line_integral_xydy(  p1, p, mesh_top%tol_dist, LI_xydy )

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
    it = 0
    DO WHILE (.NOT.Finished)
    
      it = it + 1
      IF (it > mesh_top%nV + mesh_bot%nV) THEN
        WRITE(0,*) '    Remapping - integrate around domain boundary got stuck (east)!'
        STOP
      END IF

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
      CALL line_integral_xdy(   p1, p, mesh_top%tol_dist, LI_xdy  )
      CALL line_integral_mxydx( p1, p, mesh_top%tol_dist, LI_mxydx)
      CALL line_integral_xydy(  p1, p, mesh_top%tol_dist, LI_xydy )

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
    it = 0
    DO WHILE (.NOT.Finished)
    
      it = it + 1
      IF (it > mesh_top%nV + mesh_bot%nV) THEN
        WRITE(0,*) '    Remapping - integrate around domain boundary got stuck (north)!'
        STOP
      END IF

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
      CALL line_integral_xdy(   p1, p, mesh_top%tol_dist, LI_xdy  )
      CALL line_integral_mxydx( p1, p, mesh_top%tol_dist, LI_mxydx)
      CALL line_integral_xydy(  p1, p, mesh_top%tol_dist, LI_xydy )

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
    it = 0
    DO WHILE (.NOT.Finished)
    
      it = it + 1
      IF (it > mesh_top%nV + mesh_bot%nV) THEN
        WRITE(0,*) '    Remapping - integrate around domain boundary got stuck (west)!'
        STOP
      END IF

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
      CALL line_integral_xdy(   p1, p, mesh_top%tol_dist, LI_xdy  )
      CALL line_integral_mxydx( p1, p, mesh_top%tol_dist, LI_mxydx)
      CALL line_integral_xydy(  p1, p, mesh_top%tol_dist, LI_xydy )

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
  SUBROUTINE calculate_remapping_weights_from_line_integrals( mesh_top, mesh_bot, r_top, map)
    ! Calculate the actual remapping weights from the three line integrals
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                                          INTENT(IN)    :: mesh_top
    TYPE(type_mesh),                                          INTENT(IN)    :: mesh_bot
    TYPE(type_remapping_conservative_intermediate_shared),    INTENT(INOUT) :: r_top
    TYPE(type_remapping_conservative),                        INTENT(INOUT) :: map
    
    ! Local variables:
    INTEGER                                                :: vi_top, vli_top, vi_bot
    
    ! Allocate shared memory
    CALL allocate_shared_int_0D(              map%n_tot, map%wn_tot)
    CALL allocate_shared_int_1D( mesh_top%nV, map%vli1,  map%wvli1 )
    CALL allocate_shared_int_1D( mesh_top%nV, map%vli2,  map%wvli2 )
    CALL allocate_shared_int_1D( r_top%n_tot, map%vi,    map%wvi   )
    CALL allocate_shared_dp_1D(  r_top%n_tot, map%w0,    map%ww0   )
    CALL allocate_shared_dp_1D(  r_top%n_tot, map%w1x,   map%ww1x  )
    CALL allocate_shared_dp_1D(  r_top%n_tot, map%w1y,   map%ww1y  )
    
    IF (par%master) map%n_tot = r_top%n_tot
    
    DO vi_top = mesh_top%vi1, mesh_top%vi2
    
      map%vli1( vi_top) = r_top%vli1( vi_top)
      map%vli2( vi_top) = r_top%vli2( vi_top)
      
      IF (mesh_top%edge_index( vi_top) > 0) THEN
        ! Exception: since the current remapping code does not integrate over the domain boundary,
        ! remapping weights for edge vertices are really far off. Since boundary conditions apply
        
        DO vli_top = map%vli1( vi_top), map%vli2( vi_top)
          ! there anyway, just set the weights to zero.DO vli_top = map%vli1( vi_top), map%vli2( vi_top)
          map%vi(  vli_top) = r_top%vi_opp( vli_top)
          map%w0(  vli_top) = 0._dp
          map%w1x( vli_top) = 0._dp
          map%w1y( vli_top) = 0._dp
        END DO
        
      ELSE ! IF (mesh_top%edge_index( vi_top) > 0) THEN
      
        DO vli_top = map%vli1( vi_top), map%vli2( vi_top)
          vi_bot = r_top%vi_opp( vli_top)
          map%vi(  vli_top) = vi_bot
          map%w0(  vli_top) = (r_top%LI_xdy(   vli_top) / mesh_top%A( vi_top))
          map%w1x( vli_top) = (r_top%LI_mxydx( vli_top) / mesh_top%A( vi_top)) - (mesh_bot%VorGC( vi_bot,1) * map%w0( vli_top))
          map%w1y( vli_top) = (r_top%LI_xydy(  vli_top) / mesh_top%A( vi_top)) - (mesh_bot%VorGC( vi_bot,2) * map%w0( vli_top))
        END DO
      
      END IF ! IF (mesh_top%edge_index( vi_top) > 0) THEN
      
    END DO !  DO vi_top = mesh_top%vi1, mesh_top%vi2
    
    ! Deallocate r_dst
    CALL deallocate_shared( r_top%wn_max   )
    CALL deallocate_shared( r_top%wn_tot   )
    CALL deallocate_shared( r_top%wnV      )
    CALL deallocate_shared( r_top%wvli1    )
    CALL deallocate_shared( r_top%wvli2    )
    CALL deallocate_shared( r_top%wvi_opp  )
    CALL deallocate_shared( r_top%wLI_xdy  )
    CALL deallocate_shared( r_top%wLI_mxydx)
    CALL deallocate_shared( r_top%wLI_xydy )
    
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
    
    d_src( mesh_src%vi1:mesh_src%vi2) = 1._dp
    CALL remap_cons_2nd_order_2D( mesh_src, mesh_dst, map, d_src, d_dst)
    
    ! Ignore edge vertices, as they're useless, and maybe 99% of remapping errors occur there
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
      IF (mesh_dst%edge_index( vi_dst) > 0) d_dst( vi_dst) = 1._dp
    END DO
    CALL sync
    
    ! Check for vertices where conservative remapping went wrong. If there are not too many, replace them with linear interpolation.
    n_wrong = 0
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
      IF (ABS(1._dp - d_dst( vi_dst)) > 1E-1_dp) THEN
      
        n_wrong = n_wrong + 1
        
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
        CALL find_triangle_area( pa, pb, pc, Atot)
        CALL find_triangle_area( pb, pc, p , Aa  )
        CALL find_triangle_area( pc, pa, p , Ab  )
        CALL find_triangle_area( pa, pb, p , Ac  )
        
        ! Add new weights to the list
        IF     (map%vli2( vi_dst) == map%vli1( vi_dst) + 1) THEN
          ! Conservative remapping provided only 2 contributing src vertices
        
          map%vi(  map%vli1( vi_dst):map%vli2( vi_dst)) = 0
          map%w0(  map%vli1( vi_dst):map%vli2( vi_dst)) = 0._dp
          map%w1x( map%vli1( vi_dst):map%vli2( vi_dst)) = 0._dp
          map%w1y( map%vli1( vi_dst):map%vli2( vi_dst)) = 0._dp
          
          map%vli2( vi_dst) = map%vli1( vi_dst) + 1
          
          IF     (Aa <= Ab .AND. Aa <= Ac) THEN
            ! A is the smallest; only list B and C
            map%vi( map%vli1( vi_dst):map%vli2( vi_dst)) = [vib_src, vic_src]
            map%w0( map%vli1( vi_dst):map%vli2( vi_dst)) = [Ab/(Ab+Ac), Ac/(Ab+Ac)] 
          ELSEIF (Ab <= Aa .AND. Ab <= Ac) THEN
            ! B is the smallest; only list A and C
            map%vi( map%vli1( vi_dst):map%vli2( vi_dst)) = [via_src, vic_src]
            map%w0( map%vli1( vi_dst):map%vli2( vi_dst)) = [Aa/(Aa+Ac), Ac/(Aa+Ac)] 
          ELSE
            ! C is the smallest; only list A and B
            map%vi( map%vli1( vi_dst):map%vli2( vi_dst)) = [via_src, vib_src]
            map%w0( map%vli1( vi_dst):map%vli2( vi_dst)) = [Aa/(Aa+Ab), Aa/(Aa+Ab)] 
          END IF
          
        ELSEIF (map%vli2( vi_dst) == map%vli1( vi_dst) ) THEN
          ! Conservative remapping provided only 1 contributing src vertex
        
          map%vi(  map%vli1( vi_dst)) = 0
          map%w0(  map%vli1( vi_dst)) = 0._dp
          map%w1x( map%vli1( vi_dst)) = 0._dp
          map%w1y( map%vli1( vi_dst)) = 0._dp
          
          map%vli2( vi_dst) = map%vli1( vi_dst) 
          
          IF     (Aa >= Ab .AND. Aa >= Ac) THEN
            ! A is the largest; list only that one
            map%vi( map%vli1( vi_dst)) = via_src
            map%w0( map%vli1( vi_dst)) = 1._dp 
          ELSEIF (Ab >= Aa .AND. Ab >= Ac) THEN
            ! B is the largest; list only that one
            map%vi( map%vli1( vi_dst)) = vib_src
            map%w0( map%vli1( vi_dst)) = 1._dp 
          ELSE
            ! C is the largest; list only that one
            map%vi( map%vli1( vi_dst)) = vic_src
            map%w0( map%vli1( vi_dst)) = 1._dp 
          END IF
          
        ELSE
          ! Conservative remapping provided at least 3 contributing src vertices
        
          map%vi(  map%vli1( vi_dst):map%vli2( vi_dst)) = 0
          map%w0(  map%vli1( vi_dst):map%vli2( vi_dst)) = 0._dp
          map%w1x( map%vli1( vi_dst):map%vli2( vi_dst)) = 0._dp
          map%w1y( map%vli1( vi_dst):map%vli2( vi_dst)) = 0._dp
          
          map%vli2( vi_dst) = map%vli1( vi_dst) + 2
          
          map%vi( map%vli1( vi_dst):map%vli2( vi_dst)) = [via_src, vib_src, vic_src]
          map%w0( map%vli1( vi_dst):map%vli2( vi_dst)) = [Aa/Atot, Ab/Atot, Ac/Atot] 
          
        END IF ! IF     (map%vli2( vi_dst) == map%vli1( vi_dst) + 1) THEN
                
      END IF ! IF (ABS(1._dp - d_dst( vi_dst)) > 1E-1_dp) THEN
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
  
! == The line-tracing algorithm used in conservative remapping
  SUBROUTINE calculate_line_integral_contributions( mesh_bot, Vor_lines_top, Vor_lines_bot, Vor_vi_ti_bot, aci, &
          r_sng_nS, r_sng_vi_left, r_sng_vi_right, r_sng_LI_xdy, r_sng_LI_mxydx, r_sng_LI_xydy, vi_bot, CountCoincidences, nV_max)
    ! Trace Voronoi boundary line aci through mesh_bot, calculate the three line integrals along all sections of this line.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                         INTENT(INOUT) :: mesh_bot
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_top
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: Vor_lines_bot
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: Vor_vi_ti_bot
    INTEGER,                                 INTENT(IN)    :: aci
    INTEGER,                                 INTENT(OUT)   :: r_sng_nS
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
    r_sng_nS       = 0
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

      ! If applicable, calculate the three line integrals over the line section p-p_next and add them to the lists
      DoAdd = .TRUE.
      IF (Coincides .AND. (.NOT. CountCoincidences)) DoAdd = .FALSE.
      IF (DoAdd) THEN
    
        r_sng_nS = r_sng_nS + 1
        
        IF (r_sng_nS > nV_max) THEN
          WRITE(0,*) '  Conservative remapping exceeded memory of r_sng; nS_max is too small!'
          STOP
        END IF
        
        r_sng_vi_left(  r_sng_nS) = vi_bot_left
        r_sng_vi_right( r_sng_nS) = vi_bot_right
        CALL line_integral_xdy(   p, p_next, mesh_bot%tol_dist, r_sng_LI_xdy(   r_sng_nS))
        CALL line_integral_mxydx( p, p_next, mesh_bot%tol_dist, r_sng_LI_mxydx( r_sng_nS))
        CALL line_integral_xydy(  p, p_next, mesh_bot%tol_dist, r_sng_LI_xydy(  r_sng_nS))
    
      END IF ! IF (DoAdd) THEN

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
  
! == Manipulating intermediate forms of data lists in conservative remapping
  SUBROUTINE allocate_memory_Ac_local( r_Ac_proc, n_max, nAc)
    ! Allocate process-local memory for intermediate Ac integration data
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_remapping_conservative_intermediate_Ac_local), INTENT(INOUT) :: r_Ac_proc
    INTEGER,                                                 INTENT(IN)    :: n_max
    INTEGER,                                                 INTENT(IN)    :: nAc
    
    r_Ac_proc%n_max = n_max
    r_Ac_proc%n_tot = 0
    ALLOCATE( r_Ac_proc%nS(           nAc))
    ALLOCATE( r_Ac_proc%sli1(         nAc))
    ALLOCATE( r_Ac_proc%sli2(         nAc))
    ALLOCATE( r_Ac_proc%vi_opp_left(  r_Ac_proc%n_max))
    ALLOCATE( r_Ac_proc%vi_opp_right( r_Ac_proc%n_max))
    ALLOCATE( r_Ac_proc%LI_xdy(       r_Ac_proc%n_max))
    ALLOCATE( r_Ac_proc%LI_mxydx(     r_Ac_proc%n_max))
    ALLOCATE( r_Ac_proc%LI_xydy(      r_Ac_proc%n_max))
    
    r_Ac_proc%nS           = 0
    r_Ac_proc%sli1         = 0
    r_Ac_proc%sli2         = -1
    r_Ac_proc%vi_opp_left  = 0
    r_Ac_proc%vi_opp_right = 0
    r_Ac_proc%LI_xdy       = 0._dp
    r_Ac_proc%LI_mxydx     = 0._dp
    r_Ac_proc%LI_xydy      = 0._dp
    
  END SUBROUTINE allocate_memory_Ac_local
  SUBROUTINE extend_memory_Ac_local( r_Ac_proc, n_max_new)
    ! Extend process-local memory for intermediate Ac integration data
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_remapping_conservative_intermediate_Ac_local), INTENT(INOUT) :: r_Ac_proc
    INTEGER,                                                 INTENT(IN)    :: n_max_new
    
    ! Local variables:
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: d_temp_int
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: d_temp_dp

    r_Ac_proc%n_max = n_max_new

    ! One at a time to minimise memory spikes
    
    ALLOCATE( d_temp_int( r_Ac_proc%n_tot))
    
    ! vi_opp_left
    d_temp_int = r_Ac_proc%vi_opp_left( 1:r_Ac_proc%n_tot)
    DEALLOCATE( r_Ac_proc%vi_opp_left)
    ALLOCATE( r_Ac_proc%vi_opp_left( r_Ac_proc%n_max))
    r_Ac_proc%vi_opp_left( 1:r_Ac_proc%n_tot) = d_temp_int
    r_Ac_proc%vi_opp_left( r_Ac_proc%n_tot+1:r_Ac_proc%n_max) = 0

    ! vi_opp_right
    d_temp_int = r_Ac_proc%vi_opp_right( 1:r_Ac_proc%n_tot)
    DEALLOCATE( r_Ac_proc%vi_opp_right)
    ALLOCATE( r_Ac_proc%vi_opp_right( r_Ac_proc%n_max))
    r_Ac_proc%vi_opp_right( 1:r_Ac_proc%n_tot) = d_temp_int
    r_Ac_proc%vi_opp_right( r_Ac_proc%n_tot+1:r_Ac_proc%n_max) = 0
    
    DEALLOCATE( d_temp_int)
    
    ALLOCATE( d_temp_dp( r_Ac_proc%n_tot))
    
    ! LI_xdy
    d_temp_dp = r_Ac_proc%LI_xdy( 1:r_Ac_proc%n_tot)
    DEALLOCATE( r_Ac_proc%LI_xdy)
    ALLOCATE( r_Ac_proc%LI_xdy( r_Ac_proc%n_max))
    r_Ac_proc%LI_xdy( 1:r_Ac_proc%n_tot) = d_temp_dp
    r_Ac_proc%LI_xdy( r_Ac_proc%n_tot+1:r_Ac_proc%n_max) = 0._dp
    
    ! LI_mxydx
    d_temp_dp = r_Ac_proc%LI_mxydx( 1:r_Ac_proc%n_tot)
    DEALLOCATE( r_Ac_proc%LI_mxydx)
    ALLOCATE( r_Ac_proc%LI_mxydx( r_Ac_proc%n_max))
    r_Ac_proc%LI_mxydx( 1:r_Ac_proc%n_tot) = d_temp_dp
    r_Ac_proc%LI_mxydx( r_Ac_proc%n_tot+1:r_Ac_proc%n_max) = 0._dp
    
    ! LI_xydy
    d_temp_dp = r_Ac_proc%LI_xydy( 1:r_Ac_proc%n_tot)
    DEALLOCATE( r_Ac_proc%LI_xydy)
    ALLOCATE( r_Ac_proc%LI_xydy( r_Ac_proc%n_max))
    r_Ac_proc%LI_xydy( 1:r_Ac_proc%n_tot) = d_temp_dp
    r_Ac_proc%LI_xydy( r_Ac_proc%n_tot+1:r_Ac_proc%n_max) = 0._dp
    
    DEALLOCATE( d_temp_dp)
    
  END SUBROUTINE extend_memory_Ac_local
  SUBROUTINE gather_memory_Ac( mesh, r_Ac_proc, r_Ac, proc_domain_Ac)
    ! Gather process-local intermediate Ac integration data into shared memory
    ! Deallocates the process-local memory
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                                          INTENT(IN)    :: mesh
    TYPE(type_remapping_conservative_intermediate_Ac_local),  INTENT(INOUT) :: r_Ac_proc
    TYPE(type_remapping_conservative_intermediate_Ac_shared), INTENT(INOUT) :: r_Ac
    INTEGER,  DIMENSION(:    ),                               INTENT(IN)    :: proc_domain_Ac
    
    ! Local variables:
    INTEGER                                                :: dsli, aci, p, n_tot_proc, status(MPI_STATUS_SIZE)
    
    ! Determine list index offset for each process
    IF (par%master) dsli = r_Ac_proc%n_tot
    DO p = 1, par%n-1
      IF (par%master) THEN
        CALL MPI_SEND( dsli,            1, MPI_INTEGER, p, 1337, MPI_COMM_WORLD,         ierr)
        CALL MPI_RECV( n_tot_proc,      1, MPI_INTEGER, p, 1338, MPI_COMM_WORLD, status, ierr)
        dsli = dsli + n_tot_proc
      ELSEIF (par%i == p) THEN
        CALL MPI_RECV( dsli,            1, MPI_INTEGER, 0, 1337, MPI_COMM_WORLD, status, ierr)
        CALL MPI_SEND( r_Ac_proc%n_tot, 1, MPI_INTEGER, 0, 1338, MPI_COMM_WORLD,         ierr)
      END IF
      CALL sync
    END DO
    IF (par%master) dsli = 0
    
    ! Determine how much shared memory to allocate for the lists
    CALL allocate_shared_int_0D(      r_Ac%n_tot, r_Ac%wn_tot)
    CALL allocate_shared_int_0D(      r_Ac%n_max, r_Ac%wn_max)
    CALL MPI_REDUCE( r_Ac_proc%n_tot, r_Ac%n_tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( r_Ac_proc%n_tot, r_Ac%n_max, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    ! nS
    CALL allocate_shared_int_1D( mesh%nAc, r_Ac%nS,    r_Ac%wnS   )
    DO aci = 1, mesh%nAc
      IF (proc_domain_Ac( aci) /= par%i) CYCLE
      r_Ac%nS( aci) = r_Ac_proc%nS( aci)
    END DO
    DEALLOCATE( r_Ac_proc%nS)
    
    ! sli1
    CALL allocate_shared_int_1D( mesh%nAc, r_Ac%sli1,  r_Ac%wsli1 )
    DO aci = 1, mesh%nAc
      IF (proc_domain_Ac( aci) /= par%i) CYCLE
      r_Ac%sli1( aci) = r_Ac_proc%sli1( aci) + dsli
    END DO
    DEALLOCATE( r_Ac_proc%sli1)
    
    ! sli2
    CALL allocate_shared_int_1D( mesh%nAc, r_Ac%sli2,  r_Ac%wsli2 )
    DO aci = 1, mesh%nAc
      IF (proc_domain_Ac( aci) /= par%i) CYCLE
      r_Ac%sli2( aci) = r_Ac_proc%sli2( aci) + dsli
    END DO
    DEALLOCATE( r_Ac_proc%sli2)
    
    ! vi_opp_left
    CALL allocate_shared_int_1D( r_Ac%n_max, r_Ac%vi_opp_left, r_Ac%wvi_opp_left)
    r_Ac%vi_opp_left( dsli + 1 : dsli + r_Ac_proc%n_tot) = r_Ac_proc%vi_opp_left( 1:r_Ac_proc%n_tot)
    DEALLOCATE( r_Ac_proc%vi_opp_left)
    
    ! vi_opp_right
    CALL allocate_shared_int_1D( r_Ac%n_max, r_Ac%vi_opp_right, r_Ac%wvi_opp_right)
    r_Ac%vi_opp_right( dsli + 1 : dsli + r_Ac_proc%n_tot) = r_Ac_proc%vi_opp_right( 1:r_Ac_proc%n_tot)
    DEALLOCATE( r_Ac_proc%vi_opp_right)
    
    ! LI_xdy
    CALL allocate_shared_dp_1D( r_Ac%n_max, r_Ac%LI_xdy, r_Ac%wLI_xdy)
    r_Ac%LI_xdy( dsli + 1 : dsli + r_Ac_proc%n_tot) = r_Ac_proc%LI_xdy( 1:r_Ac_proc%n_tot)
    DEALLOCATE( r_Ac_proc%LI_xdy)
    
    ! LI_mxydx
    CALL allocate_shared_dp_1D( r_Ac%n_max, r_Ac%LI_mxydx, r_Ac%wLI_mxydx)
    r_Ac%LI_mxydx( dsli + 1 : dsli + r_Ac_proc%n_tot) = r_Ac_proc%LI_mxydx( 1:r_Ac_proc%n_tot)
    DEALLOCATE( r_Ac_proc%LI_mxydx)
    
    ! LI_xydy
    CALL allocate_shared_dp_1D( r_Ac%n_max, r_Ac%LI_xydy, r_Ac%wLI_xydy)
    r_Ac%LI_xydy( dsli + 1 : dsli + r_Ac_proc%n_tot) = r_Ac_proc%LI_xydy( 1:r_Ac_proc%n_tot)
    DEALLOCATE( r_Ac_proc%LI_xydy)
    
  END SUBROUTINE gather_memory_Ac
  SUBROUTINE allocate_memory_local( r_proc, n_max, nV)
    ! Allocate process-local memory for intermediate integration data
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_remapping_conservative_intermediate_local),    INTENT(INOUT) :: r_proc
    INTEGER,                                                 INTENT(IN)    :: n_max
    INTEGER,                                                 INTENT(IN)    :: nV
    
    r_proc%n_max = n_max
    r_proc%n_tot = 0
    ALLOCATE( r_proc%nV(       nV))
    ALLOCATE( r_proc%vli1(     nV))
    ALLOCATE( r_proc%vli2(     nV))
    ALLOCATE( r_proc%vi_opp(   r_proc%n_max))
    ALLOCATE( r_proc%LI_xdy(   r_proc%n_max))
    ALLOCATE( r_proc%LI_mxydx( r_proc%n_max))
    ALLOCATE( r_proc%LI_xydy(  r_proc%n_max))
    
    r_proc%nV       = 0
    r_proc%vli1     = 0
    r_proc%vli2     = 0
    r_proc%vi_opp   = 0
    r_proc%LI_xdy   = 0._dp
    r_proc%LI_mxydx = 0._dp
    r_proc%LI_xydy  = 0._dp
    
  END SUBROUTINE allocate_memory_local
  SUBROUTINE extend_memory_local( r_proc, n_max_new)
    ! Extend process-local memory for intermediate integration data
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_remapping_conservative_intermediate_local), INTENT(INOUT) :: r_proc
    INTEGER,                                              INTENT(IN)    :: n_max_new
    
    ! Local variables:
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: d_temp_int
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: d_temp_dp

    r_proc%n_max = n_max_new

    ! One at a time to minimise memory spikes

    ALLOCATE( d_temp_int( r_proc%n_tot))
    
    ! vi_opp
    d_temp_int = r_proc%vi_opp( 1:r_proc%n_tot)
    DEALLOCATE( r_proc%vi_opp)
    ALLOCATE( r_proc%vi_opp( r_proc%n_max))
    r_proc%vi_opp( 1:r_proc%n_tot) = d_temp_int
    r_proc%vi_opp( r_proc%n_tot+1:r_proc%n_max) = 0
    
    DEALLOCATE( d_temp_int)
    
    ALLOCATE( d_temp_dp( r_proc%n_tot))
    
    ! LI_xdy
    d_temp_dp = r_proc%LI_xdy( 1:r_proc%n_tot)
    DEALLOCATE( r_proc%LI_xdy)
    ALLOCATE( r_proc%LI_xdy( r_proc%n_max))
    r_proc%LI_xdy( 1:r_proc%n_tot) = d_temp_dp
    r_proc%LI_xdy( r_proc%n_tot+1:r_proc%n_max) = 0._dp
    
    ! LI_mxydx
    d_temp_dp = r_proc%LI_mxydx( 1:r_proc%n_tot)
    DEALLOCATE( r_proc%LI_mxydx)
    ALLOCATE( r_proc%LI_mxydx( r_proc%n_max))
    r_proc%LI_mxydx( 1:r_proc%n_tot) = d_temp_dp
    r_proc%LI_mxydx( r_proc%n_tot+1:r_proc%n_max) = 0._dp
    
    ! LI_xydy
    d_temp_dp = r_proc%LI_xydy( 1:r_proc%n_tot)
    DEALLOCATE( r_proc%LI_xydy)
    ALLOCATE( r_proc%LI_xydy( r_proc%n_max))
    r_proc%LI_xydy( 1:r_proc%n_tot) = d_temp_dp
    r_proc%LI_xydy( r_proc%n_tot+1:r_proc%n_max) = 0._dp
    
    DEALLOCATE( d_temp_dp)
    
  END SUBROUTINE extend_memory_local
  SUBROUTINE gather_memory( mesh, r_proc, r)
    ! Gather process-local intermediate Ac integration data into shared memory
    ! Deallocates the process-local memory
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                                       INTENT(IN)    :: mesh
    TYPE(type_remapping_conservative_intermediate_local),  INTENT(INOUT) :: r_proc
    TYPE(type_remapping_conservative_intermediate_shared), INTENT(INOUT) :: r
    
    ! Local variables:
    INTEGER                                                :: dvli, p, n_tot_proc, status(MPI_STATUS_SIZE)
    
    ! Determine list index offset for each process
    IF (par%master) dvli = r_proc%n_tot
    DO p = 1, par%n-1
      IF (par%master) THEN
        CALL MPI_SEND( dvli,         1, MPI_INTEGER, p, 0,           MPI_COMM_WORLD,         ierr)
        CALL MPI_RECV( n_tot_proc,   1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        dvli = dvli + n_tot_proc
      ELSEIF (par%i == p) THEN
        CALL MPI_RECV( dvli,         1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        CALL MPI_SEND( r_proc%n_tot, 1, MPI_INTEGER, 0, 0,           MPI_COMM_WORLD,         ierr)
      END IF
      CALL sync
    END DO
    IF (par%master) dvli = 0
    
    CALL allocate_shared_int_0D( r%n_max, r%wn_max)
    CALL allocate_shared_int_0D( r%n_tot, r%wn_tot)
    CALL MPI_REDUCE( r_proc%n_tot, r%n_max, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( r_proc%n_tot, r%n_tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    ! Let all processes move their data (one variable at a time to keep memory spikes to a minimum)
    
    ! nV
    CALL allocate_shared_int_1D( mesh%nV, r%nV,    r%wnV   )
    r%nV( mesh%vi1:mesh%vi2) = r_proc%nV( mesh%vi1:mesh%vi2)
    DEALLOCATE( r_proc%nV)
    
    ! vli1
    CALL allocate_shared_int_1D( mesh%nV, r%vli1,  r%wvli1 )
    r%vli1( mesh%vi1:mesh%vi2) = r_proc%vli1( mesh%vi1:mesh%vi2) + dvli
    DEALLOCATE( r_proc%vli1)
    
    ! vli2
    CALL allocate_shared_int_1D( mesh%nV, r%vli2,  r%wvli2 )
    r%vli2( mesh%vi1:mesh%vi2) = r_proc%vli2( mesh%vi1:mesh%vi2) + dvli
    DEALLOCATE( r_proc%vli2)
    
    ! vi_opp
    CALL allocate_shared_int_1D( r%n_max, r%vi_opp, r%wvi_opp)
    r%vi_opp(   1 + dvli : r_proc%n_tot + dvli) = r_proc%vi_opp(   1:r_proc%n_tot)
    DEALLOCATE( r_proc%vi_opp)
    
    ! LI_xdy
    CALL allocate_shared_dp_1D( r%n_max, r%LI_xdy, r%wLI_xdy)
    r%LI_xdy(   1 + dvli : r_proc%n_tot + dvli) = r_proc%LI_xdy(   1:r_proc%n_tot)
    DEALLOCATE( r_proc%LI_xdy)
    
    ! LI_mxydx
    CALL allocate_shared_dp_1D( r%n_max, r%LI_mxydx, r%wLI_mxydx)
    r%LI_mxydx( 1 + dvli : r_proc%n_tot + dvli) = r_proc%LI_mxydx( 1:r_proc%n_tot)
    DEALLOCATE( r_proc%LI_mxydx)
    
    ! LI_xydy
    CALL allocate_shared_dp_1D( r%n_max, r%LI_xydy, r%wLI_xydy)
    r%LI_xydy(  1 + dvli : r_proc%n_tot + dvli) = r_proc%LI_xydy(  1:r_proc%n_tot)
    DEALLOCATE( r_proc%LI_xydy)
    
  END SUBROUTINE gather_memory
  
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
        
    d_dst(mesh_dst%vi1:mesh_dst%vi2) = 0._dp
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
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
        
    d_dst(mesh_dst%vi1:mesh_dst%vi2,:) = 0._dp
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
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
        
    d_dst(mesh_dst%vi1:mesh_dst%vi2,:) = 0._dp
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
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
        
    d_dst(mesh_dst%vi1:mesh_dst%vi2) = 0._dp
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
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
        
    d_dst(mesh_dst%vi1:mesh_dst%vi2,:) = 0._dp
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
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
    INTEGER                                                :: vi_dst, vli
    
    d_dst( mesh_dst%vi1:mesh_dst%vi2) = 0._dp        
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
      DO vli = map%vli1( vi_dst), map%vli2( vi_dst)
        d_dst( vi_dst) = d_dst( vi_dst) + (d_src( map%vi( vli)) * map%w0( vli))
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
    INTEGER                                                :: vi_dst, vli, k
    
    d_dst( mesh_dst%vi1:mesh_dst%vi2,:) = 0._dp        
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
      DO vli = map%vli1( vi_dst), map%vli2( vi_dst)
        DO k = 1, nz
          d_dst( vi_dst,k) = d_dst( vi_dst,k) + (d_src( map%vi( vli),k) * map%w0( vli))
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
    INTEGER                                                :: vi_dst, vli
    REAL(dp), DIMENSION(:), POINTER                        :: ddx_src, ddy_src
    INTEGER                                                :: wddx_src, wddy_src
    
    CALL allocate_shared_dp_1D(  mesh_src%nV, ddx_src, wddx_src)
    CALL allocate_shared_dp_1D(  mesh_src%nV, ddy_src, wddy_src)
    
    CALL ddx_a_to_a_2D( mesh_src, d_src, ddx_src)
    CALL ddy_a_to_a_2D( mesh_src, d_src, ddy_src)
    
    d_dst( mesh_dst%vi1:mesh_dst%vi2) = 0._dp
        
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
      DO vli = map%vli1( vi_dst), map%vli2( vi_dst)
        d_dst( vi_dst) = d_dst( vi_dst) + (d_src(   map%vi( vli)) * map%w0(  vli)) + &
                                          (ddx_src( map%vi( vli)) * map%w1x( vli)) + &
                                          (ddy_src( map%vi( vli)) * map%w1y( vli))
      END DO
    END DO
    
    CALL deallocate_shared( wddx_src)
    CALL deallocate_shared( wddy_src)
    
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
    INTEGER                                                :: vi_dst, vli, m
    REAL(dp), DIMENSION(:,:), POINTER                      :: ddx_src, ddy_src
    INTEGER                                                :: wddx_src, wddy_src
    
    CALL allocate_shared_dp_2D(  mesh_src%nV, 12, ddx_src, wddx_src)
    CALL allocate_shared_dp_2D(  mesh_src%nV, 12, ddy_src, wddy_src)
    
    CALL ddx_a_to_a_3D( mesh_src, d_src, ddx_src)
    CALL ddy_a_to_a_3D( mesh_src, d_src, ddy_src)
    
    d_dst( mesh_dst%vi1:mesh_dst%vi2,:) = 0._dp
        
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
      DO vli = map%vli1( vi_dst), map%vli2( vi_dst)
        DO m = 1, 12
          d_dst( vi_dst,m) = d_dst( vi_dst,m) + (d_src(   map%vi( vli),m) * map%w0(  vli)) + &
                                                (ddx_src( map%vi( vli),m) * map%w1x( vli)) + &
                                                (ddy_src( map%vi( vli),m) * map%w1y( vli))
        END DO
      END DO
    END DO
    
    CALL deallocate_shared( wddx_src)
    CALL deallocate_shared( wddy_src)
    
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
    INTEGER                                                :: vi_dst, vli, k
    REAL(dp), DIMENSION(:,:), POINTER                      :: ddx_src, ddy_src
    INTEGER                                                :: wddx_src, wddy_src
    
    CALL allocate_shared_dp_2D(  mesh_src%nV, C%nz, ddx_src, wddx_src)
    CALL allocate_shared_dp_2D(  mesh_src%nV, C%nz, ddy_src, wddy_src)
    
    CALL ddx_a_to_a_3D( mesh_src, d_src, ddx_src)
    CALL ddy_a_to_a_3D( mesh_src, d_src, ddy_src)
    
    d_dst( mesh_dst%vi1:mesh_dst%vi2,:) = 0._dp
        
    DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
      DO vli = map%vli1( vi_dst), map%vli2( vi_dst)
        DO k = 1, nz
          d_dst( vi_dst,k) = d_dst( vi_dst,k) + (d_src(   map%vi( vli),k) * map%w0(  vli)) + &
                                                (ddx_src( map%vi( vli),k) * map%w1x( vli)) + &
                                                (ddy_src( map%vi( vli),k) * map%w1y( vli))
        END DO
      END DO
    END DO
    
    CALL deallocate_shared( wddx_src)
    CALL deallocate_shared( wddy_src)
    
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
  
! == Smoothing operations on the mesh
  SUBROUTINE smooth_Gaussian_2D( mesh, grid, d_mesh, r)
    ! Use pseudo-conservative remapping to map the data from the mesh
    ! to the square grid. Apply the smoothing on the gridded data, then map
    ! it back to the mesh. The numerical diffusion arising from the two mapping
    ! operations is not a problem since we're smoothing the data anyway.
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),              INTENT(INOUT) :: d_mesh
    REAL(dp),                                INTENT(IN)    :: r
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER                    :: d_grid
    INTEGER                                                :: wd_grid
    
    ! Allocate shared memory for the gridded data
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid, wd_grid)
    
    ! Map data from the mesh to the grid
    CALL map_mesh2grid_2D( mesh, grid, d_mesh, d_grid)
    
    ! Smooth data on the grid
    CALL smooth_Gaussian_2D_grid( grid, d_grid, r)
    
    ! Map smoothed data back to the mesh
    CALL map_grid2mesh_2D( mesh, grid, d_grid, d_mesh)
    
    ! Deallocate gridded data
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE smooth_Gaussian_2D
  SUBROUTINE smooth_Gaussian_3D( mesh, grid, d_mesh, r, nz)
    ! Use pseudo-conservative remapping to map the data from the mesh
    ! to the square grid. Apply the smoothing on the gridded data, then map
    ! it back to the mesh. The numerical diffusion arising from the two mapping
    ! operations is not a problem since we're smoothing the data anyway.
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: d_mesh
    REAL(dp),                                INTENT(IN)    :: r
    INTEGER,                                 INTENT(IN)    :: nz
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:,:), POINTER                    :: d_grid
    INTEGER                                                :: wd_grid
    
    ! Allocate shared memory for the gridded data
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, nz, d_grid, wd_grid)
    
    ! Map data from the mesh to the grid
    CALL map_mesh2grid_3D( mesh, grid, d_mesh, d_grid)
    
    ! Smooth data on the grid
    CALL smooth_Gaussian_3D_grid( grid, d_grid, r, nz)
    
    ! Map smoothed data back to the mesh
    CALL map_grid2mesh_3D( mesh, grid, d_grid, d_mesh)
    
    ! Deallocate gridded data
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE smooth_Gaussian_3D
  
END MODULE mesh_mapping_module
