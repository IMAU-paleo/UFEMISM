MODULE mesh_rectangular_module

  USE mpi
  USE configuration_module,          ONLY: dp, C
  USE parallel_module,               ONLY: par, sync, allocate_shared_int_0D, allocate_shared_dp_0D, &
                                                      allocate_shared_int_1D, allocate_shared_dp_1D, &
                                                      allocate_shared_int_2D, allocate_shared_dp_2D, &
                                                      allocate_shared_int_3D, allocate_shared_dp_3D, &
                                                      allocate_shared_bool_1D, deallocate_shared
  USE data_types_netcdf_module,      ONLY: type_netcdf_restart, type_netcdf_help_fields
  USE data_types_module,             ONLY: type_mesh, type_model_region, type_grid
  USE mesh_help_functions_module,    ONLY: get_lat_lon_coordinates, check_mesh, find_Voronoi_cell_areas, inverse_oblique_sg_projection, &
                                           switch_vertices, find_Voronoi_cell_geometric_centres, partition_list, update_triangle_circumcenter
  USE mesh_memory_module,            ONLY: allocate_mesh_primary, allocate_mesh_secondary
  USE mesh_Delaunay_module,          ONLY: ensure_Delaunay
  USE mesh_ArakawaC_module,          ONLY: make_Ac_mesh
  USE mesh_mapping_module,           ONLY: create_remapping_arrays, remap_cons_1st_order_2D, remap_cons_1st_order_3D, remap_cons_2nd_order_2D, remap_cons_2nd_order_3D

  IMPLICIT NONE
  
  CONTAINS
  
  ! Create a regular square grid and a "semi-rectangular mesh"
  SUBROUTINE create_semi_rectangular_mesh( region, grid, dx)
    ! Create a structured "semi-rectangular" mesh. These are used for writing output on a square
    ! grid, running a GIA model (either ELRA or SELEN, both of which need a square grid), and
    ! smoothing data in the climate matrix (I've yet to figure out a conservative smoothing algorithm
    ! on an unstructured mesh). However, the size of the domain covered by the model mesh might
    ! not be an integer multiple of the desired square grid resolution. Since the conservative remapping
    ! code requires two meshes to square the same domain boundary, a little "nudging" is required.
    ! We allow the square grid to extend outside of the model domain. The square mesh corresponding
    ! to the grid is slightly deformed, so that the boundary vertices coincide with the actual model
    ! domain boundary. This means that remapping is not strictly conservative for the boundary pixels,
    ! but since most data fields we want to remap are zero on the boundary anyway, this is okay.
    
    USE parameters_module, ONLY: pi
    
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    TYPE(type_grid),            INTENT(INOUT)     :: grid
    REAL(dp),                   INTENT(IN)        :: dx
        
    ! Local variables:
    REAL(dp)                                      :: xmid, ymid
    INTEGER                                       :: nsx, nsy, nx, ny
    INTEGER                                       :: nV, nTri, i, j, nTri_per_col, til, tir, vi, v1, v2, v3, v4
    INTEGER                                       :: ia,ja,via,ib,jb,vib,ti
    INTEGER                                       :: vi_SE, vi_NE, vi_NW
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp
    
    ! Allocate shared memory
    CALL allocate_shared_int_0D( grid%nx,           grid%wnx          )
    CALL allocate_shared_int_0D( grid%ny,           grid%wny          )
    CALL allocate_shared_dp_0D(  grid%dx,           grid%wdx          )
    CALL allocate_shared_dp_0D(  grid%xmin,         grid%wxmin        )
    CALL allocate_shared_dp_0D(  grid%xmax,         grid%wxmax        )
    CALL allocate_shared_dp_0D(  grid%ymin,         grid%wymin        )
    CALL allocate_shared_dp_0D(  grid%ymax,         grid%wymax        )
    CALL allocate_shared_dp_0D(  grid%lambda_M,     grid%wlambda_M    )
    CALL allocate_shared_dp_0D(  grid%phi_M,        grid%wphi_M       )
    CALL allocate_shared_dp_0D(  grid%alpha_stereo, grid%walpha_stereo)
  
    ! Projection parameters are copied from the mesh
    IF (par%master) THEN
      grid%lambda_M     = region%mesh%lambda_M
      grid%phi_M        = region%mesh%phi_M
      grid%alpha_stereo = region%mesh%alpha_stereo
    END IF
    CALL sync
    
    ! Determine the number of grid cells
    ! (it might well be that the model domain size is not an integer multiple of the square grid resolution)
    
    ! Domain centre
    xmid = (region%mesh%xmax + region%mesh%xmin) / 2._dp
    ymid = (region%mesh%ymax + region%mesh%ymin) / 2._dp
    
    ! Half-maximum number of square grid cells that can fit inside the model domain
    ! (since we want there to be an uneven number, with the centre one lying exactly in the middle)
    nsx = FLOOR( 0.5_dp * (region%mesh%xmax - region%mesh%xmin) / dx)
    nsy = FLOOR( 0.5_dp * (region%mesh%ymax - region%mesh%ymin) / dx)
    
    ! Check if we need to extend this by one cell to fully cover the model domain.
    IF (xmid + REAL(nsx,dp) * dx < region%mesh%xmax) nsx = nsx + 1
    IF (ymid + REAL(nsy,dp) * dx < region%mesh%ymax) nsy = nsy + 1

    ! Total number of grid cells
    IF (par%master) THEN
      grid%dx = dx
      grid%nx = 1 + 2*nsx
      grid%ny = 1 + 2*nsy
    END IF
    CALL sync
    
    !IF (par%master) WRITE (0,'(A,F5.2,A,I4,A,I4,A)') '   Initialising square grid at ', dx/1000._dp, ' km resolution: [', grid%nx, ' x ', grid%ny, '] pixels'
    
    ! Assign range to each processor
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)
    
    ! Allocate shared memory for x and y
    CALL allocate_shared_dp_1D( grid%nx, grid%x, grid%wx)
    CALL allocate_shared_dp_1D( grid%ny, grid%y, grid%wy)
    
    ! Fill in x and y
    IF (par%master) THEN
      DO i = 1, grid%nx
        grid%x( i) = -nsx*grid%dx + (i-1)*grid%dx
      END DO
      DO j = 1, grid%ny
        grid%y( j) = -nsy*grid%dx + (j-1)*grid%dx
      END DO
      
      grid%xmin = MINVAL(grid%x)
      grid%xmax = MAXVAL(grid%x)
      grid%ymin = MINVAL(grid%y)
      grid%ymax = MAXVAL(grid%y)
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Lat,lon coordinates
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, grid%lat, grid%wlat)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, grid%lon, grid%wlon)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL inverse_oblique_sg_projection( grid%x( i), grid%y( j), grid%lambda_M, grid%phi_M, grid%alpha_stereo, grid%lon( j,i), grid%lat( j,i))
    END DO
    END DO
    CALL sync
    
    ! Determine number of vertices and triangles, allocate mesh memory
    nV   = grid%nx * grid%ny
    nTri = 2 * (grid%nx-1) * (grid%ny-1)
    CALL allocate_mesh_primary( grid%mesh, region%name, nV, nTri, C%nconmax)
    
    ! Maps between the square mesh and the actual 2D grid
    CALL allocate_shared_int_2D( nV,      2,       grid%map_square2grid, grid%wmap_square2grid)
    CALL allocate_shared_int_2D( grid%nx, grid%ny, grid%map_grid2square, grid%wmap_grid2square)
    
    IF (par%master) THEN
        
      ! Metadata
      ! ========
      
      grid%mesh%xmin             = region%mesh%xmin
      grid%mesh%xmax             = region%mesh%xmax
      grid%mesh%ymin             = region%mesh%ymin
      grid%mesh%ymax             = region%mesh%ymax
      grid%mesh%tol_dist         = region%mesh%tol_dist
      grid%mesh%nV               = nV
      grid%mesh%nTri             = nTri
      grid%mesh%perturb_dir      = 0
      grid%mesh%alpha_min        = 45._dp * (pi / 180._dp)
      grid%mesh%dz_max_ice       = C%dz_max_ice
      grid%mesh%res_max          = grid%dx
      grid%mesh%res_max_margin   = grid%dx
      grid%mesh%res_max_gl       = grid%dx
      grid%mesh%res_max_cf       = grid%dx
      grid%mesh%res_max_mountain = grid%dx
      grid%mesh%res_max_coast    = grid%dx
      grid%mesh%res_min          = grid%dx
      grid%mesh%resolution_min   = grid%dx
      grid%mesh%resolution_max   = grid%dx
      
      ! ===== Primary data
      ! ===================
          
      nTri_per_col = 2 * (grid%ny-1)
      
      nx = grid%nx
      ny = grid%ny
    
      DO i = 1, nx
        DO j = 1, ny
  
          ! Vertex index
          vi = (i-1)*ny + j
          
          ! Maps
          grid%map_square2grid( vi,:) = [i,j]
          grid%map_grid2square( i,j ) = vi
  
          ! Vertex coordinates
          grid%mesh%V(vi,:) = [grid%mesh%xmin + REAL(i-1)*dx, grid%mesh%ymin + REAL(j-1)*dx]
  
          ! Triangle indices
          til = (i-1)*nTri_per_col + (j-1)*2 + 1
          tir = (i-1)*nTri_per_col + (j-1)*2 + 2
  
          v1 = vi
          v2 = vi + ny
          v3 = v2 + 1
          v4 = v1 + 1
  
          IF (i < nx .AND. j < ny) THEN
            grid%mesh%Tri(til,:) = [v1, v2, v4]
            grid%mesh%Tri(tir,:) = [v2, v3, v4]
          END IF
  
          ! Determine other properties - vertex connectivity, inverse triangles,
          ! edge index, etc.
          IF (i==1) THEN
            IF (j==1) THEN
              ! SW corner
              grid%mesh%nC(              vi    ) = 2
              grid%mesh%C(               vi,1:2) = [ny+1, 2]
              grid%mesh%niTri(           vi    ) = 1
              grid%mesh%iTri(            vi,1  ) = 1
              grid%mesh%edge_index(      vi    ) = 6
  
              grid%mesh%Tricc(           til,:)  = grid%mesh%V(vi,:) + [dx/2, dx/2]
              grid%mesh%TriC(            til,:)  = [2, 0, 0]
              grid%mesh%Tri_edge_index(  til  )  = 5
  
              grid%mesh%Tricc(           tir,:)  = grid%mesh%V(vi,:) + [dx/2, dx/2]
              grid%mesh%TriC(            tir,:)  = [3, 1, nTri_per_col+1]
              grid%mesh%Tri_edge_index(  tir  )  = 0
            ELSEIF (j==ny) THEN
              ! NW corner
              vi_NW = vi
              grid%mesh%nC(              vi    ) = 3
              grid%mesh%C(               vi,1:3) = [ny-1, 2*ny-1, 2*ny]
              grid%mesh%niTri(           vi    ) = 2
              grid%mesh%iTri(            vi,1:2) = [nTri_per_col-1, nTri_per_col]
              grid%mesh%edge_index(      vi    ) = 8
            ELSE
              ! W edge
              grid%mesh%nC(              vi    ) = 4
              grid%mesh%C(               vi,1:4) = [vi-1, vi+ny-1, vi+ny, vi+1]
              grid%mesh%niTri(           vi    ) = 3
              grid%mesh%iTri(            vi,1:3) = [vi*2-3, vi*2-2, vi*2-1]
              grid%mesh%edge_index(      vi    ) = 7
  
              grid%mesh%Tricc(           til,:)  = grid%mesh%V(vi,:) + [dx/2, dx/2]
              grid%mesh%TriC(            til,:)  = [til+1, 0, til-1]
              grid%mesh%Tri_edge_index(  til  )  = 7
  
              grid%mesh%Tricc(           tir,:)  = grid%mesh%V(vi,:) + [dx/2, dx/2]
              grid%mesh%TriC(            tir,:)  = [tir+1, tir-1, tir+nTri_per_col-1]
              grid%mesh%Tri_edge_index(  tir  )  = 0
  
              IF (j==ny-1) THEN
                grid%mesh%TriC(          tir,:)  = [0, tir-1, tir+nTri_per_col-1]
                grid%mesh%Tri_edge_index(tir  )  = 1
              END IF
            END IF
          ELSEIF (i==nx) THEN
            IF (j==1) THEN
              ! SE corner
              vi_SE = vi
              grid%mesh%nC(              vi    ) = 3
              grid%mesh%C(               vi,1:3) = [vi+1, vi-ny+1, vi-ny]
              grid%mesh%niTri(           vi    ) = 2
              grid%mesh%iTri(            vi,1:2) = [nTri_per_col*(nx-2)+2, nTri_per_col*(nx-2)+1]
              grid%mesh%edge_index(      vi    ) = 4
            ELSEIF (j==ny) THEN
              ! NE corner
              vi_NE = vi
              grid%mesh%nC(              vi    ) = 2
              grid%mesh%C(               vi,1:2) = [vi-ny, vi-1]
              grid%mesh%niTri(           vi    ) = 1
              grid%mesh%iTri(            vi,1  ) = grid%mesh%nTri
              grid%mesh%edge_index(      vi    ) = 2
            ELSE
              ! E edge
              grid%mesh%nC(              vi    ) = 4
              grid%mesh%C(               vi,1:4) = [vi+1, vi-ny+1, vi-ny, vi-1]
              grid%mesh%niTri(           vi    ) = 3
              grid%mesh%iTri(            vi,1:3) = [nTri_per_col*(nx-2)+(j*2), nTri_per_col*(nx-2)+(j*2)-1, nTri_per_col*(nx-2)+(j*2)-2]
              grid%mesh%edge_index(      vi    ) = 3
            END IF
          ELSE
            if (j==1) THEN
              ! S edge
              grid%mesh%nC(              vi    ) = 4
              grid%mesh%C(               vi,1:4) = [vi+ny, vi+1, vi-ny+1, vi-ny]
              grid%mesh%niTri(           vi    ) = 3
              grid%mesh%iTri(            vi,1:3) = [nTri_per_col*(i-1)+1, nTri_per_col*(i-2)+2, nTri_per_col*(i-2)+1]
              grid%mesh%edge_index(      vi    ) = 5
  
              grid%mesh%Tricc(           til,:)  = grid%mesh%V(vi,:) + [dx/2, dx/2]
              grid%mesh%TriC(            til,:)  = [til+1, til-nTri_per_col+1, 0]
              grid%mesh%Tri_edge_index(  til  )  = 5
  
              grid%mesh%Tricc(           tir,:)  = grid%mesh%V(vi,:) + [dx/2, dx/2]
              grid%mesh%TriC(            tir,:)  = [tir+1, tir-1, tir+nTri_per_col-1]
              grid%mesh%Tri_edge_index(  tir  )  = 0
  
              IF (i==nx-1) THEN
                grid%mesh%TriC(          tir,:)  = [tir+1, tir-1, 0]
                grid%mesh%Tri_edge_index(tir  )  = 3
              END IF
            ELSEIF (j==ny) THEN
              ! N edge
              grid%mesh%nC(              vi    ) = 4
              grid%mesh%C(               vi,1:4) = [vi-ny, vi-1, vi+ny-1, vi+ny]
              grid%mesh%niTri(           vi    ) = 3
              grid%mesh%iTri(            vi,1:3) = [nTri_per_col*(i-1), nTri_per_col*i-1, nTri_per_col*i]
              grid%mesh%edge_index(      vi    ) = 1
            ELSE
              ! non-edge vertex
              grid%mesh%nC(              vi    ) = 6
              grid%mesh%C(               vi,1:6) = [vi+ny, vi+1, vi-ny+1, vi-ny, vi-1, vi+ny-1]
              grid%mesh%niTri(           vi    ) = 6
              grid%mesh%iTri(            vi,1:6) = [ nTri_per_col*(i-1)+j*2-1, &
                                                nTri_per_col*(i-2)+j*2, &
                                                nTri_per_col*(i-2)+j*2-1, &
                                                nTri_per_col*(i-2)+j*2-2, &
                                                nTri_per_col*(i-1)+j*2-3, &
                                                nTri_per_col*(i-1)+j*2-2]
              grid%mesh%edge_index(      vi    ) = 0
  
              grid%mesh%Tricc(           til,:)  = grid%mesh%V(vi,:) + [dx/2, dx/2]
              grid%mesh%TriC(            til,:)  = [til+1, til-nTri_per_col+1, til-1]
              grid%mesh%Tri_edge_index(  til  )  = 0
  
              grid%mesh%Tricc(           tir,:)  = grid%mesh%V(vi,:) + [dx/2, dx/2]
              grid%mesh%TriC(            tir,:)  = [tir+1, tir-1, tir+nTri_per_col-1]
              grid%mesh%Tri_edge_index(  tir  )  = 0
  
              IF (i<nx-1 .AND. j==ny-1) THEN
                grid%mesh%TriC(          tir,:)  = [0, tir-1, tir+nTri_per_col-1]
                grid%mesh%Tri_edge_index(tir  )  = 1
              ELSEIF (i==nx-1 .AND. j<ny-1) THEN
                grid%mesh%TriC(          tir,:)  = [tir+1, tir-1, 0]
                grid%mesh%Tri_edge_index(tir  )  = 3
              ELSEIF (i==nx-1 .AND. j==ny-1) THEN
                grid%mesh%TriC(          tir,:)  = [0, tir-1, 0]
                grid%mesh%Tri_edge_index(tir  )  = 1
              END IF
            END IF
          END IF 
  
        END DO
      END DO
      
      ! Make sure vertices 2, 3 and 4 are at the SE, NE and NW corners (some routines expect this)
      ia  = grid%map_square2grid( 2,    1)
      ja  = grid%map_square2grid( 2,    2)
      via = 2
      ib  = grid%map_square2grid( vi_SE,1)
      jb  = grid%map_square2grid( vi_SE,2)
      vib = vi_SE
      CALL switch_vertices( grid%mesh, 2, vi_SE)
      grid%map_square2grid( 2    ,:) = [ib,jb]
      grid%map_square2grid( vi_SE,:) = [ia,ja]
      grid%map_grid2square( ia,ja) = vi_SE
      grid%map_grid2square( ib,jb) = 2
      
      ia  = grid%map_square2grid( 3,    1)
      ja  = grid%map_square2grid( 3,    2)
      via = 3
      ib  = grid%map_square2grid( vi_NE,1)
      jb  = grid%map_square2grid( vi_NE,2)
      vib = vi_NE
      CALL switch_vertices( grid%mesh, 3, vi_NE)
      grid%map_square2grid( 3    ,:) = [ib,jb]
      grid%map_square2grid( vi_NE,:) = [ia,ja]
      grid%map_grid2square( ia,ja) = vi_NE
      grid%map_grid2square( ib,jb) = 3
      
      ia  = grid%map_square2grid( 4,    1)
      ja  = grid%map_square2grid( 4,    2)
      via = 4
      ib  = grid%map_square2grid( vi_NW,1)
      jb  = grid%map_square2grid( vi_NW,2)
      vib = vi_NW
      CALL switch_vertices( grid%mesh, 4, vi_NW)
      grid%map_square2grid( 4    ,:) = [ib,jb]
      grid%map_square2grid( vi_NW,:) = [ia,ja]
      grid%map_grid2square( ia,ja) = vi_NW
      grid%map_grid2square( ib,jb) = 4
      
      ! Move boundary vertices to the domain boundary
      DO i = 1, grid%nx
        ! South
        j = 1
        vi = grid%map_grid2square( i,j)
        grid%mesh%V( vi,2) = grid%mesh%ymin
        ! North
        j = grid%ny
        vi = grid%map_grid2square( i,j)
        grid%mesh%V( vi,2) = grid%mesh%ymax
      END DO
      DO j = 1, grid%ny
        ! West
        i = 1
        vi = grid%map_grid2square( i,j)
        grid%mesh%V( vi,1) = grid%mesh%xmin
        ! East
        i = grid%nx
        vi = grid%map_grid2square( i,j)
        grid%mesh%V( vi,1) = grid%mesh%xmax
      END DO
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Determine vertex and triangle domains
    CALL partition_list( grid%mesh%nV,   par%i, par%n, grid%mesh%v1, grid%mesh%v2)
    CALL partition_list( grid%mesh%nTri, par%i, par%n, grid%mesh%t1, grid%mesh%t2)
    
    ! The squeezing of the outermost rows of triangles can sometimes violate the Delaunay criterion
    CALL ensure_Delaunay( grid%mesh)
    
    ! Update triangle circumcenters
    DO ti = grid%mesh%t1, grid%mesh%t1
      CALL update_triangle_circumcenter( grid%mesh, ti)
    END DO
    CALL sync
        
    ! Check if we did everything correctly
    CALL check_mesh( grid%mesh)
    
    ! Calculate extra mesh data (lat,lon for SELEN, VorGC for conservative remapping)
    CALL allocate_mesh_secondary(             grid%mesh)
    CALL make_Ac_mesh(                        grid%mesh)
    CALL find_Voronoi_cell_areas(             grid%mesh)
    CALL get_lat_lon_coordinates(             grid%mesh)
    CALL find_Voronoi_cell_geometric_centres( grid%mesh)
    
    ! Create remapping arrays between the square mesh and the model mesh
    CALL create_remapping_arrays( region%mesh, grid%mesh,   grid%map_mesh2square)
    CALL create_remapping_arrays( grid%mesh,   region%mesh, grid%map_square2mesh)
    
  END SUBROUTINE create_semi_rectangular_mesh
  
  ! Convert data between a grid (2D array) and a square mesh (1D array)
  SUBROUTINE map_square_to_grid_2D( mesh, map_square2grid, d_square, d_grid)
    ! Map data from the square mesh to a real 2D grid
    
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,  DIMENSION(:,:  ), INTENT(IN)        :: map_square2grid
    REAL(dp), DIMENSION(:    ), INTENT(IN)        :: d_square
    REAL(dp), DIMENSION(:,:  ), INTENT(OUT)       :: d_grid
    
    ! Local variables:
    INTEGER                                       :: vi, i, j
    
    DO vi = mesh%v1, mesh%v2
      i = map_square2grid( vi,1)
      j = map_square2grid( vi,2)
      d_grid( i,j) = d_square( vi)
    END DO
    CALL sync
  
  END SUBROUTINE map_square_to_grid_2D
  SUBROUTINE map_square_to_grid_3D( mesh, map_square2grid, d_square, d_grid)
    ! Map data from the square mesh to a real 2D grid
    
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,  DIMENSION(:,:  ), INTENT(IN)        :: map_square2grid
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: d_square
    REAL(dp), DIMENSION(:,:,:), INTENT(OUT)       :: d_grid
    
    ! Local variables:
    INTEGER                                       :: vi, i, j
    
    DO vi = mesh%v1, mesh%v2
      i = map_square2grid( vi,1)
      j = map_square2grid( vi,2)
      d_grid( i,j,:) = d_square( vi,:)
    END DO
    CALL sync
  
  END SUBROUTINE map_square_to_grid_3D
  SUBROUTINE map_grid_to_square_2D( mesh, map_square2grid, d_grid, d_square)
    ! Map data from the square mesh to a real 2D grid
    
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,  DIMENSION(:,:  ), INTENT(IN)        :: map_square2grid
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: d_grid
    REAL(dp), DIMENSION(:    ), INTENT(OUT)       :: d_square
    
    ! Local variables:
    INTEGER                                       :: vi, i, j
    
    DO vi = mesh%v1, mesh%v2
      i = map_square2grid( vi,1)
      j = map_square2grid( vi,2)
      d_square( vi) = d_grid( i,j)
    END DO
    CALL sync
  
  END SUBROUTINE map_grid_to_square_2D
  SUBROUTINE map_grid_to_square_3D( mesh, map_square2grid, d_grid, d_square)
    ! Map data from the square mesh to a real 2D grid
    
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,  DIMENSION(:,:  ), INTENT(IN)        :: map_square2grid
    REAL(dp), DIMENSION(:,:,:), INTENT(IN)        :: d_grid
    REAL(dp), DIMENSION(:,:  ), INTENT(OUT)       :: d_square
    
    ! Local variables:
    INTEGER                                       :: vi, i, j
    
    DO vi = mesh%v1, mesh%v2
      i = map_square2grid( vi,1)
      j = map_square2grid( vi,2)
      d_square( vi,:) = d_grid( i,j,:)
    END DO
    CALL sync
  
  END SUBROUTINE map_grid_to_square_3D
  
  ! Write to the "grid" restart and help_fields files (routines are here instead of
  ! the netcdf_module because they require access to the remapping routines)
  SUBROUTINE write_to_restart_file_grid( region, netcdf)
  
    USE netcdf_module, ONLY: open_netcdf_file, handle_error, close_netcdf_file
    USE netcdf,        ONLY: nf90_put_var
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),   INTENT(INOUT) :: region
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf
    
    ! Open the file for writing
    IF (par%master) CALL open_netcdf_file( netcdf%filename, netcdf%ncid)
        
    ! Time
    IF (par%master) CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time,             region%time,      start=(/         netcdf%ti/)))
    
    ! Map and write data
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hi,               netcdf%id_var_Hi,               netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hb,               netcdf%id_var_Hb,               netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hs,               netcdf%id_var_Hs,               netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%U_SIA,            netcdf%id_var_U_SIA,            netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%V_SIA,            netcdf%id_var_V_SIA,            netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%U_SSA,            netcdf%id_var_U_SSA,            netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%V_SSA,            netcdf%id_var_V_SSA,            netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti,               netcdf%id_var_Ti,               netcdf%ti, C%nZ)
    CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%FirnDepth,        netcdf%id_var_FirnDepth,        netcdf%ti, 12  )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%MeltPreviousYear, netcdf%id_var_MeltPreviousYear, netcdf%ti      )
    
    ! Close the file
    IF (par%master) CALL close_netcdf_file(netcdf%ncid)
    
    ! Increase time frame counter
    IF (par%master) netcdf%ti = netcdf%ti + 1
    CALL sync
        
  END SUBROUTINE write_to_restart_file_grid
  SUBROUTINE write_to_help_fields_file_grid( region, netcdf)
  
    USE netcdf_module, ONLY: open_netcdf_file, handle_error, close_netcdf_file
    USE netcdf,        ONLY: nf90_put_var
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),       INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields), INTENT(INOUT) :: netcdf
    
    ! Open the file for writing
    IF (par%master) CALL open_netcdf_file( netcdf%filename, netcdf%ncid)
        
    ! Time
    IF (par%master) CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time, region%time, start=(/ netcdf%ti/)))
    
    ! Write data
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_01, C%help_field_01)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_02, C%help_field_02)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_03, C%help_field_03)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_04, C%help_field_04)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_05, C%help_field_05)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_06, C%help_field_06)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_07, C%help_field_07)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_08, C%help_field_08)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_09, C%help_field_09)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_10, C%help_field_10)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_11, C%help_field_11)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_12, C%help_field_12)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_13, C%help_field_13)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_14, C%help_field_14)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_15, C%help_field_15)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_16, C%help_field_16)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_17, C%help_field_17)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_18, C%help_field_18)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_19, C%help_field_19)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_20, C%help_field_20)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_21, C%help_field_21)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_22, C%help_field_22)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_23, C%help_field_23)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_24, C%help_field_24)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_25, C%help_field_25)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_26, C%help_field_26)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_27, C%help_field_27)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_28, C%help_field_28)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_29, C%help_field_29)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_30, C%help_field_30)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_31, C%help_field_31)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_32, C%help_field_32)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_33, C%help_field_33)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_34, C%help_field_34)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_35, C%help_field_35)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_36, C%help_field_36)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_37, C%help_field_37)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_38, C%help_field_38)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_39, C%help_field_39)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_40, C%help_field_40)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_41, C%help_field_41)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_42, C%help_field_42)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_43, C%help_field_43)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_44, C%help_field_44)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_45, C%help_field_45)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_46, C%help_field_46)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_47, C%help_field_47)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_48, C%help_field_48)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_49, C%help_field_49)
    CALL write_help_field_grid( region, netcdf, netcdf%id_help_field_50, C%help_field_50)
    
    ! Close the file
    IF (par%master) CALL close_netcdf_file(netcdf%ncid)
    
    ! Increase time frame counter
    IF (par%master) netcdf%ti = netcdf%ti + 1
    CALL sync
        
  END SUBROUTINE write_to_help_fields_file_grid
  SUBROUTINE write_help_field_grid( region, netcdf, id_var, field_name)
    ! Write the current model state to the existing output file
  
    USE netcdf_module, ONLY: open_netcdf_file, handle_error, close_netcdf_file
    USE netcdf,        ONLY: nf90_put_var
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),       INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields), INTENT(INOUT) :: netcdf
    INTEGER,                       INTENT(IN)    :: id_var
    CHARACTER(LEN=*),              INTENT(IN)    :: field_name
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: vi
    REAL(dp), DIMENSION(:    ), POINTER           :: d_year
    INTEGER                                       :: wd_year
    
    CALL allocate_shared_dp_1D( region%mesh%nV, d_year, wd_year)
    
    IF (field_name == 'none') THEN
      RETURN
      
    ! Fields with no time dimension
    ! =============================
      
    ! Lat/lon
    ELSEIF (field_name == 'lat') THEN
      IF (par%master) CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%grid_output%lat ))
    ELSEIF (field_name == 'lon') THEN
      IF (par%master) CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%grid_output%lon ))
    
!    ! Geothermal heat flux
!    ELSEIF (field_name == 'GHF') THEN
!      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%GHF_Aa,          (/1, 1        /))
      
    ! Fields with a time dimension
    ! ============================
      
    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hi, id_var, netcdf%ti)
    ELSEIF (field_name == 'Hb') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hb, id_var, netcdf%ti)
    ELSEIF (field_name == 'Hs') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hs, id_var, netcdf%ti)
    ELSEIF (field_name == 'SL') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%SL, id_var, netcdf%ti)
    ELSEIF (field_name == 'dHs_dx') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dHs_dx, id_var, netcdf%ti)
    ELSEIF (field_name == 'dHs_dy') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dHs_dy, id_var, netcdf%ti)
      
    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'Cpi') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Cpi, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'Ki') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ki, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti(:,C%nZ), id_var, netcdf%ti)
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti_pmp, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'A_flow') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%A_flow, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'A_flow_mean') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%A_flow_mean, id_var, netcdf%ti)
      
    ! Velocity fields
    ELSEIF (field_name == 'U_SIA') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%U_SIA, id_var, netcdf%ti)
    ELSEIF (field_name == 'V_SIA') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%V_SIA, id_var, netcdf%ti)
    ELSEIF (field_name == 'U_SSA') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%U_SSA, id_var, netcdf%ti)
    ELSEIF (field_name == 'V_SSA') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%V_SSA, id_var, netcdf%ti)
    ELSEIF (field_name == 'U_vav') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%U_vav, id_var, netcdf%ti)
    ELSEIF (field_name == 'V_vav') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%V_vav, id_var, netcdf%ti)
    ELSEIF (field_name == 'U_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%U_surf, id_var, netcdf%ti)
    ELSEIF (field_name == 'V_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%V_surf, id_var, netcdf%ti)
    ELSEIF (field_name == 'U_base') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%U_base, id_var, netcdf%ti)
    ELSEIF (field_name == 'V_base') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%V_base, id_var, netcdf%ti)
    ELSEIF (field_name == 'U_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%U_3D, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'V_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%V_3D, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'W_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%W_3D, id_var, netcdf%ti, C%nZ)
      
    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate%applied%T2m, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'T2m_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%climate%applied%T2m( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'Precip') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate%applied%Precip, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Precip_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%climate%applied%Precip( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate%applied%Wind_WE, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Wind_WE_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%climate%applied%Wind_WE( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate%applied%Wind_SN, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Wind_SN_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%climate%applied%Wind_SN( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
      
    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%SMB, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'SMB_year') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%SMB_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%BMB%BMB_sheet, id_var, netcdf%ti)
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%BMB%BMB_shelf, id_var, netcdf%ti)
    ELSEIF (field_name == 'BMB') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%BMB%BMB, id_var, netcdf%ti)
    ELSEIF (field_name == 'Snowfall') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Snowfall, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Snowfall_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%SMB%Snowfall( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'Rainfall') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Rainfall, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Rainfall_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%SMB%Rainfall( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%AddedFirn, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'AddedFirn_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%SMB%AddedFirn( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'Refreezing') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Refreezing, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Refreezing_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%SMB%Refreezing( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'Runoff') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Runoff, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Runoff_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%SMB%Runoff( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'Albedo') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Albedo, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Albedo_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%SMB%Albedo( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%FirnDepth, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'FirnDepth_year') THEN
      DO vi = region%mesh%v1, region%mesh%v2
        d_year( vi) = SUM( region%SMB%FirnDepth( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
      
      ! NOTE: masks commented out; mapping masks between grids is meaningless
      
!    ! Masks
!    ELSEIF (field_name == 'mask') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask, id_var, netcdf%ti)
!    ELSEIF (field_name == 'mask_land') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_land, id_var, netcdf%ti)
!    ELSEIF (field_name == 'mask_ocean') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_ocean, id_var, netcdf%ti)
!    ELSEIF (field_name == 'mask_lake') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_lake, id_var, netcdf%ti)
!    ELSEIF (field_name == 'mask_ice') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_ice, id_var, netcdf%ti)
!    ELSEIF (field_name == 'mask_sheet') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_sheet, id_var, netcdf%ti)
!    ELSEIF (field_name == 'mask_shelf') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_shelf, id_var, netcdf%ti)
!    ELSEIF (field_name == 'mask_coast') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_coast, id_var, netcdf%ti)
!    ELSEIF (field_name == 'mask_margin') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_margin, id_var, netcdf%ti)
!    ELSEIF (field_name == 'mask_gl') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_gl, id_var, netcdf%ti)
!    ELSEIF (field_name == 'mask_cf') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_cf, id_var, netcdf%ti)
      
    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      d_year( region%mesh%v1:region%mesh%v2) = region%ice%phi_fric_AaAc( region%mesh%v1:region%mesh%v2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    ELSEIF (field_name == 'tau_yield') THEN
      d_year( region%mesh%v1:region%mesh%v2) = region%ice%tau_c_AaAc( region%mesh%v1:region%mesh%v2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
      
!    ! Isotopes
!    ELSEIF (field_name == 'iso_ice') THEN
!      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%IsoIce,            (/1, ti /))
!    ELSEIF (field_name == 'iso_surf') THEN
!      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%IsoSurf,           (/1, ti /))
    
    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      d_year( region%mesh%v1:region%mesh%v2) = region%ice%Hb( region%mesh%v1:region%mesh%v2) - region%PD%Hb( region%mesh%v1:region%mesh%v2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, d_year, id_var, netcdf%ti)
    
    ELSE
      WRITE(0,*) ' ERROR: help field "', TRIM(field_name), '" not implemented in write_help_field!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    CALL deallocate_shared( wd_year)
    
  END SUBROUTINE write_help_field_grid
  
  ! Use the existing remapping routines to map a model data field from the model
  ! mesh to the square mesh, then to the grid, then finally write it to a NetCDF file.
  SUBROUTINE map_and_write_to_grid_netcdf_int_2D( ncid, mesh, grid, d, id_var, ti)
  
    USE netcdf_module, ONLY: open_netcdf_file, handle_error, close_netcdf_file
    USE netcdf,        ONLY: nf90_put_var
   
    IMPLICIT NONE
    
    ! Input variables:
    INTEGER,                    INTENT(IN)        :: ncid
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_grid),            INTENT(INOUT)     :: grid
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: d
    INTEGER,                    INTENT(IN)        :: id_var
    INTEGER,                    INTENT(IN)        :: ti
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER           :: d_dp
    REAL(dp), DIMENSION(:    ), POINTER           :: d_square
    REAL(dp), DIMENSION(:,:  ), POINTER           :: d_grid
    INTEGER                                       :: wd_dp, wd_square, wd_grid
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( grid%mesh%nV,     d_dp,     wd_dp    )
    CALL allocate_shared_dp_1D( grid%mesh%nV,     d_square, wd_square)
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid,   wd_grid  )
    
    ! Convert from integer to real
    d_dp( mesh%v1:mesh%v2) = REAL( d( mesh%v1:mesh%v2), dp)
    
    ! Map data from the model mesh to the square mesh
    CALL remap_cons_1st_order_2D(       grid%mesh, grid%map_mesh2square%conservative, d_dp, d_square)
   ! CALL remap_cons_2nd_order_2D( mesh, grid%mesh, grid%map_mesh2square%conservative, d_dp, d_square)
    
    ! Map data from the square mesh to the grid
    CALL map_square_to_grid_2D( grid%mesh, grid%map_square2grid, d_square, d_grid)
    
    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, ti/) ))
    
    ! Deallocate shared memory
    CALL deallocate_shared( wd_dp)
    CALL deallocate_shared( wd_square)
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE map_and_write_to_grid_netcdf_int_2D
  SUBROUTINE map_and_write_to_grid_netcdf_dp_2D( ncid, mesh, grid, d, id_var, ti)
  
    USE netcdf_module, ONLY: open_netcdf_file, handle_error, close_netcdf_file
    USE netcdf,        ONLY: nf90_put_var
   
    IMPLICIT NONE
    
    ! Input variables:
    INTEGER,                    INTENT(IN)        :: ncid
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_grid),            INTENT(INOUT)     :: grid
    REAL(dp), DIMENSION(:    ), INTENT(IN)        :: d
    INTEGER,                    INTENT(IN)        :: id_var
    INTEGER,                    INTENT(IN)        :: ti
    
    ! Local variables:
    INTEGER                                       :: int_dummy
    REAL(dp), DIMENSION(:    ), POINTER           :: d_square
    REAL(dp), DIMENSION(:,:  ), POINTER           :: d_grid
    INTEGER                                       :: wd_square, wd_grid
    
    int_dummy = mesh%nV
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( grid%mesh%nV,     d_square, wd_square)
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid,   wd_grid  )
    
    ! Map data from the model mesh to the square mesh
    CALL remap_cons_1st_order_2D(       grid%mesh, grid%map_mesh2square%conservative, d, d_square)
   ! CALL remap_cons_2nd_order_2D( mesh, grid%mesh, grid%map_mesh2square%conservative, d, d_square)
    
    ! Map data from the square mesh to the grid
    CALL map_square_to_grid_2D( grid%mesh, grid%map_square2grid, d_square, d_grid)
    
    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, ti/) ))
    
    ! Deallocate shared memory
    CALL deallocate_shared( wd_square)
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE map_and_write_to_grid_netcdf_dp_2D
  SUBROUTINE map_and_write_to_grid_netcdf_dp_3D( ncid, mesh, grid, d, id_var, ti, nz)
  
    USE netcdf_module, ONLY: open_netcdf_file, handle_error, close_netcdf_file
    USE netcdf,        ONLY: nf90_put_var
   
    IMPLICIT NONE
    
    ! Input variables:
    INTEGER,                    INTENT(IN)        :: ncid
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_grid),            INTENT(INOUT)     :: grid
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: d
    INTEGER,                    INTENT(IN)        :: id_var
    INTEGER,                    INTENT(IN)        :: ti
    INTEGER,                    INTENT(IN)        :: nz
    
    ! Local variables:
    INTEGER                                       :: int_dummy
    REAL(dp), DIMENSION(:,:  ), POINTER           :: d_square
    REAL(dp), DIMENSION(:,:,:), POINTER           :: d_grid
    INTEGER                                       :: wd_square, wd_grid
    
    int_dummy = mesh%nV
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%mesh%nV,     nz, d_square, wd_square)
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, nz, d_grid,   wd_grid  )
    
    ! Map data from the model mesh to the square mesh
    CALL remap_cons_1st_order_3D(       grid%mesh, grid%map_mesh2square%conservative, d, d_square, nz)
  !  CALL remap_cons_2nd_order_3D( mesh, grid%mesh, grid%map_mesh2square%conservative, d, d_square, nz)
    
    ! Map data from the square mesh to the grid
    CALL map_square_to_grid_3D( grid%mesh, grid%map_square2grid, d_square, d_grid)
    
    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, 1, ti/) ))
    
    ! Deallocate shared memory
    CALL deallocate_shared( wd_square)
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE map_and_write_to_grid_netcdf_dp_3D

END MODULE mesh_rectangular_module
