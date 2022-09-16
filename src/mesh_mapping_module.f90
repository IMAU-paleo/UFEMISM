MODULE mesh_mapping_module

  ! Routines for creating mapping arrays and mapping data between meshes and grids

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE petscksp
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr, MatDestroy
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
  USE data_types_module,               ONLY: type_mesh, type_remapping_mesh_mesh, type_grid, type_grid_lonlat, &
                                             type_remapping_lonlat2mesh, type_single_row_mapping_matrices, type_sparse_matrix_CSR_dp
  USE mesh_help_functions_module,      ONLY: is_in_triangle, write_mesh_to_text_file, lies_on_line_segment, segment_intersection, &
                                             find_containing_vertex, find_containing_triangle, is_in_Voronoi_cell, &
                                             find_shared_Voronoi_boundary, cross2, find_Voronoi_cell_vertices, find_triangle_area
  USE utilities_module,                ONLY: line_integral_xdy, line_integral_mxydx, line_integral_xydy, smooth_Gaussian_2D_grid, &
                                             smooth_Gaussian_3D_grid
  USE petsc_module,                    ONLY: multiply_PETSc_matrix_with_vector_1D, multiply_PETSc_matrix_with_vector_2D, mat_CSR2petsc
  USE mesh_operators_module,           ONLY: calc_matrix_operators_grid, apply_Neumann_BC_direct_a_2D, apply_Neumann_BC_direct_a_3D
  USE sparse_matrix_module,            ONLY: allocate_matrix_CSR_dist, add_entry_CSR_dist, finalise_matrix_CSR_dist, deallocate_matrix_CSR

  IMPLICIT NONE

  LOGICAL, PARAMETER :: do_check_matrices = .TRUE.

CONTAINS

! == Mapping data between a grid and a mesh
  SUBROUTINE map_grid2mesh_2D( grid, mesh, d_grid, d_mesh)
    ! Map a 2-D data field from the grid to the mesh using 2nd-order conservative remapping.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_grid2mesh_2D'
    INTEGER                                            :: n1,n2,n,i,j
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_grid_vec
    INTEGER                                            :: wd_grid_vec

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh,1) /= mesh%nV .OR. SIZE( d_grid,1) /= grid%nx .OR. SIZE( d_grid,2) /= grid%ny) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( grid%n, d_grid_vec, wd_grid_vec)

    ! Reshape data from vector form to grid form
    CALL partition_list( grid%n, par%i, par%n, n1, n2)
    DO n = n1, n2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      d_grid_vec( n) = d_grid( i,j)
    END DO
    CALL sync

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( grid%M_map_grid2mesh, d_grid_vec, d_mesh)

    ! Fix border elements because the remapping often is inaccurate there
    CALL apply_Neumann_BC_direct_a_2D( mesh, d_mesh)

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid_vec)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_grid2mesh_2D
  SUBROUTINE map_grid2mesh_3D( grid, mesh, d_grid, d_mesh)
    ! Map a 3-D data field from the grid to the mesh using 2nd-order conservative remapping.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: d_grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_grid2mesh_3D'
    INTEGER                                            :: n1,n2,n,i,j,nz
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_grid_vec
    INTEGER                                            :: wd_grid_vec

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh,1) /= mesh%nV .OR. SIZE( d_grid,1) /= grid%nx .OR. SIZE( d_grid,2) /= grid%ny .OR. SIZE( d_grid,3) /= SIZE( d_mesh,2)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    nz = SIZE( d_mesh,2)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%n, nz, d_grid_vec, wd_grid_vec)

    ! Reshape data from vector form to grid form
    CALL partition_list( grid%n, par%i, par%n, n1, n2)
    DO n = n1, n2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      d_grid_vec( n,:) = d_grid( i,j,:)
    END DO
    CALL sync

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( grid%M_map_grid2mesh, d_grid_vec, d_mesh)

    ! Fix border elements because the remapping often is inaccurate there
    CALL apply_Neumann_BC_direct_a_3D( mesh, d_mesh)

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid_vec)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_grid2mesh_3D
  SUBROUTINE map_mesh2grid_2D( mesh, grid, d_mesh, d_grid)
    ! Map a 2-D data field from the mesh to the grid using 2nd-order conservative remapping.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_mesh2grid_2D'
    INTEGER                                            :: n1,n2,n,i,j
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_grid_vec
    INTEGER                                            :: wd_grid_vec

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh,1) /= mesh%nV .OR. SIZE( d_grid,1) /= grid%nx .OR. SIZE( d_grid,2) /= grid%ny) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( grid%n, d_grid_vec, wd_grid_vec)

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_1D( grid%M_map_mesh2grid, d_mesh, d_grid_vec)

    ! Reshape data from vector form to grid form
    CALL partition_list( grid%n, par%i, par%n, n1, n2)
    DO n = n1, n2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      d_grid( i,j) = d_grid_vec( n)
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid_vec)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_mesh2grid_2D
  SUBROUTINE map_mesh2grid_3D( mesh, grid, d_mesh, d_grid)
    ! Map a 3-D data field from the mesh to the grid using 2nd-order conservative remapping.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_mesh
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_mesh2grid_3D'
    INTEGER                                            :: n1,n2,n,i,j,nz
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_grid_vec
    INTEGER                                            :: wd_grid_vec

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh,1) /= mesh%nV .OR. SIZE( d_grid,1) /= grid%nx .OR. SIZE( d_grid,2) /= grid%ny .OR. SIZE( d_mesh,2) /= SIZE( d_grid,3)) THEN
      CALL crash('data fields are the wrong size!')
    END IF

    nz = SIZE( d_grid,3)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%n, nz, d_grid_vec, wd_grid_vec)

    ! Perform the mapping operation as a matrix multiplication
    CALL multiply_PETSc_matrix_with_vector_2D( grid%M_map_mesh2grid, d_mesh, d_grid_vec)

    ! Reshape data from vector form to grid form
    CALL partition_list( grid%n, par%i, par%n, n1, n2)
    DO n = n1, n2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      d_grid( i,j,:) = d_grid_vec( n,:)
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid_vec)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_mesh2grid_3D

! == Subroutine for mapping data from a lat/lon-grid and the mesh ==
  SUBROUTINE create_remapping_arrays_lonlat_mesh( mesh, grid_lonlat, map)
    ! Create remapping arrays for remapping data from a global lon/lat-grid to the model mesh
    ! using bilinear interpolation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid_lonlat),              INTENT(IN)    :: grid_lonlat
    TYPE(type_remapping_lonlat2mesh),    INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_remapping_arrays_lonlat_mesh'
    INTEGER                                            :: vi
    INTEGER                                            :: il,iu,jl,ju
    REAL(dp)                                           :: wil,wiu,wjl,wju

    ! Add routine to path
    CALL init_routine( routine_name)

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
      il  = MAX(1, MIN( grid_lonlat%nlon-1, 1 + FLOOR((mesh%lon( vi) - MINVAL(grid_lonlat%lon)) / (grid_lonlat%lon(2) - grid_lonlat%lon(1)))))
      iu  = il + 1
      wil = (grid_lonlat%lon(iu) - mesh%lon( vi)) / (grid_lonlat%lon(2) - grid_lonlat%lon(1))
      wiu = 1._dp - wil

      ! Exception for pixels near the zero meridian
      IF (mesh%lon( vi) < MINVAL(grid_lonlat%lon)) THEN
        il  = grid_lonlat%nlon
        iu  = 1
        wil = (grid_lonlat%lon( iu) - mesh%lon( vi)) / (grid_lonlat%lon(2) - grid_lonlat%lon(1))
        wiu = 1._dp - wil
      ELSEIF (mesh%lon( vi) > MAXVAL(grid_lonlat%lon)) THEN
        il  = grid_lonlat%nlon
        iu  = 1
        wiu = (mesh%lon( vi) - grid_lonlat%lon( il)) / (grid_lonlat%lon(2) - grid_lonlat%lon(1))
        wil = 1._dp - wiu
      END IF

      jl  = MAX(1, MIN( grid_lonlat%nlat-1, 1 + FLOOR((mesh%lat( vi) - MINVAL(grid_lonlat%lat)) / (grid_lonlat%lat(2) - grid_lonlat%lat(1)))))
      ju  = jl + 1
      wjl = (grid_lonlat%lat( ju) - mesh%lat( vi)) / (grid_lonlat%lat(2) - grid_lonlat%lat(1))
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

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=8)

  END SUBROUTINE create_remapping_arrays_lonlat_mesh
  SUBROUTINE map_lonlat2mesh_2D( mesh, map, d_grid, d_mesh)
    ! Map data from a lon/lat-grid to the model mesh using bilinear interpolation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_remapping_lonlat2mesh),    INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_grid
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_lonlat2mesh_2D'
    INTEGER                                            :: vi
    INTEGER                                            :: il,iu,jl,ju
    REAL(dp)                                           :: wil,wiu,wjl,wju

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_lonlat2mesh_2D
  SUBROUTINE map_lonlat2mesh_3D( mesh, map, d_grid, d_mesh)
    ! Map data from a lon/lat-grid to the model mesh using bilinear interpolation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_remapping_lonlat2mesh),    INTENT(IN)    :: map
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: d_grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_lonlat2mesh_3D'
    INTEGER                                            :: vi
    INTEGER                                            :: il,iu,jl,ju
    REAL(dp)                                           :: wil,wiu,wjl,wju

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_lonlat2mesh_3D
  SUBROUTINE deallocate_remapping_arrays_lonlat_mesh( map)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_remapping_lonlat2mesh),    INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_remapping_arrays_lonlat_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL deallocate_shared( map%wilat1)
    CALL deallocate_shared( map%wilat2)
    CALL deallocate_shared( map%wilon1)
    CALL deallocate_shared( map%wilon2)
    CALL deallocate_shared( map%wwlat1)
    CALL deallocate_shared( map%wwlat2)
    CALL deallocate_shared( map%wwlon1)
    CALL deallocate_shared( map%wwlon2)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_remapping_arrays_lonlat_mesh

! == Medium-level remapping routines, where the memory containing the data is reallocated
  SUBROUTINE remap_field_dp_2D( mesh_src, mesh_dst, map, d, w, method)
    ! Remap a 2-D data field from mesh_src to mesh_dst using the specified remapping method. Includes memory reallocation.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map          ! Remapping matrices
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(INOUT) :: d            ! Pointer to the data
    INTEGER,                             INTENT(INOUT) :: w            ! MPI window to the shared memory space containing that data
    CHARACTER(LEN=*),                    INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_field_dp_2D'
    REAL(dp), DIMENSION(:    ), POINTER                :: d_temp
    INTEGER                                            :: w_temp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh_src%nV) THEN
      CALL crash('data field is the wrong size!')
    END IF

    ! Allocate temporary memory
    CALL allocate_shared_dp_1D( mesh_src%nV, d_temp, w_temp)

    ! Copy data to temporary memory
    d_temp( mesh_src%vi1: mesh_src%vi2) = d( mesh_src%vi1: mesh_src%vi2)

    ! Deallocate memory
    CALL deallocate_shared( w)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh_dst%nV, d, w)

    ! Map data field from source mesh to new mesh
    IF     (method == 'trilin') THEN
      CALL multiply_PETSc_matrix_with_vector_1D( map%M_trilin,            d_temp, d)
    ELSEIF (method == 'nearest_neighbour') THEN
      CALL multiply_PETSc_matrix_with_vector_1D( map%M_nearest_neighbour, d_temp, d)
    ELSEIF (method == 'cons_1st_order') THEN
      CALL multiply_PETSc_matrix_with_vector_1D( map%M_cons_1st_order,    d_temp, d)
    ELSEIF (method == 'cons_2nd_order') THEN
      CALL multiply_PETSc_matrix_with_vector_1D( map%M_cons_2nd_order,    d_temp, d)
    ELSE
      CALL crash('unknown remapping method "' // TRIM( method) // '"!')
    END IF

    ! Deallocate temporary memory
    CALL deallocate_shared( w_temp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_field_dp_2D
  SUBROUTINE remap_field_dp_3D( mesh_src, mesh_dst, map, d, w, method)
    ! Remap a 3-D data field from mesh_src to mesh_dst using the specified remapping method. Includes memory reallocation.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_src
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_dst
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map          ! Remapping matrices
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: d            ! Pointer to the data
    INTEGER,                             INTENT(INOUT) :: w            ! MPI window to the shared memory space containing that data
    CHARACTER(LEN=*),                    INTENT(IN)    :: method

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_field_dp_3D'
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_temp
    INTEGER                                            :: w_temp, nz

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh_src%nV) THEN
      CALL crash('data field is the wrong size!')
    END IF

    nz = SIZE( d,2)

    ! Allocate temporary memory
    CALL allocate_shared_dp_2D( mesh_src%nV, nz, d_temp, w_temp)

    ! Copy data to temporary memory
    d_temp( mesh_src%vi1: mesh_src%vi2,:) = d( mesh_src%vi1: mesh_src%vi2,:)

    ! Deallocate memory
    CALL deallocate_shared( w)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh_dst%nV, nz, d, w)

    ! Map data field from source mesh to new mesh
    IF     (method == 'trilin') THEN
      CALL multiply_PETSc_matrix_with_vector_2D( map%M_trilin,            d_temp, d)
    ELSEIF (method == 'nearest_neighbour') THEN
      CALL multiply_PETSc_matrix_with_vector_2D( map%M_nearest_neighbour, d_temp, d)
    ELSEIF (method == 'cons_1st_order') THEN
      CALL multiply_PETSc_matrix_with_vector_2D( map%M_cons_1st_order,    d_temp, d)
    ELSEIF (method == 'cons_2nd_order') THEN
      CALL multiply_PETSc_matrix_with_vector_2D( map%M_cons_2nd_order,    d_temp, d)
    ELSE
      CALL crash('data fields are the wrong size!')
    END IF

    ! Deallocate temporary memory
    CALL deallocate_shared( w_temp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_field_dp_3D

! == Smoothing operations on the mesh
  SUBROUTINE smooth_Gaussian_2D( mesh, grid, d_mesh, r)
    ! Use 2nd-order conservative remapping to map the 2-D data from the mesh
    ! to the square grid. Apply the smoothing on the gridded data, then map
    ! it back to the mesh. The numerical diffusion arising from the two mapping
    ! operations is not a problem since we're smoothing the data anyway.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_mesh
    REAL(dp),                            INTENT(IN)    :: r

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Gaussian_2D'
    REAL(dp), DIMENSION(:,:  ), POINTER                :: d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh,1) /= mesh%nV) THEN
      CALL crash('data field is the wrong size!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid, wd_grid)

    ! Map data to the grid
    CALL map_mesh2grid_2D( mesh, grid, d_mesh, d_grid)

    ! Apply smoothing on the gridded data
    CALL smooth_Gaussian_2D_grid( grid, d_grid, r)

    ! Map data back to the mesh
    CALL map_grid2mesh_2D( grid, mesh, d_grid, d_mesh)

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE smooth_Gaussian_2D
  SUBROUTINE smooth_Gaussian_3D( mesh, grid, d_mesh, r)
    ! Use 2nd-order conservative remapping to map the 3-D data from the mesh
    ! to the square grid. Apply the smoothing on the gridded data, then map
    ! it back to the mesh. The numerical diffusion arising from the two mapping
    ! operations is not a problem since we're smoothing the data anyway.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d_mesh
    REAL(dp),                            INTENT(IN)    :: r

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Gaussian_3D'
    REAL(dp), DIMENSION(:,:,:), POINTER                :: d_grid
    INTEGER                                            :: wd_grid, nz

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d_mesh,1) /= mesh%nV) THEN
      CALL crash('data field is the wrong size!')
    END IF

    nz = SIZE( d_mesh,2)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, nz, d_grid, wd_grid)

    ! Map data to the grid
    CALL map_mesh2grid_3D( mesh, grid, d_mesh, d_grid)

    ! Apply smoothing on the gridded data
    CALL smooth_Gaussian_3D_grid( grid, d_grid, r)

    ! Map data back to the mesh
    CALL map_grid2mesh_3D( grid, mesh, d_grid, d_mesh)

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE smooth_Gaussian_3D

! == Calculating the remapping matrices
  SUBROUTINE calc_remapping_operator_grid2mesh( grid, mesh)
    ! Calculate the remapping operators from the square grid to the mesh using 2nd-order conservative remapping
    !
    ! NOTE: the current implementation is a compromise. For "small" triangles (defined as having an area smaller
    !       than four times that of a square grid cell), a 2nd-order conservative remapping operation is calculated
    !       explicitly, using the line integrals around area of overlap. However, for "large" triangles (defined as
    !       all the rest), the result is generally very close to simply averaging over all the overlapping grid cells.
    !       Explicitly calculating the line integrals around all the grid cells is very slow, so this
    !       seems like a reasonable compromise.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_remapping_operator_grid2mesh'
    TYPE(PetscErrorCode)                               :: perr
    LOGICAL                                            :: count_coincidences
    INTEGER                                            :: nrows_A, ncols_A, nnz_est, nnz_est_proc, nnz_per_row_max
    TYPE(type_sparse_matrix_CSR_dp)                    :: A_xdy_a_g_CSR, A_mxydx_a_g_CSR, A_xydy_a_g_CSR
    TYPE(tMat)                                         :: A_xdy_a_g    , A_mxydx_a_g    , A_xydy_a_g
    INTEGER                                            :: vi1, vi2, vi
    INTEGER                                            :: nVor, vori1, vori2
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: Vor
    REAL(dp), DIMENSION(2)                             :: p, q
    INTEGER                                            :: k, n, i, j, kk, vj
    REAL(dp)                                           :: xl, xu, yl, yu
    REAL(dp), DIMENSION(2)                             :: sw, se, nw, ne
    INTEGER                                            :: vi_hint
    REAL(dp)                                           :: xmin, xmax, ymin, ymax
    INTEGER                                            :: il, iu, jl, ju
    TYPE(type_single_row_mapping_matrices)             :: single_row_Vor, single_row_grid
    TYPE(tMat)                                         :: w0, w1x, w1y
    INTEGER                                            :: ncols
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: cols
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: vals, w0_row, w1x_row, w1y_row
    REAL(dp)                                           :: A_overlap_tot
    TYPE(tMat)                                         :: grid_M_ddx, grid_M_ddy
    TYPE(tMat)                                         :: M1, M2

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
  ! ===========================================================================

    ! Matrix sise
    nrows_A         = mesh%nV  ! to
    ncols_A         = grid%n   ! from
    nnz_est         = 4 * MAX( nrows_A, ncols_A)
    nnz_est_proc    = CEILING( REAL( nnz_est, dp) / REAL( par%n, dp))
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh%A) / (grid%dx**2)), &
                                    CEILING( 2._dp * (grid%dx**2) / MINVAL( mesh%A))) )

    CALL allocate_matrix_CSR_dist( A_xdy_a_g_CSR  , nrows_A, ncols_A, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( A_mxydx_a_g_CSR, nrows_A, ncols_A, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( A_xydy_a_g_CSR , nrows_A, ncols_A, nnz_est_proc)

    ! Allocate memory for single row results
    single_row_Vor%n_max = nnz_per_row_max
    single_row_Vor%n     = 0
    ALLOCATE( single_row_Vor%index_left( single_row_Vor%n_max))
    ALLOCATE( single_row_Vor%LI_xdy(     single_row_Vor%n_max))
    ALLOCATE( single_row_Vor%LI_mxydx(   single_row_Vor%n_max))
    ALLOCATE( single_row_Vor%LI_xydy(    single_row_Vor%n_max))

    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    ALLOCATE( single_row_grid%index_left( single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_xdy(     single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_xydy(    single_row_grid%n_max))

    ALLOCATE( Vor( mesh%nC_mem+2,2))

    !IF (par%master) WRITE(0,*) 'calc_remapping_operator_grid2mesh - calculating all the line integrals...'

    ! Calculate line integrals around all Voronoi cells
    DO vi = mesh%vi1, mesh%vi2

      IF (mesh%A( vi) < 4._dp * grid%dx**2) THEN
        ! This Voronoi cell is small enough to warrant a proper line integral

        ! Clean up single row results
        single_row_Vor%n          = 0
        single_row_Vor%index_left = 0
        single_row_Vor%LI_xdy     = 0._dp
        single_row_Vor%LI_mxydx   = 0._dp
        single_row_Vor%LI_xydy    = 0._dp

        ! Integrate around the complete Voronoi cell boundary
        CALL find_Voronoi_cell_vertices( mesh, vi, Vor, nVor)
        IF (mesh%edge_index( vi) > 0) THEN
          Vor( nVor+1,:) = Vor( 1,:)
          nVor = nVor + 1
        END IF
        DO vori1 = 1, nVor-1
          vori2 = vori1 + 1
          IF (vori2 > nVor) vori2 = 1
          p = Vor( vori1,:)
          q = Vor( vori2,:)
          count_coincidences = .TRUE.
          CALL trace_line_grid( grid, p, q, single_row_Vor, count_coincidences)
        END DO

        ! Next integrate around the grid cells overlapping with this triangle
        DO k = 1, single_row_Vor%n

          ! Clean up single row results
          single_row_grid%n          = 0
          single_row_grid%index_left = 0
          single_row_grid%LI_xdy     = 0._dp
          single_row_grid%LI_mxydx   = 0._dp
          single_row_grid%LI_xydy    = 0._dp

          ! The grid cell
          n  = single_row_Vor%index_left( k)
          i  = grid%n2ij( n,1)
          j  = grid%n2ij( n,2)

          xl = grid%x( i) - grid%dx / 2._dp
          xu = grid%x( i) + grid%dx / 2._dp
          yl = grid%y( j) - grid%dx / 2._dp
          yu = grid%y( j) + grid%dx / 2._dp

          sw = [xl,yl]
          nw = [xl,yu]
          se = [xu,yl]
          ne = [xu,yu]

          ! Integrate around the grid cell
          vi_hint = vi
          count_coincidences = .FALSE.
          CALL trace_line_Vor( mesh, sw, se, single_row_grid, count_coincidences, vi_hint)
          CALL trace_line_Vor( mesh, se, ne, single_row_grid, count_coincidences, vi_hint)
          CALL trace_line_Vor( mesh, ne, nw, single_row_grid, count_coincidences, vi_hint)
          CALL trace_line_Vor( mesh, nw, sw, single_row_grid, count_coincidences, vi_hint)

          ! Add contribution for this particular triangle
          DO kk = 1, single_row_grid%n
            vj = single_row_grid%index_left( kk)
            IF (vj == vi) THEN
              ! Add contribution to this triangle
              single_row_Vor%LI_xdy(   k) = single_row_Vor%LI_xdy(   k) + single_row_grid%LI_xdy(   kk)
              single_row_Vor%LI_mxydx( k) = single_row_Vor%LI_mxydx( k) + single_row_grid%LI_mxydx( kk)
              single_row_Vor%LI_xydy(  k) = single_row_Vor%LI_xydy(  k) + single_row_grid%LI_xydy(  kk)
              EXIT
            END IF
          END DO ! DO kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          CALL add_entry_CSR_dist( A_xdy_a_g_CSR  , vi, n, single_row_Vor%LI_xdy(   k))
          CALL add_entry_CSR_dist( A_mxydx_a_g_CSR, vi, n, single_row_Vor%LI_mxydx( k))
          CALL add_entry_CSR_dist( A_xydy_a_g_CSR , vi, n, single_row_Vor%LI_xydy(  k))

        END DO ! DO k = 1, single_row_Vor%n

      ELSE ! IF (mesh%A( vi) < 4._dp * grid%dx**2) THEN
        ! This Voronoi cell is big enough that we can just average over the grid cells it contains

        ! Clean up single row results
        single_row_Vor%n = 0

        ! Find the square of grid cells enveloping this Voronoi cell
        CALL find_Voronoi_cell_vertices( mesh, vi, Vor, nVor)

        xmin = MINVAL( Vor( 1:nVor,1))
        xmax = MAXVAL( Vor( 1:nVor,1))
        ymin = MINVAL( Vor( 1:nVor,2))
        ymax = MAXVAL( Vor( 1:nVor,2))

        il = MAX( 1, MIN( grid%nx, 1 + FLOOR( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        iu = MAX( 1, MIN( grid%nx, 1 + FLOOR( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx) ))
        jl = MAX( 1, MIN( grid%ny, 1 + FLOOR( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx) ))
        ju = MAX( 1, MIN( grid%ny, 1 + FLOOR( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx) ))

        ! Check which of the grid cells in this square lie inside the triangle
        DO i = il, iu
        DO j = jl, ju

          n = grid%ij2n( i,j)
          p = [grid%x( i), grid%y( j)]

          IF (is_in_Voronoi_cell( mesh, p, vi)) THEN
            ! This grid cell lies inside the triangle; add it to the single row
            single_row_Vor%n = single_row_Vor%n + 1
            single_row_Vor%index_left( single_row_Vor%n) = n
            single_row_Vor%LI_xdy(     single_row_Vor%n) = grid%dx**2
            single_row_Vor%LI_mxydx(   single_row_Vor%n) = grid%x( i) * grid%dx**2
            single_row_Vor%LI_xydy(    single_row_Vor%n) = grid%y( j) * grid%dx**2
          END IF

        END DO
        END DO

        ! Add entries to the big matrices
        DO k = 1, single_row_Vor%n
          n = single_row_Vor%index_left( k)
          CALL add_entry_CSR_dist( A_xdy_a_g_CSR  , vi, n, single_row_Vor%LI_xdy(   k))
          CALL add_entry_CSR_dist( A_mxydx_a_g_CSR, vi, n, single_row_Vor%LI_mxydx( k))
          CALL add_entry_CSR_dist( A_xydy_a_g_CSR , vi, n, single_row_Vor%LI_xydy(  k))
        END DO

      END IF ! IF (mesh%A( vi) < 4._dp * grid%dx**2) THEN

    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( single_row_Vor%index_left )
    DEALLOCATE( single_row_Vor%LI_xdy     )
    DEALLOCATE( single_row_Vor%LI_mxydx   )
    DEALLOCATE( single_row_Vor%LI_xydy    )

    DEALLOCATE( single_row_grid%index_left )
    DEALLOCATE( single_row_grid%LI_xdy     )
    DEALLOCATE( single_row_grid%LI_mxydx   )
    DEALLOCATE( single_row_grid%LI_xydy    )

    DEALLOCATE( Vor)

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( A_xdy_a_g_CSR  , mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( A_mxydx_a_g_CSR, mesh%vi1, mesh%vi2)
    CALL finalise_matrix_CSR_dist( A_xydy_a_g_CSR , mesh%vi1, mesh%vi2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( A_xdy_a_g_CSR  , A_xdy_a_g  )
    CALL mat_CSR2petsc( A_mxydx_a_g_CSR, A_mxydx_a_g)
    CALL mat_CSR2petsc( A_xydy_a_g_CSR , A_xydy_a_g )

    ! Clean up the Fortran versions
    CALL deallocate_matrix_CSR( A_xdy_a_g_CSR  )
    CALL deallocate_matrix_CSR( A_mxydx_a_g_CSR)
    CALL deallocate_matrix_CSR( A_xydy_a_g_CSR )

  ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
  ! ==============================================================

    !IF (par%master) WRITE(0,*) 'calc_remapping_operator_grid2mesh - calculating w0, w1x, w1y...'

    CALL MatDuplicate( A_xdy_a_g, MAT_SHARE_NONZERO_PATTERN, w0 , perr)
    CALL MatDuplicate( A_xdy_a_g, MAT_SHARE_NONZERO_PATTERN, w1x, perr)
    CALL MatDuplicate( A_xdy_a_g, MAT_SHARE_NONZERO_PATTERN, w1y, perr)

    ALLOCATE( cols(    nnz_per_row_max))
    ALLOCATE( vals(    nnz_per_row_max))
    ALLOCATE( w0_row(  nnz_per_row_max))
    ALLOCATE( w1x_row( nnz_per_row_max))
    ALLOCATE( w1y_row( nnz_per_row_max))

    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( A_xdy_a_g, vi1, vi2, perr)

    DO vi = vi1+1, vi2 ! +1 because PETSc indexes from 0

      ! w0
      CALL MatGetRow( A_xdy_a_g, vi-1, ncols, cols, vals, perr)
      A_overlap_tot = SUM( vals( 1:ncols))
      DO k = 1, ncols
        w0_row( k) = vals( k) / A_overlap_tot
        CALL MatSetValues( w0, 1, vi-1, 1, cols( k), w0_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( A_xdy_a_g, vi-1, ncols, cols, vals, perr)

      ! w1x
      CALL MatGetRow( A_mxydx_a_g, vi-1, ncols, cols, vals, perr)
      DO k = 1, ncols
        n = cols( k)+1
        i = grid%n2ij( n,1)
        j = grid%n2ij( n,2)
        w1x_row( k) = (vals( k) / A_overlap_tot) - (grid%x( i) * w0_row( k))
        CALL MatSetValues( w1x, 1, vi-1, 1, cols( k), w1x_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( A_mxydx_a_g, vi-1, ncols, cols, vals, perr)

      ! w1y
      CALL MatGetRow( A_xydy_a_g, vi-1, ncols, cols, vals, perr)
      DO k = 1, ncols
        n = cols( k)+1
        i = grid%n2ij( n,1)
        j = grid%n2ij( n,2)
        w1y_row( k) = (vals( k) / A_overlap_tot) - (grid%y( j) * w0_row( k))
        CALL MatSetValues( w1y, 1, vi-1, 1, cols( k), w1y_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( A_xydy_a_g, vi-1, ncols, cols, vals, perr)

    END DO
    CALL sync

    CALL MatDestroy( A_xdy_a_g  , perr)
    CALL MatDestroy( A_mxydx_a_g, perr)
    CALL MatDestroy( A_xydy_a_g , perr)

    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.

    CALL MatAssemblyBegin( w0 , MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( w1x, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( w1y, MAT_FINAL_ASSEMBLY, perr)

    CALL MatAssemblyEnd(   w0 , MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   w1x, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   w1y, MAT_FINAL_ASSEMBLY, perr)

    ! Calculate the remapping matrix

    !IF (par%master) WRITE(0,*) 'calc_remapping_operator_grid2mesh - calculating remapping matrix...'

    CALL calc_matrix_operators_grid( grid, grid_M_ddx, grid_M_ddy)

    CALL MatDuplicate( w0, MAT_COPY_VALUES, grid%M_map_grid2mesh, perr)
    CALL MatMatMult( w1x, grid_M_ddx, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    CALL MatMatMult( w1y, grid_M_ddy, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M2, perr)

    CALL MatDestroy( grid_M_ddx    , perr)
    CALL MatDestroy( grid_M_ddy    , perr)
    CALL MatDestroy( w0            , perr)
    CALL MatDestroy( w1x           , perr)
    CALL MatDestroy( w1y           , perr)

    CALL MatAXPY( grid%M_map_grid2mesh, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( grid%M_map_grid2mesh, 1._dp, M2, DIFFERENT_NONZERO_PATTERN, perr)

    CALL MatDestroy( M1, perr)
    CALL MatDestroy( M2, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_remapping_operator_grid2mesh
  SUBROUTINE calc_remapping_operator_mesh2grid( mesh, grid)
    ! Calculate the remapping operators from the mesh to the square grid using 2nd-order conservative remapping
    !
    ! NOTE: the current implementation is a compromise. For "small" triangles (defined as being having an area smaller
    !       than four times that of a square grid cell), a 2nd-order conservative remapping operation is calculated
    !       explicitly, using the line integrals around area of overlap. However, for "large" triangles (defined as
    !       all the rest), the result is generally very close to simply averaging over all the overlapping grid cells.
    !       Explicitly calculating the line integrals around all the grid cells is prohibitively slow, so this
    !       seems like a reasonable compromise.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_grid),                     INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_remapping_operator_mesh2grid'
    TYPE(PetscErrorCode)                               :: perr
    LOGICAL                                            :: count_coincidences
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  overlaps_with_small_triangle,  containing_triangle
    INTEGER                                            :: woverlaps_with_small_triangle, wcontaining_triangle
    INTEGER                                            :: ti
    INTEGER                                            :: via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    REAL(dp)                                           :: xmin, xmax, ymin, ymax
    INTEGER                                            :: il, iu, jl, ju
    INTEGER                                            :: i, j, ii, jj
    INTEGER                                            :: nrows_A, ncols_A, nnz_est, nnz_est_proc, nnz_per_row_max
    TYPE(type_sparse_matrix_CSR_dp)                    :: A_xdy_g_b_CSR, A_mxydx_g_b_CSR, A_xydy_g_b_CSR
    TYPE(tMat)                                         :: A_xdy_g_b    , A_mxydx_g_b    , A_xydy_g_b
    TYPE(type_single_row_mapping_matrices)             :: single_row_grid, single_row_Tri
    INTEGER                                            :: n1, n2, n, ti_hint
    REAL(dp), DIMENSION(2)                             :: p
    REAL(dp)                                           :: xl, xu, yl, yu
    REAL(dp), DIMENSION(2)                             :: sw, se, nw, ne
    INTEGER                                            :: k, kk, nn
    REAL(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy
    TYPE(tMat)                                         :: w0, w1x, w1y
    INTEGER                                            :: ncols
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: cols
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: vals, w0_row, w1x_row, w1y_row
    REAL(dp)                                           :: A_overlap_tot
    TYPE(tMat)                                         :: M1, M2

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Find all grid cells that overlap with small triangles
  ! ========================================================

    !IF (par%master) WRITE(0,*) 'calc_remapping_operator_mesh2grid - finding all grid cells overlapping with small/big triangles...'

    CALL allocate_shared_int_2D( grid%nx, grid%ny, overlaps_with_small_triangle, woverlaps_with_small_triangle)
    CALL allocate_shared_int_2D( grid%nx, grid%ny, containing_triangle         , wcontaining_triangle         )

    DO ti = mesh%ti1, mesh%ti2

      ! The three vertices spanning this triangle
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      pa  = mesh%V( via,:)
      pb  = mesh%V( vib,:)
      pc  = mesh%V( vic,:)

      ! The square enveloping this triangle
      xmin = MIN( MIN( pa(1), pb(1)), pc(1))
      xmax = MAX( MAX( pa(1), pb(1)), pc(1))
      ymin = MIN( MIN( pa(2), pb(2)), pc(2))
      ymax = MAX( MAX( pa(2), pb(2)), pc(2))

      ! The square of grid cells enveloping this triangle
      il = 1 + FLOOR( (xmin - grid%xmin + grid%dx / 2._dp) / grid%dx)
      iu = 1 + FLOOR( (xmax - grid%xmin + grid%dx / 2._dp) / grid%dx)
      jl = 1 + FLOOR( (ymin - grid%ymin + grid%dx / 2._dp) / grid%dx)
      ju = 1 + FLOOR( (ymax - grid%ymin + grid%dx / 2._dp) / grid%dx)

      il = MAX( 1      , il - 1)
      iu = MIN( grid%nx, iu + 1)
      jl = MAX( 1      , jl - 1)
      ju = MIN( grid%ny, ju + 1)

      IF (mesh%TriA( ti) < 4._dp * grid%dx**2) THEN
        ! This triangle is small; mark all grid cells it overlaps with

        ! Mark all these grid cells
        DO i = il, iu
        DO j = jl, ju
          overlaps_with_small_triangle( i,j) = 1
        END DO
        END DO

      ELSE
        ! This triangle is large; mark all grid cells it contains

        ! Mark all these grid cells
        DO i = il, iu
        DO j = jl, ju
          p = [grid%x( i), grid%y( j)]
          IF (is_in_triangle( pa, pb, pc, p)) THEN
            containing_triangle( i,j) = ti
          END IF
        END DO
        END DO

      END IF ! IF (mesh%TriA( ti) < 4._dp * grid%dx**2) THEN

    END DO
    CALL sync

    ! Treat grid cells that possibly were not yet marked before
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (containing_triangle( i,j) == 0 .AND. overlaps_with_small_triangle( i,j) == 0) THEN
        ! This grid cell does not overlap with a small triangle, but was not yet marked
        ! as being contained inside a large one; find the large triangle containing it.

        ! For efficiency, find the nearest grid cell that does list which large
        ! triangle contains it; use that as a hint for the triangle search
        n = 0
        ti_hint = 0
        DO WHILE (ti_hint == 0)
          n = n+1
          ! Safety
          IF (n > MAX( grid%nx, grid%ny)) EXIT
          il = MAX( 1      , i-n)
          iu = MIN( grid%nx, i+n)
          jl = MAX( 1      , j-n)
          ju = MIN( grid%ny, j+n)
          DO ii = il, iu
          DO jj = jl, ju
            IF (containing_triangle( ii,jj) > 0) THEN
              ti_hint = containing_triangle( ii,jj)
              EXIT
            END IF
          END DO
          IF (ti_hint > 0) EXIT
          END DO
        END DO
        IF (ti_hint == 0) ti_hint = 1

        ! Find the triangle containing this grid cell
        p = [MAX( mesh%xmin, MIN( mesh%xmax, grid%x( i) )), MAX( mesh%ymin, MIN( mesh%ymax, grid%y( j) ))]
        CALL find_containing_triangle( mesh, p, ti_hint)
        containing_triangle( i,j) = ti_hint

      END IF

    END DO
    END DO
    CALL sync

  ! == Integrate around all grid cells that overlap with small triangles
  ! ====================================================================

    !IF (par%master) WRITE(0,*) 'calc_remapping_operator_mesh2grid - calculating all line integrals...'

    ! Initialise the three matrices using the native UFEMISM CSR-matrix format

    ! Matrix sise
    nrows_A         = grid%n     ! to
    ncols_A         = mesh%nTri  ! from
    nnz_est         = 4 * MAX( nrows_A, ncols_A)
    nnz_est_proc    = CEILING( REAL( nnz_est, dp) / REAL( par%n, dp))
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh%TriA) / (grid%dx**2)), &
                                    CEILING( 2._dp * (grid%dx**2) / MINVAL( mesh%TriA))) )

    CALL allocate_matrix_CSR_dist( A_xdy_g_b_CSR  , nrows_A, ncols_A, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( A_mxydx_g_b_CSR, nrows_A, ncols_A, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( A_xydy_g_b_CSR , nrows_A, ncols_A, nnz_est_proc)

    ! Allocate memory for single row results
    single_row_grid%n_max = nnz_per_row_max
    single_row_grid%n     = 0
    ALLOCATE( single_row_grid%index_left( single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_xdy(     single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_mxydx(   single_row_grid%n_max))
    ALLOCATE( single_row_grid%LI_xydy(    single_row_grid%n_max))

    single_row_Tri%n_max = nnz_per_row_max
    single_row_Tri%n     = 0
    ALLOCATE( single_row_Tri%index_left( single_row_Tri%n_max))
    ALLOCATE( single_row_Tri%LI_xdy(     single_row_Tri%n_max))
    ALLOCATE( single_row_Tri%LI_mxydx(   single_row_Tri%n_max))
    ALLOCATE( single_row_Tri%LI_xydy(    single_row_Tri%n_max))

    ti_hint = 1

    CALL partition_list( grid%n, par%i, par%n, n1, n2)

    DO n = n1, n2

      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      p = [grid%x( i), grid%y( j)]

      IF (overlaps_with_small_triangle( i,j) == 1) THEN
        ! This grid cell overlaps with a small triangle; integrate around it, and around
        ! all triangles overlapping with it

        ! The four sides of the grid cell
        xl = grid%x( i) - grid%dx / 2._dp
        xu = grid%x( i) + grid%dx / 2._dp
        yl = grid%y( j) - grid%dx / 2._dp
        yu = grid%y( j) + grid%dx / 2._dp

        sw = [xl, yl]
        nw = [xl, yu]
        se = [xu, yl]
        ne = [xu, yu]

        ! Clear the single row results
        single_row_grid%n          = 0
        single_row_grid%index_left = 0
        single_row_grid%LI_xdy     = 0._dp
        single_row_grid%LI_mxydx   = 0._dp
        single_row_grid%LI_xydy    = 0._dp

        ! Integrate over all four sides
        count_coincidences = .TRUE.
        CALL trace_line_tri( mesh, sw, se, single_row_grid, count_coincidences, ti_hint)
        CALL trace_line_tri( mesh, se, ne, single_row_grid, count_coincidences, ti_hint)
        CALL trace_line_tri( mesh, ne, nw, single_row_grid, count_coincidences, ti_hint)
        CALL trace_line_tri( mesh, nw, sw, single_row_grid, count_coincidences, ti_hint)

        ! Next, integrate around all the triangles overlapping with this grid cell
        DO k = 1, single_row_grid%n

          ti = single_row_grid%index_left( k)

          ! The three vertices spanning this triangle
          via = mesh%Tri( ti,1)
          vib = mesh%Tri( ti,2)
          vic = mesh%Tri( ti,3)

          pa  = mesh%V( via,:)
          pb  = mesh%V( vib,:)
          pc  = mesh%V( vic,:)

          ! Clear the single row results
          single_row_Tri%n = 0
          single_row_Tri%index_left = 0
          single_row_Tri%LI_xdy     = 0._dp
          single_row_Tri%LI_mxydx   = 0._dp
          single_row_Tri%LI_xydy    = 0._dp

          ! Integrate over all three triangle sides
          count_coincidences = .FALSE.
          CALL trace_line_grid( grid, pa, pb, single_row_Tri, count_coincidences)
          CALL trace_line_grid( grid, pb, pc, single_row_Tri, count_coincidences)
          CALL trace_line_grid( grid, pc, pa, single_row_Tri, count_coincidences)

          ! Add contribution for this particular grid cell
          DO kk = 1, single_row_Tri%n
            nn = single_row_Tri%index_left( kk)
            IF (nn == n) THEN
              ! Add contribution to this triangle
              single_row_grid%LI_xdy(   k) = single_row_grid%LI_xdy(   k) + single_row_Tri%LI_xdy(   kk)
              single_row_grid%LI_mxydx( k) = single_row_grid%LI_mxydx( k) + single_row_Tri%LI_mxydx( kk)
              single_row_grid%LI_xydy(  k) = single_row_grid%LI_xydy(  k) + single_row_Tri%LI_xydy(  kk)
              EXIT
            END IF
          END DO ! DO kk = 1, single_row_grid%n

          ! Add entries to the big matrices
          CALL add_entry_CSR_dist( A_xdy_g_b_CSR  , n, ti, single_row_grid%LI_xdy(   k))
          CALL add_entry_CSR_dist( A_mxydx_g_b_CSR, n, ti, single_row_grid%LI_mxydx( k))
          CALL add_entry_CSR_dist( A_xydy_g_b_CSR , n, ti, single_row_grid%LI_xydy(  k))

        END DO ! DO k = 1, single_row_grid%n

      ELSE ! IF (overlaps_with_small_triangle( i,j) == 1) THEN
        ! This grid cell does not overlap with a small triangle; use only the
        ! contribution from the nearest triangle

        ti_hint = containing_triangle( i,j)

        LI_xdy   = grid%dx**2
        LI_mxydx = grid%dx**2 * grid%x( i)
        LI_xydy  = grid%dx**2 * grid%y( j)

        CALL add_entry_CSR_dist( A_xdy_g_b_CSR  , n, ti_hint, LI_xdy  )
        CALL add_entry_CSR_dist( A_mxydx_g_b_CSR, n, ti_hint, LI_mxydx)
        CALL add_entry_CSR_dist( A_xydy_g_b_CSR , n, ti_hint, LI_xydy )

      END IF ! IF (overlaps_with_small_triangle( i,j) == 1) THEN

    END DO ! DO n = n1, n2

    ! Clean up after yourself
    CALL deallocate_shared( woverlaps_with_small_triangle)
    CALL deallocate_shared( wcontaining_triangle         )

    DEALLOCATE( single_row_grid%index_left )
    DEALLOCATE( single_row_grid%LI_xdy     )
    DEALLOCATE( single_row_grid%LI_mxydx   )
    DEALLOCATE( single_row_grid%LI_xydy    )

    DEALLOCATE( single_row_Tri%index_left )
    DEALLOCATE( single_row_Tri%LI_xdy     )
    DEALLOCATE( single_row_Tri%LI_mxydx   )
    DEALLOCATE( single_row_Tri%LI_xydy    )

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( A_xdy_g_b_CSR  , n1, n2)
    CALL finalise_matrix_CSR_dist( A_mxydx_g_b_CSR, n1, n2)
    CALL finalise_matrix_CSR_dist( A_xydy_g_b_CSR , n1, n2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( A_xdy_g_b_CSR  , A_xdy_g_b  )
    CALL mat_CSR2petsc( A_mxydx_g_b_CSR, A_mxydx_g_b)
    CALL mat_CSR2petsc( A_xydy_g_b_CSR , A_xydy_g_b )

    ! Clean up the Fortran versions
    CALL deallocate_matrix_CSR( A_xdy_g_b_CSR  )
    CALL deallocate_matrix_CSR( A_mxydx_g_b_CSR)
    CALL deallocate_matrix_CSR( A_xydy_g_b_CSR )

  ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
  ! ==============================================================

    !IF (par%master) WRITE(0,*) 'calc_remapping_operator_mesh2grid - calculating w0, w1x, w1y...'

    CALL MatDuplicate( A_xdy_g_b, MAT_SHARE_NONZERO_PATTERN, w0 , perr)
    CALL MatDuplicate( A_xdy_g_b, MAT_SHARE_NONZERO_PATTERN, w1x, perr)
    CALL MatDuplicate( A_xdy_g_b, MAT_SHARE_NONZERO_PATTERN, w1y, perr)

    ALLOCATE( cols(    nnz_per_row_max))
    ALLOCATE( vals(    nnz_per_row_max))
    ALLOCATE( w0_row(  nnz_per_row_max))
    ALLOCATE( w1x_row( nnz_per_row_max))
    ALLOCATE( w1y_row( nnz_per_row_max))

    CALL MatGetOwnershipRange( A_xdy_g_b, n1, n2, perr)

    DO n = n1+1, n2 ! +1 because PETSc indexes from 0

      ! w0
      CALL MatGetRow( A_xdy_g_b, n-1, ncols, cols, vals, perr)
      A_overlap_tot = SUM( vals( 1:ncols))
      DO k = 1, ncols
        w0_row( k) = vals( k) / A_overlap_tot
        CALL MatSetValues( w0, 1, n-1, 1, cols( k), w0_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( A_xdy_g_b, n-1, ncols, cols, vals, perr)

      ! w1x
      CALL MatGetRow( A_mxydx_g_b, n-1, ncols, cols, vals, perr)
      DO k = 1, ncols
        ti = cols( k)+1
        w1x_row( k) = (vals( k) / A_overlap_tot) - (mesh%TriGC( ti,1) * w0_row( k))
        CALL MatSetValues( w1x, 1, n-1, 1, cols( k), w1x_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( A_mxydx_g_b, n-1, ncols, cols, vals, perr)

      ! w1y
      CALL MatGetRow( A_xydy_g_b, n-1, ncols, cols, vals, perr)
      DO k = 1, ncols
        ti = cols( k)+1
        w1y_row( k) = (vals( k) / A_overlap_tot) - (mesh%TriGC( ti,2) * w0_row( k))
        CALL MatSetValues( w1y, 1, n-1, 1, cols( k), w1y_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( A_xydy_g_b, n-1, ncols, cols, vals, perr)

    END DO
    CALL sync

    CALL MatDestroy( A_xdy_g_b  , perr)
    CALL MatDestroy( A_mxydx_g_b, perr)
    CALL MatDestroy( A_xydy_g_b , perr)

    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.

    CALL MatAssemblyBegin( w0 , MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( w1x, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyBegin( w1y, MAT_FINAL_ASSEMBLY, perr)

    CALL MatAssemblyEnd(   w0 , MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   w1x, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   w1y, MAT_FINAL_ASSEMBLY, perr)

    ! Calculate the remapping matrix

    !IF (par%master) WRITE(0,*) 'calc_remapping_operator_mesh2grid - calculating remapping matrix...'

    CALL MatMatMult( w0,  mesh%M_map_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, grid%M_map_mesh2grid, perr)
    CALL MatMatMult( w1x, mesh%M_ddx_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    CALL MatMatMult( w1y, mesh%M_ddy_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M2, perr)

    CALL MatDestroy( w0            , perr)
    CALL MatDestroy( w1x           , perr)
    CALL MatDestroy( w1y           , perr)

    CALL MatAXPY( grid%M_map_mesh2grid, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( grid%M_map_mesh2grid, 1._dp, M2, DIFFERENT_NONZERO_PATTERN, perr)

    CALL MatDestroy( M1, perr)
    CALL MatDestroy( M2, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_remapping_operator_mesh2grid
  SUBROUTINE calc_remapping_operators_mesh_mesh( mesh_src, mesh_dst, map)
    ! Calculate all the remapping operators between two meshes

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_src
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_dst
    TYPE(type_remapping_mesh_mesh),      INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_remapping_operators_mesh_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL calc_remapping_operators_mesh_mesh_trilin(            mesh_src, mesh_dst, map%M_trilin)
    CALL calc_remapping_operators_mesh_mesh_nearest_neighbour( mesh_src, mesh_dst, map%M_nearest_neighbour)
    CALL calc_remapping_operators_mesh_mesh_conservative(      mesh_src, mesh_dst, map%M_cons_1st_order, map%M_cons_2nd_order)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_remapping_operators_mesh_mesh
  SUBROUTINE calc_remapping_operators_mesh_mesh_trilin( mesh_src, mesh_dst, M_trilin)
    ! Calculate the trilinear interpolation operator from mesh_src to mesh_dst

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_src
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_dst
    TYPE(tMat),                          INTENT(INOUT) :: M_trilin

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_remapping_operators_mesh_mesh_trilin'
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER                                            :: vi_dst
    REAL(dp), DIMENSION(2)                             :: p
    INTEGER                                            :: ti_src, via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    REAL(dp)                                           :: Atri_abp, Atri_bcp, Atri_cap, Atri_abc, wa, wb, wc

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Use PETSc routines to initialise the matrix object
  ! =====================================================

    ! Matrix size
    nrows           = mesh_dst%nV   ! to
    ncols           = mesh_src%nV   ! from
    nnz_per_row_max = 3

    ! Initialise the matrix object
    CALL MatCreate( PETSC_COMM_WORLD, M_trilin, perr)

    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( M_trilin, 'mpiaij', perr)

    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( M_trilin, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)

    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( M_trilin, perr)

    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( M_trilin, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)

    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( M_trilin, istart, iend, perr)

    ! For all mesh_dst vertices, find the mesh_src triangle containing them
    ti_src = 1
    DO vi_dst = istart+1, iend ! +1 because PETSc indexes from 0

      p = mesh_dst%V( vi_dst,:)
      CALL find_containing_triangle( mesh_src, p, ti_src)

      ! Calculate the trilinear interpolation weights
      via = mesh_src%Tri( ti_src,1)
      vib = mesh_src%Tri( ti_src,2)
      vic = mesh_src%Tri( ti_src,3)

      pa  = mesh_src%V( via,:)
      pb  = mesh_src%V( vib,:)
      pc  = mesh_src%V( vic,:)

      CALL find_triangle_area( pa, pb, p, Atri_abp)
      CALL find_triangle_area( pb, pc, p, Atri_bcp)
      CALL find_triangle_area( pc, pa, p, Atri_cap)
      Atri_abc = Atri_abp + Atri_bcp + Atri_cap

      wa = Atri_bcp / Atri_abc
      wb = Atri_cap / Atri_abc
      wc = Atri_abp / Atri_abc

      ! Add to the matrix
      CALL MatSetValues( M_trilin, 1, vi_dst-1, 1, via-1, wa, INSERT_VALUES, perr)
      CALL MatSetValues( M_trilin, 1, vi_dst-1, 1, vib-1, wb, INSERT_VALUES, perr)
      CALL MatSetValues( M_trilin, 1, vi_dst-1, 1, vic-1, wc, INSERT_VALUES, perr)

    END DO ! DO vi_dst = mesh_dst%vi1, mesh_dst%vi2
    CALL sync

    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.

    CALL MatAssemblyBegin( M_trilin, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_trilin, MAT_FINAL_ASSEMBLY, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_remapping_operators_mesh_mesh_trilin
  SUBROUTINE calc_remapping_operators_mesh_mesh_nearest_neighbour( mesh_src, mesh_dst, M_nearest_neighbour)
    ! Calculate the nearest-neighbour interpolation operator from mesh_src to mesh_dst

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_src
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_dst
    TYPE(tMat),                          INTENT(INOUT) :: M_nearest_neighbour

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_remapping_operators_mesh_mesh_nearest_neighbour'
    INTEGER                                            :: ncols, nrows, nnz_per_row_max, istart, iend
    INTEGER                                            :: vi_dst
    REAL(dp), DIMENSION(2)                             :: p
    INTEGER                                            :: vi_src

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Use PETSc routines to initialise the matrix object
  ! =====================================================

    ! Matrix size
    nrows           = mesh_dst%nV   ! to
    ncols           = mesh_src%nV   ! from
    nnz_per_row_max = 1

    ! Initialise the matrix object
    CALL MatCreate( PETSC_COMM_WORLD, M_nearest_neighbour, perr)

    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( M_nearest_neighbour, 'mpiaij', perr)

    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( M_nearest_neighbour, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, perr)

    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( M_nearest_neighbour, perr)

    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( M_nearest_neighbour, nnz_per_row_max+1, PETSC_NULL_INTEGER, nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)

    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( M_nearest_neighbour, istart, iend, perr)

    ! For all mesh_dst vertices, find the mesh_src triangle containing them
    vi_src = 1
    DO vi_dst = istart+1, iend ! +1 because PETSc indexes from 0

      p = mesh_dst%V( vi_dst,:)
      CALL find_containing_vertex( mesh_src, p, vi_src)

      ! Add to the matrix
      CALL MatSetValues( M_nearest_neighbour, 1, vi_dst-1, 1, vi_src-1, 1._dp, INSERT_VALUES, perr)

    END DO ! DO vi_dst = istart+1, iend ! +1 because PETSc indexes from 0
    CALL sync

    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.

    CALL MatAssemblyBegin( M_nearest_neighbour, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   M_nearest_neighbour, MAT_FINAL_ASSEMBLY, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_remapping_operators_mesh_mesh_nearest_neighbour
  SUBROUTINE calc_remapping_operators_mesh_mesh_conservative( mesh_src, mesh_dst, M_cons_1st_order, M_cons_2nd_order)
    ! Calculate the 1st- and 2nd-order conservative remapping operators from mesh_src to mesh_dst

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_src
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_dst
    TYPE(tMat),                          INTENT(INOUT) :: M_cons_1st_order, M_cons_2nd_order

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_remapping_operators_mesh_mesh_conservative'
    TYPE(PetscErrorCode)                               :: perr
    LOGICAL                                            :: count_coincidences
    INTEGER                                            :: nnz_per_row_max
    TYPE(tMat)                                         :: B_xdy_b_a  , B_mxydx_b_a  , B_xydy_b_a
    TYPE(tMat)                                         :: B_xdy_a_b  , B_mxydx_a_b  , B_xydy_a_b
    TYPE(tMat)                                         :: B_xdy_b_a_T, B_mxydx_b_a_T, B_xydy_b_a_T
    TYPE(tMat)                                         :: w0, w1x, w1y
    INTEGER                                            :: istart, iend, n, k, ti
    INTEGER                                            :: ncols
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: cols
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: vals, w0_row, w1x_row, w1y_row
    REAL(dp)                                           :: A_overlap_tot
    TYPE(tMat)                                         :: M1, M2

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Integrate around the Voronoi cells of the destination mesh through the triangles of the source mesh
    count_coincidences = .TRUE.
    !IF (par%master) WRITE(0,*) 'calc_remapping_operators_mesh_mesh_conservative - integrating dst Voronoi cells through src triangles...'
    CALL integrate_Voronoi_cells_through_triangles( mesh_dst, mesh_src, B_xdy_a_b, B_mxydx_a_b, B_xydy_a_b, count_coincidences)

    ! Integrate around the triangles of the source mesh through the Voronoi cells of the destination mesh
    count_coincidences = .FALSE.
    !IF (par%master) WRITE(0,*) 'calc_remapping_operators_mesh_mesh_conservative - integrating src triangles through dst Voronoi cells...'
    CALL integrate_triangles_through_Voronoi_cells( mesh_src, mesh_dst, B_xdy_b_a, B_mxydx_b_a, B_xydy_b_a, count_coincidences)

    ! Transpose line integral matrices
    !IF (par%master) WRITE(0,*) 'calc_remapping_operators_mesh_mesh_conservative - transposing line integral matrices...'
    CALL MatCreateTranspose( B_xdy_b_a  , B_xdy_b_a_T  , perr)
    CALL MatCreateTranspose( B_mxydx_b_a, B_mxydx_b_a_T, perr)
    CALL MatCreateTranspose( B_xydy_b_a , B_xydy_b_a_T , perr)

    ! Combine line integrals around areas of overlap to get surface integrals over areas of overlap
    !IF (par%master) WRITE(0,*) 'calc_remapping_operators_mesh_mesh_conservative - combining line integrals around areas of overlap...'
    CALL MatAXPY( B_xdy_a_b  , 1._dp, B_xdy_b_a_T  , UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( B_mxydx_a_b, 1._dp, B_mxydx_b_a_T, UNKNOWN_NONZERO_PATTERN, perr)
    CALL MatAXPY( B_xydy_a_b , 1._dp, B_xydy_b_a_T , UNKNOWN_NONZERO_PATTERN, perr)

    CALL MatDestroy( B_xdy_b_a_T  , perr)
    CALL MatDestroy( B_mxydx_b_a_T, perr)
    CALL MatDestroy( B_xydy_b_a_T , perr)

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    !IF (par%master) WRITE(0,*) 'calc_remapping_operators_mesh_mesh_conservative - calculating remapping weights...'
    CALL MatDuplicate( B_xdy_a_b, MAT_SHARE_NONZERO_PATTERN, w0 , perr)
    CALL MatDuplicate( B_xdy_a_b, MAT_SHARE_NONZERO_PATTERN, w1x, perr)
    CALL MatDuplicate( B_xdy_a_b, MAT_SHARE_NONZERO_PATTERN, w1y, perr)

    ! Estimate maximum number of non-zeros per row (i.e. maximum number of grid cells overlapping with a mesh triangle)
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh_src%TriA) / MINVAL( mesh_dst%A   )), &
                                    CEILING( 2._dp * MAXVAL( mesh_dst%A   ) / MINVAL( mesh_src%TriA))) )

    ! Allocate memory for a single matrix row
    ALLOCATE( cols(    nnz_per_row_max))
    ALLOCATE( vals(    nnz_per_row_max))
    ALLOCATE( w0_row(  nnz_per_row_max))
    ALLOCATE( w1x_row( nnz_per_row_max))
    ALLOCATE( w1y_row( nnz_per_row_max))

    CALL MatGetOwnershipRange( B_xdy_a_b  , istart, iend, perr)

    DO n = istart+1, iend ! +1 because PETSc indexes from 0

      ! w0
      CALL MatGetRow( B_xdy_a_b, n-1, ncols, cols, vals, perr)
      A_overlap_tot = SUM( vals( 1:ncols))
      DO k = 1, ncols
        w0_row( k) = vals( k) / A_overlap_tot
        CALL MatSetValues( w0, 1, n-1, 1, cols( k), w0_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( B_xdy_a_b, n-1, ncols, cols, vals, perr)

      ! w1x
      CALL MatGetRow( B_mxydx_a_b, n-1, ncols, cols, vals, perr)
      DO k = 1, ncols
        ti = cols( k)+1
        w1x_row( k) = (vals( k) / A_overlap_tot) - (mesh_src%TriGC( ti,1) * w0_row( k))
        CALL MatSetValues( w1x, 1, n-1, 1, cols( k), w1x_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( B_mxydx_a_b, n-1, ncols, cols, vals, perr)

      ! w1y
      CALL MatGetRow( B_xydy_a_b, n-1, ncols, cols, vals, perr)
      DO k = 1, ncols
        ti = cols( k)+1
        w1y_row( k) = (vals( k) / A_overlap_tot) - (mesh_src%TriGC( ti,2) * w0_row( k))
        CALL MatSetValues( w1y, 1, n-1, 1, cols( k), w1y_row( k), INSERT_VALUES, perr)
      END DO
      CALL MatRestoreRow( B_xydy_a_b, n-1, ncols, cols, vals, perr)

    END DO
    CALL sync

    CALL MatDestroy( B_xdy_a_b  , perr)
    CALL MatDestroy( B_mxydx_a_b, perr)
    CALL MatDestroy( B_xydy_a_b , perr)

    ! Calculate the remapping matrices
    !IF (par%master) WRITE(0,*) 'calc_remapping_operators_mesh_mesh_conservative - combining weights into remapping matrix...'

    ! 1st-order = w0 * map_a_b
    CALL MatMatMult( w0 , mesh_src%M_map_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_cons_1st_order, perr)

    ! 2nd-order = 1st-order + w1x * ddx_a_b + w1y * ddy_a_b
    CALL MatDuplicate( M_cons_1st_order, MAT_COPY_VALUES, M_cons_2nd_order, perr)
    CALL MatMatMult( w1x, mesh_src%M_ddx_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M1, perr)  ! This can be done more efficiently now that the non-zero structure is known...
    CALL MatMatMult( w1y, mesh_src%M_ddy_a_b, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M2, perr)

    CALL MatDestroy( w0            , perr)
    CALL MatDestroy( w1x           , perr)
    CALL MatDestroy( w1y           , perr)

    CALL MatAXPY( M_cons_2nd_order, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, perr)
    CALL MatAXPY( M_cons_2nd_order, 1._dp, M2, DIFFERENT_NONZERO_PATTERN, perr)

    CALL MatDestroy( M1, perr)
    CALL MatDestroy( M2, perr)

    !IF (par%master) WRITE(0,*) 'calc_remapping_operators_mesh_mesh_conservative - done!'

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_remapping_operators_mesh_mesh_conservative

  ! Integrate around triangles/Voronoi cells through triangles/Voronoi cells
  SUBROUTINE integrate_triangles_through_Voronoi_cells( mesh_tri, mesh_Vor, B_xdy_b_a, B_mxydx_b_a, B_xydy_b_a, count_coincidences)
    ! Integrate around the triangles of mesh_tri through the Voronoi cells of mesh_Vor

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_tri
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_Vor
    TYPE(tMat),                          INTENT(INOUT) :: B_xdy_b_a
    TYPE(tMat),                          INTENT(INOUT) :: B_mxydx_b_a
    TYPE(tMat),                          INTENT(INOUT) :: B_xydy_b_a
    LOGICAL,                             INTENT(IN)    :: count_coincidences

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'integrate_triangles_through_Voronoi_cells'
    INTEGER                                            :: nrows_A, ncols_A, nnz_est, nnz_est_proc, nnz_per_row_max
    TYPE(type_sparse_matrix_CSR_dp)                    :: B_xdy_b_a_CSR, B_mxydx_b_a_CSR, B_xydy_b_a_CSR
    TYPE(type_single_row_mapping_matrices)             :: single_row
    INTEGER                                            :: via, vib, vic, ti, vi_hint, k
    REAL(dp), DIMENSION(2)                             :: p, q

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
  ! ===========================================================================

    ! Matrix sise
    nrows_A         = mesh_tri%nTri  ! to
    ncols_A         = mesh_Vor%nV    ! from
    nnz_est         = 4 * MAX( nrows_A, ncols_A)
    nnz_est_proc    = CEILING( REAL( nnz_est, dp) / REAL( par%n, dp))
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh_tri%TriA) / MINVAL( mesh_Vor%A   )), &
                                    CEILING( 2._dp * MAXVAL( mesh_Vor%A   ) / MINVAL( mesh_tri%TriA)) ))

    CALL allocate_matrix_CSR_dist( B_xdy_b_a_CSR  , nrows_A, ncols_A, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( B_mxydx_b_a_CSR, nrows_A, ncols_A, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( B_xydy_b_a_CSR , nrows_A, ncols_A, nnz_est_proc)

    ! Initialise results from integrating a single triangle through the Voronoi cells
    single_row%n_max = 100
    single_row%n     = 0
    ALLOCATE( single_row%index_left( single_row%n_max))
    ALLOCATE( single_row%LI_xdy(     single_row%n_max))
    ALLOCATE( single_row%LI_mxydx(   single_row%n_max))
    ALLOCATE( single_row%LI_xydy(    single_row%n_max))

  ! == Trace all the line segments to fill the matrices
  ! ===================================================

    vi_hint = 1

    DO ti = mesh_tri%ti1, mesh_tri%ti2

      !WRITE(0,*) '  Process ', par%i, ' - integrating mesh triangle ', ti, '/', mesh_tri%nTri, ' sides through the Voronoi cells of the opposite mesh...'

      ! Clean up single row results
      single_row%n            = 0
      single_row%index_left   = 0
      single_row%LI_xdy       = 0
      single_row%LI_mxydx     = 0
      single_row%LI_xydy      = 0

      ! The three vertices spanning this triangle
      via = mesh_tri%Tri( ti,1)
      vib = mesh_tri%Tri( ti,2)
      vic = mesh_tri%Tri( ti,3)

      ! Integrate over the three triangle sides
      p = mesh_tri%V( via,:)
      q = mesh_tri%V( vib,:)
      CALL trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      p = mesh_tri%V( vib,:)
      q = mesh_tri%V( vic,:)
      CALL trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      p = mesh_tri%V( vic,:)
      q = mesh_tri%V( via,:)
      CALL trace_line_Vor( mesh_Vor, p, q, single_row, count_coincidences, vi_hint)

      ! Add the results for this triangle to the sparse matrix
      DO k = 1, single_row%n
        CALL add_entry_CSR_dist( B_xdy_b_a_CSR  , ti, single_row%index_left( k), single_row%LI_xdy(   k))
        CALL add_entry_CSR_dist( B_mxydx_b_a_CSR, ti, single_row%index_left( k), single_row%LI_mxydx( k))
        CALL add_entry_CSR_dist( B_xydy_b_a_CSR , ti, single_row%index_left( k), single_row%LI_xydy(  k))
      END DO

    END DO ! DO ti = mesh_tri%ti1, mesh_tri%ti2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( B_xdy_b_a_CSR  , mesh_tri%ti1, mesh_tri%ti2)
    CALL finalise_matrix_CSR_dist( B_mxydx_b_a_CSR, mesh_tri%ti1, mesh_tri%ti2)
    CALL finalise_matrix_CSR_dist( B_xydy_b_a_CSR , mesh_tri%ti1, mesh_tri%ti2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( B_xdy_b_a_CSR  , B_xdy_b_a  )
    CALL mat_CSR2petsc( B_mxydx_b_a_CSR, B_mxydx_b_a)
    CALL mat_CSR2petsc( B_xydy_b_a_CSR , B_xydy_b_a )

    ! Clean up the Fortran versions
    CALL deallocate_matrix_CSR( B_xdy_b_a_CSR  )
    CALL deallocate_matrix_CSR( B_mxydx_b_a_CSR)
    CALL deallocate_matrix_CSR( B_xydy_b_a_CSR )

    ! Clean up after yourself
    DEALLOCATE( single_row%index_left )
    DEALLOCATE( single_row%LI_xdy     )
    DEALLOCATE( single_row%LI_mxydx   )
    DEALLOCATE( single_row%LI_xydy    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_triangles_through_Voronoi_cells
  SUBROUTINE integrate_Voronoi_cells_through_triangles( mesh_Vor, mesh_tri, B_xdy_a_b, B_mxydx_a_b, B_xydy_a_b, count_coincidences)
    ! Integrate around the grid cells of the grid through the triangles of the mesh

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_Vor
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_tri
    TYPE(tMat),                          INTENT(INOUT) :: B_xdy_a_b
    TYPE(tMat),                          INTENT(INOUT) :: B_mxydx_a_b
    TYPE(tMat),                          INTENT(INOUT) :: B_xydy_a_b
    LOGICAL,                             INTENT(IN)    :: count_coincidences

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'integrate_Voronoi_cells_through_triangles'
    INTEGER                                            :: nrows_A, ncols_A, nnz_est, nnz_est_proc, nnz_per_row_max
    TYPE(type_sparse_matrix_CSR_dp)                    :: B_xdy_a_b_CSR, B_mxydx_a_b_CSR, B_xydy_a_b_CSR
    TYPE(type_single_row_mapping_matrices)             :: single_row
    INTEGER                                            :: vi, nVor, vori1, vori2, k, ti_hint
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: Vor
    REAL(dp), DIMENSION(2)                             :: p, q


    ! Add routine to path
    CALL init_routine( routine_name)
  ! == Initialise the three matrices using the native UFEMISM CSR-matrix format
  ! ===========================================================================

    ! Matrix sise
    nrows_A         = mesh_Vor%nV    ! to
    ncols_A         = mesh_tri%nTri  ! from
    nnz_est         = 4 * MAX( nrows_A, ncols_A)
    nnz_est_proc    = CEILING( REAL( nnz_est, dp) / REAL( par%n, dp))
    nnz_per_row_max = MAX( 32, MAX( CEILING( 2._dp * MAXVAL( mesh_tri%TriA) / MINVAL( mesh_vor%A   )), &
                                    CEILING( 2._dp * MAXVAL( mesh_vor%A   ) / MINVAL( mesh_tri%TriA)) ))

    CALL allocate_matrix_CSR_dist( B_xdy_a_b_CSR  , nrows_A, ncols_A, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( B_mxydx_a_b_CSR, nrows_A, ncols_A, nnz_est_proc)
    CALL allocate_matrix_CSR_dist( B_xydy_a_b_CSR , nrows_A, ncols_A, nnz_est_proc)

    ! Initialise results from integrating a single triangle through the Voronoi cells
    single_row%n_max = 100
    single_row%n     = 0
    ALLOCATE( single_row%index_left( single_row%n_max))
    ALLOCATE( single_row%LI_xdy(     single_row%n_max))
    ALLOCATE( single_row%LI_mxydx(   single_row%n_max))
    ALLOCATE( single_row%LI_xydy(    single_row%n_max))

  ! == Trace all the line segments to fill the matrices
  ! ===================================================

    ALLOCATE( Vor( mesh_Vor%nC_mem+2,2))

    ti_hint = 1

    DO vi = mesh_Vor%vi1, mesh_Vor%vi2 ! +1 because PETSc indexes from 0

      !WRITE(0,*) '  Process ', par%i, ' - integrating around Voronoi cell ', vi, '/', mesh_Vor%nV, ' through the opposite mesh triangles...'

      ! Clean up single row results
      single_row%n            = 0
      single_row%index_left   = 0
      single_row%LI_xdy       = 0
      single_row%LI_mxydx     = 0
      single_row%LI_xydy      = 0

      ! Integrate over the complete Voronoi cell boundary
      CALL find_Voronoi_cell_vertices( mesh_Vor, vi, Vor, nVor)
      IF (mesh_Vor%edge_index( vi) > 0) THEN
        Vor( nVor+1,:) = Vor( 1,:)
        nVor = nVor + 1
      END IF
      DO vori1 = 1, nVor-1
        vori2 = vori1 + 1
        IF (vori2 > nVor) vori2 = 1
        p = Vor( vori1,:)
        q = Vor( vori2,:)
        CALL trace_line_tri( mesh_tri, p, q, single_row, count_coincidences, ti_hint)
      END DO

      ! Add the results for this triangle to the sparse matrix
      DO k = 1, single_row%n
        CALL add_entry_CSR_dist( B_xdy_a_b_CSR  , vi, single_row%index_left( k), single_row%LI_xdy(   k))
        CALL add_entry_CSR_dist( B_mxydx_a_b_CSR, vi, single_row%index_left( k), single_row%LI_mxydx( k))
        CALL add_entry_CSR_dist( B_xydy_a_b_CSR , vi, single_row%index_left( k), single_row%LI_xydy(  k))
      END DO

    END DO ! DO vi = mesh_Vor%vi1, mesh_Vor%vi2
    CALL sync

    ! Assemble matrices
    CALL finalise_matrix_CSR_dist( B_xdy_a_b_CSR  , mesh_Vor%vi1, mesh_Vor%vi2)
    CALL finalise_matrix_CSR_dist( B_mxydx_a_b_CSR, mesh_Vor%vi1, mesh_Vor%vi2)
    CALL finalise_matrix_CSR_dist( B_xydy_a_b_CSR , mesh_Vor%vi1, mesh_Vor%vi2)

    ! Convert matrices from Fortran to PETSc types
    CALL mat_CSR2petsc( B_xdy_a_b_CSR  , B_xdy_a_b  )
    CALL mat_CSR2petsc( B_mxydx_a_b_CSR, B_mxydx_a_b)
    CALL mat_CSR2petsc( B_xydy_a_b_CSR , B_xydy_a_b )

    ! Clean up the Fortran versions
    CALL deallocate_matrix_CSR( B_xdy_a_b_CSR  )
    CALL deallocate_matrix_CSR( B_mxydx_a_b_CSR)
    CALL deallocate_matrix_CSR( B_xydy_a_b_CSR )

    ! Clean up after yourself
    DEALLOCATE( Vor)
    DEALLOCATE( single_row%index_left )
    DEALLOCATE( single_row%LI_xdy     )
    DEALLOCATE( single_row%LI_mxydx   )
    DEALLOCATE( single_row%LI_xydy    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_Voronoi_cells_through_triangles

  ! Add the values for a single row of the three line-integral matrices
  SUBROUTINE add_integrals_to_single_row(  single_row, index_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
    ! Add the values for a single row of the three line-integral matrices

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_single_row_mapping_matrices), INTENT(INOUT) :: single_row
    INTEGER,                                INTENT(IN)    :: index_left
    REAL(dp),                               INTENT(IN)    :: LI_xdy, LI_mxydx, LI_xydy
    LOGICAL,                                INTENT(IN)    :: coincides, count_coincidences

    ! Local variables:
    LOGICAL                                            :: do_add_integrals, is_listed
    INTEGER                                            :: i, i_add

    ! Check whether we actually need to add the line integrals
    do_add_integrals = .TRUE.
    IF (coincides .AND. (.NOT. count_coincidences)) do_add_integrals = .FALSE.

    ! Check if an entry from this left-hand vertex is already listed
    is_listed = .FALSE.
    i_add     = 0

    DO i = 1, single_row%n
      IF (single_row%index_left( i) == index_left) THEN
        is_listed = .TRUE.
        i_add     = i
        EXIT
      END IF
    END DO
    IF (.NOT. is_listed) THEN
      single_row%n = single_row%n + 1
      i_add = single_row%n
    END IF

    ! Add data
    single_row%index_left( i_add) = index_left
    IF (do_add_integrals) THEN
      single_row%LI_xdy(   i_add) = single_row%LI_xdy(   i_add) + LI_xdy
      single_row%LI_mxydx( i_add) = single_row%LI_mxydx( i_add) + LI_mxydx
      single_row%LI_xydy(  i_add) = single_row%LI_xydy(  i_add) + LI_xydy
    END IF

    ! If necessary, extend memory
    IF (single_row%n > single_row%n_max - 10) CALL extend_single_row_memory( single_row, 100)

  END SUBROUTINE add_integrals_to_single_row
  SUBROUTINE extend_single_row_memory( single_row, n_extra)
    ! Extend memory for a single row of the three line-integral matrices

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_single_row_mapping_matrices), INTENT(INOUT) :: single_row
    INTEGER,                                INTENT(IN)    :: n_extra

    ! Local variables:
    INTEGER                                            :: n
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: index_left_temp
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: LI_xdy_temp, LI_mxydx_temp, LI_xydy_temp

    n = single_row%n

    ! Allocate temporary memory
    ALLOCATE( index_left_temp( n))
    ALLOCATE( LI_xdy_temp(     n))
    ALLOCATE( LI_mxydx_temp(   n))
    ALLOCATE( LI_xydy_temp(    n))

    ! Copy data to temporary memory
    index_left_temp = single_row%index_left( 1:n)
    LI_xdy_temp     = single_row%LI_xdy(     1:n)
    LI_mxydx_temp   = single_row%LI_mxydx(   1:n)
    LI_xydy_temp    = single_row%LI_xydy(    1:n)

    ! Deallocate memory
    DEALLOCATE( single_row%index_left)
    DEALLOCATE( single_row%LI_xdy    )
    DEALLOCATE( single_row%LI_mxydx  )
    DEALLOCATE( single_row%LI_xydy   )

    ! Allocate new, extended memory
    single_row%n_max = single_row%n_max + n_extra
    ALLOCATE( single_row%index_left( single_row%n_max))
    ALLOCATE( single_row%LI_xdy(     single_row%n_max))
    ALLOCATE( single_row%LI_mxydx(   single_row%n_max))
    ALLOCATE( single_row%LI_xydy(    single_row%n_max))

    ! Copy data back from temporary memory
    single_row%index_left( 1:n) = index_left_temp
    single_row%LI_xdy(     1:n) = LI_xdy_temp
    single_row%LI_mxydx(   1:n) = LI_mxydx_temp
    single_row%LI_xydy(    1:n) = LI_xydy_temp

    ! Deallocate temporary memory
    DEALLOCATE( index_left_temp)
    DEALLOCATE( LI_xdy_temp    )
    DEALLOCATE( LI_mxydx_temp  )
    DEALLOCATE( LI_xydy_temp   )

  END SUBROUTINE extend_single_row_memory

  ! Line tracing algorithm through mesh triangles
  SUBROUTINE trace_line_tri( mesh, p, q, single_row, count_coincidences, ti_hint)
    ! Trace the line [pq] through the triangles of the mesh and calculate
    ! the three line integrals for the line segments inside the different triangles

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    TYPE(type_single_row_mapping_matrices), INTENT(INOUT) :: single_row
    LOGICAL,                             INTENT(IN)    :: count_coincidences
    INTEGER,                             INTENT(INOUT) :: ti_hint

    ! Local variables:
    REAL(dp), DIMENSION(2)                             :: pp, qq, pa, pb, pc
    LOGICAL                                            :: is_valid_line
    INTEGER                                            :: edge_index_pq
    LOGICAL                                            :: finished
    INTEGER                                            :: n_cycles
    INTEGER                                            :: ti_in, vi_on, aci_on
    REAL(dp), DIMENSION(2)                             :: p_next
    INTEGER                                            :: ti_left
    LOGICAL                                            :: coincides
    REAL(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy
    INTEGER                                            :: ti_p, ti_q, vi_p, vi_q, ti_next

    ! Crop the line [pq] so that it lies within the mesh domain
    CALL crop_line_to_domain( p, q, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp, qq, is_valid_line)

    IF (.NOT. is_valid_line) THEN
      ! [pq] doesn't pass through the mesh domain anywhere
      RETURN
    END IF

    ! Check whether [pq] lies on the domain border
    edge_index_pq = 0
    IF     (ABS( p(1) - mesh%xmin) < mesh%tol_dist .AND. ABS( q(1) - mesh%xmin) < mesh%tol_dist) THEN
      ! pq lies on the western border
      edge_index_pq = 7
    ELSEIF (ABS( p(1) - mesh%xmax) < mesh%tol_dist .AND. ABS( q(1) - mesh%xmax) < mesh%tol_dist) THEN
      ! pq lies on the eastern border
      edge_index_pq = 3
    ELSEIF (ABS( p(2) - mesh%ymin) < mesh%tol_dist .AND. ABS( q(2) - mesh%ymin) < mesh%tol_dist) THEN
      ! pq lies on the southern border
      edge_index_pq = 5
    ELSEIF (ABS( p(2) - mesh%ymax) < mesh%tol_dist .AND. ABS( q(2) - mesh%ymax) < mesh%tol_dist) THEN
      ! pq lies on the northern border
      edge_index_pq = 1
    END IF

    IF (edge_index_pq == 0) THEN
      ! [pq] lies in the mesh interior

      ! Initialise the coincidence indicators for the point p, i.e. check IF p either...
      !    - lies inside the Voronoi cell of vertex vi_in, ...
      !    - lies on the circumcentre of triangle ti_on, or...
      !    - lies on the shared Voronoi cell boundary represented by edge aci_on
      CALL trace_line_tri_start( mesh, pp, ti_hint, ti_in, vi_on, aci_on)

      ! Iteratively trace the line through the mesh
      finished = .FALSE.
      n_cycles = 0
      DO WHILE (.NOT. finished)

        ! Find the point p_next where [pq] crosses into the next Voronoi cell
        IF     (ti_in  > 0) THEN
          ! p lies inside triangle ti_in
          CALL trace_line_tri_ti(  mesh, pp, qq, p_next, ti_in, vi_on, aci_on, ti_left, coincides, finished)
        ELSEIF (vi_on  > 0) THEN
          ! p lies on vertex vi_on
          CALL trace_line_tri_vi(  mesh, pp, qq, p_next, ti_in, vi_on, aci_on, ti_left, coincides, finished)
        ELSEIF (aci_on > 0) THEN
          ! p lies on edge aci_on
          CALL trace_line_tri_aci( mesh, pp, qq, p_next, ti_in, vi_on, aci_on, ti_left, coincides, finished)
        END IF

        ! Calculate the three line integrals
        CALL line_integral_xdy(   pp, p_next, mesh%tol_dist, LI_xdy  )
        CALL line_integral_mxydx( pp, p_next, mesh%tol_dist, LI_mxydx)
        CALL line_integral_xydy(  pp, p_next, mesh%tol_dist, LI_xydy )

        ! Add them to the results structure
        CALL add_integrals_to_single_row( single_row, ti_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

        ! Cycle the pointer
        pp = p_next

        ! Safety
        n_cycles = n_cycles + 1
        IF (n_cycles > mesh%nV) THEN
          CALL crash('trace_line_tri - iterative tracer got stuck!')
        END IF

        ! Update ti_hint, for more efficiency
        ti_hint = ti_left

      END DO ! DO WHILE (.NOT. finished)

    ELSE ! IF (edge_index_pq == 0) THEN
      ! [pq] lies on the mesh domain border

      ! Safety
      IF (edge_index_pq == 1) THEN
        ! North; q should be west of p
        IF (qq(1) >= pp(1)) THEN
          CALL crash('trace_line_tri - pq is not oriented counter-clockwise!')
        END IF
      ELSEIF (edge_index_pq == 3) THEN
        ! East; q should be north of p
        IF (qq(2) <= pp(2)) THEN
          CALL crash('trace_line_tri - pq is not oriented counter-clockwise!')
        END IF
      ELSEIF (edge_index_pq == 5) THEN
        ! South; q should be eat of p
        IF (qq(1) <= pp(1)) THEN
          CALL crash('trace_line_tri - pq is not oriented counter-clockwise!')
        END IF
      ELSEIF (edge_index_pq == 7) THEN
        ! West; q should be south of p
        IF (qq(2) >= pp(2)) THEN
          CALL crash('trace_line_tri - pq is not oriented counter-clockwise!')
        END IF
      END IF

      ! Find the triangles containing p and q
      vi_p = mesh%Tri( ti_hint,1)
      CALL find_containing_vertex( mesh, pp, vi_p)
      vi_q = vi_p
      CALL find_containing_vertex( mesh, qq, vi_q)

      ! Update ti_hint, for more efficiency
      ti_hint = mesh%iTri( vi_q,1)

      IF (edge_index_pq == 1) THEN
        ! Northern border
        IF (p(1) > mesh%V( vi_p,1)) THEN
          ! p lies east of vi_p; last iTriangle
          ti_p = mesh%iTri( vi_p, mesh%niTri( vi_p))
        ELSE
          ! p lies west of vi_p; first iTriangle
          ti_p = mesh%iTri( vi_p, 1)
        END IF
        IF (q(1) > mesh%V( vi_q,1)) THEN
          ! q lies east of vi_q; last iTriangle
          ti_q = mesh%iTri( vi_q, mesh%niTri( vi_q))
        ELSE
          ! q lies west of vi_q; first iTriangle
          ti_q = mesh%iTri( vi_q, 1)
        END IF
      ELSEIF (edge_index_pq == 3) THEN
        ! Eastern border
        IF (p(2) < mesh%V( vi_p,2)) THEN
          ! p lies south of vi_p; last iTriangle
          ti_p = mesh%iTri( vi_p, mesh%niTri( vi_p))
        ELSE
          ! p lies north of vi_p; first iTriangle
          ti_p = mesh%iTri( vi_p, 1)
        END IF
        IF (q(2) < mesh%V( vi_q,2)) THEN
          ! q lies south of vi_q; last iTriangle
          ti_q = mesh%iTri( vi_q, mesh%niTri( vi_q))
        ELSE
          ! q lies north of vi_q; first iTriangle
          ti_q = mesh%iTri( vi_q, 1)
        END IF
      ELSEIF (edge_index_pq == 5) THEN
        ! Southern border
        IF (p(1) < mesh%V( vi_p,1)) THEN
          ! p lies west of vi_p; last iTriangle
          ti_p = mesh%iTri( vi_p, mesh%niTri( vi_p))
        ELSE
          ! p lies east of vi_p; first iTriangle
          ti_p = mesh%iTri( vi_p, 1)
        END IF
        IF (q(1) < mesh%V( vi_q,1)) THEN
          ! q lies west of vi_q; last iTriangle
          ti_q = mesh%iTri( vi_q, mesh%niTri( vi_q))
        ELSE
          ! q lies east of vi_q; first iTriangle
          ti_q = mesh%iTri( vi_q, 1)
        END IF
      ELSEIF (edge_index_pq == 7) THEN
        ! Western border
        IF (p(2) > mesh%V( vi_p,2)) THEN
          ! p lies north of vi_p; last iTriangle
          ti_p = mesh%iTri( vi_p, mesh%niTri( vi_p))
        ELSE
          ! p lies south of vi_p; first iTriangle
          ti_p = mesh%iTri( vi_p, 1)
        END IF
        IF (q(2) > mesh%V( vi_q,2)) THEN
          ! q lies north of vi_q; last iTriangle
          ti_q = mesh%iTri( vi_q, mesh%niTri( vi_q))
        ELSE
          ! q lies south of vi_q; first iTriangle
          ti_q = mesh%iTri( vi_q, 1)
        END IF
      END IF

      ! Safety
      pa = mesh%V( mesh%Tri( ti_p,1),:)
      pb = mesh%V( mesh%Tri( ti_p,2),:)
      pc = mesh%V( mesh%Tri( ti_p,3),:)
      IF (.NOT. is_in_triangle( pa, pb, pc, pp)) THEN
        CALL crash('trace_line_tri - border version, p is not inside ti_p!')
      END IF

      pa = mesh%V( mesh%Tri( ti_q,1),:)
      pb = mesh%V( mesh%Tri( ti_q,2),:)
      pc = mesh%V( mesh%Tri( ti_q,3),:)
      IF (.NOT. is_in_triangle( pa, pb, pc, qq)) THEN
        CALL crash('trace_line_tri - border version, p is not inside ti_q!')
      END IF

      ! Iteratively trace the line through the mesh
      finished = .FALSE.
      n_cycles = 0
      DO WHILE (.NOT. finished)

        ! Find the point p_next where [pq] crosses into the next triangle
        CALL trace_line_tri_border( mesh, pp, qq, ti_p, ti_q, p_next, ti_next, coincides, finished)

        ! Calculate the three line integrals
        CALL line_integral_xdy(   pp, p_next, mesh%tol_dist, LI_xdy  )
        CALL line_integral_mxydx( pp, p_next, mesh%tol_dist, LI_mxydx)
        CALL line_integral_xydy(  pp, p_next, mesh%tol_dist, LI_xydy )

        ! Add them to the results structure
        CALL add_integrals_to_single_row( single_row, ti_p, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

        ! Cycle the pointer
        pp   = p_next
        ti_p = ti_next

        ! Safety
        n_cycles = n_cycles + 1
        IF (n_cycles > mesh%nV) THEN
          CALL crash('trace_line_tri - iterative tracer (border version) got stuck!')
        END IF

      END DO ! DO WHILE (.NOT. finished)

    END IF ! IF (edge_index_pq == 0) THEN

  END SUBROUTINE trace_line_tri
  SUBROUTINE trace_line_tri_start( mesh, p, ti_hint, ti_in, vi_on, aci_on)
    ! Initialise the coincidence indicators for the point p, i.e. check IF p either...
    !    - lies inside triangle ti_in, ...
    !    - lies on vertex vi_on, or...
    !    - lies on edge aci_on

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    INTEGER,                             INTENT(INOUT) :: ti_hint
    INTEGER,                             INTENT(OUT)   :: ti_in
    INTEGER,                             INTENT(OUT)   :: vi_on
    INTEGER,                             INTENT(OUT)   :: aci_on

    ! Local variables:
    INTEGER                                            :: via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    INTEGER                                            :: vvi, vj, aci

    ! Initialise
    ti_in  = 0
    vi_on  = 0
    aci_on = 0

    ! Find the triangle containing p
    CALL find_containing_triangle( mesh, p, ti_hint)

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_hint,1)
    vib = mesh%Tri( ti_hint,2)
    vic = mesh%Tri( ti_hint,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

    ! Check IF p lies on any of the three vertices
    IF     (NORM2( pa - p) < mesh%tol_dist) THEN
      ! p lies on via
      vi_on = via
      RETURN
    ELSEIF (NORM2( pb - p) < mesh%tol_dist) THEN
      ! p lies on vib
      vi_on = vib
      RETURN
    ELSEIF (NORM2( pc - p) < mesh%tol_dist) THEN
      ! p lies on vic
      vi_on = vic
      RETURN
    END IF

    ! Check IF p lies on any of the three edges
    IF     (lies_on_line_segment( pa, pb, p, mesh%tol_dist)) THEN
      ! p lies on the edge connecting via and vib
      DO vvi = 1, mesh%nC( via)
        vj  = mesh%C(    via,vvi)
        aci = mesh%iAci( via,vvi)
        IF (vj == vib) THEN
          aci_on = aci
          RETURN
        END IF
      END DO
    ELSEIF (lies_on_line_segment( pb, pc, p, mesh%tol_dist)) THEN
      ! p lies on the edge connecting vib and vic
      DO vvi = 1, mesh%nC( vib)
        vj  = mesh%C(    vib,vvi)
        aci = mesh%iAci( vib,vvi)
        IF (vj == vic) THEN
          aci_on = aci
          RETURN
        END IF
      END DO
    ELSEIF (lies_on_line_segment( pc, pa, p, mesh%tol_dist)) THEN
      ! p lies on the edge connecting vic and via
      DO vvi = 1, mesh%nC( vic)
        vj  = mesh%C(    vic,vvi)
        aci = mesh%iAci( vic,vvi)
        IF (vj == via) THEN
          aci_on = aci
          RETURN
        END IF
      END DO
    END IF

    ! IF p lies not on the vertices or edges of the triangle, then it must lie inside of it
    ti_in = ti_hint

  END SUBROUTINE trace_line_tri_start
  SUBROUTINE trace_line_tri_ti(  mesh, p, q, p_next, ti_in, vi_on, aci_on, ti_left, coincides, finished)
    ! Given the line [pq], where p lies inside triangle ti_in,
    ! find the point p_next where [pq] crosses into the next triangle.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: ti_in
    INTEGER,                             INTENT(INOUT) :: vi_on
    INTEGER,                             INTENT(INOUT) :: aci_on
    INTEGER,                             INTENT(OUT)   :: ti_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc
    INTEGER                                            :: vvi, vj, aci
    REAL(dp), DIMENSION(2)                             :: llis
    LOGICAL                                            :: do_cross

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_in,1)
    vib = mesh%Tri( ti_in,2)
    vic = mesh%Tri( ti_in,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

    ! Safety
    IF (ti_in == 0 .OR. vi_on > 0 .OR. aci_on > 0 .OR. (.NOT. is_in_triangle( pa, pb, pc, p))) THEN
      CALL crash('trace_line_tri_ti - coincidence indicators dont make sense!')
    END IF

    ! Check if q lies inside the same triangle
    IF (is_in_triangle( pa, pb, pc, q)) THEN
      ! q lies inside the same triangle
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      aci_on    = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on vertex via
    IF (NORM2( pa - q) < mesh%tol_dist) THEN
      ! q lies on vertex via
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      aci_on    = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on vertex vib
    IF (NORM2( pb - q) < mesh%tol_dist) THEN
      ! q lies on vertex vib
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      aci_on    = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on vertex vic
    IF (NORM2( pc - q) < mesh%tol_dist) THEN
      ! q lies on vertex vic
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      aci_on    = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on edge via-vib
    IF (lies_on_line_segment( pa, pb, q, mesh%tol_dist)) THEN
      ! q lies on edge via-vib
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      aci_on    = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on edge vib-vic
    IF (lies_on_line_segment( pb, pc, q, mesh%tol_dist)) THEN
      ! q lies on edge vib-vic
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      aci_on    = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on edge vic-via
    IF (lies_on_line_segment( pc, pa, q, mesh%tol_dist)) THEN
      ! q lies on edge vic-via
      p_next    = q
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      aci_on    = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if [pq] passes through via
    IF (lies_on_line_segment( p, q, pa, mesh%tol_dist)) THEN
      ! [pq] passes through via
      p_next    = pa
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = via
      aci_on    = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through vib
    IF (lies_on_line_segment( p, q, pb, mesh%tol_dist)) THEN
      ! [pq] passes through vib
      p_next    = pb
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = vib
      aci_on    = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through vic
    IF (lies_on_line_segment( p, q, pc, mesh%tol_dist)) THEN
      ! [pq] passes through vic
      p_next    = pc
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = vic
      aci_on    = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] crosses edge via-vib
    CALL segment_intersection( p, q, pa, pb, llis, do_cross, mesh%tol_dist)
    IF (do_cross) THEN
      ! [pq] crosses edge [via,vib]
      ! Find the edge connecting via and vib
      DO vvi = 1, mesh%nC( via)
        vj  = mesh%C(    via,vvi)
        aci = mesh%iAci( via,vvi)
        IF (vj == vib) THEN
          aci_on = aci
          EXIT
        END IF
      END DO
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] crosses edge vib-vic
    CALL segment_intersection( p, q, pb, pc, llis, do_cross, mesh%tol_dist)
    IF (do_cross) THEN
      ! [pq] crosses edge [vib,vic]
      ! Find the edge connecting vib and vic
      DO vvi = 1, mesh%nC( vib)
        vj  = mesh%C(    vib,vvi)
        aci = mesh%iAci( vib,vvi)
        IF (vj == vic) THEN
          aci_on = aci
          EXIT
        END IF
      END DO
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] crosses edge vic-via
    CALL segment_intersection( p, q, pc, pa, llis, do_cross, mesh%tol_dist)
    IF (do_cross) THEN
      ! [pq] crosses edge [vic,via]
      ! Find the edge connecting vic and via
      DO vvi = 1, mesh%nC( vic)
        vj  = mesh%C(    vic,vvi)
        aci = mesh%iAci( vic,vvi)
        IF (vj == via) THEN
          aci_on = aci
          EXIT
        END IF
      END DO
      p_next    = llis
      ti_left   = ti_in
      ti_in     = 0
      vi_on     = 0
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_tri_ti - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_tri_ti
  SUBROUTINE trace_line_tri_vi(  mesh, p, q, p_next, ti_in, vi_on, aci_on, ti_left, coincides, finished)
    ! Given the line [pq], where p lies on vertex vi_on,
    ! find the point p_next where [pq] crosses into the next triangle.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: ti_in
    INTEGER,                             INTENT(INOUT) :: vi_on
    INTEGER,                             INTENT(INOUT) :: aci_on
    INTEGER,                             INTENT(OUT)   :: ti_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vic
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc, pv
    INTEGER                                            :: vvi, vj, aci, n1, n2, n3, vti, ti
    REAL(dp), DIMENSION(2)                             :: llis
    LOGICAL                                            :: do_cross

    ! Safety
    IF (ti_in > 0 .OR. vi_on == 0 .OR. aci_on > 0 .OR. NORM2( p - mesh%V( vi_on,:)) > mesh%tol_dist) THEN
      CALL crash('trace_line_tri_vi - coincidence indicators dont make sense!')
    END IF

    ! Check IF q lies on any of the edges originating in this vertex
    DO vvi = 1, mesh%nC( vi_on)
      vj  = mesh%C(    vi_on,vvi)
      aci = mesh%iAci( vi_on,vvi)
      pv  = mesh%V( vj,:)
      IF (NORM2( mesh%V( vj,:) - q) < mesh%tol_dist .OR. &
          lies_on_line_segment( p, pv, q, mesh%tol_dist)) THEN
        ! q lies on edge aci, connecting vi_on and vj
        IF (mesh%Aci( aci,1) == vi_on) THEN
          ti_left = mesh%Aci( aci,5)
        ELSE
          ti_left = mesh%Aci( aci,6)
        END IF
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        aci_on    = 0
        coincides = .TRUE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check IF q lies inside any of the triangles surrounding vi_on
    DO vti = 1, mesh%niTri( vi_on)
      ti  = mesh%iTri( vi_on,vti)
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      pa  = mesh%V( via,:)
      pb  = mesh%V( vib,:)
      pc  = mesh%V( vic,:)
      IF (is_in_triangle( pa, pb, pc, q) .OR. &
          lies_on_line_segment( pa, pb, q, mesh%tol_dist) .OR. &
          lies_on_line_segment( pb, pc, q, mesh%tol_dist) .OR. &
          lies_on_line_segment( pc, pa, q, mesh%tol_dist)) THEN
        ! q lies inside adjacent triangle ti
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        aci_on    = 0
        ti_left   = ti
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check IF [pq] passes through any of the neighbouring vertices
    DO vvi = 1, mesh%nC( vi_on)
      vj  = mesh%C(    vi_on,vvi)
      aci = mesh%iAci( vi_on,vvi)
      pv  = mesh%V( vj,:)
      IF (lies_on_line_segment( p, q, pv, mesh%tol_dist)) THEN
        ! [pq] passes through neighbouring vertex vj, which is connected to vi_on by edge aci
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        aci_on    = 0
        IF (mesh%Aci( aci,1) == vi_on) THEN
          ti_left = mesh%Aci( aci,5)
        ELSE
          ti_left = mesh%Aci( aci,6)
        END IF
        coincides = .TRUE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check IF [pq] exits any of the adjacent triangles
    DO vti = 1, mesh%niTri( vi_on)
      ti  = mesh%iTri( vi_on,vti)
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        n3 = n2 + 1
        IF (n3 == 4) n3 = 1
        IF (mesh%Tri( ti,n1) == vi_on) THEN
          vib = mesh%Tri( ti,n2)
          vic = mesh%Tri( ti,n3)
          pb  = mesh%V( vib,:)
          pc  = mesh%V( vic,:)
          ! Find the opposite triangle edge
          aci = 0
          DO vvi = 1, mesh%nC( vib)
            vj = mesh%C( vib,vvi)
            IF (vj == vic) THEN
              aci = mesh%iAci( vib,vvi)
              EXIT
            END IF
          END DO
          CALL segment_intersection( p, q, pb, pc, llis, do_cross, mesh%tol_dist)
          IF (do_cross) THEN
            ! [pq] exits triangle ti through the opposite edge aci
            p_next    = llis
            ti_in     = 0
            vi_on     = 0
            aci_on    = aci
            ti_left   = ti
            coincides = .FALSE.
            finished  = .FALSE.
            RETURN
          END IF
        END IF
      END DO
    END DO

    ! This point should not be reachable!
    CALL crash('trace_line_tri_vi - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_tri_vi
  SUBROUTINE trace_line_tri_aci( mesh, p, q, p_next, ti_in, vi_on, aci_on, ti_left, coincides, finished)
    ! Given the line [pq], where p lies on edge aci,
    ! find the point p_next where [pq] crosses into the next triangle.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: ti_in
    INTEGER,                             INTENT(INOUT) :: vi_on
    INTEGER,                             INTENT(INOUT) :: aci_on
    INTEGER,                             INTENT(OUT)   :: ti_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vil, vir, til, tir
    REAL(dp), DIMENSION(2)                             :: pa, pb, pl, pr
    INTEGER                                            :: vvi, vj, aci
    REAL(dp), DIMENSION(2)                             :: llis
    LOGICAL                                            :: do_cross

    ! Some more info about this edge
    via = mesh%Aci( aci_on,1)
    vib = mesh%Aci( aci_on,2)
    vil = mesh%Aci( aci_on,3)
    vir = mesh%Aci( aci_on,4)
    til = mesh%Aci( aci_on,5)
    tir = mesh%Aci( aci_on,6)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    IF (vil > 0) pl  = mesh%V( vil,:)
    IF (vir > 0) pr  = mesh%V( vir,:)

    ! Safety
    IF (ti_in > 0 .OR. vi_on > 0 .OR. aci_on == 0 .OR. (.NOT. lies_on_line_segment( pa, pb, p, mesh%tol_dist))) THEN
      CALL crash('trace_line_tri_aci - coincidence indicators dont make sense!')
    END IF

    ! Check IF q lies on the same edge in the direction of via
    IF (lies_on_line_segment( p, pa, q, mesh%tol_dist)) THEN
      ! q lies on the same edge in the direction of via
      p_next    = q
      ti_in     = 0
      vi_on     = 0
      aci_on    = 0
      ti_left   = tir
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies on the same edge in the direction of vib
    IF (lies_on_line_segment( p, pb, q, mesh%tol_dist)) THEN
      ! q lies on the same edge in the direction of vib
      p_next    = q
      ti_in     = 0
      vi_on     = 0
      aci_on    = 0
      ti_left   = til
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies inside either of the two adjacent triangles
    IF (til > 0) THEN
      IF (is_in_triangle( pa, pb, pl, q)) THEN
        ! q lies inside triangle til
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        aci_on    = 0
        ti_left   = til
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END IF
    IF (tir > 0) THEN
      IF (is_in_triangle( pa, pr, pb, q)) THEN
        ! q lies inside triangle tir
        p_next    = q
        ti_in     = 0
        vi_on     = 0
        aci_on    = 0
        ti_left   = tir
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END IF

    ! Check IF [pq] passes through pa
    IF (lies_on_line_segment( p, q, pa, mesh%tol_dist)) THEN
      ! [pq] passes through pa
      p_next    = pa
      ti_in     = 0
      vi_on     = via
      aci_on    = 0
      ti_left   = tir
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] passes through pb
    IF (lies_on_line_segment( p, q, pb, mesh%tol_dist)) THEN
      ! [pq] passes through pb
      p_next    = pb
      ti_in     = 0
      vi_on     = vib
      aci_on    = 0
      ti_left   = til
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] passes through pl
    IF (til > 0) THEN
      IF (lies_on_line_segment( p, q, pl, mesh%tol_dist)) THEN
        ! [pq] passes through pl
        p_next    = pl
        ti_in     = 0
        vi_on     = vil
        aci_on    = 0
        ti_left   = til
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check IF [pq] passes through pr
    IF (tir > 0) THEN
      IF (lies_on_line_segment( p, q, pr, mesh%tol_dist)) THEN
        ! [pq] passes through pr
        p_next    = pr
        ti_in     = 0
        vi_on     = vir
        aci_on    = 0
        ti_left   = tir
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check IF [pq] crosses edge [via,vil]
    IF (til > 0) THEN
      CALL segment_intersection( p, q, pa, pl, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses edge [via,vil]
        ! Find the edge connecting via and vil
        DO vvi = 1, mesh%nC( via)
          vj  = mesh%C(    via,vvi)
          aci = mesh%iAci( via,vvi)
          IF (vj == vil) THEN
            aci_on = aci
            EXIT
          END IF
        END DO
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = til
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check IF [pq] crosses edge [vil,vib]
    IF (til > 0) THEN
      CALL segment_intersection( p, q, pl, pb, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses edge [vil,vib]
        ! Find the edge connecting vil and vib
        DO vvi = 1, mesh%nC( vil)
          vj  = mesh%C(    vil,vvi)
          aci = mesh%iAci( vil,vvi)
          IF (vj == vib) THEN
            aci_on = aci
            EXIT
          END IF
        END DO
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = til
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check IF [pq] crosses edge [via,vir]
    IF (tir > 0) THEN
      CALL segment_intersection( p, q, pa, pr, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses edge [via,vir]
        ! Find the edge connecting via and vir
        DO vvi = 1, mesh%nC( via)
          vj  = mesh%C(    via,vvi)
          aci = mesh%iAci( via,vvi)
          IF (vj == vir) THEN
            aci_on = aci
            EXIT
          END IF
        END DO
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = tir
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check IF [pq] crosses edge [vir,vib]
    IF (tir > 0) THEN
      CALL segment_intersection( p, q, pr, pb, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses edge [vir,vib]
        ! Find the edge connecting vir and vib
        DO vvi = 1, mesh%nC( vir)
          vj  = mesh%C(    vir,vvi)
          aci = mesh%iAci( vir,vvi)
          IF (vj == vib) THEN
            aci_on = aci
            EXIT
          END IF
        END DO
        p_next    = llis
        ti_in     = 0
        vi_on     = 0
        ti_left   = tir
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_tri_aci - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_tri_aci
  SUBROUTINE trace_line_tri_border( mesh, p, q, ti_p, ti_q, p_next, ti_next, coincides, finished)
    ! Given the line [pq] that lies on the mesh domain border, where p lies inside
    ! the triangle ti_p, find the point p_next where [pq] crosses into the next triangle.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    INTEGER,                             INTENT(IN)    :: ti_p, ti_q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(OUT)   :: ti_next
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vic, vi_exit
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc

    ! Check IF q lies inside the same triangle as p
    IF (ti_q == ti_p) THEN
      ! q lies inside the same triangle as p
      p_next    = q
      ti_next   = 0
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Find the triangle vertex lying on [pq]
    via = mesh%Tri( ti_p,1)
    vib = mesh%Tri( ti_p,2)
    vic = mesh%Tri( ti_p,3)

    pa  = mesh%V( via,:)
    pb  = mesh%V( vib,:)
    pc  = mesh%V( vic,:)

    vi_exit = 0
    IF     (lies_on_line_segment( p, q, pa, mesh%tol_dist) .AND. NORM2( p - pa) > mesh%tol_dist) THEN
      ! [pq] exits triangle ti_p through vertex via
      vi_exit = via
    ELSEIF (lies_on_line_segment( p, q, pb, mesh%tol_dist) .AND. NORM2( p - pb) > mesh%tol_dist) THEN
      ! [pq] exits triangle ti_p through vertex vib
      vi_exit = vib
    ELSEIF (lies_on_line_segment( p, q, pc, mesh%tol_dist) .AND. NORM2( p - pc) > mesh%tol_dist) THEN
      ! [pq] exits triangle ti_p through vertex vic
      vi_exit = vic
    END IF

    ! Safety
    IF (vi_exit == 0) THEN
      CALL crash('trace_line_tri_border - couldnt find vertex where [pq] exits triangle ti_p!')
    END IF

    ! Answer
    p_next    = mesh%V( vi_exit,:)
    ti_next   = mesh%iTri( vi_exit,1)
    coincides = .TRUE.
    finished  = .FALSE.

  END SUBROUTINE trace_line_tri_border

  ! Line tracing algorithm through mesh Voronoi cells
  SUBROUTINE trace_line_Vor( mesh, p, q, single_row, count_coincidences, vi_hint)
  ! Trace the line [pq] through the Voronoi cells of the mesh and calculate
  ! the three line integrals for the line segments inside the different Voronoi cells

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    TYPE(type_single_row_mapping_matrices), INTENT(INOUT) :: single_row
    LOGICAL,                             INTENT(IN)    :: count_coincidences
    INTEGER,                             INTENT(INOUT) :: vi_hint

    ! Local variables:
    REAL(dp), DIMENSION(2)                             :: pp,qq
    LOGICAL                                            :: is_valid_line
    INTEGER                                            :: edge_index_pq
    LOGICAL                                            :: finished
    INTEGER                                            :: n_cycles
    INTEGER                                            :: vi_in, ti_on, aci_on
    REAL(dp), DIMENSION(2)                             :: p_next
    INTEGER                                            :: vi_left
    LOGICAL                                            :: coincides
    REAL(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy
    INTEGER                                            :: vi_p, vi_q, vi_next

    ! Crop the line [pq] so that it lies within the mesh domain
    CALL crop_line_to_domain( p, q, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, pp, qq, is_valid_line)

    IF (.NOT.is_valid_line) THEN
      ! [pq] doesn't pass through the mesh domain anywhere
      RETURN
    END IF

    ! Check whether [pq] lies on the domain border
    edge_index_pq = 0
    IF     (ABS( p(1) - mesh%xmin) < mesh%tol_dist .AND. ABS( q(1) - mesh%xmin) < mesh%tol_dist) THEN
      ! pq lies on the western border
      edge_index_pq = 7
    ELSEIF (ABS( p(1) - mesh%xmax) < mesh%tol_dist .AND. ABS( q(1) - mesh%xmax) < mesh%tol_dist) THEN
      ! pq lies on the eastern border
      edge_index_pq = 3
    ELSEIF (ABS( p(2) - mesh%ymin) < mesh%tol_dist .AND. ABS( q(2) - mesh%ymin) < mesh%tol_dist) THEN
      ! pq lies on the southern border
      edge_index_pq = 5
    ELSEIF (ABS( p(2) - mesh%ymax) < mesh%tol_dist .AND. ABS( q(2) - mesh%ymax) < mesh%tol_dist) THEN
      ! pq lies on the northern border
      edge_index_pq = 1
    END IF

    IF (edge_index_pq == 0) THEN
      ! [pq] lies in the mesh interior

      ! Initialise the coincidence indicators for the point p, i.e. check IF p either...
      !    - lies inside the Voronoi cell of vertex vi_in, ...
      !    - lies on the circumcentre of triangle ti_on, or...
      !    - lies on the shared Voronoi cell boundary represented by edge aci_on
      CALL trace_line_Vor_start( mesh, pp, vi_hint, vi_in, ti_on, aci_on)

      ! Iteratively trace the line through the mesh
      finished = .FALSE.
      n_cycles = 0
      DO WHILE (.NOT.finished)

        ! Find the point p_next where [pq] crosses into the next Voronoi cell
        IF     (vi_in  > 0) THEN
          ! p lies inside the Voronoi cell of vertex vi_in
          CALL trace_line_Vor_vi(  mesh, pp, qq, p_next, vi_in, ti_on, aci_on, vi_left, coincides, finished)
        ELSEIF (ti_on  > 0) THEN
          ! p lies on the circumcentre of triangle ti_on
          CALL trace_line_Vor_ti(  mesh, pp, qq, p_next, vi_in, ti_on, aci_on, vi_left, coincides, finished)
        ELSEIF (aci_on > 0) THEN
          ! p lies on the shared Voronoi cell boundary represented by edge aci_on
          CALL trace_line_Vor_aci( mesh, pp, qq, p_next, vi_in, ti_on, aci_on, vi_left, coincides, finished)
        END IF

        ! Calculate the three line integrals
        CALL line_integral_xdy(   pp, p_next, mesh%tol_dist, LI_xdy  )
        CALL line_integral_mxydx( pp, p_next, mesh%tol_dist, LI_mxydx)
        CALL line_integral_xydy(  pp, p_next, mesh%tol_dist, LI_xydy )

        ! Add them to the results structure
        IF (NORM2( p_next - pp) > mesh%tol_dist) THEN
          CALL add_integrals_to_single_row( single_row, vi_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
        END IF

        ! Cycle the pointer
        pp = p_next

        ! Safety
        n_cycles = n_cycles + 1
        IF (n_cycles > mesh%nV) THEN
          CALL crash('trace_line_Vor - iterative tracer got stuck!')
        END IF

        ! Update vi_hint, for more efficiency
        vi_hint = vi_left

      END DO ! DO WHILE (.NOT.finished)

    ELSE ! IF (edge_index_pq == 0) THEN
      ! [pq] lies on the mesh domain border

      ! Safety
      IF (edge_index_pq == 1) THEN
        ! North; q should be west of p
        IF (qq(1) >= pp(1)) THEN
          CALL crash('trace_line_Vor - pq is not oriented counter-clockwise!')
        END IF
      ELSEIF (edge_index_pq == 3) THEN
        ! East; q should be north of p
        IF (qq(2) <= pp(2)) THEN
          CALL crash('trace_line_Vor - pq is not oriented counter-clockwise!')
        END IF
      ELSEIF (edge_index_pq == 5) THEN
        ! South; q should be eat of p
        IF (qq(1) <= pp(1)) THEN
          CALL crash('trace_line_Vor - pq is not oriented counter-clockwise!')
        END IF
      ELSEIF (edge_index_pq == 7) THEN
        ! West; q should be south of p
        IF (qq(2) >= pp(2)) THEN
          CALL crash('trace_line_Vor - pq is not oriented counter-clockwise!')
        END IF
      END IF

      ! Find the Voronoi cells containing p and q
      vi_p = vi_hint
      CALL find_containing_vertex( mesh, pp, vi_p)
      vi_q = vi_p
      CALL find_containing_vertex( mesh, qq, vi_q)

      ! Update vi_hint, for more efficiency
      vi_hint = vi_q

      ! Iteratively trace the line through the mesh
      finished = .FALSE.
      n_cycles = 0
      DO WHILE (.NOT.finished)

        ! Find the point p_next where [pq] crosses into the next Voronoi cell
        CALL trace_line_Vor_border( mesh, qq, vi_p, vi_q, p_next, vi_next, coincides, finished)

        ! Calculate the three line integrals
        CALL line_integral_xdy(   pp, p_next, mesh%tol_dist, LI_xdy  )
        CALL line_integral_mxydx( pp, p_next, mesh%tol_dist, LI_mxydx)
        CALL line_integral_xydy(  pp, p_next, mesh%tol_dist, LI_xydy )

        ! Add them to the results structure
        IF (NORM2( p_next - pp) > mesh%tol_dist) THEN
          CALL add_integrals_to_single_row( single_row, vi_p, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)
        END IF

        ! Cycle the pointer
        pp   = p_next
        vi_p = vi_next

        ! Safety
        n_cycles = n_cycles + 1
        IF (n_cycles > mesh%nV) THEN
          CALL crash('trace_line_Vor - iterative tracer (border version) got stuck!')
        END IF

      END DO ! DO WHILE (.NOT.finished)

    END IF ! IF (edge_index_pq == 0) THEN

  END SUBROUTINE trace_line_Vor
  SUBROUTINE trace_line_Vor_start( mesh, p, vi_hint, vi_in, ti_on, aci_on)
    ! Initialise the coincidence indicators for the point p, i.e. check if p either...
    !    - lies inside the Voronoi cell of vertex vi_in, ...
    !    - lies on the circumcentre of triangle ti_on, or...
    !    - lies on the shared Voronoi cell boundary represented by edge aci_on

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    INTEGER,                             INTENT(INOUT) :: vi_hint
    INTEGER,                             INTENT(OUT)   :: vi_in
    INTEGER,                             INTENT(OUT)   :: ti_on
    INTEGER,                             INTENT(OUT)   :: aci_on

    ! Local variables:
    INTEGER                                            :: vti, ti, vvi, aci
    REAL(dp), DIMENSION(2)                             :: cc1, cc2

    ! Initialise
    vi_in  = 0
    ti_on  = 0
    aci_on = 0

    ! Find the vertex whose Voronoi cell contains p
    CALL find_containing_vertex( mesh, p, vi_hint)

    ! Check IF p lies on any of the surrounding triangles' circumcentres
    DO vti = 1, mesh%niTri( vi_hint)
      ti = mesh%iTri( vi_hint,vti)
      IF (NORM2( mesh%Tricc( ti,:) - p) < mesh%tol_dist) THEN
        ! p lies on the circumcentre of triangle ti
        ti_on = ti
        RETURN
      END IF
    END DO

    ! Check IF p lies on any of the shared Voronoi boundaries
    DO vvi = 1, mesh%nC( vi_hint)
      aci = mesh%iAci( vi_hint,vvi)
      CALL find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
      IF (lies_on_line_segment( cc1, cc2, p, mesh%tol_dist)) THEN
        ! p lies on the shared Voronoi cell boundary represented by edge aci
        aci_on = aci
        RETURN
      END IF
    END DO

    ! IF p lies not on the boundary of the Voronoi cell, then it must lie inside of it
    vi_in = vi_hint

  END SUBROUTINE trace_line_Vor_start
  SUBROUTINE trace_line_Vor_vi(  mesh, p, q, p_next, vi_in, ti_on, aci_on, vi_left, coincides, finished)
    ! Given the line [pq], where p lies inside the Voronoi cell of vertex vi_in,
    ! find the point p_next where [pq] crosses into the next Voronoi cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: vi_in
    INTEGER,                             INTENT(INOUT) :: ti_on
    INTEGER,                             INTENT(INOUT) :: aci_on
    INTEGER,                             INTENT(OUT)   :: vi_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: vti, ti, vvi, aci
    REAL(dp), DIMENSION(2)                             :: cc1, cc2, r, llis
    LOGICAL                                            :: do_cross

    ! Safety
    IF (vi_in == 0 .OR. ti_on > 0 .OR. aci_on > 0 .OR. (.NOT. is_in_Voronoi_cell( mesh, p, vi_in))) THEN
      CALL crash('trace_line_Vor_vi - coincidence indicators dont make sense!')
    END IF

    ! Check IF q lies inside the same Voronoi cell
    IF (is_in_Voronoi_cell( mesh, q, vi_in)) THEN
      ! q lies inside the same Voronoi cell
      p_next    = q
      vi_left   = vi_in
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies on any of the surrounding triangles' circumcentres
    DO vti = 1, mesh%niTri( vi_in)
      ti = mesh%iTri( vi_in,vti)
      IF (NORM2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) THEN
        ! q lies on this triangle's circumcentre
        p_next    = q
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = 0
        aci_on    = 0
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check IF q lies on the Voronoi boundary
    DO vvi = 1, mesh%nC( vi_in)
      aci = mesh%iAci( vi_in,vvi)
      CALL find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
      IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) THEN
        ! q lies on this shared Voronoi boundary
        p_next    = q
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = 0
        aci_on    = 0
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check IF [pq] passes through any of the surrounding triangles' circumcentres
    DO vti = 1, mesh%niTri( vi_in)
      ti = mesh%iTri( vi_in,vti)
      r  = mesh%Tricc( ti,:)
      IF (lies_on_line_segment( p, q, r, mesh%tol_dist)) THEN
        ! [pq] passes through this triangle's circumcentre
        p_next    = mesh%Tricc( ti,:)
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = ti
        aci_on    = 0
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check IF [pq] passes through any of the shared Voronoi boundaries
    DO vvi = 1, mesh%nC( vi_in)
      aci = mesh%iAci( vi_in,vvi)
      CALL find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] passes through this shared Voronoi boundary
        p_next    = llis
        vi_left   = vi_in
        vi_in     = 0
        ti_on     = 0
        aci_on    = aci
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! This point should not be reachable!
    CALL crash('trace_line_Vor_vi - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_Vor_vi
  SUBROUTINE trace_line_Vor_ti(  mesh, p, q, p_next, vi_in, ti_on, aci_on, vi_left, coincides, finished)
    ! Given the line [pq], where p lies on the circumcentre of triangle ti_on,
    ! find the point p_next where [pq] crosses into the next Voronoi cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: vi_in
    INTEGER,                             INTENT(INOUT) :: ti_on
    INTEGER,                             INTENT(INOUT) :: aci_on
    INTEGER,                             INTENT(OUT)   :: vi_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vic, vvi, vj, aci, acab, acbc, acca, tj
    REAL(dp), DIMENSION(2)                             :: cc, cc1, cc2, llis
    LOGICAL                                            :: do_cross

    ! Safety
    IF (vi_in > 0 .OR. ti_on == 0 .OR. aci_on > 0 .OR. NORM2( mesh%Tricc( ti_on,:) - p) > mesh%tol_dist) THEN
      CALL crash('trace_line_Vor_ti - coincidence indicators dont make sense!')
    END IF

    ! The three vertices spanning the triangle
    via = mesh%Tri( ti_on,1)
    vib = mesh%Tri( ti_on,2)
    vic = mesh%Tri( ti_on,3)

    ! Find the three Voronoi cell boundaries that meet here
    acab = 0
    DO vvi = 1, mesh%nC( via)
      vj  = mesh%C(    via,vvi)
      aci = mesh%iAci( via,vvi)
      IF (vj == vib) THEN
        acab = aci
        EXIT
      END IF
    END DO
    acbc = 0
    DO vvi = 1, mesh%nC( vib)
      vj  = mesh%C(    vib,vvi)
      aci = mesh%iAci( vib,vvi)
      IF (vj == vic) THEN
        acbc = aci
        EXIT
      END IF
    END DO
    acca = 0
    DO vvi = 1, mesh%nC( vic)
      vj  = mesh%C(    vic,vvi)
      aci = mesh%iAci( vic,vvi)
      IF (vj == via) THEN
        acca = aci
        EXIT
      END IF
    END DO

    ! Check IF q lies on the Voronoi cell boundary separating via from vib
    CALL find_shared_Voronoi_boundary( mesh, acab, cc1, cc2)
    IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .OR. &
        NORM2( cc1 - q) < mesh%tol_dist .OR. &
        NORM2( cc2 - q) < mesh%tol_dist) THEN
      ! q lies on the Voronoi cell boundary separating via from vib
      IF      (mesh%Aci( acab,5) == ti_on) THEN
        vi_left = mesh%Aci( acab,2)
      ELSE
        vi_left = mesh%Aci( acab,1)
      END IF
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies on the Voronoi cell boundary separating vib from vic
    CALL find_shared_Voronoi_boundary( mesh, acbc, cc1, cc2)
    IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .OR. &
        NORM2( cc1 - q) < mesh%tol_dist .OR. &
        NORM2( cc2 - q) < mesh%tol_dist) THEN
      ! q lies on the Voronoi cell boundary separating vib from vic
      IF (mesh%Aci( acbc,5) == ti_on) THEN
        vi_left = mesh%Aci( acbc,2)
      ELSE
        vi_left = mesh%Aci( acbc,1)
      END IF
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies on the Voronoi cell boundary separating vic from via
    CALL find_shared_Voronoi_boundary( mesh, acca, cc1, cc2)
    IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist) .OR. &
        NORM2( cc1 - q) < mesh%tol_dist .OR. &
        NORM2( cc2 - q) < mesh%tol_dist) THEN
      ! q lies on the Voronoi cell boundary separating vic from via
      IF (mesh%Aci( acca,5) == ti_on) THEN
        vi_left = mesh%Aci( acca,2)
      ELSE
        vi_left = mesh%Aci( acca,1)
      END IF
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies inside any of the three adjacent Voronoi cells
    IF (is_in_Voronoi_cell( mesh, q, via)) THEN
      ! q lies inside the Voronoi cell of via
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      vi_left   = via
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF
    IF (is_in_Voronoi_cell( mesh, q, vib)) THEN
      ! q lies inside the Voronoi cell of vib
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      vi_left   = vib
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF
    IF (is_in_Voronoi_cell( mesh, q, vic)) THEN
      ! q lies inside the Voronoi cell of vic
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      vi_left   = vic
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF [pq] passes through the circumcentre of any of the three neighbouring triangles
    tj = mesh%TriC( ti_on,1)
    IF (tj > 0) THEN
      cc = mesh%Tricc( tj,:)
      IF (lies_on_line_segment( p, q, cc, mesh%tol_dist)) THEN
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        aci_on    = 0
        vi_left   = vic
        coincides = .TRUE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    tj = mesh%TriC( ti_on,2)
    IF (tj > 0) THEN
      cc = mesh%Tricc( tj,:)
      IF (lies_on_line_segment( p, q, cc, mesh%tol_dist)) THEN
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        aci_on    = 0
        vi_left   = via
        coincides = .TRUE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    tj = mesh%TriC( ti_on,3)
    IF (tj > 0) THEN
      cc = mesh%Tricc( tj,:)
      IF (lies_on_line_segment( p, q, cc, mesh%tol_dist)) THEN
        ! [pq] passes through the circumcentre of this neighbouring triangle
        p_next    = cc
        vi_in     = 0
        ti_on     = tj
        aci_on    = 0
        vi_left   = vib
        coincides = .TRUE.
        finished  = .FALSE.
        RETURN
      END IF
    END IF

    ! Check IF [pq] crosses the boundary of the Voronoi cell of via
    DO vvi = 1, mesh%nC( via)
      vj = mesh%C( via,vvi)
      IF (vj == vib .OR. vj == vic) CYCLE
      aci = mesh%iAci( via,vvi)
      CALL find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        aci_on    = aci
        vi_left   = via
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check IF [pq] crosses the boundary of the Voronoi cell of vib
    DO vvi = 1, mesh%nC( vib)
      vj = mesh%C( vib,vvi)
      IF (vj == via .OR. vj == vic) CYCLE
      aci = mesh%iAci( vib,vvi)
      CALL find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        aci_on    = aci
        vi_left   = vib
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check IF [pq] crosses the boundary of the Voronoi cell of vic
    DO vvi = 1, mesh%nC( vic)
      vj = mesh%C( vic,vvi)
      IF (vj == via .OR. vj == vib) CYCLE
      aci = mesh%iAci( vic,vvi)
      CALL find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses this part of the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        aci_on    = aci
        vi_left   = vic
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! This point should not be reachable!
    CALL crash('trace_line_Vor_ti - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_Vor_ti
  SUBROUTINE trace_line_Vor_aci( mesh, p, q, p_next, vi_in, ti_on, aci_on, vi_left, coincides, finished)
    ! Given the line [pq], where p lies on the shared Voronoi boundary represented by edge aci_on,
    ! find the point p_next where [pq] crosses into the next Voronoi cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(INOUT) :: vi_in
    INTEGER,                             INTENT(INOUT) :: ti_on
    INTEGER,                             INTENT(INOUT) :: aci_on
    INTEGER,                             INTENT(OUT)   :: vi_left
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Local variables:
    INTEGER                                            :: via, vib, vil, vir, til, tir, vvi, aci, vti, ti
    REAL(dp), DIMENSION(2)                             :: cc1, cc2, ccl, ccr, llis
    LOGICAL                                            :: do_cross

    ! Find the endpoints of this shared Voronoi boundary
    CALL find_shared_Voronoi_boundary( mesh, aci_on, cc1, cc2)

    ! Safety
    IF (vi_in > 0 .OR. ti_on > 0 .OR. aci_on == 0 .OR. (.NOT. lies_on_line_segment( cc1, cc2, p, mesh%tol_dist))) THEN
      CALL crash('trace_line_Vor_aci - coincidence indicators dont make sense!')
    END IF

    ! A bit more detail is needed
    via = mesh%Aci( aci_on,1)
    vib = mesh%Aci( aci_on,2)
    vil = mesh%Aci( aci_on,3)
    vir = mesh%Aci( aci_on,4)
    til = mesh%Aci( aci_on,5)
    tir = mesh%Aci( aci_on,6)

    IF (til == 0) THEN
      ! Apparently aci lies on the domain border and has no triangle on its left-hand side
      ccr = cc1
      ccl = cc2
    ELSEIF (tir == 0) THEN
      ! Apparently aci lies on the domain border and has no triangle on its right-hand side
      ccl = cc1
      ccr = cc2
    ELSE
      ! aci lies in the interior and has triangles on both sides
      ccl = mesh%Tricc( til,:)
      ccr = mesh%Tricc( tir,:)
    END IF

    ! Check IF q coincides with ccl
    IF (NORM2( ccl - q) < mesh%tol_dist) THEN
      ! q coincides with ccl
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      vi_left   = via
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q coincides with ccr
    IF (NORM2( ccr - q) < mesh%tol_dist) THEN
      ! q coincides with ccr
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      vi_left   = vib
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies inside the Voronoi cell of via
    IF (is_in_Voronoi_cell( mesh, q, via)) THEN
      ! q lies inside the Voronoi cell of via
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      vi_left   = via
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies inside the Voronoi cell of vib
    IF (is_in_Voronoi_cell( mesh, q, vib)) THEN
      ! q lies inside the Voronoi cell of vib
      p_next    = q
      vi_in     = 0
      ti_on     = 0
      aci_on    = 0
      vi_left   = vib
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies on the circumcentre of any of the triangles surrounding via
    DO vti = 1, mesh%niTri( via)
      ti = mesh%iTri( via,vti)
      IF (NORM2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) THEN
        ! q lies on this triangle's circumcentre
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        aci_on    = 0
        vi_left   = via
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check IF q lies on the circumcentre of any of the triangles surrounding vib
    DO vti = 1, mesh%niTri( vib)
      ti = mesh%iTri( vib,vti)
      IF (NORM2( mesh%Tricc( ti,:) - q) < mesh%tol_dist) THEN
        ! q lies on this triangle's circumcentre
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        aci_on    = 0
        vi_left   = vib
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check IF q lies on boundary of the Voronoi cell of via
    DO vvi = 1, mesh%nC( via)
      aci = mesh%iAci( via,vvi)
      CALL find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
      IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) THEN
        ! q lies on this shared Voronoi boundary
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        aci_on    = 0
        vi_left   = via
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check IF q lies on boundary of the Voronoi cell of vib
    DO vvi = 1, mesh%nC( vib)
      aci = mesh%iAci( vib,vvi)
      CALL find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
      IF (lies_on_line_segment( cc1, cc2, q, mesh%tol_dist)) THEN
        ! q lies on this shared Voronoi boundary
        p_next    = q
        vi_in     = 0
        ti_on     = 0
        aci_on    = 0
        vi_left   = vib
        coincides = .FALSE.
        finished  = .TRUE.
        RETURN
      END IF
    END DO

    ! Check IF [pq] passes through either of this boundary's endpoints
    IF (lies_on_line_segment( p, q, ccl, mesh%tol_dist)) THEN
      ! [pq] passes through ccl
      p_next    = ccl
      vi_in     = 0
      ti_on     = til
      aci_on    = 0
      vi_left   = via
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF
    IF (lies_on_line_segment( p, q, ccr, mesh%tol_dist)) THEN
      ! [pq] passes through ccr
      p_next    = ccr
      vi_in     = 0
      ti_on     = tir
      aci_on    = 0
      vi_left   = vib
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF pq crosses the boundary of the Voronoi cell of via
    DO vvi = 1, mesh%nC( via)
      aci = mesh%iAci( via,vvi)
      IF (aci == aci_on) CYCLE
      CALL find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] passes through the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        aci_on    = aci
        vi_left   = via
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! Check IF pq crosses the boundary of the Voronoi cell of vib
    DO vvi = 1, mesh%nC( vib)
      aci = mesh%iAci( vib,vvi)
      IF (aci == aci_on) CYCLE
      CALL find_shared_Voronoi_boundary( mesh, aci, cc1, cc2)
      CALL segment_intersection( p, q, cc1, cc2, llis, do_cross, mesh%tol_dist)
      IF (do_cross) THEN
        ! [pq] passes through the boundary of the Voronoi cell of via
        p_next    = llis
        vi_in     = 0
        ti_on     = 0
        aci_on    = aci
        vi_left   = vib
        coincides = .FALSE.
        finished  = .FALSE.
        RETURN
      END IF
    END DO

    ! This point should not be reachable!
    CALL crash('trace_line_Vor_aci - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_Vor_aci
  SUBROUTINE trace_line_Vor_border( mesh, q, vi_p, vi_q, p_next, vi_next, coincides, finished)
    ! Given the line [pq] that lies on the mesh domain border, where p lies inside
    ! the Voronoi cell of vertex vi_p, find the point p_next where [pq] crosses into the next Voronoi cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: q
    INTEGER,                             INTENT(IN)    :: vi_p, vi_q
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(OUT)   :: vi_next
    LOGICAL,                             INTENT(OUT)   :: coincides
    LOGICAL,                             INTENT(OUT)   :: finished

    ! Check if q lies inside the same Voronoi cell as p
    IF (vi_q == vi_p) THEN
      ! q lies inside the same Voronoi cell as p
      p_next    = q
      vi_next   = 0
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Find the next Voronoi cell
    vi_next   = mesh%C( vi_p,1)
    p_next    = (mesh%V( vi_p,:) + mesh%V( vi_next,:)) / 2._dp
    coincides = .TRUE.
    finished  = .FALSE.

  END SUBROUTINE trace_line_Vor_border

  ! Line tracing algorithm through square grid cells
  SUBROUTINE trace_line_grid( grid, p, q, single_row, count_coincidences)
    ! Trace the line [pq] through the grid and calculate the three
    ! line integrals for the line segments inside the different grid cells.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    TYPE(type_single_row_mapping_matrices), INTENT(INOUT) :: single_row
    LOGICAL,                             INTENT(IN)    :: count_coincidences

    ! Local variables:
    REAL(dp)                                           :: xmin, xmax, ymin, ymax
    REAL(dp), DIMENSION(2)                             :: pp,qq
    LOGICAL                                            :: is_valid_line
    INTEGER                                            :: edge_index_pq
    LOGICAL                                            :: finished
    INTEGER                                            :: n_cycles
    INTEGER,  DIMENSION(2)                             :: aij_in, bij_on, cxij_on, cyij_on
    REAL(dp), DIMENSION(2)                             :: p_next
    INTEGER                                            :: n_left
    LOGICAL                                            :: coincides
    REAL(dp)                                           :: LI_xdy, LI_mxydx, LI_xydy
    INTEGER                                            :: i_p, j_p, i_q, j_q

    ! Crop the line [pq] so that it lies within the domain
    xmin = grid%xmin - grid%dx / 2
    xmax = grid%xmax + grid%dx / 2
    ymin = grid%ymin - grid%dx / 2
    ymax = grid%ymax + grid%dx / 2
    CALL crop_line_to_domain( p, q, xmin, xmax, ymin, ymax, grid%tol_dist, pp, qq, is_valid_line)

    IF (.NOT. is_valid_line) THEN
      ! [pq] doesn't pass through the domain anywhere
      RETURN
    END IF

    ! Check whether [pq] lies on the domain border
    edge_index_pq = 0
    IF     (ABS( p(1) - xmin) < grid%tol_dist .AND. ABS( q(1) - xmin) < grid%tol_dist) THEN
      ! pq lies on the western border
      edge_index_pq = 7
    ELSEIF (ABS( p(1) - xmax) < grid%tol_dist .AND. ABS( q(1) - xmax) < grid%tol_dist) THEN
      ! pq lies on the eastern border
      edge_index_pq = 3
    ELSEIF (ABS( p(2) - ymin) < grid%tol_dist .AND. ABS( q(2) - ymin) < grid%tol_dist) THEN
      ! pq lies on the southern border
      edge_index_pq = 5
    ELSEIF (ABS( p(2) - ymax) < grid%tol_dist .AND. ABS( q(2) - ymax) < grid%tol_dist) THEN
      ! pq lies on the northern border
      edge_index_pq = 1
    END IF

    IF (edge_index_pq == 0) THEN
      ! [pq] lies in the domain interior

      ! Initialise the coincidence indicators for the point p, i.e. check IF p either...
      !    - lies inside grid cell aij_in, ...
      !    - lies on the b-grid point bij_on, or...
      !    - lies on the edge cij_on
      CALL trace_line_grid_start( grid, pp, aij_in, bij_on, cxij_on, cyij_on)

      ! Iteratively trace the line through the mesh
      finished = .FALSE.
      n_cycles = 0
      DO WHILE (.NOT. finished)

        ! Find the point p_next where [pq] crosses into the next Voronoi cell
        IF     (aij_in(  1) > 0) THEN
          ! p lies inside a-grid cell aij_in
          CALL trace_line_grid_a(  grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        ELSEIF (bij_on(  1) > 0) THEN
          ! p lies on b-grid point bij_on
          CALL trace_line_grid_b(  grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        ELSEIF (cxij_on( 1) > 0) THEN
          ! p lies on cx-grid edge cxij_on
          CALL trace_line_grid_cx( grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        ELSEIF (cyij_on( 1) > 0) THEN
          ! p lies on cy-grid edge cyij_on
          CALL trace_line_grid_cy( grid, pp, qq, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
        END IF

        ! Calculate the three line integrals
        CALL line_integral_xdy(   pp, p_next, grid%tol_dist, LI_xdy  )
        CALL line_integral_mxydx( pp, p_next, grid%tol_dist, LI_mxydx)
        CALL line_integral_xydy(  pp, p_next, grid%tol_dist, LI_xydy )

        ! Add them to the results structure
        CALL add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

        ! Cycle the pointer
        pp = p_next

        ! Safety
        n_cycles = n_cycles + 1
        IF (n_cycles > grid%n) THEN
          CALL crash('trace_line_grid - iterative tracer got stuck!')
        END IF

      END DO ! DO WHILE (.NOT. finished)

    ELSE ! IF (edge_index_pq == 0) THEN
      ! [pq] lies on the mesh domain border

      ! Safety
      IF (edge_index_pq == 1) THEN
        ! North; q should be west of p
        IF (qq(1) >= pp(1)) THEN
          CALL crash('trace_line_grid - pq is not oriented counter-clockwise!')
        END IF
      ELSEIF (edge_index_pq == 3) THEN
        ! East; q should be north of p
        IF (qq(2) <= pp(2)) THEN
          CALL crash('trace_line_grid - pq is not oriented counter-clockwise!')
        END IF
      ELSEIF (edge_index_pq == 5) THEN
        ! South; q should be eat of p
        IF (qq(1) <= pp(1)) THEN
          CALL crash('trace_line_grid - pq is not oriented counter-clockwise!')
        END IF
      ELSEIF (edge_index_pq == 7) THEN
        ! West; q should be south of p
        IF (qq(2) >= pp(2)) THEN
          CALL crash('trace_line_grid - pq is not oriented counter-clockwise!')
        END IF
      END IF

      ! Find the grid cells containing p and q
      IF (edge_index_pq == 1) THEN
        ! [pq] lies on the northern boundary

        ! Find the grid cells containing p and q
        i_p = 1 + FLOOR( (pp(1) - xmin + grid%dx / 2._dp) / grid%dx)
        i_q = 1 + FLOOR( (qq(1) - xmin + grid%dx / 2._dp) / grid%dx)
        j_p = grid%ny
        j_q = grid%ny

        ! Iteratively trace the line through the mesh
        n_cycles = 0
        DO WHILE (i_p > i_q)

          ! Find the point where [pq] crosses into the next grid cell (very easy on a square grid, hurray!)
          p_next    = [grid%x( i_p) - grid%dx / 2._dp, ymax + grid%dx / 2._dp]
          n_left    = grid%ij2n( i_p, j_p)
          coincides = .TRUE.

          ! Calculate the three line integrals
          CALL line_integral_xdy(   pp, p_next, grid%tol_dist, LI_xdy  )
          CALL line_integral_mxydx( pp, p_next, grid%tol_dist, LI_mxydx)
          CALL line_integral_xydy(  pp, p_next, grid%tol_dist, LI_xydy )

          ! Add them to the results structure
          CALL add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

          ! Cycle the pointer
          pp  = p_next
          i_p = i_p - 1

          ! Safety
          n_cycles = n_cycles + 1
          IF (n_cycles > grid%nx) THEN
            CALL crash('trace_line_grid - iterative tracer (north border version) got stuck!')
          END IF

        END DO ! DO WHILE (i_p > i_q)

        ! Add the last section (inside the grid cell containing q)
        CALL line_integral_xdy(   pp, qq, grid%tol_dist, LI_xdy  )
        CALL line_integral_mxydx( pp, qq, grid%tol_dist, LI_mxydx)
        CALL line_integral_xydy(  pp, qq, grid%tol_dist, LI_xydy )

        ! Add them to the results structure
        n_left = grid%ij2n( i_q, j_q)
        CALL add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

      ELSEIF (edge_index_pq == 3) THEN
        ! [pq] lies on the eastern boundary

        ! Find the grid cells containing p and q
        j_p = 1 + FLOOR( (pp(2) - ymin + grid%dx / 2._dp) / grid%dx)
        j_q = 1 + FLOOR( (qq(2) - ymin + grid%dx / 2._dp) / grid%dx)
        i_p = grid%nx
        i_q = grid%nx

        ! Iteratively trace the line through the mesh
        n_cycles = 0
        DO WHILE (j_p < j_q)

          ! Find the point where [pq] crosses into the next grid cell (very easy on a square grid, hurray!)
          p_next    = [xmax + grid%dx / 2._dp, grid%y( j_p) + grid%dx / 2._dp]
          n_left    = grid%ij2n( i_p, j_p)
          coincides = .TRUE.

          ! Calculate the three line integrals
          CALL line_integral_xdy(   pp, p_next, grid%tol_dist, LI_xdy  )
          CALL line_integral_mxydx( pp, p_next, grid%tol_dist, LI_mxydx)
          CALL line_integral_xydy(  pp, p_next, grid%tol_dist, LI_xydy )

          ! Add them to the results structure
          CALL add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

          ! Cycle the pointer
          pp  = p_next
          j_p = j_p + 1

          ! Safety
          n_cycles = n_cycles + 1
          IF (n_cycles > grid%nx) THEN
            CALL crash('trace_line_grid - iterative tracer (east border version) got stuck!')
          END IF

        END DO ! DO WHILE (j_p < j_q)

        ! Add the last section (inside the grid cell containing q)
        CALL line_integral_xdy(   pp, qq, grid%tol_dist, LI_xdy  )
        CALL line_integral_mxydx( pp, qq, grid%tol_dist, LI_mxydx)
        CALL line_integral_xydy(  pp, qq, grid%tol_dist, LI_xydy )

        ! Add them to the results structure
        n_left = grid%ij2n( i_q, j_q)
        CALL add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

      ELSEIF (edge_index_pq == 5) THEN
        ! [pq] lies on the southern boundary

        ! Find the grid cells containing p and q
        i_p = 1 + FLOOR( (pp(1) - xmin + grid%dx / 2._dp) / grid%dx)
        i_q = 1 + FLOOR( (qq(1) - xmin + grid%dx / 2._dp) / grid%dx)
        j_p = 1
        j_q = 1

        ! Iteratively trace the line through the mesh
        n_cycles = 0
        DO WHILE (i_p < i_q)

          ! Find the point where [pq] crosses into the next grid cell (very easy on a square grid, hurray!)
          p_next    = [grid%x( i_p) + grid%dx / 2._dp, ymin - grid%dx / 2._dp]
          n_left    = grid%ij2n( i_p, j_p)
          coincides = .TRUE.

          ! Calculate the three line integrals
          CALL line_integral_xdy(   pp, p_next, grid%tol_dist, LI_xdy  )
          CALL line_integral_mxydx( pp, p_next, grid%tol_dist, LI_mxydx)
          CALL line_integral_xydy(  pp, p_next, grid%tol_dist, LI_xydy )

          ! Add them to the results structure
          CALL add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

          ! Cycle the pointer
          pp  = p_next
          i_p = i_p + 1

          ! Safety
          n_cycles = n_cycles + 1
          IF (n_cycles > grid%nx) THEN
            CALL crash('trace_line_grid - iterative tracer (south border version) got stuck!')
          END IF

        END DO ! DO WHILE (i_p < i_q)

        ! Add the last section (inside the grid cell containing q)
        CALL line_integral_xdy(   pp, qq, grid%tol_dist, LI_xdy  )
        CALL line_integral_mxydx( pp, qq, grid%tol_dist, LI_mxydx)
        CALL line_integral_xydy(  pp, qq, grid%tol_dist, LI_xydy )

        ! Add them to the results structure
        n_left = grid%ij2n( i_q, j_q)
        CALL add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

      ELSEIF (edge_index_pq == 7) THEN
        ! [pq] lies on the western boundary

        ! Find the grid cells containing p and q
        j_p = 1 + FLOOR( (pp(2) - ymin + grid%dx / 2._dp) / grid%dx)
        j_q = 1 + FLOOR( (qq(2) - ymin + grid%dx / 2._dp) / grid%dx)
        i_p = 1
        i_q = 1

        ! Iteratively trace the line through the mesh
        n_cycles = 0
        DO WHILE (j_p > j_q)

          ! Find the point where [pq] crosses into the next grid cell (very easy on a square grid, hurray!)
          p_next    = [xmin - grid%dx / 2._dp, grid%y( j_p) - grid%dx / 2._dp]
          n_left    = grid%ij2n( i_p, j_p)
          coincides = .TRUE.

          ! Calculate the three line integrals
          CALL line_integral_xdy(   pp, p_next, grid%tol_dist, LI_xdy  )
          CALL line_integral_mxydx( pp, p_next, grid%tol_dist, LI_mxydx)
          CALL line_integral_xydy(  pp, p_next, grid%tol_dist, LI_xydy )

          ! Add them to the results structure
          CALL add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

          ! Cycle the pointer
          pp  = p_next
          j_p = j_p - 1

          ! Safety
          n_cycles = n_cycles + 1
          IF (n_cycles > grid%nx) THEN
            CALL crash('trace_line_grid - iterative tracer (west border version) got stuck!')
          END IF

        END DO ! DO WHILE (j_p > j_q)

        ! Add the last section (inside the grid cell containing q)
        CALL line_integral_xdy(   pp, qq, grid%tol_dist, LI_xdy  )
        CALL line_integral_mxydx( pp, qq, grid%tol_dist, LI_mxydx)
        CALL line_integral_xydy(  pp, qq, grid%tol_dist, LI_xydy )

        ! Add them to the results structure
        n_left = grid%ij2n( i_q, j_q)
        CALL add_integrals_to_single_row( single_row, n_left, LI_xdy, LI_mxydx, LI_xydy, coincides, count_coincidences)

      END IF ! IF (edge_index_pq == 1) THEN

    END IF ! IF (edge_index_pq == 0) THEN

  END SUBROUTINE trace_line_grid
  SUBROUTINE trace_line_grid_start( grid, p,    aij_in, bij_on, cxij_on, cyij_on)
    ! Initialise the coincidence indicators for the point p, i.e. check IF p either...
    !    - lies inside grid cell aij_in, ...
    !    - lies on the b-grid point bij_on, or...
    !    - lies on the edge cij_on

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    INTEGER,  DIMENSION(2),              INTENT(OUT)   :: aij_in, bij_on, cxij_on, cyij_on

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: xl,xu,yl,yu

    ! Initialise
    aij_in  = [0,0]
    bij_on  = [0,0]
    cxij_on = [0,0]
    cyij_on = [0,0]

    ! Find the grid cell containing p
    i = 1 + FLOOR( (p(1) - grid%xmin + grid%dx / 2._dp) / grid%dx)
    j = 1 + FLOOR( (p(2) - grid%ymin + grid%dx / 2._dp) / grid%dx)

    ! This grid cell's boundary
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! Check IF p lies on either of the four surrounding b-grid points
    IF     (i > 1       .AND. j > 1       .AND. abs( p(1) - xl) < grid%tol_dist .AND. abs( p(2) - yl) < grid%tol_dist) THEN
      ! p coincides with the southwest corner
      bij_on = [i-1,j-1]
      RETURN
    ELSEIF (i > 1       .AND. j < grid%ny .AND. abs( p(1) - xl) < grid%tol_dist .AND. abs( p(2) - yu) < grid%tol_dist) THEN
      ! p coincides with the northwest corner
      bij_on = [i-1,j  ]
      RETURN
    ELSEIF (i < grid%nx .AND. j < 1       .AND. abs( p(1) - xu) < grid%tol_dist .AND. abs( p(2) - yl) < grid%tol_dist) THEN
      ! p coincides with the southeast corner
      bij_on = [i  ,j-1]
      RETURN
    ELSEIF (i < grid%nx .AND. j < grid%ny .AND. abs( p(1) - xu) < grid%tol_dist .AND. abs( p(2) - yu) < grid%tol_dist) THEN
      ! p coincides with the northeast corner
      bij_on = [i  ,j  ]
      RETURN
    END IF

    ! Check IF p lies on any of the four borders
    IF     (i > 1       .AND. abs( p(1) - xl) < grid%tol_dist) THEN
      ! p coincides with the western border
      cxij_on = [i-1,j  ]
      RETURN
    ELSEIF (i < grid%nx .AND. abs( p(1) - xu) < grid%tol_dist) THEN
      ! p coincides with the eastern border
      cxij_on = [i  ,j  ]
      RETURN
    ELSEIF (j > 1       .AND. abs( p(2) - yl) < grid%tol_dist) THEN
      ! p coincides with the southern border
      cyij_on = [i  ,j-1]
      RETURN
    ELSEIF (j < grid%ny .AND. abs( p(2) - yu) < grid%tol_dist) THEN
      ! p coincides with the northern border
      cyij_on = [i  ,j  ]
      RETURN
    END IF

    ! p doesn't lie on the corners or borders, so it must lie inside the grid cell
    aij_in = [i,j]

  END SUBROUTINE trace_line_grid_start
  SUBROUTINE trace_line_grid_a(     grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies inside grid cell aij_in,
    ! find the point p_next where [pq] crosses into the next grid cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    INTEGER,  DIMENSION(2),              INTENT(INOUT) :: aij_in, bij_on, cxij_on, cyij_on
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(OUT)   :: n_left
    LOGICAL,                             INTENT(OUT)   :: coincides, finished

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: xl,xu,yl,yu
    REAL(dp), DIMENSION(2)                             :: sw,nw,se,ne
    LOGICAL                                            :: do_cross
    REAL(dp), DIMENSION(2)                             :: llis

    ! Safety
    IF (aij_in(1) == 0 .OR. bij_on(1) > 0 .OR. cxij_on(1) > 0 .OR. cyij_on(1) > 0) THEN
      CALL crash('trace_line_grid_a - coincidence indicators dont make sense!')
    END IF

    i = aij_in( 1)
    j = aij_in( 2)

    ! This grid cell's boundary
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! More safety
    IF (p(1) < xl .OR. p(1) > xu .OR. p(2) < yl .OR. p(2) > yu) THEN
      CALL crash('trace_line_grid_a - coincidence indicators dont make sense!')
    END IF

    ! Check if q lies inside the same grid cell
    IF (q(1) >= xl - grid%tol_dist .AND. &
        q(1) <= xu + grid%tol_dist .AND. &
        q(2) >= yl - grid%tol_dist .AND. &
        q(2) <= yu + grid%tol_dist) THEN
      ! q lies inside the same grid cell
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if pq passes through any of the four corners
    sw = [xl,yl]
    nw = [xl,yu]
    se = [xu,yl]
    ne = [xu,yu]

    IF (lies_on_line_segment( p, q, sw, grid%tol_dist)) THEN
      ! [pq] exits this grid cell through the southwest corner
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    IF (lies_on_line_segment( p, q, nw, grid%tol_dist)) THEN
      ! [pq] exits this grid cell through the northwest corner
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    IF (lies_on_line_segment( p, q, se, grid%tol_dist)) THEN
      ! [pq] exits this grid cell through the southeast corner
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    IF (lies_on_line_segment( p, q, ne, grid%tol_dist)) THEN
      ! [pq] exits this grid cell through the northeast corner
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through any of the four boundaries
    CALL segment_intersection( p, q, sw, nw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits this grid cell through the western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    CALL segment_intersection( p, q, se, ne, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits this grid cell through the eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    CALL segment_intersection( p, q, sw, se, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits this grid cell through the southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    CALL segment_intersection( p, q, nw, ne, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits this grid cell through the northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j  ]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_grid_a - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_grid_a
  SUBROUTINE trace_line_grid_b(     grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies on b-grid point bij_on
    ! find the point p_next where [pq] crosses into the next grid cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    INTEGER,  DIMENSION(2),              INTENT(INOUT) :: aij_in, bij_on, cxij_on, cyij_on
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(OUT)   :: n_left
    LOGICAL,                             INTENT(OUT)   :: coincides, finished

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: x,y,xl,xu,yl,yu
    REAL(dp), DIMENSION(2)                             :: sw,nw,se,ne,ww,ee,ss,nn
    LOGICAL                                            :: do_cross
    REAL(dp), DIMENSION(2)                             :: llis

    ! Safety
    IF (aij_in(1) > 0 .OR. bij_on(1) == 0 .OR. cxij_on(1) > 0 .OR. cyij_on(1) > 0) THEN
      CALL crash('trace_line_grid_b - coincidence indicators dont make sense!')
    END IF

    i = bij_on( 1)
    j = bij_on( 2)

    ! The eight surrounding b-grid points spanning the four surrounding a-grid cells
    x  = grid%x( i) + grid%dx / 2._dp
    y  = grid%y( j) + grid%dx / 2._dp
    xl = x - grid%dx
    xu = x + grid%dx
    yl = y - grid%dx
    yu = y + grid%dx

    sw = [xl,yl]
    ww = [xl,y ]
    nw = [xl,yu]
    ss = [x ,yl]
    nn = [x ,yu]
    se = [xu,yl]
    ee = [xu,y ]
    ne = [xu,yu]

    ! More safety
    IF (abs( p(1) - x) > grid%tol_dist .OR. abs( p(2) - y) > grid%tol_dist) THEN
      CALL crash('trace_line_grid_b - coincidence indicators dont make sense!')
    END IF

    ! Check if q lies on the cy-grid edge to the west
    IF (q(1) < x + grid%tol_dist .AND. q(1) > xl - grid%tol_dist .AND. abs( q(2) - y) < grid%tol_dist) THEN
      ! q lies on the cy-grid edge to the west
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the cy-grid edge to the east
    IF (q(1) > x - grid%tol_dist .AND. q(1) < xu + grid%tol_dist .AND. abs( q(2) - y) < grid%tol_dist) THEN
      ! q lies on the cy-grid edge to the east
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the cx-grid edge to the south
    IF (q(2) < y + grid%tol_dist .AND. q(2) > yl - grid%tol_dist .AND. abs( q(1) - x) < grid%tol_dist) THEN
      ! q lies on the cx-grid edge to the south
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the cx-grid edge to the north
    IF (q(2) > y - grid%tol_dist .AND. q(2) < yu + grid%tol_dist .AND. abs( q(1) - x) < grid%tol_dist) THEN
      ! q lies on the cx-grid edge to the north
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the a-grid cell to the northwest
    IF (q(1) > xl - grid%tol_dist .AND. q(1) < x  + grid%tol_dist .AND. &
        q(2) > y  - grid%tol_dist .AND. q(2) < yu + grid%tol_dist) THEN
      ! q lies inside the a-grid cell to the northwest
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the a-grid cell to the northeast
    IF (q(1) > x  - grid%tol_dist .AND. q(1) < xu + grid%tol_dist .AND. &
        q(2) > y  - grid%tol_dist .AND. q(2) < yu + grid%tol_dist) THEN
      ! q lies inside the a-grid cell to the northeast
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the a-grid cell to the southeast
    IF (q(1) > x  - grid%tol_dist .AND. q(1) < xu + grid%tol_dist .AND. &
        q(2) > yl - grid%tol_dist .AND. q(2) < y  + grid%tol_dist) THEN
      ! q lies inside the a-grid cell to the southeast
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j  )
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the a-grid cell to the southwest
    IF (q(1) > xl - grid%tol_dist .AND. q(1) < x  + grid%tol_dist .AND. &
        q(2) > yl - grid%tol_dist .AND. q(2) < y  + grid%tol_dist) THEN
      ! q lies inside the a-grid cell to the southwest
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i ,j  )
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the west
    IF (lies_on_line_segment( p, q, ww, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the west
      p_next    = ww
      aij_in    = [0,0]
      bij_on    = [i-1,j]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i ,j  )
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the north
    IF (lies_on_line_segment( p, q, nn, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the west
      p_next    = nn
      aij_in    = [0,0]
      bij_on    = [i,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i  ,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the east
    IF (lies_on_line_segment( p, q, ee, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the east
      p_next    = ee
      aij_in    = [0,0]
      bij_on    = [i+1,j]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the south
    IF (lies_on_line_segment( p, q, ss, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the south
      p_next    = ss
      aij_in    = [0,0]
      bij_on    = [i,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j  )
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the northwest through its western boundary
    CALL segment_intersection( p, q, ww, nw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the northwest through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the northwest through its northern boundary
    CALL segment_intersection( p, q, nw, nn, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the northwest through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j+1]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the northeast through its northern boundary
    CALL segment_intersection( p, q, nn, ne, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the northeast through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j+1]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the northeast through its eastern boundary
    CALL segment_intersection( p, q, ne, ee, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the northeast through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the southeast through its eastern boundary
    CALL segment_intersection( p, q, ee, se, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the southeast through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the southeast through its southern boundary
    CALL segment_intersection( p, q, se, ss, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the southeast through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j-1]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the southwest through its southern boundary
    CALL segment_intersection( p, q, ss, sw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the southwest through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the southwest through its western boundary
    CALL segment_intersection( p, q, sw, ww, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the southwest through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_grid_b - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_grid_b
  SUBROUTINE trace_line_grid_cx(    grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies on cx-grid edge cxij_on
    ! find the point p_next where [pq] crosses into the next grid cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    INTEGER,  DIMENSION(2),              INTENT(INOUT) :: aij_in, bij_on, cxij_on, cyij_on
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(OUT)   :: n_left
    LOGICAL,                             INTENT(OUT)   :: coincides, finished

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: x,yl,yu
    REAL(dp), DIMENSION(2)                             :: sw,nw,se,ne,ss,nn
    LOGICAL                                            :: do_cross
    REAL(dp), DIMENSION(2)                             :: llis

    ! Safety
    IF (aij_in(1) > 0 .OR. bij_on(1) > 0 .OR. cxij_on(1) == 0 .OR. cyij_on(1) > 0) THEN
      CALL crash('trace_line_grid_cx - coincidence indicators dont make sense!')
    END IF

    i = cxij_on( 1)
    j = cxij_on( 2)

    ! This c-grid edge
    x  = grid%x( i) + grid%dx / 2._dp
    yl = grid%y( j) - grid%dx / 2._dp
    yu = grid%y( j) + grid%dx / 2._dp

    ! The b-grid points spanning the two adjacent a-grid cells
    sw = [x - grid%dx, yl]
    nw = [x - grid%dx, yu]
    ss = [x          , yl]
    nn = [x          , yu]
    se = [x + grid%dx, yl]
    ne = [x + grid%dx, yu]

    ! More safety
    IF (p(2) < yl .OR. p(2) > yu .OR. ABS( p(1) - x) > grid%tol_dist) THEN
      CALL crash('trace_line_grid_cx - coincidence indicators dont make sense!')
    END IF

    ! Check IF q lies on the same cx-grid cell in the southern direction
    IF (q(2) < p(2) .AND. q(2) >= yl - grid%tol_dist .AND. ABS( q(1) - x) < grid%tol_dist) THEN
      ! q lies on the same cx-grid cell in the southern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies on the same cx-grid cell in the northern direction
    IF (q(2) > p(2) .AND. q(2) <= yu + grid%tol_dist .AND. ABS( q(1) - x) < grid%tol_dist) THEN
      ! q lies on the same cx-grid cell in the northern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies inside the grid cell to the west
    IF (q(2) >= yl - grid%tol_dist .AND. q(2) <= yu + grid%tol_dist .AND. &
        q(1) >= x - grid%dx - grid%tol_dist .AND. q(1) <= x + grid%tol_dist) THEN
      ! q lies inside the grid cell to the west
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF q lies inside the grid cell to the east
    IF (q(2) >= yl - grid%tol_dist .AND. q(2) <= yu + grid%tol_dist .AND. &
        q(1) <= x + grid%dx + grid%tol_dist .AND. q(1) >= x - grid%tol_dist) THEN
      ! q lies inside the grid cell to the east
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check IF [pq] passes through the b-grid point to the south
    IF (lies_on_line_segment( p, q, ss, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the south
      p_next    = ss
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] passes through the b-grid point to the north
    IF (lies_on_line_segment( p, q, nn, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the north
      p_next    = nn
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the west through the b-grid point to the northwest
    IF (lies_on_line_segment( p, q, nw, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the west through the b-grid point to the northwest
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the west through the b-grid point to the southwest
    IF (lies_on_line_segment( p, q, sw, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the west through the b-grid point to the southwest
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the east through the b-grid point to the northeast
    IF (lies_on_line_segment( p, q, ne, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the west through the b-grid point to the northeast
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i+1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the east through the b-grid point to the southeast
    IF (lies_on_line_segment( p, q, se, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the west through the b-grid point to the southeast
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i+1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the west through its southern boundary
    CALL segment_intersection( p, q, ss, sw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the west through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i  ,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the west through its western boundary
    CALL segment_intersection( p, q, sw, nw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the west through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the west through its northern boundary
    CALL segment_intersection( p, q, nw, nn, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the west through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the east through its northern boundary
    CALL segment_intersection( p, q, nn, ne, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the east through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the east through its eastern boundary
    CALL segment_intersection( p, q, ne, se, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the east through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i+1,j]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check IF [pq] exits the a-grid cell to the east through its southern boundary
    CALL segment_intersection( p, q, se, ss, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the east through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i+1,j-1]
      n_left    = grid%ij2n( i+1,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_grid_cx - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_grid_cx
  SUBROUTINE trace_line_grid_cy(    grid, p, q, aij_in, bij_on, cxij_on, cyij_on, p_next, n_left, coincides, finished)
    ! Given the line [pq], where p lies on cy-grid edge cyij_on
    ! find the point p_next where [pq] crosses into the next grid cell.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p,q
    INTEGER,  DIMENSION(2),              INTENT(INOUT) :: aij_in, bij_on, cxij_on, cyij_on
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: p_next
    INTEGER,                             INTENT(OUT)   :: n_left
    LOGICAL,                             INTENT(OUT)   :: coincides, finished

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: xl,xu,y
    REAL(dp), DIMENSION(2)                             :: sw,nw,se,ne,ww,ee
    LOGICAL                                            :: do_cross
    REAL(dp), DIMENSION(2)                             :: llis

    ! Safety
    IF (aij_in(1) > 0 .OR. bij_on(1) > 0 .OR. cxij_on(1) > 0 .OR. cyij_on(1) == 0) THEN
      CALL crash('trace_line_grid_cy - coincidence indicators dont make sense!')
    END IF

    i = cyij_on( 1)
    j = cyij_on( 2)

    ! This c-grid edge
    xl = grid%x( i) - grid%dx / 2._dp
    xu = grid%x( i) + grid%dx / 2._dp
    y  = grid%y( j) + grid%dx / 2._dp

    ! The b-grid points spanning the two adjacent a-grid cells
    sw = [xl, y - grid%dx]
    se = [xu, y - grid%dx]
    ww = [xl, y          ]
    ee = [xu, y          ]
    nw = [xl, y + grid%dx]
    ne = [xu, y + grid%dx]

    ! More safety
    IF (p(1) < xl .OR. p(1) > xu .OR. abs( p(2) - y) > grid%tol_dist) THEN
      CALL crash('trace_line_grid_cy - coincidence indicators dont make sense!')
    END IF

    ! Check if q lies on the same cy-grid cell in the western direction
    IF (q(1) < p(1) .AND. q(1) >= xl - grid%tol_dist .AND. abs( q(2) - y) < grid%tol_dist) THEN
      ! q lies on the same cy-grid cell in the western direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies on the same cy-grid cell in the eastern direction
    IF (q(1) > p(1) .AND. q(1) <= xu + grid%tol_dist .AND. abs( q(2) - y) < grid%tol_dist) THEN
      ! q lies on the same cy-grid cell in the eastern direction
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .TRUE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the grid cell to the south
    IF (q(1) >= xl - grid%tol_dist .AND. q(1) <= xu + grid%tol_dist .AND. &
        q(2) >= y - grid%dx - grid%tol_dist .AND. q(2) <= y + grid%tol_dist) THEN
      ! q lies inside the grid cell to the south
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if q lies inside the grid cell to the north
    IF (q(1) >= xl - grid%tol_dist .AND. q(1) <= xu + grid%tol_dist .AND. &
        q(2) <= y + grid%dx + grid%tol_dist .AND. q(2) >= y - grid%tol_dist) THEN
      ! q lies inside the grid cell to the north
      p_next    = q
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .TRUE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the west
    IF (lies_on_line_segment( p, q, ww, grid%tol_dist))  THEN
      ! [pq] passes through the b-grid point to the west
      p_next    = ww
      aij_in    = [0,0]
      bij_on    = [i-1,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] passes through the b-grid point to the east
    IF (lies_on_line_segment( p, q, ee, grid%tol_dist)) THEN
      ! [pq] passes through the b-grid point to the east
      p_next    = ee
      aij_in    = [0,0]
      bij_on    = [i  ,j  ]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .TRUE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the north through the b-grid point to the northwest
    IF (lies_on_line_segment( p, q, nw, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the north through the b-grid point to the northwest
      p_next    = nw
      aij_in    = [0,0]
      bij_on    = [i-1,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the north through the b-grid point to the northeast
    IF (lies_on_line_segment( p, q, ne, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the north through the b-grid point to the northeast
      p_next    = ne
      aij_in    = [0,0]
      bij_on    = [i  ,j+1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the south through the b-grid point to the southwest
    IF (lies_on_line_segment( p, q, sw, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the north through the b-grid point to the southwest
      p_next    = sw
      aij_in    = [0,0]
      bij_on    = [i-1,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the south through the b-grid point to the southeast
    IF (lies_on_line_segment( p, q, se, grid%tol_dist)) THEN
      ! [pq] exits the a-grid cell to the north through the b-grid point to the southeast
      p_next    = se
      aij_in    = [0,0]
      bij_on    = [i  ,j-1]
      cxij_on   = [0,0]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the north through its western boundary
    CALL segment_intersection( p, q, ww, nw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the north through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the north through its northern boundary
    CALL segment_intersection( p, q, nw, ne, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the north through its northern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j+1]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the north through its eastern boundary
    CALL segment_intersection( p, q, ne, ee, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the north through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j+1]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j+1)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the south through its eastern boundary
    CALL segment_intersection( p, q, ee, se, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the south through its eastern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i  ,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the south through its southern boundary
    CALL segment_intersection( p, q, se, sw, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the south through its southern boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [0,0]
      cyij_on   = [i,j-1]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! Check if [pq] exits the a-grid cell to the south through its western boundary
    CALL segment_intersection( p, q, sw, ww, llis, do_cross, grid%tol_dist)
    IF (do_cross) THEN
      ! [pq] exits the a-grid cell to the south through its western boundary
      p_next    = llis
      aij_in    = [0,0]
      bij_on    = [0,0]
      cxij_on   = [i-1,j  ]
      cyij_on   = [0,0]
      n_left    = grid%ij2n( i,j)
      coincides = .FALSE.
      finished  = .FALSE.
      RETURN
    END IF

    ! This point should not be reachable!
    CALL crash('trace_line_grid_cy - reached the unreachable end of the subroutine!')

  END SUBROUTINE trace_line_grid_cy

  SUBROUTINE crop_line_to_domain( p, q, xmin, xmax, ymin, ymax, tol_dist, pp, qq, is_valid_line)
    ! Crop the line [pq] so that it lies within the specified domain;
    ! if [pq] doesn't pass through the domain at all, return is_valid_line = .FALSE.

    IMPLICIT NONE

    ! In/output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p, q
    REAL(dp),                            INTENT(IN)    :: xmin, xmax, ymin, ymax, tol_dist
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: pp, qq
    LOGICAL,                             INTENT(OUT)   :: is_valid_line

    ! Local variables:
    REAL(dp), DIMENSION(2)                             :: sw, se, nw, ne, llis
    INTEGER                                            :: edge_index_p, edge_index_q
    LOGICAL                                            :: do_cross

    pp = p
    qq = q
    is_valid_line = .TRUE.

    sw = [xmin,ymin]
    se = [xmax,ymin]
    nw = [xmin,ymax]
    ne = [xmax,ymax]

    ! Determine in which quadrants p and q lie
    ! (same as with edge_index; 1-8 clockwise starting north, 0 means inside

    IF     (pp(1) >= xmin .AND. pp(1) <= xmax .AND. pp(2) > ymax) THEN
      ! North
      edge_index_p = 1
    ELSEIF (pp(1) > xmax .AND. pp(2) > ymax) THEN
      ! Northeast
      edge_index_p = 2
    ELSEIF (pp(1) > xmax .AND. pp(2) >= ymin .AND. pp(2) <= ymax) THEN
      ! East
      edge_index_p = 3
    ELSEIF (pp(1) > xmax .AND. pp(2) < ymin) THEN
      ! Southeast
      edge_index_p = 4
    ELSEIF (pp(1) >= xmin .AND. pp(1) <= xmax .AND. pp(2) < ymin) THEN
      ! South
      edge_index_p = 5
    ELSEIF (pp(1) < xmin .AND. pp(2) < ymin) THEN
      ! Southwest
      edge_index_p = 6
    ELSEIF (pp(1) < xmin .AND. pp(2) >= ymin .AND. pp(2) <= ymax) THEN
      ! West
      edge_index_p = 7
    ELSEIF (pp(1) < xmin .AND. pp(2) > ymax) THEN
      ! Northwest
      edge_index_p = 8
    ELSE
      ! Inside the mesh domain
      edge_index_p = 0
    END IF

    IF     (qq(1) >= xmin .AND. qq(1) <= xmax .AND. qq(2) > ymax) THEN
      ! North
      edge_index_q = 1
    ELSEIF (qq(1) > xmax .AND. qq(2) > ymax) THEN
      ! Northeast
      edge_index_q = 2
    ELSEIF (qq(1) > xmax .AND. qq(2) >= ymin .AND. qq(2) <= ymax) THEN
      ! East
      edge_index_q = 3
    ELSEIF (qq(1) > xmax .AND. qq(2) < ymin) THEN
      ! Southeast
      edge_index_q = 4
    ELSEIF (qq(1) >= xmin .AND. qq(1) <= xmax .AND. qq(2) < ymin) THEN
      ! South
      edge_index_q = 5
    ELSEIF (qq(1) < xmin .AND. qq(2) < ymin) THEN
      ! Southwest
      edge_index_q = 6
    ELSEIF (qq(1) < xmin .AND. qq(2) >= ymin .AND. qq(2) <= ymax) THEN
      ! West
      edge_index_q = 7
    ELSEIF (qq(1) < xmin .AND. qq(2) > ymax) THEN
      ! Northwest
      edge_index_q = 8
    ELSE
      ! Inside the mesh domain
      edge_index_q = 0
    END IF

    IF (edge_index_p == 0 .AND. edge_index_q == 0) THEN
      ! Both p and q lie inside the mesh domain
      RETURN
    END IF

    IF (edge_index_p == 0 .AND. edge_index_q > 0) THEN
      ! p lies inside the mesh domain, q lies outside

      ! Check IF [pq] passes through any of the four corners
      IF     (lies_on_line_segment( pp, qq, sw, tol_dist)) THEN
        ! [pq] passes through the southwest corner of the mesh
        qq = sw
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, se, tol_dist)) THEN
        ! [pq] passes through the southeast corner of the mesh
        qq = se
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, nw, tol_dist)) THEN
        ! [pq] passes through the northwest corner of the mesh
        qq = nw
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, ne, tol_dist)) THEN
        ! [pq] passes through the northeast corner of the mesh
        qq = ne
        RETURN
      END IF

      ! Check IF [pq] crosses any of the four borders

      ! South
      CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the southern border
        qq = llis
        RETURN
      END IF

      ! West
      CALL segment_intersection( pp, qq, sw, nw, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the western border
        qq = llis
        RETURN
      END IF

      ! North
      CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the northern border
        qq = llis
        RETURN
      END IF

      ! East
      CALL segment_intersection( pp, qq, se, ne, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the eastern border
        qq = llis
        RETURN
      END IF

      ! This point should be unreachable
      CALL crash('crop_line_to_mesh_domain - reached the unreachable point (p inside, q outside)!')

    END IF ! IF (edge_index_p == 0 .AND. edge_index_q > 0)

    IF (edge_index_q == 0 .AND. edge_index_p > 0) THEN
      ! q lies inside the mesh domain, p lies outside

      ! Check IF [pq] passes through any of the four corners
      IF     (lies_on_line_segment( pp, qq, sw, tol_dist)) THEN
        ! [pq] passes through the southwest corner of the mesh
        pp = sw
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, se, tol_dist)) THEN
        ! [pq] passes through the southeast corner of the mesh
        pp = se
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, nw, tol_dist)) THEN
        ! [pq] passes through the northwest corner of the mesh
        pp = nw
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, ne, tol_dist)) THEN
        ! [pq] passes through the northeast corner of the mesh
        pp = ne
        RETURN
      END IF

      ! Check IF [pq] crosses any of the four borders

      ! South
      CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the southern border
        pp = llis
        RETURN
      END IF

      ! West
      CALL segment_intersection( pp, qq, sw, nw, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the western border
        pp = llis
        RETURN
      END IF

      ! North
      CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the northern border
        pp = llis
        RETURN
      END IF

      ! East
      CALL segment_intersection( pp, qq, se, ne, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the eastern border
        pp = llis
        RETURN
      END IF

      ! This point should be unreachable
      CALL crash('crop_line_to_mesh_domain - reached the unreachable point (q inside, p outside)!')

    END IF ! IF (edge_index_q == 0 .AND. edge_index_p > 0)

    ! Both p and q lie outside the mesh domain

    IF     (pp(1) < xmin .AND. qq(1) < xmin) THEN
      ! Both p and q lie west of the western mesh border; [pq] cannot pass through the mesh domain
      is_valid_line = .FALSE.
      RETURN
    ELSEIF (pp(1) > xmax .AND. qq(1) > xmax) THEN
      ! Both p and q lie east of the eastern mesh border; [pq] cannot pass through the mesh domain
      is_valid_line = .FALSE.
      RETURN
    ELSEIF (pp(2) < ymin .AND. qq(2) < ymin) THEN
      ! Both p and q lie south of the southern mesh border; [pq] cannot pass through the mesh domain
      is_valid_line = .FALSE.
      RETURN
    ELSEIF (pp(2) > ymax .AND. qq(2) > ymax) THEN
      ! Both p and q lie north of the northern mesh border; [pq] cannot pass through the mesh domain
      is_valid_line = .FALSE.
      RETURN
    END IF

    IF (edge_index_p == 1) THEN
      ! p lies in the northern quadrant

      IF (edge_index_q == 3) THEN
        ! q lies in the eastern quadrant; check IF [pq] cuts through the northeast corner

        IF (cross2( (ne - qq), (pp - qq)) > 0) THEN
          ! [pq] cuts through the northeast corner
          CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect nw-ne, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, ne, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect ne-se, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSEIF (edge_index_q == 7) THEN
        ! q lies in the western quadrant; check IF [pq] cuts through the northwest corner

        IF (cross2( (pp - qq), (ne - qq)) > 0) THEN
          ! [pq] cuts through the northeast corner
          CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect nw-ne, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, nw, sw, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect nw-sw, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSE
        CALL crash('crop_line_to_mesh_domain - edge_index_p = {int_01}, edge_index_q = {int_02}', int_01 = edge_index_p, int_02 = edge_index_q)
      END IF

    ELSEIF (edge_index_p == 3) THEN
      ! p lies in the eastern quadrant

      IF (edge_index_q == 1) THEN
        ! q lies in the northern quadrant; check IF [pq] cuts through the northeast corner

        IF (cross2( (ne - pp), (qq - pp)) > 0) THEN
          ! [pq] cuts through the northeast corner
          CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect nw-ne, but it doesnt!')
          END IF
          qq = llis
          CALL segment_intersection( pp, qq, ne, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect ne-se, but it doesnt!')
          END IF
          pp = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSEIF (edge_index_q == 5) THEN
        ! q lies in the southern quadrant; check IF [pq] cuts through the southeast corner

        IF (cross2( (se - qq), (pp - qq)) > 0) THEN
          ! [pq] cuts through the southeast corner
          CALL segment_intersection( pp, qq, se, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect se-ne, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-se, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSE
        CALL crash('crop_line_to_mesh_domain - edge_index_p = {int_01}, edge_index_q = {int_02}', int_01 = edge_index_p, int_02 = edge_index_q)
      END IF

    ELSEIF (edge_index_p == 5) THEN
      ! p lies in the southern quadrant

      IF (edge_index_q == 3) THEN
        ! q lies in the eastern quadrant; check IF [pq] cuts through the southeast corner

        IF (cross2( (se - pp), (qq - pp)) > 0) THEN
          ! [pq] cuts through the southwest corner
          CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-se, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, se, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect se-ne, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSEIF (edge_index_q == 7) THEN
        ! q lies in the western quadrant; check IF [pq] cuts through the southwest corner

        IF (cross2( (qq - pp), (sw - pp)) > 0) THEN
          ! [pq] cuts through the southwest corner
          CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-se, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, sw, nw, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-nw, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSE
        CALL crash('crop_line_to_mesh_domain - edge_index_p = {int_01}, edge_index_q = {int_02}', int_01 = edge_index_p, int_02 = edge_index_q)
      END IF

    ELSEIF (edge_index_p == 7) THEN
      ! p lies in the western quadrant

      IF (edge_index_q == 5) THEN
        ! q lies in the southern quadrant; check IF [pq] cuts through the southwest corner

        IF (cross2( (pp - qq), (sw - qq)) > 0) THEN
          ! [pq] cuts through the southwest corner
          CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-se, but it doesnt!')
          END IF
          qq = llis
          CALL segment_intersection( pp, qq, sw, nw, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-nw, but it doesnt!')
          END IF
          pp = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSEIF (edge_index_q == 1) THEN
        ! q lies in the northern quadrant; check IF [pq] cuts through the northwest corner

        IF (cross2( (qq - pp), (nw - pp)) > 0) THEN
          ! [pq] cuts through the northwest corner
          CALL segment_intersection( pp, qq, sw, nw, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-nw, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect nw-ne, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSE
        CALL crash('crop_line_to_mesh_domain - edge_index_p = {int_01}, edge_index_q = {int_02}', int_01 = edge_index_p, int_02 = edge_index_q)
      END IF

    END IF ! IF (edge_index_p == 1)

    ! This point should be unreachable
    CALL crash('crop_line_to_mesh_domain - reached the unreachable end of the subroutine!')

  END SUBROUTINE crop_line_to_domain

! == Clean up after yourself
  SUBROUTINE deallocate_remapping_operators_mesh2grid( grid)
    ! Deallocate the remapping operators between the mesh and the square grid

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_remapping_operators_mesh2grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL MatDestroy( grid%M_map_mesh2grid, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_remapping_operators_mesh2grid
  SUBROUTINE deallocate_remapping_operators_grid2mesh( grid)
    ! Deallocate the remapping operators between the mesh and the square grid

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_remapping_operators_grid2mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL MatDestroy( grid%M_map_grid2mesh, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_remapping_operators_grid2mesh
  SUBROUTINE deallocate_remapping_operators_mesh_mesh( map)
    ! Deallocate the remapping operators between two meshes

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_remapping_mesh_mesh),      INTENT(INOUT) :: map

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_remapping_operators_mesh_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL MatDestroy( map%M_trilin           , perr)
    CALL MatDestroy( map%M_nearest_neighbour, perr)
    CALL MatDestroy( map%M_cons_1st_order   , perr)
    CALL MatDestroy( map%M_cons_2nd_order   , perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_remapping_operators_mesh_mesh

END MODULE mesh_mapping_module
