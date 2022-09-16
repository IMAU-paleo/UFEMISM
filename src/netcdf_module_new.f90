MODULE netcdf_module_new

  ! Contains all the subroutines for reading, creating, and writing to NetCDF files.

! ===== Preamble =====
! ====================

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
  USE data_types_module,               ONLY: type_debug_fields, type_mesh, type_grid, type_grid_lonlat, type_remapping_lonlat2mesh, &
                                             type_remapping_mesh_mesh
  USE data_types_netcdf_module_new
  USE netcdf,                          ONLY: NF90_NOERR, NF90_OPEN, NF90_CLOSE, NF90_NOWRITE, NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
                                             NF90_INQ_VARID, NF90_INQUIRE_VARIABLE, NF90_MAX_VAR_DIMS, NF90_GET_VAR, &
                                             NF90_CREATE, NF90_NOCLOBBER, NF90_NETCDF4
  USE mesh_memory_module,              ONLY: allocate_mesh_primary, allocate_mesh_secondary, deallocate_mesh_all
  USE mesh_help_functions_module,      ONLY: calc_triangle_geometric_centres, find_Voronoi_cell_areas, calc_lat_lon_coordinates, &
                                             find_triangle_areas, find_connection_widths, determine_mesh_resolution, check_mesh, &
                                             find_Voronoi_cell_geometric_centres
  USE mesh_ArakawaC_module,            ONLY: make_Ac_mesh
  USE mesh_operators_module,           ONLY: calc_matrix_operators_mesh
  USE mesh_mapping_module,             ONLY: calc_remapping_operator_grid2mesh, deallocate_remapping_operators_grid2mesh, &
                                             map_grid2mesh_2D, create_remapping_arrays_lonlat_mesh, map_lonlat2mesh_2D, &
                                             deallocate_remapping_arrays_lonlat_mesh, calc_remapping_operators_mesh_mesh, &
                                             remap_field_dp_2D, deallocate_remapping_operators_mesh_mesh
  USE utilities_module,                ONLY: transpose_dp_2D, flip_x_dp_2D, flip_y_dp_2D

  USE netcdf_module, ONLY: debug, write_to_debug_file

  IMPLICIT NONE

  INTEGER                 :: nerr
!  TYPE(type_debug_fields) :: debug_NAM, debug_EAS, debug_GRL, debug_ANT, debug

  ! Possible names for different dimensions and variables
  ! =====================================================

  ! Different options for the name of a dimension or variable can now be tried.
  ! They are separated by a double vertical bar ||

  ! Dimensions
  CHARACTER(LEN=256)      :: dim_name_options_x              = 'x||X||x1||X1||nx||NX||x-coordinate||X-coordinate||easting||Easting'
  CHARACTER(LEN=256)      :: dim_name_options_y              = 'y||Y||y1||Y1||ny||NY||y-coordinate||Y-coordinate||northing||Northing'
  CHARACTER(LEN=256)      :: dim_name_options_lon            = 'lon||Lon||long||Long||longitude||Longitude'
  CHARACTER(LEN=256)      :: dim_name_options_lat            = 'lat||Lat||latitude||Latitude'
  CHARACTER(LEN=256)      :: dim_name_options_time           = 'time||Time||t||nt'

  CHARACTER(LEN=256)      :: var_name_options_x              = 'x||X||x1||X1||nx||NX||x-coordinate||X-coordinate||easting||Easting'
  CHARACTER(LEN=256)      :: var_name_options_y              = 'y||Y||y1||Y1||ny||NY||y-coordinate||Y-coordinate||northing||Northing'
  CHARACTER(LEN=256)      :: var_name_options_lon            = 'lon||Lon||long||Long||longitude||Longitude'
  CHARACTER(LEN=256)      :: var_name_options_lat            = 'lat||Lat||latitude||Latitude'
  CHARACTER(LEN=256)      :: var_name_options_time           = 'time||Time||t||nt'

  ! Variables
  CHARACTER(LEN=256)      :: var_name_options_Hi             = 'Hi||thickness||lithk'
  CHARACTER(LEN=256)      :: var_name_options_Hb             = 'Hb||bed||topg'
  CHARACTER(LEN=256)      :: var_name_options_Hs             = 'Hs||surface||orog'

CONTAINS

  SUBROUTINE NetCDF_test_01( mesh, region_name)
    ! A simple test of the new NetCDF functionality

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_test_01'
    CHARACTER(LEN=256)                                 :: filename, field_name_options
    REAL(dp), DIMENSION(:    ), POINTER                ::  d
    INTEGER                                            :: wd

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name) // '...')

    CALL allocate_shared_dp_1D( mesh%nV, d, wd)

    filename = '/Users/berends/Documents/Models/UFEMISM/tools/matlab/tijn/BedMachine_Greenland_v4_20km_xflip.nc'
    field_name_options = var_name_options_Hi
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name)
    IF (par%master) THEN
      debug%dp_2D_a_01 = d
      CALL write_to_debug_file
    END IF
    CALL sync

    filename = '/Users/berends/Documents/Models/UFEMISM/results_test/restart_GRL_00001.nc'
    field_name_options = var_name_options_Hi
    CALL read_field_from_file_2D( filename, field_name_options, mesh, d, region_name, time_to_read = -10._dp)
    IF (par%master) THEN
      debug%dp_2D_a_02 = d
      CALL write_to_debug_file
    END IF
    CALL sync

    CALL crash('finished!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_test_01

! ===== Top-level functions =====
! ===============================

  SUBROUTINE read_field_from_file_2D( filename, field_name_options, mesh, d, region_name, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model mesh.
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid, or a regional mesh - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_file_2D'
    LOGICAL                                            :: file_exists
    LOGICAL                                            :: has_xy_grid, has_lonlat_grid, has_mesh, has_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Find out on what kind of grid the file is defined
    CALL inquire_xy_grid(     filename, has_xy_grid    )
    CALL inquire_lonlat_grid( filename, has_lonlat_grid)
    CALL inquire_mesh(        filename, has_mesh       )

    ! Files with more than one grid are not recognised
    IF (has_xy_grid     .AND. has_lonlat_grid) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a lon/lat-grid!')
    IF (has_xy_grid     .AND. has_mesh       ) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a mesh!')
    IF (has_lonlat_grid .AND. has_mesh       ) CALL crash('file "' // TRIM( filename) // '" contains both a lon/lat-grid and a mesh!')

    ! If a want data from a specific timeframe, check if the file actually has a time dimension
    IF (PRESENT( time_to_read)) THEN
      CALL inquire_time( filename, has_time)
      IF (.NOT. has_time) CALL crash('file "' // TRIM( filename) // '" does not have a (recognisable) time dimension!')
    END IF

    ! Choose the appropriate subroutine
    IF (has_xy_grid) THEN
      CALL read_field_from_xy_file_2D(     filename, field_name_options, mesh, d             , time_to_read)
    ELSEIF (has_lonlat_grid) THEN
      CALL read_field_from_lonlat_file_2D( filename, field_name_options, mesh, d             , time_to_read)
    ELSEIF (has_mesh) THEN
      CALL read_field_from_mesh_file_2D(   filename, field_name_options, mesh, d, region_name, time_to_read)
    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid, or mesh!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_2D

! ===== Separate versions for reading from x/y-grid, lon/lat-grid, or mesh files =====
! ====================================================================================

  SUBROUTINE read_field_from_xy_file_2D( filename, field_name_options, mesh, d, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model mesh.
    !
    ! Assumes the file contains data on an x/y-grid.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_xy_file_2D'
    TYPE(type_grid)                                    :: grid
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_grid
    INTEGER                                            :: wd_grid
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS)              ::  dims_of_var
    INTEGER                                            :: nx, ny, nt, id_dim_x, id_dim_y, id_dim_time
    CHARACTER(LEN=256)                                 :: indexing, xdir, ydir
    INTEGER                                            :: i,j,iopp,jopp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, grid)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var)

    ! If we couldn't find it, crash the model
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Inquire file dimensions
    CALL inquire_dim_multiple_options( filename, dim_name_options_x, nx, id_dim_x)
    CALL inquire_dim_multiple_options( filename, dim_name_options_y, ny, id_dim_y)
    IF (PRESENT( time_to_read)) CALL inquire_dim_multiple_options( filename, dim_name_options_time, nt, id_dim_time)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Inquire dimension info of this variable
    IF (par%master) THEN
      nerr = NF90_INQUIRE_VARIABLE( ncid, id_var, dimids = dims_of_var, ndims = ndims_of_var)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_INQUIRE_VARIABLE failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL MPI_BCAST( ndims_of_var, 1                , MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( dims_of_var , NF90_MAX_VAR_DIMS, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Check if the dimensions are correct
    IF (.NOT. PRESENT( time_to_read)) THEN
      ! We expect a field without a time dimension; check if this is the case
      IF (.NOT. (ndims_of_var == 2 .AND. ANY( dims_of_var == id_dim_x) .AND. ANY( dims_of_var == id_dim_y))) THEN
        CALL crash('variable "' // TRIM( field_name_options) // '" in file "' // TRIM( filename) // '" does not have the right dimensions!')
      END IF
    ELSE
      ! We expect a field with a time dimension; check if this is the case
      IF (.NOT. (ndims_of_var == 3 .AND. ANY( dims_of_var == id_dim_x) .AND. ANY( dims_of_var == id_dim_y) .AND. ANY( dims_of_var == id_dim_time))) THEN
        CALL crash('variable "' // TRIM( field_name_options) // '" in file "' // TRIM( filename) // '" does not have the right dimensions!')
      END IF
    END IF

    ! Determine the indexing of this field
    IF     (dims_of_var( 1) == id_dim_x .AND. dims_of_var( 2) == id_dim_y) THEN
      indexing = 'xy'
    ELSEIF (dims_of_var( 1) == id_dim_y .AND. dims_of_var( 2) == id_dim_x) THEN
      indexing = 'yx'
    ELSE
      CALL crash('variable "' // TRIM( field_name_options) // '" in file "' // TRIM( filename) // '" does not have x and y as dimensions!')
    END IF

    ! Determine dimension directions
    IF (grid%x( 2) > grid%x( 1)) THEN
      xdir = 'normal'
    ELSE
      xdir = 'reverse'
    END IF
    IF (grid%y( 2) > grid%y( 1)) THEN
      ydir = 'normal'
    ELSE
      ydir = 'reverse'
    END IF

    ! Allocate shared memory
    IF     (indexing == 'xy') THEN
      CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid, wd_grid)
    ELSEIF (indexing == 'yx') THEN
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_grid, wd_grid)
    END IF

    ! Read the data from the file
    CALL read_var_multiple_options_dp_2D( filename, field_name_options, d_grid, time_to_read)

    ! Perform necessary corrections to the gridded data

    IF (indexing == 'yx') THEN
      CALL transpose_dp_2D( d_grid, wd_grid)
    END IF

    IF (xdir == 'reverse') THEN
      CALL flip_x_dp_2D( d_grid)
      IF (par%master) THEN
        DO i = 1, grid%nx
          iopp = grid%nx + 1 - i
          IF (iopp <= i) EXIT
          grid%x( i   ) = grid%x( i) + grid%x( iopp)
          grid%x( iopp) = grid%x( i) - grid%x( iopp)
          grid%x( i   ) = grid%x( i) - grid%x( iopp)
        END DO
      END IF ! IF (par%master) THEN
      CALL sync
    END IF ! IF (xdir == 'reverse') THEN

    IF (ydir == 'reverse') THEN
      CALL flip_y_dp_2D( d_grid)
      IF (par%master) THEN
        DO j = 1, grid%ny
          jopp = grid%ny + 1 - j
          IF (jopp <= j) EXIT
          grid%y( j   ) = grid%y( j) + grid%y( jopp)
          grid%y( jopp) = grid%y( j) - grid%y( jopp)
          grid%y( j   ) = grid%y( j) - grid%y( jopp)
        END DO
      END IF ! IF (par%master) THEN
      CALL sync
    END IF ! IF (ydir == 'reverse') THEN

    ! Calculate remapping operator for this grid
    CALL calc_remapping_operator_grid2mesh( grid, mesh)

    ! Map the data from the grid to the mesh
    CALL map_grid2mesh_2D( grid, mesh, d_grid, d)

    ! Clean up after yourself
    CALL deallocate_shared( grid%wnx      )
    CALL deallocate_shared( grid%wny      )
    CALL deallocate_shared( grid%wx       )
    CALL deallocate_shared( grid%wy       )
    CALL deallocate_shared( grid%wxmin    )
    CALL deallocate_shared( grid%wxmax    )
    CALL deallocate_shared( grid%wymin    )
    CALL deallocate_shared( grid%wymax    )
    CALL deallocate_shared( grid%wdx      )
    CALL deallocate_shared( grid%wtol_dist)
    CALL deallocate_shared( grid%wn       )
    CALL deallocate_shared( grid%wij2n    )
    CALL deallocate_shared( grid%wn2ij    )
    CALL deallocate_remapping_operators_grid2mesh( grid)
    CALL deallocate_shared( wd_grid       )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_xy_file_2D

  SUBROUTINE read_field_from_lonlat_file_2D( filename, field_name_options, mesh, d, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model mesh.
    !
    ! Assumes the file contains data on a lon/lat-grid.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_lonlat_file_2D'
    TYPE(type_grid_lonlat)                             :: grid_lonlat
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_grid
    INTEGER                                            :: wd_grid
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS)              ::  dims_of_var
    INTEGER                                            :: nlon, nlat, nt, id_dim_lon, id_dim_lat, id_dim_time
    CHARACTER(LEN=256)                                 :: indexing, londir, latdir
    INTEGER                                            :: i,j,iopp,jopp
    TYPE(type_remapping_lonlat2mesh)                   :: map

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, grid_lonlat)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var)

    ! If we couldn't find it, crash the model
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Inquire file dimensions
    CALL inquire_dim_multiple_options( filename, dim_name_options_lon, nlon, id_dim_lon)
    CALL inquire_dim_multiple_options( filename, dim_name_options_lat, nlat, id_dim_lat)
    IF (PRESENT( time_to_read)) CALL inquire_dim_multiple_options( filename, dim_name_options_time, nt, id_dim_time)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Inquire dimension info of this variable
    IF (par%master) THEN
      nerr = NF90_INQUIRE_VARIABLE( ncid, id_var, dimids = dims_of_var, ndims = ndims_of_var)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_INQUIRE_VARIABLE failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL MPI_BCAST( ndims_of_var, 1                , MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( dims_of_var , NF90_MAX_VAR_DIMS, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Check if the dimensions are correct
    IF (.NOT. PRESENT( time_to_read)) THEN
      ! We expect a field without a time dimension; check if this is the case
      IF (.NOT. (ndims_of_var == 2 .AND. ANY( dims_of_var == id_dim_lon) .AND. ANY( dims_of_var == id_dim_lat))) THEN
        CALL crash('variable "' // TRIM( field_name_options) // '" in file "' // TRIM( filename) // '" does not have the right dimensions!')
      END IF
    ELSE
      ! We expect a field with a time dimension; check if this is the case
      IF (.NOT. (ndims_of_var == 3 .AND. ANY( dims_of_var == id_dim_lon) .AND. ANY( dims_of_var == id_dim_lat) .AND. ANY( dims_of_var == id_dim_time))) THEN
        CALL crash('variable "' // TRIM( field_name_options) // '" in file "' // TRIM( filename) // '" does not have the right dimensions!')
      END IF
    END IF

    ! Determine the indexing of this field
    IF     (dims_of_var( 1) == id_dim_lon .AND. dims_of_var( 2) == id_dim_lat) THEN
      indexing = 'lonlat'
    ELSEIF (dims_of_var( 1) == id_dim_lat .AND. dims_of_var( 2) == id_dim_lon) THEN
      indexing = 'latlon'
    ELSE
      CALL crash('variable "' // TRIM( field_name_options) // '" in file "' // TRIM( filename) // '" does not have lon and lat as dimensions!')
    END IF

    ! Determine dimension directions
    IF (grid_lonlat%lon( 2) > grid_lonlat%lon( 1)) THEN
      londir = 'normal'
    ELSE
      londir = 'reverse'
    END IF
    IF (grid_lonlat%lat( 2) > grid_lonlat%lat( 1)) THEN
      latdir = 'normal'
    ELSE
      latdir = 'reverse'
    END IF

    ! Allocate shared memory
    IF     (indexing == 'lonlat') THEN
      CALL allocate_shared_dp_2D( grid_lonlat%nlon, grid_lonlat%nlat, d_grid, wd_grid)
    ELSEIF (indexing == 'latlon') THEN
      CALL allocate_shared_dp_2D( grid_lonlat%nlat, grid_lonlat%nlon, d_grid, wd_grid)
    END IF

    ! Read the data from the file
    CALL read_var_multiple_options_dp_2D( filename, field_name_options, d_grid, time_to_read)

    ! Perform necessary corrections to the gridded data

    IF (indexing == 'latlon') THEN
      CALL transpose_dp_2D( d_grid, wd_grid)
    END IF

    IF (londir == 'reverse') THEN
      CALL flip_x_dp_2D( d_grid)
      IF (par%master) THEN
        DO i = 1, grid_lonlat%nlon
          iopp = grid_lonlat%nlon + 1 - i
          IF (iopp <= i) EXIT
          grid_lonlat%lon( i   ) = grid_lonlat%lon( i) + grid_lonlat%lon( iopp)
          grid_lonlat%lon( iopp) = grid_lonlat%lon( i) - grid_lonlat%lon( iopp)
          grid_lonlat%lon( i   ) = grid_lonlat%lon( i) - grid_lonlat%lon( iopp)
        END DO
      END IF ! IF (par%master) THEN
      CALL sync
    END IF ! IF (xdir == 'reverse') THEN

    IF (latdir == 'reverse') THEN
      CALL flip_y_dp_2D( d_grid)
      IF (par%master) THEN
        DO j = 1, grid_lonlat%nlat
          jopp = grid_lonlat%nlat + 1 - j
          IF (jopp <= j) EXIT
          grid_lonlat%lat( j   ) = grid_lonlat%lat( j) + grid_lonlat%lat( jopp)
          grid_lonlat%lat( jopp) = grid_lonlat%lat( j) - grid_lonlat%lat( jopp)
          grid_lonlat%lat( j   ) = grid_lonlat%lat( j) - grid_lonlat%lat( jopp)
        END DO
      END IF ! IF (par%master) THEN
      CALL sync
    END IF ! IF (xdir == 'reverse') THEN

    CALL correct_longitude_shifts_and_range( filename, grid_lonlat, d_grid)

    ! Calculate remapping operator for this grid
    CALL create_remapping_arrays_lonlat_mesh( mesh, grid_lonlat, map)

    ! Map the data from the grid to the mesh
    CALL map_lonlat2mesh_2D( mesh, map, d_grid, d)

    ! Clean up after yourself
    CALL deallocate_shared( grid_lonlat%wnlon  )
    CALL deallocate_shared( grid_lonlat%wnlat  )
    CALL deallocate_shared( grid_lonlat%wlon   )
    CALL deallocate_shared( grid_lonlat%wlat   )
    CALL deallocate_shared( grid_lonlat%wlonmin)
    CALL deallocate_shared( grid_lonlat%wlonmax)
    CALL deallocate_shared( grid_lonlat%wlatmin)
    CALL deallocate_shared( grid_lonlat%wlatmax)
    CALL deallocate_shared( grid_lonlat%wdlon  )
    CALL deallocate_shared( grid_lonlat%wdlat  )
    CALL deallocate_remapping_arrays_lonlat_mesh( map)
    CALL deallocate_shared( wd_grid            )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_lonlat_file_2D

  SUBROUTINE correct_longitude_shifts_and_range( filename, grid_lonlat, d_grid)
    ! Make sure longitude is bounded between 0 and 360, and increases monotonically

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid_lonlat
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'correct_longitude_shifts_and_range'
    INTEGER                                            :: i,j
    LOGICAL                                            :: is_correct
    INTEGER                                            :: n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the grid is already correct, do nothing
    is_correct = .TRUE.
    IF (ANY( grid_lonlat%lon < 0._dp) .OR. ANY( grid_lonlat%lon > 360._dp)) is_correct = .FALSE.
    DO i = 1, grid_lonlat%nlon - 1
      IF (grid_lonlat%lon( i+1) <= grid_lonlat%lon( i)) is_correct = .FALSE.
    END DO
    CALL sync
    IF (is_correct) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Limit values to [0,360]
    IF (par%master) THEN
      DO i = 1, grid_lonlat%nlon
        IF     (grid_lonlat%lon( i) <   0._dp) THEN
          grid_lonlat%lon( i) = grid_lonlat%lon( i) + 360._dp
        ELSEIF (grid_lonlat%lon( i) > 360._dp) THEN
          grid_lonlat%lon( i) = grid_lonlat%lon( i) - 360._dp
        END IF
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Fix shifts
    n = 0
    DO i = 1, grid_lonlat%nlon-1
      IF (grid_lonlat%lon( i) > grid_lonlat%lon( i+1)) THEN
        n = i
        EXIT
      END IF
    END DO

    IF (n > 0) THEN

      ! Fix lon
      IF (par%master) THEN
        grid_lonlat%lon = [grid_lonlat%lon( n+1:grid_lonlat%nlon), grid_lonlat%lon( 1:n)]
      END IF
      CALL sync

      ! Fix data field
      DO j = grid_lonlat%j1, grid_lonlat%j2
        d_grid( :,j) = [d_grid( n+1:grid_lonlat%nlon,j), d_grid( 1:n,j)]
      END DO
      CALL sync

    END IF ! IF (n > 0) THEN

    ! The grid should now be correct
    is_correct = .TRUE.
    IF (ANY( grid_lonlat%lon < 0._dp) .OR. ANY( grid_lonlat%lon > 360._dp)) is_correct = .FALSE.
    DO i = 1, grid_lonlat%nlon - 1
      IF (grid_lonlat%lon( i+1) <= grid_lonlat%lon( i)) is_correct = .FALSE.
    END DO
    IF (.NOT. is_correct) CALL crash('something is seriously wrong with the longitude of file "' // TRIM( filename) // '"!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE correct_longitude_shifts_and_range

  SUBROUTINE read_field_from_mesh_file_2D( filename, field_name_options, mesh, d, region_name, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model mesh.
    !
    ! Assumes the file contains data on a mesh.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_mesh_file_2D'
    TYPE(type_mesh)                                    :: mesh_from_file
    INTEGER                                            :: ncid
    INTEGER                                            :: id_var
    TYPE(type_netcdf_mesh)                             :: netcdf
    INTEGER                                            :: nV, nt, id_dim_time
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION(NF90_MAX_VAR_DIMS)              ::  dims_of_var
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_mesh_from_file
    INTEGER                                            :: wd_mesh_from_file
    TYPE(type_remapping_mesh_mesh)                     :: map

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the mesh from the file
    CALL setup_mesh_from_file( filename, mesh_from_file, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var)

    ! If we couldn't find it, crash the model
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Inquire file dimensions
    CALL inquire_dim( ncid, netcdf%name_dim_vi, nV, netcdf%id_dim_vi)
    IF (PRESENT( time_to_read)) CALL inquire_dim_multiple_options( filename, dim_name_options_time, nt, id_dim_time)

    ! Inquire dimension info of this variable
    IF (par%master) THEN
      nerr = NF90_INQUIRE_VARIABLE( ncid, id_var, dimids = dims_of_var, ndims = ndims_of_var)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_INQUIRE_VARIABLE failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL MPI_BCAST( ndims_of_var, 1                , MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( dims_of_var , NF90_MAX_VAR_DIMS, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Check if the dimensions are correct
    IF (.NOT. PRESENT( time_to_read)) THEN
      ! We expect a field without a time dimension; check if this is the case
      IF (.NOT. (ndims_of_var == 1 .AND. ANY( dims_of_var == netcdf%id_dim_vi))) THEN
        CALL crash('variable "' // TRIM( field_name_options) // '" in file "' // TRIM( filename) // '" does not have the right dimensions!')
      END IF
    ELSE
      ! We expect a field with a time dimension; check if this is the case
      IF (.NOT. (ndims_of_var == 2 .AND. ANY( dims_of_var == netcdf%id_dim_vi) .AND. ANY( dims_of_var == id_dim_time))) THEN
        CALL crash('variable "' // TRIM( field_name_options) // '" in file "' // TRIM( filename) // '" does not have the right dimensions!')
      END IF
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh_from_file%nV, d_mesh_from_file, wd_mesh_from_file)

    ! Read the data from the file
    CALL read_var_multiple_options_dp_1D( filename, field_name_options, d_mesh_from_file, time_to_read)

    ! Calculate remapping operator for this mesh
    CALL calc_remapping_operators_mesh_mesh( mesh_from_file, mesh, map)

    ! Map the data from the file mesh to the model mesh
    CALL remap_field_dp_2D( mesh_from_file, mesh, map, d_mesh_from_file, wd_mesh_from_file, 'cons_2nd_order')

    ! Copy to the output data memory
    d( mesh%vi1:mesh%vi2) = d_mesh_from_file( mesh%vi1:mesh%vi2)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_mesh_all( mesh_from_file)
    CALL deallocate_remapping_operators_mesh_mesh( map)
    CALL deallocate_shared( wd_mesh_from_file)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_mesh_file_2D

! ===== Set up grids/mesh from a NetCDF file =====
! ================================================

  SUBROUTINE setup_xy_grid_from_file( filename, grid)
    ! Set up an x/y-grid from a NetCDF file
    !
    ! Assumes no memory has yet been allocated for the grid at all

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid),                     INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_xy_grid_from_file'
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    INTEGER                                            :: id_dim_x, id_dim_y, id_var_x, id_var_y
    INTEGER                                            :: i,j,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for the grid size
    CALL allocate_shared_int_0D( grid%nx, grid%wnx)
    CALL allocate_shared_int_0D( grid%ny, grid%wny)

    ! Read the x and y grid size
    CALL inquire_dim_multiple_options( filename, dim_name_options_x, grid%nx, id_dim_x)
    CALL inquire_dim_multiple_options( filename, dim_name_options_y, grid%ny, id_dim_y)

    ! Allocate memory for x and y
    CALL allocate_shared_dp_1D( grid%nx, grid%x , grid%wx )
    CALL allocate_shared_dp_1D( grid%ny, grid%y , grid%wy )

    ! Inquire variable IDs
    CALL read_var_multiple_options_dp_1D( filename, var_name_options_x, grid%x)
    CALL read_var_multiple_options_dp_1D( filename, var_name_options_y, grid%y)

    ! Allocate memory for, and calculate, some secondary grid properties
    CALL allocate_shared_dp_0D(                    grid%xmin    , grid%wxmin    )
    CALL allocate_shared_dp_0D(                    grid%xmax    , grid%wxmax    )
    CALL allocate_shared_dp_0D(                    grid%ymin    , grid%wymin    )
    CALL allocate_shared_dp_0D(                    grid%ymax    , grid%wymax    )
    CALL allocate_shared_dp_0D(                    grid%dx      , grid%wdx      )
    CALL allocate_shared_dp_0D(                    grid%tol_dist, grid%wtol_dist)
    CALL allocate_shared_int_0D(                   grid%n       , grid%wn       )
    grid%n = grid%nx * grid%ny
    CALL allocate_shared_int_2D( grid%nx, grid%ny, grid%ij2n    , grid%wij2n    )
    CALL allocate_shared_int_2D( grid%n , 2,       grid%n2ij    , grid%wn2ij    )

    ! Calculate secondary grid data
    IF (par%master) THEN

      ! Resolution
      grid%dx   = ABS( grid%x( 2) - grid%x( 1))

      ! Safety
      DO i = 1, grid%nx-1
        IF (ABS( 1._dp - ABS(grid%x( i+1) - grid%x( i)) / grid%dx) > 1E-6_dp) CALL crash('file "' // TRIM( filename) // '" has an irregular x-dimension!')
      END DO
      DO j = 1, grid%ny-1
        IF (ABS( 1._dp - ABS(grid%y( j+1) - grid%y( j)) / grid%dx) > 1E-6_dp) CALL crash('file "' // TRIM( filename) // '" has an irregular y-dimension!')
      END DO

      ! Domain size
      grid%xmin = MINVAL( grid%x)
      grid%xmax = MAXVAL( grid%x)
      grid%ymin = MINVAL( grid%y)
      grid%ymax = MAXVAL( grid%y)

      ! Tolerance; points lying within this distance of each other are treated as identical
      grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

      ! Conversion tables for grid-form vs. vector-form data
      n = 0
      DO i = 1, grid%nx
        IF (MOD(i,2) == 1) THEN
          DO j = 1, grid%ny
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = grid%ny, 1, -1
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO

    END IF ! IF (par%master) THEN
    CALL sync

    ! Set up parallelisation domains
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 13)

  END SUBROUTINE setup_xy_grid_from_file

  SUBROUTINE setup_lonlat_grid_from_file( filename, grid_lonlat)
    ! Set up a lon/lat-grid from a NetCDF file
    !
    ! Assumes no memory has yet been allocated for the grid at all

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid_lonlat

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_lonlat_grid_from_file'
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_var_lon, id_var_lat
    INTEGER                                            :: i,j,n
    REAL(dp)                                           :: dlon

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for the grid size
    CALL allocate_shared_int_0D( grid_lonlat%nlon, grid_lonlat%wnlon)
    CALL allocate_shared_int_0D( grid_lonlat%nlat, grid_lonlat%wnlat)

    ! Read the lon and lat grid size
    CALL inquire_dim_multiple_options( filename, dim_name_options_lon, grid_lonlat%nlon, id_dim_lon)
    CALL inquire_dim_multiple_options( filename, dim_name_options_lat, grid_lonlat%nlat, id_dim_lat)

    ! Allocate memory for lon and lat
    CALL allocate_shared_dp_1D( grid_lonlat%nlon, grid_lonlat%lon, grid_lonlat%wlon)
    CALL allocate_shared_dp_1D( grid_lonlat%nlat, grid_lonlat%lat, grid_lonlat%wlat)

    ! Inquire variable IDs
    CALL read_var_multiple_options_dp_1D( filename, var_name_options_lon, grid_lonlat%lon)
    CALL read_var_multiple_options_dp_1D( filename, var_name_options_lat, grid_lonlat%lat)

    ! Allocate memory for, and calculate, some secondary grid properties
    CALL allocate_shared_dp_0D(                                  grid_lonlat%lonmin    , grid_lonlat%wlonmin    )
    CALL allocate_shared_dp_0D(                                  grid_lonlat%lonmax    , grid_lonlat%wlonmax    )
    CALL allocate_shared_dp_0D(                                  grid_lonlat%latmin    , grid_lonlat%wlatmin    )
    CALL allocate_shared_dp_0D(                                  grid_lonlat%latmax    , grid_lonlat%wlatmax    )
    CALL allocate_shared_dp_0D(                                  grid_lonlat%dlon      , grid_lonlat%wdlon      )
    CALL allocate_shared_dp_0D(                                  grid_lonlat%dlat      , grid_lonlat%wdlat      )

    ! Calculate secondary grid data
    IF (par%master) THEN

      ! Resolution
      grid_lonlat%dlon = ABS( grid_lonlat%lon( 2) - grid_lonlat%lon( 1))
      grid_lonlat%dlat = ABS( grid_lonlat%lat( 2) - grid_lonlat%lat( 1))

      ! Safety
      DO i = 1, grid_lonlat%nlon - 1
        ! Check for regularity in longitude, but allow for 360-degree jumps
        dlon = MIN( ABS( grid_lonlat%lon( i+1) - grid_lonlat%lon( i)), ABS( grid_lonlat%lon( i+1) + 360._dp - grid_lonlat%lon( i)))
        IF (ABS( 1._dp - dlon / grid_lonlat%dlon) > 1E-6_dp) CALL crash('file "' // TRIM( filename) // '" has an irregular longitude dimension!')
      END DO
      DO j = 1, grid_lonlat%nlat - 1
        IF (ABS( 1._dp - ABS( grid_lonlat%lat( j+1) - grid_lonlat%lat( j)) / grid_lonlat%dlat) > 1E-6_dp) &
          CALL crash('file "' // TRIM( filename) // '" has an irregular latitude dimension!')
      END DO

      ! Domain size
      grid_lonlat%lonmin = MINVAL( grid_lonlat%lon)
      grid_lonlat%lonmax = MAXVAL( grid_lonlat%lon)
      grid_lonlat%latmin = MINVAL( grid_lonlat%lat)
      grid_lonlat%latmax = MAXVAL( grid_lonlat%lat)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Set up parallelisation domains
    CALL partition_list( grid_lonlat%nlon, par%i, par%n, grid_lonlat%i1, grid_lonlat%i2)
    CALL partition_list( grid_lonlat%nlat, par%i, par%n, grid_lonlat%j1, grid_lonlat%j2)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE setup_lonlat_grid_from_file

  SUBROUTINE setup_mesh_from_file( filename, mesh, region_name)
    ! Set up a mesh from a NetCDF file
    !
    ! Assumes no memory has yet been allocated for the mesh at all

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_mesh_from_file'
    TYPE(type_netcdf_mesh)                             :: netcdf
    INTEGER                                            :: ncid
    INTEGER                                            :: nV_mem, nTri_mem, nC_mem, n_two, n_three
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Get the mesh size from the file
    CALL inquire_dim( ncid, netcdf%name_dim_vi, nV_mem  , netcdf%id_dim_vi)
    CALL inquire_dim( ncid, netcdf%name_dim_ti, nTri_mem, netcdf%id_dim_ti)
    CALL inquire_dim( ncid, netcdf%name_dim_ci, nC_mem  , netcdf%id_dim_ci)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Allocate memory for the mesh
    CALL allocate_mesh_primary( mesh, region_name, nV_mem, nTri_mem, nC_mem)

    ! Read primary mesh data
    CALL read_var_multiple_options_dp_2D(  filename, netcdf%name_var_V             , mesh%V             )
    CALL read_var_multiple_options_int_1D( filename, netcdf%name_var_nC            , mesh%nC            )
    CALL read_var_multiple_options_int_2D( filename, netcdf%name_var_C             , mesh%C             )
    CALL read_var_multiple_options_int_1D( filename, netcdf%name_var_niTri         , mesh%niTri         )
    CALL read_var_multiple_options_int_2D( filename, netcdf%name_var_iTri          , mesh%iTri          )
    CALL read_var_multiple_options_int_1D( filename, netcdf%name_var_edge_index    , mesh%edge_index    )
    CALL read_var_multiple_options_int_2D( filename, netcdf%name_var_Tri           , mesh%Tri           )
    CALL read_var_multiple_options_dp_2D(  filename, netcdf%name_var_Tricc         , mesh%Tricc         )
    CALL read_var_multiple_options_int_2D( filename, netcdf%name_var_TriC          , mesh%TriC          )
    CALL read_var_multiple_options_int_1D( filename, netcdf%name_var_Tri_edge_index, mesh%Tri_edge_index)

    ! Calculate secondary mesh data

    IF (par%master) THEN
      mesh%nV       = mesh%nV_mem
      mesh%nTri     = mesh%nTri_mem
      mesh%xmin     = MINVAL( mesh%V( :,1))
      mesh%xmax     = MAXVAL( mesh%V( :,1))
      mesh%ymin     = MINVAL( mesh%V( :,2))
      mesh%ymax     = MAXVAL( mesh%V( :,2))
      mesh%tol_dist = ((mesh%xmax - mesh%xmin) + (mesh%ymax-mesh%ymin)) * tol / 2._dp
    END IF
    CALL sync

    ! Determine vertex and triangle domains
    CALL partition_list( mesh%nV,   par%i, par%n, mesh%vi1, mesh%vi2)
    CALL partition_list( mesh%nTri, par%i, par%n, mesh%ti1, mesh%ti2)

    ! Calculate extra mesh data
    CALL allocate_mesh_secondary(             mesh)    ! Adds  9 MPI windows
    CALL calc_triangle_geometric_centres(     mesh)
    CALL find_Voronoi_cell_areas(             mesh)
    CALL calc_lat_lon_coordinates(            mesh)
    CALL find_triangle_areas(                 mesh)
    CALL find_connection_widths(              mesh)
    CALL make_Ac_mesh(                        mesh)    ! Adds  5 MPI windows
    CALL calc_matrix_operators_mesh(          mesh)    ! Adds 42 MPI windows (6 CSR matrices, 7 windows each)
    CALL determine_mesh_resolution(           mesh)
    CALL find_Voronoi_cell_geometric_centres( mesh)

    CALL check_mesh( mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE setup_mesh_from_file

! ===== Flexible looking for dimensions and variables =====
! =========================================================

  ! Look for grids and meshes
  SUBROUTINE inquire_xy_grid( filename, has_xy_grid)
    ! Inquire if a NetCDF file contains all the dimensions and variables
    ! describing a regular x/y-grid.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_xy_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_xy_grid'
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: nx, ny
    INTEGER                                            :: id_var_x, id_var_y

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Look for x and y dimensions and variables
    CALL inquire_dim_multiple_options( filename, dim_name_options_x, nx, id_dim_x)
    CALL inquire_dim_multiple_options( filename, dim_name_options_y, ny, id_dim_y)
    CALL inquire_var_multiple_options( filename, var_name_options_x,     id_var_x)
    CALL inquire_var_multiple_options( filename, var_name_options_y,     id_var_y)

    ! Check if everything is there
    has_xy_grid = .TRUE.
    IF (id_dim_x == -1 .OR. id_dim_y == -1 .OR. id_var_x == -1 .OR. id_var_y == -1) has_xy_grid = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_xy_grid

  SUBROUTINE inquire_lonlat_grid( filename, has_lonlat_grid)
    ! Inquire if a NetCDF file contains all the dimensions and variables
    ! describing a regular lon/lat-grid.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_lonlat_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_lonlat_grid'
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_lon, id_dim_lat
    INTEGER                                            :: nlon, nlat
    INTEGER                                            :: id_var_lon, id_var_lat

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Look for x and y dimensions and variables
    CALL inquire_dim_multiple_options( filename, dim_name_options_lon, nlon, id_dim_lon)
    CALL inquire_dim_multiple_options( filename, dim_name_options_lat, nlat, id_dim_lat)
    CALL inquire_var_multiple_options( filename, var_name_options_lon,       id_var_lon)
    CALL inquire_var_multiple_options( filename, var_name_options_lat,       id_var_lat)

    ! Check if everything is there
    has_lonlat_grid = .TRUE.
    IF (id_dim_lon == -1 .OR. id_dim_lat == -1 .OR. id_var_lon == -1 .OR. id_var_lat == -1) has_lonlat_grid = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_lonlat_grid

  SUBROUTINE inquire_mesh( filename, has_mesh)
    ! Inquire if a NetCDF file contains all the dimensions and variables
    ! describing a mesh.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_mesh'
    TYPE(type_netcdf_mesh)                             :: netcdf
    INTEGER                                            :: nV, nTri, nC_mem, n_two, n_three

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, netcdf%ncid)

    ! Look for mesh dimensions and variables
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_vi   , nV     , netcdf%id_dim_vi   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_ti   , nTri   , netcdf%id_dim_ti   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_ci   , nC_mem , netcdf%id_dim_ci   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_two  , n_two  , netcdf%id_dim_two  )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_three, n_three, netcdf%id_dim_three)

    CALL inquire_var( netcdf%ncid, netcdf%name_var_V             , netcdf%id_var_V             )
    CALL inquire_var( netcdf%ncid, netcdf%name_var_Tri           , netcdf%id_var_Tri           )
    CALL inquire_var( netcdf%ncid, netcdf%name_var_nC            , netcdf%id_var_nC            )
    CALL inquire_var( netcdf%ncid, netcdf%name_var_C             , netcdf%id_var_C             )
    CALL inquire_var( netcdf%ncid, netcdf%name_var_niTri         , netcdf%id_var_niTri         )
    CALL inquire_var( netcdf%ncid, netcdf%name_var_iTri          , netcdf%id_var_iTri          )
    CALL inquire_var( netcdf%ncid, netcdf%name_var_edge_index    , netcdf%id_var_edge_index    )
    CALL inquire_var( netcdf%ncid, netcdf%name_var_Tricc         , netcdf%id_var_Tricc         )
    CALL inquire_var( netcdf%ncid, netcdf%name_var_TriC          , netcdf%id_var_TriC          )
    CALL inquire_var( netcdf%ncid, netcdf%name_var_Tri_edge_index, netcdf%id_var_Tri_edge_index)

    ! Check if everything is there
    has_mesh = .TRUE.

    IF (netcdf%id_dim_vi             == -1) has_mesh = .FALSE.
    IF (netcdf%id_dim_ti             == -1) has_mesh = .FALSE.
    IF (netcdf%id_dim_ci             == -1) has_mesh = .FALSE.
    IF (netcdf%id_dim_two            == -1) has_mesh = .FALSE.
    IF (netcdf%id_dim_three          == -1) has_mesh = .FALSE.

    IF (netcdf%id_var_V              == -1) has_mesh = .FALSE.
    IF (netcdf%id_var_Tri            == -1) has_mesh = .FALSE.
    IF (netcdf%id_var_nC             == -1) has_mesh = .FALSE.
    IF (netcdf%id_var_C              == -1) has_mesh = .FALSE.
    IF (netcdf%id_var_niTri          == -1) has_mesh = .FALSE.
    IF (netcdf%id_var_iTri           == -1) has_mesh = .FALSE.
    IF (netcdf%id_var_edge_index     == -1) has_mesh = .FALSE.
    IF (netcdf%id_var_Tricc          == -1) has_mesh = .FALSE.
    IF (netcdf%id_var_TriC           == -1) has_mesh = .FALSE.
    IF (netcdf%id_var_Tri_edge_index == -1) has_mesh = .FALSE.

    ! Close the file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_mesh

  SUBROUTINE inquire_time( filename, has_time)
    ! Inquire if a NetCDF file contains a time dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_time'
    INTEGER                                            :: ncid, nt, id_dim_time, id_var_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    has_time = .TRUE.

    ! Check if the file contains a time dimension
    CALL inquire_dim_multiple_options( filename, dim_name_options_time, nt, id_dim_time)
    IF (nt < 1 .OR. id_dim_time == -1) has_time = .FALSE.

    ! Check if the file contains a time variable
    CALL inquire_var_multiple_options( filename, var_name_options_time,     id_var_time)
    IF (id_var_time == -1) has_time = .FALSE.

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_time

  ! Look for dimensions
  SUBROUTINE inquire_dim_multiple_options( filename, dim_name_options, dim_length, id_dim)
    ! Inquire if the (open) NetCDF file indicated by port ncid contains
    ! a dimension by name of dim_name. If so, return its length and identifier.
    ! If not, return -1 for both.
    !
    ! Supports providing multiple options for the dimension name, separated by two
    ! vertical bars || e.g. if we're looking for an X-dimension, we could do something like:
    !
    ! CALL inquire_dim_multiple_options( ncid, dim_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', dim_length, id_dim)
    !
    ! IF more than one match is found, crash the model.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: dim_name_options
    INTEGER,                             INTENT(OUT)   :: dim_length
    INTEGER,                             INTENT(OUT)   :: id_dim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_dim_multiple_options'
    INTEGER                                            :: ncid
    CHARACTER(LEN=LEN( dim_name_options))              :: dim_name_options_redux, dim_name
    INTEGER                                            :: i, n_matches, dim_length_loc, id_dim_loc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Try all options provided in dim_name_options

    dim_name_options_redux = TRIM( dim_name_options)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( dim_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        dim_name = dim_name_options_redux( 1:i-1)
        dim_name_options_redux = dim_name_options_redux( i+2:LEN_TRIM( dim_name_options_redux))

      ELSE
        ! Only one option is left over

        dim_name = dim_name_options_redux
        dim_name_options_redux( 1:LEN( dim_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_dim( ncid, dim_name, dim_length_loc, id_dim_loc)

      IF (id_dim_loc == -1) THEN
        ! No dimension by this name was found; try the next option
      ELSE
        ! A dimension by this name was found; hurray!
        n_matches  = n_matches + 1
        dim_length = dim_length_loc
        id_dim     = id_dim_loc
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( dim_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional dimension names were found in the NetCDF file
      dim_length = -1
      id_dim     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided dimension names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match; hurray!
    END IF

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_dim_multiple_options

  ! Look for variables
  SUBROUTINE inquire_var_multiple_options( filename, var_name_options, id_var)
    ! Inquire if the (open) NetCDF file indicated by port ncid contains
    ! a variable by name of var_name. If so, return its identifier.
    ! If not, return -1.
    !
    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! CALL inquire_var_multiple_options( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', id_var)
    !
    ! IF more than one match is found, crash the model.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name_options
    INTEGER,                             INTENT(OUT)   :: id_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_var_multiple_options'
    INTEGER                                            :: ncid
    CHARACTER(LEN=LEN( var_name_options))              :: var_name_options_redux, var_name
    INTEGER                                            :: i, n_matches, id_var_loc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Try all options provided in var_name_options

    var_name_options_redux = TRIM( var_name_options)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( var_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        var_name = var_name_options_redux( 1:i-1)
        var_name_options_redux = var_name_options_redux( i+2:LEN_TRIM( var_name_options_redux))

      ELSE
        ! Only one option is left over

        var_name = TRIM( var_name_options_redux)
        var_name_options_redux( 1:LEN( var_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_var( ncid, var_name, id_var_loc)

      IF (id_var_loc == -1) THEN
        ! No variable by this name was found; try the next option
      ELSE
        ! A variable by this name was found; hurray!
        n_matches  = n_matches + 1
        id_var     = id_var_loc
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( var_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional variable names were found in the NetCDF file
      id_var     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided variable names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match; hurray!
    END IF

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_var_multiple_options

  SUBROUTINE read_var_multiple_options_dp_1D( filename, var_name_options, d, time_to_read)
    ! Read a variable by any of the names in var_name_options from the file.
    !
    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! CALL read_var_multiple_options_dp_1D( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', d)
    !
    ! IF more than one match is found, crash the model.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_multiple_options_dp_1D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=LEN( var_name_options))              :: var_name_options_redux, var_name
    INTEGER                                            :: i, n_matches, id_var_loc, id_var
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Try all options provided in var_name_options

    var_name_options_redux = TRIM( var_name_options)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( var_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        var_name = var_name_options_redux( 1:i-1)
        var_name_options_redux = var_name_options_redux( i+2:LEN_TRIM( var_name_options_redux))

      ELSE
        ! Only one option is left over

        var_name = TRIM( var_name_options_redux)
        var_name_options_redux( 1:LEN( var_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_var( ncid, var_name, id_var_loc)

      IF (id_var_loc == -1) THEN
        ! No variable by this name was found; try the next option
      ELSE
        ! A variable by this name was found; hurray!
        n_matches  = n_matches + 1
        id_var     = id_var_loc
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( var_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional variable names were found in the NetCDF file
      id_var     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided variable names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match; hurray!
    END IF

    ! Safety
    IF (id_var == -1) CALL crash('no match found for variable "' // TRIM( var_name_options) // '" for file "' // TRIM( filename) // '"!')

    ! If needed, determine which timeframe to read
    IF (PRESENT( time_to_read)) CALL find_timeframe( filename, time_to_read, ti)

    ! Read the variable
    IF (par%master) THEN
      IF (.NOT. PRESENT( time_to_read)) THEN
        ! Read timeless data
        nerr = NF90_GET_VAR( ncid, id_var, d)
      ELSE
        ! Read data from the determined timeframe
        nerr = NF90_GET_VAR( ncid, id_var, d, start = (/ 1, ti /) )
      END IF
      ! Safety
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_multiple_options_dp_1D

  SUBROUTINE read_var_multiple_options_dp_2D( filename, var_name_options, d, time_to_read)
    ! Read a variable by any of the names in var_name_options from the file.
    !
    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! CALL read_var_multiple_options_dp_2D( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', d)
    !
    ! IF more than one match is found, crash the model.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_multiple_options_dp_2D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=LEN( var_name_options))              :: var_name_options_redux, var_name
    INTEGER                                            :: i, n_matches, id_var_loc, id_var
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Try all options provided in var_name_options

    var_name_options_redux = TRIM( var_name_options)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( var_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        var_name = var_name_options_redux( 1:i-1)
        var_name_options_redux = var_name_options_redux( i+2:LEN_TRIM( var_name_options_redux))

      ELSE
        ! Only one option is left over

        var_name = TRIM( var_name_options_redux)
        var_name_options_redux( 1:LEN( var_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_var( ncid, var_name, id_var_loc)

      IF (id_var_loc == -1) THEN
        ! No variable by this name was found; try the next option
      ELSE
        ! A variable by this name was found; hurray!
        n_matches  = n_matches + 1
        id_var     = id_var_loc
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( var_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional variable names were found in the NetCDF file
      id_var     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided variable names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match; hurray!
    END IF

    ! Safety
    IF (id_var == -1) CALL crash('no match found for variable "' // TRIM( var_name_options) // '" for file "' // TRIM( filename) // '"!')

    ! If needed, determine which timeframe to read
    IF (PRESENT( time_to_read)) CALL find_timeframe( filename, time_to_read, ti)

    ! Read the variable
    IF (par%master) THEN
      IF (.NOT. PRESENT( time_to_read)) THEN
        ! Read timeless data
        nerr = NF90_GET_VAR( ncid, id_var, d)
      ELSE
        ! Read data from the determined timeframe
        nerr = NF90_GET_VAR( ncid, id_var, d, start = (/ 1, 1, ti /) )
      END IF
      ! Safety
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_multiple_options_dp_2D

  SUBROUTINE read_var_multiple_options_dp_3D( filename, var_name_options, d, time_to_read)
    ! Read a variable by any of the names in var_name_options from the file.
    !
    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! CALL read_var_multiple_options_dp_3D( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', d)
    !
    ! IF more than one match is found, crash the model.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name_options
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_multiple_options_dp_3D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=LEN( var_name_options))              :: var_name_options_redux, var_name
    INTEGER                                            :: i, n_matches, id_var_loc, id_var
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Try all options provided in var_name_options

    var_name_options_redux = TRIM( var_name_options)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( var_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        var_name = var_name_options_redux( 1:i-1)
        var_name_options_redux = var_name_options_redux( i+2:LEN_TRIM( var_name_options_redux))

      ELSE
        ! Only one option is left over

        var_name = TRIM( var_name_options_redux)
        var_name_options_redux( 1:LEN( var_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_var( ncid, var_name, id_var_loc)

      IF (id_var_loc == -1) THEN
        ! No variable by this name was found; try the next option
      ELSE
        ! A variable by this name was found; hurray!
        n_matches  = n_matches + 1
        id_var     = id_var_loc
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( var_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional variable names were found in the NetCDF file
      id_var     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided variable names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match; hurray!
    END IF

    ! Safety
    IF (id_var == -1) CALL crash('no match found for variable "' // TRIM( var_name_options) // '" for file "' // TRIM( filename) // '"!')

    ! If needed, determine which timeframe to read
    IF (PRESENT( time_to_read)) CALL find_timeframe( filename, time_to_read, ti)

    ! Read the variable
    IF (par%master) THEN
      IF (.NOT. PRESENT( time_to_read)) THEN
        ! Read timeless data
        nerr = NF90_GET_VAR( ncid, id_var, d)
      ELSE
        ! Read data from the determined timeframe
        nerr = NF90_GET_VAR( ncid, id_var, d, start = (/ 1, 1, 1, ti /) )
      END IF
      ! Safety
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_multiple_options_dp_3D

  SUBROUTINE read_var_multiple_options_int_1D( filename, var_name_options, d, time_to_read)
    ! Read a variable by any of the names in var_name_options from the file.
    !
    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! CALL read_var_multiple_options_int_1D( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', d)
    !
    ! IF more than one match is found, crash the model.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name_options
    INTEGER,  DIMENSION(:    ),          INTENT(INOUT) :: d
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_multiple_options_int_1D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=LEN( var_name_options))              :: var_name_options_redux, var_name
    INTEGER                                            :: i, n_matches, id_var_loc, id_var
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Try all options provided in var_name_options

    var_name_options_redux = TRIM( var_name_options)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( var_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        var_name = var_name_options_redux( 1:i-1)
        var_name_options_redux = var_name_options_redux( i+2:LEN_TRIM( var_name_options_redux))

      ELSE
        ! Only one option is left over

        var_name = TRIM( var_name_options_redux)
        var_name_options_redux( 1:LEN( var_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_var( ncid, var_name, id_var_loc)

      IF (id_var_loc == -1) THEN
        ! No variable by this name was found; try the next option
      ELSE
        ! A variable by this name was found; hurray!
        n_matches  = n_matches + 1
        id_var     = id_var_loc
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( var_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional variable names were found in the NetCDF file
      id_var     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided variable names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match; hurray!
    END IF

    ! Safety
    IF (id_var == -1) CALL crash('no match found for variable "' // TRIM( var_name_options) // '" for file "' // TRIM( filename) // '"!')

    ! If needed, determine which timeframe to read
    IF (PRESENT( time_to_read)) CALL find_timeframe( filename, time_to_read, ti)

    ! Read the variable
    IF (par%master) THEN
      IF (.NOT. PRESENT( time_to_read)) THEN
        ! Read timeless data
        nerr = NF90_GET_VAR( ncid, id_var, d)
      ELSE
        ! Read data from the determined timeframe
        nerr = NF90_GET_VAR( ncid, id_var, d, start = (/ 1, ti /) )
      END IF
      ! Safety
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_multiple_options_int_1D

  SUBROUTINE read_var_multiple_options_int_2D( filename, var_name_options, d, time_to_read)
    ! Read a variable by any of the names in var_name_options from the file.
    !
    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! CALL read_var_multiple_options_int_2D( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', d)
    !
    ! IF more than one match is found, crash the model.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name_options
    INTEGER,  DIMENSION(:,:  ),          INTENT(INOUT) :: d
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_multiple_options_int_2D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=LEN( var_name_options))              :: var_name_options_redux, var_name
    INTEGER                                            :: i, n_matches, id_var_loc, id_var
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Try all options provided in var_name_options

    var_name_options_redux = TRIM( var_name_options)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( var_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        var_name = var_name_options_redux( 1:i-1)
        var_name_options_redux = var_name_options_redux( i+2:LEN_TRIM( var_name_options_redux))

      ELSE
        ! Only one option is left over

        var_name = TRIM( var_name_options_redux)
        var_name_options_redux( 1:LEN( var_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_var( ncid, var_name, id_var_loc)

      IF (id_var_loc == -1) THEN
        ! No variable by this name was found; try the next option
      ELSE
        ! A variable by this name was found; hurray!
        n_matches  = n_matches + 1
        id_var     = id_var_loc
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( var_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional variable names were found in the NetCDF file
      id_var     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided variable names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match; hurray!
    END IF

    ! Safety
    IF (id_var == -1) CALL crash('no match found for variable "' // TRIM( var_name_options) // '" for file "' // TRIM( filename) // '"!')

    ! If needed, determine which timeframe to read
    IF (PRESENT( time_to_read)) CALL find_timeframe( filename, time_to_read, ti)

    ! Read the variable
    IF (par%master) THEN
      IF (.NOT. PRESENT( time_to_read)) THEN
        ! Read timeless data
        nerr = NF90_GET_VAR( ncid, id_var, d)
      ELSE
        ! Read data from the determined timeframe
        nerr = NF90_GET_VAR( ncid, id_var, d, start = (/ 1, 1, ti /) )
      END IF
      ! Safety
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_multiple_options_int_2D

  SUBROUTINE read_var_multiple_options_int_3D( filename, var_name_options, d, time_to_read)
    ! Read a variable by any of the names in var_name_options from the file.
    !
    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! CALL read_var_multiple_options_int_3D( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', d)
    !
    ! IF more than one match is found, crash the model.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name_options
    INTEGER,  DIMENSION(:,:,:),          INTENT(INOUT) :: d
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_multiple_options_int_3D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=LEN( var_name_options))              :: var_name_options_redux, var_name
    INTEGER                                            :: i, n_matches, id_var_loc, id_var
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Try all options provided in var_name_options

    var_name_options_redux = TRIM( var_name_options)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( var_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        var_name = var_name_options_redux( 1:i-1)
        var_name_options_redux = var_name_options_redux( i+2:LEN_TRIM( var_name_options_redux))

      ELSE
        ! Only one option is left over

        var_name = TRIM( var_name_options_redux)
        var_name_options_redux( 1:LEN( var_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_var( ncid, var_name, id_var_loc)

      IF (id_var_loc == -1) THEN
        ! No variable by this name was found; try the next option
      ELSE
        ! A variable by this name was found; hurray!
        n_matches  = n_matches + 1
        id_var     = id_var_loc
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( var_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional variable names were found in the NetCDF file
      id_var     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided variable names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match; hurray!
    END IF

    ! Safety
    IF (id_var == -1) CALL crash('no match found for variable "' // TRIM( var_name_options) // '" for file "' // TRIM( filename) // '"!')

    ! If needed, determine which timeframe to read
    IF (PRESENT( time_to_read)) CALL find_timeframe( filename, time_to_read, ti)

    ! Read the variable
    IF (par%master) THEN
      IF (.NOT. PRESENT( time_to_read)) THEN
        ! Read timeless data
        nerr = NF90_GET_VAR( ncid, id_var, d)
      ELSE
        ! Read data from the determined timeframe
        nerr = NF90_GET_VAR( ncid, id_var, d, start = (/ 1, 1, 1, ti /) )
      END IF
      ! Safety
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_multiple_options_int_3D

  ! Find timeframe
  SUBROUTINE find_timeframe( filename, time_to_read, ti)
    ! Find the timeframe in the file that is closest to the desired time.
    ! If the file has no time dimension or variable, throw an error.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    REAL(dp),                            INTENT(IN)    :: time_to_read
    INTEGER,                             INTENT(OUT)   :: ti

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_timeframe'
    INTEGER                                            :: ncid
    INTEGER                                            :: nt, id_dim_time, id_var_time
    REAL(dp), DIMENSION(:    ), POINTER                ::  time
    INTEGER                                            :: wtime
    INTEGER                                            :: tii
    REAL(dp)                                           :: dt_min

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Check if the file contains a time dimension
    CALL inquire_dim_multiple_options( filename, dim_name_options_time, nt, id_dim_time)
    IF (nt < 1 .OR. id_dim_time == -1) CALL crash('no valid time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Check if the file contains a time variable
    CALL inquire_var_multiple_options( filename, var_name_options_time, id_var_time)
    IF (id_var_time == -1) CALL crash('no valid time variable could be found in file "' // TRIM( filename) // '"!')

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( nt, time, wtime)

    ! Read time from file (technically recursive, hopefully this works!)
    CALL read_var_multiple_options_dp_1D( filename, var_name_options_time, time)

    ! Find timeframe closest to desired time
    IF (time( 1) > time_to_read) THEN
      ! Desired time beyond lower limit
      CALL warning('desired timeframe at t = {dp_01} before start of file time for file "' // TRIM( filename) // '"; reading data from t = {dp_02} instead!', dp_01 = time_to_read, dp_02 = time( 1))
      ti = 1
    ELSEIF (time( nt) < time_to_read) THEN
      ! Desired time beyond upper limit
      CALL warning('desired timeframe at t = {dp_01} after end of file time for file "' // TRIM( filename) // '"; reading data from t = {dp_02} instead!', dp_01 = time_to_read, dp_02 = time( nt))
      ti = nt
    ELSE
      ! Desired time is within the file time
      dt_min = HUGE( 1._dp)
      DO tii = 1, nt
        IF (ABS( time( tii) - time_to_read) < dt_min) THEN
          ti = tii
          dt_min = ABS( time( tii) - time_to_read)
        END IF
      END DO
      IF (dt_min > 0._dp) THEN
        CALL warning('desired timeframe at t = {dp_01} not present in file "' // TRIM( filename) // '"; reading data from closest match at t = {dp_02} instead!', dp_01 = time_to_read, dp_02 = time( ti))
      END IF
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wtime)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE find_timeframe

! ===== Basic NetCDF wrapper functions =====
! ==========================================

  ! Routines for reading existing NetCDF files
  SUBROUTINE open_existing_netcdf_file_for_reading( filename, ncid)
    ! Open the NetCDF file in the specified location for reading only,
    ! and return its identifier.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(OUT)   :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'open_netcdf_file_for_reading'
    LOGICAL                                            :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Open the NetCDF file with read-only access
    IF (par%master) THEN
      nerr = NF90_OPEN( TRIM( filename), NF90_NOWRITE, ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_OPEN failed!')
    END IF ! IF (par%master) THEN

    CALL MPI_BCAST( ncid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE open_existing_netcdf_file_for_reading

  SUBROUTINE inquire_dim( ncid, dim_name, dim_length, id_dim)
    ! Inquire if the (open) NetCDF file indicated by port ncid contains
    ! a dimension by name of dim_name. If so, return its length and identifier.
    ! If not, return -1 for both.

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: dim_name
    INTEGER,                             INTENT(OUT)   :: dim_length
    INTEGER,                             INTENT(OUT)   :: id_dim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_dim'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      ! Check if a dimension of this name exists in the file
      nerr = NF90_INQ_DIMID( ncid, dim_name, id_dim)

      IF (nerr /= NF90_NOERR) THEN
        ! If a dimension by this name does not exist, return -1 for the length and ID
        id_dim     = -1
        dim_length = -1
      ELSE
        ! If a dimension by this name exists, find its length
        nerr = NF90_INQUIRE_DIMENSION( ncid, id_dim, len = dim_length)
        IF (nerr /= NF90_NOERR) CALL crash('NF90_INQUIRE_DIMENSION failed!')
      END IF

    END IF ! IF (par%master) THEN
    CALL sync

    CALL MPI_BCAST( id_dim    , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( dim_length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_dim

  SUBROUTINE inquire_var( ncid, var_name, id_var)
    ! Inquire if the (open) NetCDF file indicated by port ncid contains
    ! a variable by name of var_name. If so, return its identifier.
    ! If not, return -1 for the identifier.

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                             INTENT(IN)    :: ncid
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    INTEGER,                             INTENT(OUT)   :: id_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_var'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      ! Check if a variable of this name exists in the file
      nerr = NF90_INQ_VARID( ncid, var_name, id_var)

      IF (nerr /= NF90_NOERR) THEN
        ! If a variable by this name does not exist, return -1 for the ID
        id_var     = -1
      END IF

    END IF ! IF (par%master) THEN
    CALL sync

    CALL MPI_BCAST( id_var, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_var

  SUBROUTINE close_netcdf_file( ncid)
    ! Close an opened NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                             INTENT(OUT)   :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'close_netcdf_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Close netCDF file:
    IF (par%master) THEN
      nerr = NF90_CLOSE( ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_CLOSE failed!')
    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE close_netcdf_file

  ! Routines for creating new NetCDF files
  SUBROUTINE create_new_netcdf_file_for_writing( filename, ncid)
    ! Create a new NetCDF file in the specified location for writing.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(OUT)   :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_new_netcdf_file_for_writing'
    LOGICAL                                            :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this file already exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" already exists!')
    END IF

    ! Create the NetCDF file
    IF (par%master) THEN
      nerr = NF90_CREATE( filename, IOR( NF90_NOCLOBBER, NF90_NETCDF4), ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_CREATE failed!')
    END IF ! IF (par%master) THEN

    CALL MPI_BCAST( ncid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_new_netcdf_file_for_writing

END MODULE netcdf_module_new