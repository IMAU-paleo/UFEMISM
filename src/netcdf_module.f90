MODULE netcdf_module

  ! Contains all the subroutines for reading, creating, and writing to NetCDF files.

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr, mat_petsc2CSR
  USE petscksp
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
  USE data_types_netcdf_module
  USE data_types_module,               ONLY: type_model_region, type_mesh, type_grid, type_reference_geometry, type_forcing_data, &
                                             type_debug_fields, &
                                             type_climate_snapshot_global, type_sparse_matrix_CSR_dp, &
                                             type_ocean_snapshot_global, type_highres_ocean_data, &
                                             type_restart_data, type_netcdf_resource_tracker, &
                                             type_direct_SMB_forcing_global, type_direct_climate_forcing_global, &
                                             type_direct_SMB_forcing_regional, type_direct_climate_forcing_regional, &
                                             type_SELEN_global
  USE netcdf,                          ONLY: nf90_max_var_dims, nf90_create, nf90_close, nf90_clobber, nf90_share, nf90_unlimited , &
                                             nf90_enddef, nf90_put_var, nf90_sync, nf90_def_var, nf90_int, nf90_put_att, nf90_def_dim, &
                                             nf90_open, nf90_write, nf90_inq_dimid, nf90_inquire_dimension, nf90_inquire, nf90_double, &
                                             nf90_inq_varid, nf90_inquire_variable, nf90_get_var, nf90_noerr, nf90_strerror, nf90_float
  USE mesh_mapping_module,             ONLY: map_mesh2grid_2D, map_mesh2grid_3D
  USE sparse_matrix_module,            ONLY: deallocate_matrix_CSR

  IMPLICIT NONE

  TYPE(type_debug_fields) :: debug_NAM, debug_EAS, debug_GRL, debug_ANT, debug

CONTAINS

! ===== Basic NetCDF wrapper functions =====
! ==========================================

  SUBROUTINE open_netcdf_file( filename, ncid)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: filename
    INTEGER,          INTENT(OUT) :: ncid

    ! Open netCDF file:
    CALL handle_error(nf90_open(filename, IOR(nf90_write,nf90_share), ncid))

  END SUBROUTINE open_netcdf_file

  SUBROUTINE close_netcdf_file( ncid)
    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: ncid

    ! Close netCDF file:
    CALL handle_error(nf90_close(ncid))

  END SUBROUTINE close_netcdf_file

  SUBROUTINE create_dim( ncid, dim_name, length, id_dim)
    ! Subroutine for creating netCDF dimensions more convenient:
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN) :: ncid
    CHARACTER(LEN=*),           INTENT(IN) :: dim_name
    INTEGER,                    INTENT(IN) :: length

    ! Output variables:
    INTEGER, INTENT(OUT)               :: id_dim

    CALL handle_error(nf90_def_dim(ncid,dim_name,length,id_dim))

  END SUBROUTINE create_dim

  SUBROUTINE create_int_var( ncid, var_name, id_dims, id_var, long_name, units, missing_value)
    ! Subroutine for creating netCDF variables of type nf90_int more convenient:

    ! Input variables:
    INTEGER,                      INTENT(IN)  :: ncid
    CHARACTER(LEN=*),             INTENT(IN)  :: var_name
    INTEGER, DIMENSION(:),        INTENT(IN)  :: id_dims
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: long_name
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: units
    REAL(dp),           OPTIONAL, INTENT(IN)  :: missing_value

    ! Output variables:
    INTEGER,                      INTENT(OUT) :: id_var

    CALL handle_error(nf90_def_var(ncid,var_name,nf90_int,id_dims,id_var))
    IF(PRESENT(long_name))     CALL handle_error(nf90_put_att(ncid,id_var,'long_name',long_name))
    IF(PRESENT(units))         CALL handle_error(nf90_put_att(ncid,id_var,'units',units))
    IF(PRESENT(missing_value)) CALL handle_error(nf90_put_att(ncid,id_var,'missing_value',missing_value))

  END SUBROUTINE create_int_var

  SUBROUTINE create_single_var( ncid, var_name, id_dims, id_var, long_name, units, missing_value)
    ! Subroutine for creating netCDF variables of type nf90_FLOAT more convenient:

    ! Input variables:
    INTEGER,                      INTENT(IN)  :: ncid
    CHARACTER(LEN=*),             INTENT(IN)  :: var_name
    INTEGER, DIMENSION(:),        INTENT(IN)  :: id_dims
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: long_name
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: units
    REAL(dp),           OPTIONAL, INTENT(IN)  :: missing_value

    ! Output variables:
    INTEGER,                      INTENT(OUT) :: id_var

    CALL handle_error(nf90_def_var(ncid,var_name,nf90_float,id_dims,id_var))
    IF(PRESENT(long_name))     CALL handle_error(nf90_put_att(ncid,id_var,'long_name',long_name))
    IF(PRESENT(units))         CALL handle_error(nf90_put_att(ncid,id_var,'units',units))
    IF(PRESENT(missing_value)) CALL handle_error(nf90_put_att(ncid,id_var,'missing_value',missing_value))

  END SUBROUTINE create_single_var

  SUBROUTINE create_double_var( ncid, var_name, id_dims, id_var, long_name, units, missing_value)
    ! Subroutine for creating netCDF variables of type nf90_DOUBLE more convenient:

    ! Input variables:
    INTEGER,                      INTENT(IN)  :: ncid
    CHARACTER(LEN=*),             INTENT(IN)  :: var_name
    INTEGER, DIMENSION(:),        INTENT(IN)  :: id_dims
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: long_name
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: units
    REAL(dp),           OPTIONAL, INTENT(IN)  :: missing_value

    ! Output variables:
    INTEGER,                      INTENT(OUT) :: id_var

    CALL handle_error(nf90_def_var(ncid,var_name,nf90_double,id_dims,id_var))
    IF(PRESENT(long_name))     CALL handle_error(nf90_put_att(ncid,id_var,'long_name',long_name))
    IF(PRESENT(units))         CALL handle_error(nf90_put_att(ncid,id_var,'units',units))
    IF(PRESENT(missing_value)) CALL handle_error(nf90_put_att(ncid,id_var,'missing_value',missing_value))

  END SUBROUTINE create_double_var

  SUBROUTINE inquire_dim( ncid, dim_name, dim_length, id_dim)
    ! Inquire the id of a dimension and return its length.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)  :: ncid
    CHARACTER(LEN=*),           INTENT(IN)  :: dim_name

    ! Output variables:
    INTEGER,                    INTENT(OUT) :: dim_length
    INTEGER,                    INTENT(OUT) :: id_dim

    CALL handle_error(nf90_inq_dimid(ncid,dim_name,id_dim))
    CALL handle_error(nf90_inquire_dimension(ncid, id_dim, len=dim_length))

  END SUBROUTINE inquire_dim

  SUBROUTINE inquire_int_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_int.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF (xtype /= nf90_int) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is not nf90_int!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_int_var

  SUBROUTINE inquire_single_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_DOUBLE.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF (xtype /= nf90_float) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is not nf90_float!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_single_var

  SUBROUTINE inquire_double_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_DOUBLE.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF(xtype /= nf90_double) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is not nf90_double!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_double_var

  SUBROUTINE inquire_single_or_double_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_FLOAT or nf90_DOUBLE.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF (xtype /= nf90_double .AND. xtype /= nf90_float) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is neither nf90_float nor nf90_double!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_single_or_double_var

  SUBROUTINE handle_error( stat, message)
    USE netcdf, ONLY: nf90_noerr, nf90_strerror
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN) :: stat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: message

    IF (stat /= nf90_noerr) THEN
      IF (PRESENT( message)) THEN
        CALL crash( message)
      ELSE
        CALL crash( 'netcdf error')
      END IF
    END IF

  END SUBROUTINE handle_error

! ===== Useful tools =====
! ========================

  SUBROUTINE get_grid_from_file( filename, grid)
    ! Take an unallocated grid object and fill it with data from a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename
    TYPE(type_grid),                     INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'get_grid_from_file'
    TYPE(type_netcdf_grid)                             :: netcdf
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    INTEGER                                            :: i,j,n, var_type

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire if everything we need is present in the file, and obtain the grid size
    CALL allocate_shared_int_0D( grid%nx, grid%wnx)
    CALL allocate_shared_int_0D( grid%ny, grid%wny)

    IF (par%master) THEN

      ! Open the netcdf file
      netcdf%filename = filename
      CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

      ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
      CALL inquire_dim( netcdf%ncid, netcdf%name_dim_x, grid%nx, netcdf%id_dim_x)
      CALL inquire_dim( netcdf%ncid, netcdf%name_dim_y, grid%ny, netcdf%id_dim_y)

      ! Inquire variable id's. Make sure that each variable has the correct dimensions:
      CALL inquire_single_or_double_var( netcdf%ncid, netcdf%name_var_x, (/ netcdf%id_dim_x /), netcdf%id_var_x)
      CALL inquire_single_or_double_var( netcdf%ncid, netcdf%name_var_y, (/ netcdf%id_dim_y /), netcdf%id_var_y)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Allocate memory
    CALL allocate_shared_dp_1D( grid%nx,          grid%x, grid%wx)
    CALL allocate_shared_dp_1D(          grid%ny, grid%y, grid%wy)

    ! Read x and y
    IF (par%master) THEN

      ! Read the grid data
      CALL handle_error( nf90_get_var( netcdf%ncid, netcdf%id_var_x, grid%x, start = (/ 1 /) ))
      CALL handle_error( nf90_get_var( netcdf%ncid, netcdf%id_var_y, grid%y, start = (/ 1 /) ))

      ! Close the netcdf file
      CALL close_netcdf_file( netcdf%ncid)

    END IF
    CALL sync

    ! Allocate memory for secondary data
    CALL allocate_shared_dp_0D(                    grid%dx      , grid%wdx      )
    CALL allocate_shared_dp_0D(                    grid%xmin    , grid%wxmin    )
    CALL allocate_shared_dp_0D(                    grid%xmax    , grid%wxmax    )
    CALL allocate_shared_dp_0D(                    grid%ymin    , grid%wymin    )
    CALL allocate_shared_dp_0D(                    grid%ymax    , grid%wymax    )
    CALL allocate_shared_dp_0D(                    grid%tol_dist, grid%wtol_dist)
    CALL allocate_shared_int_0D(                   grid%n       , grid%wn       )
    grid%n = grid%nx * grid%ny
    CALL allocate_shared_int_2D( grid%nx, grid%ny, grid%ij2n    , grid%wij2n    )
    CALL allocate_shared_int_2D( grid%n , 2,       grid%n2ij    , grid%wn2ij    )

    ! Calculate secondary grid data
    IF (par%master) THEN

      ! Resolution
      grid%dx   = grid%x( 2) - grid%x( 1)

      ! Domain size
      grid%xmin = grid%x( 1)
      grid%xmax = grid%x( grid%nx)
      grid%ymin = grid%y( 1)
      grid%ymax = grid%y( grid%ny)

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

    END IF
    CALL sync

    ! Parallelisation domains
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 13)

  END SUBROUTINE get_grid_from_file

! ===== Main output functions =====
! =================================

  SUBROUTINE create_output_files( region)
    ! Create a new set of output NetCDF files (restart + help_fields + debug)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_output_files'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (par%master) WRITE(0,*) '  Creating output files...'

    ! Get output file names
    CALL get_output_filenames( region)

    ! Create the files
    CALL create_restart_file_mesh(     region, region%restart_mesh)
    CALL create_restart_file_grid(     region, region%restart_grid)
    CALL create_help_fields_file_mesh( region, region%help_fields_mesh)
    CALL create_help_fields_file_grid( region, region%help_fields_grid)
    CALL create_debug_file(            region)

    ! ISMIP6 output
    IF (C%do_write_ISMIP_output) THEN
      CALL create_ISMIP_output_files( region)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_output_files

  SUBROUTINE write_to_output_files( region)
    ! Write the current model state to the existing output files

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_output_files'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,'(A,F8.3,A)') '   t = ', region%time/1e3, ' kyr - writing output...'

    CALL write_to_restart_file_mesh(     region, region%restart_mesh)
    CALL write_to_restart_file_grid(     region, region%restart_grid)
    CALL write_to_help_fields_file_mesh( region, region%help_fields_mesh)
    CALL write_to_help_fields_file_grid( region, region%help_fields_grid)
    IF (C%do_write_ISMIP_output) CALL write_to_ISMIP_output_files( region)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_output_files

  SUBROUTINE get_output_filenames( region)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'get_output_filenames'
    CHARACTER(LEN=256)          :: short_filename
    LOGICAL                     :: ex
    INTEGER                     :: n
    CHARACTER(LEN=256)          :: ns

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! restart file (mesh)
    ! ===================

    short_filename = 'restart_NAM_00001.nc'
    short_filename(9:11) = region%name
    n = 1

    INQUIRE( FILE=(TRIM(C%output_dir) // TRIM(short_filename)), EXIST=ex )

    DO WHILE (ex)

     n=n+1

     WRITE(ns,*) n
     ns = ADJUSTL(ns)

     IF (n<10) THEN
       short_filename = short_filename(1:12) // '0000' // TRIM(ns) // '.nc'
     ELSEIF (n<100) THEN
       short_filename = short_filename(1:12) // '000' // TRIM(ns) // '.nc'
     ELSEIF (n<1000) THEN
       short_filename = short_filename(1:12) // '00' // TRIM(ns) // '.nc'
     ELSEIF (n<10000) THEN
       short_filename = short_filename(1:12) // '0' // TRIM(ns) // '.nc'
     ELSE
       short_filename = short_filename(1:12) // TRIM(ns) // '.nc'
     END IF

     INQUIRE( FILE=(TRIM(C%output_dir) // TRIM(short_filename)), EXIST=ex )

    END DO

    DO n = 1, 256
      region%restart_mesh%filename(n:n) = ' '
    END DO
    region%restart_mesh%filename = TRIM(C%output_dir)//TRIM(short_filename)

    ! help_fields file (mesh)
    ! =======================

    short_filename = 'help_fields_NAM_00001.nc'
    short_filename(13:15) = region%name
    n = 1

    INQUIRE( FILE=(TRIM(C%output_dir) // TRIM(short_filename)), EXIST=ex )

    DO WHILE (ex)

     n=n+1

     WRITE(ns,*) n
     ns = ADJUSTL(ns)

     IF (n<10) THEN
       short_filename = short_filename(1:16) // '0000' // TRIM(ns) // '.nc'
     ELSEIF (n<100) THEN
       short_filename = short_filename(1:16) // '000' // TRIM(ns) // '.nc'
     ELSEIF (n<1000) THEN
       short_filename = short_filename(1:16) // '00' // TRIM(ns) // '.nc'
     ELSEIF (n<10000) THEN
       short_filename = short_filename(1:16) // '0' // TRIM(ns) // '.nc'
     ELSE
       short_filename = short_filename(1:16) // TRIM(ns) // '.nc'
     END IF

     INQUIRE( FILE=(TRIM(C%output_dir) // TRIM(short_filename)), EXIST=ex )

    END DO

    DO n = 1, 256
      region%help_fields_mesh%filename(n:n) = ' '
    END DO
    region%help_fields_mesh%filename = TRIM(C%output_dir)//TRIM(short_filename)

    ! restart file (grid)
    ! ===================

    short_filename = 'restart_grid_NAM.nc'
    short_filename(14:16) = region%name
    DO n = 1, 256
      region%restart_grid%filename(n:n) = ' '
    END DO
    region%restart_grid%filename = TRIM(C%output_dir)//TRIM(short_filename)

    ! help_fields file (grid)
    ! =======================

    short_filename = 'help_fields_grid_NAM.nc'
    short_filename(18:20) = region%name
    DO n = 1, 256
      region%help_fields_grid%filename(n:n) = ' '
    END DO
    region%help_fields_grid%filename = TRIM(C%output_dir)//TRIM(short_filename)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_output_filenames

! ===== Create and write to output NetCDF files (mesh versions) =====
! ===================================================================

  SUBROUTINE write_to_restart_file_mesh( region, netcdf)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),   INTENT(INOUT) :: region
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_restart_file_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the file for writing
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Time
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time,             region%time,                    start = (/        netcdf%ti/)))

    ! Write data

    ! Geometry
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Hi,               region%ice%Hi_a,                start = (/ 1,     netcdf%ti/)))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Hb,               region%ice%Hb_a,                start = (/ 1,     netcdf%ti/)))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Hs,               region%ice%Hs_a,                start = (/ 1,     netcdf%ti/)))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_SL,               region%ice%SL_a,                start = (/ 1,     netcdf%ti/)))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_dHb,              region%ice%dHb_a,               start = (/ 1,     netcdf%ti/)))

    ! Bed roughness
    IF (C%choice_sliding_law == 'Weertman' .OR. &
        C%choice_sliding_law == 'Tsai2015' .OR. &
        C%choice_sliding_law == 'Schoof2005') THEN

      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_beta_sq,          region%ice%beta_sq_a,           start = (/ 1,     netcdf%ti/)))

    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised' .OR. &
            C%choice_sliding_law == 'Zoet-Iverson') THEN

      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_phi_fric,         region%ice%phi_fric_a,          start = (/ 1,     netcdf%ti/)))

    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Temperature
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Ti,                 region%ice%Ti_a,                start = (/ 1, 1,  netcdf%ti/)))

    ! SMB
    IF (C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_FirnDepth,        region%SMB%FirnDepth,           start = (/ 1, 1,  netcdf%ti/)))
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_MeltPreviousYear, region%SMB%MeltPreviousYear,    start = (/ 1,     netcdf%ti/)))
      IF (C%do_SMB_IMAUITM_inversion) THEN
        CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_C_abl_constant_inv, region%SMB%C_abl_constant_inv, start = (/ 1,  netcdf%ti/)))
        CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_C_abl_Ts_inv,       region%SMB%C_abl_Ts_inv,       start = (/ 1,  netcdf%ti/)))
        CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_C_abl_Q_inv,        region%SMB%C_abl_Q_inv,        start = (/ 1,  netcdf%ti/)))
        CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_C_refr_inv,         region%SMB%C_refr_inv,         start = (/ 1,  netcdf%ti/)))
      END IF
    END IF

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

    ! Increase time frame counter
    netcdf%ti = netcdf%ti + 1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_mesh

  SUBROUTINE write_to_help_fields_file_mesh( region, netcdf)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),       INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_help_fields_file_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the file for writing
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Time
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time, region%time, start=(/ netcdf%ti/)))

    ! Write data
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_01, C%help_field_01)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_02, C%help_field_02)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_03, C%help_field_03)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_04, C%help_field_04)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_05, C%help_field_05)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_06, C%help_field_06)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_07, C%help_field_07)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_08, C%help_field_08)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_09, C%help_field_09)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_10, C%help_field_10)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_11, C%help_field_11)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_12, C%help_field_12)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_13, C%help_field_13)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_14, C%help_field_14)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_15, C%help_field_15)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_16, C%help_field_16)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_17, C%help_field_17)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_18, C%help_field_18)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_19, C%help_field_19)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_20, C%help_field_20)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_21, C%help_field_21)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_22, C%help_field_22)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_23, C%help_field_23)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_24, C%help_field_24)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_25, C%help_field_25)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_26, C%help_field_26)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_27, C%help_field_27)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_28, C%help_field_28)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_29, C%help_field_29)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_30, C%help_field_30)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_31, C%help_field_31)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_32, C%help_field_32)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_33, C%help_field_33)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_34, C%help_field_34)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_35, C%help_field_35)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_36, C%help_field_36)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_37, C%help_field_37)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_38, C%help_field_38)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_39, C%help_field_39)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_40, C%help_field_40)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_41, C%help_field_41)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_42, C%help_field_42)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_43, C%help_field_43)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_44, C%help_field_44)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_45, C%help_field_45)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_46, C%help_field_46)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_47, C%help_field_47)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_48, C%help_field_48)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_49, C%help_field_49)
    CALL write_help_field_mesh( region, netcdf, netcdf%id_help_field_50, C%help_field_50)

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

    ! Increase time frame counter
    netcdf%ti = netcdf%ti + 1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_help_fields_file_mesh

  SUBROUTINE write_help_field_mesh( region, netcdf, id_var, field_name)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(IN)    :: region
    TYPE(type_netcdf_help_fields),  INTENT(IN)    :: netcdf
    INTEGER,                        INTENT(IN)    :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_help_field_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    IF     (field_name == 'lat') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%mesh%lat, start=(/1 /) ))
    ELSEIF (field_name == 'lon') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%mesh%lon, start=(/1 /) ))

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%GHF_a, start=(/1 /) ))

    ! Fields with a time dimension
    ! ============================

    ! Mesh
    ELSEIF (field_name == 'resolution') THEN
      ! Not needed, this is already part of regular mesh data

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%Hi_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'Hb') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%Hb_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'Hs') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%Hs_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'SL') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%SL_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'dHi') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%dHi_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'dHs') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%dHs_a, start=(/1, netcdf%ti /) ))

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%Ti_a, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Cpi') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%Cpi_a, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Ki') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%Ki_a, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%Ti_a(:,C%nZ), start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%Ti_pmp_a, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%A_flow_3D_a, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%A_flow_vav_a, start=(/1, netcdf%ti /) ))

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%u_3D_a, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'v_3D') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%v_3D_a, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'u_3D_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%u_3D_b, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'v_3D_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%v_3D_b, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'w_3D') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%w_3D_a, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'u_vav') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%u_vav_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'v_vav') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%v_vav_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'u_vav_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%u_vav_b, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'v_vav_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%v_vav_b, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%uabs_vav_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'uabs_vav_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%uabs_vav_b, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'u_surf') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%u_surf_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'v_surf') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%v_surf_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'u_surf_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%u_surf_b, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'v_surf_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%v_surf_b, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%uabs_surf_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'uabs_surf_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%uabs_surf_b, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'u_base') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%u_base_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'v_base') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%v_base_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'u_base_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%u_base_b, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'v_base_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%v_base_b, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'uabs_base') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%uabs_base_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'uabs_base_b') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%uabs_base_b, start=(/1, netcdf%ti /) ))

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%climate_matrix%applied%T2m, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'T2m_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%climate_matrix%applied%T2m,2)/12._dp, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'Precip') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%climate_matrix%applied%Precip, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Precip_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%climate_matrix%applied%Precip,2)/12._dp, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%climate_matrix%applied%Wind_WE, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%climate_matrix%applied%Wind_WE,2)/12._dp, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%climate_matrix%applied%Wind_SN, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%climate_matrix%applied%Wind_SN,2)/12._dp, start=(/1, netcdf%ti /) ))

    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%SMB, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'SMB_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%SMB_year, start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%BMB%BMB_sheet, start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%BMB%BMB_shelf, start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'BMB') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%BMB%BMB, start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'Snowfall') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%Snowfall, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%SMB%Snowfall,2), start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'Rainfall') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%Rainfall, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%SMB%Rainfall,2), start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%AddedFirn, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%SMB%AddedFirn,2), start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'Refreezing') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%Refreezing, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%Refreezing_year, start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'Melt') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%Melt, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Melt_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%SMB%Melt,2), start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'Runoff') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%Runoff, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%SMB%AddedFirn,2), start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'Albedo') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%Albedo, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%SMB%Albedo,2)/12._dp, start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%FirnDepth, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%SMB%FirnDepth,2)/12._dp, start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'C_abl_constant_inv') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%C_abl_constant_inv, start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'C_abl_Ts_inv') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%C_abl_Ts_inv, start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'C_abl_Q_inv') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%C_abl_Q_inv, start=(/1,  netcdf%ti /) ))
    ELSEIF (field_name == 'C_refr_inv') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%SMB%C_refr_inv, start=(/1,  netcdf%ti /) ))

    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'mask_land') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_land_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_ocean_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'mask_lake') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_lake_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'mask_ice') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_ice_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_sheet_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_shelf_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'mask_coast') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_coast_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'mask_margin') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_margin_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'mask_gl') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_gl_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'mask_cf') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%mask_cf_a, start=(/1, netcdf%ti /) ))

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%phi_fric_a(1:region%mesh%nV), start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'tau_yield') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%tauc_a(1:region%mesh%nV), start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'beta_sq') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%beta_sq_a(1:region%mesh%nV), start=(/1, netcdf%ti /) ))

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%IsoIce, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'iso_surf') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%IsoSurf, start=(/1, netcdf%ti /) ))

    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%Hb_a - region%refgeo_PD%Hb, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%dHb_dt_a, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'dSL_dt') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%dSL_dt_a, start=(/1, netcdf%ti /) ))
    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_help_field_mesh

  SUBROUTINE create_restart_file_mesh( region, netcdf)
    ! Create a new restart NetCDF file, write the current mesh data to it.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_restart),      INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_restart_file_mesh'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, six, vii, time, zeta, month

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set time frame index to 1
    netcdf%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( netcdf%filename) // '" already exists!')
    END IF

    ! Create netCDF file
    ! WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error(nf90_create(netcdf%filename,IOR(nf90_clobber,nf90_share),netcdf%ncid))

    ! Mesh data
    ! =========

    ! Define dimensions
    CALL create_dim( netcdf%ncid, netcdf%name_dim_vi,           region%mesh%nV,          netcdf%id_dim_vi          ) ! Vertex indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ti,           region%mesh%nTri,        netcdf%id_dim_ti          ) ! Triangle indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ci,           region%mesh%nC_mem,      netcdf%id_dim_ci          ) ! Connection indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_aci,          region%mesh%nAc,         netcdf%id_dim_aci         ) ! Staggered vertex indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ciplusone,    region%mesh%nC_mem+1,    netcdf%id_dim_ciplusone   ) ! connection indices plus one (neighbour function arrays have one more column)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_two,          2,                       netcdf%id_dim_two         ) ! 2 (each vertex has an X and Y coordinates)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_three,        3,                       netcdf%id_dim_three       ) ! 3 (each triangle has three vertices)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_six,          6,                       netcdf%id_dim_six         ) ! 4 (each staggered vertex lists four regular vertices and two triangles)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_vii_transect, region%mesh%nV_transect, netcdf%id_dim_vii_transect) ! Number of vertex pairs in the transect

    ! Placeholders for the dimension ID's, for shorter code
    vi        = netcdf%id_dim_vi
    ti        = netcdf%id_dim_ti
    ci        = netcdf%id_dim_ci
    aci       = netcdf%id_dim_aci
    ciplusone = netcdf%id_dim_ciplusone
    two       = netcdf%id_dim_two
    three     = netcdf%id_dim_three
    six       = netcdf%id_dim_six
    vii       = netcdf%id_dim_vii_transect

    ! Define variables
    CALL create_double_var( netcdf%ncid, netcdf%name_var_V,                [vi,  two  ], netcdf%id_var_V,                long_name='Vertex coordinates', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Tri,              [ti,  three], netcdf%id_var_Tri,              long_name='Vertex indices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_nC,               [vi        ], netcdf%id_var_nC,               long_name='Number of connected vertices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_C,                [vi,  ci   ], netcdf%id_var_C,                long_name='Indices of connected vertices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_niTri,            [vi        ], netcdf%id_var_niTri,            long_name='Number of inverse triangles')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_iTri,             [vi,  ci   ], netcdf%id_var_iTri,             long_name='Indices of inverse triangles')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_edge_index,       [vi        ], netcdf%id_var_edge_index,       long_name='Edge index')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Tricc,            [ti,  two  ], netcdf%id_var_Tricc,            long_name='Triangle circumcenter', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_TriC,             [ti,  three], netcdf%id_var_TriC,             long_name='Triangle neighbours')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Tri_edge_index,   [ti        ], netcdf%id_var_Tri_edge_index,   long_name='Triangle edge index')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_VAc,              [aci, two  ], netcdf%id_var_VAc,              long_name='Staggered vertex coordinates', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Aci,              [aci, six  ], netcdf%id_var_Aci,              long_name='Staggered to regular vertex indices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_iAci,             [vi,  ci   ], netcdf%id_var_iAci,             long_name='Regular to staggered vertex indices')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_A,                [vi        ], netcdf%id_var_A,                long_name='Vertex Voronoi cell area', units='m^2')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_R,                [vi        ], netcdf%id_var_R,                long_name='Vertex resolution', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_vi_transect,      [vii, two  ], netcdf%id_var_vi_transect,      long_name='Transect vertex pairs')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_w_transect,       [vii, two  ], netcdf%id_var_w_transect,       long_name='Transect interpolation weights')

    ! Model output
    ! ============

    ! Define dimensions
    CALL create_dim( netcdf%ncid, netcdf%name_dim_zeta,  C%nZ,           netcdf%id_dim_zeta ) ! Scaled vertical coordinate
    CALL create_dim( netcdf%ncid, netcdf%name_dim_month, 12,             netcdf%id_dim_month) ! Months (for monthly data)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time,  nf90_unlimited, netcdf%id_dim_time ) ! Time frames

    ! Placeholders for the dimension ID's, for shorter code
    time  = netcdf%id_dim_time
    zeta  = netcdf%id_dim_zeta
    month = netcdf%id_dim_month

    ! Define dimension variables
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time,  [time  ], netcdf%id_var_time,  long_name='Time', units='years')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_zeta,  [zeta  ], netcdf%id_var_zeta,  long_name='Vertical scaled coordinate', units='unitless (0 = ice surface, 1 = bedrock)')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_month, [month ], netcdf%id_var_month, long_name='Month', units='1-12')

    ! Define model data variables

    ! Geometry
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hi,               [vi,        time], netcdf%id_var_Hi,               long_name='Ice thickness', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hb,               [vi,        time], netcdf%id_var_Hb,               long_name='Bedrock elevation', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hs,               [vi,        time], netcdf%id_var_Hs,               long_name='Surface elevation', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_SL,               [vi,        time], netcdf%id_var_SL,               long_name='Sea surface change', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_dHb,              [vi,        time], netcdf%id_var_dHb,              long_name='Bedrock deformation', units='m')

    ! Bed roughness
    CALL create_double_var( netcdf%ncid, netcdf%name_var_beta_sq,          [vi,        time], netcdf%id_var_beta_sq,          long_name='Bed roughness', units='?')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_phi_fric,         [vi,        time], netcdf%id_var_phi_fric,         long_name='Bed roughness', units='?')

    ! Temperature
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Ti,               [vi, zeta,  time], netcdf%id_var_Ti,               long_name='Ice temperature', units='K')

    ! SMB
    IF (C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_FirnDepth,        [vi, month, time], netcdf%id_var_FirnDepth,          long_name='Firn depth', units='m')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_MeltPreviousYear, [vi,        time], netcdf%id_var_MeltPreviousYear,   long_name='Melt during previous year', units='mie')
      IF (C%do_SMB_IMAUITM_inversion) THEN
        CALL create_double_var( netcdf%ncid, netcdf%name_var_C_abl_constant_inv, [vi,    time], netcdf%id_var_C_abl_constant_inv, long_name='Threshold ablation factor', units='-')
        CALL create_double_var( netcdf%ncid, netcdf%name_var_C_abl_Ts_inv,       [vi,    time], netcdf%id_var_C_abl_Ts_inv,       long_name='Temperature ablation factor', units='-')
        CALL create_double_var( netcdf%ncid, netcdf%name_var_C_abl_Q_inv,        [vi,    time], netcdf%id_var_C_abl_Q_inv,        long_name='Insolation ablation factor', units='-')
        CALL create_double_var( netcdf%ncid, netcdf%name_var_C_refr_inv,         [vi,    time], netcdf%id_var_C_refr_inv,         long_name='Refreezing factor', units='-')
      END IF
    END IF

    ! Leave definition mode
    CALL handle_error(nf90_enddef( netcdf%ncid))

    ! Write mesh data
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_V,               region%mesh%V             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tri,             region%mesh%Tri           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_nC,              region%mesh%nC            ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_C,               region%mesh%C             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_niTri,           region%mesh%niTri         ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_iTri,            region%mesh%iTri          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_edge_index,      region%mesh%edge_index    ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tricc,           region%mesh%Tricc         ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_TriC,            region%mesh%TriC          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tri_edge_index,  region%mesh%Tri_edge_index))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_VAc,             region%mesh%VAc           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Aci,             region%mesh%Aci           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_iAci,            region%mesh%iAci          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_A,               region%mesh%A             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_R,               region%mesh%R             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_vi_transect,     region%mesh%vi_transect   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_w_transect,      region%mesh%w_transect    ))

    ! Write zeta and month dimension variables
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_mesh

  SUBROUTINE create_help_fields_file_mesh( region, netcdf)
    ! Create a new help_fields NetCDF file, write the current mesh data to it.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_fields_file_mesh'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, six, vii, time, zeta, month

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set time frame index to 1
    netcdf%ti = 1

    ! Create a new help file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( netcdf%filename) // '" already exists!')
    END IF

    ! Create netCDF file
    ! WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error(nf90_create(netcdf%filename,IOR(nf90_clobber,nf90_share),netcdf%ncid))

    ! Mesh data
    ! =========

    ! Define dimensions
    CALL create_dim( netcdf%ncid, netcdf%name_dim_vi,           region%mesh%nV,          netcdf%id_dim_vi          ) ! Vertex indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ti,           region%mesh%nTri,        netcdf%id_dim_ti          ) ! Triangle indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ci,           region%mesh%nC_mem,      netcdf%id_dim_ci          ) ! Connection indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_aci,          region%mesh%nAc,         netcdf%id_dim_aci         ) ! Staggered vertex indices
    CALL create_dim( netcdf%ncid, netcdf%name_dim_ciplusone,    region%mesh%nC_mem+1,    netcdf%id_dim_ciplusone   ) ! connection indices plus one (neighbour function arrays have one more column)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_two,          2,                       netcdf%id_dim_two         ) ! 2 (each vertex has an X and Y coordinates)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_three,        3,                       netcdf%id_dim_three       ) ! 3 (each triangle has three vertices)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_six,          6,                       netcdf%id_dim_six         ) ! 4 (each staggered vertex lists four regular vertices and two triangles)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_vii_transect, region%mesh%nV_transect, netcdf%id_dim_vii_transect) ! Number of vertex pairs in the transect

    ! Placeholders for the dimension ID's, for shorter code
    vi        = netcdf%id_dim_vi
    ti        = netcdf%id_dim_ti
    ci        = netcdf%id_dim_ci
    aci       = netcdf%id_dim_aci
    ciplusone = netcdf%id_dim_ciplusone
    two       = netcdf%id_dim_two
    three     = netcdf%id_dim_three
    six       = netcdf%id_dim_six
    vii       = netcdf%id_dim_vii_transect

    ! Define variables
    CALL create_double_var( netcdf%ncid, netcdf%name_var_V,                [vi,  two  ], netcdf%id_var_V,                long_name='Vertex coordinates', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Tri,              [ti,  three], netcdf%id_var_Tri,              long_name='Vertex indices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_nC,               [vi        ], netcdf%id_var_nC,               long_name='Number of connected vertices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_C,                [vi,  ci   ], netcdf%id_var_C,                long_name='Indices of connected vertices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_niTri,            [vi        ], netcdf%id_var_niTri,            long_name='Number of inverse triangles')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_iTri,             [vi,  ci   ], netcdf%id_var_iTri,             long_name='Indices of inverse triangles')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_edge_index,       [vi        ], netcdf%id_var_edge_index,       long_name='Edge index')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Tricc,            [ti,  two  ], netcdf%id_var_Tricc,            long_name='Triangle circumcenter', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_TriC,             [ti,  three], netcdf%id_var_TriC,             long_name='Triangle neighbours')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Tri_edge_index,   [ti        ], netcdf%id_var_Tri_edge_index,   long_name='Triangle edge index')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_VAc,              [aci, two  ], netcdf%id_var_VAc,              long_name='Staggered vertex coordinates', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Aci,              [aci, six  ], netcdf%id_var_Aci,              long_name='Staggered to regular vertex indices')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_iAci,             [vi,  ci   ], netcdf%id_var_iAci,             long_name='Regular to staggered vertex indices')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_A,                [vi        ], netcdf%id_var_A,                long_name='Vertex Voronoi cell area', units='m^2')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_R,                [vi        ], netcdf%id_var_R,                long_name='Vertex resolution', units='m')
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_vi_transect,      [vii, two  ], netcdf%id_var_vi_transect,      long_name='Transect vertex pairs')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_w_transect,       [vii, two  ], netcdf%id_var_w_transect,       long_name='Transect interpolation weights')

    ! Model output
    ! ============

    ! Define dimensions
    CALL create_dim( netcdf%ncid, netcdf%name_dim_zeta,  C%nZ,           netcdf%id_dim_zeta ) ! Scaled vertical coordinate
    CALL create_dim( netcdf%ncid, netcdf%name_dim_month, 12,             netcdf%id_dim_month) ! Months (for monthly data)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time,  nf90_unlimited, netcdf%id_dim_time ) ! Time frames

    ! Placeholders for the dimension ID's, for shorter code
    time  = netcdf%id_dim_time
    zeta  = netcdf%id_dim_zeta
    month = netcdf%id_dim_month

    ! Define dimension variables
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time,  [time  ], netcdf%id_var_time,  long_name='Time', units='years')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_zeta,  [zeta  ], netcdf%id_var_zeta,  long_name='Vertical scaled coordinate', units='unitless (0 = ice surface, 1 = bedrock)')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_month, [month ], netcdf%id_var_month, long_name='Month', units='1-12')

    ! Define model data variables

    ! Define data variables
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_01, C%help_field_01)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_02, C%help_field_02)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_03, C%help_field_03)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_04, C%help_field_04)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_05, C%help_field_05)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_06, C%help_field_06)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_07, C%help_field_07)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_08, C%help_field_08)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_09, C%help_field_09)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_10, C%help_field_10)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_11, C%help_field_11)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_12, C%help_field_12)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_13, C%help_field_13)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_14, C%help_field_14)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_15, C%help_field_15)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_16, C%help_field_16)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_17, C%help_field_17)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_18, C%help_field_18)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_19, C%help_field_19)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_20, C%help_field_20)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_21, C%help_field_21)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_22, C%help_field_22)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_23, C%help_field_23)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_24, C%help_field_24)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_25, C%help_field_25)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_26, C%help_field_26)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_27, C%help_field_27)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_28, C%help_field_28)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_29, C%help_field_29)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_30, C%help_field_30)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_31, C%help_field_31)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_32, C%help_field_32)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_33, C%help_field_33)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_34, C%help_field_34)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_35, C%help_field_35)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_36, C%help_field_36)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_37, C%help_field_37)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_38, C%help_field_38)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_39, C%help_field_39)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_40, C%help_field_40)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_41, C%help_field_41)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_42, C%help_field_42)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_43, C%help_field_43)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_44, C%help_field_44)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_45, C%help_field_45)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_46, C%help_field_46)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_47, C%help_field_47)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_48, C%help_field_48)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_49, C%help_field_49)
    CALL create_help_field_mesh( netcdf, netcdf%id_help_field_50, C%help_field_50)

    ! Leave definition mode
    CALL handle_error(nf90_enddef( netcdf%ncid))

    ! Write mesh data
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_V,               region%mesh%V             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tri,             region%mesh%Tri           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_nC,              region%mesh%nC            ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_C,               region%mesh%C             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_niTri,           region%mesh%niTri         ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_iTri,            region%mesh%iTri          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_edge_index,      region%mesh%edge_index    ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tricc,           region%mesh%Tricc         ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_TriC,            region%mesh%TriC          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Tri_edge_index,  region%mesh%Tri_edge_index))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_VAc,             region%mesh%VAc           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Aci,             region%mesh%Aci           ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_iAci,            region%mesh%iAci          ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_A,               region%mesh%A             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_R,               region%mesh%R             ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_vi_transect,     region%mesh%vi_transect   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_w_transect,      region%mesh%w_transect    ))

    ! Write zeta and month dimension variables
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_fields_file_mesh

  SUBROUTINE create_help_field_mesh( netcdf, id_var, field_name)
    ! Add a data field to the help_fields file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf
    INTEGER,                        INTENT(INOUT) :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_help_field_mesh'
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, six, vii, ai, tai, t, z, m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Placeholders for the dimension ID's, for shorter code
    vi        = netcdf%id_dim_vi
    ti        = netcdf%id_dim_ti
    ci        = netcdf%id_dim_ci
    aci       = netcdf%id_dim_aci
    ciplusone = netcdf%id_dim_ciplusone
    two       = netcdf%id_dim_two
    three     = netcdf%id_dim_three
    six       = netcdf%id_dim_six
    vii       = netcdf%id_dim_vii_transect
    ai        = netcdf%id_dim_ai
    tai       = netcdf%id_dim_tai
    t         = netcdf%id_dim_time
    z         = netcdf%id_dim_zeta
    m         = netcdf%id_dim_month

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    ELSEIF (field_name == 'lat') THEN
      CALL create_double_var( netcdf%ncid, 'lat',                      [vi      ], id_var, long_name='Latitude',  units='degrees north')
    ELSEIF (field_name == 'lon') THEN
      CALL create_double_var( netcdf%ncid, 'lon',                      [vi      ], id_var, long_name='Longitude', units='degrees east')

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL create_double_var( netcdf%ncid, 'GHF',                      [vi      ], id_var, long_name='Geothermal heat flux', units='J m^-2 yr^-1')

    ! Fields with a time dimension
    ! ============================

    ! Mesh
    ELSEIF (field_name == 'resolution') THEN
      ! Not needed, this is already part of regular mesh data

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL create_double_var( netcdf%ncid, 'Hi',                       [vi,    t], id_var, long_name='Ice thickness', units='m')
    ELSEIF (field_name == 'Hb') THEN
      CALL create_double_var( netcdf%ncid, 'Hb',                       [vi,    t], id_var, long_name='Bedrock elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'Hs') THEN
      CALL create_double_var( netcdf%ncid, 'Hs',                       [vi,    t], id_var, long_name='Surface elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'SL') THEN
      CALL create_double_var( netcdf%ncid, 'SL',                       [vi,    t], id_var, long_name='Geoid elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'dHi') THEN
      CALL create_double_var( netcdf%ncid, 'dHi',                      [vi,    t], id_var, long_name='Ice thickness difference w.r.t PD', units='m')
    ELSEIF (field_name == 'dHs') THEN
      CALL create_double_var( netcdf%ncid, 'dHs',                      [vi,    t], id_var, long_name='Ice elevation difference w.r.t PD', units='m')

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL create_double_var( netcdf%ncid, 'Ti',                       [vi, z, t], id_var, long_name='Englacial temperature', units='K')
    ELSEIF (field_name == 'Cpi') THEN
      CALL create_double_var( netcdf%ncid, 'Cpi',                      [vi, z, t], id_var, long_name='Ice heat capacity', units='J kg^-1 K^-1')
    ELSEIF (field_name == 'Ki') THEN
      CALL create_double_var( netcdf%ncid, 'Ki',                       [vi, z, t], id_var, long_name='Ice thermal conductivity', units='J m^-1 K^-1 yr^-1')
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL create_double_var( netcdf%ncid, 'Ti_basal',                 [vi,    t], id_var, long_name='Ice basal temperature', units='K')
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL create_double_var( netcdf%ncid, 'Ti_pmp',                   [vi, z, t], id_var, long_name='Ice pressure melting point temperature', units='K')
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL create_double_var( netcdf%ncid, 'A_flow_3D',                [vi, z, t], id_var, long_name='Ice flow factor', units='Pa^-3 y^-1')
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL create_double_var( netcdf%ncid, 'A_flow_vav',               [vi,    t], id_var, long_name='Vertically averaged ice flow factor', units='Pa^-3 y^-1')

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL create_double_var( netcdf%ncid, 'u_3D',                     [vi, z, t], id_var, long_name='3D ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_3D') THEN
      CALL create_double_var( netcdf%ncid, 'v_3D',                     [vi, z, t], id_var, long_name='3D ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_3D_b') THEN
      CALL create_double_var( netcdf%ncid, 'u_3D_b',                   [ti, z, t], id_var, long_name='3D ice x-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'v_3D_b') THEN
      CALL create_double_var( netcdf%ncid, 'v_3D_b',                   [ti, z, t], id_var, long_name='3D ice y-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'w_3D') THEN
      CALL create_double_var( netcdf%ncid, 'w_3D',                     [vi, z, t], id_var, long_name='3D ice z-velocity', units='m/yr')
    ELSEIF (field_name == 'u_vav') THEN
      CALL create_double_var( netcdf%ncid, 'u_vav',                    [vi,    t], id_var, long_name='Vertically averaged ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_vav') THEN
      CALL create_double_var( netcdf%ncid, 'v_vav',                    [vi,    t], id_var, long_name='Vertically averaged ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'u_vav_b') THEN
      CALL create_double_var( netcdf%ncid, 'u_vav_b',                  [ti,    t], id_var, long_name='Vertically averaged ice x-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'v_vav_b') THEN
      CALL create_double_var( netcdf%ncid, 'v_vav_b',                  [ti,    t], id_var, long_name='Vertically averaged ice y-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_vav',                 [vi,    t], id_var, long_name='Vertically averaged ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_vav_b') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_vav_b',               [ti,    t], id_var, long_name='Vertically averaged ice velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'u_surf') THEN
      CALL create_double_var( netcdf%ncid, 'u_surf',                   [vi,    t], id_var, long_name='Surface ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_surf') THEN
      CALL create_double_var( netcdf%ncid, 'v_surf',                   [vi,    t], id_var, long_name='Surface ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_surf_b') THEN
      CALL create_double_var( netcdf%ncid, 'u_surf_b',                 [ti,    t], id_var, long_name='Surface ice x-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'v_surf_b') THEN
      CALL create_double_var( netcdf%ncid, 'v_surf_b',                 [ti,    t], id_var, long_name='Surface ice y-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_surf',                [vi,    t], id_var, long_name='Surface ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_surf_b') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_surf_b',              [ti,    t], id_var, long_name='Surface ice velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'u_base') THEN
      CALL create_double_var( netcdf%ncid, 'u_base',                   [vi,    t], id_var, long_name='Basal ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_base') THEN
      CALL create_double_var( netcdf%ncid, 'v_base',                   [vi,    t], id_var, long_name='Basal ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_base_b') THEN
      CALL create_double_var( netcdf%ncid, 'u_base_b',                 [ti,    t], id_var, long_name='Basal ice x-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'v_base_b') THEN
      CALL create_double_var( netcdf%ncid, 'v_base_b',                 [ti,    t], id_var, long_name='Basal ice y-velocity (b-grid)', units='m/yr')
    ELSEIF (field_name == 'uabs_base') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_base',                [vi,    t], id_var, long_name='Basal ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_base_b') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_base_b',              [ti,    t], id_var, long_name='Basal ice velocity (b-grid)', units='m/yr')

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL create_double_var( netcdf%ncid, 'T2m',                      [vi, m, t], id_var, long_name='Monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'T2m_year') THEN
      CALL create_double_var( netcdf%ncid, 'T2m_year',                 [vi,    t], id_var, long_name='Annual mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'Precip') THEN
      CALL create_double_var( netcdf%ncid, 'Precip',                   [vi, m, t], id_var, long_name='Monthly total precipitation', units='mm')
    ELSEIF (field_name == 'Precip_year') THEN
      CALL create_double_var( netcdf%ncid, 'Precip_year',              [vi,    t], id_var, long_name='Annual total precipitation', units='mm')
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_WE',                  [vi, m, t], id_var, long_name='Monthly mean zonal wind', units='m/s')
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_WE_year',             [vi,    t], id_var, long_name='Annual mean zonal wind', units='m/s')
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_SN',                  [vi, m, t], id_var, long_name='Monthly mean meridional wind', units='m/s')
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_SN_year',             [vi,    t], id_var, long_name='Annual mean meridional wind', units='m/s')

    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL create_double_var( netcdf%ncid, 'SMB',                      [vi, m, t], id_var, long_name='Monthly surface mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'SMB_year') THEN
      CALL create_double_var( netcdf%ncid, 'SMB_year',                 [vi,    t], id_var, long_name='Annual surface mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL create_double_var( netcdf%ncid, 'BMB_sheet',                [vi,    t], id_var, long_name='Annual basal mass balance for grounded ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL create_double_var( netcdf%ncid, 'BMB_shelf',                [vi,    t], id_var, long_name='Annual basal mass balance for floating ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB') THEN
      CALL create_double_var( netcdf%ncid, 'BMB',                      [vi,    t], id_var, long_name='Annual basal mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'Snowfall') THEN
      CALL create_double_var( netcdf%ncid, 'Snowfall',                 [vi, m, t], id_var, long_name='Monthly total snowfall', units='m water equivalent')
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL create_double_var( netcdf%ncid, 'Snowfall_year',            [vi,    t], id_var, long_name='Annual total snowfall', units='m water equivalent')
    ELSEIF (field_name == 'Rainfall') THEN
      CALL create_double_var( netcdf%ncid, 'Rainfall',                 [vi, m, t], id_var, long_name='Monthly total rainfall', units='m water equivalent')
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL create_double_var( netcdf%ncid, 'Rainfall_year',            [vi,    t], id_var, long_name='Annual total rainfall', units='m water equivalent')
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL create_double_var( netcdf%ncid, 'AddedFirn',                [vi, m, t], id_var, long_name='Monthly total added firn', units='m water equivalent')
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL create_double_var( netcdf%ncid, 'AddedFirn_year',           [vi,    t], id_var, long_name='Annual total added firn', units='m water equivalent')
    ELSEIF (field_name == 'Refreezing') THEN
      CALL create_double_var( netcdf%ncid, 'Refreezing',               [vi, m, t], id_var, long_name='Monthly total refreezing', units='m water equivalent')
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL create_double_var( netcdf%ncid, 'Refreezing_year',          [vi,    t], id_var, long_name='Annual total refreezing', units='m water equivalent')
    ELSEIF (field_name == 'Melt') THEN
      CALL create_double_var( netcdf%ncid, 'Melt',                     [vi, m, t], id_var, long_name='Monthly total melting', units='m water equivalent')
    ELSEIF (field_name == 'Melt_year') THEN
      CALL create_double_var( netcdf%ncid, 'Melt_year',                [vi,    t], id_var, long_name='Annual total melt', units='m water equivalent')
    ELSEIF (field_name == 'Runoff') THEN
      CALL create_double_var( netcdf%ncid, 'Runoff',                   [vi, m, t], id_var, long_name='Monthly total runoff', units='m water equivalent')
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL create_double_var( netcdf%ncid, 'Runoff_year',              [vi,    t], id_var, long_name='Annual total runoff', units='m water equivalent')
    ELSEIF (field_name == 'Albedo') THEN
      CALL create_double_var( netcdf%ncid, 'Albedo',                   [vi, m, t], id_var, long_name='Monthly mean albedo')
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL create_double_var( netcdf%ncid, 'Albedo_year',              [vi,    t], id_var, long_name='Annual mean albedo')
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL create_double_var( netcdf%ncid, 'FirnDepth',                [vi, m, t], id_var, long_name='Monthly mean firn layer depth', units='m water equivalent')
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL create_double_var( netcdf%ncid, 'FirnDepth_year',           [vi,    t], id_var, long_name='Annual mean firn layer depth', units='m water equivalent')
    ELSEIF (field_name == 'C_abl_constant_inv') THEN
      CALL create_double_var( netcdf%ncid, 'C_abl_constant_inv',       [vi,    t], id_var, long_name='Constant ablation', units='-')
    ELSEIF (field_name == 'C_abl_Ts_inv') THEN
      CALL create_double_var( netcdf%ncid, 'C_abl_Ts_inv',             [vi,    t], id_var, long_name='Temperature ablation factor', units='-')
    ELSEIF (field_name == 'C_abl_Q_inv') THEN
      CALL create_double_var( netcdf%ncid, 'C_abl_Q_inv',              [vi,    t], id_var, long_name='Insolation ablation factor', units='-')
    ELSEIF (field_name == 'C_refr_inv') THEN
      CALL create_double_var( netcdf%ncid, 'C_refr_inv',               [vi,    t], id_var, long_name='Refreezing factor', units='-')

    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL create_int_var(    netcdf%ncid, 'mask',                     [vi,    t], id_var, long_name='Mask')
    ELSEIF (field_name == 'mask_land') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_land',                [vi,    t], id_var, long_name='Land mask')
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_ocean',               [vi,    t], id_var, long_name='Ocean mask')
    ELSEIF (field_name == 'mask_lake') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_lake',                [vi,    t], id_var, long_name='Lake mask')
    ELSEIF (field_name == 'mask_ice') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_ice',                 [vi,    t], id_var, long_name='Ice mask')
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_sheet',               [vi,    t], id_var, long_name='Sheet mask')
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_shelf',               [vi,    t], id_var, long_name='Shelf mask')
    ELSEIF (field_name == 'mask_coast') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_coast',               [vi,    t], id_var, long_name='Coast mask')
    ELSEIF (field_name == 'mask_margin') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_margin',              [vi,    t], id_var, long_name='Margin mask')
    ELSEIF (field_name == 'mask_gl') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_gl',                  [vi,    t], id_var, long_name='Grounding-line mask')
    ELSEIF (field_name == 'mask_cf') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_cf',                  [vi,    t], id_var, long_name='Calving-front mask')

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL create_double_var( netcdf%ncid, 'phi_fric',                 [vi,    t], id_var, long_name='Till friction angle', units='degrees')
    ELSEIF (field_name == 'tau_yield') THEN
      CALL create_double_var( netcdf%ncid, 'tau_yield',                [vi,    t], id_var, long_name='Basal yield stress', units='Pa')
    ELSEIF (field_name == 'beta_sq') THEN
      CALL create_double_var( netcdf%ncid, 'beta_sq',                  [vi,    t], id_var, long_name='Sliding coefficient', units='?')

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL create_double_var( netcdf%ncid, 'iso_ice',                  [vi,    t], id_var, long_name='Vertically averaged ice d18O', units='per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL create_double_var( netcdf%ncid, 'iso_surf',                 [vi,    t], id_var, long_name='d18O of precipitation', units='per mille')

    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      CALL create_double_var( netcdf%ncid, 'dHb',                      [vi,    t], id_var, long_name='Change in bedrock elevation w.r.t. PD', units='m')
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL create_double_var( netcdf%ncid, 'dHb_dt',                 [vi,    t], id_var, long_name='Bedrock deformation rates', units='m/yr')
    ELSEIF (field_name == 'dSL_dt') THEN
      CALL create_double_var( netcdf%ncid, 'dSL_dt',                 [vi,    t], id_var, long_name='Geoid deformation rates', units='m/yr')
    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_field_mesh

! ===== Create and write to output NetCDF files (grid versions) =====
! ===================================================================

  SUBROUTINE write_to_restart_file_grid( region, netcdf)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),   INTENT(INOUT) :: region
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_restart_file_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the file for writing
    IF (par%master) CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Time
    IF (par%master) CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time,             region%time,      start=(/         netcdf%ti/)))

    ! Map and write data

    ! Geometry
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hi_a,             netcdf%id_var_Hi,               netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hb_a,             netcdf%id_var_Hb,               netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hs_a,             netcdf%id_var_Hs,               netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%SL_a,             netcdf%id_var_SL,               netcdf%ti      )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dHb_a,            netcdf%id_var_dHb,              netcdf%ti      )

    ! Bed roughness
    IF (C%choice_sliding_law == 'Weertman' .OR. C%choice_sliding_law == 'Tsai2015' .OR. C%choice_sliding_law == 'Schoof2005') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%beta_sq_a,      netcdf%id_var_beta_sq,          netcdf%ti      )

    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. C%choice_sliding_law == 'Coulomb_regularised' .OR. C%choice_sliding_law == 'Zoet-Iverson') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%phi_fric_a,     netcdf%id_var_phi_fric,         netcdf%ti      )

    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Temperature
    CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti_a,             netcdf%id_var_Ti,               netcdf%ti, C%nZ)

    ! SMB
    IF (C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%FirnDepth,        netcdf%id_var_FirnDepth,        netcdf%ti, 12  )
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%MeltPreviousYear, netcdf%id_var_MeltPreviousYear, netcdf%ti      )
      IF (C%do_SMB_IMAUITM_inversion) THEN
        CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%C_abl_constant_inv, netcdf%id_var_C_abl_constant_inv, netcdf%ti)
        CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%C_abl_Ts_inv,       netcdf%id_var_C_abl_Ts_inv,       netcdf%ti)
        CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%C_abl_Q_inv,        netcdf%id_var_C_abl_Q_inv,        netcdf%ti)
        CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%C_refr_inv,         netcdf%id_var_C_refr_inv,         netcdf%ti)
      END IF
    END IF

    ! Close the file
    IF (par%master) CALL close_netcdf_file(netcdf%ncid)

    ! Increase time frame counter
    IF (par%master) netcdf%ti = netcdf%ti + 1
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_grid

  SUBROUTINE write_to_help_fields_file_grid( region, netcdf)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),       INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_help_fields_file_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_help_fields_file_grid

  SUBROUTINE write_help_field_grid( region, netcdf, id_var, field_name)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),       INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields), INTENT(INOUT) :: netcdf
    INTEGER,                       INTENT(IN)    :: id_var
    CHARACTER(LEN=*),              INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'write_help_field_grid'
    INTEGER                                      :: vi
    REAL(dp), DIMENSION(:    ), POINTER          ::  dp_2D_a
    INTEGER                                      :: wdp_2D_a

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( region%mesh%nV, dp_2D_a, wdp_2D_a)

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    IF     (field_name == 'lat') THEN
      IF (par%master) CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%grid_output%lat ))
    ELSEIF (field_name == 'lon') THEN
      IF (par%master) CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%grid_output%lon ))

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%GHF_a, id_var, netcdf%ti)

    ! Fields with a time dimension
    ! ============================

    ! Mesh
    ELSEIF (field_name == 'resolution') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%mesh%R, id_var, netcdf%ti)

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hi_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Hb') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hb_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Hs') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hs_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'SL') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%SL_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'dHi') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dHi_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'dHs') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dHs_a, id_var, netcdf%ti)

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'Cpi') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Cpi_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'Ki') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ki_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti_a(:,C%nZ), id_var, netcdf%ti)
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti_pmp_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%A_flow_3D_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%A_flow_vav_a, id_var, netcdf%ti)

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%u_3D_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'v_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%v_3D_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'u_3D_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_3D_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'w_3D') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%w_3D_a, id_var, netcdf%ti, C%nZ)
    ELSEIF (field_name == 'u_vav') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%u_vav_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'v_vav') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%v_vav_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'u_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%uabs_vav_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'uabs_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'u_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%u_surf_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'v_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%v_surf_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'u_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%uabs_surf_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'uabs_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'u_base') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%u_base_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'v_base') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%v_base_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'u_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_base') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%uabs_base_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'uabs_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate_matrix%applied%T2m, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'T2m_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate_matrix%applied%T2m( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Precip') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate_matrix%applied%Precip, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Precip_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate_matrix%applied%Precip( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate_matrix%applied%Wind_WE, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Wind_WE_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate_matrix%applied%Wind_WE( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate_matrix%applied%Wind_SN, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Wind_SN_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate_matrix%applied%Wind_SN( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)

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
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Snowfall( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Rainfall') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Rainfall, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Rainfall_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Rainfall( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%AddedFirn, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'AddedFirn_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%AddedFirn( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Refreezing') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Refreezing, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Refreezing_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Refreezing( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Melt') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Melt, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Melt_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Melt( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Runoff') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Runoff, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Runoff_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Runoff( vi,:))
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Albedo') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%Albedo, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Albedo_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%Albedo( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%FirnDepth, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'FirnDepth_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%SMB%FirnDepth( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'C_abl_constant_inv') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%C_abl_constant_inv, id_var, netcdf%ti)
    ELSEIF (field_name == 'C_abl_Ts_inv') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%C_abl_Ts_inv, id_var, netcdf%ti)
    ELSEIF (field_name == 'C_abl_Q_inv') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%C_abl_Q_inv, id_var, netcdf%ti)
    ELSEIF (field_name == 'C_refr_inv') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%C_refr_inv, id_var, netcdf%ti)

    ! Masks
    ! NOTE: not meant to be included, as mapping masks between grids is meaningless. These lines are needed so the model can output
    !       the mesh version of these variables without crashing at the ELSE below.
    ELSEIF (field_name == 'mask'       ) THEN
    ELSEIF (field_name == 'mask_land'  ) THEN
    ELSEIF (field_name == 'mask_ocean' ) THEN
    ELSEIF (field_name == 'mask_lake'  ) THEN
    ELSEIF (field_name == 'mask_ice'   ) THEN
    ELSEIF (field_name == 'mask_sheet' ) THEN
    ELSEIF (field_name == 'mask_shelf' ) THEN
    ELSEIF (field_name == 'mask_coast' ) THEN
    ELSEIF (field_name == 'mask_margin') THEN
    ELSEIF (field_name == 'mask_gl'    ) THEN
    ELSEIF (field_name == 'mask_cf'    ) THEN

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%phi_fric_a( region%mesh%vi1:region%mesh%vi2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'tau_yield') THEN
      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%tauc_a( region%mesh%vi1:region%mesh%vi2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'beta_sq') THEN
      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%beta_sq_a( region%mesh%vi1:region%mesh%vi2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%IsoIce, id_var, netcdf%ti)
    ELSEIF (field_name == 'iso_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%IsoSurf, id_var, netcdf%ti)

    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%Hb_a( region%mesh%vi1:region%mesh%vi2) - region%refgeo_PD%Hb( region%mesh%vi1:region%mesh%vi2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dHb_dt_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'dSL_dt') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%dSL_dt_a, id_var, netcdf%ti)
    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wdp_2D_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_help_field_grid

  SUBROUTINE create_restart_file_grid( region, netcdf)
    ! Create a new restart NetCDF file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_restart),      INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_restart_file_grid'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x,y,z,m,t

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! If the file already exists, return.
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Create netCDF file

    ! WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error(nf90_create(netcdf%filename,IOR(nf90_clobber,nf90_share),netcdf%ncid))

    ! Set time frame index to 1
    netcdf%ti = 1

    ! Define dimensions:
    CALL create_dim( netcdf%ncid, netcdf%name_dim_x,         region%grid_output%nx, netcdf%id_dim_x    )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_y,         region%grid_output%ny, netcdf%id_dim_y    )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_zeta,      C%nZ,                  netcdf%id_dim_zeta ) ! Scaled vertical coordinate
    CALL create_dim( netcdf%ncid, netcdf%name_dim_month,     12,                    netcdf%id_dim_month) ! Months (for monthly data)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time,      nf90_unlimited,        netcdf%id_dim_time ) ! Time frames

    ! Placeholders for the dimension ID's, for shorter code
    x = netcdf%id_dim_x
    y = netcdf%id_dim_y
    z = netcdf%id_dim_zeta
    m = netcdf%id_dim_month
    t = netcdf%id_dim_time

    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.

    ! Dimension variables: zeta, month, time
    CALL create_double_var( netcdf%ncid, netcdf%name_var_x,                  [x            ], netcdf%id_var_x,                  long_name='X-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_y,                  [   y         ], netcdf%id_var_y,                  long_name='Y-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_zeta,               [      z      ], netcdf%id_var_zeta,               long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_month,              [         m   ], netcdf%id_var_month,              long_name='Month', units='1-12'    )
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time,               [            t], netcdf%id_var_time,               long_name='Time', units='years'   )

    ! Ice model data

    ! Geometry
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hi,                 [x, y,       t], netcdf%id_var_Hi,                 long_name='Ice thickness', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hb,                 [x, y,       t], netcdf%id_var_Hb,                 long_name='Bedrock elevation', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hs,                 [x, y,       t], netcdf%id_var_Hs,                 long_name='Surface elevation', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_SL,                 [x, y,       t], netcdf%id_var_SL,                 long_name='Sea surface change', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_dHb,                [x, y,       t], netcdf%id_var_dHb,                long_name='Bedrock deformation', units='m')

    ! Bed roughness
    CALL create_double_var( netcdf%ncid, netcdf%name_var_beta_sq,            [x, y,       t], netcdf%id_var_beta_sq,            long_name='Bed roughness', units='?')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_phi_fric,           [x, y,       t], netcdf%id_var_phi_fric,           long_name='Bed roughness', units='?')

    ! Temperature
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Ti,                 [x, y, z,    t], netcdf%id_var_Ti,                 long_name='Ice temperature', units='K')

    ! SMB
    IF (C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_FirnDepth,        [x, y,    m, t], netcdf%id_var_FirnDepth,          long_name='Firn depth', units='m')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_MeltPreviousYear, [x, y,       t], netcdf%id_var_MeltPreviousYear,   long_name='Melt during previous year', units='mie')
      IF (C%do_SMB_IMAUITM_inversion) THEN
        CALL create_double_var( netcdf%ncid, netcdf%name_var_C_abl_constant_inv,   [x, y, t], netcdf%id_var_C_abl_constant_inv, long_name='Threshold ablation factor', units='-')
        CALL create_double_var( netcdf%ncid, netcdf%name_var_C_abl_Ts_inv,         [x, y, t], netcdf%id_var_C_abl_Ts_inv,       long_name='Temperature ablation factor', units='-')
        CALL create_double_var( netcdf%ncid, netcdf%name_var_C_abl_Q_inv,          [x, y, t], netcdf%id_var_C_abl_Q_inv,        long_name='Insolation ablation factor', units='-')
        CALL create_double_var( netcdf%ncid, netcdf%name_var_C_refr_inv,           [x, y, t], netcdf%id_var_C_refr_inv,         long_name='Refreezing factor', units='-')
      END IF
    END IF

    ! Leave definition mode
    CALL handle_error(nf90_enddef( netcdf%ncid))

    ! Write zeta and month dimension variables
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_x,        region%grid_output%x                     ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_y,        region%grid_output%y                     ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_grid

  SUBROUTINE create_help_fields_file_grid( region, netcdf)
    ! Create a new help fields file, containing secondary model output (not needed for a restart, but interesting to look at)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_fields_file_grid'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x, y, z, m, t

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! If the file already exists, return.
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error(nf90_create(netcdf%filename,IOR(nf90_clobber,nf90_share),netcdf%ncid))

    ! Set time frame index to 1
    netcdf%ti = 1

    ! Define dimensions:
    CALL create_dim( netcdf%ncid, netcdf%name_dim_x,     region%grid_output%nx, netcdf%id_dim_x    )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_y,     region%grid_output%ny, netcdf%id_dim_y    )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_zeta,  C%nZ,                  netcdf%id_dim_zeta )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_month, 12,                    netcdf%id_dim_month)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time,  nf90_unlimited,        netcdf%id_dim_time )

    ! Placeholders for the dimension ID's, for shorter code
    x = netcdf%id_dim_x
    y = netcdf%id_dim_y
    z = netcdf%id_dim_zeta
    m = netcdf%id_dim_month
    t = netcdf%id_dim_time

    ! Dimension variables: zeta, month, time
    CALL create_double_var( netcdf%ncid, netcdf%name_var_x,     [netcdf%id_dim_x    ], netcdf%id_var_x,     long_name='X-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_y,     [netcdf%id_dim_y    ], netcdf%id_var_y,     long_name='Y-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_zeta,  [netcdf%id_dim_zeta ], netcdf%id_var_zeta,  long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_month, [netcdf%id_dim_month], netcdf%id_var_month, long_name='Month', units='1-12'    )
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time,  [netcdf%id_dim_time ], netcdf%id_var_time,  long_name='Time', units='years'   )

    ! Define data variables
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_01, C%help_field_01)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_02, C%help_field_02)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_03, C%help_field_03)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_04, C%help_field_04)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_05, C%help_field_05)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_06, C%help_field_06)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_07, C%help_field_07)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_08, C%help_field_08)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_09, C%help_field_09)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_10, C%help_field_10)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_11, C%help_field_11)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_12, C%help_field_12)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_13, C%help_field_13)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_14, C%help_field_14)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_15, C%help_field_15)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_16, C%help_field_16)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_17, C%help_field_17)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_18, C%help_field_18)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_19, C%help_field_19)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_20, C%help_field_20)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_21, C%help_field_21)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_22, C%help_field_22)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_23, C%help_field_23)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_24, C%help_field_24)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_25, C%help_field_25)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_26, C%help_field_26)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_27, C%help_field_27)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_28, C%help_field_28)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_29, C%help_field_29)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_30, C%help_field_30)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_31, C%help_field_31)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_32, C%help_field_32)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_33, C%help_field_33)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_34, C%help_field_34)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_35, C%help_field_35)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_36, C%help_field_36)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_37, C%help_field_37)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_38, C%help_field_38)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_39, C%help_field_39)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_40, C%help_field_40)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_41, C%help_field_41)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_42, C%help_field_42)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_43, C%help_field_43)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_44, C%help_field_44)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_45, C%help_field_45)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_46, C%help_field_46)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_47, C%help_field_47)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_48, C%help_field_48)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_49, C%help_field_49)
    CALL create_help_field_grid( netcdf, netcdf%id_help_field_50, C%help_field_50)

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( netcdf%ncid))

    ! Write the x, y, zeta, months, and lat/lon variable data
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_x,        region%grid_output%x                     ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_y,        region%grid_output%y                     ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_fields_file_grid

  SUBROUTINE create_help_field_grid( netcdf, id_var, field_name)
    ! Add a data field to the help_fields file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf
    INTEGER,                        INTENT(INOUT) :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_help_field_grid'
    INTEGER                                       :: x, y, t, z, m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Placeholders for the dimension ID's, for shorter code
    x         = netcdf%id_dim_x
    y         = netcdf%id_dim_y
    t         = netcdf%id_dim_time
    z         = netcdf%id_dim_zeta
    m         = netcdf%id_dim_month

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    ELSEIF (field_name == 'lat') THEN
      CALL create_double_var( netcdf%ncid, 'lat',                      [x, y      ], id_var, long_name='Latitude',  units='degrees north')
    ELSEIF (field_name == 'lon') THEN
      CALL create_double_var( netcdf%ncid, 'lon',                      [x, y      ], id_var, long_name='Longitude', units='degrees east')

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL create_double_var( netcdf%ncid, 'GHF',                      [x, y      ], id_var, long_name='Geothermal heat flux', units='J m^-2 yr^-1')

    ! Fields with a time dimension
    ! ============================

    ! Mesh
    ELSEIF (field_name == 'resolution') THEN
      CALL create_double_var( netcdf%ncid, 'resolution',               [x, y,    t], id_var, long_name='Mesh resolution', units='m')

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL create_double_var( netcdf%ncid, 'Hi',                       [x, y,    t], id_var, long_name='Ice thickness', units='m')
    ELSEIF (field_name == 'Hb') THEN
      CALL create_double_var( netcdf%ncid, 'Hb',                       [x, y,    t], id_var, long_name='Bedrock elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'Hs') THEN
      CALL create_double_var( netcdf%ncid, 'Hs',                       [x, y,    t], id_var, long_name='Surface elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'SL') THEN
      CALL create_double_var( netcdf%ncid, 'SL',                       [x, y,    t], id_var, long_name='Geoid elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'dHi') THEN
      CALL create_double_var( netcdf%ncid, 'dHi',                      [x, y,    t], id_var, long_name='Ice thickness difference w.r.t. PD', units='m')
    ELSEIF (field_name == 'dHs') THEN
      CALL create_double_var( netcdf%ncid, 'dHs',                      [x, y,    t], id_var, long_name='Ice elevation difference w.r.t. PD', units='m')

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL create_double_var( netcdf%ncid, 'Ti',                       [x, y, z, t], id_var, long_name='Englacial temperature', units='K')
    ELSEIF (field_name == 'Cpi') THEN
      CALL create_double_var( netcdf%ncid, 'Cpi',                      [x, y, z, t], id_var, long_name='Ice heat capacity', units='J kg^-1 K^-1')
    ELSEIF (field_name == 'Ki') THEN
      CALL create_double_var( netcdf%ncid, 'Ki',                       [x, y, z, t], id_var, long_name='Ice thermal conductivity', units='J m^-1 K^-1 yr^-1')
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL create_double_var( netcdf%ncid, 'Ti_basal',                 [x, y,    t], id_var, long_name='Ice basal temperature', units='K')
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL create_double_var( netcdf%ncid, 'Ti_pmp',                   [x, y, z, t], id_var, long_name='Ice pressure melting point temperature', units='K')
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL create_double_var( netcdf%ncid, 'A_flow_3D',                [x, y, z, t], id_var, long_name='Ice flow factor', units='Pa^-3 y^-1')
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL create_double_var( netcdf%ncid, 'A_flow_vav',               [x, y,    t], id_var, long_name='Vertically averaged ice flow factor', units='Pa^-3 y^-1')

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL create_double_var( netcdf%ncid, 'u_3D',                     [x, y, z, t], id_var, long_name='3D ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_3D') THEN
      CALL create_double_var( netcdf%ncid, 'v_3D',                     [x, y, z, t], id_var, long_name='3D ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_3D_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_3D_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'w_3D') THEN
      CALL create_double_var( netcdf%ncid, 'w_3D',                     [x, y, z, t], id_var, long_name='3D ice z-velocity', units='m/yr')
    ELSEIF (field_name == 'u_vav') THEN
      CALL create_double_var( netcdf%ncid, 'u_vav',                    [x, y,    t], id_var, long_name='Vertically averaged ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_vav') THEN
      CALL create_double_var( netcdf%ncid, 'v_vav',                    [x, y,    t], id_var, long_name='Vertically averaged ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'u_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_vav',                 [x, y,    t], id_var, long_name='Vertically averaged ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_vav_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'u_surf') THEN
      CALL create_double_var( netcdf%ncid, 'u_surf',                   [x, y,    t], id_var, long_name='Surface ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_surf') THEN
      CALL create_double_var( netcdf%ncid, 'v_surf',                   [x, y,    t], id_var, long_name='Surface ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_surf',                [x, y,    t], id_var, long_name='Surface ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_surf_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'u_base') THEN
      CALL create_double_var( netcdf%ncid, 'u_base',                   [x, y,    t], id_var, long_name='Basal ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_base') THEN
      CALL create_double_var( netcdf%ncid, 'v_base',                   [x, y,    t], id_var, long_name='Basal ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'u_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'v_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file
    ELSEIF (field_name == 'uabs_base') THEN
      CALL create_double_var( netcdf%ncid, 'uabs_base',                [x, y,    t], id_var, long_name='Basal ice velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_base_b') THEN
      ! Do nothing; b-grid variables are only written to the mesh-version of the file

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL create_double_var( netcdf%ncid, 'T2m',                      [x, y, m, t], id_var, long_name='Monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'T2m_year') THEN
      CALL create_double_var( netcdf%ncid, 'T2m_year',                 [x, y,    t], id_var, long_name='Annual mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'Precip') THEN
      CALL create_double_var( netcdf%ncid, 'Precip',                   [x, y, m, t], id_var, long_name='Monthly total precipitation', units='mm')
    ELSEIF (field_name == 'Precip_year') THEN
      CALL create_double_var( netcdf%ncid, 'Precip_year',              [x, y,    t], id_var, long_name='Annual total precipitation', units='mm')
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_WE',                  [x, y, m, t], id_var, long_name='Monthly mean zonal wind', units='m/s')
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_WE_year',             [x, y,    t], id_var, long_name='Annual mean zonal wind', units='m/s')
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_SN',                  [x, y, m, t], id_var, long_name='Monthly mean meridional wind', units='m/s')
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL create_double_var( netcdf%ncid, 'Wind_SN_year',             [x, y,    t], id_var, long_name='Annual mean meridional wind', units='m/s')

    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL create_double_var( netcdf%ncid, 'SMB',                      [x, y, m, t], id_var, long_name='Monthly surface mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'SMB_year') THEN
      CALL create_double_var( netcdf%ncid, 'SMB_year',                 [x, y,    t], id_var, long_name='Annual surface mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL create_double_var( netcdf%ncid, 'BMB_sheet',                [x, y,    t], id_var, long_name='Annual basal mass balance for grounded ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL create_double_var( netcdf%ncid, 'BMB_shelf',                [x, y,    t], id_var, long_name='Annual basal mass balance for floating ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB') THEN
      CALL create_double_var( netcdf%ncid, 'BMB',                      [x, y,    t], id_var, long_name='Annual basal mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'Snowfall') THEN
      CALL create_double_var( netcdf%ncid, 'Snowfall',                 [x, y, m, t], id_var, long_name='Monthly total snowfall', units='m water equivalent')
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL create_double_var( netcdf%ncid, 'Snowfall_year',            [x, y,    t], id_var, long_name='Annual total snowfall', units='m water equivalent')
    ELSEIF (field_name == 'Rainfall') THEN
      CALL create_double_var( netcdf%ncid, 'Rainfall',                 [x, y, m, t], id_var, long_name='Monthly total rainfall', units='m water equivalent')
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL create_double_var( netcdf%ncid, 'Rainfall_year',            [x, y,    t], id_var, long_name='Annual total rainfall', units='m water equivalent')
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL create_double_var( netcdf%ncid, 'AddedFirn',                [x, y, m, t], id_var, long_name='Monthly total added firn', units='m water equivalent')
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL create_double_var( netcdf%ncid, 'AddedFirn_year',           [x, y,    t], id_var, long_name='Annual total added firn', units='m water equivalent')
    ELSEIF (field_name == 'Refreezing') THEN
      CALL create_double_var( netcdf%ncid, 'Refreezing',               [x, y, m, t], id_var, long_name='Monthly total refreezing', units='m water equivalent')
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL create_double_var( netcdf%ncid, 'Refreezing_year',          [x, y,    t], id_var, long_name='Annual total refreezing', units='m water equivalent')
    ELSEIF (field_name == 'Melt') THEN
      CALL create_double_var( netcdf%ncid, 'Melt',                     [x, y, m, t], id_var, long_name='Monthly total melt', units='m water equivalent')
    ELSEIF (field_name == 'Melt_year') THEN
      CALL create_double_var( netcdf%ncid, 'Melt_year',                [x, y,    t], id_var, long_name='Annual total melt', units='m water equivalent')
    ELSEIF (field_name == 'Runoff') THEN
      CALL create_double_var( netcdf%ncid, 'Runoff',                   [x, y, m, t], id_var, long_name='Monthly total runoff', units='m water equivalent')
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL create_double_var( netcdf%ncid, 'Runoff_year',              [x, y,    t], id_var, long_name='Annual total runoff', units='m water equivalent')
    ELSEIF (field_name == 'Albedo') THEN
      CALL create_double_var( netcdf%ncid, 'Albedo',                   [x, y, m, t], id_var, long_name='Monthly mean albedo')
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL create_double_var( netcdf%ncid, 'Albedo_year',              [x, y,    t], id_var, long_name='Annual mean albedo')
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL create_double_var( netcdf%ncid, 'FirnDepth',                [x, y, m, t], id_var, long_name='Monthly mean firn layer depth', units='m water equivalent')
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL create_double_var( netcdf%ncid, 'FirnDepth_year',           [x, y,    t], id_var, long_name='Annual mean firn layer depth', units='m water equivalent')
    ELSEIF (field_name == 'C_abl_constant_inv') THEN
      CALL create_double_var( netcdf%ncid, 'C_abl_constant_inv',       [x, y,    t], id_var, long_name='Constant ablation', units='-')
    ELSEIF (field_name == 'C_abl_Ts_inv') THEN
      CALL create_double_var( netcdf%ncid, 'C_abl_Ts_inv',             [x, y,    t], id_var, long_name='Temperature ablation factor', units='-')
    ELSEIF (field_name == 'C_abl_Q_inv') THEN
      CALL create_double_var( netcdf%ncid, 'C_abl_Q_inv',              [x, y,    t], id_var, long_name='Insolation ablation factor', units='-')
    ELSEIF (field_name == 'C_refr_inv') THEN
      CALL create_double_var( netcdf%ncid, 'C_refr_inv',               [x, y,    t], id_var, long_name='Refreezing factor', units='-')

    ! NOTE: masks commented out; mapping masks between grids is meaningless

    ! Masks
    ELSEIF (field_name == 'mask') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask',                     [x, y,    t], id_var, long_name='mask')
    ELSEIF (field_name == 'mask_land') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_land',                [x, y,    t], id_var, long_name='land mask')
    ELSEIF (field_name == 'mask_ocean') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_ocean',               [x, y,    t], id_var, long_name='ocean mask')
    ELSEIF (field_name == 'mask_lake') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_lake',                [x, y,    t], id_var, long_name='lake mask')
    ELSEIF (field_name == 'mask_ice') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_ice',                 [x, y,    t], id_var, long_name='ice mask')
    ELSEIF (field_name == 'mask_sheet') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_sheet',               [x, y,    t], id_var, long_name='sheet mask')
    ELSEIF (field_name == 'mask_shelf') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_shelf',               [x, y,    t], id_var, long_name='shelf mask')
    ELSEIF (field_name == 'mask_coast') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_coast',               [x, y,    t], id_var, long_name='coast mask')
    ELSEIF (field_name == 'mask_margin') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_margin',              [x, y,    t], id_var, long_name='margin mask')
    ELSEIF (field_name == 'mask_gl') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_gl',                  [x, y,    t], id_var, long_name='grounding-line mask')
    ELSEIF (field_name == 'mask_cf') THEN
      ! CALL create_int_var(    netcdf%ncid, 'mask_cf',                  [x, y,    t], id_var, long_name='calving-front mask')

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL create_double_var( netcdf%ncid, 'phi_fric',                 [x, y,    t], id_var, long_name='till friction angle', units='degrees')
    ELSEIF (field_name == 'tau_yield') THEN
      CALL create_double_var( netcdf%ncid, 'tau_yield',                [x, y,    t], id_var, long_name='basal yield stress', units='Pa')
    ELSEIF (field_name == 'beta_sq') THEN
      CALL create_double_var( netcdf%ncid, 'beta_sq',                  [x, y,    t], id_var, long_name='sliding coefficient', units='?')

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL create_double_var( netcdf%ncid, 'iso_ice',                  [x, y,    t], id_var, long_name='Vertically averaged ice d18O', units='per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL create_double_var( netcdf%ncid, 'iso_surf',                 [x, y,    t], id_var, long_name='d18O of precipitation', units='per mille')

    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      CALL create_double_var( netcdf%ncid, 'dHb',                      [x, y,    t], id_var, long_name='Change in bedrock elevation w.r.t. PD', units='m')
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL create_double_var( netcdf%ncid, 'dHb_dt',                 [x, y,    t], id_var, long_name='Bedrock deformation rates', units='m/yr')
    ELSEIF (field_name == 'dSL_dt') THEN
      CALL create_double_var( netcdf%ncid, 'dSL_dt',                 [x, y,    t], id_var, long_name='Geoid deformation rates', units='m/yr')
    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_field_grid

  ! Map a model data field from the model mesh to the output grid, and write it to a NetCDF file.
  SUBROUTINE map_and_write_to_grid_netcdf_dp_2D(  ncid, mesh, grid, d_mesh, id_var, ti)

    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)        :: ncid
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_grid),            INTENT(IN)        :: grid
    REAL(dp), DIMENSION(:    ), INTENT(IN)        :: d_mesh
    INTEGER,                    INTENT(IN)        :: id_var
    INTEGER,                    INTENT(IN)        :: ti

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_and_write_to_grid_netcdf_dp_2D'
    REAL(dp), DIMENSION(:,:  ), POINTER           :: d_grid
    INTEGER                                       :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid, wd_grid)

    ! Map data from the model mesh to the square grid
    CALL map_mesh2grid_2D( mesh, grid, d_mesh, d_grid)

    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, ti/) ))

    ! Deallocate shared memory
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_and_write_to_grid_netcdf_dp_2D

  SUBROUTINE map_and_write_to_grid_netcdf_dp_2D_notime(  ncid, mesh, grid, d_mesh, id_var)

    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)        :: ncid
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_grid),            INTENT(IN)        :: grid
    REAL(dp), DIMENSION(:    ), INTENT(IN)        :: d_mesh
    INTEGER,                    INTENT(IN)        :: id_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_and_write_to_grid_netcdf_dp_2D_notime'
    REAL(dp), DIMENSION(:,:  ), POINTER           :: d_grid
    INTEGER                                       :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid, wd_grid)

    ! Map data from the model mesh to the square grid
    CALL map_mesh2grid_2D( mesh, grid, d_mesh, d_grid)

    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/ 1, 1/) ))

    ! Deallocate shared memory
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_and_write_to_grid_netcdf_dp_2D_notime

  SUBROUTINE map_and_write_to_grid_netcdf_dp_3D(  ncid, mesh, grid, d_mesh, id_var, ti, nz)

    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)        :: ncid
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_grid),            INTENT(IN)        :: grid
    REAL(dp), DIMENSION(:,:  ), INTENT(IN)        :: d_mesh
    INTEGER,                    INTENT(IN)        :: id_var
    INTEGER,                    INTENT(IN)        :: ti
    INTEGER,                    INTENT(IN)        :: nz

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_and_write_to_grid_netcdf_dp_3D'
    REAL(dp), DIMENSION(:,:,:), POINTER           :: d_grid
    INTEGER                                       :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, nz, d_grid, wd_grid)

    ! Map data from the model mesh to the square grid
    CALL map_mesh2grid_3D( mesh, grid, d_mesh, d_grid)

    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, 1, ti/) ))

    ! Deallocate shared memory
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_and_write_to_grid_netcdf_dp_3D

! ===== Read all kinds of input files =====
! =========================================

  ! A restart file produced by an earlier run
  SUBROUTINE inquire_restart_file_mesh( netcdf, nV, nTri, nC_mem)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf
    INTEGER,                   INTENT(OUT)   :: nV, nTri, nC_mem

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER            :: routine_name = 'inquire_restart_file_mesh'
    INTEGER                                  :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_vi,    nV,        netcdf%id_dim_vi   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_ti,    nTri,      netcdf%id_dim_ti   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_ci,    nC_mem,    netcdf%id_dim_ci   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_two,   int_dummy, netcdf%id_dim_two  )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_three, int_dummy, netcdf%id_dim_three)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_V,              (/ netcdf%id_dim_vi, netcdf%id_dim_two  /), netcdf%id_var_V             )
    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_nC,             (/ netcdf%id_dim_vi                     /), netcdf%id_var_nC            )
    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_C,              (/ netcdf%id_dim_vi, netcdf%id_dim_ci   /), netcdf%id_var_C             )
    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_niTri,          (/ netcdf%id_dim_vi                     /), netcdf%id_var_niTri         )
    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_iTri,           (/ netcdf%id_dim_vi, netcdf%id_dim_ci   /), netcdf%id_var_iTri          )
    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_edge_index,     (/ netcdf%id_dim_vi                     /), netcdf%id_var_edge_index    )
    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_Tri,            (/ netcdf%id_dim_ti, netcdf%id_dim_three/), netcdf%id_var_Tri           )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Tricc,          (/ netcdf%id_dim_ti, netcdf%id_dim_two  /), netcdf%id_var_Tricc         )
    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_TriC,           (/ netcdf%id_dim_ti, netcdf%id_dim_three/), netcdf%id_var_TriC          )
    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_Tri_edge_index, (/ netcdf%id_dim_ti                     /), netcdf%id_var_Tri_edge_index)

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_restart_file_mesh

  SUBROUTINE inquire_restart_file_init( netcdf)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER            :: routine_name = 'inquire_restart_file_init'
    INTEGER                                  :: nZ, nt, nm, k
    REAL(dp), DIMENSION(:    ), ALLOCATABLE  :: zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_zeta,  nZ, netcdf%id_dim_zeta )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_time,  nt, netcdf%id_dim_time )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_month, nm, netcdf%id_dim_month)

    IF (nZ /= C%nZ) THEN
      CALL crash('nZ in restart file doesnt match nZ in config!')
    END IF

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_zeta,             (/ netcdf%id_dim_zeta                    /), netcdf%id_var_zeta            )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_time,             (/ netcdf%id_dim_time                    /), netcdf%id_var_time            )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_month,            (/ netcdf%id_dim_month                   /), netcdf%id_var_month           )

    ! Read zeta, check if it matches the config zeta levels
    ALLOCATE( zeta( nZ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_zeta, zeta, start = (/ 1 /) ))
    DO k = 1, C%nz
      IF (ABS(C%zeta(k) - zeta(k)) > 0.0001_dp) THEN
        CALL warning('Vertical coordinate zeta in restart file doesnt match zeta in config!')
      END IF
      END DO
    DEALLOCATE( zeta)

    ! Inquire model data

    ! Geometry
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Hi,               (/ netcdf%id_dim_vi,                      netcdf%id_dim_time /), netcdf%id_var_Hi              )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Hb,               (/ netcdf%id_dim_vi,                      netcdf%id_dim_time /), netcdf%id_var_Hb              )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Hs,               (/ netcdf%id_dim_vi,                      netcdf%id_dim_time /), netcdf%id_var_Hs              )

    ! Bed roughness

    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_beta_sq,          (/ netcdf%id_dim_vi,                      netcdf%id_dim_time /), netcdf%id_var_beta_sq         )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_phi_fric,         (/ netcdf%id_dim_vi,                      netcdf%id_dim_time /), netcdf%id_var_phi_fric        )

    ! Temperature
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Ti,               (/ netcdf%id_dim_vi, netcdf%id_dim_zeta,  netcdf%id_dim_time /), netcdf%id_var_Ti              )

    ! SMB
    IF (C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL inquire_double_var( netcdf%ncid, netcdf%name_var_MeltPreviousYear, (/ netcdf%id_dim_vi,                      netcdf%id_dim_time /), netcdf%id_var_MeltPreviousYear  )
      CALL inquire_double_var( netcdf%ncid, netcdf%name_var_FirnDepth,        (/ netcdf%id_dim_vi, netcdf%id_dim_month, netcdf%id_dim_time /), netcdf%id_var_FirnDepth         )
      IF (C%do_SMB_IMAUITM_inversion .AND. C%SMB_IMAUITM_inv_choice_init_C == 'restart') THEN
        CALL inquire_double_var( netcdf%ncid, netcdf%name_var_C_abl_constant_inv, (/ netcdf%id_dim_vi,                  netcdf%id_dim_time /), netcdf%id_var_C_abl_constant_inv)
        CALL inquire_double_var( netcdf%ncid, netcdf%name_var_C_abl_Ts_inv,       (/ netcdf%id_dim_vi,                  netcdf%id_dim_time /), netcdf%id_var_C_abl_Ts_inv      )
        CALL inquire_double_var( netcdf%ncid, netcdf%name_var_C_abl_Q_inv,        (/ netcdf%id_dim_vi,                  netcdf%id_dim_time /), netcdf%id_var_C_abl_Q_inv       )
        CALL inquire_double_var( netcdf%ncid, netcdf%name_var_C_refr_inv,         (/ netcdf%id_dim_vi,                  netcdf%id_dim_time /), netcdf%id_var_C_refr_inv        )
      END IF
    END IF

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_restart_file_init

  SUBROUTINE read_restart_file_mesh( mesh, netcdf)
    ! Read mesh data from a restart file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),           INTENT(INOUT) :: mesh
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER            :: routine_name = 'read_restart_file_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_V,              mesh%V,              start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_nC,             mesh%nC,             start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_C,              mesh%C,              start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_niTri,          mesh%niTri,          start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_iTri,           mesh%iTri,           start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_edge_index,     mesh%edge_index,     start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Tri,            mesh%Tri,            start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Tricc,          mesh%Tricc,          start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_TriC,           mesh%TriC,           start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Tri_edge_index, mesh%Tri_edge_index, start = (/ 1    /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_restart_file_mesh

  SUBROUTINE read_restart_file_init( name, restart, netcdf)
    ! Read mesh data from a restart file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),              INTENT(IN)    :: name
    TYPE(type_restart_data),       INTENT(INOUT) :: restart
    TYPE(type_netcdf_restart),     INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'read_restart_file_init'
    INTEGER                                      :: nt, ti, ti_min
    REAL(dp), DIMENSION(:    ), ALLOCATABLE      :: time
    REAL(dp)                                     :: dt_min, dt, time_to_restart_from

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set the restarting time
   IF     (name == 'NAM') THEN
     time_to_restart_from = C%time_to_restart_from_NAM
   ELSEIF (name == 'EAS') THEN
     time_to_restart_from = C%time_to_restart_from_EAS
   ELSEIF (name == 'GRL') THEN
     time_to_restart_from = C%time_to_restart_from_GRL
   ELSEIF (name == 'ANT') THEN
     time_to_restart_from = C%time_to_restart_from_ANT
   END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Read time, determine which time frame to read
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_time,  nt, netcdf%id_dim_time )

    ALLOCATE( time( nt))

    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_time, time, start = (/ 1 /) ))

    IF (time_to_restart_from < MINVAL(time) .OR. time_to_restart_from > MAXVAL(time)) THEN
      CALL crash('time_to_restart_from outside range of restart file!')
    END IF

    ti_min = 0
    dt_min = 1E8_dp
    DO ti = 1, nt
      dt = ABS(time( ti) - time_to_restart_from)
      IF (dt < dt_min) THEN
        ti_min = ti
        dt_min = dt
      END IF
    END DO
    ti = ti_min

    DEALLOCATE( time)

    ! Read the data

    ! Geometry
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Hi,               restart%Hi,               start = (/ 1,    ti /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Hb,               restart%Hb,               start = (/ 1,    ti /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Hs,               restart%Hs,               start = (/ 1,    ti /) ))

    ! Bed roughness
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_beta_sq,          restart%beta_sq,          start = (/ 1,    ti /) ))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_phi_fric,         restart%phi_fric,         start = (/ 1,    ti /) ))

    ! Temperature
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Ti,               restart%Ti,               start = (/ 1, 1, ti /) ))

    ! SMB
    IF (C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_MeltPreviousYear, restart%MeltPreviousYear, start = (/ 1,    ti /) ))
      CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_FirnDepth,        restart%FirnDepth,        start = (/ 1, 1, ti /) ))
      IF (C%do_SMB_IMAUITM_inversion .AND. C%SMB_IMAUITM_inv_choice_init_C == 'restart') THEN
        CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_C_abl_constant_inv, restart%C_abl_constant_inv, start = (/ 1,    ti /) ))
        CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_C_abl_Ts_inv,       restart%C_abl_Ts_inv,       start = (/ 1,    ti /) ))
        CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_C_abl_Q_inv,        restart%C_abl_Q_inv,        start = (/ 1,    ti /) ))
        CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_C_refr_inv,         restart%C_refr_inv,         start = (/ 1,    ti /) ))
      END IF
    END IF

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_restart_file_init

  ! Reference ice-sheet geometry (ice thickness, bed topography, and surface elevation)
  SUBROUTINE inquire_reference_geometry_file( refgeo)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_reference_geometry), INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_reference_geometry_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the netcdf file
    CALL open_netcdf_file( refgeo%netcdf%filename, refgeo%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( refgeo%netcdf%ncid, refgeo%netcdf%name_dim_x, refgeo%grid%nx, refgeo%netcdf%id_dim_x)
    CALL inquire_dim( refgeo%netcdf%ncid, refgeo%netcdf%name_dim_y, refgeo%grid%ny, refgeo%netcdf%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_x,  (/ refgeo%netcdf%id_dim_x                         /), refgeo%netcdf%id_var_x )
    CALL inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_y,  (/                         refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_y )

    CALL inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_Hi, (/ refgeo%netcdf%id_dim_x, refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_Hi)
    CALL inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_Hb, (/ refgeo%netcdf%id_dim_x, refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_Hb)
    CALL inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_Hs, (/ refgeo%netcdf%id_dim_x, refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_Hs)

    ! Close the netcdf file
    CALL close_netcdf_file( refgeo%netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_reference_geometry_file

  SUBROUTINE read_reference_geometry_file(    refgeo)
    ! Read reference geometry data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_reference_geometry), INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_reference_geometry_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the netcdf file
    CALL open_netcdf_file( refgeo%netcdf%filename, refgeo%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_x,      refgeo%grid%x,  start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_y,      refgeo%grid%y,  start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_Hi,     refgeo%Hi_grid, start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_Hb,     refgeo%Hb_grid, start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_Hs,     refgeo%Hs_grid, start = (/ 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( refgeo%netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_reference_geometry_file

  ! Insolation solution (e.g. Laskar 2004)
  SUBROUTINE inquire_insolation_data_file( forcing)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_insolation_data_file'
    INTEGER                                :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_time,     forcing%ins_nyears,        forcing%netcdf_ins%id_dim_time)
    CALL inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_month,    int_dummy,                 forcing%netcdf_ins%id_dim_month)
    CALL inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_lat,      forcing%ins_nlat,          forcing%netcdf_ins%id_dim_lat)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_time,  (/ forcing%netcdf_ins%id_dim_time                                                                 /), forcing%netcdf_ins%id_var_time)
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_month, (/ forcing%netcdf_ins%id_dim_month                                                                /), forcing%netcdf_ins%id_var_month)
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_lat,   (/ forcing%netcdf_ins%id_dim_lat                                                                  /), forcing%netcdf_ins%id_var_lat)
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_Q_TOA, (/ forcing%netcdf_ins%id_dim_time, forcing%netcdf_ins%id_dim_month, forcing%netcdf_ins%id_dim_lat /), forcing%netcdf_ins%id_var_Q_TOA)

    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ins%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_insolation_data_file

  SUBROUTINE read_insolation_data_file( forcing, ti0, ti1, ins_Q_TOA0, ins_Q_TOA1)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    INTEGER,                        INTENT(IN)    :: ti0, ti1
    REAL(dp), DIMENSION(:,:),       INTENT(OUT)   :: ins_Q_TOA0, ins_Q_TOA1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_insolation_data_file'
    INTEGER                                       :: mi, li
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Q_temp0, Q_temp1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE( Q_temp0(1, 12, forcing%ins_nlat))
    ALLOCATE( Q_temp1(1, 12, forcing%ins_nlat))

    ! Read data
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_Q_TOA, Q_temp0, start = (/ ti0, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_Q_TOA, Q_temp1, start = (/ ti1, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /), stride = (/ 1, 1, 1 /) ))
    CALL close_netcdf_file(forcing%netcdf_ins%ncid)

    ! Store the data in the shared memory structure
    DO mi = 1, 12
    DO li = 1, forcing%ins_nlat
      ins_Q_TOA0( li,mi) = Q_temp0( 1,mi,li)
      ins_Q_TOA1( li,mi) = Q_temp1( 1,mi,li)
    END DO
    END DO

    ! Clean up temporary memory
    DEALLOCATE(Q_temp0)
    DEALLOCATE(Q_temp1)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_insolation_data_file

  SUBROUTINE read_insolation_data_file_time_lat( forcing)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_insolation_data_file_time_lat'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_time,    forcing%ins_time,    start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_lat,     forcing%ins_lat,     start = (/ 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ins%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_insolation_data_file_time_lat

  ! Geothermal heat flux
  SUBROUTINE inquire_geothermal_heat_flux_file( forcing)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_geothermal_heat_flux_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ghf%filename, forcing%netcdf_ghf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_dim_lon, forcing%grid_ghf%nlon, forcing%netcdf_ghf%id_dim_lon)
    CALL inquire_dim( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_dim_lat, forcing%grid_ghf%nlat, forcing%netcdf_ghf%id_dim_lat)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_lon, (/ forcing%netcdf_ghf%id_dim_lon                                /), forcing%netcdf_ghf%id_var_lon)
    CALL inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_lat, (/ forcing%netcdf_ghf%id_dim_lat                                /), forcing%netcdf_ghf%id_var_lat)
    CALL inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_ghf, (/ forcing%netcdf_ghf%id_dim_lon, forcing%netcdf_ghf%id_dim_lat /), forcing%netcdf_ghf%id_var_ghf)

    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ghf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_geothermal_heat_flux_file

  SUBROUTINE read_geothermal_heat_flux_file( forcing)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_geothermal_heat_flux_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ghf%filename, forcing%netcdf_ghf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_lon, forcing%grid_ghf%lon, start=(/1   /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_lat, forcing%grid_ghf%lat, start=(/1   /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_ghf, forcing%ghf_ghf,      start=(/1, 1/) ))

    ! Close the NetCDF file
    CALL close_netcdf_file(forcing%netcdf_ghf%ncid)

    ! Convert from W m-2 (J m-2 s-1) to J m-2 yr-1
    forcing%ghf_ghf = forcing%ghf_ghf * sec_per_year

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_geothermal_heat_flux_file

  ! Write a sparse matrix to a NetCDF file
  SUBROUTINE write_PETSc_matrix_to_NetCDF( A, filename)
    ! Write a PETSc matrix to a NetCDF file

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(tMat),                          INTENT(IN)    :: A
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_PETSc_matrix_to_NetCDF'
    TYPE(type_sparse_matrix_CSR_dp)                    :: A_CSR

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get matrix in CSR format using native Fortran arrays
    CALL mat_petsc2CSR( A, A_CSR)

    ! Write the CSR matrix to a file
    CALL write_CSR_matrix_to_NetCDF( A_CSR, filename)

    ! Clean up after yourself
    CALL deallocate_matrix_CSR( A_CSR)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_PETSc_matrix_to_NetCDF

  SUBROUTINE write_CSR_matrix_to_NetCDF( A_CSR, filename)
    ! Write a CSR matrix to a NetCDF file

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR_dp),     INTENT(IN)    :: A_CSR
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_CSR_matrix_to_NetCDF'
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_m, id_dim_mp1, id_dim_n, id_dim_nnz
    INTEGER                                            :: id_var_ptr, id_var_index, id_var_val

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      ! Safety
      INQUIRE(EXIST=file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL crash('file "' // TRIM( filename) // '" already exists!')
      END IF

      WRITE(0,*) '   NOTE: writing CSR matrix to file "', TRIM(filename), '"'

      ! Create netCDF file
      CALL handle_error( nf90_create( TRIM(C%output_dir)//TRIM(filename), IOR(nf90_clobber,nf90_share), ncid))

      ! Define dimensions:
      CALL create_dim( ncid, 'm',      A_CSR%m  , id_dim_m  )
      CALL create_dim( ncid, 'mplus1', A_CSR%m+1, id_dim_mp1)
      CALL create_dim( ncid, 'n',      A_CSR%n  , id_dim_n  )
      CALL create_dim( ncid, 'nnz',    A_CSR%nnz, id_dim_nnz)

      ! Define variables:
      ! The order of the CALL statements for the different variables determines their
      ! order of appearence in the netcdf file.

      CALL create_int_var(    ncid, 'ptr',   [id_dim_mp1], id_var_ptr  , long_name = 'ptr'  )
      CALL create_int_var(    ncid, 'index', [id_dim_nnz], id_var_index, long_name = 'index')
      CALL create_double_var( ncid, 'val',   [id_dim_nnz], id_var_val  , long_name = 'val'  )

      ! Leave definition mode
      CALL handle_error( nf90_enddef( ncid))

      ! Write data
      CALL handle_error( nf90_put_var( ncid, id_var_ptr  , A_CSR%ptr  ))
      CALL handle_error( nf90_put_var( ncid, id_var_index, A_CSR%index))
      CALL handle_error( nf90_put_var( ncid, id_var_val  , A_CSR%val  ))

      ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
      CALL handle_error(nf90_sync( ncid))

      ! Close the file
      CALL close_netcdf_file( ncid)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_CSR_matrix_to_NetCDF

  ! Present-day observed global climate (e.g. ERA-40)
  SUBROUTINE inquire_PD_obs_global_climate_file( clim)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_climate_snapshot_global), INTENT(INOUT) :: clim

    ! Local variables:
    INTEGER                               :: int_dummy

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lat,     clim%nlat,  clim%netcdf%id_dim_lat)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lon,     clim%nlon,  clim%netcdf%id_dim_lon)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_month,   int_dummy,  clim%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lat,      (/ clim%netcdf%id_dim_lat                                                   /), clim%netcdf%id_var_lat)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lon,      (/ clim%netcdf%id_dim_lon                                                   /), clim%netcdf%id_var_lon)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Hs,       (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat                           /), clim%netcdf%id_var_Hs)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m,      (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_T2m)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Precip,   (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_Precip)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Wind_WE,  (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_Wind_WE)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Wind_SN,  (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_Wind_SN)

    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

  END SUBROUTINE inquire_PD_obs_global_climate_file

  SUBROUTINE read_PD_obs_global_climate_file(    clim)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_climate_snapshot_global), INTENT(INOUT) :: clim

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lon,     clim%lon,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lat,     clim%lat,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Hs,      clim%Hs,      start = (/ 1, 1    /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,     clim%T2m,     start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip,  clim%Precip,  start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Wind_WE, clim%Wind_WE, start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Wind_SN, clim%Wind_SN, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

  END SUBROUTINE read_PD_obs_global_climate_file

  ! GCM global climate (climate matrix snapshots)
  SUBROUTINE inquire_GCM_global_climate_file( clim)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_climate_snapshot_global), INTENT(INOUT) :: clim

    ! Local variables:
    INTEGER                                     :: int_dummy

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lat,     clim%nlat,  clim%netcdf%id_dim_lat)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lon,     clim%nlon,  clim%netcdf%id_dim_lon)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_month,   int_dummy,  clim%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lat,      (/ clim%netcdf%id_dim_lat                                                   /),  clim%netcdf%id_var_lat)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lon,      (/ clim%netcdf%id_dim_lon                                                   /),  clim%netcdf%id_var_lon)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Hs,       (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat                           /),  clim%netcdf%id_var_Hs)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m,      (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /),  clim%netcdf%id_var_T2m)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Precip,   (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /),  clim%netcdf%id_var_Precip)
    !CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Wind_WE,  (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /),  clim%netcdf%id_var_Wind_WE)
    !CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Wind_SN,  (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /),  clim%netcdf%id_var_Wind_SN)

    ! Close the netcdf file
    CALL close_netcdf_file(clim%netcdf%ncid)

  END SUBROUTINE inquire_GCM_global_climate_file

  SUBROUTINE read_GCM_global_climate_file(    clim)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_climate_snapshot_global), INTENT(INOUT) :: clim

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file(clim%netcdf%filename, clim%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lon,     clim%lon,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lat,     clim%lat,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Hs,      clim%Hs,      start = (/ 1, 1    /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,     clim%T2m,     start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip,  clim%Precip,  start = (/ 1, 1, 1 /) ))
    !CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Wind_WE, clim%Wind_WE, start = (/ 1, 1, 1 /) ))
    !CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Wind_SN, clim%Wind_SN, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file(clim%netcdf%ncid)

  END SUBROUTINE read_GCM_global_climate_file

  ! Insolation solution (e.g. Laskar 2004)
  SUBROUTINE inquire_insolation_file( forcing)
    IMPLICIT NONE

    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing

    ! Local variables:
    INTEGER                                :: int_dummy

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_time,     forcing%ins_nyears,        forcing%netcdf_ins%id_dim_time)
    CALL inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_month,    int_dummy,                 forcing%netcdf_ins%id_dim_month)
    CALL inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_lat,      forcing%ins_nlat,          forcing%netcdf_ins%id_dim_lat)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_time,  (/ forcing%netcdf_ins%id_dim_time                                                                 /), forcing%netcdf_ins%id_var_time)
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_month, (/ forcing%netcdf_ins%id_dim_month                                                                /), forcing%netcdf_ins%id_var_month)
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_lat,   (/ forcing%netcdf_ins%id_dim_lat                                                                  /), forcing%netcdf_ins%id_var_lat)
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_Q_TOA, (/ forcing%netcdf_ins%id_dim_time, forcing%netcdf_ins%id_dim_month, forcing%netcdf_ins%id_dim_lat /), forcing%netcdf_ins%id_var_Q_TOA)

    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ins%ncid)

  END SUBROUTINE inquire_insolation_file

  SUBROUTINE read_insolation_file_timeframes( forcing, ti0, ti1)
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    INTEGER,                        INTENT(IN)    :: ti0, ti1

    ! Local variables:
    INTEGER                                       :: mi, li
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Q_temp0, Q_temp1

    IF (.NOT. par%master) RETURN

    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE( Q_temp0(1, 12, forcing%ins_nlat))
    ALLOCATE( Q_temp1(1, 12, forcing%ins_nlat))

    ! Read data
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_Q_TOA, Q_temp0, start = (/ ti0, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_Q_TOA, Q_temp1, start = (/ ti1, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /), stride = (/ 1, 1, 1 /) ))
    CALL close_netcdf_file(forcing%netcdf_ins%ncid)

    ! Store the data in the shared memory structure
    DO mi = 1, 12
    DO li = 1, forcing%ins_nlat
      forcing%ins_Q_TOA0( li,mi) = Q_temp0( 1,mi,li)
      forcing%ins_Q_TOA1( li,mi) = Q_temp1( 1,mi,li)
    END DO
    END DO

    ! Clean up temporary memory
    DEALLOCATE(Q_temp0)
    DEALLOCATE(Q_temp1)

  END SUBROUTINE read_insolation_file_timeframes

  SUBROUTINE read_insolation_file_time_lat( forcing)
    IMPLICIT NONE

    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_time,    forcing%ins_time,    start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_lat,     forcing%ins_lat,     start = (/ 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ins%ncid)

  END SUBROUTINE read_insolation_file_time_lat

  ! Present-day observed global ocean (e.g. WOA18)
  SUBROUTINE inquire_PD_obs_global_ocean_file( ocn)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_ocean_snapshot_global), INTENT(INOUT) :: ocn

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( ocn%netcdf%filename, ocn%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_lat,     ocn%nlat,         ocn%netcdf%id_dim_lat    )
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_lon,     ocn%nlon,         ocn%netcdf%id_dim_lon    )
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_z_ocean, ocn%nz_ocean_raw, ocn%netcdf%id_dim_z_ocean)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_lat,     (/ ocn%netcdf%id_dim_lat     /), ocn%netcdf%id_var_lat    )
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_lon,     (/ ocn%netcdf%id_dim_lon     /), ocn%netcdf%id_var_lon    )
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_z_ocean, (/ ocn%netcdf%id_dim_z_ocean /), ocn%netcdf%id_var_z_ocean)

    CALL inquire_double_var( ocn%netcdf%ncid, TRIM(C%name_ocean_temperature), (/ ocn%netcdf%id_dim_lon, ocn%netcdf%id_dim_lat, ocn%netcdf%id_dim_z_ocean /),  ocn%netcdf%id_var_T_ocean)
    CALL inquire_double_var( ocn%netcdf%ncid, TRIM(C%name_ocean_salinity)   , (/ ocn%netcdf%id_dim_lon, ocn%netcdf%id_dim_lat, ocn%netcdf%id_dim_z_ocean /),  ocn%netcdf%id_var_S_ocean)

    ! Close the netcdf file
    CALL close_netcdf_file( ocn%netcdf%ncid)

  END SUBROUTINE inquire_PD_obs_global_ocean_file

  SUBROUTINE read_PD_obs_global_ocean_file(    ocn)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_ocean_snapshot_global), INTENT(INOUT) :: ocn

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( ocn%netcdf%filename, ocn%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_lon,     ocn%lon,         start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_lat,     ocn%lat,         start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_z_ocean, ocn%z_ocean_raw, start = (/ 1       /) ))

    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_T_ocean, ocn%T_ocean_raw, start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_S_ocean, ocn%S_ocean_raw, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( ocn%netcdf%ncid)

  END SUBROUTINE read_PD_obs_global_ocean_file

  ! GCM global ocean (ocean matrix snapshots)
  SUBROUTINE inquire_GCM_global_ocean_file( ocn)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_ocean_snapshot_global), INTENT(INOUT) :: ocn

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( ocn%netcdf%filename, ocn%netcdf%ncid)

     ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_lat,     ocn%nlat,         ocn%netcdf%id_dim_lat    )
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_lon,     ocn%nlon,         ocn%netcdf%id_dim_lon    )
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_z_ocean, ocn%nz_ocean_raw, ocn%netcdf%id_dim_z_ocean)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_lat,     (/ ocn%netcdf%id_dim_lat     /), ocn%netcdf%id_var_lat    )
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_lon,     (/ ocn%netcdf%id_dim_lon     /), ocn%netcdf%id_var_lon    )
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_z_ocean, (/ ocn%netcdf%id_dim_z_ocean /), ocn%netcdf%id_var_z_ocean)

    CALL inquire_double_var( ocn%netcdf%ncid, TRIM(C%name_ocean_temperature), (/ ocn%netcdf%id_dim_lon, ocn%netcdf%id_dim_lat, ocn%netcdf%id_dim_z_ocean /), ocn%netcdf%id_var_T_ocean)
    CALL inquire_double_var( ocn%netcdf%ncid, TRIM(C%name_ocean_salinity)   , (/ ocn%netcdf%id_dim_lon, ocn%netcdf%id_dim_lat, ocn%netcdf%id_dim_z_ocean /), ocn%netcdf%id_var_S_ocean)

    ! Close the netcdf file
    CALL close_netcdf_file( ocn%netcdf%ncid)

  END SUBROUTINE inquire_GCM_global_ocean_file

  SUBROUTINE read_GCM_global_ocean_file(    ocn)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_ocean_snapshot_global), INTENT(INOUT) :: ocn

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( ocn%netcdf%filename, ocn%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_lon,     ocn%lon,         start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_lat,     ocn%lat,         start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_z_ocean, ocn%z_ocean_raw, start = (/ 1       /) ))

    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_T_ocean, ocn%T_ocean_raw, start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_S_ocean, ocn%S_ocean_raw, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( ocn%netcdf%ncid)

  END SUBROUTINE read_GCM_global_ocean_file

  ! High-resolution geometry used for extrapolating ocean data
  SUBROUTINE inquire_hires_geometry_file( hires)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_highres_ocean_data), INTENT(INOUT) :: hires

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( hires%netcdf_geo%filename, hires%netcdf_geo%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( hires%netcdf_geo%ncid, hires%netcdf_geo%name_dim_x, hires%grid%nx, hires%netcdf_geo%id_dim_x)
    CALL inquire_dim( hires%netcdf_geo%ncid, hires%netcdf_geo%name_dim_y, hires%grid%ny, hires%netcdf_geo%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_x,  (/ hires%netcdf_geo%id_dim_x                        /), hires%netcdf_geo%id_var_x )
    CALL inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_y,  (/ hires%netcdf_geo%id_dim_y                        /), hires%netcdf_geo%id_var_y )
    CALL inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_Hi, (/ hires%netcdf_geo%id_dim_x, hires%netcdf_geo%id_dim_y /), hires%netcdf_geo%id_var_Hi)
    CALL inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_Hb, (/ hires%netcdf_geo%id_dim_x, hires%netcdf_geo%id_dim_y /), hires%netcdf_geo%id_var_Hb)

    ! Close the netcdf file
    CALL close_netcdf_file( hires%netcdf_geo%ncid)

  END SUBROUTINE inquire_hires_geometry_file

  SUBROUTINE read_hires_geometry_file(    hires)
    ! Read the high-resolution geometry netcdf file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_highres_ocean_data), INTENT(INOUT) :: hires

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( hires%netcdf_geo%filename, hires%netcdf_geo%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_x,      hires%grid%x, start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_y,      hires%grid%y, start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_Hi,     hires%Hi,     start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_Hb,     hires%Hb,     start = (/ 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( hires%netcdf_geo%ncid)

  END SUBROUTINE read_hires_geometry_file

  ! Create/read an extrapolated ocean data file
  SUBROUTINE create_extrapolated_ocean_file(  hires, hires_ocean_filename)
    ! Create a new folder extrapolated ocean data file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    CHARACTER(LEN=256),                  INTENT(IN)    :: hires_ocean_filename

    ! Local variables:
    LOGICAL                                            :: file_exists
    INTEGER                                            :: x, y, z

    IF (.NOT. par%master) RETURN

    ! Create a new file and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):

    hires%netcdf%filename = hires_ocean_filename
    INQUIRE(EXIST=file_exists, FILE = TRIM( hires%netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM(hires%netcdf%filename) // '" already exists!')
    END IF

    ! Create hires%netcdf file
    ! WRITE(0,*) '    Creating new hires%netcdf file at ', TRIM( hires%netcdf%filename)
    CALL handle_error(nf90_create( hires%netcdf%filename, IOR(nf90_clobber,nf90_share), hires%netcdf%ncid))

    ! Define dimensions:
    CALL create_dim( hires%netcdf%ncid, hires%netcdf%name_dim_x,       hires%grid%nx, hires%netcdf%id_dim_x      )
    CALL create_dim( hires%netcdf%ncid, hires%netcdf%name_dim_y,       hires%grid%ny, hires%netcdf%id_dim_y      )
    CALL create_dim( hires%netcdf%ncid, hires%netcdf%name_dim_z_ocean, C%nz_ocean,    hires%netcdf%id_dim_z_ocean)

    ! Placeholders for the dimension ID's, for shorter code
    x = hires%netcdf%id_dim_x
    y = hires%netcdf%id_dim_y
    z = hires%netcdf%id_dim_z_ocean

    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the hires%netcdf file.

    ! Dimension variables
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_x,       [x      ], hires%netcdf%id_var_x,       long_name='X-coordinate', units='m')
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_y,       [   y   ], hires%netcdf%id_var_y,       long_name='Y-coordinate', units='m')
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_z_ocean, [      z], hires%netcdf%id_var_z_ocean, long_name='Depth in ocean', units='m')

    ! Extrapolated ocean data
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_T_ocean, [x, y, z], hires%netcdf%id_var_T_ocean, long_name='3-D ocean temperature', units='K')
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_S_ocean, [x, y, z], hires%netcdf%id_var_S_ocean, long_name='3-D ocean salinity', units='PSU')

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( hires%netcdf%ncid))

    ! Write the data
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_x,        hires%grid%x ))
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_y,        hires%grid%y ))
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_z_ocean,  C%z_ocean    ))

    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_T_ocean, hires%T_ocean, start=(/ 1,1,1 /)))
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_S_ocean, hires%S_ocean, start=(/ 1,1,1 /)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( hires%netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file( hires%netcdf%ncid)

  END SUBROUTINE create_extrapolated_ocean_file

  SUBROUTINE inquire_extrapolated_ocean_file( hires)
    ! Check if the right dimensions and variables are present in the file.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_highres_ocean_data), INTENT(INOUT) :: hires

    ! Local variables:
    INTEGER                                      :: int_dummy

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( hires%netcdf%filename, hires%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( hires%netcdf%ncid, hires%netcdf%name_dim_x,       hires%grid%nx,   hires%netcdf%id_dim_x      )
    CALL inquire_dim( hires%netcdf%ncid, hires%netcdf%name_dim_y,       hires%grid%ny,   hires%netcdf%id_dim_y      )
    CALL inquire_dim( hires%netcdf%ncid, hires%netcdf%name_dim_z_ocean, int_dummy,  hires%netcdf%id_dim_z_ocean)

    ! Safety
    IF (int_dummy /= C%nz_ocean) THEN
      WRITE(0,*) 'inquire_extrapolated_ocean_file - ERROR: nz_ocean in file "', TRIM( hires%netcdf%filename), '" doesnt match ice model settings!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_x,       (/ hires%netcdf%id_dim_x                                                     /), hires%netcdf%id_var_x      )
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_y,       (/                        hires%netcdf%id_dim_y                              /), hires%netcdf%id_var_y      )
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_z_ocean, (/                                               hires%netcdf%id_dim_z_ocean /), hires%netcdf%id_var_z_ocean)

    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_T_ocean, (/ hires%netcdf%id_dim_x, hires%netcdf%id_dim_y, hires%netcdf%id_dim_z_ocean /), hires%netcdf%id_var_T_ocean)
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_S_ocean, (/ hires%netcdf%id_dim_x, hires%netcdf%id_dim_y, hires%netcdf%id_dim_z_ocean /), hires%netcdf%id_var_S_ocean)

    ! Close the netcdf file
    CALL close_netcdf_file( hires%netcdf%ncid)

  END SUBROUTINE inquire_extrapolated_ocean_file

  SUBROUTINE read_extrapolated_ocean_file(    hires)
    ! Read the extrapolated ocean data netcdf file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_highres_ocean_data), INTENT(INOUT) :: hires

    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file( hires%netcdf%filename, hires%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_x,       hires%grid%x,  start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_y,       hires%grid%y,  start = (/ 1       /) ))

    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_T_ocean, hires%T_ocean, start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_S_ocean, hires%S_ocean, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( hires%netcdf%ncid)

  END SUBROUTINE read_extrapolated_ocean_file

  ! Direct global SMB forcing
  SUBROUTINE inquire_direct_global_SMB_forcing_file( clim)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_direct_SMB_forcing_global), INTENT(INOUT) :: clim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'inquire_direct_global_SMB_forcing_file'
    INTEGER                                             :: lon,lat,t

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_global') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_global"!')
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lon,   clim%nlon,   clim%netcdf%id_dim_lon  )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lat,   clim%nlat,   clim%netcdf%id_dim_lat  )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_time,  clim%nyears, clim%netcdf%id_dim_time )

    ! Abbreviate dimension ID's for more readable code
    lon = clim%netcdf%id_dim_lon
    lat = clim%netcdf%id_dim_lat
    t   = clim%netcdf%id_dim_time

    ! Inquire variable id's. Make sure that each variable has the correct dimensions.
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lon,      (/ lon         /), clim%netcdf%id_var_lon     )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lat,      (/      lat    /), clim%netcdf%id_var_lat     )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_time,     (/           t /), clim%netcdf%id_var_time    )

    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m_year, (/ lon, lat, t /), clim%netcdf%id_var_T2m_year)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_SMB_year, (/ lon, lat, t /), clim%netcdf%id_var_SMB_year)

    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_direct_global_SMB_forcing_file

  SUBROUTINE read_direct_global_SMB_file_timeframes( clim, ti0, ti1)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_direct_SMB_forcing_global), INTENT(INOUT) :: clim
    INTEGER,                        INTENT(IN)          :: ti0, ti1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'read_direct_global_SMB_file_timeframes'
    INTEGER                                             :: loni,lati
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE           :: T2m_temp0, T2m_temp1, SMB_temp0, SMB_temp1

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_global') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_global"!')
    END IF

    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE( T2m_temp0( clim%nlon, clim%nlat, 1))
    ALLOCATE( T2m_temp1( clim%nlon, clim%nlat, 1))
    ALLOCATE( SMB_temp0( clim%nlon, clim%nlat, 1))
    ALLOCATE( SMB_temp1( clim%nlon, clim%nlat, 1))

    ! Open netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m_year, T2m_temp0, start = (/ 1, 1, ti0 /), count = (/ clim%nlon, clim%nlat, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m_year, T2m_temp1, start = (/ 1, 1, ti1 /), count = (/ clim%nlon, clim%nlat, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_SMB_year, SMB_temp0, start = (/ 1, 1, ti0 /), count = (/ clim%nlon, clim%nlat, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_SMB_year, SMB_temp1, start = (/ 1, 1, ti1 /), count = (/ clim%nlon, clim%nlat, 1 /), stride = (/ 1, 1, 1 /) ))

     ! Close netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Store the data in the shared memory structure
    DO loni = 1, clim%nlon
    DO lati = 1, clim%nlat
      clim%T2m_year0( loni,lati) = T2m_temp0( loni,lati,1)
      clim%T2m_year1( loni,lati) = T2m_temp1( loni,lati,1)
      clim%SMB_year0( loni,lati) = SMB_temp0( loni,lati,1)
      clim%SMB_year1( loni,lati) = SMB_temp1( loni,lati,1)
    END DO
    END DO

    ! Clean up after yourself
    DEALLOCATE( T2m_temp0)
    DEALLOCATE( T2m_temp1)
    DEALLOCATE( SMB_temp0)
    DEALLOCATE( SMB_temp1)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_direct_global_SMB_file_timeframes

  SUBROUTINE read_direct_global_SMB_file_time_latlon( clim)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_direct_SMB_forcing_global), INTENT(INOUT) :: clim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_global_SMB_file_time_latlon'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_global') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_global"!')
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_time, clim%time, start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lat,  clim%lat,  start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lon,  clim%lon,  start = (/ 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_direct_global_SMB_file_time_latlon

  ! Direct global climate forcing
  SUBROUTINE inquire_direct_global_climate_forcing_file( clim)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_direct_climate_forcing_global), INTENT(INOUT) :: clim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER          :: routine_name = 'inquire_direct_global_climate_forcing_file'
    INTEGER                                :: lon,lat,t,m
    INTEGER                                :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_global') THEN
      CALL crash('should only be called when choice_climate_model = "direct_global"!')
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lon,   clim%nlon,   clim%netcdf%id_dim_lon  )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lat,   clim%nlat,   clim%netcdf%id_dim_lat  )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_month, int_dummy,   clim%netcdf%id_dim_month)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_time,  clim%nyears, clim%netcdf%id_dim_time )

    ! Abbreviate dimension ID's for more readable code
    lon = clim%netcdf%id_dim_lon
    lat = clim%netcdf%id_dim_lat
    t   = clim%netcdf%id_dim_time
    m   = clim%netcdf%id_dim_month

    ! Inquire variable id's. Make sure that each variable has the correct dimensions.
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lon,    (/ lon            /), clim%netcdf%id_var_lon   )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lat,    (/      lat       /), clim%netcdf%id_var_lat   )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_month,  (/           m    /), clim%netcdf%id_var_month )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_time,   (/              t /), clim%netcdf%id_var_time  )

    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m,    (/ lon, lat, m, t /), clim%netcdf%id_var_T2m   )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Precip, (/ lon, lat, m, t /), clim%netcdf%id_var_Precip)

    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_direct_global_climate_forcing_file

  SUBROUTINE read_direct_global_climate_file_timeframes( clim, ti0, ti1)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_direct_climate_forcing_global), INTENT(INOUT) :: clim
    INTEGER,                        INTENT(IN)              :: ti0, ti1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'read_direct_global_climate_file_timeframes'
    INTEGER                                                 :: loni,lati,m
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE               :: T2m_temp0, T2m_temp1, Precip_temp0, Precip_temp1

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_global') THEN
      CALL crash('should only be called when choice_climate_model = "direct_global"!')
    END IF

    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE(    T2m_temp0( clim%nlon, clim%nlat, 12, 1))
    ALLOCATE(    T2m_temp1( clim%nlon, clim%nlat, 12, 1))
    ALLOCATE( Precip_temp0( clim%nlon, clim%nlat, 12, 1))
    ALLOCATE( Precip_temp1( clim%nlon, clim%nlat, 12, 1))

    ! Open netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,    T2m_temp0,    start = (/ 1, 1, 1, ti0 /), count = (/ clim%nlon, clim%nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,    T2m_temp1,    start = (/ 1, 1, 1, ti1 /), count = (/ clim%nlon, clim%nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip, Precip_temp0, start = (/ 1, 1, 1, ti0 /), count = (/ clim%nlon, clim%nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip, Precip_temp1, start = (/ 1, 1, 1, ti1 /), count = (/ clim%nlon, clim%nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))

     ! Close netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Store the data in the shared memory structure
    DO m    = 1, 12
    DO loni = 1, clim%nlon
    DO lati = 1, clim%nlat
      clim%T2m0(    loni,lati,m) =    T2m_temp0( loni,lati,m,1)
      clim%T2m1(    loni,lati,m) =    T2m_temp1( loni,lati,m,1)
      clim%Precip0( loni,lati,m) = Precip_temp0( loni,lati,m,1)
      clim%Precip1( loni,lati,m) = Precip_temp1( loni,lati,m,1)
    END DO
    END DO
    END DO

    ! Clean up after yourself
    DEALLOCATE(    T2m_temp0)
    DEALLOCATE(    T2m_temp1)
    DEALLOCATE( Precip_temp0)
    DEALLOCATE( Precip_temp1)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_direct_global_climate_file_timeframes

  SUBROUTINE read_direct_global_climate_file_time_latlon( clim)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_direct_climate_forcing_global), INTENT(INOUT) :: clim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_global_climate_file_time_latlon'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_global') THEN
      CALL crash('should only be called when choice_climate_model = "direct_global"!')
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_time, clim%time, start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lat,  clim%lat,  start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lon,  clim%lon,  start = (/ 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_direct_global_climate_file_time_latlon

  ! Direct regional SMB forcing
  SUBROUTINE inquire_direct_regional_SMB_forcing_file( clim)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_direct_SMB_forcing_regional), INTENT(INOUT) :: clim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'inquire_direct_regional_SMB_forcing_file'
    INTEGER                                               :: x,y,t

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_regional"!')
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_x,     clim%grid%nx, clim%netcdf%id_dim_x   )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_y,     clim%grid%ny, clim%netcdf%id_dim_y   )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_time,  clim%nyears,  clim%netcdf%id_dim_time)

    ! Abbreviate dimension ID's for more readable code
    x = clim%netcdf%id_dim_x
    y = clim%netcdf%id_dim_y
    t = clim%netcdf%id_dim_time

    ! Inquire variable id's. Make sure that each variable has the correct dimensions.
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_x,        (/ x       /), clim%netcdf%id_var_x       )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_y,        (/    y    /), clim%netcdf%id_var_y       )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_time,     (/       t /), clim%netcdf%id_var_time    )

    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m_year, (/ x, y, t /), clim%netcdf%id_var_T2m_year)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_SMB_year, (/ x, y, t /), clim%netcdf%id_var_SMB_year)

    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_direct_regional_SMB_forcing_file

  SUBROUTINE read_direct_regional_SMB_file_timeframes( clim, ti0, ti1)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_direct_SMB_forcing_regional), INTENT(INOUT) :: clim
    INTEGER,                                INTENT(IN)    :: ti0, ti1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'read_direct_regional_SMB_file_timeframes'
    INTEGER                                               :: i,j
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE             :: T2m_temp0, T2m_temp1, SMB_temp0, SMB_temp1

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_regional"!')
    END IF

    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE( T2m_temp0( clim%grid%nx, clim%grid%ny, 1))
    ALLOCATE( T2m_temp1( clim%grid%nx, clim%grid%ny, 1))
    ALLOCATE( SMB_temp0( clim%grid%nx, clim%grid%ny, 1))
    ALLOCATE( SMB_temp1( clim%grid%nx, clim%grid%ny, 1))

    ! Open netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m_year, T2m_temp0, start = (/ 1, 1, ti0 /), count = (/ clim%grid%nx, clim%grid%ny, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m_year, T2m_temp1, start = (/ 1, 1, ti1 /), count = (/ clim%grid%nx, clim%grid%nx, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_SMB_year, SMB_temp0, start = (/ 1, 1, ti0 /), count = (/ clim%grid%nx, clim%grid%ny, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_SMB_year, SMB_temp1, start = (/ 1, 1, ti1 /), count = (/ clim%grid%nx, clim%grid%nx, 1 /), stride = (/ 1, 1, 1 /) ))

     ! Close netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Store the data in the shared memory structure
    DO i = 1, clim%grid%nx
    DO j = 1, clim%grid%ny
      clim%T2m_year0_raw( j,i) = T2m_temp0( i,j,1)
      clim%T2m_year1_raw( j,i) = T2m_temp1( i,j,1)
      clim%SMB_year0_raw( j,i) = SMB_temp0( i,j,1)
      clim%SMB_year1_raw( j,i) = SMB_temp1( i,j,1)
    END DO
    END DO

    ! Clean up after yourself
    DEALLOCATE( T2m_temp0)
    DEALLOCATE( T2m_temp1)
    DEALLOCATE( SMB_temp0)
    DEALLOCATE( SMB_temp1)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_direct_regional_SMB_file_timeframes

  SUBROUTINE read_direct_regional_SMB_file_time_xy( clim)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_direct_SMB_forcing_regional), INTENT(INOUT) :: clim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'read_direct_regional_SMB_file_time_xy'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_regional"!')
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_time, clim%time,   start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_x,    clim%grid%x, start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_y,    clim%grid%y, start = (/ 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_direct_regional_SMB_file_time_xy

  ! Direct regional climate forcing
  SUBROUTINE inquire_direct_regional_climate_forcing_file( clim)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_direct_climate_forcing_regional), INTENT(INOUT) :: clim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                             :: routine_name = 'inquire_direct_regional_climate_forcing_file'
    INTEGER                                                   :: x,y,t,m
    INTEGER                                                   :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_climate_model = "direct_regional"!')
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_x,     clim%grid%nx, clim%netcdf%id_dim_x    )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_y,     clim%grid%ny, clim%netcdf%id_dim_y    )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_month, int_dummy,    clim%netcdf%id_dim_month)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_time,  clim%nyears,  clim%netcdf%id_dim_time )

    ! Abbreviate dimension ID's for more readable code
    x = clim%netcdf%id_dim_x
    y = clim%netcdf%id_dim_y
    t = clim%netcdf%id_dim_time
    m = clim%netcdf%id_dim_month

    ! Inquire variable id's. Make sure that each variable has the correct dimensions.
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_x,      (/ x          /), clim%netcdf%id_var_x     )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_y,      (/    y       /), clim%netcdf%id_var_y     )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_month,  (/       m    /), clim%netcdf%id_var_month )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_time,   (/          t /), clim%netcdf%id_var_time  )

    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m,    (/ x, y, m, t /), clim%netcdf%id_var_T2m   )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Precip, (/ x, y, m, t /), clim%netcdf%id_var_Precip)

    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_direct_regional_climate_forcing_file

  SUBROUTINE read_direct_regional_climate_file_timeframes( clim, ti0, ti1)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_direct_climate_forcing_regional), INTENT(INOUT) :: clim
    INTEGER,                        INTENT(IN)    :: ti0, ti1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_regional_climate_file_timeframes'
    INTEGER                                       :: i,j,m
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: T2m_temp0, T2m_temp1, Precip_temp0, Precip_temp1

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_climate_model = "direct_regional"!')
    END IF

    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE(    T2m_temp0( clim%grid%nx, clim%grid%ny, 12, 1))
    ALLOCATE(    T2m_temp1( clim%grid%nx, clim%grid%ny, 12, 1))
    ALLOCATE( Precip_temp0( clim%grid%nx, clim%grid%ny, 12, 1))
    ALLOCATE( Precip_temp1( clim%grid%nx, clim%grid%ny, 12, 1))

    ! Open netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,    T2m_temp0,    start = (/ 1, 1, 1, ti0 /), count = (/ clim%grid%nx, clim%grid%ny, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,    T2m_temp1,    start = (/ 1, 1, 1, ti1 /), count = (/ clim%grid%nx, clim%grid%nx, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip, Precip_temp0, start = (/ 1, 1, 1, ti0 /), count = (/ clim%grid%nx, clim%grid%nx, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip, Precip_temp1, start = (/ 1, 1, 1, ti1 /), count = (/ clim%grid%nx, clim%grid%nx, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))

     ! Close netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Store the data in the shared memory structure
    DO m = 1, 12
    DO i = 1, clim%grid%nx
    DO j = 1, clim%grid%ny
      clim%T2m0_raw(    m,j,i) =    T2m_temp0( i,j,m,1)
      clim%T2m1_raw(    m,j,i) =    T2m_temp1( i,j,m,1)
      clim%Precip0_raw( m,j,i) = Precip_temp0( i,j,m,1)
      clim%Precip1_raw( m,j,i) = Precip_temp1( i,j,m,1)
    END DO
    END DO
    END DO

    ! Clean up after yourself
    DEALLOCATE(    T2m_temp0)
    DEALLOCATE(    T2m_temp1)
    DEALLOCATE( Precip_temp0)
    DEALLOCATE( Precip_temp1)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_direct_regional_climate_file_timeframes

  SUBROUTINE read_direct_regional_climate_file_time_xy( clim)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_direct_climate_forcing_regional), INTENT(INOUT) :: clim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_regional_climate_file_time_xy'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_climate_model = "direct_regional"!')
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_time, clim%time,   start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_x,    clim%grid%x, start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_y,    clim%grid%y, start = (/ 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_direct_regional_climate_file_time_xy

! ===== SELEN =====
! =================

  ! Global topography for SELEN
  SUBROUTINE inquire_SELEN_global_topo_file( SELEN)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_SELEN_global_topo_file'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    INQUIRE(EXIST=file_exists, FILE = TRIM( SELEN%netcdf_topo%filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( SELEN%netcdf_topo%filename) // '" does not exist!')
    ELSE
      CALL open_netcdf_file( SELEN%netcdf_topo%filename, SELEN%netcdf_topo%ncid)
    END IF

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_vi,    SELEN%mesh%nV,     SELEN%netcdf_topo%id_dim_vi)
    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_ti,    SELEN%mesh%nTri,   SELEN%netcdf_topo%id_dim_ti)
    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_ci,    SELEN%mesh%nC_mem, SELEN%netcdf_topo%id_dim_ci)
    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_three, int_dummy,         SELEN%netcdf_topo%id_dim_three)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_V,     (/ SELEN%netcdf_topo%id_dim_vi, SELEN%netcdf_topo%id_dim_three /),  SELEN%netcdf_topo%id_var_V       )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_Tri,   (/ SELEN%netcdf_topo%id_dim_ti, SELEN%netcdf_topo%id_dim_three /),  SELEN%netcdf_topo%id_var_Tri     )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_nC,    (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_nC      )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_C,     (/ SELEN%netcdf_topo%id_dim_vi, SELEN%netcdf_topo%id_dim_ci    /),  SELEN%netcdf_topo%id_var_C       )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_niTri, (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_niTri   )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_iTri,  (/ SELEN%netcdf_topo%id_dim_vi, SELEN%netcdf_topo%id_dim_ci    /),  SELEN%netcdf_topo%id_var_iTri    )
    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_lat,   (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_lat     )
    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_lon,   (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_lon     )
    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_Hb,    (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_Hb      )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_ianc,  (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_ianc    )

    ! Close the netcdf file
    CALL close_netcdf_file(SELEN%netcdf_topo%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_SELEN_global_topo_file

  SUBROUTINE read_SELEN_global_topo_file( SELEN)
    ! Read the init netcdf file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_SELEN_global_topo_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file(SELEN%netcdf_topo%filename, SELEN%netcdf_topo%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_V,     SELEN%mesh%V,     start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_Tri,   SELEN%mesh%Tri,   start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_nC,    SELEN%mesh%nC,    start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_C,     SELEN%mesh%C,     start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_niTri, SELEN%mesh%niTri, start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_iTri,  SELEN%mesh%iTri,  start = (/ 1, 1 /) ))

    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_lat,   SELEN%mesh%lat,   start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_lon,   SELEN%mesh%lon,   start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_Hb,    SELEN%topo_ref,   start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_ianc,  SELEN%mesh%ianc,  start = (/ 1    /) ))

    ! Close the netcdf file
    CALL close_netcdf_file(SELEN%netcdf_topo%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_SELEN_global_topo_file

  ! SELEN output file
  SUBROUTINE create_SELEN_output_file( SELEN)
    ! Create a new NetCDF output file for SELEN (on the irregular global mesh)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_SELEN_output_file'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, three, time, ki

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set time frame index to 1
    SELEN%output%ti = 1

    ! Set output filename
    SELEN%output%filename = TRIM(C%output_dir) // 'SELEN_output.nc'

    ! Create a new output file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(SELEN%output%filename))
    IF(file_exists) THEN
      CALL crash('file "' // TRIM( SELEN%output%filename) // '" already exists!')
    END IF

    ! Create netCDF file
    CALL handle_error(nf90_create(SELEN%output%filename,IOR(nf90_clobber,nf90_share),SELEN%output%ncid))

    ! Mesh data
    ! =========

    ! Define dimensions
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_vi,           SELEN%mesh%nV,           SELEN%output%id_dim_vi          ) ! Vertex indices
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_ti,           SELEN%mesh%nTri,         SELEN%output%id_dim_ti          ) ! Triangle indices
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_ci,           SELEN%mesh%nC_mem,       SELEN%output%id_dim_ci          ) ! Connection indices
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_three,        3,                       SELEN%output%id_dim_three       ) ! 3 (each vertex has three coordinates, each triangle has three vertices)

    ! Placeholders for the dimension ID's, for shorter code
    vi        = SELEN%output%id_dim_vi
    ti        = SELEN%output%id_dim_ti
    ci        = SELEN%output%id_dim_ci
    three     = SELEN%output%id_dim_three

    ! Define variables
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_V,                [vi,  three], SELEN%output%id_var_V,                long_name='Vertex coordinates', units='m')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_Tri,              [ti,  three], SELEN%output%id_var_Tri,              long_name='Vertex indices')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_nC,               [vi        ], SELEN%output%id_var_nC,               long_name='Number of connected vertices')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_C,                [vi,  ci   ], SELEN%output%id_var_C,                long_name='Indices of connected vertices')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_niTri,            [vi        ], SELEN%output%id_var_niTri,            long_name='Number of inverse triangles')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_iTri,             [vi,  ci   ], SELEN%output%id_var_iTri,             long_name='Indices of inverse triangles')

    ! Model output
    ! ============

    ! Define dimensions
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_time,  nf90_unlimited,         SELEN%output%id_dim_time ) ! Time frames
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_ki,    C%SELEN_irreg_time_n+1, SELEN%output%id_dim_ki   ) ! Window frames

    ! Placeholders for the dimension ID's, for shorter code
    time  = SELEN%output%id_dim_time
    ki    = SELEN%output%id_dim_ki

    ! Define dimension variables
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_time,  [time  ], SELEN%output%id_var_time,  long_name='Time', units='years')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_ki,    [ki    ], SELEN%output%id_var_ki,    long_name='Window frames', units='years')

    ! Define model data variables
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_lat,              [vi             ], SELEN%output%id_var_lat,              long_name='Latitude', units='degrees north')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_lon,              [vi             ], SELEN%output%id_var_lon,              long_name='Longtitude', units='degrees east')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_Hi,               [vi,        time], SELEN%output%id_var_Hi,               long_name='Surface load', units='mie')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_Hi_rel,           [vi,        time], SELEN%output%id_var_Hi_rel,           long_name='Relative surface load', units='mie')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_U,                [vi,        time], SELEN%output%id_var_U,                long_name='Land surface change', units='m')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_N,                [vi,        time], SELEN%output%id_var_N,                long_name='Sea surface change', units='m')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_ocean_function,   [vi,        time], SELEN%output%id_var_ocean_function,   long_name='Ocean function (1 = ocean)')

    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_load_history,     [vi,   ki,  time], SELEN%output%id_var_load_history,     long_name='Load history', units='mie')

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( SELEN%output%ncid))

    ! Write mesh data
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_V,               SELEN%mesh%V             ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_Tri,             SELEN%mesh%Tri           ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_nC,              SELEN%mesh%nC            ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_C,               SELEN%mesh%C             ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_niTri,           SELEN%mesh%niTri         ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_iTri,            SELEN%mesh%iTri          ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_lat,             SELEN%mesh%lat           ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_lon,             SELEN%mesh%lon           ))

    ! Window frames
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_ki, (/0._dp, C%SELEN_irreg_time_window/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( SELEN%output%ncid))

    ! Close the file
    CALL close_netcdf_file(SELEN%output%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_SELEN_output_file

  SUBROUTINE write_to_SELEN_output_file( SELEN, time)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_to_SELEN_output_file'
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN
    REAL(dp),                       INTENT(IN)    :: time

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the file for writing
    CALL open_netcdf_file( SELEN%output%filename, SELEN%output%ncid)

    ! Time
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_time, time, start = (/ SELEN%output%ti /)))

    ! Model data
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_Hi,             SELEN%Hi_glob,                        start = (/ 1,    SELEN%output%ti/) ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_Hi_rel,         SELEN%Hi_rel_glob,                    start = (/ 1,    SELEN%output%ti/) ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_U,              SELEN%U_glob,                         start = (/ 1,    SELEN%output%ti/) ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_N,              SELEN%N_glob,                         start = (/ 1,    SELEN%output%ti/) ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_ocean_function, SELEN%of_glob,                        start = (/ 1,    SELEN%output%ti/) ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_load_history,   SELEN%ice_loading_history_irreg_glob, start = (/ 1, 1, SELEN%output%ti/) ))

    ! Close the file
    CALL close_netcdf_file(SELEN%output%ncid)

    ! Increase time frame counter
    SELEN%output%ti = SELEN%output%ti + 1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_SELEN_output_file

! ===== ISMIP-style output =====
! ==============================

  SUBROUTINE create_ISMIP_output_files( region)
    ! Create all the ISMIP output files

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region), INTENT(IN)    :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_ISMIP_output_files'
    CHARACTER(LEN=256)                                 :: icesheet_code, foldername

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Code for the ice sheet name in the ISMIP output file names
    IF     (region%name == 'NAM') THEN
      icesheet_code = 'NAIS'
    ELSEIF (region%name == 'EAS') THEN
      icesheet_code = 'EUIS'
    ELSEIF (region%name == 'GRL') THEN
      icesheet_code = 'GIS'
    ELSEIF (region%name == 'ANT') THEN
      icesheet_code = 'AIS'
    ELSEIF (region%name == 'PAT') THEN
      icesheet_code = 'PIS'
    ELSE
      icesheet_code = 'beep'
      CALL crash('unknown region "' // TRIM( region%name) // '"!')
    END IF

    ! Create a subdirectory within the output directory
    foldername = TRIM( C%output_dir) // TRIM(                icesheet_code  ) // '_' // &
                                        TRIM( C%ISMIP_output_group_code     ) // '_' // &
                                        TRIM( C%ISMIP_output_model_code     ) // '_' // &
                                        TRIM( C%ISMIP_output_experiment_code)
    CALL system('mkdir ' // foldername)

    ! Create all the ISMIP output files
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'lithk'                , 'land_ice_thickness'                                               , 'm'             )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'orog'                 , 'surface_altitude'                                                 , 'm'             )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'topg'                 , 'bedrock_altitude'                                                 , 'm'             )
    CALL create_ISMIP_output_file_field_notime( foldername, icesheet_code, region%grid_output, 'hfgeoubed'            , 'upward_geothermal_heat_flux_at_ground_level'                      , 'W m-2'         )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'acabf'                , 'land_ice_surface_specific_mass_balance_flux'                      , 'kg m-2 s-1'    )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'libmassbfgr'          , 'land_ice_basal_specific_mass_balance_flux'                        , 'kg m-2 s-1'    )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'libmassbffl'          , 'land_ice_basal_specific_mass_balance_flux'                        , 'kg m-2 s-1'    )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'dlithkdt'             , 'tendency_of_land_ice_thickness'                                   , 'm s-1'         )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'xvelsurf'             , 'land_ice_surface_x_velocity'                                      , 'm s-1'         )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'yvelsurf'             , 'land_ice_surface_y_velocity'                                      , 'm s-1'         )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'zvelsurf'             , 'land_ice_surface_upward_velocity'                                 , 'm s-1'         )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'xvelbase'             , 'land_ice_basal_x_velocity'                                        , 'm s-1'         )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'yvelbase'             , 'land_ice_basal_y_velocity'                                        , 'm s-1'         )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'zvelbase'             , 'land_ice_basal_upward_velocity'                                   , 'm s-1'         )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'xvelmean'             , 'land_ice_vertical_mean_x_velocity'                                , 'm s-1'         )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'yvelmean'             , 'land_ice_vertical_mean_y_velocity'                                , 'm s-1'         )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'litemptop'            , 'temperature_at_top_of_ice_sheet_model'                            , 'K'             )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'litempbotgr'          , 'temperature_at_base_of_ice_sheet_model'                           , 'K'             )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'litempbotfl'          , 'temperature_at_base_of_ice_sheet_model'                           , 'K'             )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'strbasemag'           , 'land_ice_basal_drag'                                              , 'Pa'            )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'licalvf'              , 'land_ice_specific_mass_flux_due_to_calving'                       , 'kg m-2 s-1'    )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'lifmassbf'            , 'land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting' , 'kg m-2 s-1'    )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'sftgif'               , 'land_ice_area_fraction'                                           , '1'             )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'sftgrf'               , 'grounded_ice_sheet_area_fraction'                                 , '1'             )
    CALL create_ISMIP_output_file_field(        foldername, icesheet_code, region%grid_output, 'sftflf'               , 'floating_ice_shelf_area_fraction'                                 , '1'             )
    CALL create_ISMIP_output_file_scalar(       foldername, icesheet_code,              'lim'                  , 'land_ice_mass'                                                    , 'kg'            )
    CALL create_ISMIP_output_file_scalar(       foldername, icesheet_code,              'limnsw'               , 'land_ice_mass_not_displacing_sea_water'                           , 'kg'            )
    CALL create_ISMIP_output_file_scalar(       foldername, icesheet_code,              'iareagr'              , 'grounded_ice_sheet_area'                                          , 'm2'            )
    CALL create_ISMIP_output_file_scalar(       foldername, icesheet_code,              'iareafl'              , 'floating_ice_sheet_area'                                          , 'm2'            )
    CALL create_ISMIP_output_file_scalar(       foldername, icesheet_code,              'tendacabf'            , 'tendency_of_land_ice_mass_due_to_surface_mass_balance'            , 'kg s-1'        )
    CALL create_ISMIP_output_file_scalar(       foldername, icesheet_code,              'tendlibmassbf'        , 'tendency_of_land_ice_mass_due_to_basal_mass_balance'              , 'kg s-1'        )
    CALL create_ISMIP_output_file_scalar(       foldername, icesheet_code,              'tendlibmassbffl'      , 'tendency_of_land_ice_mass_due_to_basal_mass_balance'              , 'kg s-1'        )
    CALL create_ISMIP_output_file_scalar(       foldername, icesheet_code,              'tendlicalvf'          , 'tendency_of_land_ice_mass_due_to_calving'                         , 'kg s-1'        )
    CALL create_ISMIP_output_file_scalar(       foldername, icesheet_code,              'tendlifmassbf'        , 'tendency_of_land_ice_mass_due_to_calving_and_ice_front_melting'   , 'kg s-1'        )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_ISMIP_output_files

  SUBROUTINE create_ISMIP_output_file_scalar( foldername, icesheet_code, variable_name, standard_name, units)
    ! Create a single ISMIP output file for a scalar

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name, standard_name, units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_ISMIP_output_file_scalar'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_t
    INTEGER                                            :: id_var_t, id_var_scalar

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                 ) // '_' // &
                                          TRIM(                icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_experiment_code) // '.nc'

    ! If the file already exists, return.
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (file_exists) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Create netcdf file
    CALL handle_error( nf90_create( filename, IOR( nf90_clobber,nf90_share), ncid))

    ! Define dimensions:
    CALL create_dim( ncid, 'time',      nf90_unlimited,     id_dim_t      )

    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.

    ! time variable (needs some attributes that are not in the standard subroutine)
    CALL handle_error( nf90_def_var( ncid, 'time', nf90_float, [id_dim_t], id_var_t))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'standard_name', 'time'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'long_name'    , 'time'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'units'        , 'days since ' // TRIM( C%ISMIP_output_basetime)))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'calendar'     , '360_day'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'axis'         , 'T'))

    ! Field variable (needs some attributes that are not in the standard subroutine)
    CALL handle_error( nf90_def_var( ncid, variable_name, nf90_float, [id_dim_t], id_var_scalar))
    CALL handle_error( nf90_put_att( ncid, id_var_scalar, 'standard_name', standard_name))
    CALL handle_error( nf90_put_att( ncid, id_var_scalar, 'units'        , units))
    CALL handle_error( nf90_put_att( ncid, id_var_scalar, 'missing_value', 1.E20))

    ! Leave definition mode:
    CALL handle_error( nf90_enddef( ncid))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( ncid))

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_ISMIP_output_file_scalar

  SUBROUTINE create_ISMIP_output_file_field( foldername, icesheet_code, grid, variable_name, standard_name, units)
    ! Create a single ISMIP output file for an [x,y,t] data field

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name, standard_name, units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_ISMIP_output_file_field'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_t
    INTEGER                                            :: id_var_x, id_var_y, id_var_t, id_var_field

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                 ) // '_' // &
                                          TRIM(                icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_experiment_code) // '.nc'

    ! If the file already exists, return.
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (file_exists) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Create netcdf file
    CALL handle_error( nf90_create( filename, IOR( nf90_clobber,nf90_share), ncid))

    ! Define dimensions:
    CALL create_dim( ncid, 'x',         grid%nx,            id_dim_x      )
    CALL create_dim( ncid, 'y',         grid%ny,            id_dim_y      )
    CALL create_dim( ncid, 'time',      nf90_unlimited,     id_dim_t      )

    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.

    ! x,y variables
    CALL create_single_var( ncid, 'x', [id_dim_x], id_var_x, long_name = 'x-coordinate', units = 'm')
    CALL create_single_var( ncid, 'y', [id_dim_y], id_var_y, long_name = 'y-coordinate', units = 'm')

    ! time variable (needs some attributes that are not in the standard subroutine)
    CALL handle_error( nf90_def_var( ncid, 'time', nf90_float, [id_dim_t], id_var_t))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'standard_name', 'time'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'long_name'    , 'time'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'units'        , 'days since ' // TRIM( C%ISMIP_output_basetime)))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'calendar'     , '360_day'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'axis'         , 'T'))

    ! Field variable (needs some attributes that are not in the standard subroutine)
    CALL handle_error( nf90_def_var( ncid, variable_name, nf90_float, [id_dim_x, id_dim_y, id_dim_t], id_var_field))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'standard_name', standard_name))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'units'        , units))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'missing_value', 1.E20))

    ! Leave definition mode:
    CALL handle_error( nf90_enddef( ncid))

    ! Write the x, y variable data
    CALL handle_error( nf90_put_var( ncid, id_var_x, grid%x))
    CALL handle_error( nf90_put_var( ncid, id_var_y, grid%y))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( ncid))

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_ISMIP_output_file_field

  SUBROUTINE create_ISMIP_output_file_field_notime( foldername, icesheet_code, grid, variable_name, standard_name, units)
    ! Create a single ISMIP output file for an [x,y] data field

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name, standard_name, units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_ISMIP_output_file_field_notime'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: id_var_x, id_var_y, id_var_field

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                 ) // '_' // &
                                          TRIM(                icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_experiment_code) // '.nc'

    ! If the file already exists, return.
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (file_exists) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Create netcdf file
    CALL handle_error( nf90_create( filename, IOR( nf90_clobber,nf90_share), ncid))

    ! Define dimensions:
    CALL create_dim( ncid, 'x',         grid%nx,            id_dim_x      )
    CALL create_dim( ncid, 'y',         grid%ny,            id_dim_y      )

    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.

    ! x,y variables
    CALL create_single_var( ncid, 'x', [id_dim_x], id_var_x, long_name = 'x-coordinate', units = 'm')
    CALL create_single_var( ncid, 'y', [id_dim_y], id_var_y, long_name = 'y-coordinate', units = 'm')

    ! Field variable (needs some attributes that are not in the standard subroutine)
    CALL handle_error( nf90_def_var( ncid, variable_name, nf90_float, [id_dim_x, id_dim_y], id_var_field))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'standard_name', standard_name))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'units'        , units))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'missing_value', 1.E20))

    ! Leave definition mode:
    CALL handle_error( nf90_enddef( ncid))

    ! Write the x, y variable data
    CALL handle_error( nf90_put_var( ncid, id_var_x, grid%x))
    CALL handle_error( nf90_put_var( ncid, id_var_y, grid%y))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( ncid))

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_ISMIP_output_file_field_notime

  SUBROUTINE write_to_ISMIP_output_files( region)
    ! Write to all the ISMIP output files

    USE parameters_module, ONLY: ice_density, sec_per_year

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region), INTENT(IN)    :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                        :: routine_name = 'write_to_ISMIP_output_files'
    INTEGER                                              :: vi
    CHARACTER(LEN=256)                                   :: icesheet_code, foldername
    REAL(dp), PARAMETER                                  :: missing_value = 1.E20
    REAL(dp), DIMENSION( :    ), POINTER                 :: Ti_base_gr
    REAL(dp), DIMENSION( :    ), POINTER                 :: Ti_base_fl
    REAL(dp), DIMENSION( :    ), POINTER                 :: basal_drag
    REAL(dp), DIMENSION( :    ), POINTER                 :: calving_flux
    REAL(dp), DIMENSION( :    ), POINTER                 :: calving_and_front_melt_flux
    REAL(dp), DIMENSION( :    ), POINTER                 :: land_ice_area_fraction
    REAL(dp), DIMENSION( :    ), POINTER                 :: grounded_ice_sheet_area_fraction
    REAL(dp), DIMENSION( :    ), POINTER                 :: floating_ice_shelf_area_fraction
    INTEGER :: wTi_base_gr, wTi_base_fl, wbasal_drag, wcalving_flux, wcalving_and_front_melt_flux
    INTEGER :: wland_ice_area_fraction, wgrounded_ice_sheet_area_fraction, wfloating_ice_shelf_area_fraction
    REAL(dp)                                             :: land_ice_mass
    REAL(dp)                                             :: mass_above_floatation
    REAL(dp)                                             :: grounded_ice_sheet_area
    REAL(dp)                                             :: floating_ice_sheet_area
    REAL(dp)                                             :: total_SMB
    REAL(dp)                                             :: total_BMB
    REAL(dp)                                             :: total_BMB_shelf
    REAL(dp)                                             :: total_calving_flux
    REAL(dp)                                             :: total_calving_and_front_melt_flux

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Code for the ice sheet name in the ISMIP output file names
    IF     (region%name == 'NAM') THEN
      icesheet_code = 'NAIS'
    ELSEIF (region%name == 'EAS') THEN
      icesheet_code = 'EUIS'
    ELSEIF (region%name == 'GRL') THEN
      icesheet_code = 'GIS'
    ELSEIF (region%name == 'ANT') THEN
      icesheet_code = 'AIS'
    ELSEIF (region%name == 'PAT') THEN
      icesheet_code = 'PIS'
    ELSE
      icesheet_code = 'beep'
      CALL crash('unknown region "' // TRIM( region%name) // '"!')
    END IF

    ! The folder where the ISMIp6 output files are located
    foldername = TRIM( C%output_dir) // TRIM(                icesheet_code  ) // '_' // &
                                        TRIM( C%ISMIP_output_group_code     ) // '_' // &
                                        TRIM( C%ISMIP_output_model_code     ) // '_' // &
                                        TRIM( C%ISMIP_output_experiment_code)

    ! Calculate some quantities that are not natively in the ice model

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( region%mesh%nV, Ti_base_gr                      , wTi_base_gr                      )
    CALL allocate_shared_dp_1D( region%mesh%nV, Ti_base_fl                      , wTi_base_fl                      )
    CALL allocate_shared_dp_1D( region%mesh%nV, basal_drag                      , wbasal_drag                      )
    CALL allocate_shared_dp_1D( region%mesh%nV, calving_flux                    , wcalving_flux                    )
    CALL allocate_shared_dp_1D( region%mesh%nV, calving_and_front_melt_flux     , wcalving_and_front_melt_flux     )
    CALL allocate_shared_dp_1D( region%mesh%nV, land_ice_area_fraction          , wland_ice_area_fraction          )
    CALL allocate_shared_dp_1D( region%mesh%nV, grounded_ice_sheet_area_fraction, wgrounded_ice_sheet_area_fraction)
    CALL allocate_shared_dp_1D( region%mesh%nV, floating_ice_shelf_area_fraction, wfloating_ice_shelf_area_fraction)

    ! 2-D fields
    Ti_base_gr(                       region%mesh%vi1:region%mesh%vi2) = missing_value
    Ti_base_fl(                       region%mesh%vi1:region%mesh%vi2) = missing_value
    basal_drag(                       region%mesh%vi1:region%mesh%vi2) = missing_value
    calving_flux(                     region%mesh%vi1:region%mesh%vi2) = missing_value
    calving_and_front_melt_flux(      region%mesh%vi1:region%mesh%vi2) = missing_value
    land_ice_area_fraction(           region%mesh%vi1:region%mesh%vi2) = missing_value
    grounded_ice_sheet_area_fraction( region%mesh%vi1:region%mesh%vi2) = missing_value
    floating_ice_shelf_area_fraction( region%mesh%vi1:region%mesh%vi2) = missing_value

    ! Scalars (integrated values)
    land_ice_mass                     = 0._dp
    mass_above_floatation             = 0._dp
    grounded_ice_sheet_area           = 0._dp
    floating_ice_sheet_area           = 0._dp
    total_SMB                         = 0._dp
    total_BMB                         = 0._dp
    total_BMB_shelf                   = 0._dp
    total_calving_flux                = 0._dp
    total_calving_and_front_melt_flux = 0._dp

    DO vi = 1, region%mesh%vi1, region%mesh%vi2

      ! Ice base temperature separate for sheet and shelf
      ! =================================================

      IF (region%ice%mask_sheet_a( vi) == 1) THEN
        Ti_base_gr( vi) = region%ice%Ti_a( vi,C%nz)
      END IF

      IF (region%ice%mask_shelf_a( vi) == 1) THEN
        Ti_base_fl( vi) = region%ice%Ti_a( vi,C%nz)
      END IF

      ! Basal drag
      ! ==========

      IF (region%ice%mask_ice_a( vi) == 1 .AND. region%ice%f_grnd_a( vi) > 0._dp) THEN
        basal_drag( vi) = region%ice%uabs_base_a( vi) * region%ice%beta_a( vi) * region%ice%f_grnd_a( vi)**2
      END IF

      ! Calving and front melting fluxes
      ! ================================

      calving_flux(                vi) = 0._dp ! FIXME
      calving_and_front_melt_flux( vi) = 0._dp ! FIXME

      ! Ice fractions
      ! =============

      IF (region%ice%mask_cf_a( vi) == 0) THEN
        land_ice_area_fraction( vi) = REAL( region%ice%mask_ice_a( vi), dp)
      ELSE
        land_ice_area_fraction( vi) = region%ice%float_margin_frac_a( vi)
      END IF

      IF (region%ice%mask_ice_a( vi) == 1) THEN
        grounded_ice_sheet_area_fraction( vi) = region%ice%f_grnd_a( vi)
      ELSE
        grounded_ice_sheet_area_fraction( vi) = 0._dp
      END IF

      floating_ice_shelf_area_fraction( vi) = REAL( region%ice%mask_ice_a( vi), dp) * MAX( (1._dp - region%ice%f_grnd_a( vi)), region%ice%float_margin_frac_a( vi))

      ! Integrated values
      ! =================

      land_ice_mass                     = land_ice_mass                     + (region%ice%Hi_a(                  vi) * region%mesh%A( vi) * ice_density)    ! kg
      mass_above_floatation             = mass_above_floatation             + (region%ice%TAF_a(                 vi) * region%mesh%A( vi) * ice_density)    ! kg
      grounded_ice_sheet_area           = grounded_ice_sheet_area           + (grounded_ice_sheet_area_fraction( vi) * region%mesh%A( vi))                  ! m2
      floating_ice_sheet_area           = floating_ice_sheet_area           + (floating_ice_shelf_area_fraction( vi) * region%mesh%A( vi))                  ! m2
      total_SMB                         = total_SMB                         + (land_ice_area_fraction( vi) * region%SMB%SMB_year(  vi) * region%mesh%A( vi) * ice_density / sec_per_year) ! kg s-1
      total_BMB                         = total_BMB                         + (land_ice_area_fraction( vi) * region%BMB%BMB(       vi) * region%mesh%A( vi) * ice_density / sec_per_year) ! kg s-1
      total_BMB_shelf                   = total_BMB_shelf                   + (land_ice_area_fraction( vi) * region%BMB%BMB_shelf( vi) * region%mesh%A( vi) * ice_density / sec_per_year) ! kg s-1
      total_calving_flux                = total_calving_flux                + (calving_flux( vi)                                       * region%mesh%A( vi) * ice_density / sec_per_year) ! kg s-1
      total_calving_and_front_melt_flux = total_calving_and_front_melt_flux + (calving_and_front_melt_flux( vi)                        * region%mesh%A( vi) * ice_density / sec_per_year) ! kg s-1

    END DO
    CALL sync

    ! Write to all the ISMIP output files
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%Hi_a                    , 'lithk'                    )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%Hs_a                    , 'orog'                     )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%Hb_a                    , 'topg'                     )
    CALL write_to_ISMIP_output_file_field_notime(  region%mesh, region%grid_output, foldername, icesheet_code,              region%ice%GHF_a                   , 'hfgeoubed'                )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%SMB%SMB_year                , 'acabf'                    )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%BMB%BMB_sheet               , 'libmassbfgr'              )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%BMB%BMB_shelf               , 'libmassbffl'              )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%dHs_dt_a                , 'dlithkdt'                 )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%u_surf_a                , 'xvelsurf'                 )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%v_surf_a                , 'yvelsurf'                 )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%w_surf_a                , 'zvelsurf'                 )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%u_base_a                , 'xvelbase'                 )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%v_base_a                , 'yvelbase'                 )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%w_base_a                , 'zvelbase'                 )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%u_vav_a                 , 'xvelmean'                 )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%v_vav_a                 , 'yvelmean'                 )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, region%ice%Ti_a( :,1)              , 'litemptop'                )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, Ti_base_gr                         , 'litempbotgr'              )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, Ti_base_fl                         , 'litempbotfl'              )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, basal_drag                         , 'strbasemag'               )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, calving_flux                       , 'licalvf'                  )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, calving_and_front_melt_flux        , 'lifmassbf'                )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, land_ice_area_fraction             , 'sftgif'                   )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, grounded_ice_sheet_area_fraction   , 'sftgrf'                   )
    CALL write_to_ISMIP_output_file_field(         region%mesh, region%grid_output, foldername, icesheet_code, region%time, floating_ice_shelf_area_fraction   , 'sftflf'                   )
    CALL write_to_ISMIP_output_file_scalar(                                         foldername, icesheet_code, region%time, land_ice_mass                      , 'lim'                      )
    CALL write_to_ISMIP_output_file_scalar(                                         foldername, icesheet_code, region%time, mass_above_floatation              , 'limnsw'                   )
    CALL write_to_ISMIP_output_file_scalar(                                         foldername, icesheet_code, region%time, grounded_ice_sheet_area            , 'iareagr'                  )
    CALL write_to_ISMIP_output_file_scalar(                                         foldername, icesheet_code, region%time, floating_ice_sheet_area            , 'iareafl'                  )
    CALL write_to_ISMIP_output_file_scalar(                                         foldername, icesheet_code, region%time, total_SMB                          , 'tendacabf'                )
    CALL write_to_ISMIP_output_file_scalar(                                         foldername, icesheet_code, region%time, total_BMB                          , 'tendlibmassbf'            )
    CALL write_to_ISMIP_output_file_scalar(                                         foldername, icesheet_code, region%time, total_BMB_shelf                    , 'tendlibmassbffl'          )
    CALL write_to_ISMIP_output_file_scalar(                                         foldername, icesheet_code, region%time, total_calving_flux                 , 'tendlicalvf'              )
    CALL write_to_ISMIP_output_file_scalar(                                         foldername, icesheet_code, region%time, total_calving_and_front_melt_flux  , 'tendlifmassbf'            )

    ! Clean up after yourself
    CALL deallocate_shared( wTi_base_gr                      )
    CALL deallocate_shared( wTi_base_fl                      )
    CALL deallocate_shared( wbasal_drag                      )
    CALL deallocate_shared( wcalving_flux                    )
    CALL deallocate_shared( wcalving_and_front_melt_flux     )
    CALL deallocate_shared( wland_ice_area_fraction          )
    CALL deallocate_shared( wgrounded_ice_sheet_area_fraction)
    CALL deallocate_shared( wfloating_ice_shelf_area_fraction)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_ISMIP_output_files

  SUBROUTINE write_to_ISMIP_output_file_scalar( foldername, icesheet_code, time, d, variable_name)
    ! Write a single scalar to the corresponding ISMIP output file

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time
    REAL(dp),                            INTENT(IN)    :: d
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_ISMIP_output_file_scalar'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_t
    INTEGER                                            :: id_var_t, id_var
    INTEGER                                            :: nt

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                 ) // '_' // &
                                          TRIM(                icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_experiment_code) // '.nc'

    ! Open the netcdf file
    CALL open_netcdf_file( filename, ncid)

    ! Inquire for dimension IDs
    CALL inquire_dim( ncid, 'time', nt, id_dim_t)

    ! Inquire for variable IDs
    CALL inquire_single_var( ncid, 'time'       , (/ id_dim_t /), id_var_t)
    CALL inquire_single_var( ncid, variable_name, (/ id_dim_t /), id_var  )

    ! Write time
    CALL handle_error( nf90_put_var( ncid, id_var_t, days_since_ISMIP_basetime( time), start = (/ nt+1 /)))

    ! Write data to the NetCDF file
    CALL handle_error( nf90_put_var( ncid, id_var, d, start = (/ nt+1 /)))

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_ISMIP_output_file_scalar

  SUBROUTINE write_to_ISMIP_output_file_field( mesh, grid, foldername, icesheet_code, time, d, variable_name)
    ! Write a single [x,y,t] data field to the corresponding ISMIP output file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp),                            INTENT(IN)    :: time
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_ISMIP_output_file_field'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_t
    INTEGER                                            :: id_var_t, id_var
    INTEGER                                            :: nx, ny, nt

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                 ) // '_' // &
                                          TRIM(                icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_experiment_code) // '.nc'

    ! Open the netcdf file
    IF (par%master) CALL open_netcdf_file( filename, ncid)

    ! Inquire for dimension IDs
    IF (par%master) CALL inquire_dim( ncid, 'x'   , nx, id_dim_x)
    IF (par%master) CALL inquire_dim( ncid, 'y'   , ny, id_dim_y)
    IF (par%master) CALL inquire_dim( ncid, 'time', nt, id_dim_t)

    ! Inquire for variable IDs
    IF (par%master) CALL inquire_single_var( ncid, 'time'       , (/                     id_dim_t /), id_var_t)
    IF (par%master) CALL inquire_single_var( ncid, variable_name, (/ id_dim_x, id_dim_y, id_dim_t /), id_var  )

    ! Write time
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var_t, days_since_ISMIP_basetime( time), start = (/ nt+1 /)))

    ! Write data to the NetCDF file
    CALL map_and_write_to_grid_netcdf_dp_2D( ncid, mesh, grid, d, id_var, nt+1)

    ! Close the file
    IF (par%master) CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_ISMIP_output_file_field

  SUBROUTINE write_to_ISMIP_output_file_field_notime( mesh, grid, foldername, icesheet_code, d, variable_name)
    ! Write a single [x,y,t] data field to the corresponding ISMIP output file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_ISMIP_output_file_field_notime'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: id_var
    INTEGER                                            :: nx, ny

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                 ) // '_' // &
                                          TRIM(                icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP_output_experiment_code) // '.nc'

    ! Open the netcdf file
    IF (par%master) CALL open_netcdf_file( filename, ncid)

    ! Inquire for dimension IDs
    IF (par%master) CALL inquire_dim( ncid, 'x'   , nx, id_dim_x)
    IF (par%master) CALL inquire_dim( ncid, 'y'   , ny, id_dim_y)

    ! Inquire for variable IDs
    IF (par%master) CALL inquire_single_var( ncid, variable_name, (/ id_dim_x, id_dim_y /), id_var  )

    ! Write data to the NetCDF file
    CALL map_and_write_to_grid_netcdf_dp_2D_notime( ncid, mesh, grid, d, id_var)

    ! Close the file
    IF (par%master) CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_ISMIP_output_file_field_notime

  FUNCTION ISMIP_unit_conversion_field( d, variable_name) RESULT( d_conv)
    ! Convert data fields from IMAU-ICE units to SI units

    USE parameters_module, ONLY: sec_per_year, ice_density

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d
    CHARACTER(LEN=*),                    INTENT(IN)    :: variable_name
    REAL(dp), DIMENSION(SIZE(d,1))                     :: d_conv

    IF     (variable_name == 'lithk') THEN
      ! land_ice_thickness
      ! Ice model units: m
      ! SI units: m
      d_conv = d
    ELSEIF (variable_name == 'orog') THEN
      ! surface_altitude
      ! Ice model units: m
      ! SI units: m
      d_conv = d
    ELSEIF (variable_name == 'topg') THEN
      ! bedrock_altitude
      ! Ice model units: m
      ! SI units: m
      d_conv = d
    ELSEIF (variable_name == 'hfgeoubed') THEN
      ! upward_geothermal_heat_flux_at_ground_level
      ! Ice model units: J m-2 yr-1
      ! SI units: J m-2 s-1
      d_conv = d / sec_per_year
    ELSEIF (variable_name == 'acabf') THEN
      ! land_ice_surface_specific_mass_balance_flux
      ! Ice model units: m.i.e./yr
      ! SI units: kg m-2 s-1
      d_conv = d * ice_density / sec_per_year
    ELSEIF (variable_name == 'libmassbfgr') THEN
      ! land_ice_basal_specific_mass_balance_flux
      ! Ice model units: m.i.e./yr
      ! SI units: kg m-2 s-1
      d_conv = d * ice_density / sec_per_year
    ELSEIF (variable_name == 'libmassbffl') THEN
      ! land_ice_basal_specific_mass_balance_flux
      ! Ice model units: m.i.e./yr
      ! SI units: kg m-2 s-1
      d_conv = d * ice_density / sec_per_year
    ELSEIF (variable_name == 'dlithkdt') THEN
      ! tendency_of_land_ice_thickness
      ! Ice model units: m/yr
      ! SI units: m s-1
      d_conv = d / sec_per_year
    ELSEIF (variable_name == 'xvelsurf' .OR. &
            variable_name == 'yvelsurf' .OR. &
            variable_name == 'zvelsurf' .OR. &
            variable_name == 'xvelbase' .OR. &
            variable_name == 'yvelbase' .OR. &
            variable_name == 'zvelbase' .OR. &
            variable_name == 'xvelmean' .OR. &
            variable_name == 'yvelmean') THEN
      ! different ice velocities
      ! Ice model units: m/yr
      ! SI units: m s-1
      d_conv = d / sec_per_year
    ELSEIF (variable_name == 'litemptop') THEN
      ! temperature_at_top_of_ice_sheet_model
      ! Ice model units: K
      ! SI units: K
      d_conv = d
    ELSEIF (variable_name == 'litempbotgr') THEN
      ! temperature_at_base_of_ice_sheet_model
      ! Ice model units: K
      ! SI units: K
      d_conv = d
    ELSEIF (variable_name == 'litempbotfl') THEN
      ! temperature_at_base_of_ice_sheet_model
      ! Ice model units: K
      ! SI units: K
      d_conv = d
    ELSEIF (variable_name == 'strbasemag') THEN
      ! land_ice_basal_drag
      ! Ice model units: Pa
      ! SI units: Pa
      d_conv = d
    ELSEIF (variable_name == 'licalvf' .OR. &
            variable_name == 'lifmassbf') THEN
      ! land_ice_specific_mass_flux_due_to_calving (possibly also front melting)
      ! Ice model units: m/yr
      ! SI units: kg m-2 s-1
      d_conv = d * ice_density / sec_per_year
    ELSEIF (variable_name == 'sftgif' .OR. &
            variable_name == 'sftgrf' .OR. &
            variable_name == 'sftflf') THEN
      ! ice fractions
      ! Ice model units:
      ! SI units:
      d_conv = d
    ELSE
      ! Unknown variable name
      CALL crash('ISMIP_unit_conversion: unknown variable name "' // TRIM( variable_name) // '"!')
    END IF

  END FUNCTION ISMIP_unit_conversion_field

  FUNCTION days_since_ISMIP_basetime( time) RESULT( ndays)
    ! Calculate the number of days since ISMIP basetime
    !
    ! Assume basetime equals t = 0

    USE parameters_module, ONLY: sec_per_year

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time
    REAL(dp)                                           :: ndays

    ndays = time * 360._dp

  END FUNCTION days_since_ISMIP_basetime

! ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
! ============================================================

  SUBROUTINE inquire_ISMIP_forcing_SMB_baseline_file( netcdf)
    ! Check if the right dimensions and variables are present in the file.

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_ISMIP_forcing_SMB_baseline_file'
    INTEGER                                       :: nx,ny

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_x, nx, netcdf%id_dim_x)
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_y, ny, netcdf%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_x,   (/ netcdf%id_dim_x                  /), netcdf%id_var_x  )
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_y,   (/                  netcdf%id_dim_y /), netcdf%id_var_y  )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_SMB, (/ netcdf%id_dim_x, netcdf%id_dim_y /), netcdf%id_var_SMB)

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_ISMIP_forcing_SMB_baseline_file

  SUBROUTINE read_ISMIP_forcing_SMB_baseline_file( netcdf, SMB)
    ! Read grid and data from the NetCDF file

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing),   INTENT(INOUT) :: netcdf
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_ISMIP_forcing_SMB_baseline_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Read the field data
    CALL handle_error( nf90_get_var( netcdf%ncid, netcdf%id_var_SMB, SMB, start = (/ 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_ISMIP_forcing_SMB_baseline_file

  SUBROUTINE inquire_ISMIP_forcing_ST_baseline_file( netcdf)
    ! Check if the right dimensions and variables are present in the file.

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_ISMIP_forcing_ST_baseline_file'
    INTEGER                                       :: nx,ny

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_x, nx, netcdf%id_dim_x)
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_y, ny, netcdf%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_x,   (/ netcdf%id_dim_x                  /), netcdf%id_var_x  )
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_y,   (/                  netcdf%id_dim_y /), netcdf%id_var_y  )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_ST,  (/ netcdf%id_dim_x, netcdf%id_dim_y /), netcdf%id_var_ST )

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_ISMIP_forcing_ST_baseline_file

  SUBROUTINE read_ISMIP_forcing_ST_baseline_file( netcdf, ST)
    ! Read grid and data from the NetCDF file

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing),   INTENT(INOUT) :: netcdf
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: ST

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_ISMIP_forcing_ST_baseline_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Read the field data
    CALL handle_error( nf90_get_var( netcdf%ncid, netcdf%id_var_ST, ST, start = (/ 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_ISMIP_forcing_ST_baseline_file

  SUBROUTINE inquire_ISMIP_forcing_aSMB_file( netcdf)
    ! Check if the right dimensions and variables are present in the file.

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_ISMIP_forcing_aSMB_file'
    INTEGER                                       :: nx, ny, nt

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_x   , nx, netcdf%id_dim_x   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_y   , ny, netcdf%id_dim_y   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_time, nt, netcdf%id_dim_time)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_x,    (/ netcdf%id_dim_x                                      /), netcdf%id_var_x   )
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_y,    (/                  netcdf%id_dim_y                     /), netcdf%id_var_y   )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_aSMB, (/ netcdf%id_dim_x, netcdf%id_dim_y, netcdf%id_dim_time /), netcdf%id_var_aSMB)

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_ISMIP_forcing_aSMB_file

  SUBROUTINE read_ISMIP_forcing_aSMB_file( netcdf, aSMB)
    ! Read grid and data from the NetCDF file

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing), INTENT(INOUT) :: netcdf
    REAL(dp), DIMENSION(:,:  ),            INTENT(INOUT) :: aSMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_ISMIP_forcing_aSMB_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Read the field data
    CALL handle_error( nf90_get_var( netcdf%ncid, netcdf%id_var_aSMB, aSMB, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_ISMIP_forcing_aSMB_file

  SUBROUTINE inquire_ISMIP_forcing_dSMBdz_file( netcdf)
    ! Check if the right dimensions and variables are present in the file.

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_ISMIP_forcing_dSMBdz_file'
    INTEGER                                       :: nx, ny, nt

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_x   , nx, netcdf%id_dim_x   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_y   , ny, netcdf%id_dim_y   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_time, nt, netcdf%id_dim_time)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_x,      (/ netcdf%id_dim_x                                      /), netcdf%id_var_x   )
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_y,      (/                  netcdf%id_dim_y                     /), netcdf%id_var_y   )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_dSMBdz, (/ netcdf%id_dim_x, netcdf%id_dim_y, netcdf%id_dim_time /), netcdf%id_var_dSMBdz)

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_ISMIP_forcing_dSMBdz_file

  SUBROUTINE read_ISMIP_forcing_dSMBdz_file( netcdf, dSMBdz)
    ! Read grid and data from the NetCDF file

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing), INTENT(INOUT) :: netcdf
    REAL(dp), DIMENSION(:,:  ),            INTENT(INOUT) :: dSMBdz

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_ISMIP_forcing_dSMBdz_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Read the field data
    CALL handle_error( nf90_get_var( netcdf%ncid, netcdf%id_var_dSMBdz, dSMBdz, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_ISMIP_forcing_dSMBdz_file

  SUBROUTINE inquire_ISMIP_forcing_aST_file( netcdf)
    ! Check if the right dimensions and variables are present in the file.

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_ISMIP_forcing_aST_file'
    INTEGER                                       :: nx, ny, nt

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_x   , nx, netcdf%id_dim_x   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_y   , ny, netcdf%id_dim_y   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_time, nt, netcdf%id_dim_time)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_x,    (/ netcdf%id_dim_x                                      /), netcdf%id_var_x   )
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_y,    (/                  netcdf%id_dim_y                     /), netcdf%id_var_y   )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_aST, (/ netcdf%id_dim_x, netcdf%id_dim_y, netcdf%id_dim_time /), netcdf%id_var_aST)

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_ISMIP_forcing_aST_file

  SUBROUTINE read_ISMIP_forcing_aST_file( netcdf, aST)
    ! Read grid and data from the NetCDF file

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing), INTENT(INOUT) :: netcdf
    REAL(dp), DIMENSION(:,:  ),            INTENT(INOUT) :: aST

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_ISMIP_forcing_aST_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Read the field data
    CALL handle_error( nf90_get_var( netcdf%ncid, netcdf%id_var_aST, aST, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_ISMIP_forcing_aST_file

  SUBROUTINE inquire_ISMIP_forcing_dSTdz_file( netcdf)
    ! Check if the right dimensions and variables are present in the file.

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_ISMIP_forcing_dSTdz_file'
    INTEGER                                       :: nx, ny, nt

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_x   , nx, netcdf%id_dim_x   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_y   , ny, netcdf%id_dim_y   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_time, nt, netcdf%id_dim_time)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_x,     (/ netcdf%id_dim_x                                      /), netcdf%id_var_x   )
    CALL inquire_single_var( netcdf%ncid, netcdf%name_var_y,     (/                  netcdf%id_dim_y                     /), netcdf%id_var_y   )
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_dSTdz, (/ netcdf%id_dim_x, netcdf%id_dim_y, netcdf%id_dim_time /), netcdf%id_var_dSTdz)

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_ISMIP_forcing_dSTdz_file

  SUBROUTINE read_ISMIP_forcing_dSTdz_file( netcdf, dSTdz)
    ! Read grid and data from the NetCDF file

    ! In/output variables:
    TYPE(type_netcdf_ISMIP_style_forcing), INTENT(INOUT) :: netcdf
    REAL(dp), DIMENSION(:,:  ),            INTENT(INOUT) :: dSTdz

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_ISMIP_forcing_dSTdz_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Read the field data
    CALL handle_error( nf90_get_var( netcdf%ncid, netcdf%id_var_dSTdz, dSTdz, start = (/ 1, 1, 1 /) ))

    ! Close the netcdf file
    CALL close_netcdf_file( netcdf%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_ISMIP_forcing_dSTdz_file

! ===== Create and write to debug NetCDF file =====
! =================================================

  SUBROUTINE write_to_debug_file
    ! Write the current set of debug data fields to the debug NetCDF file

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER          :: routine_name = 'write_to_debug_file'
    INTEGER                                :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (.NOT. C%do_write_debug_data) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the file for writing
    CALL open_netcdf_file( debug%netcdf%filename, ncid)

    ! Write data
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_01, debug%int_2D_a_01, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_02, debug%int_2D_a_02, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_03, debug%int_2D_a_03, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_04, debug%int_2D_a_04, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_05, debug%int_2D_a_05, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_06, debug%int_2D_a_06, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_07, debug%int_2D_a_07, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_08, debug%int_2D_a_08, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_09, debug%int_2D_a_09, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_a_10, debug%int_2D_a_10, start = (/ 1 /) ))

    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_01, debug%int_2D_b_01, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_02, debug%int_2D_b_02, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_03, debug%int_2D_b_03, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_04, debug%int_2D_b_04, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_05, debug%int_2D_b_05, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_06, debug%int_2D_b_06, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_07, debug%int_2D_b_07, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_08, debug%int_2D_b_08, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_09, debug%int_2D_b_09, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_b_10, debug%int_2D_b_10, start = (/ 1 /) ))

    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_01, debug%int_2D_c_01, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_02, debug%int_2D_c_02, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_03, debug%int_2D_c_03, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_04, debug%int_2D_c_04, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_05, debug%int_2D_c_05, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_06, debug%int_2D_c_06, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_07, debug%int_2D_c_07, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_08, debug%int_2D_c_08, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_09, debug%int_2D_c_09, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_int_2D_c_10, debug%int_2D_c_10, start = (/ 1 /) ))

    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_01, debug%dp_2D_a_01, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_02, debug%dp_2D_a_02, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_03, debug%dp_2D_a_03, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_04, debug%dp_2D_a_04, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_05, debug%dp_2D_a_05, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_06, debug%dp_2D_a_06, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_07, debug%dp_2D_a_07, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_08, debug%dp_2D_a_08, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_09, debug%dp_2D_a_09, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_a_10, debug%dp_2D_a_10, start = (/ 1 /) ))

    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_01, debug%dp_2D_b_01, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_02, debug%dp_2D_b_02, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_03, debug%dp_2D_b_03, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_04, debug%dp_2D_b_04, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_05, debug%dp_2D_b_05, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_06, debug%dp_2D_b_06, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_07, debug%dp_2D_b_07, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_08, debug%dp_2D_b_08, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_09, debug%dp_2D_b_09, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_b_10, debug%dp_2D_b_10, start = (/ 1 /) ))

    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_01, debug%dp_2D_c_01, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_02, debug%dp_2D_c_02, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_03, debug%dp_2D_c_03, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_04, debug%dp_2D_c_04, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_05, debug%dp_2D_c_05, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_06, debug%dp_2D_c_06, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_07, debug%dp_2D_c_07, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_08, debug%dp_2D_c_08, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_09, debug%dp_2D_c_09, start = (/ 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_c_10, debug%dp_2D_c_10, start = (/ 1 /) ))

    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_01, debug%dp_3D_a_01, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_02, debug%dp_3D_a_02, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_03, debug%dp_3D_a_03, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_04, debug%dp_3D_a_04, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_05, debug%dp_3D_a_05, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_06, debug%dp_3D_a_06, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_07, debug%dp_3D_a_07, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_08, debug%dp_3D_a_08, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_09, debug%dp_3D_a_09, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_3D_a_10, debug%dp_3D_a_10, start = (/ 1, 1 /) ))

    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_01, debug%dp_2D_monthly_a_01, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_02, debug%dp_2D_monthly_a_02, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_03, debug%dp_2D_monthly_a_03, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_04, debug%dp_2D_monthly_a_04, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_05, debug%dp_2D_monthly_a_05, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_06, debug%dp_2D_monthly_a_06, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_07, debug%dp_2D_monthly_a_07, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_08, debug%dp_2D_monthly_a_08, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_09, debug%dp_2D_monthly_a_09, start = (/ 1, 1 /) ))
    CALL handle_error( nf90_put_var( ncid, debug%netcdf%id_var_dp_2D_monthly_a_10, debug%dp_2D_monthly_a_10, start = (/ 1, 1 /) ))

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_debug_file

  SUBROUTINE create_debug_file( region)
    ! Create the debug NetCDF file; a lot of data fields but no time dimension.

    USE data_types_netcdf_module, ONLY: type_netcdf_debug

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_debug_file'
    TYPE(type_netcdf_debug)                       :: debug_temp
    CHARACTER(LEN=20)                             :: short_filename
    INTEGER                                       :: n
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, six, vii, zeta, month

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (.NOT. C%do_write_debug_data) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Determine debug NetCDF filename for this model region
    short_filename = 'debug_NAM.nc'
    short_filename(7:9) = region%name
    DO n = 1, 256
      debug_temp%filename(n:n) = ' '
    END DO
    debug_temp%filename = TRIM(C%output_dir)//TRIM(short_filename)

    ! Delete existing debug file
    INQUIRE(EXIST=file_exists, FILE = TRIM(debug_temp%filename))
    IF (file_exists) THEN
      CALL system('rm -f ' // debug_temp%filename)
    END IF

    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( debug_temp%filename)
    CALL handle_error(nf90_create(debug_temp%filename,IOR(nf90_clobber,nf90_share),debug_temp%ncid))

    ! Mesh data
    ! =========

    ! Define dimensions
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_vi,           region%mesh%nV,          debug_temp%id_dim_vi          ) ! Vertex indices
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_ti,           region%mesh%nTri,        debug_temp%id_dim_ti          ) ! Triangle indices
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_ci,           region%mesh%nC_mem,      debug_temp%id_dim_ci          ) ! Connection indices
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_aci,          region%mesh%nAc,         debug_temp%id_dim_aci         ) ! Staggered vertex indices
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_ciplusone,    region%mesh%nC_mem+1,    debug_temp%id_dim_ciplusone   ) ! connection indices plus one (neighbour function arrays have one more column)
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_two,          2,                       debug_temp%id_dim_two         ) ! 2 (each vertex has an X and Y coordinates)
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_three,        3,                       debug_temp%id_dim_three       ) ! 3 (each triangle has three vertices)
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_six,          6,                       debug_temp%id_dim_six         ) ! 4 (each staggered vertex lists four regular vertices and two triangles)
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_vii_transect, region%mesh%nV_transect, debug_temp%id_dim_vii_transect) ! Number of vertex pairs in the transect

    ! Placeholders for the dimension ID's, for shorter code
    vi        = debug_temp%id_dim_vi
    ti        = debug_temp%id_dim_ti
    ci        = debug_temp%id_dim_ci
    aci       = debug_temp%id_dim_aci
    ciplusone = debug_temp%id_dim_ciplusone
    two       = debug_temp%id_dim_two
    three     = debug_temp%id_dim_three
    six       = debug_temp%id_dim_six
    vii       = debug_temp%id_dim_vii_transect

    ! Define variables
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_V,                [vi,  two  ], debug_temp%id_var_V,                long_name='Vertex coordinates', units='m')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_Tri,              [ti,  three], debug_temp%id_var_Tri,              long_name='Vertex indices')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_nC,               [vi        ], debug_temp%id_var_nC,               long_name='Number of connected vertices')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_C,                [vi,  ci   ], debug_temp%id_var_C,                long_name='Indices of connected vertices')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_niTri,            [vi        ], debug_temp%id_var_niTri,            long_name='Number of inverse triangles')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_iTri,             [vi,  ci   ], debug_temp%id_var_iTri,             long_name='Indices of inverse triangles')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_edge_index,       [vi        ], debug_temp%id_var_edge_index,       long_name='Edge index')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_Tricc,            [ti,  two  ], debug_temp%id_var_Tricc,            long_name='Triangle circumcenter', units='m')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_TriC,             [ti,  three], debug_temp%id_var_TriC,             long_name='Triangle neighbours')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_Tri_edge_index,   [ti        ], debug_temp%id_var_Tri_edge_index,   long_name='Triangle edge index')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_VAc,              [aci, two  ], debug_temp%id_var_VAc,              long_name='Staggered vertex coordinates', units='m')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_Aci,              [aci, six  ], debug_temp%id_var_Aci,              long_name='Staggered to regular vertex indices')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_iAci,             [vi,  ci   ], debug_temp%id_var_iAci,             long_name='Regular to staggered vertex indices')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_A,                [vi        ], debug_temp%id_var_A,                long_name='Vertex Voronoi cell area', units='m^2')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_R,                [vi        ], debug_temp%id_var_R,                long_name='Vertex resolution', units='m')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_vi_transect,      [vii, two  ], debug_temp%id_var_vi_transect,      long_name='Transect vertex pairs')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_w_transect,       [vii, two  ], debug_temp%id_var_w_transect,       long_name='Transect interpolation weights')

    ! Model output
    ! ============

    ! Define dimensions
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_zeta,  C%nZ, debug_temp%id_dim_zeta ) ! Scaled vertical coordinate
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_month, 12,   debug_temp%id_dim_month) ! Months (for monthly data)

    ! Placeholders for the dimension ID's, for shorter code
    zeta  = debug_temp%id_dim_zeta
    month = debug_temp%id_dim_month

    ! Define dimension variables
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_zeta,  [zeta  ], debug_temp%id_var_zeta,  long_name='Vertical scaled coordinate', units='unitless (0 = ice surface, 1 = bedrock)')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_month, [month ], debug_temp%id_var_month, long_name='Month', units='1-12')

    ! Data
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_01,  [vi], debug_temp%id_var_int_2D_a_01,  long_name='2D int a-grid (vertex) variable 01')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_02,  [vi], debug_temp%id_var_int_2D_a_02,  long_name='2D int a-grid (vertex) variable 02')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_03,  [vi], debug_temp%id_var_int_2D_a_03,  long_name='2D int a-grid (vertex) variable 03')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_04,  [vi], debug_temp%id_var_int_2D_a_04,  long_name='2D int a-grid (vertex) variable 04')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_05,  [vi], debug_temp%id_var_int_2D_a_05,  long_name='2D int a-grid (vertex) variable 05')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_06,  [vi], debug_temp%id_var_int_2D_a_06,  long_name='2D int a-grid (vertex) variable 06')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_07,  [vi], debug_temp%id_var_int_2D_a_07,  long_name='2D int a-grid (vertex) variable 07')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_08,  [vi], debug_temp%id_var_int_2D_a_08,  long_name='2D int a-grid (vertex) variable 08')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_09,  [vi], debug_temp%id_var_int_2D_a_09,  long_name='2D int a-grid (vertex) variable 09')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_a_10,  [vi], debug_temp%id_var_int_2D_a_10,  long_name='2D int a-grid (vertex) variable 10')

    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_01,  [ti], debug_temp%id_var_int_2D_b_01,  long_name='2D int b-grid (triangle) variable 01')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_02,  [ti], debug_temp%id_var_int_2D_b_02,  long_name='2D int b-grid (triangle) variable 02')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_03,  [ti], debug_temp%id_var_int_2D_b_03,  long_name='2D int b-grid (triangle) variable 03')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_04,  [ti], debug_temp%id_var_int_2D_b_04,  long_name='2D int b-grid (triangle) variable 04')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_05,  [ti], debug_temp%id_var_int_2D_b_05,  long_name='2D int b-grid (triangle) variable 05')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_06,  [ti], debug_temp%id_var_int_2D_b_06,  long_name='2D int b-grid (triangle) variable 06')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_07,  [ti], debug_temp%id_var_int_2D_b_07,  long_name='2D int b-grid (triangle) variable 07')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_08,  [ti], debug_temp%id_var_int_2D_b_08,  long_name='2D int b-grid (triangle) variable 08')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_09,  [ti], debug_temp%id_var_int_2D_b_09,  long_name='2D int b-grid (triangle) variable 09')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_b_10,  [ti], debug_temp%id_var_int_2D_b_10,  long_name='2D int b-grid (triangle) variable 10')

    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_01, [aci], debug_temp%id_var_int_2D_c_01,  long_name='2D int c-grid (edge) variable 01')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_02, [aci], debug_temp%id_var_int_2D_c_02,  long_name='2D int c-grid (edge) variable 02')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_03, [aci], debug_temp%id_var_int_2D_c_03,  long_name='2D int c-grid (edge) variable 03')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_04, [aci], debug_temp%id_var_int_2D_c_04,  long_name='2D int c-grid (edge) variable 04')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_05, [aci], debug_temp%id_var_int_2D_c_05,  long_name='2D int c-grid (edge) variable 05')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_06, [aci], debug_temp%id_var_int_2D_c_06,  long_name='2D int c-grid (edge) variable 06')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_07, [aci], debug_temp%id_var_int_2D_c_07,  long_name='2D int c-grid (edge) variable 07')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_08, [aci], debug_temp%id_var_int_2D_c_08,  long_name='2D int c-grid (edge) variable 08')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_09, [aci], debug_temp%id_var_int_2D_c_09,  long_name='2D int c-grid (edge) variable 09')
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_int_2D_c_10, [aci], debug_temp%id_var_int_2D_c_10,  long_name='2D int c-grid (edge) variable 10')

    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_01,  [vi], debug_temp%id_var_dp_2D_a_01,  long_name='2D dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_02,  [vi], debug_temp%id_var_dp_2D_a_02,  long_name='2D dp a-grid (vertex) variable 02')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_03,  [vi], debug_temp%id_var_dp_2D_a_03,  long_name='2D dp a-grid (vertex) variable 03')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_04,  [vi], debug_temp%id_var_dp_2D_a_04,  long_name='2D dp a-grid (vertex) variable 04')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_05,  [vi], debug_temp%id_var_dp_2D_a_05,  long_name='2D dp a-grid (vertex) variable 05')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_06,  [vi], debug_temp%id_var_dp_2D_a_06,  long_name='2D dp a-grid (vertex) variable 06')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_07,  [vi], debug_temp%id_var_dp_2D_a_07,  long_name='2D dp a-grid (vertex) variable 07')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_08,  [vi], debug_temp%id_var_dp_2D_a_08,  long_name='2D dp a-grid (vertex) variable 08')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_09,  [vi], debug_temp%id_var_dp_2D_a_09,  long_name='2D dp a-grid (vertex) variable 09')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_a_10,  [vi], debug_temp%id_var_dp_2D_a_10,  long_name='2D dp a-grid (vertex) variable 10')

    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_01,  [ti], debug_temp%id_var_dp_2D_b_01,  long_name='2D dp b-grid (triangle) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_02,  [ti], debug_temp%id_var_dp_2D_b_02,  long_name='2D dp b-grid (triangle) variable 02')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_03,  [ti], debug_temp%id_var_dp_2D_b_03,  long_name='2D dp b-grid (triangle) variable 03')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_04,  [ti], debug_temp%id_var_dp_2D_b_04,  long_name='2D dp b-grid (triangle) variable 04')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_05,  [ti], debug_temp%id_var_dp_2D_b_05,  long_name='2D dp b-grid (triangle) variable 05')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_06,  [ti], debug_temp%id_var_dp_2D_b_06,  long_name='2D dp b-grid (triangle) variable 06')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_07,  [ti], debug_temp%id_var_dp_2D_b_07,  long_name='2D dp b-grid (triangle) variable 07')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_08,  [ti], debug_temp%id_var_dp_2D_b_08,  long_name='2D dp b-grid (triangle) variable 08')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_09,  [ti], debug_temp%id_var_dp_2D_b_09,  long_name='2D dp b-grid (triangle) variable 09')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_b_10,  [ti], debug_temp%id_var_dp_2D_b_10,  long_name='2D dp b-grid (triangle) variable 10')

    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_01, [aci], debug_temp%id_var_dp_2D_c_01,  long_name='2D dp c-grid (edge) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_02, [aci], debug_temp%id_var_dp_2D_c_02,  long_name='2D dp c-grid (edge) variable 02')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_03, [aci], debug_temp%id_var_dp_2D_c_03,  long_name='2D dp c-grid (edge) variable 03')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_04, [aci], debug_temp%id_var_dp_2D_c_04,  long_name='2D dp c-grid (edge) variable 04')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_05, [aci], debug_temp%id_var_dp_2D_c_05,  long_name='2D dp c-grid (edge) variable 05')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_06, [aci], debug_temp%id_var_dp_2D_c_06,  long_name='2D dp c-grid (edge) variable 06')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_07, [aci], debug_temp%id_var_dp_2D_c_07,  long_name='2D dp c-grid (edge) variable 07')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_08, [aci], debug_temp%id_var_dp_2D_c_08,  long_name='2D dp c-grid (edge) variable 08')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_09, [aci], debug_temp%id_var_dp_2D_c_09,  long_name='2D dp c-grid (edge) variable 09')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_c_10, [aci], debug_temp%id_var_dp_2D_c_10,  long_name='2D dp c-grid (edge) variable 10')

    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_01,  [vi, zeta], debug_temp%id_var_dp_3D_a_01,  long_name='3D dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_02,  [vi, zeta], debug_temp%id_var_dp_3D_a_02,  long_name='3D dp a-grid (vertex) variable 02')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_03,  [vi, zeta], debug_temp%id_var_dp_3D_a_03,  long_name='3D dp a-grid (vertex) variable 03')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_04,  [vi, zeta], debug_temp%id_var_dp_3D_a_04,  long_name='3D dp a-grid (vertex) variable 04')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_05,  [vi, zeta], debug_temp%id_var_dp_3D_a_05,  long_name='3D dp a-grid (vertex) variable 05')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_06,  [vi, zeta], debug_temp%id_var_dp_3D_a_06,  long_name='3D dp a-grid (vertex) variable 06')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_07,  [vi, zeta], debug_temp%id_var_dp_3D_a_07,  long_name='3D dp a-grid (vertex) variable 07')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_08,  [vi, zeta], debug_temp%id_var_dp_3D_a_08,  long_name='3D dp a-grid (vertex) variable 08')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_09,  [vi, zeta], debug_temp%id_var_dp_3D_a_09,  long_name='3D dp a-grid (vertex) variable 09')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_3D_a_10,  [vi, zeta], debug_temp%id_var_dp_3D_a_10,  long_name='3D dp a-grid (vertex) variable 10')

    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_01,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_01,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_02,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_02,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_03,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_03,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_04,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_04,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_05,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_05,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_06,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_06,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_07,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_07,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_08,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_08,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_09,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_09,  long_name='2D-monthly dp a-grid (vertex) variable 01')
    CALL create_double_var( debug_temp%ncid, debug_temp%name_var_dp_2D_monthly_a_10,  [vi, month], debug_temp%id_var_dp_2D_monthly_a_10,  long_name='2D-monthly dp a-grid (vertex) variable 01')

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( debug_temp%ncid))

    ! Write mesh data
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_V,               region%mesh%V             ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_Tri,             region%mesh%Tri           ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_nC,              region%mesh%nC            ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_C,               region%mesh%C             ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_niTri,           region%mesh%niTri         ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_iTri,            region%mesh%iTri          ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_edge_index,      region%mesh%edge_index    ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_Tricc,           region%mesh%Tricc         ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_TriC,            region%mesh%TriC          ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_Tri_edge_index,  region%mesh%Tri_edge_index))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_VAc,             region%mesh%VAc           ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_Aci,             region%mesh%Aci           ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_iAci,            region%mesh%iAci          ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_A,               region%mesh%A             ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_R,               region%mesh%R             ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_vi_transect,     region%mesh%vi_transect   ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_w_transect,      region%mesh%w_transect    ))

    ! Write zeta and month dimension variables
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( debug_temp%ncid, debug_temp%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( debug_temp%ncid))

    ! Close the file
    CALL close_netcdf_file( debug_temp%ncid)

    ! Copy NetCDF data to the relevant debug structure
    IF     (region%name == 'NAM') THEN
      debug_NAM%netcdf = debug_temp
    ELSEIF (region%name == 'EAS') THEN
      debug_EAS%netcdf = debug_temp
    ELSEIF (region%name == 'GRL') THEN
      debug_GRL%netcdf = debug_temp
    ELSEIF (region%name == 'ANT') THEN
      debug_ANT%netcdf = debug_temp
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_debug_file

  ! Manage memory for the debug data fields
  SUBROUTINE associate_debug_fields( region)
    ! Since the dimensions vary, each region needs its own set of debug fields. However, if
    ! we make them part of the "region" TYPE, they need to be passed to every subroutine as an
    ! argument before they can be used, which is a lot of hassle. So instead they are saved as
    ! global variables of this module, where they can be accessed from anywhere. This is done
    ! via the "intermediary" set of pointers, which are bound to the region-specific debug structure
    ! with this here subroutine.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),        INTENT(IN)    :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'associate_debug_fields'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Copy the netcdf ID's
    IF (par%master) THEN
      IF     (region%name == 'NAM') THEN
        debug%netcdf = debug_NAM%netcdf
      ELSEIF (region%name == 'EAS') THEN
        debug%netcdf = debug_EAS%netcdf
      ELSEIF (region%name == 'GRL') THEN
        debug%netcdf = debug_GRL%netcdf
      ELSEIF (region%name == 'ANT') THEN
        debug%netcdf = debug_ANT%netcdf
      END IF
    END IF
    CALL sync

    ! If necessary (i.e. every time except the first ever time this subroutine is called), de-associate the intermediary pointers first.
    IF (ASSOCIATED(debug%dp_2D_a_01)) THEN

      NULLIFY( debug%int_2D_a_01)
      NULLIFY( debug%int_2D_a_02)
      NULLIFY( debug%int_2D_a_03)
      NULLIFY( debug%int_2D_a_04)
      NULLIFY( debug%int_2D_a_05)
      NULLIFY( debug%int_2D_a_06)
      NULLIFY( debug%int_2D_a_07)
      NULLIFY( debug%int_2D_a_08)
      NULLIFY( debug%int_2D_a_09)
      NULLIFY( debug%int_2D_a_10)

      NULLIFY( debug%int_2D_b_01)
      NULLIFY( debug%int_2D_b_02)
      NULLIFY( debug%int_2D_b_03)
      NULLIFY( debug%int_2D_b_04)
      NULLIFY( debug%int_2D_b_05)
      NULLIFY( debug%int_2D_b_06)
      NULLIFY( debug%int_2D_b_07)
      NULLIFY( debug%int_2D_b_08)
      NULLIFY( debug%int_2D_b_09)
      NULLIFY( debug%int_2D_b_10)

      NULLIFY( debug%int_2D_c_01)
      NULLIFY( debug%int_2D_c_02)
      NULLIFY( debug%int_2D_c_03)
      NULLIFY( debug%int_2D_c_04)
      NULLIFY( debug%int_2D_c_05)
      NULLIFY( debug%int_2D_c_06)
      NULLIFY( debug%int_2D_c_07)
      NULLIFY( debug%int_2D_c_08)
      NULLIFY( debug%int_2D_c_09)
      NULLIFY( debug%int_2D_c_10)

      NULLIFY( debug%dp_2D_a_01)
      NULLIFY( debug%dp_2D_a_02)
      NULLIFY( debug%dp_2D_a_03)
      NULLIFY( debug%dp_2D_a_04)
      NULLIFY( debug%dp_2D_a_05)
      NULLIFY( debug%dp_2D_a_06)
      NULLIFY( debug%dp_2D_a_07)
      NULLIFY( debug%dp_2D_a_08)
      NULLIFY( debug%dp_2D_a_09)
      NULLIFY( debug%dp_2D_a_10)

      NULLIFY( debug%dp_2D_b_01)
      NULLIFY( debug%dp_2D_b_02)
      NULLIFY( debug%dp_2D_b_03)
      NULLIFY( debug%dp_2D_b_04)
      NULLIFY( debug%dp_2D_b_05)
      NULLIFY( debug%dp_2D_b_06)
      NULLIFY( debug%dp_2D_b_07)
      NULLIFY( debug%dp_2D_b_08)
      NULLIFY( debug%dp_2D_b_09)
      NULLIFY( debug%dp_2D_b_10)

      NULLIFY( debug%dp_2D_c_01)
      NULLIFY( debug%dp_2D_c_02)
      NULLIFY( debug%dp_2D_c_03)
      NULLIFY( debug%dp_2D_c_04)
      NULLIFY( debug%dp_2D_c_05)
      NULLIFY( debug%dp_2D_c_06)
      NULLIFY( debug%dp_2D_c_07)
      NULLIFY( debug%dp_2D_c_08)
      NULLIFY( debug%dp_2D_c_09)
      NULLIFY( debug%dp_2D_c_10)

      NULLIFY( debug%dp_3D_a_01)
      NULLIFY( debug%dp_3D_a_02)
      NULLIFY( debug%dp_3D_a_03)
      NULLIFY( debug%dp_3D_a_04)
      NULLIFY( debug%dp_3D_a_05)
      NULLIFY( debug%dp_3D_a_06)
      NULLIFY( debug%dp_3D_a_07)
      NULLIFY( debug%dp_3D_a_08)
      NULLIFY( debug%dp_3D_a_09)
      NULLIFY( debug%dp_3D_a_10)

      NULLIFY( debug%dp_2D_monthly_a_01)
      NULLIFY( debug%dp_2D_monthly_a_02)
      NULLIFY( debug%dp_2D_monthly_a_03)
      NULLIFY( debug%dp_2D_monthly_a_04)
      NULLIFY( debug%dp_2D_monthly_a_05)
      NULLIFY( debug%dp_2D_monthly_a_06)
      NULLIFY( debug%dp_2D_monthly_a_07)
      NULLIFY( debug%dp_2D_monthly_a_08)
      NULLIFY( debug%dp_2D_monthly_a_09)
      NULLIFY( debug%dp_2D_monthly_a_10)

    END IF

    ! Bind to the actual memory for this region
    IF (region%name == 'NAM') THEN

      debug%int_2D_a_01 => debug_NAM%int_2D_a_01
      debug%int_2D_a_02 => debug_NAM%int_2D_a_02
      debug%int_2D_a_03 => debug_NAM%int_2D_a_03
      debug%int_2D_a_04 => debug_NAM%int_2D_a_04
      debug%int_2D_a_05 => debug_NAM%int_2D_a_05
      debug%int_2D_a_06 => debug_NAM%int_2D_a_06
      debug%int_2D_a_07 => debug_NAM%int_2D_a_07
      debug%int_2D_a_08 => debug_NAM%int_2D_a_08
      debug%int_2D_a_09 => debug_NAM%int_2D_a_09
      debug%int_2D_a_10 => debug_NAM%int_2D_a_10

      debug%int_2D_b_01 => debug_NAM%int_2D_b_01
      debug%int_2D_b_02 => debug_NAM%int_2D_b_02
      debug%int_2D_b_03 => debug_NAM%int_2D_b_03
      debug%int_2D_b_04 => debug_NAM%int_2D_b_04
      debug%int_2D_b_05 => debug_NAM%int_2D_b_05
      debug%int_2D_b_06 => debug_NAM%int_2D_b_06
      debug%int_2D_b_07 => debug_NAM%int_2D_b_07
      debug%int_2D_b_08 => debug_NAM%int_2D_b_08
      debug%int_2D_b_09 => debug_NAM%int_2D_b_09
      debug%int_2D_b_10 => debug_NAM%int_2D_b_10

      debug%int_2D_c_01 => debug_NAM%int_2D_c_01
      debug%int_2D_c_02 => debug_NAM%int_2D_c_02
      debug%int_2D_c_03 => debug_NAM%int_2D_c_03
      debug%int_2D_c_04 => debug_NAM%int_2D_c_04
      debug%int_2D_c_05 => debug_NAM%int_2D_c_05
      debug%int_2D_c_06 => debug_NAM%int_2D_c_06
      debug%int_2D_c_07 => debug_NAM%int_2D_c_07
      debug%int_2D_c_08 => debug_NAM%int_2D_c_08
      debug%int_2D_c_09 => debug_NAM%int_2D_c_09
      debug%int_2D_c_10 => debug_NAM%int_2D_c_10

      debug%dp_2D_a_01 => debug_NAM%dp_2D_a_01
      debug%dp_2D_a_02 => debug_NAM%dp_2D_a_02
      debug%dp_2D_a_03 => debug_NAM%dp_2D_a_03
      debug%dp_2D_a_04 => debug_NAM%dp_2D_a_04
      debug%dp_2D_a_05 => debug_NAM%dp_2D_a_05
      debug%dp_2D_a_06 => debug_NAM%dp_2D_a_06
      debug%dp_2D_a_07 => debug_NAM%dp_2D_a_07
      debug%dp_2D_a_08 => debug_NAM%dp_2D_a_08
      debug%dp_2D_a_09 => debug_NAM%dp_2D_a_09
      debug%dp_2D_a_10 => debug_NAM%dp_2D_a_10

      debug%dp_2D_b_01 => debug_NAM%dp_2D_b_01
      debug%dp_2D_b_02 => debug_NAM%dp_2D_b_02
      debug%dp_2D_b_03 => debug_NAM%dp_2D_b_03
      debug%dp_2D_b_04 => debug_NAM%dp_2D_b_04
      debug%dp_2D_b_05 => debug_NAM%dp_2D_b_05
      debug%dp_2D_b_06 => debug_NAM%dp_2D_b_06
      debug%dp_2D_b_07 => debug_NAM%dp_2D_b_07
      debug%dp_2D_b_08 => debug_NAM%dp_2D_b_08
      debug%dp_2D_b_09 => debug_NAM%dp_2D_b_09
      debug%dp_2D_b_10 => debug_NAM%dp_2D_b_10

      debug%dp_2D_c_01 => debug_NAM%dp_2D_c_01
      debug%dp_2D_c_02 => debug_NAM%dp_2D_c_02
      debug%dp_2D_c_03 => debug_NAM%dp_2D_c_03
      debug%dp_2D_c_04 => debug_NAM%dp_2D_c_04
      debug%dp_2D_c_05 => debug_NAM%dp_2D_c_05
      debug%dp_2D_c_06 => debug_NAM%dp_2D_c_06
      debug%dp_2D_c_07 => debug_NAM%dp_2D_c_07
      debug%dp_2D_c_08 => debug_NAM%dp_2D_c_08
      debug%dp_2D_c_09 => debug_NAM%dp_2D_c_09
      debug%dp_2D_c_10 => debug_NAM%dp_2D_c_10

      debug%dp_3D_a_01 => debug_NAM%dp_3D_a_01
      debug%dp_3D_a_02 => debug_NAM%dp_3D_a_02
      debug%dp_3D_a_03 => debug_NAM%dp_3D_a_03
      debug%dp_3D_a_04 => debug_NAM%dp_3D_a_04
      debug%dp_3D_a_05 => debug_NAM%dp_3D_a_05
      debug%dp_3D_a_06 => debug_NAM%dp_3D_a_06
      debug%dp_3D_a_07 => debug_NAM%dp_3D_a_07
      debug%dp_3D_a_08 => debug_NAM%dp_3D_a_08
      debug%dp_3D_a_09 => debug_NAM%dp_3D_a_09
      debug%dp_3D_a_10 => debug_NAM%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_NAM%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_NAM%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_NAM%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_NAM%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_NAM%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_NAM%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_NAM%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_NAM%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_NAM%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_NAM%dp_2D_monthly_a_10

    ELSEIF (region%name == 'EAS') THEN

      debug%int_2D_a_01 => debug_EAS%int_2D_a_01
      debug%int_2D_a_02 => debug_EAS%int_2D_a_02
      debug%int_2D_a_03 => debug_EAS%int_2D_a_03
      debug%int_2D_a_04 => debug_EAS%int_2D_a_04
      debug%int_2D_a_05 => debug_EAS%int_2D_a_05
      debug%int_2D_a_06 => debug_EAS%int_2D_a_06
      debug%int_2D_a_07 => debug_EAS%int_2D_a_07
      debug%int_2D_a_08 => debug_EAS%int_2D_a_08
      debug%int_2D_a_09 => debug_EAS%int_2D_a_09
      debug%int_2D_a_10 => debug_EAS%int_2D_a_10

      debug%int_2D_b_01 => debug_EAS%int_2D_b_01
      debug%int_2D_b_02 => debug_EAS%int_2D_b_02
      debug%int_2D_b_03 => debug_EAS%int_2D_b_03
      debug%int_2D_b_04 => debug_EAS%int_2D_b_04
      debug%int_2D_b_05 => debug_EAS%int_2D_b_05
      debug%int_2D_b_06 => debug_EAS%int_2D_b_06
      debug%int_2D_b_07 => debug_EAS%int_2D_b_07
      debug%int_2D_b_08 => debug_EAS%int_2D_b_08
      debug%int_2D_b_09 => debug_EAS%int_2D_b_09
      debug%int_2D_b_10 => debug_EAS%int_2D_b_10

      debug%int_2D_c_01 => debug_EAS%int_2D_c_01
      debug%int_2D_c_02 => debug_EAS%int_2D_c_02
      debug%int_2D_c_03 => debug_EAS%int_2D_c_03
      debug%int_2D_c_04 => debug_EAS%int_2D_c_04
      debug%int_2D_c_05 => debug_EAS%int_2D_c_05
      debug%int_2D_c_06 => debug_EAS%int_2D_c_06
      debug%int_2D_c_07 => debug_EAS%int_2D_c_07
      debug%int_2D_c_08 => debug_EAS%int_2D_c_08
      debug%int_2D_c_09 => debug_EAS%int_2D_c_09
      debug%int_2D_c_10 => debug_EAS%int_2D_c_10

      debug%dp_2D_a_01 => debug_EAS%dp_2D_a_01
      debug%dp_2D_a_02 => debug_EAS%dp_2D_a_02
      debug%dp_2D_a_03 => debug_EAS%dp_2D_a_03
      debug%dp_2D_a_04 => debug_EAS%dp_2D_a_04
      debug%dp_2D_a_05 => debug_EAS%dp_2D_a_05
      debug%dp_2D_a_06 => debug_EAS%dp_2D_a_06
      debug%dp_2D_a_07 => debug_EAS%dp_2D_a_07
      debug%dp_2D_a_08 => debug_EAS%dp_2D_a_08
      debug%dp_2D_a_09 => debug_EAS%dp_2D_a_09
      debug%dp_2D_a_10 => debug_EAS%dp_2D_a_10

      debug%dp_2D_b_01 => debug_EAS%dp_2D_b_01
      debug%dp_2D_b_02 => debug_EAS%dp_2D_b_02
      debug%dp_2D_b_03 => debug_EAS%dp_2D_b_03
      debug%dp_2D_b_04 => debug_EAS%dp_2D_b_04
      debug%dp_2D_b_05 => debug_EAS%dp_2D_b_05
      debug%dp_2D_b_06 => debug_EAS%dp_2D_b_06
      debug%dp_2D_b_07 => debug_EAS%dp_2D_b_07
      debug%dp_2D_b_08 => debug_EAS%dp_2D_b_08
      debug%dp_2D_b_09 => debug_EAS%dp_2D_b_09
      debug%dp_2D_b_10 => debug_EAS%dp_2D_b_10

      debug%dp_2D_c_01 => debug_EAS%dp_2D_c_01
      debug%dp_2D_c_02 => debug_EAS%dp_2D_c_02
      debug%dp_2D_c_03 => debug_EAS%dp_2D_c_03
      debug%dp_2D_c_04 => debug_EAS%dp_2D_c_04
      debug%dp_2D_c_05 => debug_EAS%dp_2D_c_05
      debug%dp_2D_c_06 => debug_EAS%dp_2D_c_06
      debug%dp_2D_c_07 => debug_EAS%dp_2D_c_07
      debug%dp_2D_c_08 => debug_EAS%dp_2D_c_08
      debug%dp_2D_c_09 => debug_EAS%dp_2D_c_09
      debug%dp_2D_c_10 => debug_EAS%dp_2D_c_10

      debug%dp_3D_a_01 => debug_EAS%dp_3D_a_01
      debug%dp_3D_a_02 => debug_EAS%dp_3D_a_02
      debug%dp_3D_a_03 => debug_EAS%dp_3D_a_03
      debug%dp_3D_a_04 => debug_EAS%dp_3D_a_04
      debug%dp_3D_a_05 => debug_EAS%dp_3D_a_05
      debug%dp_3D_a_06 => debug_EAS%dp_3D_a_06
      debug%dp_3D_a_07 => debug_EAS%dp_3D_a_07
      debug%dp_3D_a_08 => debug_EAS%dp_3D_a_08
      debug%dp_3D_a_09 => debug_EAS%dp_3D_a_09
      debug%dp_3D_a_10 => debug_EAS%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_EAS%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_EAS%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_EAS%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_EAS%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_EAS%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_EAS%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_EAS%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_EAS%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_EAS%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_EAS%dp_2D_monthly_a_10

    ELSEIF (region%name == 'GRL') THEN

      debug%int_2D_a_01 => debug_GRL%int_2D_a_01
      debug%int_2D_a_02 => debug_GRL%int_2D_a_02
      debug%int_2D_a_03 => debug_GRL%int_2D_a_03
      debug%int_2D_a_04 => debug_GRL%int_2D_a_04
      debug%int_2D_a_05 => debug_GRL%int_2D_a_05
      debug%int_2D_a_06 => debug_GRL%int_2D_a_06
      debug%int_2D_a_07 => debug_GRL%int_2D_a_07
      debug%int_2D_a_08 => debug_GRL%int_2D_a_08
      debug%int_2D_a_09 => debug_GRL%int_2D_a_09
      debug%int_2D_a_10 => debug_GRL%int_2D_a_10

      debug%int_2D_b_01 => debug_GRL%int_2D_b_01
      debug%int_2D_b_02 => debug_GRL%int_2D_b_02
      debug%int_2D_b_03 => debug_GRL%int_2D_b_03
      debug%int_2D_b_04 => debug_GRL%int_2D_b_04
      debug%int_2D_b_05 => debug_GRL%int_2D_b_05
      debug%int_2D_b_06 => debug_GRL%int_2D_b_06
      debug%int_2D_b_07 => debug_GRL%int_2D_b_07
      debug%int_2D_b_08 => debug_GRL%int_2D_b_08
      debug%int_2D_b_09 => debug_GRL%int_2D_b_09
      debug%int_2D_b_10 => debug_GRL%int_2D_b_10

      debug%int_2D_c_01 => debug_GRL%int_2D_c_01
      debug%int_2D_c_02 => debug_GRL%int_2D_c_02
      debug%int_2D_c_03 => debug_GRL%int_2D_c_03
      debug%int_2D_c_04 => debug_GRL%int_2D_c_04
      debug%int_2D_c_05 => debug_GRL%int_2D_c_05
      debug%int_2D_c_06 => debug_GRL%int_2D_c_06
      debug%int_2D_c_07 => debug_GRL%int_2D_c_07
      debug%int_2D_c_08 => debug_GRL%int_2D_c_08
      debug%int_2D_c_09 => debug_GRL%int_2D_c_09
      debug%int_2D_c_10 => debug_GRL%int_2D_c_10

      debug%dp_2D_a_01 => debug_GRL%dp_2D_a_01
      debug%dp_2D_a_02 => debug_GRL%dp_2D_a_02
      debug%dp_2D_a_03 => debug_GRL%dp_2D_a_03
      debug%dp_2D_a_04 => debug_GRL%dp_2D_a_04
      debug%dp_2D_a_05 => debug_GRL%dp_2D_a_05
      debug%dp_2D_a_06 => debug_GRL%dp_2D_a_06
      debug%dp_2D_a_07 => debug_GRL%dp_2D_a_07
      debug%dp_2D_a_08 => debug_GRL%dp_2D_a_08
      debug%dp_2D_a_09 => debug_GRL%dp_2D_a_09
      debug%dp_2D_a_10 => debug_GRL%dp_2D_a_10

      debug%dp_2D_b_01 => debug_GRL%dp_2D_b_01
      debug%dp_2D_b_02 => debug_GRL%dp_2D_b_02
      debug%dp_2D_b_03 => debug_GRL%dp_2D_b_03
      debug%dp_2D_b_04 => debug_GRL%dp_2D_b_04
      debug%dp_2D_b_05 => debug_GRL%dp_2D_b_05
      debug%dp_2D_b_06 => debug_GRL%dp_2D_b_06
      debug%dp_2D_b_07 => debug_GRL%dp_2D_b_07
      debug%dp_2D_b_08 => debug_GRL%dp_2D_b_08
      debug%dp_2D_b_09 => debug_GRL%dp_2D_b_09
      debug%dp_2D_b_10 => debug_GRL%dp_2D_b_10

      debug%dp_2D_c_01 => debug_GRL%dp_2D_c_01
      debug%dp_2D_c_02 => debug_GRL%dp_2D_c_02
      debug%dp_2D_c_03 => debug_GRL%dp_2D_c_03
      debug%dp_2D_c_04 => debug_GRL%dp_2D_c_04
      debug%dp_2D_c_05 => debug_GRL%dp_2D_c_05
      debug%dp_2D_c_06 => debug_GRL%dp_2D_c_06
      debug%dp_2D_c_07 => debug_GRL%dp_2D_c_07
      debug%dp_2D_c_08 => debug_GRL%dp_2D_c_08
      debug%dp_2D_c_09 => debug_GRL%dp_2D_c_09
      debug%dp_2D_c_10 => debug_GRL%dp_2D_c_10

      debug%dp_3D_a_01 => debug_GRL%dp_3D_a_01
      debug%dp_3D_a_02 => debug_GRL%dp_3D_a_02
      debug%dp_3D_a_03 => debug_GRL%dp_3D_a_03
      debug%dp_3D_a_04 => debug_GRL%dp_3D_a_04
      debug%dp_3D_a_05 => debug_GRL%dp_3D_a_05
      debug%dp_3D_a_06 => debug_GRL%dp_3D_a_06
      debug%dp_3D_a_07 => debug_GRL%dp_3D_a_07
      debug%dp_3D_a_08 => debug_GRL%dp_3D_a_08
      debug%dp_3D_a_09 => debug_GRL%dp_3D_a_09
      debug%dp_3D_a_10 => debug_GRL%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_GRL%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_GRL%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_GRL%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_GRL%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_GRL%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_GRL%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_GRL%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_GRL%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_GRL%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_GRL%dp_2D_monthly_a_10

    ELSEIF (region%name == 'ANT') THEN

      debug%int_2D_a_01 => debug_ANT%int_2D_a_01
      debug%int_2D_a_02 => debug_ANT%int_2D_a_02
      debug%int_2D_a_03 => debug_ANT%int_2D_a_03
      debug%int_2D_a_04 => debug_ANT%int_2D_a_04
      debug%int_2D_a_05 => debug_ANT%int_2D_a_05
      debug%int_2D_a_06 => debug_ANT%int_2D_a_06
      debug%int_2D_a_07 => debug_ANT%int_2D_a_07
      debug%int_2D_a_08 => debug_ANT%int_2D_a_08
      debug%int_2D_a_09 => debug_ANT%int_2D_a_09
      debug%int_2D_a_10 => debug_ANT%int_2D_a_10

      debug%int_2D_b_01 => debug_ANT%int_2D_b_01
      debug%int_2D_b_02 => debug_ANT%int_2D_b_02
      debug%int_2D_b_03 => debug_ANT%int_2D_b_03
      debug%int_2D_b_04 => debug_ANT%int_2D_b_04
      debug%int_2D_b_05 => debug_ANT%int_2D_b_05
      debug%int_2D_b_06 => debug_ANT%int_2D_b_06
      debug%int_2D_b_07 => debug_ANT%int_2D_b_07
      debug%int_2D_b_08 => debug_ANT%int_2D_b_08
      debug%int_2D_b_09 => debug_ANT%int_2D_b_09
      debug%int_2D_b_10 => debug_ANT%int_2D_b_10

      debug%int_2D_c_01 => debug_ANT%int_2D_c_01
      debug%int_2D_c_02 => debug_ANT%int_2D_c_02
      debug%int_2D_c_03 => debug_ANT%int_2D_c_03
      debug%int_2D_c_04 => debug_ANT%int_2D_c_04
      debug%int_2D_c_05 => debug_ANT%int_2D_c_05
      debug%int_2D_c_06 => debug_ANT%int_2D_c_06
      debug%int_2D_c_07 => debug_ANT%int_2D_c_07
      debug%int_2D_c_08 => debug_ANT%int_2D_c_08
      debug%int_2D_c_09 => debug_ANT%int_2D_c_09
      debug%int_2D_c_10 => debug_ANT%int_2D_c_10

      debug%dp_2D_a_01 => debug_ANT%dp_2D_a_01
      debug%dp_2D_a_02 => debug_ANT%dp_2D_a_02
      debug%dp_2D_a_03 => debug_ANT%dp_2D_a_03
      debug%dp_2D_a_04 => debug_ANT%dp_2D_a_04
      debug%dp_2D_a_05 => debug_ANT%dp_2D_a_05
      debug%dp_2D_a_06 => debug_ANT%dp_2D_a_06
      debug%dp_2D_a_07 => debug_ANT%dp_2D_a_07
      debug%dp_2D_a_08 => debug_ANT%dp_2D_a_08
      debug%dp_2D_a_09 => debug_ANT%dp_2D_a_09
      debug%dp_2D_a_10 => debug_ANT%dp_2D_a_10

      debug%dp_2D_b_01 => debug_ANT%dp_2D_b_01
      debug%dp_2D_b_02 => debug_ANT%dp_2D_b_02
      debug%dp_2D_b_03 => debug_ANT%dp_2D_b_03
      debug%dp_2D_b_04 => debug_ANT%dp_2D_b_04
      debug%dp_2D_b_05 => debug_ANT%dp_2D_b_05
      debug%dp_2D_b_06 => debug_ANT%dp_2D_b_06
      debug%dp_2D_b_07 => debug_ANT%dp_2D_b_07
      debug%dp_2D_b_08 => debug_ANT%dp_2D_b_08
      debug%dp_2D_b_09 => debug_ANT%dp_2D_b_09
      debug%dp_2D_b_10 => debug_ANT%dp_2D_b_10

      debug%dp_2D_c_01 => debug_ANT%dp_2D_c_01
      debug%dp_2D_c_02 => debug_ANT%dp_2D_c_02
      debug%dp_2D_c_03 => debug_ANT%dp_2D_c_03
      debug%dp_2D_c_04 => debug_ANT%dp_2D_c_04
      debug%dp_2D_c_05 => debug_ANT%dp_2D_c_05
      debug%dp_2D_c_06 => debug_ANT%dp_2D_c_06
      debug%dp_2D_c_07 => debug_ANT%dp_2D_c_07
      debug%dp_2D_c_08 => debug_ANT%dp_2D_c_08
      debug%dp_2D_c_09 => debug_ANT%dp_2D_c_09
      debug%dp_2D_c_10 => debug_ANT%dp_2D_c_10

      debug%dp_3D_a_01 => debug_ANT%dp_3D_a_01
      debug%dp_3D_a_02 => debug_ANT%dp_3D_a_02
      debug%dp_3D_a_03 => debug_ANT%dp_3D_a_03
      debug%dp_3D_a_04 => debug_ANT%dp_3D_a_04
      debug%dp_3D_a_05 => debug_ANT%dp_3D_a_05
      debug%dp_3D_a_06 => debug_ANT%dp_3D_a_06
      debug%dp_3D_a_07 => debug_ANT%dp_3D_a_07
      debug%dp_3D_a_08 => debug_ANT%dp_3D_a_08
      debug%dp_3D_a_09 => debug_ANT%dp_3D_a_09
      debug%dp_3D_a_10 => debug_ANT%dp_3D_a_10

      debug%dp_2D_monthly_a_01 => debug_ANT%dp_2D_monthly_a_01
      debug%dp_2D_monthly_a_02 => debug_ANT%dp_2D_monthly_a_02
      debug%dp_2D_monthly_a_03 => debug_ANT%dp_2D_monthly_a_03
      debug%dp_2D_monthly_a_04 => debug_ANT%dp_2D_monthly_a_04
      debug%dp_2D_monthly_a_05 => debug_ANT%dp_2D_monthly_a_05
      debug%dp_2D_monthly_a_06 => debug_ANT%dp_2D_monthly_a_06
      debug%dp_2D_monthly_a_07 => debug_ANT%dp_2D_monthly_a_07
      debug%dp_2D_monthly_a_08 => debug_ANT%dp_2D_monthly_a_08
      debug%dp_2D_monthly_a_09 => debug_ANT%dp_2D_monthly_a_09
      debug%dp_2D_monthly_a_10 => debug_ANT%dp_2D_monthly_a_10

    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE associate_debug_fields

  SUBROUTINE initialise_debug_fields( region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_debug_fields'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (region%name == 'NAM') THEN
      CALL initialise_debug_fields_region( debug_NAM, region%mesh)
    ELSEIF (region%name == 'EAS') THEN
      CALL initialise_debug_fields_region( debug_EAS, region%mesh)
    ELSEIF (region%name == 'GRL') THEN
      CALL initialise_debug_fields_region( debug_GRL, region%mesh)
    ELSEIF (region%name == 'ANT') THEN
      CALL initialise_debug_fields_region( debug_ANT, region%mesh)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 80)

  END SUBROUTINE initialise_debug_fields

  SUBROUTINE initialise_debug_fields_region( debug, mesh)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_debug_fields),         INTENT(INOUT)     :: debug
    TYPE(type_mesh),                 INTENT(IN)        :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_debug_fields_region'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_01, debug%wint_2D_a_01)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_02, debug%wint_2D_a_02)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_03, debug%wint_2D_a_03)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_04, debug%wint_2D_a_04)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_05, debug%wint_2D_a_05)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_06, debug%wint_2D_a_06)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_07, debug%wint_2D_a_07)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_08, debug%wint_2D_a_08)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_09, debug%wint_2D_a_09)
    CALL allocate_shared_int_1D( mesh%nV, debug%int_2D_a_10, debug%wint_2D_a_10)

    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_01, debug%wint_2D_b_01)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_02, debug%wint_2D_b_02)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_03, debug%wint_2D_b_03)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_04, debug%wint_2D_b_04)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_05, debug%wint_2D_b_05)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_06, debug%wint_2D_b_06)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_07, debug%wint_2D_b_07)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_08, debug%wint_2D_b_08)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_09, debug%wint_2D_b_09)
    CALL allocate_shared_int_1D( mesh%nTri, debug%int_2D_b_10, debug%wint_2D_b_10)

    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_01, debug%wint_2D_c_01)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_02, debug%wint_2D_c_02)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_03, debug%wint_2D_c_03)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_04, debug%wint_2D_c_04)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_05, debug%wint_2D_c_05)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_06, debug%wint_2D_c_06)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_07, debug%wint_2D_c_07)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_08, debug%wint_2D_c_08)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_09, debug%wint_2D_c_09)
    CALL allocate_shared_int_1D( mesh%nAc, debug%int_2D_c_10, debug%wint_2D_c_10)

    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_01, debug%wdp_2D_a_01)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_02, debug%wdp_2D_a_02)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_03, debug%wdp_2D_a_03)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_04, debug%wdp_2D_a_04)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_05, debug%wdp_2D_a_05)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_06, debug%wdp_2D_a_06)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_07, debug%wdp_2D_a_07)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_08, debug%wdp_2D_a_08)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_09, debug%wdp_2D_a_09)
    CALL allocate_shared_dp_1D( mesh%nV, debug%dp_2D_a_10, debug%wdp_2D_a_10)

    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_01, debug%wdp_2D_b_01)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_02, debug%wdp_2D_b_02)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_03, debug%wdp_2D_b_03)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_04, debug%wdp_2D_b_04)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_05, debug%wdp_2D_b_05)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_06, debug%wdp_2D_b_06)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_07, debug%wdp_2D_b_07)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_08, debug%wdp_2D_b_08)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_09, debug%wdp_2D_b_09)
    CALL allocate_shared_dp_1D( mesh%nTri, debug%dp_2D_b_10, debug%wdp_2D_b_10)

    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_01, debug%wdp_2D_c_01)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_02, debug%wdp_2D_c_02)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_03, debug%wdp_2D_c_03)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_04, debug%wdp_2D_c_04)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_05, debug%wdp_2D_c_05)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_06, debug%wdp_2D_c_06)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_07, debug%wdp_2D_c_07)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_08, debug%wdp_2D_c_08)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_09, debug%wdp_2D_c_09)
    CALL allocate_shared_dp_1D( mesh%nAc, debug%dp_2D_c_10, debug%wdp_2D_c_10)

    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_01, debug%wdp_3D_a_01)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_02, debug%wdp_3D_a_02)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_03, debug%wdp_3D_a_03)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_04, debug%wdp_3D_a_04)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_05, debug%wdp_3D_a_05)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_06, debug%wdp_3D_a_06)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_07, debug%wdp_3D_a_07)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_08, debug%wdp_3D_a_08)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_09, debug%wdp_3D_a_09)
    CALL allocate_shared_dp_2D( mesh%nV, C%nz, debug%dp_3D_a_10, debug%wdp_3D_a_10)

    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_01, debug%wdp_2D_monthly_a_01)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_02, debug%wdp_2D_monthly_a_02)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_03, debug%wdp_2D_monthly_a_03)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_04, debug%wdp_2D_monthly_a_04)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_05, debug%wdp_2D_monthly_a_05)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_06, debug%wdp_2D_monthly_a_06)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_07, debug%wdp_2D_monthly_a_07)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_08, debug%wdp_2D_monthly_a_08)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_09, debug%wdp_2D_monthly_a_09)
    CALL allocate_shared_dp_2D( mesh%nV, 12, debug%dp_2D_monthly_a_10, debug%wdp_2D_monthly_a_10)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 80)

  END SUBROUTINE initialise_debug_fields_region

  SUBROUTINE reallocate_debug_fields( region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'reallocate_debug_fields'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (region%name == 'NAM') THEN
      CALL deallocate_debug_fields_region( debug_NAM)
      CALL initialise_debug_fields_region( debug_NAM, region%mesh)
    ELSEIF (region%name == 'EAS') THEN
      CALL deallocate_debug_fields_region( debug_EAS)
      CALL initialise_debug_fields_region( debug_EAS, region%mesh)
    ELSEIF (region%name == 'GRL') THEN
      CALL deallocate_debug_fields_region( debug_GRL)
      CALL initialise_debug_fields_region( debug_GRL, region%mesh)
    ELSEIF (region%name == 'ANT') THEN
      CALL deallocate_debug_fields_region( debug_ANT)
      CALL initialise_debug_fields_region( debug_ANT, region%mesh)
    END IF
    CALL associate_debug_fields( region)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE reallocate_debug_fields

  SUBROUTINE deallocate_debug_fields_region( debug)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_debug_fields),         INTENT(INOUT)     :: debug

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_debug_fields_region'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL deallocate_shared( debug%wint_2D_a_01)
    CALL deallocate_shared( debug%wint_2D_a_02)
    CALL deallocate_shared( debug%wint_2D_a_03)
    CALL deallocate_shared( debug%wint_2D_a_04)
    CALL deallocate_shared( debug%wint_2D_a_05)
    CALL deallocate_shared( debug%wint_2D_a_06)
    CALL deallocate_shared( debug%wint_2D_a_07)
    CALL deallocate_shared( debug%wint_2D_a_08)
    CALL deallocate_shared( debug%wint_2D_a_09)
    CALL deallocate_shared( debug%wint_2D_a_10)

    CALL deallocate_shared( debug%wint_2D_b_01)
    CALL deallocate_shared( debug%wint_2D_b_02)
    CALL deallocate_shared( debug%wint_2D_b_03)
    CALL deallocate_shared( debug%wint_2D_b_04)
    CALL deallocate_shared( debug%wint_2D_b_05)
    CALL deallocate_shared( debug%wint_2D_b_06)
    CALL deallocate_shared( debug%wint_2D_b_07)
    CALL deallocate_shared( debug%wint_2D_b_08)
    CALL deallocate_shared( debug%wint_2D_b_09)
    CALL deallocate_shared( debug%wint_2D_b_10)

    CALL deallocate_shared( debug%wint_2D_c_01)
    CALL deallocate_shared( debug%wint_2D_c_02)
    CALL deallocate_shared( debug%wint_2D_c_03)
    CALL deallocate_shared( debug%wint_2D_c_04)
    CALL deallocate_shared( debug%wint_2D_c_05)
    CALL deallocate_shared( debug%wint_2D_c_06)
    CALL deallocate_shared( debug%wint_2D_c_07)
    CALL deallocate_shared( debug%wint_2D_c_08)
    CALL deallocate_shared( debug%wint_2D_c_09)
    CALL deallocate_shared( debug%wint_2D_c_10)

    CALL deallocate_shared( debug%wdp_2D_a_01)
    CALL deallocate_shared( debug%wdp_2D_a_02)
    CALL deallocate_shared( debug%wdp_2D_a_03)
    CALL deallocate_shared( debug%wdp_2D_a_04)
    CALL deallocate_shared( debug%wdp_2D_a_05)
    CALL deallocate_shared( debug%wdp_2D_a_06)
    CALL deallocate_shared( debug%wdp_2D_a_07)
    CALL deallocate_shared( debug%wdp_2D_a_08)
    CALL deallocate_shared( debug%wdp_2D_a_09)
    CALL deallocate_shared( debug%wdp_2D_a_10)

    CALL deallocate_shared( debug%wdp_2D_b_01)
    CALL deallocate_shared( debug%wdp_2D_b_02)
    CALL deallocate_shared( debug%wdp_2D_b_03)
    CALL deallocate_shared( debug%wdp_2D_b_04)
    CALL deallocate_shared( debug%wdp_2D_b_05)
    CALL deallocate_shared( debug%wdp_2D_b_06)
    CALL deallocate_shared( debug%wdp_2D_b_07)
    CALL deallocate_shared( debug%wdp_2D_b_08)
    CALL deallocate_shared( debug%wdp_2D_b_09)
    CALL deallocate_shared( debug%wdp_2D_b_10)

    CALL deallocate_shared( debug%wdp_2D_c_01)
    CALL deallocate_shared( debug%wdp_2D_c_02)
    CALL deallocate_shared( debug%wdp_2D_c_03)
    CALL deallocate_shared( debug%wdp_2D_c_04)
    CALL deallocate_shared( debug%wdp_2D_c_05)
    CALL deallocate_shared( debug%wdp_2D_c_06)
    CALL deallocate_shared( debug%wdp_2D_c_07)
    CALL deallocate_shared( debug%wdp_2D_c_08)
    CALL deallocate_shared( debug%wdp_2D_c_09)
    CALL deallocate_shared( debug%wdp_2D_c_10)

    CALL deallocate_shared( debug%wdp_3D_a_01)
    CALL deallocate_shared( debug%wdp_3D_a_02)
    CALL deallocate_shared( debug%wdp_3D_a_03)
    CALL deallocate_shared( debug%wdp_3D_a_04)
    CALL deallocate_shared( debug%wdp_3D_a_05)
    CALL deallocate_shared( debug%wdp_3D_a_06)
    CALL deallocate_shared( debug%wdp_3D_a_07)
    CALL deallocate_shared( debug%wdp_3D_a_08)
    CALL deallocate_shared( debug%wdp_3D_a_09)
    CALL deallocate_shared( debug%wdp_3D_a_10)

    CALL deallocate_shared( debug%wdp_2D_monthly_a_01)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_02)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_03)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_04)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_05)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_06)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_07)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_08)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_09)
    CALL deallocate_shared( debug%wdp_2D_monthly_a_10)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE deallocate_debug_fields_region

! ===== Create and write to resource tracking file =====
! ======================================================

  SUBROUTINE write_to_resource_tracking_file( netcdf, time, tcomp_tot)
    ! Write to the resource tracking output file

    USE configuration_module, ONLY: resource_tracker, mem_use_tot_max

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_resource_tracker), INTENT(INOUT) :: netcdf
    REAL(dp),                           INTENT(IN)    :: time
    REAL(dp),                           INTENT(IN)    :: tcomp_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'write_to_resource_tracking_file'
    INTEGER                                           :: i,n
    INTEGER,  DIMENSION(1024)                         :: path_int_enc

    IF (.NOT. par%master) RETURN

    ! Open the file for writing
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Time
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time, time, start = (/netcdf%ti/)))

    ! Actual variables
    ! ================

    ! Total model resource use
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_tcomp_tot, tcomp_tot      , start = (/ netcdf%ti /) ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_mem_tot  , mem_use_tot_max, start = (/ netcdf%ti /) ))

    ! Per-subroutine resource use

    n = SIZE( resource_tracker)

    DO i = 1, n

      ! Subroutine name
      CALL encode_subroutine_path_as_integer( resource_tracker( i)%routine_path, path_int_enc)
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_names( i), path_int_enc ))

      ! Computation time
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_tcomp( i), resource_tracker( i)%tcomp      , start = (/ netcdf%ti /) ))

      ! Memory use (defined as maximum over the preceding coupling interval)
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_mem(   i), resource_tracker( i)%mem_use_max, start = (/ netcdf%ti /) ))

    END DO

    ! Close the file
    CALL close_netcdf_file( netcdf%ncid)

    ! Increase time frame counter
    netcdf%ti = netcdf%ti + 1

  END SUBROUTINE write_to_resource_tracking_file

  SUBROUTINE create_resource_tracking_file( netcdf)
    ! Create the resource tracking output file

    USE configuration_module, ONLY: resource_tracker

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_resource_tracker), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'create_resource_tracking_file'
    LOGICAL                                           :: file_exists
    INTEGER                                           :: t,nl
    INTEGER                                           :: i,n
    CHARACTER(LEN=256)                                :: var_name, long_name

    IF (.NOT. par%master) RETURN

    ! Set time frame index to 1
    netcdf%ti = 1

    ! Create a new file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    netcdf%filename = TRIM(C%output_dir) // '/resource_tracking.nc'
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( netcdf%filename) // '" already exists!')
    END IF

    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error( nf90_create( netcdf%filename, IOR( nf90_clobber, nf90_share), netcdf%ncid))

    ! Define dimensions:
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time       , nf90_unlimited, netcdf%id_dim_time       )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_name_length, 1024          , netcdf%id_dim_name_length)

    ! Placeholders for the dimension ID's, for shorter code
    t  = netcdf%id_dim_time
    nl = netcdf%id_dim_name_length

    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.

    ! Dimension variables: time
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time , [t], netcdf%id_var_time, long_name='Time', units='years'   )

    ! Actual variables
    ! ================

    ! Total model resource use
    CALL create_double_var( netcdf%ncid, 'tcomp_tot', [t], netcdf%id_var_tcomp_tot, long_name='Computation time', units='s'    )
    CALL create_double_var( netcdf%ncid, 'mem_tot'  , [t], netcdf%id_var_mem_tot  , long_name='Memory use'      , units='bytes')

    ! Per-subroutine resource use

    n = SIZE( resource_tracker)

    ALLOCATE( netcdf%id_var_names( n))
    ALLOCATE( netcdf%id_var_tcomp( n))
    ALLOCATE( netcdf%id_var_mem(   n))

    DO i = 1, n

      ! Subroutine name
      ! ===============

      ! Generate variable name (name_00001, name_00002, etc.)
      var_name(  1:256) = ' '
      long_name( 1:256) = ' '
      IF     (i < 10) THEN
        WRITE( var_name ,'(A,I1)') 'name_0000', i
      ELSEIF (i < 100) THEN
        WRITE( var_name,'(A,I2)') 'name_000', i
      ELSEIF (i < 1000) THEN
        WRITE( var_name,'(A,I3)') 'name_00', i
      ELSEIF (i < 10000) THEN
        WRITE( var_name,'(A,I4)') 'name_0', i
      ELSEIF (i < 100000) THEN
        WRITE( var_name,'(A,I5)') 'name_', i
      END IF

      WRITE( long_name,'(A,I1)') 'Full name of subroutine #', i

      ! Create the variable in the NetCDF file
      CALL create_int_var( netcdf%ncid, var_name, [nl], netcdf%id_var_names( i),  long_name = long_name)

      ! Computation time
      ! ================

      ! Generate variable name (tcomp_00001, tcomp_00002, etc.)
      var_name(  1:256) = ' '
      long_name( 1:256) = ' '
      IF     (i < 10) THEN
        WRITE( var_name ,'(A,I1)') 'tcomp_0000', i
      ELSEIF (i < 100) THEN
        WRITE( var_name,'(A,I2)') 'tcomp_000', i
      ELSEIF (i < 1000) THEN
        WRITE( var_name,'(A,I3)') 'tcomp_00', i
      ELSEIF (i < 10000) THEN
        WRITE( var_name,'(A,I4)') 'tcomp_0', i
      ELSEIF (i < 100000) THEN
        WRITE( var_name,'(A,I5)') 'tcomp_', i
      END IF

      WRITE( long_name,'(A,I5)') 'Computation time for subroutine #', i

      ! Create the variable in the NetCDF file
      CALL create_double_var( netcdf%ncid, var_name, [t], netcdf%id_var_tcomp( i),  long_name = long_name, units = 's', missing_value = 0._dp)

      ! Memory use
      ! ==========

      ! Generate variable name (mem_00001, mem_00002, etc.)
      var_name(  1:256) = ' '
      long_name( 1:256) = ' '
      IF     (i < 10) THEN
        WRITE( var_name ,'(A,I1)') 'mem_0000', i
      ELSEIF (i < 100) THEN
        WRITE( var_name,'(A,I2)') 'mem_000', i
      ELSEIF (i < 1000) THEN
        WRITE( var_name,'(A,I3)') 'mem_00', i
      ELSEIF (i < 10000) THEN
        WRITE( var_name,'(A,I4)') 'mem_0', i
      ELSEIF (i < 100000) THEN
        WRITE( var_name,'(A,I5)') 'mem_', i
      END IF

      WRITE( long_name,'(A,I5)') 'Memory use for subroutine #', i

      ! Create the variable in the NetCDF file
      CALL create_double_var( netcdf%ncid, var_name, [t], netcdf%id_var_mem( i),  long_name = long_name, units = 'bytes', missing_value = 0._dp)

    END DO

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( netcdf%ncid))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

  END SUBROUTINE create_resource_tracking_file

  SUBROUTINE encode_subroutine_path_as_integer( subroutine_path, path_int_enc)
    ! Encode the current subroutine path as an integer array so it can be saved as a NetCDF variable
    !
    ! Use the simplest possible encoding:
    !
    !  ' ' = -1 (empty character)
    !
    !    0 = 0
    !    1 = 1
    !    ...
    !    9 = 9
    !
    !    a = 10
    !    b = 11
    !    c = 12
    !    ...
    !    z = 36
    !
    !    A = 37
    !    B = 38
    !    C = 39
    !    ...
    !    Z = 62
    !
    !    _ = 63 (underscore)
    !    / = 64 (forward slash)
    !    ( = 65 (left  bracket)
    !    ) = 66 (right bracket)

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=1024),                INTENT(IN)    :: subroutine_path
    INTEGER,  DIMENSION(1024),          INTENT(OUT)   :: path_int_enc

    ! Local variables:
    INTEGER                                           :: i

    path_int_enc = 0

    DO i = 1, 1024

      SELECT CASE ( subroutine_path( i:i))
      CASE( ' ')
        path_int_enc( i) = -1
      CASE( '0')
        path_int_enc( i) = 0
      CASE( '1')
        path_int_enc( i) = 1
      CASE( '2')
        path_int_enc( i) = 2
      CASE( '3')
        path_int_enc( i) = 3
      CASE( '4')
        path_int_enc( i) = 4
      CASE( '5')
        path_int_enc( i) = 5
      CASE( '6')
        path_int_enc( i) = 6
      CASE( '7')
        path_int_enc( i) = 7
      CASE( '8')
        path_int_enc( i) = 8
      CASE( '9')
        path_int_enc( i) = 9
      CASE( 'a')
        path_int_enc( i) = 11
      CASE( 'b')
        path_int_enc( i) = 12
      CASE( 'c')
        path_int_enc( i) = 13
      CASE( 'd')
        path_int_enc( i) = 14
      CASE( 'e')
        path_int_enc( i) = 15
      CASE( 'f')
        path_int_enc( i) = 16
      CASE( 'g')
        path_int_enc( i) = 17
      CASE( 'h')
        path_int_enc( i) = 18
      CASE( 'i')
        path_int_enc( i) = 19
      CASE( 'j')
        path_int_enc( i) = 20
      CASE( 'k')
        path_int_enc( i) = 21
      CASE( 'l')
        path_int_enc( i) = 22
      CASE( 'm')
        path_int_enc( i) = 23
      CASE( 'n')
        path_int_enc( i) = 24
      CASE( 'o')
        path_int_enc( i) = 25
      CASE( 'p')
        path_int_enc( i) = 26
      CASE( 'q')
        path_int_enc( i) = 27
      CASE( 'r')
        path_int_enc( i) = 28
      CASE( 's')
        path_int_enc( i) = 29
      CASE( 't')
        path_int_enc( i) = 30
      CASE( 'u')
        path_int_enc( i) = 31
      CASE( 'v')
        path_int_enc( i) = 32
      CASE( 'w')
        path_int_enc( i) = 33
      CASE( 'x')
        path_int_enc( i) = 34
      CASE( 'y')
        path_int_enc( i) = 35
      CASE( 'z')
        path_int_enc( i) = 36
      CASE( 'A')
        path_int_enc( i) = 37
      CASE( 'B')
        path_int_enc( i) = 38
      CASE( 'C')
        path_int_enc( i) = 39
      CASE( 'D')
        path_int_enc( i) = 40
      CASE( 'E')
        path_int_enc( i) = 41
      CASE( 'F')
        path_int_enc( i) = 42
      CASE( 'G')
        path_int_enc( i) = 43
      CASE( 'H')
        path_int_enc( i) = 44
      CASE( 'I')
        path_int_enc( i) = 45
      CASE( 'J')
        path_int_enc( i) = 46
      CASE( 'K')
        path_int_enc( i) = 47
      CASE( 'L')
        path_int_enc( i) = 48
      CASE( 'M')
        path_int_enc( i) = 49
      CASE( 'N')
        path_int_enc( i) = 50
      CASE( 'O')
        path_int_enc( i) = 51
      CASE( 'P')
        path_int_enc( i) = 52
      CASE( 'Q')
        path_int_enc( i) = 53
      CASE( 'R')
        path_int_enc( i) = 54
      CASE( 'S')
        path_int_enc( i) = 55
      CASE( 'T')
        path_int_enc( i) = 56
      CASE( 'U')
        path_int_enc( i) = 57
      CASE( 'V')
        path_int_enc( i) = 58
      CASE( 'W')
        path_int_enc( i) = 59
      CASE( 'X')
        path_int_enc( i) = 60
      CASE( 'Y')
        path_int_enc( i) = 61
      CASE( 'Z')
        path_int_enc( i) = 62
      CASE( '_')
        path_int_enc( i) = 63
      CASE( '/')
        path_int_enc( i) = 64
      CASE( '(')
        path_int_enc( i) = 65
      CASE( ')')
        path_int_enc( i) = 66
      CASE DEFAULT
        CALL crash('unknown character in routine_path "' // TRIM( subroutine_path) // '"!')
      END SELECT

    END DO

  END SUBROUTINE encode_subroutine_path_as_integer

END MODULE netcdf_module