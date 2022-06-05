MODULE restart_module
  ! Routines for restarting the model from output of an earlier run.
  ! Read primary mesh data from a NetCDF file, calculate secondary mesh data,
  ! and read ice model data.

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
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  
  ! Import specific functionality
  USE data_types_netcdf_module,        ONLY: type_netcdf_restart
  USE data_types_module,               ONLY: type_model_region
  USE netcdf_module,                   ONLY: inquire_restart_file_mesh, read_restart_file_mesh, inquire_restart_file_init, &
                                             read_restart_file_init
  USE mesh_memory_module,              ONLY: allocate_mesh_primary, allocate_mesh_secondary
  USE mesh_help_functions_module,      ONLY: find_Voronoi_cell_areas, get_lat_lon_coordinates, find_triangle_areas, &
                                             find_connection_widths, determine_mesh_resolution, find_POI_xy_coordinates, &
                                             find_POI_vertices_and_weights, find_Voronoi_cell_geometric_centres, check_mesh, &
                                             calc_triangle_geometric_centres
  USE mesh_ArakawaC_module,            ONLY: make_Ac_mesh
  USE mesh_operators_module,           ONLY: calc_matrix_operators_mesh
  USE mesh_creation_module,            ONLY: create_transect
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE read_mesh_from_restart_file( region)
    ! Read primary mesh data from the restart file of a previous run, calculate secondary mesh data.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_mesh_from_restart_file'

    INTEGER                                            :: nV, nTri, nC_mem
    REAL(dp)                                           :: xmin, xmax, ymin, ymax
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

   ! Set the filename
   IF     (region%name == 'NAM') THEN
     region%restart%netcdf%filename = C%filename_restart_NAM
   ELSEIF (region%name == 'EAS') THEN
     region%restart%netcdf%filename = C%filename_restart_EAS
   ELSEIF (region%name == 'GRL') THEN
     region%restart%netcdf%filename = C%filename_restart_GRL
   ELSEIF (region%name == 'ANT') THEN
     region%restart%netcdf%filename = C%filename_restart_ANT
   END IF

   IF (par%master) WRITE(0,*) '  Reading mesh from restart file "', TRIM(region%restart%netcdf%filename), '"...'

   ! Get the mesh size from the restart file
   IF (par%master) CALL inquire_restart_file_mesh( region%restart%netcdf, nV, nTri, nC_mem)
   CALL MPI_BCAST( nV,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST( nTri,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST( nC_mem, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


   ! Allocate memory for primary mesh data
   CALL allocate_mesh_primary( region%mesh, region%name, nV, nTri, nC_mem)

   IF (par%master) region%mesh%nV   = nV
   IF (par%master) region%mesh%nTri = nTri

   ! Read primary mesh data
   IF (par%master) CALL read_restart_file_mesh( region%mesh, region%restart%netcdf)
   CALL sync

   ! Determine vertex and triangle domains
   CALL partition_list( region%mesh%nV,   par%i, par%n, region%mesh%vi1, region%mesh%vi2)
   CALL partition_list( region%mesh%nTri, par%i, par%n, region%mesh%ti1, region%mesh%ti2)

   ! Calculate some mesh metadata
   xmin = MINVAL( region%mesh%V( region%mesh%vi1:region%mesh%vi2,1) )
   xmax = MAXVAL( region%mesh%V( region%mesh%vi1:region%mesh%vi2,1) )
   ymin = MINVAL( region%mesh%V( region%mesh%vi1:region%mesh%vi2,2) )
   ymax = MAXVAL( region%mesh%V( region%mesh%vi1:region%mesh%vi2,2) )
   CALL MPI_REDUCE( xmin, region%mesh%xmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_REDUCE( xmax, region%mesh%xmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_REDUCE( ymin, region%mesh%ymin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_REDUCE( ymax, region%mesh%ymax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
   IF (par%master) region%mesh%tol_dist = ((region%mesh%xmax - region%mesh%xmin) + (region%mesh%ymax - region%mesh%ymin)) * tol / 2._dp

   ! Calculate extra mesh data
    CALL allocate_mesh_secondary(                 region%mesh)    ! Adds  9 MPI windows
    CALL calc_triangle_geometric_centres(         region%mesh)
    CALL find_Voronoi_cell_areas(                 region%mesh)
    CALL get_lat_lon_coordinates(                 region%mesh)
    CALL find_triangle_areas(                     region%mesh)
    CALL find_connection_widths(                  region%mesh)
    CALL make_Ac_mesh(                            region%mesh)    ! Adds  5 MPI windows
    CALL calc_matrix_operators_mesh(              region%mesh)    ! Adds 42 MPI windows (6 CSR matrices, 7 windows each)
    CALL determine_mesh_resolution(               region%mesh)
    IF (par%master) CALL find_POI_xy_coordinates( region%mesh)
    CALL sync
    CALL find_POI_vertices_and_weights(           region%mesh)
    CALL find_Voronoi_cell_geometric_centres(     region%mesh)
    CALL create_transect(                         region%mesh)

    CALL check_mesh( region%mesh)

    IF (par%master) THEN
      WRITE(0,'(A)')                '    Finished restarting mesh.'
      WRITE(0,'(A,I6)')             '     Vertices  : ', region%mesh%nV
      WRITE(0,'(A,I6)')             '     Triangles : ', region%mesh%nTri
      WRITE(0,'(A,F7.1,A,F7.1,A)')  '     Resolution: ', region%mesh%resolution_min/1000._dp, ' - ', region%mesh%resolution_max/1000._dp, ' km'
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 110)

  END SUBROUTINE read_mesh_from_restart_file
  SUBROUTINE read_init_data_from_restart_file( region)
    ! Read initial model data from the restart file of a previous run

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_init_data_from_restart_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  Reading initial data from restart file "', TRIM(region%restart%netcdf%filename), '"...'

    ! Check if all the data is there
    IF (par%master) CALL inquire_restart_file_init( region%restart%netcdf)
    CALL sync

    ! Allocate memory
    CALL allocate_shared_dp_1D( region%mesh%nV,       region%restart%Hi,               region%restart%wHi              )
    CALL allocate_shared_dp_1D( region%mesh%nV,       region%restart%Hb,               region%restart%wHb              )
    CALL allocate_shared_dp_1D( region%mesh%nV,       region%restart%Hs,               region%restart%wHs              )
    CALL allocate_shared_dp_2D( region%mesh%nV, C%nZ, region%restart%Ti,               region%restart%wTi              )
    CALL allocate_shared_dp_1D( region%mesh%nV,       region%restart%MeltPreviousYear, region%restart%wMeltPreviousYear)
    CALL allocate_shared_dp_2D( region%mesh%nV, 12,   region%restart%FirnDepth,        region%restart%wFirnDepth       )

    ! Read data from the restart file
    IF (par%master) CALL read_restart_file_init( region%restart, region%restart%netcdf)
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 6)

  END SUBROUTINE read_init_data_from_restart_file

END MODULE restart_module
