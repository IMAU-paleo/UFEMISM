MODULE restart_module
  ! Routines for restarting the model from output of an earlier run.
  ! Read primary mesh data from a NetCDF file, calculate secondary mesh data,
  ! read ice model data, and remap some key data after a mesh update.

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
                                             allocate_shared_dp_1D, allocate_shared_dp_2D
  USE netcdf_module,                   ONLY: inquire_restart_file_mesh, read_restart_file_mesh, &
                                             inquire_restart_file_init, read_restart_file_init
  USE data_types_netcdf_module,        ONLY: type_netcdf_restart
  USE data_types_module,               ONLY: type_model_region, type_mesh, type_restart_data, type_remapping_mesh_mesh
  USE mesh_memory_module,              ONLY: allocate_mesh_primary, allocate_mesh_secondary, deallocate_mesh_all
  USE mesh_help_functions_module,      ONLY: find_Voronoi_cell_areas, get_lat_lon_coordinates, find_triangle_areas, &
                                             find_connection_widths, determine_mesh_resolution, find_POI_xy_coordinates, &
                                             find_POI_vertices_and_weights, find_Voronoi_cell_geometric_centres, check_mesh, &
                                             calc_triangle_geometric_centres
  USE mesh_ArakawaC_module,            ONLY: make_Ac_mesh
  USE mesh_operators_module,           ONLY: calc_matrix_operators_mesh
  USE mesh_creation_module,            ONLY: create_transect
  USE mesh_mapping_module,             ONLY: calc_remapping_operators_mesh_mesh, map_mesh2mesh_2D, &
                                             deallocate_remapping_operators_mesh_mesh

  IMPLICIT NONE

CONTAINS

! ===== Restart mesh =====
! ========================

  SUBROUTINE read_mesh_from_restart_file( mesh, restart, region_name, region_time)
    ! Read primary mesh data from the restart file of a previous run, calculate secondary mesh data.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),               INTENT(INOUT) :: mesh
    TYPE(type_restart_data),       INTENT(INOUT) :: restart
    CHARACTER(LEN=3),              INTENT(IN)    :: region_name
    REAL(dp),                      INTENT(IN)    :: region_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'read_mesh_from_restart_file'

    INTEGER                                      :: nV, nTri, nC_mem
    REAL(dp)                                     :: xmin, xmax, ymin, ymax
    REAL(dp), PARAMETER                          :: tol = 1E-9_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set the filename
    IF     (region_name == 'NAM') THEN
      restart%netcdf%filename = C%filename_restart_NAM
    ELSEIF (region_name == 'EAS') THEN
      restart%netcdf%filename = C%filename_restart_EAS
    ELSEIF (region_name == 'GRL') THEN
      restart%netcdf%filename = C%filename_restart_GRL
    ELSEIF (region_name == 'ANT') THEN
      restart%netcdf%filename = C%filename_restart_ANT
    END IF

    IF (par%master .AND. region_time == C%start_time_of_run) THEN
     WRITE(0,*) '  Reading mesh from restart file "', TRIM(restart%netcdf%filename), '"...'
    END IF

    ! Get the mesh size from the restart file
    IF (par%master) THEN
      CALL inquire_restart_file_mesh( restart%netcdf, nV, nTri, nC_mem)
    END IF

    CALL MPI_BCAST( nV,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( nTri,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( nC_mem, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Allocate memory for primary mesh data
    CALL allocate_mesh_primary( mesh, region_name, nV, nTri, nC_mem)

    IF (par%master) THEN
      mesh%nV   = nV
      mesh%nTri = nTri
    END IF

    ! Read primary mesh data
    IF (par%master) THEN
      CALL read_restart_file_mesh( mesh, restart%netcdf)
    END IF
    CALL sync

    ! Determine vertex and triangle domains
    CALL partition_list( mesh%nV,   par%i, par%n, mesh%vi1, mesh%vi2)
    CALL partition_list( mesh%nTri, par%i, par%n, mesh%ti1, mesh%ti2)

    ! Calculate some mesh metadata
    xmin = MINVAL( mesh%V( mesh%vi1:mesh%vi2,1) )
    xmax = MAXVAL( mesh%V( mesh%vi1:mesh%vi2,1) )
    ymin = MINVAL( mesh%V( mesh%vi1:mesh%vi2,2) )
    ymax = MAXVAL( mesh%V( mesh%vi1:mesh%vi2,2) )

    CALL MPI_REDUCE( xmin, mesh%xmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( xmax, mesh%xmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ymin, mesh%ymin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ymax, mesh%ymax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    IF (par%master) THEN
      mesh%tol_dist = ((mesh%xmax - mesh%xmin) + &
                       (mesh%ymax - mesh%ymin)) * tol / 2._dp
    END IF

    ! Calculate extra mesh data
    CALL allocate_mesh_secondary(                 mesh)    ! Adds  9 MPI windows
    CALL calc_triangle_geometric_centres(         mesh)
    CALL find_Voronoi_cell_areas(                 mesh)
    CALL get_lat_lon_coordinates(                 mesh)
    CALL find_triangle_areas(                     mesh)
    CALL find_connection_widths(                  mesh)
    CALL make_Ac_mesh(                            mesh)    ! Adds  5 MPI windows
    CALL calc_matrix_operators_mesh(              mesh)    ! Adds 42 MPI windows (6 CSR matrices, 7 windows each)
    CALL determine_mesh_resolution(               mesh)
    IF (par%master) THEN
      CALL find_POI_xy_coordinates( mesh)
    END IF
    CALL sync
    CALL find_POI_vertices_and_weights(           mesh)
    CALL find_Voronoi_cell_geometric_centres(     mesh)
    CALL create_transect(                         mesh)

    CALL check_mesh( mesh)

    IF (par%master .AND. region_time == C%start_time_of_run) THEN
      WRITE(0,'(A,I6)')             '    Vertices  : ', mesh%nV
      WRITE(0,'(A,I6)')             '    Triangles : ', mesh%nTri
      WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', mesh%resolution_min/1000._dp, ' - ', &
                                                        mesh%resolution_max/1000._dp, ' km'
      WRITE(0,'(A)')                '   Finished restarting mesh.'
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 110)

  END SUBROUTINE read_mesh_from_restart_file

! ===== Restart data =====
! ========================

  SUBROUTINE read_init_data_from_restart_file( restart, region_name)
    ! Read initial model data from the restart file of a previous run

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_restart_data),       INTENT(INOUT) :: restart
    CHARACTER(LEN=3),              INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'read_init_data_from_restart_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  Reading data from restart file "', TRIM(restart%netcdf%filename), '"...'

    ! Allocate and read restart mesh
    CALL read_mesh_from_restart_file( restart%mesh, restart, region_name, C%start_time_of_run-9e99_dp)

    ! Check if all the data is there
    IF (par%master) CALL inquire_restart_file_init( restart%netcdf)
    CALL sync

    ! Allocate memory
    CALL allocate_shared_dp_1D( restart%mesh%nV, restart%Hi, restart%wHi)
    CALL allocate_shared_dp_1D( restart%mesh%nV, restart%Hb, restart%wHb)
    CALL allocate_shared_dp_1D( restart%mesh%nV, restart%Hs, restart%wHs)

    CALL allocate_shared_dp_1D( restart%mesh%nV, restart%dHi_dt, restart%wdHi_dt)
    CALL allocate_shared_dp_1D( restart%mesh%nV, restart%dHb_dt, restart%wdHb_dt)

    CALL allocate_shared_dp_1D( restart%mesh%nV, restart%dHi_dt_ave, restart%wdHi_dt_ave)

    CALL allocate_shared_dp_2D( restart%mesh%nTri, C%nz, restart%u_3D, restart%wu_3D)
    CALL allocate_shared_dp_2D( restart%mesh%nTri, C%nz, restart%v_3D, restart%wv_3D)

    CALL allocate_shared_dp_1D( restart%mesh%nV, restart%SL,  restart%wSL )
    CALL allocate_shared_dp_1D( restart%mesh%nV, restart%dHb, restart%wdHb)

    CALL allocate_shared_dp_1D( restart%mesh%nV, restart%beta_sq,  restart%wbeta_sq )
    CALL allocate_shared_dp_1D( restart%mesh%nV, restart%phi_fric, restart%wphi_fric)

    IF (C%basal_roughness_restart_type == 'average') THEN
      CALL allocate_shared_dp_1D( restart%mesh%nV, restart%phi_fric_ave, restart%wphi_fric_ave)
    END IF

    CALL allocate_shared_dp_2D( restart%mesh%nV, C%nz, restart%Ti, restart%wTi)

    IF (C%choice_SMB_model == 'IMAU-ITM') THEN

      CALL allocate_shared_dp_1D( restart%mesh%nV,     restart%MeltPreviousYear, restart%wMeltPreviousYear)
      CALL allocate_shared_dp_2D( restart%mesh%nV, 12, restart%FirnDepth,        restart%wFirnDepth       )

      IF (C%do_SMB_IMAUITM_inversion .AND. C%SMB_IMAUITM_inv_choice_init_C == 'restart') THEN
        CALL allocate_shared_dp_1D( restart%mesh%nV, restart%C_abl_constant_inv, restart%wC_abl_constant_inv)
        CALL allocate_shared_dp_1D( restart%mesh%nV, restart%C_abl_Ts_inv,       restart%wC_abl_Ts_inv      )
        CALL allocate_shared_dp_1D( restart%mesh%nV, restart%C_abl_Q_inv,        restart%wC_abl_Q_inv       )
        CALL allocate_shared_dp_1D( restart%mesh%nV, restart%C_refr_inv,         restart%wC_refr_inv        )
      END IF

    END IF

    ! Read data from the restart file
    IF (par%master) CALL read_restart_file_init( region_name, restart, restart%netcdf)
    CALL sync

    ! Deallocate restart mesh
    CALL deallocate_mesh_all( restart%mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 12)

  END SUBROUTINE read_init_data_from_restart_file

! ===== Remapping of restart data =====
! =====================================

  SUBROUTINE remap_restart_data( region)
    ! Remap key data from the restart mesh to a newly updated model mesh

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                  :: routine_name = 'remap_restart_data'
    INTEGER                                        :: int_dummy
    TYPE(type_remapping_mesh_mesh)                 :: map

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. C%do_basal_sliding_inversion) THEN

      ! To prevent compiler warnings for unused variables
      int_dummy = map%int_dummy

      IF (par%master) WRITE(0,*) '   Remapping key restart data...'

      ! Allocate and read restart mesh
      CALL read_mesh_from_restart_file( region%restart%mesh, region%restart, region%name, region%time)

      ! Calculate the mapping array
      CALL calc_remapping_operators_mesh_mesh( region%restart%mesh, region%mesh, map)

      ! Map data field from source mesh to new mesh
      IF (C%basal_roughness_restart_type == 'last') THEN
        CALL map_mesh2mesh_2D( region%restart%mesh, region%mesh, map, &
                               region%restart%phi_fric, region%ice%phi_fric_a, &
                               'nearest_neighbour')

      ELSEIF (C%basal_roughness_restart_type == 'average') THEN
        CALL map_mesh2mesh_2D( region%restart%mesh, region%mesh, map, &
                               region%restart%phi_fric_ave, region%ice%phi_fric_a, &
                               'nearest_neighbour')
      ELSE
        CALL crash('unknown basal_roughness_restart_type "' // TRIM( C%basal_roughness_restart_type) // '"!')
      END IF

      ! Deallocate shared memory for the mapping array
      CALL deallocate_remapping_operators_mesh_mesh( map)

      ! Deallocate restart mesh
      CALL deallocate_mesh_all( region%restart%mesh)

    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_restart_data

END MODULE restart_module