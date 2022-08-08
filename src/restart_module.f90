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
                                             allocate_shared_dp_1D, allocate_shared_dp_2D, &
                                             deallocate_shared
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
  USE utilities_module,                ONLY: time_display
  USE ice_dynamics_module,             ONLY: run_ice_model

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
    CALL finalise_routine( routine_name, n_extra_windows_expected = 111)

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
    TYPE(type_remapping_mesh_mesh)                 :: map
    REAL(dp), DIMENSION(:    ), POINTER            ::  bed_method1,  bed_method2
    INTEGER                                        :: wbed_method1, wbed_method2
    INTEGER                                        :: vi, int_dummy

    ! == Initialisation
    ! =================

    ! Add routine to path
    CALL init_routine( routine_name)

    IF ( (.NOT. C%choice_basal_roughness == 'restart') .AND. &
         (.NOT. C%do_use_hi_memory) ) THEN
      CALL finalise_routine(routine_name)
      RETURN
    END IF

    ! To prevent compiler warnings for unused variables
    int_dummy = map%int_dummy

    ! == Mapping basics
    ! =================

    IF (par%master) WRITE(0,*) '   Remapping key restart data...'

    ! Allocate and read restart mesh
    CALL read_mesh_from_restart_file( region%restart%mesh, region%restart, region%name, region%time)

    ! Calculate the mapping array
    CALL calc_remapping_operators_mesh_mesh( region%restart%mesh, region%mesh, map)

    ! == Basal roughness
    ! ==================

    IF (C%choice_basal_roughness == 'restart') THEN

      CALL allocate_shared_dp_1D( region%mesh%nV, bed_method1, wbed_method1)
      CALL allocate_shared_dp_1D( region%mesh%nV, bed_method2, wbed_method2)

      ! Map bed roughness from restart mesh to new mesh
      ! Use a combination of two methods to balance out interpolation errors
      ! Use exp(log()) trick to account for different orders of magnitude
      IF (C%basal_roughness_restart_type == 'last') THEN
        ! Use log of last output from previous run
        ! Use a mapping method that keeps sharpness
        CALL map_mesh2mesh_2D( region%restart%mesh, region%mesh, map, &
                               LOG(region%restart%phi_fric), bed_method1, &
                               'trilin')

        ! Use log of last output from previous run
        ! Use a mapping method that smooths the field
        CALL map_mesh2mesh_2D( region%restart%mesh, region%mesh, map, &
                               LOG(region%restart%phi_fric), bed_method2, &
                               'cons_2nd_order')

      ELSEIF (C%basal_roughness_restart_type == 'average') THEN
        ! Use log of averaged bed roughness from previous run
        ! Use a mapping method that keeps sharpness
        CALL map_mesh2mesh_2D( region%restart%mesh, region%mesh, map, &
                               LOG(region%restart%phi_fric_ave), bed_method1, &
                               'trilin')

        ! Use log of last output from previous run
        ! Use a mapping method that smooths the field
        CALL map_mesh2mesh_2D( region%restart%mesh, region%mesh, map, &
                               LOG(region%restart%phi_fric_ave), bed_method2, &
                               'cons_2nd_order')
      ELSE
        CALL crash('unknown basal_roughness_restart_type "' // TRIM( C%basal_roughness_restart_type) // '"!')
      END IF

      ! Use a combination of two methods to balance out interpolation errors
      ! Use exp(log()) trick to account for different orders of magnitude
      DO vi = region%mesh%vi1, region%mesh%vi2

        ! Keep very high friction values
        IF ( EXP(bed_method1( vi)) >= 20.0_dp .OR. EXP(bed_method2( vi)) >= 20.0_dp ) THEN
          region%ice%phi_fric_a( vi) = EXP( MAX(bed_method1( vi), bed_method2( vi)) )

        ! Keep very high sliding values
        ELSEIF ( EXP(bed_method1( vi)) < 0.2_dp .OR. EXP(bed_method2( vi)) < 0.2_dp ) THEN
          region%ice%phi_fric_a( vi) = EXP( MIN(bed_method1( vi), bed_method2( vi)) )

        ! Compute average of both methods for mid-range values
        ELSE
          region%ice%phi_fric_a( vi) = EXP( (bed_method1( vi) + bed_method2( vi)) / 2._dp )

        END IF

        ! Make sure bed roughness stays within the prescribed limits
        region%ice%phi_fric_a( vi) = MIN(MAX(region%ice%phi_fric_a( vi), &
                                             C%basal_sliding_inv_phi_min), &
                                             C%basal_sliding_inv_phi_max)

      END DO
      CALL sync

      ! Bring values from both methods back from the log domain to the real world
      bed_method1( region%mesh%vi1:region%mesh%vi2) = EXP(bed_method1( region%mesh%vi1:region%mesh%vi2))
      bed_method2( region%mesh%vi1:region%mesh%vi2) = EXP(bed_method2( region%mesh%vi1:region%mesh%vi2))

      ! Calibrate the resulting bed roughness within the range of the two mapping methods
      CALL adjust_remapped_bed_roughness(region, bed_method1, bed_method2)

      ! Clean up after yourself
      CALL deallocate_shared( wbed_method1)
      CALL deallocate_shared( wbed_method2)

    END IF

    ! == Ice thickness rate of change
    ! ===============================

    IF (C%do_use_hi_memory ) THEN

      ! Map data field from source mesh to new mesh
      ! Already re-allocated by remap_ice_model
      CALL map_mesh2mesh_2D( region%restart%mesh, region%mesh, map, &
                             region%restart%dHi_dt_ave, region%ice%dHi_dt_past_a, &
                             'cons_2nd_order')
    END IF

    ! == Finalisation
    ! ===============

    ! Deallocate shared memory for the mapping array
    CALL deallocate_remapping_operators_mesh_mesh( map)

    ! Deallocate restart mesh
    CALL deallocate_mesh_all( region%restart%mesh)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_restart_data

  SUBROUTINE adjust_remapped_bed_roughness( region, bed_method1, bed_method2)
    ! After a bed roughness remap, go back in time a set amount of years and
    ! run only the ice dynamics model (without ice thickness updates) until
    ! coming back to present, adjusting the remapped bed roughness within the
    ! range of values generated by the two mapping methods. The adjustment uses
    ! the averaged past rates of ice thickness change, with the aim of reducing
    ! shock increases/decreases of ice thickness due to interpolation errors.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT)  :: region
    REAL(dp), DIMENSION(region%mesh%nV), INTENT(IN)     :: bed_method1,  bed_method2

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'adjust_remapped_bed_roughness'
    INTEGER                                             :: it, vi
    REAL(dp)                                            :: dt_ave, t_end
    REAL(dp)                                            :: h_scale, h_delta

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%windup_total_years <= 0._dp) THEN
      ! Finalise routine path
      CALL finalise_routine( routine_name)
      ! Exit rutine
      RETURN
    END IF

    ! === Wind up the model ===
    ! =========================

    IF (par%master) THEN
      WRITE (0,*) '    Adjusting the remapped bed roughness field to the new mesh...'
    END IF

    ! Save current time as end-time for wind-up
    t_end = region%time

    ! Bring the timer 1000 years back in time
    region%time = region%time - 1000._dp

    ! Let the model know we want to run velocities from this point on
    region%t_last_SIA       = region%time
    region%t_next_SIA       = region%time
    region%do_SIA           = .TRUE.

    region%t_last_SSA       = region%time
    region%t_next_SSA       = region%time
    region%do_SSA           = .TRUE.

    region%t_last_DIVA      = region%time
    region%t_next_DIVA      = region%time
    region%do_DIVA          = .TRUE.

    ! Run the ice model until coming back to present
    DO WHILE (region%time < t_end)

      ! Make sure memory/senility method is off
      C%do_use_hi_memory = .FALSE.

      ! Calculate ice velocities
      CALL run_ice_model( region, t_end)

      ! Define the dHi_dt factor for scaling of inversion
      h_scale = 1.0_dp/1000._dp

      DO vi = region%mesh%vi1, region%mesh%vi2

        ! Compute difference between the current and reference dHi_dt
        h_delta = region%ice%dHi_dt_a( vi) - region%ice%dHi_dt_ave_a( vi)

        ! Invert only where the model has grounded ice
        IF (region%ice%mask_sheet_a( vi) == 1) THEN

          ! Scale the difference and restrict it to the [-1.5 1.5] range
          h_delta = MAX(-1.5_dp, MIN(1.5_dp, h_delta * h_scale))

          ! Further adjust only where needed
          IF ( (h_delta > 0._dp) .OR. &
               (h_delta < 0._dp) ) THEN

            ! Adjust based on scaled dHi_dt difference
            region%ice%phi_fric_a( vi) = region%ice%phi_fric_a( vi) * (10._dp ** (-h_delta))
            ! Constrain adjusted value to within the range of the two mapping methods
            region%ice%phi_fric_a( vi) = MAX(region%ice%phi_fric_a( vi), MIN(bed_method1( vi), bed_method2( vi)))
            region%ice%phi_fric_a( vi) = MIN(region%ice%phi_fric_a( vi), MAX(bed_method1( vi), bed_method2( vi)))

          END IF

        END IF

        ! Make sure bed roughness stays within the prescribed limits
        region%ice%phi_fric_a( vi) = MIN(MAX(region%ice%phi_fric_a( vi), &
                                             C%basal_sliding_inv_phi_min), &
                                             C%basal_sliding_inv_phi_max)

      END DO
      CALL sync

      ! Advance region time and repeat
      IF (par%master) THEN
        region%time = region%time + region%dt
      END IF
      CALL sync

    END DO

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE adjust_remapped_bed_roughness

END MODULE restart_module
