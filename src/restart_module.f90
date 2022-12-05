MODULE restart_module
  ! Routines for restarting the model from output of an earlier run.
  ! Read primary mesh data from a NetCDF file, calculate secondary mesh data,
  ! and read ice model data.

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

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
  USE netcdf_basic_module,             ONLY: open_existing_netcdf_file_for_reading, close_netcdf_file
  USE netcdf_input_module,             ONLY: setup_mesh_from_file, read_field_from_file_2D, read_field_from_file_2D_monthly, &
                                             read_field_from_file_3D
  USE data_types_netcdf_module,        ONLY: type_netcdf_restart
  USE data_types_module,               ONLY: type_model_region

  IMPLICIT NONE

CONTAINS

! == Mesh
! =======

  SUBROUTINE read_mesh_from_restart_file( region)
    ! Read primary mesh data from the restart file of a previous run, calculate secondary mesh data.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_mesh_from_restart_file'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set the filename
    IF     (region%name == 'NAM') THEN
      filename = C%filename_restart_NAM
    ELSEIF (region%name == 'EAS') THEN
      filename = C%filename_restart_EAS
    ELSEIF (region%name == 'GRL') THEN
      filename = C%filename_restart_GRL
    ELSEIF (region%name == 'ANT') THEN
      filename = C%filename_restart_ANT
    END IF

    IF (par%master) WRITE(0,*) '  Reading mesh from restart file "', TRIM( filename), '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Set up the mesh from the restart file
    CALL setup_mesh_from_file( filename, ncid, region%mesh, region%name)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    IF (par%master) THEN
      WRITE(0,'(A)')                '    Finished restarting mesh.'
      WRITE(0,'(A,I6)')             '     Vertices  : ', region%mesh%nV
      WRITE(0,'(A,I6)')             '     Triangles : ', region%mesh%nTri
      WRITE(0,'(A,F7.1,A,F7.1,A)')  '     Resolution: ', region%mesh%resolution_min/1000._dp, ' - ', region%mesh%resolution_max/1000._dp, ' km'
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 110)

  END SUBROUTINE read_mesh_from_restart_file

! == Data
! =======

  SUBROUTINE read_init_data_from_restart_file( region)
    ! Read initial model data from the restart file of a previous run

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_init_data_from_restart_file'
    CHARACTER(LEN=256)                                 :: filename
    REAL(dp)                                           :: time_to_restart_from
    TYPE(type_netcdf_restart)                          :: restart

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set the filename and time to restart from
    IF     (region%name == 'NAM') THEN
      filename = C%filename_restart_NAM
      time_to_restart_from = C%time_to_restart_from_NAM
    ELSEIF (region%name == 'EAS') THEN
      filename = C%filename_restart_EAS
      time_to_restart_from = C%time_to_restart_from_EAS
    ELSEIF (region%name == 'GRL') THEN
      filename = C%filename_restart_GRL
      time_to_restart_from = C%time_to_restart_from_GRL
    ELSEIF (region%name == 'ANT') THEN
      filename = C%filename_restart_ANT
      time_to_restart_from = C%time_to_restart_from_ANT
    END IF

    IF (par%master) WRITE(0,*) '  Reading data from restart file "', TRIM( filename), '"...'

    ! Allocate memory
    CALL allocate_shared_dp_1D( region%mesh%nV,       region%restart%Hi,                 region%restart%wHi                )
    CALL allocate_shared_dp_1D( region%mesh%nV,       region%restart%Hb,                 region%restart%wHb                )
    CALL allocate_shared_dp_1D( region%mesh%nV,       region%restart%Hs,                 region%restart%wHs                )
    CALL allocate_shared_dp_1D( region%mesh%nV,       region%restart%beta_sq,            region%restart%wbeta_sq           )
    CALL allocate_shared_dp_1D( region%mesh%nV,       region%restart%phi_fric,           region%restart%wphi_fric          )
    CALL allocate_shared_dp_2D( region%mesh%nV, C%nz, region%restart%Ti,                 region%restart%wTi                )
    CALL allocate_shared_dp_1D( region%mesh%nV,       region%restart%MeltPreviousYear,   region%restart%wMeltPreviousYear  )
    CALL allocate_shared_dp_2D( region%mesh%nV, 12,   region%restart%FirnDepth,          region%restart%wFirnDepth         )

    ! Read data from the restart file
    CALL read_field_from_file_2D(         filename, restart%name_var_Hi              , region%mesh, region%restart%Hi              , region%name, time_to_restart_from)
    CALL read_field_from_file_2D(         filename, restart%name_var_Hb              , region%mesh, region%restart%Hb              , region%name, time_to_restart_from)
    CALL read_field_from_file_2D(         filename, restart%name_var_Hs              , region%mesh, region%restart%Hs              , region%name, time_to_restart_from)
!    CALL read_field_from_file_2D(         filename, restart%name_var_beta_sq         , region%mesh, region%restart%beta_sq         , region%name, time_to_restart_from)
!    CALL read_field_from_file_2D(         filename, restart%name_var_phi_fric        , region%mesh, region%restart%phi_fric        , region%name, time_to_restart_from)
    CALL read_field_from_file_3D(         filename, restart%name_var_Ti              , region%mesh, region%restart%Ti              , region%name, time_to_restart_from)
!    CALL read_field_from_file_2D(         filename, restart%name_var_MeltPreviousYear, region%mesh, region%restart%MeltPreviousYear, region%name, time_to_restart_from)
!    CALL read_field_from_file_2D_monthly( filename, restart%name_var_FirnDepth       , region%mesh, region%restart%FirnDepth       , region%name, time_to_restart_from)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 12)

  END SUBROUTINE read_init_data_from_restart_file

END MODULE restart_module