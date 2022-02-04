MODULE netcdf_module

  ! Contains all the subroutines for reading, creating, and writing to NetCDF files.

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
  USE data_types_netcdf_module,      ONLY: type_netcdf_restart, type_netcdf_help_fields
  USE data_types_module,             ONLY: type_model_region, type_mesh, type_grid, type_reference_geometry, type_forcing_data, &
                                           type_subclimate_global, type_debug_fields, type_ICE5G_timeframe, &
                                           type_sparse_matrix_CSR, type_climate_snapshot_global
  USE netcdf,                        ONLY: nf90_max_var_dims, nf90_create, nf90_close, nf90_clobber, nf90_share, nf90_unlimited , &
                                           nf90_enddef, nf90_put_var, nf90_sync, nf90_def_var, nf90_int, nf90_put_att, nf90_def_dim, &
                                           nf90_open, nf90_write, nf90_inq_dimid, nf90_inquire_dimension, nf90_inquire, nf90_double, &
                                           nf90_inq_varid, nf90_inquire_variable, nf90_get_var, nf90_noerr, nf90_strerror, nf90_float
  USE mesh_mapping_module,           ONLY: map_mesh2grid_2D, map_mesh2grid_3D, map_mesh2grid_2D_min
  USE mesh_operators_module,         ONLY: ddx_a_to_a_2D, ddy_a_to_a_2D
  
  IMPLICIT NONE
  
  TYPE(type_debug_fields) :: debug_NAM, debug_EAS, debug_GRL, debug_ANT, debug

CONTAINS

! Main output functions
! =====================

  SUBROUTINE create_output_files( region)
    ! Create a new set of output NetCDF files (restart + help_fields + debug)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    
    IF (.NOT. par%master) RETURN
    
    IF (par%master) WRITE(0,*) '  Creating output files...'
    
    ! Get output file names
    CALL get_output_filenames( region)
    
    ! Create the files
    CALL create_restart_file_mesh(     region, region%restart_mesh)
    CALL create_restart_file_grid(     region, region%restart_grid)
    CALL create_help_fields_file_mesh( region, region%help_fields_mesh)
    CALL create_help_fields_file_grid( region, region%help_fields_grid)
    CALL create_debug_file(            region)
    
  END SUBROUTINE create_output_files
  SUBROUTINE write_to_output_files( region)
    ! Write the current model state to the existing output files
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region
    
    IF (par%master) WRITE(0,'(A,F8.2,A)') '   t = ', region%time/1e3, ' kyr - writing output...'
    
    CALL write_to_restart_file_mesh(     region, region%restart_mesh)
    CALL write_to_restart_file_grid(     region, region%restart_grid)
    CALL write_to_help_fields_file_mesh( region, region%help_fields_mesh)
    CALL write_to_help_fields_file_grid( region, region%help_fields_grid)
    
  END SUBROUTINE write_to_output_files
  SUBROUTINE get_output_filenames( region)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region

    CHARACTER(LEN=256)          :: short_filename
    LOGICAL                     :: ex
    INTEGER                     :: n    
    CHARACTER(LEN=256)          :: ns
    
    IF (.NOT. par%master) RETURN

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

  END SUBROUTINE get_output_filenames

! Create and write to output NetCDF files (mesh versions)
! =======================================================

  SUBROUTINE write_to_restart_file_mesh( region, netcdf)
    ! Write the current model state to the existing output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),   INTENT(INOUT) :: region
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf
    
    IF (.NOT. par%master) RETURN
    
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
    
    ! Temperature
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_Ti,               region%ice%Ti_a,                start = (/ 1, 1,  netcdf%ti/)))
    
    ! SMB
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_FirnDepth,        region%SMB%FirnDepth,           start = (/ 1, 1,  netcdf%ti/)))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_MeltPreviousYear, region%SMB%MeltPreviousYear,    start = (/ 1,     netcdf%ti/)))
    
    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)
    
    ! Increase time frame counter
    netcdf%ti = netcdf%ti + 1
        
  END SUBROUTINE write_to_restart_file_mesh
  SUBROUTINE write_to_help_fields_file_mesh( region, netcdf)
    ! Write the current model state to the existing output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),       INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields), INTENT(INOUT) :: netcdf
    
    IF (.NOT. par%master) RETURN

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
        
  END SUBROUTINE write_to_help_fields_file_mesh
  SUBROUTINE write_help_field_mesh( region, netcdf, id_var, field_name)
    ! Write the current model state to the existing output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(IN)    :: region
    TYPE(type_netcdf_help_fields),  INTENT(IN)    :: netcdf
    INTEGER,                        INTENT(IN)    :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name
    
    IF (field_name == 'none') RETURN
      
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
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%climate%applied%T2m, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'T2m_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%climate%applied%T2m,2)/12._dp, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'Precip') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%climate%applied%Precip, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Precip_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%climate%applied%Precip,2)/12._dp, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%climate%applied%Wind_WE, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%climate%applied%Wind_WE,2)/12._dp, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%climate%applied%Wind_SN, start=(/1, 1, netcdf%ti /) ))
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, SUM(region%climate%applied%Wind_SN,2)/12._dp, start=(/1, netcdf%ti /) ))
      
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
      
!    ! Basal conditions
!    ELSEIF (field_name == 'phi_fric') THEN
!      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%phi_fric_AaAc(1:region%mesh%nV), start=(/1, netcdf%ti /) ))
!    ELSEIF (field_name == 'tau_yield') THEN
!      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%tau_c_AaAc(1:region%mesh%nV), start=(/1, netcdf%ti /) ))
      
    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%IsoIce, start=(/1, netcdf%ti /) ))
    ELSEIF (field_name == 'iso_surf') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%IsoSurf, start=(/1, netcdf%ti /) ))
    
    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%Hb_a - region%refgeo_PD%Hb, start=(/1, netcdf%ti /) ))
    
    ELSE
      WRITE(0,*) ' ERROR: help field "', TRIM(field_name), '" not implemented in write_help_field_mesh!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE write_help_field_mesh
  
  SUBROUTINE create_restart_file_mesh( region, netcdf)
    ! Create a new restart NetCDF file, write the current mesh data to it.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_restart),      INTENT(INOUT) :: netcdf

    ! Local variables:
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, four, vii, time, zeta, month
    
    IF (.NOT. par%master) RETURN
    
    ! Set time frame index to 1
    netcdf%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      WRITE(0,*) 'ERROR: ', TRIM(netcdf%filename), ' already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Create netCDF file
!    WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
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
    CALL create_dim( netcdf%ncid, netcdf%name_dim_four,         4,                       netcdf%id_dim_four        ) ! 4 (each staggered vertex has three "neighbouring" regular vertices)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_vii_transect, region%mesh%nV_transect, netcdf%id_dim_vii_transect) ! Number of vertex pairs in the transect
    
    ! Placeholders for the dimension ID's, for shorter code
    vi        = netcdf%id_dim_vi
    ti        = netcdf%id_dim_ti
    ci        = netcdf%id_dim_ci
    aci       = netcdf%id_dim_aci
    ciplusone = netcdf%id_dim_ciplusone
    two       = netcdf%id_dim_two
    three     = netcdf%id_dim_three
    four      = netcdf%id_dim_four
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
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Aci,              [aci, four ], netcdf%id_var_Aci,              long_name='Staggered to regular vertex indices')
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
    
    ! Temperature
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Ti,               [vi, zeta,  time], netcdf%id_var_Ti,               long_name='Ice temperature', units='K')
    
    ! SMB
    CALL create_double_var( netcdf%ncid, netcdf%name_var_FirnDepth,        [vi, month, time], netcdf%id_var_FirnDepth,        long_name='Firn depth', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_MeltPreviousYear, [vi,        time], netcdf%id_var_MeltPreviousYear, long_name='Melt during previous year', units='mie')
       
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
    
  END SUBROUTINE create_restart_file_mesh
  SUBROUTINE create_help_fields_file_mesh( region, netcdf)
    ! Create a new help_fields NetCDF file, write the current mesh data to it.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf

    ! Local variables:
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, four, vii, time, zeta, month
    
    IF (.NOT. par%master) RETURN
    
    ! Set time frame index to 1
    netcdf%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      WRITE(0,*) 'ERROR: ', TRIM(netcdf%filename), ' already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Create netCDF file
!    WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
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
    CALL create_dim( netcdf%ncid, netcdf%name_dim_four,         4,                       netcdf%id_dim_four        ) ! 4 (each staggered vertex has three "neighbouring" regular vertices)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_vii_transect, region%mesh%nV_transect, netcdf%id_dim_vii_transect) ! Number of vertex pairs in the transect
    
    ! Placeholders for the dimension ID's, for shorter code
    vi        = netcdf%id_dim_vi
    ti        = netcdf%id_dim_ti
    ci        = netcdf%id_dim_ci
    aci       = netcdf%id_dim_aci
    ciplusone = netcdf%id_dim_ciplusone
    two       = netcdf%id_dim_two
    three     = netcdf%id_dim_three
    four      = netcdf%id_dim_four
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
    CALL create_int_var(    netcdf%ncid, netcdf%name_var_Aci,              [aci, four ], netcdf%id_var_Aci,              long_name='Staggered to regular vertex indices')
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
    
  END SUBROUTINE create_help_fields_file_mesh
  SUBROUTINE create_help_field_mesh( netcdf, id_var, field_name)
    ! Add a data field to the help_fields file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf
    INTEGER,                        INTENT(INOUT) :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name
    
    ! Local variables:
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, four, vii, ai, tai, t, z, m
    
    ! Placeholders for the dimension ID's, for shorter code
    vi        = netcdf%id_dim_vi
    ti        = netcdf%id_dim_ti
    ci        = netcdf%id_dim_ci
    aci       = netcdf%id_dim_aci
    ciplusone = netcdf%id_dim_ciplusone
    two       = netcdf%id_dim_two
    three     = netcdf%id_dim_three
    four      = netcdf%id_dim_four
    vii       = netcdf%id_dim_vii_transect
    ai        = netcdf%id_dim_ai
    tai       = netcdf%id_dim_tai
    t         = netcdf%id_dim_time
    z         = netcdf%id_dim_zeta
    m         = netcdf%id_dim_month
    
    IF (field_name == 'none') THEN
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
      
    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL create_int_var(    netcdf%ncid, 'mask',                     [vi,    t], id_var, long_name='mask')
    ELSEIF (field_name == 'mask_land') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_land',                [vi,    t], id_var, long_name='land mask')
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_ocean',               [vi,    t], id_var, long_name='ocean mask')
    ELSEIF (field_name == 'mask_lake') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_lake',                [vi,    t], id_var, long_name='lake mask')
    ELSEIF (field_name == 'mask_ice') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_ice',                 [vi,    t], id_var, long_name='ice mask')
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_sheet',               [vi,    t], id_var, long_name='sheet mask')
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_shelf',               [vi,    t], id_var, long_name='shelf mask')
    ELSEIF (field_name == 'mask_coast') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_coast',               [vi,    t], id_var, long_name='coast mask')
    ELSEIF (field_name == 'mask_margin') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_margin',              [vi,    t], id_var, long_name='margin mask')
    ELSEIF (field_name == 'mask_gl') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_gl',                  [vi,    t], id_var, long_name='grounding-line mask')
    ELSEIF (field_name == 'mask_cf') THEN
      CALL create_int_var(    netcdf%ncid, 'mask_cf',                  [vi,    t], id_var, long_name='calving-front mask')
      
!    ! Basal conditions
!    ELSEIF (field_name == 'phi_fric') THEN
!      CALL create_double_var( netcdf%ncid, 'phi_fric',                 [vi,    t], id_var, long_name='till friction angle', units='degrees')
!    ELSEIF (field_name == 'tau_yield') THEN
!      CALL create_double_var( netcdf%ncid, 'tau_yield',                [vi,    t], id_var, long_name='basal yield stress', units='Pa')
      
    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL create_double_var( netcdf%ncid, 'iso_ice',                  [vi,    t], id_var, long_name='Vertically averaged ice d18O', units='per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL create_double_var( netcdf%ncid, 'iso_surf',                 [vi,    t], id_var, long_name='d18O of precipitation', units='per mille')
    
    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      CALL create_double_var( netcdf%ncid, 'dHb',                      [vi,    t], id_var, long_name='Change in bedrock elevation w.r.t. PD', units='m')
      
    ELSE
      WRITE(0,*) ' ERROR: help field "', TRIM(field_name), '" not implemented in create_help_field_mesh!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE create_help_field_mesh
  
! Create and write to output NetCDF files (grid versions)
! =======================================================

  SUBROUTINE write_to_restart_file_grid( region, netcdf)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),   INTENT(INOUT) :: region
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf
    
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
    
    ! Temperature
    CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Ti_a,             netcdf%id_var_Ti,               netcdf%ti, C%nZ)
    
    ! SMB
    CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%FirnDepth,        netcdf%id_var_FirnDepth,        netcdf%ti, 12  )
    CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%SMB%MeltPreviousYear, netcdf%id_var_MeltPreviousYear, netcdf%ti      )
    
    ! Close the file
    IF (par%master) CALL close_netcdf_file(netcdf%ncid)
    
    ! Increase time frame counter
    IF (par%master) netcdf%ti = netcdf%ti + 1
    CALL sync
        
  END SUBROUTINE write_to_restart_file_grid
  SUBROUTINE write_to_help_fields_file_grid( region, netcdf)
   
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
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),       INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields), INTENT(INOUT) :: netcdf
    INTEGER,                       INTENT(IN)    :: id_var
    CHARACTER(LEN=*),              INTENT(IN)    :: field_name
    
    ! Local variables:
    INTEGER                                       :: vi
    REAL(dp), DIMENSION(:    ), POINTER           ::  dp_2D_a
    INTEGER                                       :: wdp_2D_a
    
    IF (field_name == 'none') RETURN
    
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
      IF (par%master) CALL handle_error( nf90_put_var( netcdf%ncid, id_var, region%ice%GHF_a ))
      
    ! Fields with a time dimension
    ! ============================
    
    ! Mesh
    ELSEIF (field_name == 'resolution') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D_min( netcdf%ncid, region%mesh, region%grid_output, region%mesh%R, id_var, netcdf%ti)
      
    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hi_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Hb') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hb_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Hs') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%Hs_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'SL') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%SL_a, id_var, netcdf%ti)
      
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
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate%applied%T2m, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'T2m_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate%applied%T2m( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Precip') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate%applied%Precip, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Precip_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate%applied%Precip( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate%applied%Wind_WE, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Wind_WE_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate%applied%Wind_WE( vi,:)) / 12._dp
      END DO
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL map_and_write_to_grid_netcdf_dp_3D( netcdf%ncid, region%mesh, region%grid_output, region%climate%applied%Wind_SN, id_var, netcdf%ti, 12)
    ELSEIF (field_name == 'Wind_SN_year') THEN
      DO vi = region%mesh%vi1, region%mesh%vi2
        dp_2D_a( vi) = SUM( region%climate%applied%Wind_SN( vi,:)) / 12._dp
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
      
      ! NOTE: masks commented out; mapping masks between grids is meaningless
      
    ! Masks
    ELSEIF (field_name == 'mask') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask, id_var, netcdf%ti)
    ELSEIF (field_name == 'mask_land') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_land, id_var, netcdf%ti)
    ELSEIF (field_name == 'mask_ocean') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_ocean, id_var, netcdf%ti)
    ELSEIF (field_name == 'mask_lake') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_lake, id_var, netcdf%ti)
    ELSEIF (field_name == 'mask_ice') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_ice, id_var, netcdf%ti)
    ELSEIF (field_name == 'mask_sheet') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_sheet, id_var, netcdf%ti)
    ELSEIF (field_name == 'mask_shelf') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_shelf, id_var, netcdf%ti)
    ELSEIF (field_name == 'mask_coast') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_coast, id_var, netcdf%ti)
    ELSEIF (field_name == 'mask_margin') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_margin, id_var, netcdf%ti)
    ELSEIF (field_name == 'mask_gl') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_gl, id_var, netcdf%ti)
    ELSEIF (field_name == 'mask_cf') THEN
!      CALL map_and_write_to_grid_netcdf_int_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%mask_cf, id_var, netcdf%ti)
      
!    ! Basal conditions
!    ELSEIF (field_name == 'phi_fric') THEN
!      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%phi_fric_AaAc( region%mesh%vi1:region%mesh%vi2)
!      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
!    ELSEIF (field_name == 'tau_yield') THEN
!      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%tau_c_AaAc( region%mesh%vi1:region%mesh%vi2)
!      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
      
    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%IsoIce, id_var, netcdf%ti)
    ELSEIF (field_name == 'iso_surf') THEN
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, region%ice%IsoSurf, id_var, netcdf%ti)
    
    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      dp_2D_a( region%mesh%vi1:region%mesh%vi2) = region%ice%Hb_a( region%mesh%vi1:region%mesh%vi2) - region%refgeo_PD%Hb( region%mesh%vi1:region%mesh%vi2)
      CALL map_and_write_to_grid_netcdf_dp_2D( netcdf%ncid, region%mesh, region%grid_output, dp_2D_a, id_var, netcdf%ti)
    
    ELSE
      WRITE(0,*) ' ERROR: help field "', TRIM(field_name), '" not implemented in write_help_field_grid!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Clean up after yourself
    CALL deallocate_shared( wdp_2D_a)
    
  END SUBROUTINE write_help_field_grid

  SUBROUTINE create_restart_file_grid( region, netcdf)
    ! Create a new restart NetCDF file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_restart),      INTENT(INOUT) :: netcdf

    ! Local variables:
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x,y,z,m,t
    
    IF (.NOT. par%master) RETURN

    ! If the file already exists, return.
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) RETURN
    
    ! Create netCDF file
!    WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
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
    CALL create_double_var( netcdf%ncid, netcdf%name_var_x,                [x            ], netcdf%id_var_x,                long_name='X-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_y,                [   y         ], netcdf%id_var_y,                long_name='Y-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_zeta,             [      z      ], netcdf%id_var_zeta,             long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_month,            [         m   ], netcdf%id_var_month,            long_name='Month', units='1-12'    )
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time,             [            t], netcdf%id_var_time,             long_name='Time', units='years'   )
    
    ! Ice model data
    
    ! Geometry
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hi,               [x, y,       t], netcdf%id_var_Hi,               long_name='Ice thickness', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hb,               [x, y,       t], netcdf%id_var_Hb,               long_name='Bedrock elevation', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Hs,               [x, y,       t], netcdf%id_var_Hs,               long_name='Surface elevation', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_SL,               [x, y,       t], netcdf%id_var_SL,               long_name='Sea surface change', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_dHb,              [x, y,       t], netcdf%id_var_dHb,              long_name='Bedrock deformation', units='m')
    
    ! Temperature
    CALL create_double_var( netcdf%ncid, netcdf%name_var_Ti,               [x, y, z,    t], netcdf%id_var_Ti,               long_name='Ice temperature', units='K')
    
    ! SMB
    CALL create_double_var( netcdf%ncid, netcdf%name_var_FirnDepth,        [x, y,    m, t], netcdf%id_var_FirnDepth,        long_name='Firn depth', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_MeltPreviousYear, [x, y,       t], netcdf%id_var_MeltPreviousYear, long_name='Melt during previous year', units='mie')
       
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
    
  END SUBROUTINE create_restart_file_grid
  SUBROUTINE create_help_fields_file_grid( region, netcdf)
    ! Create a new help fields file, containing secondary model output (not needed for a restart, but interesting to look at)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf

    ! Local variables:
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x, y, z, m, t
    
    IF (.NOT. par%master) RETURN

    ! If the file already exists, return.
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) RETURN
    
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
    
  END SUBROUTINE create_help_fields_file_grid
  SUBROUTINE create_help_field_grid( netcdf, id_var, field_name)
    ! Add a data field to the help_fields file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_netcdf_help_fields),  INTENT(INOUT) :: netcdf
    INTEGER,                        INTENT(INOUT) :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name
    
    ! Local variables:
    INTEGER                                       :: x, y, t, z, m
    
    ! Placeholders for the dimension ID's, for shorter code
    x         = netcdf%id_dim_x
    y         = netcdf%id_dim_y
    t         = netcdf%id_dim_time
    z         = netcdf%id_dim_zeta
    m         = netcdf%id_dim_month
    
    IF (field_name == 'none') THEN
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
      
      
      ! NOTE: masks commented out; mapping masks between grids is meaningless
      
    ! Masks
    ELSEIF (field_name == 'mask') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask',                     [x, y,    t], id_var, long_name='mask')
    ELSEIF (field_name == 'mask_land') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask_land',                [x, y,    t], id_var, long_name='land mask')
    ELSEIF (field_name == 'mask_ocean') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask_ocean',               [x, y,    t], id_var, long_name='ocean mask')
    ELSEIF (field_name == 'mask_lake') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask_lake',                [x, y,    t], id_var, long_name='lake mask')
    ELSEIF (field_name == 'mask_ice') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask_ice',                 [x, y,    t], id_var, long_name='ice mask')
    ELSEIF (field_name == 'mask_sheet') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask_sheet',               [x, y,    t], id_var, long_name='sheet mask')
    ELSEIF (field_name == 'mask_shelf') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask_shelf',               [x, y,    t], id_var, long_name='shelf mask')
    ELSEIF (field_name == 'mask_coast') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask_coast',               [x, y,    t], id_var, long_name='coast mask')
    ELSEIF (field_name == 'mask_margin') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask_margin',              [x, y,    t], id_var, long_name='margin mask')
    ELSEIF (field_name == 'mask_gl') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask_gl',                  [x, y,    t], id_var, long_name='grounding-line mask')
    ELSEIF (field_name == 'mask_cf') THEN
!      CALL create_int_var(    netcdf%ncid, 'mask_cf',                  [x, y,    t], id_var, long_name='calving-front mask')
      
    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL create_double_var( netcdf%ncid, 'phi_fric',                 [x, y,    t], id_var, long_name='till friction angle', units='degrees')
    ELSEIF (field_name == 'tau_yield') THEN
      CALL create_double_var( netcdf%ncid, 'tau_yield',                [x, y,    t], id_var, long_name='basal yield stress', units='Pa')
      
    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL create_double_var( netcdf%ncid, 'iso_ice',                  [x, y,    t], id_var, long_name='Vertically averaged ice d18O', units='per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL create_double_var( netcdf%ncid, 'iso_surf',                 [x, y,    t], id_var, long_name='d18O of precipitation', units='per mille')
    
    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      CALL create_double_var( netcdf%ncid, 'dHb',                      [x, y,    t], id_var, long_name='Change in bedrock elevation w.r.t. PD', units='m')
      
    ELSE
      WRITE(0,*) ' ERROR: help field "', TRIM(field_name), '" not implemented in create_help_field_grid!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE create_help_field_grid
  
  ! Map a model data field from the model mesh to the output grid, and write it to a NetCDF file.
  SUBROUTINE map_and_write_to_grid_netcdf_int_2D( ncid, mesh, grid, d_mesh, id_var, ti)
   
    IMPLICIT NONE
    
    ! Input variables:
    INTEGER,                    INTENT(IN)        :: ncid
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_grid),            INTENT(IN)        :: grid
    INTEGER,  DIMENSION(:    ), INTENT(IN)        :: d_mesh
    INTEGER,                    INTENT(IN)        :: id_var
    INTEGER,                    INTENT(IN)        :: ti
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER           :: d_mesh_dp
    REAL(dp), DIMENSION(:,:  ), POINTER           :: d_grid
    INTEGER                                       :: wd_mesh_dp, wd_grid
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV,          d_mesh_dp, wd_mesh_dp)
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid,    wd_grid   )
    
    ! Convert from integer to real
    d_mesh_dp( mesh%vi1:mesh%vi2) = REAL( d_mesh( mesh%vi1:mesh%vi2), dp)
    
    ! Map data from the model mesh to the square grid
    CALL map_mesh2grid_2D( mesh, grid, d_mesh_dp, d_grid)
    
    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, ti/) ))
    
    ! Deallocate shared memory
    CALL deallocate_shared( wd_mesh_dp)
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE map_and_write_to_grid_netcdf_int_2D
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
    REAL(dp), DIMENSION(:,:  ), POINTER           :: d_grid
    INTEGER                                       :: wd_grid
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid, wd_grid)
    
    ! Map data from the model mesh to the square grid
    CALL map_mesh2grid_2D( mesh, grid, d_mesh, d_grid)
    
    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, ti/) ))
    
    ! Deallocate shared memory
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE map_and_write_to_grid_netcdf_dp_2D
  SUBROUTINE map_and_write_to_grid_netcdf_dp_2D_min(  ncid, mesh, grid, d_mesh, id_var, ti)
   
    IMPLICIT NONE
    
    ! Input variables:
    INTEGER,                    INTENT(IN)        :: ncid
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_grid),            INTENT(IN)        :: grid
    REAL(dp), DIMENSION(:    ), INTENT(IN)        :: d_mesh
    INTEGER,                    INTENT(IN)        :: id_var
    INTEGER,                    INTENT(IN)        :: ti
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER           :: d_grid
    INTEGER                                       :: wd_grid
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid, wd_grid)
    
    ! Map data from the model mesh to the square grid
    CALL map_mesh2grid_2D_min( mesh, grid, d_mesh, d_grid)
    
    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, ti/) ))
    
    ! Deallocate shared memory
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE map_and_write_to_grid_netcdf_dp_2D_min
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
    REAL(dp), DIMENSION(:,:,:), POINTER           :: d_grid
    INTEGER                                       :: wd_grid
    
    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, nz, d_grid, wd_grid)
    
    ! Map data from the model mesh to the square grid
    CALL map_mesh2grid_3D( mesh, grid, d_mesh, d_grid)
    
    ! Write grid data to NetCDF
    IF (par%master) CALL handle_error( nf90_put_var( ncid, id_var, d_grid, start=(/1, 1, 1, ti/) ))
    
    ! Deallocate shared memory
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE map_and_write_to_grid_netcdf_dp_3D
  
! Create and write to debug NetCDF file
! =====================================

  SUBROUTINE write_to_debug_file
    ! Write the current set of debug data fields to the debug NetCDF file
   
    IMPLICIT NONE
    
    ! Local variables:
    INTEGER                                :: ncid
    
    IF (.NOT. par%master) RETURN
    
    IF (.NOT. C%do_write_debug_data) RETURN
    
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
        
  END SUBROUTINE write_to_debug_file
  SUBROUTINE create_debug_file( region)
    ! Create the debug NetCDF file; a lot of data fields but no time dimension.
    
    USE data_types_netcdf_module, ONLY: type_netcdf_debug
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    TYPE(type_netcdf_debug)                       :: debug_temp
    CHARACTER(LEN=20)                             :: short_filename
    INTEGER                                       :: n
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, aci, ciplusone, two, three, four, vii, zeta, month
    
    IF (.NOT. par%master) RETURN
    
    IF (.NOT. C%do_write_debug_data) RETURN

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
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_four,         4,                       debug_temp%id_dim_four        ) ! 4 (each staggered vertex has three "neighbouring" regular vertices)
    CALL create_dim( debug_temp%ncid, debug_temp%name_dim_vii_transect, region%mesh%nV_transect, debug_temp%id_dim_vii_transect) ! Number of vertex pairs in the transect
    
    ! Placeholders for the dimension ID's, for shorter code
    vi        = debug_temp%id_dim_vi
    ti        = debug_temp%id_dim_ti
    ci        = debug_temp%id_dim_ci
    aci       = debug_temp%id_dim_aci
    ciplusone = debug_temp%id_dim_ciplusone
    two       = debug_temp%id_dim_two
    three     = debug_temp%id_dim_three
    four      = debug_temp%id_dim_four
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
    CALL create_int_var(    debug_temp%ncid, debug_temp%name_var_Aci,              [aci, four ], debug_temp%id_var_Aci,              long_name='Staggered to regular vertex indices')
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
    
  END SUBROUTINE associate_debug_fields
  SUBROUTINE initialise_debug_fields( region)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    ! Local variables:
    CHARACTER(LEN=64), PARAMETER                  :: routine_name = 'initialise_debug_fields'
    INTEGER                                       :: n1, n2
    
    n1 = par%mem%n
    
    IF     (region%name == 'NAM') THEN
      CALL initialise_debug_fields_region( debug_NAM, region%mesh)
    ELSEIF (region%name == 'EAS') THEN
      CALL initialise_debug_fields_region( debug_EAS, region%mesh)
    ELSEIF (region%name == 'GRL') THEN
      CALL initialise_debug_fields_region( debug_GRL, region%mesh)
    ELSEIF (region%name == 'ANT') THEN
      CALL initialise_debug_fields_region( debug_ANT, region%mesh)
    END IF
    
    n2 = par%mem%n
    CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE initialise_debug_fields
  SUBROUTINE initialise_debug_fields_region( debug, mesh)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_debug_fields),         INTENT(INOUT)     :: debug
    TYPE(type_mesh),                 INTENT(IN)        :: mesh
    
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
    
  END SUBROUTINE initialise_debug_fields_region
  SUBROUTINE reallocate_debug_fields( region)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
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
    
  END SUBROUTINE reallocate_debug_fields
  SUBROUTINE deallocate_debug_fields_region( debug)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_debug_fields),         INTENT(INOUT)     :: debug
    
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
    
  END SUBROUTINE deallocate_debug_fields_region
  
! Read all kinds of input files
! =============================

  ! A restart file produced by an earlier run
  SUBROUTINE inquire_restart_file_mesh( netcdf, nV, nTri, nC_mem)
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf
    INTEGER,                   INTENT(OUT)   :: nV, nTri, nC_mem
    
    ! Local variables:
!    INTEGER                                  :: int_dummy
    
    WRITE(0,*) 'inquire_restart_file_mesh - FIXME!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        
!    ! Open the netcdf file
!    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)
!    
!    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
!    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_vi,    nV,        netcdf%id_dim_vi   )
!    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_ti,    nTri,      netcdf%id_dim_ti   )
!    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_ci,    nC_mem,    netcdf%id_dim_ci   )
!    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_two,   int_dummy, netcdf%id_dim_two  )
!    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_three, int_dummy, netcdf%id_dim_three)
!
!    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_V,              (/ netcdf%id_dim_vi, netcdf%id_dim_two  /), netcdf%id_var_V             )
!    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_nC,             (/ netcdf%id_dim_vi                     /), netcdf%id_var_nC            )
!    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_C,              (/ netcdf%id_dim_vi, netcdf%id_dim_ci   /), netcdf%id_var_C             )
!    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_niTri,          (/ netcdf%id_dim_vi                     /), netcdf%id_var_niTri         )
!    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_iTri,           (/ netcdf%id_dim_vi, netcdf%id_dim_ci   /), netcdf%id_var_iTri          )
!    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_edge_index,     (/ netcdf%id_dim_vi                     /), netcdf%id_var_edge_index    )
!    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_Tri,            (/ netcdf%id_dim_ti, netcdf%id_dim_three/), netcdf%id_var_Tri           )
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Tricc,          (/ netcdf%id_dim_ti, netcdf%id_dim_two  /), netcdf%id_var_Tricc         )
!    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_TriC,           (/ netcdf%id_dim_ti, netcdf%id_dim_three/), netcdf%id_var_TriC          )
!    CALL inquire_int_var(    netcdf%ncid, netcdf%name_var_Tri_edge_index, (/ netcdf%id_dim_ti                     /), netcdf%id_var_Tri_edge_index)
!        
!    ! Close the netcdf file
!    CALL close_netcdf_file( netcdf%ncid)
    
  END SUBROUTINE inquire_restart_file_mesh
  SUBROUTINE inquire_restart_file_init( netcdf)
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf
    
    ! Local variables:
!    INTEGER                                  :: nZ, nt, nm, k
!    REAL(dp), DIMENSION(:    ), ALLOCATABLE  :: zeta
    
    WRITE(0,*) 'inquire_restart_file_init - FIXME!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        
!    ! Open the netcdf file
!    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)
!    
!    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
!    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_zeta,  nZ, netcdf%id_dim_zeta )
!    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_time,  nt, netcdf%id_dim_time )
!    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_month, nm, netcdf%id_dim_month)
!    
!    IF (nZ /= C%nZ) THEN
!      WRITE(0,*) '   ERROR: nZ in restart file doesnt match nZ in config!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!
!    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_zeta,             (/ netcdf%id_dim_zeta                    /), netcdf%id_var_zeta            )
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_time,             (/ netcdf%id_dim_time                    /), netcdf%id_var_time            )
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_month,            (/ netcdf%id_dim_month                   /), netcdf%id_var_month           )
!    
!    ! Read zeta, check if it matches the config zeta
!    ALLOCATE( zeta( nZ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_zeta, zeta, start = (/ 1 /) ))
!    DO k = 1, C%nz
!      IF (ABS(C%zeta(k) - zeta(k)) > 0.0001_dp) THEN
!        WRITE(0,*) '  WARNING - vertical coordinate zeta in restart file doesnt match zeta in config!'
!      END IF
!      END DO
!    DEALLOCATE( zeta)
!    
!    ! Inquire model data
!    
!    ! Geometry
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Hi,               (/ netcdf%id_dim_vi,                      netcdf%id_dim_time /), netcdf%id_var_Hi              )
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Hb,               (/ netcdf%id_dim_vi,                      netcdf%id_dim_time /), netcdf%id_var_Hb              )
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Hs,               (/ netcdf%id_dim_vi,                      netcdf%id_dim_time /), netcdf%id_var_Hs              )
!    
!    ! Temperature
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Ti,               (/ netcdf%id_dim_vi, netcdf%id_dim_zeta,  netcdf%id_dim_time /), netcdf%id_var_Ti              )
!    
!    ! SMB
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_MeltPreviousYear, (/ netcdf%id_dim_vi,                      netcdf%id_dim_time /), netcdf%id_var_MeltPreviousYear)
!    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_FirnDepth,        (/ netcdf%id_dim_vi, netcdf%id_dim_month, netcdf%id_dim_time /), netcdf%id_var_FirnDepth       )
!        
!    ! Close the netcdf file
!    CALL close_netcdf_file( netcdf%ncid)
    
  END SUBROUTINE inquire_restart_file_init
  SUBROUTINE read_restart_file_mesh( mesh, netcdf)
    ! Read mesh data from a restart file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),           INTENT(INOUT) :: mesh
    TYPE(type_netcdf_restart), INTENT(INOUT) :: netcdf
    
    WRITE(0,*) 'read_restart_file_mesh - FIXME!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    
!    ! Open the netcdf file
!    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)
!    
!    ! Read the data
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_V,              mesh%V,              start = (/ 1, 1 /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_nC,             mesh%nC,             start = (/ 1    /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_C,              mesh%C,              start = (/ 1, 1 /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_niTri,          mesh%niTri,          start = (/ 1    /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_iTri,           mesh%iTri,           start = (/ 1, 1 /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_edge_index,     mesh%edge_index,     start = (/ 1    /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Tri,            mesh%Tri,            start = (/ 1, 1 /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Tricc,          mesh%Tricc,          start = (/ 1, 1 /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_TriC,           mesh%TriC,           start = (/ 1, 1 /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Tri_edge_index, mesh%Tri_edge_index, start = (/ 1    /) ))
!        
!    ! Close the netcdf file
!    CALL close_netcdf_file( netcdf%ncid)
    
  END SUBROUTINE read_restart_file_mesh
  SUBROUTINE read_restart_file_init( refgeo_init, netcdf)
    ! Read mesh data from a restart file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_reference_geometry), INTENT(INOUT) :: refgeo_init
    TYPE(type_netcdf_restart),   INTENT(INOUT) :: netcdf
    
    ! Local variables:
!    INTEGER                                    :: nt, ti, ti_min
!    REAL(dp), DIMENSION(:    ), ALLOCATABLE    :: time
!    REAL(dp)                                   :: dt_min, dt
    
    WRITE(0,*) 'read_restart_file_init - FIXME!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    
!    ! Open the netcdf file
!    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)
!    
!    ! Read time, determine which time frame to read
!    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_time,  nt, netcdf%id_dim_time )
!    ALLOCATE( time( nt))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_time, time, start = (/ 1 /) ))
!    
!    IF (C%time_to_restart_from < MINVAL(time) .OR. C%time_to_restart_from > MAXVAL(time)) THEN
!      WRITE(0,*) '  ERROR - time_to_restart_from ', C%time_to_restart_from, ' outside range of restart file! (range = [', MINVAL( time), ' - ', MAXVAL(time), '])'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF    
!    
!    ti_min = 0
!    dt_min = 1E8_dp
!    DO ti = 1, nt
!      dt = ABS(time( ti) - C%time_to_restart_from)
!      IF (dt < dt_min) THEN
!        ti_min = ti
!        dt_min = dt
!      END IF
!    END DO
!    ti = ti_min
!    
!    DEALLOCATE( time)
!    
!    ! Read the data
!    
!    ! Geometry
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Hi,               init%Hi,               start = (/ 1,    ti /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Hb,               init%Hb,               start = (/ 1,    ti /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Hs,               init%Hs,               start = (/ 1,    ti /) ))
!    
!    ! Temperature
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Ti,               init%Ti,               start = (/ 1, 1, ti /) ))
!    
!    ! SMB
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_MeltPreviousYear, init%MeltPreviousYear, start = (/ 1,    ti /) ))
!    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_FirnDepth,        init%FirnDepth,        start = (/ 1, 1, ti /) ))
!        
!    ! Close the netcdf file
!    CALL close_netcdf_file( netcdf%ncid)
    
  END SUBROUTINE read_restart_file_init 
  
  ! Reference ice-sheet geometry (ice thickness, bed topography, and surface elevation)
  SUBROUTINE inquire_reference_geometry_file( refgeo)
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_reference_geometry), INTENT(INOUT) :: refgeo
        
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
    
  END SUBROUTINE inquire_reference_geometry_file
  SUBROUTINE read_reference_geometry_file(    refgeo)
    ! Read reference geometry data from a NetCDF file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_reference_geometry), INTENT(INOUT) :: refgeo
    
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
    
  END SUBROUTINE read_reference_geometry_file
  
  ! Present-day observed global climate (e.g. ERA-40)
  SUBROUTINE inquire_PD_obs_data_file( PD_obs) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: PD_obs
 
    ! Local variables:
    INTEGER                               :: int_dummy
        
    ! Open the netcdf file
    CALL open_netcdf_file(PD_obs%netcdf%filename, PD_obs%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_dim_lat,     PD_obs%grid%nlat,  PD_obs%netcdf%id_dim_lat)
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_dim_lon,     PD_obs%grid%nlon,  PD_obs%netcdf%id_dim_lon)
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_dim_month,   int_dummy,         PD_obs%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_lat,      (/ PD_obs%netcdf%id_dim_lat                                                       /),  PD_obs%netcdf%id_var_lat)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_lon,      (/ PD_obs%netcdf%id_dim_lon                                                       /),  PD_obs%netcdf%id_var_lon)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_T2m,      (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_T2m)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Precip,   (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_Precip)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Hs,       (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat                             /),  PD_obs%netcdf%id_var_Hs)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Wind_WE,  (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_Wind_WE)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Wind_SN,  (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_Wind_SN)
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD_obs%netcdf%ncid)
    
  END SUBROUTINE inquire_PD_obs_data_file
  SUBROUTINE read_PD_obs_data_file(    PD_obs)
    ! Read the PD_obs0 netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: PD_obs
    
    ! Open the netcdf file
    CALL open_netcdf_file(PD_obs%netcdf%filename, PD_obs%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_lon,     PD_obs%grid%lon, start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_lat,     PD_obs%grid%lat, start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_T2m,     PD_obs%T2m,      start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Precip,  PD_obs%Precip,   start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Hs,      PD_obs%Hs_ref,   start = (/ 1, 1    /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Wind_WE, PD_obs%Wind_WE,  start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Wind_SN, PD_obs%Wind_SN,  start = (/ 1, 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD_obs%netcdf%ncid)
    
  END SUBROUTINE read_PD_obs_data_file
  
  ! GCM global climate (climate matrix snapshots)
  SUBROUTINE inquire_GCM_snapshot( snapshot) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: snapshot
 
    ! Local variables:
    INTEGER                               :: int_dummy
        
    ! Open the netcdf file
    CALL open_netcdf_file( snapshot%netcdf%filename, snapshot%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( snapshot%netcdf%ncid, snapshot%netcdf%name_dim_lat,     snapshot%grid%nlat,  snapshot%netcdf%id_dim_lat)
    CALL inquire_dim( snapshot%netcdf%ncid, snapshot%netcdf%name_dim_lon,     snapshot%grid%nlon,  snapshot%netcdf%id_dim_lon)
    CALL inquire_dim( snapshot%netcdf%ncid, snapshot%netcdf%name_dim_month,   int_dummy,           snapshot%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_lat,      (/ snapshot%netcdf%id_dim_lat                                                           /),  snapshot%netcdf%id_var_lat)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_lon,      (/ snapshot%netcdf%id_dim_lon                                                           /),  snapshot%netcdf%id_var_lon)
    !CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Hi,       (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat                               /),  snapshot%netcdf%id_var_Hi)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Hs,       (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat                               /),  snapshot%netcdf%id_var_Hs)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_T2m,      (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_T2m)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Precip,   (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_Precip)
    !CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Wind_WE,  (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_Wind_WE)
    !CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Wind_SN,  (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_Wind_SN)
        
    ! Close the netcdf file
    CALL close_netcdf_file(snapshot%netcdf%ncid)
    
  END SUBROUTINE inquire_GCM_snapshot
  SUBROUTINE read_GCM_snapshot(    snapshot)
    ! Read the PD_obs0 netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: snapshot
    
    ! Open the netcdf file
    CALL open_netcdf_file(snapshot%netcdf%filename, snapshot%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_lon,     snapshot%grid%lon, start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_lat,     snapshot%grid%lat, start = (/ 1       /) ))
   !CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Hi,      snapshot%Hi,       start = (/ 1, 1    /) ))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Hs,      snapshot%Hs_ref,   start = (/ 1, 1    /) ))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_T2m,     snapshot%T2m,      start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Precip,  snapshot%Precip,   start = (/ 1, 1, 1 /) ))
   !CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Wind_WE, snapshot%Wind_WE,  start = (/ 1, 1, 1 /) ))
   !CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Wind_SN, snapshot%Wind_SN,  start = (/ 1, 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(snapshot%netcdf%ncid)
    
  END SUBROUTINE read_GCM_snapshot
  
  ! ICE5G ice geometry (needed for GCM snapshots)
  SUBROUTINE inquire_ICE5G_data( ICE5G)
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_ICE5G_timeframe), INTENT(INOUT) :: ICE5G
        
    ! Open the netcdf file
    CALL open_netcdf_file( ICE5G%netcdf%filename, ICE5G%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( ICE5G%netcdf%ncid, ICE5G%netcdf%name_dim_lat, ICE5G%grid%nlat, ICE5G%netcdf%id_dim_lat)
    CALL inquire_dim( ICE5G%netcdf%ncid, ICE5G%netcdf%name_dim_lon, ICE5G%grid%nlon, ICE5G%netcdf%id_dim_lon)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_single_var( ICE5G%netcdf%ncid, ICE5G%netcdf%name_var_lat,      (/ ICE5G%netcdf%id_dim_lat                         /),  ICE5G%netcdf%id_var_lat     )
    CALL inquire_single_var( ICE5G%netcdf%ncid, ICE5G%netcdf%name_var_lon,      (/ ICE5G%netcdf%id_dim_lon                         /),  ICE5G%netcdf%id_var_lon     )
    CALL inquire_single_var( ICE5G%netcdf%ncid, ICE5G%netcdf%name_var_Hi,       (/ ICE5G%netcdf%id_dim_lon, ICE5G%netcdf%id_dim_lat/),  ICE5G%netcdf%id_var_Hi      )
    CALL inquire_single_var( ICE5G%netcdf%ncid, ICE5G%netcdf%name_var_Hb,       (/ ICE5G%netcdf%id_dim_lon, ICE5G%netcdf%id_dim_lat/),  ICE5G%netcdf%id_var_Hb      )
    CALL inquire_single_var( ICE5G%netcdf%ncid, ICE5G%netcdf%name_var_mask_ice, (/ ICE5G%netcdf%id_dim_lon, ICE5G%netcdf%id_dim_lat/),  ICE5G%netcdf%id_var_mask_ice)
        
    ! Close the netcdf file
    CALL close_netcdf_file(ICE5G%netcdf%ncid)
    
  END SUBROUTINE inquire_ICE5G_data
  SUBROUTINE read_ICE5G_data(    ICE5G)
    ! Read the PD_obs0 netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_ICE5G_timeframe), INTENT(INOUT) :: ICE5G
    
    ! Open the netcdf file
    CALL open_netcdf_file( ICE5G%netcdf%filename, ICE5G%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( ICE5G%netcdf%ncid, ICE5G%netcdf%id_var_lon,      ICE5G%grid%lon, start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( ICE5G%netcdf%ncid, ICE5G%netcdf%id_var_lat,      ICE5G%grid%lat, start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( ICE5G%netcdf%ncid, ICE5G%netcdf%id_var_Hi,       ICE5G%Hi,       start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( ICE5G%netcdf%ncid, ICE5G%netcdf%id_var_Hb,       ICE5G%Hb,       start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( ICE5G%netcdf%ncid, ICE5G%netcdf%id_var_mask_ice, ICE5G%mask_ice, start = (/ 1, 1 /) ))
    
    ! For some reason, "orog" in ICE5G is Hb+Hi (so like Hs without any oceans)
    ICE5G%Hb = ICE5G%Hb - ICE5G%Hi
        
    ! Close the netcdf file
    CALL close_netcdf_file(ICE5G%netcdf%ncid)
    
  END SUBROUTINE read_ICE5G_data
  
  ! Insolation solution (e.g. Laskar 2004)
  SUBROUTINE inquire_insolation_data_file( forcing)
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
 
    ! Local variables:  
    INTEGER                                :: int_dummy
            
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
    
  END SUBROUTINE inquire_insolation_data_file
  SUBROUTINE read_insolation_data_file( forcing, ti0, ti1, ins_Q_TOA0, ins_Q_TOA1) 
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    INTEGER,                        INTENT(IN)    :: ti0, ti1
    REAL(dp), DIMENSION(:,:),       INTENT(OUT)   :: ins_Q_TOA0, ins_Q_TOA1
    
    ! Local variables:
    INTEGER                                       :: mi, li
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Q_temp0, Q_temp1
    
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
   
  END SUBROUTINE read_insolation_data_file
  SUBROUTINE read_insolation_data_file_time_lat( forcing) 
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
    
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_time,    forcing%ins_time,    start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_lat,     forcing%ins_lat,     start = (/ 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ins%ncid)
    
  END SUBROUTINE read_insolation_data_file_time_lat
  
  ! Geothermal heat flux
  SUBROUTINE inquire_geothermal_heat_flux_file( forcing)
    IMPLICIT NONE

    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing

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

  END SUBROUTINE inquire_geothermal_heat_flux_file
  SUBROUTINE read_geothermal_heat_flux_file( forcing)
  
    USE parameters_module, ONLY: sec_per_year
  
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing

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

  END SUBROUTINE read_geothermal_heat_flux_file
  
  ! Write a CSR matrix to a NetCDF file
  SUBROUTINE write_CSR_matrix_to_NetCDF( AA, filename)
    ! Write a CSR matrix to a NetCDF file
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    
    ! Local variables:
    LOGICAL                                       :: file_exists
    INTEGER                                       :: ncid
    INTEGER                                       :: m, mp1, n, nnz, ptr, index, val
    
    IF (.NOT. par%master) RETURN

    ! Safety
    INQUIRE(EXIST=file_exists, FILE = TRIM( filename))
    IF (file_exists) THEN
      WRITE(0,*) 'write_CSR_matrix_to_NetCDF - ERROR: file "', TRIM(filename), '" already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    WRITE(0,*) '   NOTE: writing CSR matrix to file "', TRIM(filename), '"'
    
    ! Create netCDF file
    CALL handle_error( nf90_create( TRIM(C%output_dir)//TRIM(filename), IOR(nf90_clobber,nf90_share), ncid))
        
    ! Define dimensions:
    CALL create_dim( ncid, 'm',      AA%m  , m  )
    CALL create_dim( ncid, 'mplus1', AA%m+1, mp1)
    CALL create_dim( ncid, 'n',      AA%n  , n  )
    CALL create_dim( ncid, 'nnz',    AA%nnz, nnz)
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    CALL create_int_var(    ncid, 'ptr',   [mp1], ptr  , long_name = 'ptr'  )
    CALL create_int_var(    ncid, 'index', [nnz], index, long_name = 'index')
    CALL create_double_var( ncid, 'val',   [nnz], val  , long_name = 'val'  )
       
    ! Leave definition mode
    CALL handle_error( nf90_enddef( ncid))    
    
    ! Write data
    CALL handle_error( nf90_put_var( ncid, ptr  , AA%ptr  ))
    CALL handle_error( nf90_put_var( ncid, index, AA%index))
    CALL handle_error( nf90_put_var( ncid, val  , AA%val  ))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( ncid))
    
    ! Close the file
    CALL close_netcdf_file(ncid)
    
  END SUBROUTINE write_CSR_matrix_to_NetCDF

! Basic NetCDF wrapper functions
! ==============================

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
  SUBROUTINE create_int_var(    ncid, var_name, id_dims, id_var, long_name, units, missing_value)
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
     WRITE(0,'(3A)') 'ERROR: Actual type of variable "',var_name,'" is not nf90_DOUBLE.'
     STOP
    END IF
    IF(ndims /= SIZE(id_dims)) THEN
     WRITE(0,'(A,I5,3A,I5,A)') 'ERROR: Actual number of dimensions(', &
            ndims,') of variable "',var_name,'": does not match required number of dimensions (',SIZE(id_dims),').'
     STOP
    END IF
    IF(ANY(actual_id_dims(1:ndims) /= id_dims)) THEN
     WRITE(0,'(3A)') 'ERROR: Actual dimensions of variable "',var_name,'" does not match required dimensions.'
     STOP
    END IF
    
  END SUBROUTINE inquire_double_var
  SUBROUTINE inquire_single_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_float.
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
    IF(xtype /= nf90_float) THEN
     WRITE(0,'(3A)') 'ERROR: Actual type of variable "',var_name,'" is not nf90_float.'
     STOP
    END IF
    IF(ndims /= SIZE(id_dims)) THEN
     WRITE(0,'(A,I5,3A,I5,A)') 'ERROR: Actual number of dimensions(', &
            ndims,') of variable "',var_name,'": does not match required number of dimensions (',SIZE(id_dims),').'
     STOP
    END IF
    IF(ANY(actual_id_dims(1:ndims) /= id_dims)) THEN
     WRITE(0,'(3A)') 'ERROR: Actual dimensions of variable "',var_name,'" does not match required dimensions.'
     STOP
    END IF
    
  END SUBROUTINE inquire_single_var
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
    IF(xtype /= nf90_int) THEN
     WRITE(0,'(3A)') 'ERROR: Actual type of variable "',var_name,'" is not nf90_int.'
     STOP
    END IF
    IF(ndims /= SIZE(id_dims)) THEN
     WRITE(0,'(A,I5,3A,I5,A)') 'ERROR: Actual number of dimensions(', &
            ndims,') of variable "',var_name,'": does not match required number of dimensions (',SIZE(id_dims),').'
     STOP
    END IF
    IF(ANY(actual_id_dims(1:ndims) /= id_dims)) THEN
     WRITE(0,'(3A)') 'ERROR: Actual dimensions of variable "',var_name,'" does not match required dimensions.'
     STOP
    END IF
    
  END SUBROUTINE inquire_int_var
  SUBROUTINE handle_error( stat, message)
    USE netcdf, ONLY: nf90_noerr, nf90_strerror
    IMPLICIT NONE
    
    ! Input variables:
    INTEGER,                    INTENT(IN) :: stat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: message

    IF(stat /= nf90_noerr) THEN
     IF(PRESENT(message)) THEN
      WRITE(0,'(A,A,A,A)') 'ERROR: ', TRIM(nf90_strerror(stat)), ' concerning: ', message
     ELSE
      WRITE(0,'(A,A)')     'ERROR: ', TRIM(nf90_strerror(stat))
     END IF
     STOP
    END IF
    
  END SUBROUTINE handle_error

  !==========================

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

END MODULE netcdf_module
