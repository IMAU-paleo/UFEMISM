MODULE netcdf_output_module

! ===== Creating and writing to output files =====
! ================================================
!
! These routines create and write data to output files, both on meshes and on grids. For
! grid files, remapping of data from the provided model mesh to the output grid is done
! automatically.

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
                                             allocate_shared_int_4D,   allocate_shared_dp_4D, &
                                             allocate_shared_bool_0D,  allocate_shared_bool_1D, &
                                             reallocate_shared_int_0D, reallocate_shared_dp_0D, &
                                             reallocate_shared_int_1D, reallocate_shared_dp_1D, &
                                             reallocate_shared_int_2D, reallocate_shared_dp_2D, &
                                             reallocate_shared_int_3D, reallocate_shared_dp_3D, &
                                             deallocate_shared
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh, type_grid, type_grid_lonlat, type_remapping_lonlat2mesh, &
                                             type_remapping_mesh_mesh, type_ice_model, type_model_region
  USE netcdf,                          ONLY: NF90_UNLIMITED, NF90_INT, NF90_FLOAT, NF90_DOUBLE
  USE netcdf_basic_module,             ONLY: nerr, field_name_options_x, field_name_options_y, field_name_options_zeta, field_name_options_z_ocean, &
                                             field_name_options_lon, field_name_options_lat, field_name_options_time, field_name_options_month, &
                                             field_name_options_dim_nV, field_name_options_dim_nTri, field_name_options_dim_nC_mem, &
                                             field_name_options_dim_nAc, field_name_options_dim_two, field_name_options_dim_three, &
                                             field_name_options_dim_six, field_name_options_V, field_name_options_Tri, field_name_options_nC, &
                                             field_name_options_C, field_name_options_niTri, field_name_options_iTri, &
                                             field_name_options_edge_index, field_name_options_Tricc, field_name_options_TriC, &
                                             field_name_options_Tri_edge_index, field_name_options_VAc, field_name_options_Aci, &
                                             field_name_options_iAci, field_name_options_A, field_name_options_R, &
                                             field_name_options_Hi, field_name_options_Hb, field_name_options_Hs, field_name_options_dHb, &
                                             field_name_options_SL, field_name_options_Ti, &
                                             inquire_dim_multiple_options, inquire_var_multiple_options, get_first_option_from_list, &
                                             check_x, check_y, check_lon, check_lat, check_mesh_dimensions, check_zeta, check_month, check_time, &
                                             create_new_netcdf_file_for_writing, create_dimension, create_variable, add_attribute_char, &
                                             write_var_int_0D, write_var_int_1D, write_var_int_2D, write_var_int_3D, write_var_int_4D, &
                                             write_var_dp_0D , write_var_dp_1D , write_var_dp_2D , write_var_dp_3D , write_var_dp_4D, &
                                             check_xy_grid_field_int_2D, check_xy_grid_field_dp_2D, check_xy_grid_field_dp_2D_monthly, check_xy_grid_field_dp_3D, &
                                             check_lonlat_grid_field_int_2D, check_lonlat_grid_field_dp_2D, check_lonlat_grid_field_dp_2D_monthly, check_lonlat_grid_field_dp_3D, &
                                             check_mesh_field_int_2D, check_mesh_field_dp_2D, check_mesh_field_dp_2D_monthly, check_mesh_field_dp_3D
  USE netcdf_input_module,             ONLY: setup_xy_grid_from_file
  USE mesh_mapping_module,             ONLY: calc_remapping_operator_mesh2grid, deallocate_remapping_operators_mesh2grid, map_mesh2grid_2D, map_mesh2grid_3D
  USE utilities_module,                ONLY: deallocate_grid

  USE netcdf_module, ONLY: debug, write_to_debug_file

  IMPLICIT NONE

CONTAINS

  ! ===== Top-level functions =====
  ! ===============================

  ! Create restart and help fields files
  SUBROUTINE create_restart_file_mesh( filename, mesh)
    ! Create an empty restart file on the model mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_restart_file_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Set up the mesh in this file
    CALL setup_mesh_in_netcdf_file( filename, mesh)

    ! Add time, zeta, and month dimensions
    CALL add_time_dimension_to_file(  filename)
    CALL add_zeta_dimension_to_file(  filename)
    CALL add_month_dimension_to_file( filename)

    ! Create variables
    ! ================

    ! Geometry
    CALL add_field_mesh_dp_2D( filename, get_first_option_from_list( field_name_options_Hi ), long_name = 'Ice thickness'      , units = 'm')
    CALL add_field_mesh_dp_2D( filename, get_first_option_from_list( field_name_options_Hb ), long_name = 'Bedrock elevation'  , units = 'm w.r.t. PD sea level')
    CALL add_field_mesh_dp_2D( filename, get_first_option_from_list( field_name_options_Hs ), long_name = 'Surface elevation'  , units = 'm w.r.t. PD sea level')
    CALL add_field_mesh_dp_2D( filename, get_first_option_from_list( field_name_options_SL ), long_name = 'Geoid elevation'    , units = 'm w.r.t. PD sea level')
    CALL add_field_mesh_dp_2D( filename, get_first_option_from_list( field_name_options_dHB), long_name = 'Bedrock deformation', units = 'm')

    ! Englacial temperature
    CALL add_field_mesh_dp_3D( filename, get_first_option_from_list( field_name_options_Ti ), long_name = 'Englacial temperature'    , units = 'K')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_mesh

  SUBROUTINE create_restart_file_grid( filename, grid)
    ! Create an empty restart file with the specified grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid),                     INTENT(IN)    :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_restart_file_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Set up the mesh in this file
    CALL setup_xy_grid_in_netcdf_file( filename, grid)

    ! Add time, zeta, and month dimensions
    CALL add_time_dimension_to_file(  filename)
    CALL add_zeta_dimension_to_file(  filename)
    CALL add_month_dimension_to_file( filename)

    ! Create variables
    ! ================

    ! Geometry
    CALL add_field_grid_dp_2D( filename, get_first_option_from_list( field_name_options_Hi ), long_name = 'Ice thickness'      , units = 'm')
    CALL add_field_grid_dp_2D( filename, get_first_option_from_list( field_name_options_Hb ), long_name = 'Bedrock elevation'  , units = 'm w.r.t. PD sea level')
    CALL add_field_grid_dp_2D( filename, get_first_option_from_list( field_name_options_Hs ), long_name = 'Surface elevation'  , units = 'm w.r.t. PD sea level')
    CALL add_field_grid_dp_2D( filename, get_first_option_from_list( field_name_options_SL ), long_name = 'Geoid elevation'    , units = 'm w.r.t. PD sea level')
    CALL add_field_grid_dp_2D( filename, get_first_option_from_list( field_name_options_dHB), long_name = 'Bedrock deformation', units = 'm')

    ! Englacial temperature
    CALL add_field_grid_dp_3D( filename, get_first_option_from_list( field_name_options_Ti ), long_name = 'Englacial temperature'    , units = 'K')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_grid

  SUBROUTINE create_help_fields_file_mesh( filename, mesh)
    ! Create an empty help fields file on the model mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_fields_file_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Set up the mesh in this file
    CALL setup_mesh_in_netcdf_file( filename, mesh)

    ! Add time, zeta, and month dimensions
    CALL add_time_dimension_to_file(  filename)
    CALL add_zeta_dimension_to_file(  filename)
    CALL add_month_dimension_to_file( filename)

    ! Create variables
    ! ================

    CALL create_help_field_mesh( filename, C%help_field_01)
    CALL create_help_field_mesh( filename, C%help_field_02)
    CALL create_help_field_mesh( filename, C%help_field_03)
    CALL create_help_field_mesh( filename, C%help_field_04)
    CALL create_help_field_mesh( filename, C%help_field_05)
    CALL create_help_field_mesh( filename, C%help_field_06)
    CALL create_help_field_mesh( filename, C%help_field_07)
    CALL create_help_field_mesh( filename, C%help_field_08)
    CALL create_help_field_mesh( filename, C%help_field_09)
    CALL create_help_field_mesh( filename, C%help_field_10)
    CALL create_help_field_mesh( filename, C%help_field_11)
    CALL create_help_field_mesh( filename, C%help_field_12)
    CALL create_help_field_mesh( filename, C%help_field_13)
    CALL create_help_field_mesh( filename, C%help_field_14)
    CALL create_help_field_mesh( filename, C%help_field_15)
    CALL create_help_field_mesh( filename, C%help_field_16)
    CALL create_help_field_mesh( filename, C%help_field_17)
    CALL create_help_field_mesh( filename, C%help_field_18)
    CALL create_help_field_mesh( filename, C%help_field_19)
    CALL create_help_field_mesh( filename, C%help_field_20)
    CALL create_help_field_mesh( filename, C%help_field_21)
    CALL create_help_field_mesh( filename, C%help_field_22)
    CALL create_help_field_mesh( filename, C%help_field_23)
    CALL create_help_field_mesh( filename, C%help_field_24)
    CALL create_help_field_mesh( filename, C%help_field_25)
    CALL create_help_field_mesh( filename, C%help_field_26)
    CALL create_help_field_mesh( filename, C%help_field_27)
    CALL create_help_field_mesh( filename, C%help_field_28)
    CALL create_help_field_mesh( filename, C%help_field_29)
    CALL create_help_field_mesh( filename, C%help_field_30)
    CALL create_help_field_mesh( filename, C%help_field_31)
    CALL create_help_field_mesh( filename, C%help_field_32)
    CALL create_help_field_mesh( filename, C%help_field_33)
    CALL create_help_field_mesh( filename, C%help_field_34)
    CALL create_help_field_mesh( filename, C%help_field_35)
    CALL create_help_field_mesh( filename, C%help_field_36)
    CALL create_help_field_mesh( filename, C%help_field_37)
    CALL create_help_field_mesh( filename, C%help_field_38)
    CALL create_help_field_mesh( filename, C%help_field_39)
    CALL create_help_field_mesh( filename, C%help_field_40)
    CALL create_help_field_mesh( filename, C%help_field_41)
    CALL create_help_field_mesh( filename, C%help_field_42)
    CALL create_help_field_mesh( filename, C%help_field_43)
    CALL create_help_field_mesh( filename, C%help_field_44)
    CALL create_help_field_mesh( filename, C%help_field_45)
    CALL create_help_field_mesh( filename, C%help_field_46)
    CALL create_help_field_mesh( filename, C%help_field_47)
    CALL create_help_field_mesh( filename, C%help_field_48)
    CALL create_help_field_mesh( filename, C%help_field_49)
    CALL create_help_field_mesh( filename, C%help_field_50)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_fields_file_mesh

  SUBROUTINE create_help_fields_file_grid( filename, grid)
    ! Create an empty help fields file with the specified grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid),                     INTENT(IN)    :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_fields_file_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Set up the mesh in this file
    CALL setup_xy_grid_in_netcdf_file( filename, grid)

    ! Add time, zeta, and month dimensions
    CALL add_time_dimension_to_file(  filename)
    CALL add_zeta_dimension_to_file(  filename)
    CALL add_month_dimension_to_file( filename)

    ! Create variables
    ! ================

    CALL create_help_field_grid( filename, C%help_field_01)
    CALL create_help_field_grid( filename, C%help_field_02)
    CALL create_help_field_grid( filename, C%help_field_03)
    CALL create_help_field_grid( filename, C%help_field_04)
    CALL create_help_field_grid( filename, C%help_field_05)
    CALL create_help_field_grid( filename, C%help_field_06)
    CALL create_help_field_grid( filename, C%help_field_07)
    CALL create_help_field_grid( filename, C%help_field_08)
    CALL create_help_field_grid( filename, C%help_field_09)
    CALL create_help_field_grid( filename, C%help_field_10)
    CALL create_help_field_grid( filename, C%help_field_11)
    CALL create_help_field_grid( filename, C%help_field_12)
    CALL create_help_field_grid( filename, C%help_field_13)
    CALL create_help_field_grid( filename, C%help_field_14)
    CALL create_help_field_grid( filename, C%help_field_15)
    CALL create_help_field_grid( filename, C%help_field_16)
    CALL create_help_field_grid( filename, C%help_field_17)
    CALL create_help_field_grid( filename, C%help_field_18)
    CALL create_help_field_grid( filename, C%help_field_19)
    CALL create_help_field_grid( filename, C%help_field_20)
    CALL create_help_field_grid( filename, C%help_field_21)
    CALL create_help_field_grid( filename, C%help_field_22)
    CALL create_help_field_grid( filename, C%help_field_23)
    CALL create_help_field_grid( filename, C%help_field_24)
    CALL create_help_field_grid( filename, C%help_field_25)
    CALL create_help_field_grid( filename, C%help_field_26)
    CALL create_help_field_grid( filename, C%help_field_27)
    CALL create_help_field_grid( filename, C%help_field_28)
    CALL create_help_field_grid( filename, C%help_field_29)
    CALL create_help_field_grid( filename, C%help_field_30)
    CALL create_help_field_grid( filename, C%help_field_31)
    CALL create_help_field_grid( filename, C%help_field_32)
    CALL create_help_field_grid( filename, C%help_field_33)
    CALL create_help_field_grid( filename, C%help_field_34)
    CALL create_help_field_grid( filename, C%help_field_35)
    CALL create_help_field_grid( filename, C%help_field_36)
    CALL create_help_field_grid( filename, C%help_field_37)
    CALL create_help_field_grid( filename, C%help_field_38)
    CALL create_help_field_grid( filename, C%help_field_39)
    CALL create_help_field_grid( filename, C%help_field_40)
    CALL create_help_field_grid( filename, C%help_field_41)
    CALL create_help_field_grid( filename, C%help_field_42)
    CALL create_help_field_grid( filename, C%help_field_43)
    CALL create_help_field_grid( filename, C%help_field_44)
    CALL create_help_field_grid( filename, C%help_field_45)
    CALL create_help_field_grid( filename, C%help_field_46)
    CALL create_help_field_grid( filename, C%help_field_47)
    CALL create_help_field_grid( filename, C%help_field_48)
    CALL create_help_field_grid( filename, C%help_field_49)
    CALL create_help_field_grid( filename, C%help_field_50)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_fields_file_grid

  SUBROUTINE create_help_field_mesh( filename, field_name)
    ! Add a data field to the help_fields file

    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_field_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    ELSEIF (field_name == 'lon') THEN
      CALL warning('longitude is already by default written to output!')
    ELSEIF (field_name == 'lat') THEN
      CALL warning('latitude is already by default written to output!')

    ! Mesh data
    ELSEIF (field_name == 'resolution') THEN
      CALL warning('mesh resolution is already by default written to output!')

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL add_field_mesh_dp_2D_notime( filename, 'GHF', long_name = 'Geothermal heat flux', units = 'J m^-2 yr^-1')

    ! Fields with a time dimension
    ! ============================

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL add_field_mesh_dp_2D( filename, 'Hi', long_name = 'Ice thickness', units = 'm')
    ELSEIF (field_name == 'Hb') THEN
      CALL add_field_mesh_dp_2D( filename, 'Hb', long_name = 'Bedrock elevation', units = 'm w.r.t. PD sea-level')
    ELSEIF (field_name == 'Hs') THEN
      CALL add_field_mesh_dp_2D( filename, 'Hs', long_name = 'Surface elevation', units = 'm w.r.t. PD sea-level')
    ELSEIF (field_name == 'SL') THEN
      CALL add_field_mesh_dp_2D( filename, 'SL', long_name = 'Geoid elevation', units = 'm w.r.t. PD sea-level')
    ELSEIF (field_name == 'dHi') THEN
      CALL add_field_mesh_dp_2D( filename, 'dHi', long_name = 'Ice thickness difference w.r.t. PD', units = 'm')
    ELSEIF (field_name == 'dHb') THEN
      CALL add_field_mesh_dp_2D( filename, 'dHb', long_name = 'Bedrock elevation difference w.r.t. PD', units = 'm')
    ELSEIF (field_name == 'dHs') THEN
      CALL add_field_mesh_dp_2D( filename, 'dHs', long_name = 'Surface elevation difference w.r.t. PD', units = 'm')

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL add_field_mesh_dp_3D( filename, 'Ti', long_name = 'Englacial temperature', units = 'K')
    ELSEIF (field_name == 'Cpi') THEN
      CALL add_field_mesh_dp_3D( filename, 'Cpi', long_name = 'Ice heat capacity', units = 'J kg^-1 K^-1')
    ELSEIF (field_name == 'Ki') THEN
      CALL add_field_mesh_dp_3D( filename, 'Ki', long_name = 'Ice thermal conductivity', units = 'J m^-1 K^-1 yr^-1')
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL add_field_mesh_dp_2D( filename, 'Ti_basal', long_name = 'Ice basal temperature', units = 'K')
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL add_field_mesh_dp_3D( filename, 'Ti_pmp', long_name = 'Englacial pressure melting point temperature', units = 'K')
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL add_field_mesh_dp_3D( filename, 'A_flow_3D', long_name = 'Ice flow factor', units = 'Pa^-3 y^-1')
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL add_field_mesh_dp_2D( filename, 'A_flow_vav', long_name = 'Vertically averaged ice flow factor', units = 'Pa^-3 y^-1')

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL add_field_mesh_dp_3D( filename, 'u_3D', long_name = '3D ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_3D') THEN
      CALL add_field_mesh_dp_3D( filename, 'v_3D', long_name = '3D ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'w_3D') THEN
      CALL add_field_mesh_dp_3D( filename, 'w_3D', long_name = '3D ice z-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'u_vav') THEN
      CALL add_field_mesh_dp_2D( filename, 'u_vav', long_name = 'Vertically averaged ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_vav') THEN
      CALL add_field_mesh_dp_2D( filename, 'v_vav', long_name = 'Vertically averaged ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL add_field_mesh_dp_2D( filename, 'uabs_vav', long_name = 'Vertically averaged ice velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'u_surf') THEN
      CALL add_field_mesh_dp_2D( filename, 'u_surf', long_name = 'Surface ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_surf') THEN
      CALL add_field_mesh_dp_2D( filename, 'v_surf', long_name = 'Surface ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL add_field_mesh_dp_2D( filename, 'uabs_surf', long_name = 'Surface ice velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'u_base') THEN
      CALL add_field_mesh_dp_2D( filename, 'u_base', long_name = 'Basal ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_base') THEN
      CALL add_field_mesh_dp_2D( filename, 'v_base', long_name = 'Basal ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'uabs_base') THEN
      CALL add_field_mesh_dp_2D( filename, 'uabs_base', long_name = 'Basal ice velocity', units = 'm yr^-1')

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'T2m', long_name = 'Monthly mean 2-m air temperature', units = 'K')
    ELSEIF (field_name == 'T2m_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'T2m_year', long_name = 'Annual mean 2-m air temperature', units = 'K')
    ELSEIF (field_name == 'Precip') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Precip', long_name = 'Monthly total precipitation', units = 'mm')
    ELSEIF (field_name == 'Precip_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Precip_year', long_name = 'Annual total precipitation', units = 'mm')
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Wind_WE', long_name = 'Monthly mean zonal wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Wind_WE_year', long_name = 'Annual mean zonal wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Wind_SN', long_name = 'Monthly mean meridional wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Wind_SN_year', long_name = 'Annual mean meridional wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_LR') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Wind_LR', long_name = 'Monthly mean x-wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_LR_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Wind_LR_year', long_name = 'Annual mean x-wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_DU') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Wind_DU', long_name = 'Monthly mean y-wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_DU_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Wind_DU_year', long_name = 'Annual mean y-wind', units = 'm s^-1')

    ! Surface mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'SMB', long_name = 'Monthly surface mass balance', units = 'm ice equivalent')
    ELSEIF (field_name == 'SMB_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'SMB_year', long_name = 'Annual surface mass balance', units = 'm ice equivalent')
    ELSEIF (field_name == 'Snowfall') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Snowfall', long_name = 'Monthly total snowfall', units = 'm water equivalent')
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Snowfall_year', long_name = 'Annual total snowfall', units = 'm water equivalent')
    ELSEIF (field_name == 'Rainfall') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Rainfall', long_name = 'Monthly total rainfall', units = 'm water equivalent')
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Rainfall_year', long_name = 'Annual total rainfall', units = 'm water equivalent')
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'AddedFirn', long_name = 'Monthly total added firn', units = 'm water equivalent')
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'AddedFirn_year', long_name = 'Annual total added firn', units = 'm water equivalent')
    ELSEIF (field_name == 'Refreezing') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Refreezing', long_name = 'Monthly total refreezing', units = 'm water equivalent')
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Refreezing_year', long_name = 'Annual total refreezing', units = 'm water equivalent')
    ELSEIF (field_name == 'Melt') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Melt', long_name = 'Monthly total melt', units = 'm water equivalent')
    ELSEIF (field_name == 'Melt_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Melt_year', long_name = 'Annual total melt', units = 'm water equivalent')
    ELSEIF (field_name == 'Runoff') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Runoff', long_name = 'Monthly total runoff', units = 'm water equivalent')
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Runoff_year', long_name = 'Annual total runoff', units = 'm water equivalent')
    ELSEIF (field_name == 'Albedo') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'Albedo', long_name = 'Monthly mean albedo')
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'Albedo_year', long_name = 'Annual mean albedo')
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL add_field_mesh_dp_2D_monthly( filename, 'FirnDepth', long_name = 'Monthly mean firn layer depth', units='m water equivalent')
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL add_field_mesh_dp_2D( filename, 'FirnDepth_year', long_name = 'Annual mean firn layer depth', units='m water equivalent')

    ! Basal mass balance
    ELSEIF (field_name == 'BMB') THEN
      CALL add_field_mesh_dp_2D( filename, 'BMB', long_name = 'Annual basal mass balance', units = 'm ice equivalent')
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL add_field_mesh_dp_2D( filename, 'BMB_sheet', long_name = 'Annual basal mass balance for grounded ice', units = 'm ice equivalent')
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL add_field_mesh_dp_2D( filename, 'BMB_shelf', long_name = 'Annual basal mass balance for floating ice', units = 'm ice equivalent')

    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL add_field_mesh_int_2D( filename, 'mask', long_name = 'Mask')
    ELSEIF (field_name == 'mask_land') THEN
      CALL add_field_mesh_int_2D( filename, 'mask_land', long_name = 'Land mask')
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL add_field_mesh_int_2D( filename, 'mask_ocean', long_name = 'Ocean mask')
    ELSEIF (field_name == 'mask_lake') THEN
      CALL add_field_mesh_int_2D( filename, 'mask_lake', long_name = 'Lake mask')
    ELSEIF (field_name == 'mask_ice') THEN
      CALL add_field_mesh_int_2D( filename, 'mask_ice', long_name = 'Ice mask')
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL add_field_mesh_int_2D( filename, 'mask_sheet', long_name = 'Sheet mask')
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL add_field_mesh_int_2D( filename, 'mask_shelf', long_name = 'Shelf mask')
    ELSEIF (field_name == 'mask_coast') THEN
      CALL add_field_mesh_int_2D( filename, 'mask_coast', long_name = 'Coast mask')
    ELSEIF (field_name == 'mask_margin') THEN
      CALL add_field_mesh_int_2D( filename, 'mask_margin', long_name = 'Margin mask')
    ELSEIF (field_name == 'mask_gl') THEN
      CALL add_field_mesh_int_2D( filename, 'mask_gl', long_name = 'Grounding-line mask')
    ELSEIF (field_name == 'mask_cf') THEN
      CALL add_field_mesh_int_2D( filename, 'mask_cf', long_name = 'Calving-front mask')

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL add_field_mesh_dp_2D( filename, 'phi_fric', long_name = 'Till friction angle', units = 'degrees')
    ELSEIF (field_name == 'tau_yield') THEN
      CALL add_field_mesh_dp_2D( filename, 'tau_yield', long_name = 'Basal yield stress', units = 'Pa')

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL add_field_mesh_dp_2D( filename, 'iso_ice', long_name = 'Vertically averaged englacial d18O', units = 'per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL add_field_mesh_dp_2D( filename, 'iso_surf', long_name = 'd18O of precipitation', units = 'per mille')

    ! GIA
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL add_field_mesh_dp_2D( filename, 'dHb_dt', long_name = 'Bedrock deformation rate', units = 'm yr^-1')
    ELSEIF (field_name == 'dSL_dt') THEN
      CALL add_field_mesh_dp_2D( filename, 'dSL_dt', long_name = 'Geoid deformation rate', units = 'm yr^-1')

    ! Safety
    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_field_mesh

  SUBROUTINE create_help_field_grid( filename, field_name)
    ! Add a data field to the help_fields file

    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_field_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    ELSEIF (field_name == 'lon') THEN
      CALL warning('longitude is already by default written to output!')
    ELSEIF (field_name == 'lat') THEN
      CALL warning('latitude is already by default written to output!')

    ! Mesh data
    ELSEIF (field_name == 'resolution') THEN
      CALL warning('mesh resolution is already by default written to output!')

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL add_field_grid_dp_2D_notime( filename, 'GHF', long_name = 'Geothermal heat flux', units = 'J m^-2 yr^-1')

    ! Fields with a time dimension
    ! ============================

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL add_field_grid_dp_2D( filename, 'Hi', long_name = 'Ice thickness', units = 'm')
    ELSEIF (field_name == 'Hb') THEN
      CALL add_field_grid_dp_2D( filename, 'Hb', long_name = 'Bedrock elevation', units = 'm w.r.t. PD sea-level')
    ELSEIF (field_name == 'Hs') THEN
      CALL add_field_grid_dp_2D( filename, 'Hs', long_name = 'Surface elevation', units = 'm w.r.t. PD sea-level')
    ELSEIF (field_name == 'SL') THEN
      CALL add_field_grid_dp_2D( filename, 'SL', long_name = 'Geoid elevation', units = 'm w.r.t. PD sea-level')
    ELSEIF (field_name == 'dHi') THEN
      CALL add_field_grid_dp_2D( filename, 'dHi', long_name = 'Ice thickness difference w.r.t. PD', units = 'm')
    ELSEIF (field_name == 'dHb') THEN
      CALL add_field_grid_dp_2D( filename, 'dHb', long_name = 'Bedrock elevation difference w.r.t. PD', units = 'm')
    ELSEIF (field_name == 'dHs') THEN
      CALL add_field_grid_dp_2D( filename, 'dHs', long_name = 'Surface elevation difference w.r.t. PD', units = 'm')

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL add_field_grid_dp_3D( filename, 'Ti', long_name = 'Englacial temperature', units = 'K')
    ELSEIF (field_name == 'Cpi') THEN
      CALL add_field_grid_dp_3D( filename, 'Cpi', long_name = 'Ice heat capacity', units = 'J kg^-1 K^-1')
    ELSEIF (field_name == 'Ki') THEN
      CALL add_field_grid_dp_3D( filename, 'Ki', long_name = 'Ice thermal conductivity', units = 'J m^-1 K^-1 yr^-1')
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL add_field_grid_dp_2D( filename, 'Ti_basal', long_name = 'Ice basal temperature', units = 'K')
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL add_field_grid_dp_3D( filename, 'Ti_pmp', long_name = 'Englacial pressure melting point temperature', units = 'K')
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL add_field_grid_dp_3D( filename, 'A_flow_3D', long_name = 'Ice flow factor', units = 'Pa^-3 y^-1')
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL add_field_grid_dp_2D( filename, 'A_flow_vav', long_name = 'Vertically averaged ice flow factor', units = 'Pa^-3 y^-1')

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL add_field_grid_dp_3D( filename, 'u_3D', long_name = '3D ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_3D') THEN
      CALL add_field_grid_dp_3D( filename, 'v_3D', long_name = '3D ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'w_3D') THEN
      CALL add_field_grid_dp_3D( filename, 'w_3D', long_name = '3D ice z-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'u_vav') THEN
      CALL add_field_grid_dp_2D( filename, 'u_vav', long_name = 'Vertically averaged ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_vav') THEN
      CALL add_field_grid_dp_2D( filename, 'v_vav', long_name = 'Vertically averaged ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL add_field_grid_dp_2D( filename, 'uabs_vav', long_name = 'Vertically averaged ice velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'u_surf') THEN
      CALL add_field_grid_dp_2D( filename, 'u_surf', long_name = 'Surface ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_surf') THEN
      CALL add_field_grid_dp_2D( filename, 'v_surf', long_name = 'Surface ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL add_field_grid_dp_2D( filename, 'uabs_surf', long_name = 'Surface ice velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'u_base') THEN
      CALL add_field_grid_dp_2D( filename, 'u_base', long_name = 'Basal ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_base') THEN
      CALL add_field_grid_dp_2D( filename, 'v_base', long_name = 'Basal ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'uabs_base') THEN
      CALL add_field_grid_dp_2D( filename, 'uabs_base', long_name = 'Basal ice velocity', units = 'm yr^-1')

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'T2m', long_name = 'Monthly mean 2-m air temperature', units = 'K')
    ELSEIF (field_name == 'T2m_year') THEN
      CALL add_field_grid_dp_2D( filename, 'T2m_year', long_name = 'Annual mean 2-m air temperature', units = 'K')
    ELSEIF (field_name == 'Precip') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Precip', long_name = 'Monthly total precipitation', units = 'mm')
    ELSEIF (field_name == 'Precip_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Precip_year', long_name = 'Annual total precipitation', units = 'mm')
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Wind_WE', long_name = 'Monthly mean zonal wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Wind_WE_year', long_name = 'Annual mean zonal wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Wind_SN', long_name = 'Monthly mean meridional wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Wind_SN_year', long_name = 'Annual mean meridional wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_LR') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Wind_LR', long_name = 'Monthly mean x-wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_LR_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Wind_LR_year', long_name = 'Annual mean x-wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_DU') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Wind_DU', long_name = 'Monthly mean y-wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_DU_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Wind_DU_year', long_name = 'Annual mean y-wind', units = 'm s^-1')

    ! Surface mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'SMB', long_name = 'Monthly surface mass balance', units = 'm ice equivalent')
    ELSEIF (field_name == 'SMB_year') THEN
      CALL add_field_grid_dp_2D( filename, 'SMB_year', long_name = 'Annual surface mass balance', units = 'm ice equivalent')
    ELSEIF (field_name == 'Snowfall') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Snowfall', long_name = 'Monthly total snowfall', units = 'm water equivalent')
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Snowfall_year', long_name = 'Annual total snowfall', units = 'm water equivalent')
    ELSEIF (field_name == 'Rainfall') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Rainfall', long_name = 'Monthly total rainfall', units = 'm water equivalent')
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Rainfall_year', long_name = 'Annual total rainfall', units = 'm water equivalent')
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'AddedFirn', long_name = 'Monthly total added firn', units = 'm water equivalent')
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL add_field_grid_dp_2D( filename, 'AddedFirn_year', long_name = 'Annual total added firn', units = 'm water equivalent')
    ELSEIF (field_name == 'Refreezing') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Refreezing', long_name = 'Monthly total refreezing', units = 'm water equivalent')
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Refreezing_year', long_name = 'Annual total refreezing', units = 'm water equivalent')
    ELSEIF (field_name == 'Melt') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Melt', long_name = 'Monthly total melt', units = 'm water equivalent')
    ELSEIF (field_name == 'Melt_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Melt_year', long_name = 'Annual total melt', units = 'm water equivalent')
    ELSEIF (field_name == 'Runoff') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Runoff', long_name = 'Monthly total runoff', units = 'm water equivalent')
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Runoff_year', long_name = 'Annual total runoff', units = 'm water equivalent')
    ELSEIF (field_name == 'Albedo') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Albedo', long_name = 'Monthly mean albedo')
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Albedo_year', long_name = 'Annual mean albedo')
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'FirnDepth', long_name = 'Monthly mean firn layer depth', units='m water equivalent')
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL add_field_grid_dp_2D( filename, 'FirnDepth_year', long_name = 'Annual mean firn layer depth', units='m water equivalent')

    ! Basal mass balance
    ELSEIF (field_name == 'BMB') THEN
      CALL add_field_grid_dp_2D( filename, 'BMB', long_name = 'Annual basal mass balance', units = 'm ice equivalent')
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL add_field_grid_dp_2D( filename, 'BMB_sheet', long_name = 'Annual basal mass balance for grounded ice', units = 'm ice equivalent')
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL add_field_grid_dp_2D( filename, 'BMB_shelf', long_name = 'Annual basal mass balance for floating ice', units = 'm ice equivalent')

    ! Masks
    ! Not used, as mapping masks from the mesh to the grid is meaningless
    ELSEIF (field_name == 'mask') THEN
!      CALL add_field_grid_int_2D( filename, 'mask', long_name = 'Mask')
    ELSEIF (field_name == 'mask_land') THEN
!      CALL add_field_grid_int_2D( filename, 'mask_land', long_name = 'Land mask')
    ELSEIF (field_name == 'mask_ocean') THEN
!      CALL add_field_grid_int_2D( filename, 'mask_ocean', long_name = 'Ocean mask')
    ELSEIF (field_name == 'mask_lake') THEN
!      CALL add_field_grid_int_2D( filename, 'mask_lake', long_name = 'Lake mask')
    ELSEIF (field_name == 'mask_ice') THEN
!      CALL add_field_grid_int_2D( filename, 'mask_ice', long_name = 'Ice mask')
    ELSEIF (field_name == 'mask_sheet') THEN
!      CALL add_field_grid_int_2D( filename, 'mask_sheet', long_name = 'Sheet mask')
    ELSEIF (field_name == 'mask_shelf') THEN
!      CALL add_field_grid_int_2D( filename, 'mask_shelf', long_name = 'Shelf mask')
    ELSEIF (field_name == 'mask_coast') THEN
!      CALL add_field_grid_int_2D( filename, 'mask_coast', long_name = 'Coast mask')
    ELSEIF (field_name == 'mask_margin') THEN
!      CALL add_field_grid_int_2D( filename, 'mask_margin', long_name = 'Margin mask')
    ELSEIF (field_name == 'mask_gl') THEN
!      CALL add_field_grid_int_2D( filename, 'mask_gl', long_name = 'Grounding-line mask')
    ELSEIF (field_name == 'mask_cf') THEN
!      CALL add_field_grid_int_2D( filename, 'mask_cf', long_name = 'Calving-front mask')

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL add_field_grid_dp_2D( filename, 'phi_fric', long_name = 'Till friction angle', units = 'degrees')
    ELSEIF (field_name == 'tau_yield') THEN
      CALL add_field_grid_dp_2D( filename, 'tau_yield', long_name = 'Basal yield stress', units = 'Pa')

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL add_field_grid_dp_2D( filename, 'iso_ice', long_name = 'Vertically averaged englacial d18O', units = 'per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL add_field_grid_dp_2D( filename, 'iso_surf', long_name = 'd18O of precipitation', units = 'per mille')

    ! GIA
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL add_field_grid_dp_2D( filename, 'dHb_dt', long_name = 'Bedrock deformation rate', units = 'm yr^-1')
    ELSEIF (field_name == 'dSL_dt') THEN
      CALL add_field_grid_dp_2D( filename, 'dSL_dt', long_name = 'Geoid deformation rate', units = 'm yr^-1')

    ! Safety
    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_field_grid

  ! Write to restart and help fields files
  SUBROUTINE write_to_restart_file_mesh( filename, region)
    ! Write model output to the mesh restart file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_restart_file_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write new time to file (thus extending the time dimension by one frame, making room for the new model data)
    CALL write_time_to_file( filename, region%time)

    ! Write model data
    ! ================

    ! Geometry
    CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, field_name_options_Hi , region%ice%Hi_a )
    CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, field_name_options_Hb , region%ice%Hb_a )
    CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, field_name_options_Hs , region%ice%Hs_a )
    CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, field_name_options_dHb, region%ice%dHb_a)
    CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, field_name_options_SL , region%ice%SL_a )

    ! Englacial temperature
    CALL write_to_field_multiple_options_mesh_dp_3D( filename, region%mesh, field_name_options_Ti , region%ice%Ti_a )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_mesh

  SUBROUTINE write_to_restart_file_grid( filename, region)
    ! Write model output to the grid restart file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_restart_file_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write new time to file (thus extending the time dimension by one frame, making room for the new model data)
    CALL write_time_to_file( filename, region%time)

    ! Write model data
    ! ================

    ! Geometry
    CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, field_name_options_Hi , region%ice%Hi_a )
    CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, field_name_options_Hb , region%ice%Hb_a )
    CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, field_name_options_Hs , region%ice%Hs_a )
    CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, field_name_options_dHb, region%ice%dHb_a)
    CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, field_name_options_SL , region%ice%SL_a )

    ! Englacial temperature
    CALL write_to_field_multiple_options_grid_dp_3D( filename, region%mesh, region%name, field_name_options_Ti , region%ice%Ti_a )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_grid

  SUBROUTINE write_to_help_fields_file_mesh( filename, region)
    ! Write model output to the mesh help fields file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_help_fields_file_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write new time to file (thus extending the time dimension by one frame, making room for the new model data)
    CALL write_time_to_file( filename, region%time)

    ! Write model data
    ! ================

    CALL write_help_field_mesh( filename, region, C%help_field_01)
    CALL write_help_field_mesh( filename, region, C%help_field_02)
    CALL write_help_field_mesh( filename, region, C%help_field_03)
    CALL write_help_field_mesh( filename, region, C%help_field_04)
    CALL write_help_field_mesh( filename, region, C%help_field_05)
    CALL write_help_field_mesh( filename, region, C%help_field_06)
    CALL write_help_field_mesh( filename, region, C%help_field_07)
    CALL write_help_field_mesh( filename, region, C%help_field_08)
    CALL write_help_field_mesh( filename, region, C%help_field_09)
    CALL write_help_field_mesh( filename, region, C%help_field_10)
    CALL write_help_field_mesh( filename, region, C%help_field_11)
    CALL write_help_field_mesh( filename, region, C%help_field_12)
    CALL write_help_field_mesh( filename, region, C%help_field_13)
    CALL write_help_field_mesh( filename, region, C%help_field_14)
    CALL write_help_field_mesh( filename, region, C%help_field_15)
    CALL write_help_field_mesh( filename, region, C%help_field_16)
    CALL write_help_field_mesh( filename, region, C%help_field_17)
    CALL write_help_field_mesh( filename, region, C%help_field_18)
    CALL write_help_field_mesh( filename, region, C%help_field_19)
    CALL write_help_field_mesh( filename, region, C%help_field_20)
    CALL write_help_field_mesh( filename, region, C%help_field_21)
    CALL write_help_field_mesh( filename, region, C%help_field_22)
    CALL write_help_field_mesh( filename, region, C%help_field_23)
    CALL write_help_field_mesh( filename, region, C%help_field_24)
    CALL write_help_field_mesh( filename, region, C%help_field_25)
    CALL write_help_field_mesh( filename, region, C%help_field_26)
    CALL write_help_field_mesh( filename, region, C%help_field_27)
    CALL write_help_field_mesh( filename, region, C%help_field_28)
    CALL write_help_field_mesh( filename, region, C%help_field_29)
    CALL write_help_field_mesh( filename, region, C%help_field_30)
    CALL write_help_field_mesh( filename, region, C%help_field_31)
    CALL write_help_field_mesh( filename, region, C%help_field_32)
    CALL write_help_field_mesh( filename, region, C%help_field_33)
    CALL write_help_field_mesh( filename, region, C%help_field_34)
    CALL write_help_field_mesh( filename, region, C%help_field_35)
    CALL write_help_field_mesh( filename, region, C%help_field_36)
    CALL write_help_field_mesh( filename, region, C%help_field_37)
    CALL write_help_field_mesh( filename, region, C%help_field_38)
    CALL write_help_field_mesh( filename, region, C%help_field_39)
    CALL write_help_field_mesh( filename, region, C%help_field_40)
    CALL write_help_field_mesh( filename, region, C%help_field_41)
    CALL write_help_field_mesh( filename, region, C%help_field_42)
    CALL write_help_field_mesh( filename, region, C%help_field_43)
    CALL write_help_field_mesh( filename, region, C%help_field_44)
    CALL write_help_field_mesh( filename, region, C%help_field_45)
    CALL write_help_field_mesh( filename, region, C%help_field_46)
    CALL write_help_field_mesh( filename, region, C%help_field_47)
    CALL write_help_field_mesh( filename, region, C%help_field_48)
    CALL write_help_field_mesh( filename, region, C%help_field_49)
    CALL write_help_field_mesh( filename, region, C%help_field_50)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_help_fields_file_mesh

  SUBROUTINE write_to_help_fields_file_grid( filename, region)
    ! Write model output to the mesh help fields file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_help_fields_file_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write new time to file (thus extending the time dimension by one frame, making room for the new model data)
    CALL write_time_to_file( filename, region%time)

    ! Write model data
    ! ================

    CALL write_help_field_grid( filename, region, C%help_field_01)
    CALL write_help_field_grid( filename, region, C%help_field_02)
    CALL write_help_field_grid( filename, region, C%help_field_03)
    CALL write_help_field_grid( filename, region, C%help_field_04)
    CALL write_help_field_grid( filename, region, C%help_field_05)
    CALL write_help_field_grid( filename, region, C%help_field_06)
    CALL write_help_field_grid( filename, region, C%help_field_07)
    CALL write_help_field_grid( filename, region, C%help_field_08)
    CALL write_help_field_grid( filename, region, C%help_field_09)
    CALL write_help_field_grid( filename, region, C%help_field_10)
    CALL write_help_field_grid( filename, region, C%help_field_11)
    CALL write_help_field_grid( filename, region, C%help_field_12)
    CALL write_help_field_grid( filename, region, C%help_field_13)
    CALL write_help_field_grid( filename, region, C%help_field_14)
    CALL write_help_field_grid( filename, region, C%help_field_15)
    CALL write_help_field_grid( filename, region, C%help_field_16)
    CALL write_help_field_grid( filename, region, C%help_field_17)
    CALL write_help_field_grid( filename, region, C%help_field_18)
    CALL write_help_field_grid( filename, region, C%help_field_19)
    CALL write_help_field_grid( filename, region, C%help_field_20)
    CALL write_help_field_grid( filename, region, C%help_field_21)
    CALL write_help_field_grid( filename, region, C%help_field_22)
    CALL write_help_field_grid( filename, region, C%help_field_23)
    CALL write_help_field_grid( filename, region, C%help_field_24)
    CALL write_help_field_grid( filename, region, C%help_field_25)
    CALL write_help_field_grid( filename, region, C%help_field_26)
    CALL write_help_field_grid( filename, region, C%help_field_27)
    CALL write_help_field_grid( filename, region, C%help_field_28)
    CALL write_help_field_grid( filename, region, C%help_field_29)
    CALL write_help_field_grid( filename, region, C%help_field_30)
    CALL write_help_field_grid( filename, region, C%help_field_31)
    CALL write_help_field_grid( filename, region, C%help_field_32)
    CALL write_help_field_grid( filename, region, C%help_field_33)
    CALL write_help_field_grid( filename, region, C%help_field_34)
    CALL write_help_field_grid( filename, region, C%help_field_35)
    CALL write_help_field_grid( filename, region, C%help_field_36)
    CALL write_help_field_grid( filename, region, C%help_field_37)
    CALL write_help_field_grid( filename, region, C%help_field_38)
    CALL write_help_field_grid( filename, region, C%help_field_39)
    CALL write_help_field_grid( filename, region, C%help_field_40)
    CALL write_help_field_grid( filename, region, C%help_field_41)
    CALL write_help_field_grid( filename, region, C%help_field_42)
    CALL write_help_field_grid( filename, region, C%help_field_43)
    CALL write_help_field_grid( filename, region, C%help_field_44)
    CALL write_help_field_grid( filename, region, C%help_field_45)
    CALL write_help_field_grid( filename, region, C%help_field_46)
    CALL write_help_field_grid( filename, region, C%help_field_47)
    CALL write_help_field_grid( filename, region, C%help_field_48)
    CALL write_help_field_grid( filename, region, C%help_field_49)
    CALL write_help_field_grid( filename, region, C%help_field_50)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_help_fields_file_grid

  SUBROUTINE write_help_field_mesh( filename, region, field_name)
    ! Write a single data field to the help fields file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_model_region),             INTENT(INOUT) :: region
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_help_field_mesh'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    ELSEIF (field_name == 'lon') THEN
      CALL warning('longitude is already by default written to output!')
    ELSEIF (field_name == 'lat') THEN
      CALL warning('latitude is already by default written to output!')

    ! Mesh data
    ELSEIF (field_name == 'resolution') THEN
      CALL warning('mesh resolution is already by default written to output!')

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_notime( filename, region%mesh, 'GHF', region%ice%GHF_a)

    ! Fields with a time dimension
    ! ============================

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Hi', region%ice%Hi_a)
    ELSEIF (field_name == 'Hb') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Hb', region%ice%Hb_a)
    ELSEIF (field_name == 'Hs') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Hs', region%ice%Hs_a)
    ELSEIF (field_name == 'SL') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'SL', region%ice%SL_a)
    ELSEIF (field_name == 'dHi') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'dHi', region%ice%dHi_a)
    ELSEIF (field_name == 'dHb') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'dHb', region%ice%dHb_a)
    ELSEIF (field_name == 'dHs') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'dHs', region%ice%dHs_a)

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL write_to_field_multiple_options_mesh_dp_3D( filename, region%mesh, 'Ti', region%ice%Ti_a)
    ELSEIF (field_name == 'Cpi') THEN
      CALL write_to_field_multiple_options_mesh_dp_3D( filename, region%mesh, 'Cpi', region%ice%Cpi_a)
    ELSEIF (field_name == 'Ki') THEN
      CALL write_to_field_multiple_options_mesh_dp_3D( filename, region%mesh, 'Ki', region%ice%Ki_a)
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Ti_basal', region%ice%Ti_a( :,C%nz))
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL write_to_field_multiple_options_mesh_dp_3D( filename, region%mesh, 'Ti_pmp', region%ice%Ti_pmp_a)
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL write_to_field_multiple_options_mesh_dp_3D( filename, region%mesh, 'A_flow_3D', region%ice%A_flow_3D_a)
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'A_flow_vav', region%ice%A_flow_vav_a)

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL write_to_field_multiple_options_mesh_dp_3D( filename, region%mesh, 'u_3D', region%ice%u_3D_a)
    ELSEIF (field_name == 'v_3D') THEN
      CALL write_to_field_multiple_options_mesh_dp_3D( filename, region%mesh, 'v_3D', region%ice%v_3D_a)
    ELSEIF (field_name == 'w_3D') THEN
      CALL write_to_field_multiple_options_mesh_dp_3D( filename, region%mesh, 'w_3D', region%ice%w_3D_a)
    ELSEIF (field_name == 'u_vav') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'u_vav', region%ice%u_vav_a)
    ELSEIF (field_name == 'v_vav') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'v_vav', region%ice%v_vav_a)
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'uabs_vav', region%ice%uabs_vav_a)
    ELSEIF (field_name == 'u_surf') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'u_surf', region%ice%u_surf_a)
    ELSEIF (field_name == 'v_surf') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'v_surf', region%ice%v_surf_a)
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'uabs_surf', region%ice%uabs_surf_a)
    ELSEIF (field_name == 'u_base') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'u_base', region%ice%u_base_a)
    ELSEIF (field_name == 'v_base') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'v_base', region%ice%v_base_a)
    ELSEIF (field_name == 'uabs_base') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'uabs_base', region%ice%uabs_base_a)

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'T2m', region%climate_matrix%applied%T2m)
    ELSEIF (field_name == 'T2m_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'T2m_year', SUM( region%climate_matrix%applied%T2m,2) / 12._dp)
    ELSEIF (field_name == 'Precip') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Precip', region%climate_matrix%applied%Precip)
    ELSEIF (field_name == 'Precip_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Precip_year', SUM( region%climate_matrix%applied%Precip,2) )
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Wind_WE', region%climate_matrix%applied%Wind_WE)
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Wind_WE_year', SUM( region%climate_matrix%applied%Wind_WE,2) / 12._dp)
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Wind_SN', region%climate_matrix%applied%Wind_SN)
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Wind_SN_year', SUM( region%climate_matrix%applied%Wind_SN,2) / 12._dp)
    ELSEIF (field_name == 'Wind_LR') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Wind_LR', region%climate_matrix%applied%Wind_LR)
    ELSEIF (field_name == 'Wind_LR_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Wind_LR_year', SUM( region%climate_matrix%applied%Wind_LR,2) / 12._dp)
    ELSEIF (field_name == 'Wind_DU') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Wind_DU', region%climate_matrix%applied%Wind_DU)
    ELSEIF (field_name == 'Wind_DU_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Wind_DU_year', SUM( region%climate_matrix%applied%Wind_DU,2) / 12._dp)

    ! Surface mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'SMB', region%SMB%SMB)
    ELSEIF (field_name == 'SMB_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'SMB_year', SUM( region%SMB%SMB,2))
    ELSEIF (field_name == 'Snowfall') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Snowfall', region%SMB%Snowfall)
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Snowfall_year', SUM( region%SMB%Snowfall,2))
    ELSEIF (field_name == 'Rainfall') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Rainfall', region%SMB%Rainfall)
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Rainfall_year', SUM( region%SMB%Rainfall,2))
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'AddedFirn', region%SMB%AddedFirn)
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'AddedFirn_year', SUM( region%SMB%AddedFirn,2))
    ELSEIF (field_name == 'Refreezing') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Refreezing', region%SMB%Refreezing)
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Refreezing_year', SUM( region%SMB%Refreezing,2))
    ELSEIF (field_name == 'Melt') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Melt', region%SMB%Melt)
    ELSEIF (field_name == 'Melt_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Melt_year', SUM( region%SMB%Melt,2))
    ELSEIF (field_name == 'Runoff') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Runoff', region%SMB%Runoff)
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Runoff_year', SUM( region%SMB%Runoff,2))
    ELSEIF (field_name == 'Albedo') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'Albedo', region%SMB%Albedo)
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'Albedo_year', SUM( region%SMB%Albedo,2) / 12._dp)
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D_monthly( filename, region%mesh, 'FirnDepth', region%SMB%FirnDepth)
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'FirnDepth_year', SUM( region%SMB%FirnDepth,2) / 12._dp)

    ! Basal mass balance
    ELSEIF (field_name == 'BMB') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'BMB', region%BMB%BMB)
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'BMB_sheet', region%BMB%BMB_sheet)
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'BMB_shelf', region%BMB%BMB_shelf)

    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask', region%ice%mask_a)
    ELSEIF (field_name == 'mask_land') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask_land', region%ice%mask_land_a)
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask_ocean', region%ice%mask_ocean_a)
    ELSEIF (field_name == 'mask_lake') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask_lake', region%ice%mask_lake_a)
    ELSEIF (field_name == 'mask_ice') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask_ice', region%ice%mask_ice_a)
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask_sheet', region%ice%mask_sheet_a)
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask_shelf', region%ice%mask_shelf_a)
    ELSEIF (field_name == 'mask_coast') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask_coast', region%ice%mask_coast_a)
    ELSEIF (field_name == 'mask_margin') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask_margin', region%ice%mask_margin_a)
    ELSEIF (field_name == 'mask_gl') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask_gl', region%ice%mask_gl_a)
    ELSEIF (field_name == 'mask_cf') THEN
      CALL write_to_field_multiple_options_mesh_int_2D( filename, region%mesh, 'mask_cf', region%ice%mask_cf_a)

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'phi_fric', region%ice%phi_fric_a)
    ELSEIF (field_name == 'tau_yield') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'tau_yield', region%ice%tauc_a)

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'iso_ice', region%ice%IsoIce)
    ELSEIF (field_name == 'iso_surf') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'iso_surf', region%ice%IsoSurf)

    ! GIA
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'dHb_dt', region%ice%dHb_dt_a)
    ELSEIF (field_name == 'dSL_dt') THEN
      CALL write_to_field_multiple_options_mesh_dp_2D( filename, region%mesh, 'dSL_dt', region%ice%dSL_dt_a)

    ! Safety
    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_help_field_mesh

  SUBROUTINE write_help_field_grid( filename, region, field_name)
    ! Write a single data field to the help fields file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_model_region),             INTENT(INOUT) :: region
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_help_field_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    ELSEIF (field_name == 'lon') THEN
      CALL warning('longitude is already by default written to output!')
    ELSEIF (field_name == 'lat') THEN
      CALL warning('latitude is already by default written to output!')

    ! Mesh data
    ELSEIF (field_name == 'resolution') THEN
      CALL warning('mesh resolution is already by default written to output!')

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_notime( filename, region%mesh, region%name, 'GHF', region%ice%GHF_a)

    ! Fields with a time dimension
    ! ============================

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Hi', region%ice%Hi_a)
    ELSEIF (field_name == 'Hb') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Hb', region%ice%Hb_a)
    ELSEIF (field_name == 'Hs') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Hs', region%ice%Hs_a)
    ELSEIF (field_name == 'SL') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'SL', region%ice%SL_a)
    ELSEIF (field_name == 'dHi') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'dHi', region%ice%dHi_a)
    ELSEIF (field_name == 'dHb') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'dHb', region%ice%dHb_a)
    ELSEIF (field_name == 'dHs') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'dHs', region%ice%dHs_a)

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%mesh, region%name, 'Ti', region%ice%Ti_a)
    ELSEIF (field_name == 'Cpi') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%mesh, region%name, 'Cpi', region%ice%Cpi_a)
    ELSEIF (field_name == 'Ki') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%mesh, region%name, 'Ki', region%ice%Ki_a)
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Ti_basal', region%ice%Ti_a( :,C%nz))
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%mesh, region%name, 'Ti_pmp', region%ice%Ti_pmp_a)
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%mesh, region%name, 'A_flow_3D', region%ice%A_flow_3D_a)
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'A_flow_vav', region%ice%A_flow_vav_a)

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%mesh, region%name, 'u_3D', region%ice%u_3D_a)
    ELSEIF (field_name == 'v_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%mesh, region%name, 'v_3D', region%ice%v_3D_a)
    ELSEIF (field_name == 'w_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%mesh, region%name, 'w_3D', region%ice%w_3D_a)
    ELSEIF (field_name == 'u_vav') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'u_vav', region%ice%u_vav_a)
    ELSEIF (field_name == 'v_vav') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'v_vav', region%ice%v_vav_a)
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'uabs_vav', region%ice%uabs_vav_a)
    ELSEIF (field_name == 'u_surf') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'u_surf', region%ice%u_surf_a)
    ELSEIF (field_name == 'v_surf') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'v_surf', region%ice%v_surf_a)
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'uabs_surf', region%ice%uabs_surf_a)
    ELSEIF (field_name == 'u_base') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'u_base', region%ice%u_base_a)
    ELSEIF (field_name == 'v_base') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'v_base', region%ice%v_base_a)
    ELSEIF (field_name == 'uabs_base') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'uabs_base', region%ice%uabs_base_a)

    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'T2m', region%climate_matrix%applied%T2m)
    ELSEIF (field_name == 'T2m_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'T2m_year', SUM( region%climate_matrix%applied%T2m,2) / 12._dp)
    ELSEIF (field_name == 'Precip') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Precip', region%climate_matrix%applied%Precip)
    ELSEIF (field_name == 'Precip_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Precip_year', SUM( region%climate_matrix%applied%Precip,2) )
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Wind_WE', region%climate_matrix%applied%Wind_WE)
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Wind_WE_year', SUM( region%climate_matrix%applied%Wind_WE,2) / 12._dp)
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Wind_SN', region%climate_matrix%applied%Wind_SN)
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Wind_SN_year', SUM( region%climate_matrix%applied%Wind_SN,2) / 12._dp)
    ELSEIF (field_name == 'Wind_LR') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Wind_LR', region%climate_matrix%applied%Wind_LR)
    ELSEIF (field_name == 'Wind_LR_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Wind_LR_year', SUM( region%climate_matrix%applied%Wind_LR,2) / 12._dp)
    ELSEIF (field_name == 'Wind_DU') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Wind_DU', region%climate_matrix%applied%Wind_DU)
    ELSEIF (field_name == 'Wind_DU_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Wind_DU_year', SUM( region%climate_matrix%applied%Wind_DU,2) / 12._dp)

    ! Surface mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'SMB', region%SMB%SMB)
    ELSEIF (field_name == 'SMB_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'SMB_year', SUM( region%SMB%SMB,2))
    ELSEIF (field_name == 'Snowfall') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Snowfall', region%SMB%Snowfall)
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Snowfall_year', SUM( region%SMB%Snowfall,2))
    ELSEIF (field_name == 'Rainfall') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Rainfall', region%SMB%Rainfall)
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Rainfall_year', SUM( region%SMB%Rainfall,2))
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'AddedFirn', region%SMB%AddedFirn)
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'AddedFirn_year', SUM( region%SMB%AddedFirn,2))
    ELSEIF (field_name == 'Refreezing') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Refreezing', region%SMB%Refreezing)
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Refreezing_year', SUM( region%SMB%Refreezing,2))
    ELSEIF (field_name == 'Melt') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Melt', region%SMB%Melt)
    ELSEIF (field_name == 'Melt_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Melt_year', SUM( region%SMB%Melt,2))
    ELSEIF (field_name == 'Runoff') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Runoff', region%SMB%Runoff)
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Runoff_year', SUM( region%SMB%Runoff,2))
    ELSEIF (field_name == 'Albedo') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'Albedo', region%SMB%Albedo)
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'Albedo_year', SUM( region%SMB%Albedo,2) / 12._dp)
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%mesh, region%name, 'FirnDepth', region%SMB%FirnDepth)
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'FirnDepth_year', SUM( region%SMB%FirnDepth,2) / 12._dp)

    ! Basal mass balance
    ELSEIF (field_name == 'BMB') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'BMB', region%BMB%BMB)
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'BMB_sheet', region%BMB%BMB_sheet)
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'BMB_shelf', region%BMB%BMB_shelf)

    ! Masks
    ! Not used, as mapping masks from the mesh to the grid is meaningless
    ELSEIF (field_name == 'mask') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask', region%ice%mask_a)
    ELSEIF (field_name == 'mask_land') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask_land', region%ice%mask_land_a)
    ELSEIF (field_name == 'mask_ocean') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask_ocean', region%ice%mask_ocean_a)
    ELSEIF (field_name == 'mask_lake') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask_lake', region%ice%mask_lake_a)
    ELSEIF (field_name == 'mask_ice') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask_ice', region%ice%mask_ice_a)
    ELSEIF (field_name == 'mask_sheet') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask_sheet', region%ice%mask_sheet_a)
    ELSEIF (field_name == 'mask_shelf') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask_shelf', region%ice%mask_shelf_a)
    ELSEIF (field_name == 'mask_coast') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask_coast', region%ice%mask_coast_a)
    ELSEIF (field_name == 'mask_margin') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask_margin', region%ice%mask_margin_a)
    ELSEIF (field_name == 'mask_gl') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask_gl', region%ice%mask_gl_a)
    ELSEIF (field_name == 'mask_cf') THEN
!      CALL write_to_field_multiple_options_grid_int_2D( filename, region%mesh, region%name, 'mask_cf', region%ice%mask_cf_a)

    ! Basal conditions
    ELSEIF (field_name == 'phi_fric') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'phi_fric', region%ice%phi_fric_a)
    ELSEIF (field_name == 'tau_yield') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'tau_yield', region%ice%tauc_a)

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'iso_ice', region%ice%IsoIce)
    ELSEIF (field_name == 'iso_surf') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'iso_surf', region%ice%IsoSurf)

    ! GIA
    ELSEIF (field_name == 'dHb_dt') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'dHb_dt', region%ice%dHb_dt_a)
    ELSEIF (field_name == 'dSL_dt') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%mesh, region%name, 'dSL_dt', region%ice%dSL_dt_a)

    ! Safety
    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_help_field_grid

  ! ==== Write data to flexibly-defined fields =====
  ! ================================================

  ! Write data to a grid output file
  SUBROUTINE write_to_field_multiple_options_grid_dp_2D( filename, mesh, region_name, field_name_options, d)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on an x/y-grid
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    TYPE(type_grid)                                    :: grid
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_grid
    INTEGER                                            :: wd_grid
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_with_time
    INTEGER                                            :: wd_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D( filename, var_name, should_have_time = .TRUE.)

    ! Set up grid from file
    CALL setup_xy_grid_from_file( filename, grid, region_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid, wd_grid)

    ! Calculate mapping operator between the mesh and the grid
    CALL calc_remapping_operator_mesh2grid( mesh, grid)

    ! Map data from mesh to grid
    CALL map_mesh2grid_2D( mesh, grid, d, d_grid)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, 1, d_grid_with_time, wd_grid_with_time)

    ! Copy data
    d_grid_with_time( grid%i1:grid%i2,:,1) = d_grid( grid%i1:grid%i2,:)
    CALL sync

    ! Write data to the variable
    CALL write_var_dp_3D( filename, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid%nx, grid%ny, 1 /) )

    ! Clean up after yourself
    CALL deallocate_grid( grid)
    CALL deallocate_remapping_operators_mesh2grid( grid)
    CALL deallocate_shared( wd_grid)
    CALL deallocate_shared( wd_grid_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_2D

  SUBROUTINE write_to_field_multiple_options_grid_dp_2D_monthly( filename, mesh, region_name, field_name_options, d)
    ! Write a 2-D monthly data field defined on a mesh to a NetCDF file variable on an x/y-grid
    ! (Mind you, that's 2-D monthly in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_2D_monthly'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    TYPE(type_grid)                                    :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid
    INTEGER                                            :: wd_grid
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_grid_with_time
    INTEGER                                            :: wd_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D_monthly( filename, var_name, should_have_time = .TRUE.)

    ! Set up grid from file
    CALL setup_xy_grid_from_file( filename, grid, region_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, 12, d_grid, wd_grid)

    ! Calculate mapping operator between the mesh and the grid
    CALL calc_remapping_operator_mesh2grid( mesh, grid)

    ! Map data from mesh to grid
    CALL map_mesh2grid_3D( mesh, grid, d, d_grid)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_4D( grid%nx, grid%ny, 12, 1, d_grid_with_time, wd_grid_with_time)

    ! Copy data
    d_grid_with_time( grid%i1:grid%i2,:,:,1) = d_grid( grid%i1:grid%i2,:,:)
    CALL sync

    ! Write data to the variable
    CALL write_var_dp_4D( filename, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, 12, 1 /) )

    ! Clean up after yourself
    CALL deallocate_grid( grid)
    CALL deallocate_remapping_operators_mesh2grid( grid)
    CALL deallocate_shared( wd_grid)
    CALL deallocate_shared( wd_grid_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_2D_monthly

  SUBROUTINE write_to_field_multiple_options_grid_dp_3D( filename, mesh, region_name, field_name_options, d)
    ! Write a 3-D data field defined on a mesh to a NetCDF file variable on an x/y-grid
    ! (Mind you, that's 3-D in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_3D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    TYPE(type_grid)                                    :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid
    INTEGER                                            :: wd_grid
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_grid_with_time
    INTEGER                                            :: wd_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_3D( filename, var_name, should_have_time = .TRUE.)

    ! Set up grid from file
    CALL setup_xy_grid_from_file( filename, grid, region_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, C%nz, d_grid, wd_grid)

    ! Calculate mapping operator between the mesh and the grid
    CALL calc_remapping_operator_mesh2grid( mesh, grid)

    ! Map data from mesh to grid
    CALL map_mesh2grid_3D( mesh, grid, d, d_grid)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_4D( grid%nx, grid%ny, C%nz, 1, d_grid_with_time, wd_grid_with_time)

    ! Copy data
    d_grid_with_time( grid%i1:grid%i2,:,:,1) = d_grid( grid%i1:grid%i2,:,:)
    CALL sync

    ! Write data to the variable
    CALL write_var_dp_4D( filename, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, C%nz, 1 /) )

    ! Clean up after yourself
    CALL deallocate_grid( grid)
    CALL deallocate_remapping_operators_mesh2grid( grid)
    CALL deallocate_shared( wd_grid)
    CALL deallocate_shared( wd_grid_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_3D

  SUBROUTINE write_to_field_multiple_options_grid_dp_2D_notime( filename, mesh, region_name, field_name_options, d)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on an x/y-grid
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_2D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    TYPE(type_grid)                                    :: grid
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D( filename, var_name, should_have_time = .FALSE.)

    ! Set up grid from file
    CALL setup_xy_grid_from_file( filename, grid, region_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid, wd_grid)

    ! Calculate mapping operator between the mesh and the grid
    CALL calc_remapping_operator_mesh2grid( mesh, grid)

    ! Map data from mesh to grid
    CALL map_mesh2grid_2D( mesh, grid, d, d_grid)

    ! Write data to the variable
    CALL write_var_dp_2D( filename, id_var, d_grid)

    ! Clean up after yourself
    CALL deallocate_grid( grid)
    CALL deallocate_remapping_operators_mesh2grid( grid)
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_2D_notime

  SUBROUTINE write_to_field_multiple_options_grid_dp_2D_monthly_notime( filename, mesh, region_name, field_name_options, d)
    ! Write a 2-D monthly data field defined on a mesh to a NetCDF file variable on an x/y-grid
    ! (Mind you, that's 2-D monthly in the physical sense, so a 1-D array!)
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_2D_monthly_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    TYPE(type_grid)                                    :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D_monthly( filename, var_name, should_have_time = .FALSE.)

    ! Set up grid from file
    CALL setup_xy_grid_from_file( filename, grid, region_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, 12, d_grid, wd_grid)

    ! Calculate mapping operator between the mesh and the grid
    CALL calc_remapping_operator_mesh2grid( mesh, grid)

    ! Map data from mesh to grid
    CALL map_mesh2grid_3D( mesh, grid, d, d_grid)

    ! Write data to the variable
    CALL write_var_dp_3D( filename, id_var, d_grid)

    ! Clean up after yourself
    CALL deallocate_grid( grid)
    CALL deallocate_remapping_operators_mesh2grid( grid)
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_2D_monthly_notime

  SUBROUTINE write_to_field_multiple_options_grid_dp_3D_notime( filename, mesh, region_name, field_name_options, d)
    ! Write a 3-D data field defined on a mesh to a NetCDF file variable on an x/y-grid
    ! (Mind you, that's 3-D in the physical sense, so a 1-D array!)
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_3D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    TYPE(type_grid)                                    :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_3D( filename, var_name, should_have_time = .FALSE.)

    ! Set up grid from file
    CALL setup_xy_grid_from_file( filename, grid, region_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, C%nz, d_grid, wd_grid)

    ! Calculate mapping operator between the mesh and the grid
    CALL calc_remapping_operator_mesh2grid( mesh, grid)

    ! Map data from mesh to grid
    CALL map_mesh2grid_3D( mesh, grid, d, d_grid)

    ! Write data to the variable
    CALL write_var_dp_3D( filename, id_var, d_grid)

    ! Clean up after yourself
    CALL deallocate_grid( grid)
    CALL deallocate_remapping_operators_mesh2grid( grid)
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_3D_notime

  ! Write data to a mesh output file
  SUBROUTINE write_to_field_multiple_options_mesh_int_2D( filename, mesh, field_name_options, d)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_mesh_int_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_int_2D( filename, var_name, should_have_time = .TRUE.)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_int_2D( mesh%nV, 1, d_with_time, wd_with_time)

    ! Copy data
    d_with_time( mesh%vi1:mesh%vi2,1) = d( mesh%vi1:mesh%vi2)
    CALL sync

    ! Write data to the variable
    CALL write_var_int_2D( filename, id_var, d_with_time, start = (/ 1, ti /), count = (/ mesh%nV, 1 /) )

    ! Clean up after yourself
    CALL deallocate_shared( wd_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_mesh_int_2D

  SUBROUTINE write_to_field_multiple_options_mesh_dp_2D( filename, mesh, field_name_options, d)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_mesh_dp_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D( filename, var_name, should_have_time = .TRUE.)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, 1, d_with_time, wd_with_time)

    ! Copy data
    d_with_time( mesh%vi1:mesh%vi2,1) = d( mesh%vi1:mesh%vi2)
    CALL sync

    ! Write data to the variable
    CALL write_var_dp_2D( filename, id_var, d_with_time, start = (/ 1, ti /), count = (/ mesh%nV, 1 /) )

    ! Clean up after yourself
    CALL deallocate_shared( wd_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_mesh_dp_2D

  SUBROUTINE write_to_field_multiple_options_mesh_dp_2D_monthly( filename, mesh, field_name_options, d)
    ! Write a 2-D monthly data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D monthly in the physical sense, so a 2-D array!)
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_mesh_dp_2D_monthly'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D_monthly( filename, var_name, should_have_time = .TRUE.)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( mesh%nV, 12, 1, d_with_time, wd_with_time)

    ! Copy data
    d_with_time( mesh%vi1:mesh%vi2,:,1) = d( mesh%vi1:mesh%vi2,:)
    CALL sync

    ! Write data to the variable
    CALL write_var_dp_3D( filename, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, 12, 1 /) )

    ! Clean up after yourself
    CALL deallocate_shared( wd_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_mesh_dp_2D_monthly

  SUBROUTINE write_to_field_multiple_options_mesh_dp_3D( filename, mesh, field_name_options, d)
    ! Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 3-D in the physical sense, so a 2-D array!)
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_mesh_dp_3D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_3D( filename, var_name, should_have_time = .TRUE.)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( mesh%nV, C%nz, 1, d_with_time, wd_with_time)

    ! Copy data
    d_with_time( mesh%vi1:mesh%vi2,:,1) = d( mesh%vi1:mesh%vi2,:)
    CALL sync

    ! Write data to the variable
    CALL write_var_dp_3D( filename, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ mesh%nV, C%nz, 1 /) )

    ! Clean up after yourself
    CALL deallocate_shared( wd_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_mesh_dp_3D

  SUBROUTINE write_to_field_multiple_options_mesh_int_2D_notime( filename, mesh, field_name_options, d)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_mesh_int_2D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_int_2D( filename, var_name, should_have_time = .FALSE.)

    ! Write data to the variable
    CALL write_var_int_1D( filename, id_var, d)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_mesh_int_2D_notime

  SUBROUTINE write_to_field_multiple_options_mesh_dp_2D_notime( filename, mesh, field_name_options, d)
    ! Write a 2-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_mesh_dp_2D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D( filename, var_name, should_have_time = .FALSE.)

    ! Write data to the variable
    CALL write_var_dp_1D( filename, id_var, d)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_mesh_dp_2D_notime

  SUBROUTINE write_to_field_multiple_options_mesh_dp_2D_monthly_notime( filename, mesh, field_name_options, d)
    ! Write a 2-D monthly data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 2-D monthly in the physical sense, so a 2-D array!)
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_mesh_dp_2D_monthly_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_2D_monthly( filename, var_name, should_have_time = .FALSE.)

    ! Write data to the variable
    CALL write_var_dp_2D( filename, id_var, d)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_mesh_dp_2D_monthly_notime

  SUBROUTINE write_to_field_multiple_options_mesh_dp_3D_notime( filename, mesh, field_name_options, d)
    ! Write a 3-D data field defined on a mesh to a NetCDF file variable on the same mesh
    ! (Mind you, that's 3-D in the physical sense, so a 2-D array!)
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_mesh_dp_3D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_mesh_field_dp_3D( filename, var_name, should_have_time = .FALSE.)

    ! Write data to the variable
    CALL write_var_dp_2D( filename, id_var, d)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_mesh_dp_3D_notime

  ! Write new time value to file
  SUBROUTINE write_time_to_file( filename, time)
    ! Write new time value to file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_time_to_file'
    INTEGER                                            :: id_dim_time
    INTEGER                                            :: id_var_time
    INTEGER                                            :: nt

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check time dimension and variable
    CALL check_time( filename)

    ! Determine current length of time dimension in file
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = nt)

    ! Inquire variable id
    CALL inquire_var_multiple_options( filename, field_name_options_time, id_var_time)

    ! Write time
    nt = nt + 1
    CALL write_var_dp_1D( filename, id_var_time, (/ time /), start = (/ nt /), count = (/ 1 /) )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_time_to_file

  ! ===== Set up mesh/grid, time, and field variables in a NetCDF file =====
  ! ========================================================================

  ! Set up x/y-grid and gridded variables
  SUBROUTINE setup_xy_grid_in_netcdf_file( filename, grid)
    ! Set up a regular x/y-grid in an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid),                     INTENT(IN)    :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_xy_grid_in_netcdf_file'
    INTEGER                                            :: id_dim_x
    INTEGER                                            :: id_dim_y
    INTEGER                                            :: id_var_x
    INTEGER                                            :: id_var_y
    INTEGER                                            :: id_var_lon
    INTEGER                                            :: id_var_lat

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create x/y dimensions
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_x), grid%nx, id_dim_x)
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_y), grid%ny, id_dim_y)

    ! Create and write x/y variables

    ! x
    CALL create_variable( filename, get_first_option_from_list( field_name_options_x), NF90_DOUBLE, (/ id_dim_x /), id_var_x)
    CALL add_attribute_char( filename, id_var_x, 'long_name', 'x-coordinate')
    CALL add_attribute_char( filename, id_var_x, 'units'    , 'm'           )
    CALL write_var_dp_1D( filename, id_var_x, grid%x)
    ! y
    CALL create_variable( filename, get_first_option_from_list( field_name_options_y), NF90_DOUBLE, (/ id_dim_y /), id_var_y)
    CALL add_attribute_char( filename, id_var_y, 'long_name', 'y-coordinate')
    CALL add_attribute_char( filename, id_var_y, 'units'    , 'm'           )
    CALL write_var_dp_1D( filename, id_var_y, grid%y)

    ! Create and write lon/lat variables

    ! lon
    CALL add_field_grid_dp_2D_notime( filename, get_first_option_from_list( field_name_options_lon), long_name = 'Longitude', units = 'degrees east')
    CALL inquire_var_multiple_options( filename, field_name_options_lon, id_var_lon)
    CALL write_var_dp_2D( filename, id_var_lon, grid%lon)
    ! lat
    CALL add_field_grid_dp_2D_notime( filename, get_first_option_from_list( field_name_options_lat), long_name = 'Latitude', units = 'degrees north')
    CALL inquire_var_multiple_options( filename, field_name_options_lat, id_var_lat)
    CALL write_var_dp_2D( filename, id_var_lat, grid%lat)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_xy_grid_in_netcdf_file

  SUBROUTINE add_field_grid_int_2D( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_int_2D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename)
    CALL check_y(    filename)
    CALL check_time( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multiple_options( filename, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_INT, (/ id_dim_x, id_dim_y, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_int_2D( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_int_2D

  SUBROUTINE add_field_grid_dp_2D( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename)
    CALL check_y(    filename)
    CALL check_time( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multiple_options( filename, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D

  SUBROUTINE add_field_grid_dp_2D_monthly( filename, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D_monthly'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_month, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(     filename)
    CALL check_y(     filename)
    CALL check_month( filename)
    CALL check_time(  filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x    , id_dim_x    )
    CALL inquire_dim_multiple_options( filename, field_name_options_y    , id_dim_y    )
    CALL inquire_dim_multiple_options( filename, field_name_options_month, id_dim_month)
    CALL inquire_dim_multiple_options( filename, field_name_options_time , id_dim_time )

    ! Safety
    IF (id_dim_x     == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y     == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time  == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_month, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D_monthly( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D_monthly

  SUBROUTINE add_field_grid_dp_3D( filename, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_3D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_zeta, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename)
    CALL check_y(    filename)
    CALL check_zeta( filename)
    CALL check_time( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multiple_options( filename, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta, id_dim_zeta)
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_zeta, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_3D( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_3D

  SUBROUTINE add_field_grid_int_2D_notime( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_int_2D_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y dimensions and variables are there
    CALL check_x( filename)
    CALL check_y( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x, id_dim_x)
    CALL inquire_dim_multiple_options( filename, field_name_options_y, id_dim_y)

    ! Safety
    IF (id_dim_x == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_INT, (/ id_dim_x, id_dim_y /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_int_2D( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_int_2D_notime

  SUBROUTINE add_field_grid_dp_2D_notime( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x( filename)
    CALL check_y( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x, id_dim_x)
    CALL inquire_dim_multiple_options( filename, field_name_options_y, id_dim_y)

    ! Safety
    IF (id_dim_x == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D_notime

  SUBROUTINE add_field_grid_dp_2D_monthly_notime( filename, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D_monthly_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_month, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(     filename)
    CALL check_y(     filename)
    CALL check_month( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x    , id_dim_x    )
    CALL inquire_dim_multiple_options( filename, field_name_options_y    , id_dim_y    )
    CALL inquire_dim_multiple_options( filename, field_name_options_month, id_dim_month)

    ! Safety
    IF (id_dim_x     == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y     == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_month /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D_monthly( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D_monthly_notime

  SUBROUTINE add_field_grid_dp_3D_notime( filename, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_3D_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_zeta, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename)
    CALL check_y(    filename)
    CALL check_zeta( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multiple_options( filename, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta, id_dim_zeta)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_zeta /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_3D( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_3D_notime

  ! Set up mesh and meshed variables
  SUBROUTINE setup_mesh_in_netcdf_file( filename, mesh)
    ! Set up a mesh in an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_mesh_in_netcdf_file'
    INTEGER                                            :: id_dim_vi
    INTEGER                                            :: id_dim_ti
    INTEGER                                            :: id_dim_ci
    INTEGER                                            :: id_dim_aci
    INTEGER                                            :: id_dim_two
    INTEGER                                            :: id_dim_three
    INTEGER                                            :: id_dim_six
    INTEGER                                            :: id_var_V
    INTEGER                                            :: id_var_nC
    INTEGER                                            :: id_var_C
    INTEGER                                            :: id_var_niTri
    INTEGER                                            :: id_var_iTri
    INTEGER                                            :: id_var_edge_index
    INTEGER                                            :: id_var_Tri
    INTEGER                                            :: id_var_Tricc
    INTEGER                                            :: id_var_TriC
    INTEGER                                            :: id_var_Tri_edge_index
    INTEGER                                            :: id_var_VAc
    INTEGER                                            :: id_var_Aci
    INTEGER                                            :: id_var_iAci
    INTEGER                                            :: id_var_A
    INTEGER                                            :: id_var_R
    INTEGER                                            :: id_var_lon
    INTEGER                                            :: id_var_lat

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create mesh dimensions
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_dim_nV    ), mesh%nV    , id_dim_vi   )
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_dim_nTri  ), mesh%nTri  , id_dim_ti   )
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_dim_nC_mem), mesh%nC_mem, id_dim_ci   )
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_dim_nAc   ), mesh%nAc   , id_dim_aci  )
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_dim_two   ), 2          , id_dim_two  )
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_dim_three ), 3          , id_dim_three)
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_dim_six   ), 6          , id_dim_six  )

    ! Create and write mesh variables - vertex data

    ! V
    CALL create_variable( filename, get_first_option_from_list( field_name_options_V             ), NF90_DOUBLE, (/ id_dim_vi, id_dim_two   /), id_var_V             )
    CALL add_attribute_char( filename, id_var_V             , 'long_name'  , 'Vertex coordinates'         )
    CALL add_attribute_char( filename, id_var_V             , 'units'      , 'm'                          )
    CALL write_var_dp_2D(    filename, id_var_V             , mesh%V             )
    ! nC
    CALL create_variable( filename, get_first_option_from_list( field_name_options_nC            ), NF90_INT   , (/ id_dim_vi               /), id_var_nC            )
    CALL add_attribute_char( filename, id_var_nC            , 'long_name'  , 'Number of vertex-vertex connections')
    CALL write_var_int_1D(   filename, id_var_nC            , mesh%nC            )
    ! C
    CALL create_variable( filename, get_first_option_from_list( field_name_options_C             ), NF90_INT   , (/ id_dim_vi, id_dim_ci    /), id_var_C             )
    CALL add_attribute_char( filename, id_var_C             , 'long_name'  , 'Vertex-vertex connections')
    CALL add_attribute_char( filename, id_var_C             , 'orientation', 'counter-clockwise'          )
    CALL write_var_int_2D(   filename, id_var_C             , mesh%C             )
    ! niTri
    CALL create_variable( filename, get_first_option_from_list( field_name_options_niTri         ), NF90_INT   , (/ id_dim_vi               /), id_var_niTri         )
    CALL add_attribute_char( filename, id_var_niTri         , 'long_name'  , 'Number of vertex-triangle connections')
    CALL write_var_int_1D(   filename, id_var_niTri         , mesh%niTri         )
    ! iTri
    CALL create_variable( filename, get_first_option_from_list( field_name_options_iTri          ), NF90_INT   , (/ id_dim_vi, id_dim_ci    /), id_var_iTri          )
    CALL add_attribute_char( filename, id_var_iTri          , 'long_name'  , 'Vertex-triangle connections')
    CALL add_attribute_char( filename, id_var_iTri          , 'orientation', 'counter-clockwise'          )
    CALL write_var_int_2D(   filename, id_var_iTri          , mesh%iTri          )
    ! edge_index
    CALL create_variable( filename, get_first_option_from_list( field_name_options_edge_index    ), NF90_INT   , (/ id_dim_vi               /), id_var_edge_index    )
    CALL add_attribute_char( filename, id_var_edge_index    , 'long_name'  , 'Vertex edge index')
    CALL add_attribute_char( filename, id_var_edge_index    , 'orientation', '1 = N, 2 = NE, 3 = E, 4 = SE, 5 = S, 6 = SW, 7 = W, 8 = NW')
    CALL write_var_int_1D(   filename, id_var_edge_index    , mesh%edge_index    )

    ! Create and write mesh variables - triangle data

    ! Tri
    CALL create_variable( filename, get_first_option_from_list( field_name_options_Tri           ), NF90_INT   , (/ id_dim_ti, id_dim_three /), id_var_Tri           )
    CALL add_attribute_char( filename, id_var_Tri           , 'long_name'  , 'Vertex indices per triangle')
    CALL add_attribute_char( filename, id_var_Tri           , 'orientation', 'counter-clockwise'          )
    CALL write_var_int_2D(   filename, id_var_Tri           , mesh%Tri           )
    ! Tricc
    CALL create_variable( filename, get_first_option_from_list( field_name_options_Tricc         ), NF90_DOUBLE, (/ id_dim_ti, id_dim_two   /), id_var_Tricc         )
    CALL add_attribute_char( filename, id_var_Tricc         , 'long_name'  , 'Triangle circumcentre coordinates')
    CALL add_attribute_char( filename, id_var_Tricc         , 'units'      , 'm'                          )
    CALL write_var_dp_2D(    filename, id_var_Tricc         , mesh%Tricc         )
    ! TriC
    CALL create_variable( filename, get_first_option_from_list( field_name_options_TriC          ), NF90_INT   , (/ id_dim_ti, id_dim_three /), id_var_TriC          )
    CALL add_attribute_char( filename, id_var_TriC          , 'long_name'  , 'Triangle-triangle connections')
    CALL add_attribute_char( filename, id_var_TriC          , 'orientation', 'counter-clockwise, opposite from constituent vertices (i.e. first entry is opposite from first vertex)')
    CALL write_var_int_2D(   filename, id_var_TriC          , mesh%TriC          )
    ! Tri_edge_index
    CALL create_variable( filename, get_first_option_from_list( field_name_options_Tri_edge_index), NF90_INT   , (/ id_dim_ti               /), id_var_Tri_edge_index)
    CALL add_attribute_char( filename, id_var_Tri_edge_index, 'long_name'  , 'Triangle edge index')
    CALL add_attribute_char( filename, id_var_Tri_edge_index, 'orientation', '1 = N, 2 = NE, 3 = E, 4 = SE, 5 = S, 6 = SW, 7 = W, 8 = NW')
    CALL write_var_int_1D(   filename, id_var_Tri_edge_index, mesh%Tri_edge_index)

    ! Create and write mesh variables - edge data

    ! VAc
    CALL create_variable( filename, get_first_option_from_list( field_name_options_VAc           ), NF90_DOUBLE, (/ id_dim_aci, id_dim_two  /), id_var_VAc           )
    CALL add_attribute_char( filename, id_var_VAc           , 'long_name'  , 'Staggered vertex coordinates')
    CALL add_attribute_char( filename, id_var_VAc           , 'units'      , 'm'                           )
    CALL write_var_dp_2D(    filename, id_var_VAc           , mesh%VAc           )
    ! Aci
    CALL create_variable( filename, get_first_option_from_list( field_name_options_Aci           ), NF90_INT   , (/ id_dim_aci, id_dim_six  /), id_var_Aci           )
    CALL add_attribute_char( filename, id_var_Aci           , 'long_name'  , 'Staggered vertex indices')
    CALL add_attribute_char( filename, id_var_Aci           , 'orientation', 'vertex-from,vertex-to,vertex-left,vertex-right,triangle-left,triangle-right')
    CALL write_var_int_2D(   filename, id_var_Aci           , mesh%Aci           )
    ! iAci
    CALL create_variable( filename, get_first_option_from_list( field_name_options_iAci          ), NF90_INT   , (/ id_dim_vi,  id_dim_ci   /), id_var_iAci          )
    CALL add_attribute_char( filename, id_var_iAci          , 'long_name'  , 'Vertex-edge connections')
    CALL add_attribute_char( filename, id_var_iAci          , 'orientation', 'same as C')
    CALL write_var_int_2D( filename, id_var_iAci          , mesh%iAci          )

    ! Create and write mesh variables - resolution, Voronoi cell area, longitude, and latitude

    ! R
    CALL add_field_mesh_dp_2D_notime( filename, get_first_option_from_list( field_name_options_R             ), long_name = 'Resolution', units = 'm')
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( filename, mesh, field_name_options_R  , mesh%R  )
    ! A
    CALL add_field_mesh_dp_2D_notime( filename, get_first_option_from_list( field_name_options_A             ), long_name = 'Voronoi cell area', units = 'm^2')
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( filename, mesh, field_name_options_A  , mesh%A  )
    ! lon
    CALL add_field_mesh_dp_2D_notime( filename, get_first_option_from_list( field_name_options_lon           ), long_name = 'Longitude', units = 'degrees east')
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( filename, mesh, field_name_options_lon, mesh%lon)
    ! lat
    CALL add_field_mesh_dp_2D_notime( filename, get_first_option_from_list( field_name_options_lat           ), long_name = 'Latitude' , units = 'degrees north')
    CALL write_to_field_multiple_options_mesh_dp_2D_notime( filename, mesh, field_name_options_lat, mesh%lat)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_mesh_in_netcdf_file

  SUBROUTINE add_field_mesh_int_2D( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_int_2D'
    INTEGER                                            :: id_dim_vi, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_dim_nV, id_dim_vi  )
    CALL inquire_dim_multiple_options( filename, field_name_options_time  , id_dim_time)

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_INT, (/ id_dim_vi, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_int_2D( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_int_2D

  SUBROUTINE add_field_mesh_dp_2D( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_2D'
    INTEGER                                            :: id_dim_vi, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_dim_nV, id_dim_vi  )
    CALL inquire_dim_multiple_options( filename, field_name_options_time  , id_dim_time)

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D

  SUBROUTINE add_field_mesh_dp_2D_monthly( filename, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_2D_monthly'
    INTEGER                                            :: id_dim_vi, id_dim_month, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_month(           filename)
    CALL check_mesh_dimensions( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_dim_nV, id_dim_vi   )
    CALL inquire_dim_multiple_options( filename, field_name_options_month , id_dim_month)
    CALL inquire_dim_multiple_options( filename, field_name_options_time  , id_dim_time )

    ! Safety
    IF (id_dim_vi    == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time  == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_month, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D_monthly( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D_monthly

  SUBROUTINE add_field_mesh_dp_3D( filename, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_3D'
    INTEGER                                            :: id_dim_vi, id_dim_zeta, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_zeta(            filename)
    CALL check_mesh_dimensions( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_dim_nV, id_dim_vi  )
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta  , id_dim_zeta)
    CALL inquire_dim_multiple_options( filename, field_name_options_time  , id_dim_time)

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_zeta, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_3D( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_3D

  SUBROUTINE add_field_mesh_int_2D_notime( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_int_2D_notime'
    INTEGER                                            :: id_dim_vi, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_dim_nV, id_dim_vi  )

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_INT, (/ id_dim_vi /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_int_2D( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_int_2D_notime

  SUBROUTINE add_field_mesh_dp_2D_notime( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_2D_notime'
    INTEGER                                            :: id_dim_vi, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_mesh_dimensions( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_dim_nV, id_dim_vi  )

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_vi /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D_notime

  SUBROUTINE add_field_mesh_dp_2D_monthly_notime( filename, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_2D_monthly_notime'
    INTEGER                                            :: id_dim_vi, id_dim_month, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_month(           filename)
    CALL check_mesh_dimensions( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_dim_nV, id_dim_vi   )
    CALL inquire_dim_multiple_options( filename, field_name_options_month , id_dim_month)

    ! Safety
    IF (id_dim_vi    == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_month /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_2D_monthly( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_2D_monthly_notime

  SUBROUTINE add_field_mesh_dp_3D_notime( filename, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with a mesh

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_mesh_dp_3D_notime'
    INTEGER                                            :: id_dim_vi, id_dim_zeta, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if all mesh dimensions and variables are there
    CALL check_zeta(            filename)
    CALL check_mesh_dimensions( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_dim_nV, id_dim_vi  )
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta  , id_dim_zeta)

    ! Safety
    IF (id_dim_vi   == -1) CALL crash('no vi dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_vi, id_dim_zeta /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_mesh_field_dp_3D( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_mesh_dp_3D_notime

  ! Add extra dimensions
  SUBROUTINE add_time_dimension_to_file( filename)
    ! Add a time dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_time_dimension_to_file'
    INTEGER                                            :: id_dim_time
    INTEGER                                            :: id_var_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create time dimension
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_time), NF90_UNLIMITED, id_dim_time)

    ! Create time variable
    CALL create_variable(  filename, get_first_option_from_list( field_name_options_time), NF90_DOUBLE, (/ id_dim_time /), id_var_time)
    CALL add_attribute_char( filename, id_var_time, 'long_name', 'Time')
    CALL add_attribute_char( filename, id_var_time, 'units', 'years')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_time_dimension_to_file

  SUBROUTINE add_month_dimension_to_file( filename)
    ! Add a month dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_month_dimension_to_file'
    INTEGER                                            :: id_dim_month
    INTEGER                                            :: id_var_month

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create month dimension
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_month), 12, id_dim_month)

    ! Create month variable
    CALL create_variable(  filename, get_first_option_from_list( field_name_options_month), NF90_INT, (/ id_dim_month /), id_var_month)
    CALL add_attribute_char( filename, id_var_month, 'long_name', 'Month')
    CALL add_attribute_char( filename, id_var_month, 'units', '1-12')
    CALL add_attribute_char( filename, id_var_month, 'description', '1 = Jan, 2 = Feb, ..., 12 = Dec')

    ! Write month variable
    CALL write_var_int_1D( filename, id_var_month, (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /) )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_month_dimension_to_file

  SUBROUTINE add_zeta_dimension_to_file( filename)
    ! Add a zeta dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_zeta_dimension_to_file'
    INTEGER                                            :: id_dim_zeta
    INTEGER                                            :: id_var_zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create month dimension
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_zeta), C%nz, id_dim_zeta)

    ! Create month variable
    CALL create_variable(  filename, get_first_option_from_list( field_name_options_zeta), NF90_DOUBLE, (/ id_dim_zeta /), id_var_zeta)
    CALL add_attribute_char( filename, id_var_zeta, 'long_name', 'Scaled vertical coordinate')
    CALL add_attribute_char( filename, id_var_zeta, 'units', '0-1')
    CALL add_attribute_char( filename, id_var_zeta, 'transformation', 'zeta = (h - z) / H; zeta = 0 at the ice surface; zeta = 1 at the ice base')

    ! Write month variable
    CALL write_var_dp_1D( filename, id_var_zeta, C%zeta)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_zeta_dimension_to_file

  ! ===== Tests =====
  ! =================

  SUBROUTINE NetCDF_output_test( region)
    ! A simple test of the new NetCDF functionality

    USE netcdf_input_module, ONLY: setup_xy_grid_from_file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'NetCDF_output_test'
    CHARACTER(LEN=256)                                 :: filename_restart_mesh
    CHARACTER(LEN=256)                                 :: filename_restart_grid
    CHARACTER(LEN=256)                                 :: filename_help_fields_mesh
    CHARACTER(LEN=256)                                 :: filename_help_fields_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('running ' // TRIM( routine_name))

    ! Create files
    filename_restart_mesh     = TRIM( C%output_dir) // 'test_restart_mesh.nc'
    filename_restart_grid     = TRIM( C%output_dir) // 'test_restart_grid.nc'
    filename_help_fields_mesh = TRIM( C%output_dir) // 'test_help_fields_mesh.nc'
    filename_help_fields_grid = TRIM( C%output_dir) // 'test_help_fields_grid.nc'

    CALL create_restart_file_mesh(     TRIM( filename_restart_mesh    ), region%mesh       )
    CALL create_restart_file_grid(     TRIM( filename_restart_grid    ), region%grid_output)
    CALL create_help_fields_file_mesh( TRIM( filename_help_fields_mesh), region%mesh       )
    CALL create_help_fields_file_grid( TRIM( filename_help_fields_grid), region%grid_output)

    ! Write to files
    CALL write_to_restart_file_mesh(     TRIM( filename_restart_mesh    ), region)
    CALL write_to_restart_file_grid(     TRIM( filename_restart_grid    ), region)
    CALL write_to_help_fields_file_mesh( TRIM( filename_help_fields_mesh), region)
    CALL write_to_help_fields_file_grid( TRIM( filename_help_fields_grid), region)

    CALL warning('finished ' // TRIM( routine_name))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE NetCDF_output_test

END MODULE netcdf_output_module