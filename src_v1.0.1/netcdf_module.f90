MODULE netcdf_module
  ! This module defines methods to read and write the current state of the model from or to a NetCDF file.

  USE mpi
  USE netcdf,                      ONLY: nf90_max_var_dims, nf90_create, nf90_close, nf90_clobber, nf90_share, nf90_unlimited , &
                                         nf90_enddef, nf90_put_var, nf90_sync, nf90_def_var, nf90_int, nf90_put_att, nf90_def_dim, &
                                         nf90_open, nf90_write, nf90_inq_dimid, nf90_inquire_dimension, nf90_inquire, nf90_double, &
                                         nf90_inq_varid, nf90_inquire_variable, nf90_get_var, nf90_noerr, nf90_strerror
  USE parallel_module,             ONLY: par, sync
  USE configuration_module,        ONLY: dp, C
  USE data_types_module,           ONLY: type_model_region, type_mesh, type_PD_data_fields, type_forcing_data, &
                                         type_init_data_fields, type_subclimate_global, type_netcdf_output
  IMPLICIT NONE

CONTAINS

  ! Create and write to an output NetCDF file
  SUBROUTINE write_to_output_file( region)
    ! Write the current model state to the existing output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region
    
    WRITE(0,'(A,F8.2,A)') '   t = ', region%time/1e3, ' kyr - writing output...'
    
    ! Open the file for writing
    CALL open_netcdf_file( region%output%filename, region%output%ncid)
        
    ! Time
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_time,             region%time,                    start=(/       region%output%ti/)))
    
    ! Key output
    ! ==========
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Hi,               region%ice%Hi,                  start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Hb,               region%ice%Hb,                  start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Hs,               region%ice%Hs,                  start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_mask,             region%ice%mask,                start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_U_SIA,            region%ice%U_SIA,               start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_V_SIA,            region%ice%V_SIA,               start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_U_SSA,            region%ice%U_SSA,               start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_V_SSA,            region%ice%V_SSA,               start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Ti,               region%ice%Ti,                  start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_T2m,              region%climate%applied%T2m,     start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Precip,           region%climate%applied%Precip,  start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Albedo,           region%SMB%Albedo,              start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_SMB,              region%SMB%SMB,                 start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_BMB,              region%BMB%BMB,                 start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_FirnDepth,        region%SMB%FirnDepth,           start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_MeltPreviousYear, region%SMB%MeltPreviousYear,    start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_dHs_dx,           region%ice%dHs_dx,              start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_dHs_dy,           region%ice%dHs_dy,              start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_D_SIA,            region%ice%D_SIA,               start=(/1,     region%output%ti/)))
    
    ! Optional output
    ! ===============
    
    ! SMB components
    IF (C%save_secondary_data_SMB_components) THEN
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Wind_WE,          region%climate%applied%Wind_WE, start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Wind_SN,          region%climate%applied%Wind_SN, start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Snowfall,         region%SMB%Snowfall,            start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Melt,             region%SMB%Melt,                start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Refreezing,       region%SMB%Refreezing_year,     start=(/1,     region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Runoff,           region%SMB%Runoff,              start=(/1, 1,  region%output%ti/)))
    END IF
    
    ! Ice dynamics
    IF (C%save_secondary_data_ice_dynamics) THEN
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_U_3D,             region%ice%U_3D,                start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_V_3D,             region%ice%V_3D,                start=(/1, 1,  region%output%ti/)))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_W_3D,             region%ice%W_3D,                start=(/1, 1,  region%output%ti/)))
    END IF
    
    ! Close the file
    CALL close_netcdf_file(region%output%ncid)
    
    ! Increase time frame counter
    region%output%ti = region%output%ti + 1
        
  END SUBROUTINE write_to_output_file
  SUBROUTINE create_output_file(   region)
    ! Create a new output NetCDF file, write the current mesh data to it.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    INTEGER                                       :: cerr, ierr
    LOGICAL                                       :: file_exists
    
    ! Get output file name
    CALL GetOutputFilename(region)
    
    ! Set time frame index to 1
    region%output%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(region%output%filename))
    IF(file_exists) THEN
      WRITE(0,'(5A)') 'ERROR: ', TRIM(region%output%filename), ' already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Create netCDF file
!    WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( region%output%filename)
    CALL handle_error(nf90_create(region%output%filename,IOR(nf90_clobber,nf90_share),region%output%ncid))
        
    ! Define dimensions:
    CALL create_dim( region%output%ncid, region%output%name_dim_vi,        region%mesh%nV,     region%output%id_dim_vi        ) ! Vertex indices
    CALL create_dim( region%output%ncid, region%output%name_dim_two,       2,                  region%output%id_dim_two       ) ! 2 (each vertex has an X and Y coordinates)
    CALL create_dim( region%output%ncid, region%output%name_dim_ti,        region%mesh%nTri,   region%output%id_dim_ti        ) ! Triangle indices
    CALL create_dim( region%output%ncid, region%output%name_dim_three,     3,                  region%output%id_dim_three     ) ! 3 (each triangle has three vertices)
    CALL create_dim( region%output%ncid, region%output%name_dim_ci,        region%mesh%nC_mem, region%output%id_dim_ci        ) ! Connection indices
    CALL create_dim( region%output%ncid, region%output%name_dim_ciplusone, region%mesh%nC_mem, region%output%id_dim_ciplusone ) ! connection indices plus one (neighbour function arrays have one more column)
    CALL create_dim( region%output%ncid, region%output%name_dim_zeta,      C%nZ,               region%output%id_dim_zeta      ) ! Scaled vertical coordinate
    CALL create_dim( region%output%ncid, region%output%name_dim_month,     12,                 region%output%id_dim_month     ) ! Months (for monthly data)
    CALL create_dim( region%output%ncid, region%output%name_dim_time,      nf90_unlimited,     region%output%id_dim_time      ) ! Time frames
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.    

    ! Primary output - everything that's needed to restart a new run, and to plot/analyse data
    ! ========================================================================================
    
    ! Mesh data
    CALL create_double_var( region%output%ncid, region%output%name_var_V,                [region%output%id_dim_vi, region%output%id_dim_two    ], region%output%id_var_V,                long_name='Vertex coordinates', units='m' )
    CALL create_double_var( region%output%ncid, region%output%name_var_Tri,              [region%output%id_dim_ti, region%output%id_dim_three  ], region%output%id_var_Tri,              long_name='Vertex indices'                  )
    CALL create_double_var( region%output%ncid, region%output%name_var_nC,               [region%output%id_dim_vi                              ], region%output%id_var_nC,               long_name='Number of connected vertices'    )
    CALL create_double_var( region%output%ncid, region%output%name_var_C,                [region%output%id_dim_vi, region%output%id_dim_ci     ], region%output%id_var_C,                long_name='Indices of connected vertices'   )
    CALL create_double_var( region%output%ncid, region%output%name_var_niTri,            [region%output%id_dim_vi                              ], region%output%id_var_niTri,            long_name='Number of inverse triangles'     )
    CALL create_double_var( region%output%ncid, region%output%name_var_iTri,             [region%output%id_dim_vi, region%output%id_dim_ci     ], region%output%id_var_iTri,             long_name='Indices of inverse triangles'    )
    CALL create_double_var( region%output%ncid, region%output%name_var_edge_index,       [region%output%id_dim_vi                              ], region%output%id_var_edge_index,       long_name='Edge index'                      )
    CALL create_double_var( region%output%ncid, region%output%name_var_Tricc,            [region%output%id_dim_ti, region%output%id_dim_two    ], region%output%id_var_Tricc,            long_name='Triangle circumcenter', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_var_TriC,             [region%output%id_dim_ti, region%output%id_dim_three  ], region%output%id_var_TriC,             long_name='Triangle neighbours'             )
    CALL create_double_var( region%output%ncid, region%output%name_var_Tri_edge_index,   [region%output%id_dim_ti                              ], region%output%id_var_Tri_edge_index,   long_name='Triangle edge index'             )
    CALL create_double_var( region%output%ncid, region%output%name_var_R,                [region%output%id_dim_vi                              ], region%output%id_var_R,                long_name='Resolution',            units='m')
    CALL create_double_var( region%output%ncid, region%output%name_var_lat,              [region%output%id_dim_vi                              ], region%output%id_var_lat,              long_name='Latitude',  units='degrees north')
    CALL create_double_var( region%output%ncid, region%output%name_var_lon,              [region%output%id_dim_vi                              ], region%output%id_var_lon,              long_name='Longitude', units='degrees east' )
    
    ! Dimension variables: zeta, month, time
    CALL create_double_var( region%output%ncid, region%output%name_var_time,             [region%output%id_dim_time                            ], region%output%id_var_time,             long_name='Time', units='years'   )
    CALL create_double_var( region%output%ncid, region%output%name_var_zeta,             [region%output%id_dim_zeta                            ], region%output%id_var_zeta,             long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( region%output%ncid, region%output%name_var_month,            [region%output%id_dim_month                           ], region%output%id_var_month,            long_name='Month', units='1-12'    )
    
    ! Ice model data
    CALL create_double_var( region%output%ncid, region%output%name_var_Hi,               [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_Hi,               long_name='Ice Thickness', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_var_Hb,               [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_Hb,               long_name='Bedrock Height', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_var_Hs,               [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_Hs,               long_name='Surface Height', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_var_mask,             [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_mask,             long_name='mask', units='unitless')
    CALL create_double_var( region%output%ncid, region%output%name_var_U_SIA,            [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_U_SIA,            long_name='SIA ice x-velocity', units='m/yr')
    CALL create_double_var( region%output%ncid, region%output%name_var_V_SIA,            [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_V_SIA,            long_name='SIA ice y-velocity', units='m/yr')    
    CALL create_double_var( region%output%ncid, region%output%name_var_U_SSA,            [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_U_SSA,            long_name='SSA ice x-velocity', units='m/yr')
    CALL create_double_var( region%output%ncid, region%output%name_var_V_SSA,            [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_V_SSA,            long_name='SSA ice y-velocity', units='m/yr')              
    CALL create_double_var( region%output%ncid, region%output%name_var_Ti,               [region%output%id_dim_vi, region%output%id_dim_zeta,  region%output%id_dim_time], region%output%id_var_Ti,               long_name='Ice temperature', units='K')
    CALL create_double_var( region%output%ncid, region%output%name_var_T2m,              [region%output%id_dim_vi, region%output%id_dim_month, region%output%id_dim_time], region%output%id_var_T2m,              long_name='Monthly mean 2-, air temperature', units='K')
    CALL create_double_var( region%output%ncid, region%output%name_var_Precip,           [region%output%id_dim_vi, region%output%id_dim_month, region%output%id_dim_time], region%output%id_var_Precip,           long_name='Monthly total precipitation', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_var_Albedo,           [region%output%id_dim_vi, region%output%id_dim_month, region%output%id_dim_time], region%output%id_var_Albedo,           long_name='Surface albedo', units='unitless')
    CALL create_double_var( region%output%ncid, region%output%name_var_SMB,              [region%output%id_dim_vi, region%output%id_dim_month, region%output%id_dim_time], region%output%id_var_SMB,              long_name='Surface mass balance', units='mie/yr')
    CALL create_double_var( region%output%ncid, region%output%name_var_BMB,              [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_BMB,              long_name='Basal mass balance', units='mie/yr')
    CALL create_double_var( region%output%ncid, region%output%name_var_FirnDepth,        [region%output%id_dim_vi, region%output%id_dim_month, region%output%id_dim_time], region%output%id_var_FirnDepth,        long_name='Firn depth', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_var_MeltPreviousYear, [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_MeltPreviousYear, long_name='Melt during previous year', units='mie')
    CALL create_double_var( region%output%ncid, region%output%name_var_dHs_dx,           [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_dHs_dx,           long_name='Surface slope in x direction', units='m/m')
    CALL create_double_var( region%output%ncid, region%output%name_var_dHs_dy,           [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_dHs_dy,           long_name='Surface slope in y direction', units='m/m') 
    CALL create_double_var( region%output%ncid, region%output%name_var_D_SIA,            [region%output%id_dim_vi,                             region%output%id_dim_time], region%output%id_var_D_SIA,            long_name='SIA ice diffusivity', units='m/yr')             
                           
    ! Secondary (optional) output data
    ! ================================
    
    ! Ice dynamics
    IF (C%save_secondary_data_ice_dynamics) THEN
      CALL create_double_var( region%output%ncid, region%output%name_var_U_3D,        [region%output%id_dim_vi,  region%output%id_dim_zeta,  region%output%id_dim_time], region%output%id_var_U_3D,      long_name='3D SIA ice x-velocity', units='m/yr')
      CALL create_double_var( region%output%ncid, region%output%name_var_V_3D,        [region%output%id_dim_vi,  region%output%id_dim_zeta,  region%output%id_dim_time], region%output%id_var_V_3D,      long_name='3D SIA ice y-velocity', units='m/yr')
      CALL create_double_var( region%output%ncid, region%output%name_var_W_3D,        [region%output%id_dim_vi,  region%output%id_dim_zeta,  region%output%id_dim_time], region%output%id_var_W_3D,      long_name='3D SIA ice z-velocity', units='m/yr')                             
    END IF
                           
    ! Climate and SMB components
    IF (C%save_secondary_data_SMB_components) THEN
      CALL create_double_var( region%output%ncid, region%output%name_var_Wind_WE,     [region%output%id_dim_vi,  region%output%id_dim_month, region%output%id_dim_time], region%output%id_var_Wind_WE,    long_name='Monthly 10 m west-east wind component', units='m/s')
      CALL create_double_var( region%output%ncid, region%output%name_var_Wind_SN,     [region%output%id_dim_vi,  region%output%id_dim_month, region%output%id_dim_time], region%output%id_var_Wind_SN,    long_name='Monthly 10 m south-north wind component', units='m/s')
      CALL create_double_var( region%output%ncid, region%output%name_var_Snowfall,    [region%output%id_dim_vi,  region%output%id_dim_month, region%output%id_dim_time], region%output%id_var_Snowfall,   long_name='Monthly total snowfall', units='m')
      CALL create_double_var( region%output%ncid, region%output%name_var_Melt,        [region%output%id_dim_vi,  region%output%id_dim_month, region%output%id_dim_time], region%output%id_var_Melt,       long_name='Monthly total melt', units='m')
      CALL create_double_var( region%output%ncid, region%output%name_var_Refreezing,  [region%output%id_dim_vi,                              region%output%id_dim_time], region%output%id_var_Refreezing, long_name='Yearly total refreezing', units='m')
      CALL create_double_var( region%output%ncid, region%output%name_var_Runoff,      [region%output%id_dim_vi,  region%output%id_dim_month, region%output%id_dim_time], region%output%id_var_Runoff,     long_name='Monthly total runoff', units='m')
    END IF

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( region%output%ncid))    
    
    ! Write general data and key mesh data:
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_V,               region%mesh%V             ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Tri,             region%mesh%Tri           ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_nC,              region%mesh%nC            ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_C,               region%mesh%C             ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_niTri,           region%mesh%niTri         ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_iTri,            region%mesh%iTri          ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_edge_index,      region%mesh%edge_index    ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Tricc,           region%mesh%Tricc         ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_TriC,            region%mesh%TriC          ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_Tri_edge_index,  region%mesh%Tri_edge_index))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_R,               region%mesh%R             ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_lat,             region%mesh%lat           ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_lon,             region%mesh%lon           ))
    
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( region%output%ncid, region%output%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( region%output%ncid))
    
    ! Close the file
    CALL close_netcdf_file(region%output%ncid)
    
  END SUBROUTINE create_output_file
  
  SUBROUTINE read_PD_data_file( PD)
    ! Read the PD netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_PD_data_Fields), INTENT(INOUT) :: PD
    
    ! Open the netcdf file
    WRITE(0,*) '   Reading PD    data from file ', TRIM(PD%netcdf%filename), '...'
    CALL open_netcdf_file(PD%netcdf%filename, PD%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_x,      PD%x,           start=(/1   /)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_y,      PD%y,           start=(/1   /)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_Hi,     PD%Hi_cart,     start=(/1, 1/)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_Hb,     PD%Hb_cart,     start=(/1, 1/)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_Hs,     PD%Hs_cart,     start=(/1, 1/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD%netcdf%ncid)
    
  END SUBROUTINE read_PD_data_file 
  SUBROUTINE read_init_data_file_cart( init)
    ! Read the init netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_Fields), INTENT(INOUT) :: init
    
    ! Open the netcdf file
    WRITE(0,*) '   Reading init  data from file ', TRIM(init%netcdf%filename), '...'
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_x,      init%x,           start=(/1   /)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_y,      init%y,           start=(/1   /)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hi,     init%Hi_cart,     start=(/1, 1/)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hb,     init%Hb_cart,     start=(/1, 1/)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hs,     init%Hs_cart,     start=(/1, 1/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE read_init_data_file_cart
  SUBROUTINE read_init_data_file_mesh( init, mesh, netcdf, nt)
    ! Read the init netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    TYPE(type_mesh),                INTENT(INOUT) :: mesh
    TYPE(type_netcdf_output),       INTENT(INOUT) :: netcdf
    INTEGER,                        INTENT(IN)    :: nt
    
    ! Local variables:
    REAL(dp)                                      :: t_file
    
    ! Open the netcdf file
    WRITE(0,*) '   Reading init  data from file ', TRIM(init%netcdf%filename), '...'
    CALL open_netcdf_file(init%netcdf%filename, netcdf%ncid)
    
    ! Throw a warning when the current simulation doesn't start where the one we're initialising from ended    
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_time, t_file, start=(/nt/)))
    IF (t_file /= C%start_time_of_run) THEN
      WRITE(0,'(A,F8.2,A,F8.2,A)') '     WARNING: the init file ends at t = ', t_file, ', but the current simulation starts at t = ', C%start_time_of_run, '!'
    END IF
    
    ! Read mesh data
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_V,              mesh%V,                start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Tri,            mesh%Tri,              start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_nC,             mesh%nC,               start=(/1   /)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_C,              mesh%C,                start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_niTri,          mesh%niTri,            start=(/1   /)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_iTri,           mesh%iTri,             start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_edge_index,     mesh%edge_index,       start=(/1   /)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_TriC,           mesh%TriC,             start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Tricc,          mesh%Tricc,            start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Tri_edge_index, mesh%Tri_edge_index,   start=(/1   /)))
    
    mesh%xmin = mesh%V(1,1)
    mesh%xmax = mesh%V(2,1)
    mesh%ymin = mesh%V(2,2)
    mesh%ymax = mesh%V(3,2)
    
    ! Read init data
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Hi,                  init%Hi,               start=(/1,    nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Hb,                  init%Hb,               start=(/1,    nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Hs,                  init%Hs,               start=(/1,    nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_U_SSA,               init%U_SSA,            start=(/1,    nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_V_SSA,               init%V_SSA,            start=(/1,    nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Ti,                  init%Ti,               start=(/1, 1, nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_Albedo,              init%Albedo,           start=(/1, 1, nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_FirnDepth,           init%FirnDepth,        start=(/1, 1, nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%id_var_MeltPreviousYear,    init%MeltPreviousYear, start=(/1,    nt/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(netcdf%ncid)
    
  END SUBROUTINE read_init_data_file_mesh
  SUBROUTINE read_PD_obs_data_file( PD_obs)
    ! Read the PD_obs0 netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: PD_obs
    
    ! Open the netcdf file
    WRITE(0,*) '   Reading PD observed climate data from file ', TRIM(PD_obs%netcdf%filename), '...'
    CALL open_netcdf_file(PD_obs%netcdf%filename, PD_obs%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_lon,     PD_obs%lon,     start=(/1      /)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_lat,     PD_obs%lat,     start=(/1      /)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_T2m,     PD_obs%T2m,     start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Precip,  PD_obs%Precip,  start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Hs,      PD_obs%Hs,      start=(/1, 1   /)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Albedo,  PD_obs%Albedo,  start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Wind_WE, PD_obs%Wind_WE, start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Wind_SN, PD_obs%Wind_SN, start=(/1, 1, 1/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD_obs%netcdf%ncid)
    
  END SUBROUTINE read_PD_obs_data_file
  SUBROUTINE read_insolation_data_file( forcing, ti0, ti1) 
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    INTEGER,                        INTENT(IN)    :: ti0, ti1
    
    ! Local variables:
    INTEGER                                       :: mi, li
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Q_temp0, Q_temp1
    
    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE( Q_temp0(1, 12, forcing%ins_nlat))
    ALLOCATE( Q_temp1(1, 12, forcing%ins_nlat))
        
    ! Read data
    CALL open_netcdf_file(forcing%netcdf%filename, forcing%netcdf%ncid)
    CALL handle_error(nf90_get_var( forcing%netcdf%ncid, forcing%netcdf%id_var_Q_TOA, Q_temp0, start=(/ti0, 1, 1/), count=(/ 1, 12, forcing%ins_nlat/), stride=(/1, 1, 1/) ))    
    CALL handle_error(nf90_get_var( forcing%netcdf%ncid, forcing%netcdf%id_var_Q_TOA, Q_temp1, start=(/ti1, 1, 1/), count=(/ 1, 12, forcing%ins_nlat/), stride=(/1, 1, 1/) ))
    CALL close_netcdf_file(forcing%netcdf%ncid) 
    
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
   
  END SUBROUTINE read_insolation_data_file
  SUBROUTINE read_insolation_data_file_time_lat( forcing) 
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
    
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf%filename, forcing%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( forcing%netcdf%ncid, forcing%netcdf%id_var_time,    forcing%ins_time,    start=(/1      /)))
    CALL handle_error(nf90_get_var( forcing%netcdf%ncid, forcing%netcdf%id_var_lat,     forcing%ins_lat,     start=(/1      /)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf%ncid)
    
  END SUBROUTINE read_insolation_data_file_time_lat

  SUBROUTINE inquire_PD_data_file( PD) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_PD_data_fields), INTENT(INOUT) :: PD
        
    ! Open the netcdf file
    CALL open_netcdf_file(PD%netcdf%filename, PD%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( PD%netcdf%ncid, PD%netcdf%name_dim_x, PD%nx, PD%netcdf%id_dim_x)
    CALL inquire_dim( PD%netcdf%ncid, PD%netcdf%name_dim_y, PD%ny, PD%netcdf%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_var_x,  (/ PD%netcdf%id_dim_x                     /), PD%netcdf%id_var_x)
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_var_y,  (/ PD%netcdf%id_dim_y                     /), PD%netcdf%id_var_y)
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_var_Hi, (/ PD%netcdf%id_dim_x, PD%netcdf%id_dim_y /), PD%netcdf%id_var_Hi)
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_var_Hb, (/ PD%netcdf%id_dim_x, PD%netcdf%id_dim_y /), PD%netcdf%id_var_Hb)
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_var_Hs, (/ PD%netcdf%id_dim_x, PD%netcdf%id_dim_y /), PD%netcdf%id_var_Hs)
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD%netcdf%ncid)
    
  END SUBROUTINE inquire_PD_data_file
  SUBROUTINE inquire_init_data_file_cart( init) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields), INTENT(INOUT) :: init
        
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_x, init%nx, init%netcdf%id_dim_x)
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_y, init%ny, init%netcdf%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_x,  (/ init%netcdf%id_dim_x                       /), init%netcdf%id_var_x)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_y,  (/ init%netcdf%id_dim_y                       /), init%netcdf%id_var_y)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hi, (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y /), init%netcdf%id_var_Hi)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hb, (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y /), init%netcdf%id_var_Hb)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hs, (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y /), init%netcdf%id_var_Hs)
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE inquire_init_data_file_cart
  SUBROUTINE inquire_init_data_file_mesh( init, mesh, netcdf, nt)
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    TYPE(type_mesh),                INTENT(INOUT) :: mesh
    TYPE(type_netcdf_output),       INTENT(INOUT) :: netcdf
    INTEGER,                        INTENT(OUT)   :: nt
 
    ! Local variables:
    INTEGER                               :: int_dummy
        
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_vi,      mesh%nV,        netcdf%id_dim_vi    )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_two,     int_dummy,      netcdf%id_dim_two   )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_ti,      mesh%nTri,      netcdf%id_dim_ti    )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_three,   int_dummy,      netcdf%id_dim_three )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_ci,      mesh%nC_mem,    netcdf%id_dim_ci    )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_zeta,    int_dummy,      netcdf%id_dim_zeta  )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_month,   int_dummy,      netcdf%id_dim_month )
    CALL inquire_dim( netcdf%ncid, netcdf%name_dim_time,    nt,             netcdf%id_dim_time  )
    
    mesh%nV_mem   = mesh%nV
    mesh%nTri_mem = mesh%nTri

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_zeta,                (/ netcdf%id_dim_zeta                                          /), netcdf%id_var_zeta)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_month,               (/ netcdf%id_dim_month                                         /), netcdf%id_var_month)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_time,                (/ netcdf%id_dim_time                                          /), netcdf%id_var_time)
    
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_V,                   (/ netcdf%id_dim_vi,   netcdf%id_dim_two                       /), netcdf%id_var_V)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_nC,                  (/ netcdf%id_dim_vi                                            /), netcdf%id_var_nC)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_C,                   (/ netcdf%id_dim_vi,   netcdf%id_dim_ci                        /), netcdf%id_var_C)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_niTri,               (/ netcdf%id_dim_vi                                            /), netcdf%id_var_niTri)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_iTri,                (/ netcdf%id_dim_vi,   netcdf%id_dim_ci                        /), netcdf%id_var_iTri)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_edge_index,          (/ netcdf%id_dim_vi                                            /), netcdf%id_var_edge_index)
    
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Tri,                 (/ netcdf%id_dim_ti,   netcdf%id_dim_three                     /), netcdf%id_var_Tri)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Tricc,               (/ netcdf%id_dim_ti,   netcdf%id_dim_two                       /), netcdf%id_var_Tricc)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_TriC,                (/ netcdf%id_dim_ti,   netcdf%id_dim_three                     /), netcdf%id_var_TriC)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Tri_edge_index,      (/ netcdf%id_dim_ti                                            /), netcdf%id_var_Tri_edge_index)
        
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Hi,                  (/ netcdf%id_dim_vi,   netcdf%id_dim_time                      /), netcdf%id_var_Hi)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Hb,                  (/ netcdf%id_dim_vi,   netcdf%id_dim_time                      /), netcdf%id_var_Hb)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Hs,                  (/ netcdf%id_dim_vi,   netcdf%id_dim_time                      /), netcdf%id_var_Hs)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_U_SSA,               (/ netcdf%id_dim_vi,   netcdf%id_dim_time                      /), netcdf%id_var_U_SSA)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_V_SSA,               (/ netcdf%id_dim_vi,   netcdf%id_dim_time                      /), netcdf%id_var_V_SSA)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Ti,                  (/ netcdf%id_dim_vi,   netcdf%id_dim_zeta,  netcdf%id_dim_time /), netcdf%id_var_Ti)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_Albedo,              (/ netcdf%id_dim_vi,   netcdf%id_dim_month, netcdf%id_dim_time /), netcdf%id_var_Albedo)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_FirnDepth,           (/ netcdf%id_dim_vi,   netcdf%id_dim_month, netcdf%id_dim_time /), netcdf%id_var_FirnDepth)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_var_MeltPreviousYear,    (/ netcdf%id_dim_vi,   netcdf%id_dim_time                      /), netcdf%id_var_MeltPreviousYear)
        
    ! Close the netcdf file
    CALL close_netcdf_file(netcdf%ncid)
    
  END SUBROUTINE inquire_init_data_file_mesh
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
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_dim_lat,     PD_obs%nlat,  PD_obs%netcdf%id_dim_lat)
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_dim_lon,     PD_obs%nlon,  PD_obs%netcdf%id_dim_lon)
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_dim_month,   int_dummy,    PD_obs%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_lat,      (/ PD_obs%netcdf%id_dim_lat                                                       /),  PD_obs%netcdf%id_var_lat)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_lon,      (/ PD_obs%netcdf%id_dim_lon                                                       /),  PD_obs%netcdf%id_var_lon)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_T2m,      (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_T2m)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Precip,   (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_Precip)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Hs,       (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat                             /),  PD_obs%netcdf%id_var_Hs)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Albedo,   (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_Albedo)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Wind_WE,  (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_Wind_WE)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Wind_SN,  (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_Wind_SN)
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD_obs%netcdf%ncid)
    
  END SUBROUTINE inquire_PD_obs_data_file 
  SUBROUTINE inquire_insolation_data_file( forcing)
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
 
    ! Local variables:  
    INTEGER                                :: int_dummy
            
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf%filename, forcing%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( forcing%netcdf%ncid, forcing%netcdf%name_dim_time,     forcing%ins_nyears,        forcing%netcdf%id_dim_time)
    CALL inquire_dim( forcing%netcdf%ncid, forcing%netcdf%name_dim_month,    int_dummy,                 forcing%netcdf%id_dim_month)  
    CALL inquire_dim( forcing%netcdf%ncid, forcing%netcdf%name_dim_lat,      forcing%ins_nlat,          forcing%netcdf%id_dim_lat)
    
    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( forcing%netcdf%ncid, forcing%netcdf%name_var_time,  (/ forcing%netcdf%id_dim_time                                                         /), forcing%netcdf%id_var_time)
    CALL inquire_double_var( forcing%netcdf%ncid, forcing%netcdf%name_var_month, (/ forcing%netcdf%id_dim_month                                                        /), forcing%netcdf%id_var_month)
    CALL inquire_double_var( forcing%netcdf%ncid, forcing%netcdf%name_var_lat,   (/ forcing%netcdf%id_dim_lat                                                          /), forcing%netcdf%id_var_lat)
    CALL inquire_double_var( forcing%netcdf%ncid, forcing%netcdf%name_var_Q_TOA, (/ forcing%netcdf%id_dim_time, forcing%netcdf%id_dim_month, forcing%netcdf%id_dim_lat /), forcing%netcdf%id_var_Q_TOA)
        
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf%ncid)
    
  END SUBROUTINE inquire_insolation_data_file

  SUBROUTINE GetOutputFilename( region)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region

    CHARACTER(LEN=20)          :: short_filename
    LOGICAL                    :: ex
    INTEGER                    :: n    
    CHARACTER(LEN=256)         :: ns

    short_filename = 'results_NAM_00001.nc'
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
      region%output%filename(n:n) = ' '
    END DO
    region%output%filename = TRIM(C%output_dir)//TRIM(short_filename)

  END SUBROUTINE GetOutputFilename

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
    
END MODULE netcdf_module
