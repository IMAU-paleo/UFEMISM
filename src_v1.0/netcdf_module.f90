MODULE netcdf_module
  ! This module defines methods to read and write the current state of the model from or to a NetCDF file.

  USE netcdf,                      ONLY: nf90_max_var_dims, nf90_create, nf90_close, nf90_clobber, nf90_share, nf90_unlimited , &
                                         nf90_enddef, nf90_put_var, nf90_sync, nf90_def_var, nf90_int, nf90_put_att, nf90_def_dim, &
                                         nf90_open, nf90_write, nf90_inq_dimid, nf90_inquire_dimension, nf90_inquire, nf90_double, &
                                         nf90_inq_varid, nf90_inquire_variable, nf90_get_var, nf90_noerr, nf90_strerror
  USE parallel_module,             ONLY: par, sync
  USE configuration_module,        ONLY: dp, C
  USE data_types_module,           ONLY: type_model_region, type_mesh, type_ice_model, type_climate_model, type_PD_data_Fields, &
                                         type_init_data_Fields, type_subclimate_global, type_forcing_data, type_netcdf_file_output
  IMPLICIT NONE

CONTAINS

  ! Create and write to an output NetCDF file
  SUBROUTINE write_to_output_file(region)
    ! Write the current model state to the existing output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region
    
    WRITE(0,'(A,F8.2,A)') '   t = ', region%time/1e3, ' kyr - writing output...'
    
    ! Open the file for writing
    CALL open_netcdf_file(region%output%filename, region%output%ncid)
        
    ! Time
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_time,        region%time,                start=(/       region%output%ti/)))
    
    ! Key output
    ! ==========
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Hi,             region%ice%Hi,              start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Hb,             region%ice%Hb,              start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Hs,             region%ice%Hs,              start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Us,             region%ice%Us,              start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Vs,             region%ice%Vs,              start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_U_SSA,          region%ice%U_SSA,           start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_V_SSA,          region%ice%V_SSA,           start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Ti,             region%ice%Ti,              start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_mask,           region%ice%mask,            start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_T2m,            region%climate%applied%T2m,         start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Precip,         region%climate%applied%Precip,      start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Albedo,         region%SMB%Albedo,          start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_SMB,            region%SMB%SMB,             start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_BMB,            region%BMB%BMB,             start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_FirnDepth,      region%SMB%FirnDepth,       start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_MeltPreviousYear, region%SMB%MeltPreviousYear, start=(/1,  region%output%ti/)))
    
    ! Optional output
    ! ===============
    
    ! SMB components
    IF (C%save_SMB_components) THEN
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Wind_WE,        region%climate%applied%Wind_WE,         start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Wind_SN,        region%climate%applied%Wind_SN,         start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Snowfall,       region%SMB%Snowfall,        start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Melt,           region%SMB%Melt,            start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Refreezing,     region%SMB%Refreezing_year, start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Runoff,         region%SMB%Runoff,          start=(/1, 1,  region%output%ti/)))
    END IF
    
    ! Ice dynamics
    IF (C%save_ice_dynamics) THEN
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_dHs_dx,         region%ice%dHs_dx,          start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_dHs_dy,         region%ice%dHs_dy,          start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_U_SIA,          region%ice%U_SIA,           start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_V_SIA,          region%ice%V_SIA,           start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Ux_ac,          region%ice%Ux_ac,           start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_Uy_ac,          region%ice%Uy_ac,           start=(/1,     region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_U_3D,           region%ice%U_3D,            start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_V_3D,           region%ice%V_3D,            start=(/1, 1,  region%output%ti/)))
    CALL handle_error(nf90_put_var(region%output%ncid, region%output%var_W_3D,           region%ice%W_3D,            start=(/1, 1,  region%output%ti/)))
    END IF
    
    ! Close the file
    CALL close_netcdf_file(region%output%ncid)
    
    ! Increase time frame counter
    region%output%ti = region%output%ti + 1
        
  END SUBROUTINE write_to_output_file
  SUBROUTINE create_output_file(region)
    ! Create a new output NetCDF file, write the current mesh data to it.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region

    ! Local variables:
    INTEGER                             :: dim_nV, dim_nV2, dim_nTri, dim_nTri2, dim_AcV, dim_AcV2, dim_nconmax, dim_nconmax2, dim_zeta, dim_t, dim_month
    LOGICAL                             :: file_exists
    
    ! Get output file name
    CALL GetOutputFilename(region)
    
    ! Set time frame index to 1
    region%output%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(region%output%filename))
    IF(file_exists) THEN
     WRITE(0,'(5A)') 'ERROR: ', TRIM(region%output%filename), ' already exists!'
     STOP
    END IF
    
    ! Create netCDF file
    CALL handle_error(nf90_create(region%output%filename,IOR(nf90_clobber,nf90_share),region%output%ncid))
        
    ! Define dimensions:
    CALL create_dim( region%output%ncid, region%output%name_mesh_V,    region%mesh%nV,        dim_nV       ) ! Vertices
    CALL create_dim( region%output%ncid, region%output%name_mesh_V2,   2,                     dim_nV2      ) ! The vertices themselves have both x and y coordinates
    CALL create_dim( region%output%ncid, region%output%name_mesh_Tri,  region%mesh%nTri,      dim_nTri     ) ! Triangles
    CALL create_dim( region%output%ncid, region%output%name_mesh_Tri2, 3,                     dim_nTri2    ) ! Each triangle has three vertices
    CALL create_dim( region%output%ncid, region%output%name_mesh_AcV,  region%mesh%nAC,       dim_AcV      ) ! Arakawa C vertices
    CALL create_dim( region%output%ncid, 'Aa neighbour',               4,                     dim_AcV2     ) ! Each Arakawa C vertex has 4 regular vertex neighbours
    CALL create_dim( region%output%ncid, 'neighbour',                  region%mesh%nconmax,   dim_nconmax  ) ! Vertex connections
    CALL create_dim( region%output%ncid, 'neighbour+1',                region%mesh%nconmax+1, dim_nconmax2 ) ! Vertex connections + vertex itself (for neighbour functions)
    CALL create_dim( region%output%ncid, region%output%name_zeta,      C%NZ,                  dim_zeta     ) ! Scaled vertical coordinate
    CALL create_dim( region%output%ncid, region%output%name_month,     12,                    dim_month    ) ! Months (for monthly data)
    CALL create_dim( region%output%ncid, region%output%name_time,      nf90_unlimited,        dim_t        ) ! Time frames
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    ! Mesh data:
    CALL create_double_var( region%output%ncid, region%output%name_mesh_V,     [dim_nV,   dim_nV2],   region%output%var_mesh_V, &
                            units='m')
    CALL create_double_var( region%output%ncid, region%output%name_mesh_Tri,   [dim_nTri, dim_nTri2], region%output%var_mesh_Tri, &
                            units='vertex number')
    CALL create_double_var( region%output%ncid, region%output%name_mesh_nC,    [dim_nV],              region%output%var_mesh_nC, &
                           long_name='number of neighbours')
    CALL create_double_var( region%output%ncid, region%output%name_mesh_C,     [dim_nV, dim_nconmax], region%output%var_mesh_C, &
                           long_name='neighbours')
    CALL create_double_var( region%output%ncid, region%output%name_mesh_niTri, [dim_nV],              region%output%var_mesh_niTri, &
                           long_name='number of inverse triangles')
    CALL create_double_var( region%output%ncid, region%output%name_mesh_iTri,  [dim_nV, dim_nconmax],  region%output%var_mesh_iTri, &
                           long_name='inverse triangles')
    CALL create_double_var( region%output%ncid, region%output%name_mesh_edge_index, [dim_nV],  region%output%var_mesh_edge_index, &
                           long_name='edge index')
    CALL create_double_var( region%output%ncid, region%output%name_mesh_Tricc, [dim_nTri, dim_nV2],  region%output%var_mesh_Tricc, &
                           long_name='triangle circumcenter', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_mesh_TriC, [dim_nTri, dim_nTri2],  region%output%var_mesh_TriC, &
                           long_name='triangle neighbours')
    CALL create_double_var( region%output%ncid, region%output%name_mesh_Tri_edge_index, [dim_nTri],  region%output%var_mesh_Tri_edge_index, &
                           long_name='triangle edge index')
    CALL create_double_var( region%output%ncid, region%output%name_mesh_R,   [dim_nV], region%output%var_mesh_R, &
                           long_name = 'Resolution', units='m')

    ! zeta, time, month, lat, lon:
    CALL create_double_var( region%output%ncid, region%output%name_zeta,  [dim_zeta],  region%output%var_zeta,  long_name='Vertical scaled coordinate',   units='unitless')
    CALL create_double_var( region%output%ncid, region%output%name_month, [dim_month], region%output%var_month, long_name='Month',                        units='1-12')
    CALL create_double_var( region%output%ncid, region%output%name_lat,   [dim_nV],    region%output%var_lat,   long_name='Latitude',                     units='degrees north')
    CALL create_double_var( region%output%ncid, region%output%name_lon,   [dim_nV],    region%output%var_lon,   long_name='Longitude',                    units='degrees east')
    CALL create_double_var( region%output%ncid, region%output%name_time,  [dim_t],     region%output%var_time,  long_name='Time',                         units='years')

    ! Key output
    CALL create_double_var( region%output%ncid, region%output%name_Hi, [dim_nV, dim_t],  region%output%var_Hi, &
                           long_name='Ice Thickness', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_Hb, [dim_nV, dim_t],  region%output%var_Hb, &
                           long_name='Bedrock Height', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_Hs, [dim_nV, dim_t],  region%output%var_Hs, &
                           long_name='Surface Height', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_mask, [dim_nV, dim_t],  region%output%var_mask, &
                           long_name='mask', units='unitless')
    CALL create_double_var( region%output%ncid, region%output%name_Us, [dim_nV, dim_t],  region%output%var_Us, &
                           long_name='Ice surface velocity in x-direction', units='m/yr')
    CALL create_double_var( region%output%ncid, region%output%name_Vs, [dim_nV, dim_t],  region%output%var_Vs, &
                           long_name='Ice surface velocity in y-direction', units='m/yr')
    CALL create_double_var( region%output%ncid, region%output%name_U_SSA, [dim_nV, dim_t],  region%output%var_U_SSA, &
                           long_name='SSA ice x-velocity', units='m/yr')
    CALL create_double_var( region%output%ncid, region%output%name_V_SSA, [dim_nV, dim_t],  region%output%var_V_SSA, &
                           long_name='SSA ice y-velocity', units='m/yr')              
    CALL create_double_var( region%output%ncid, region%output%name_Ti, [dim_nV, dim_zeta, dim_t],  region%output%var_Ti, &
                           long_name='Ice temperature', units='K')
    CALL create_double_var( region%output%ncid, region%output%name_T2m, [dim_nV, dim_month, dim_t],  region%output%var_T2m, &
                           long_name='Monthly mean 2 meter air temperature', units='K')
    CALL create_double_var( region%output%ncid, region%output%name_Precip, [dim_nV, dim_month, dim_t],  region%output%var_Precip, &
                           long_name='Monthly total precipitation', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_Albedo, [dim_nV, dim_month, dim_t],  region%output%var_Albedo, &
                           long_name='Surface albedo', units='unitless')
    CALL create_double_var( region%output%ncid, region%output%name_SMB, [dim_nV, dim_month, dim_t],  region%output%var_SMB, &
                           long_name='Surface mass balance', units='mie/yr')
    CALL create_double_var( region%output%ncid, region%output%name_BMB, [dim_nV, dim_t],  region%output%var_BMB, &
                           long_name='Basal mass balance', units='mie/yr')
    CALL create_double_var( region%output%ncid, region%output%name_FirnDepth, [dim_nV, dim_month, dim_t],  region%output%var_FirnDepth, &
                           long_name='Firn depth', units='m')
    CALL create_double_var( region%output%ncid, region%output%name_MeltPreviousYear, [dim_nV, dim_t],  region%output%var_MeltPreviousYear, &
                           long_name='Melt during previous year', units='mie')
                           
    ! Optional output data
    ! ====================
    
    ! Mesh data
    IF (C%save_mesh_data) THEN
      CALL create_double_var( region%output%ncid, region%output%name_mesh_A,   [dim_nV], region%output%var_mesh_A, &
                             long_name = 'Voronoi cell area', units='m^2')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_Cw, [dim_nV, dim_nconmax],  region%output%var_mesh_Cw, &
                             long_name='connection width', units='m')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_TriA, [dim_nTri],  region%output%var_mesh_TriA, &
                             long_name='triangle area', units='m^2')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_NxTri, [dim_nTri, dim_nTri2],  region%output%var_mesh_NxTri, &
                             long_name='Triangle neighbour function NxTri')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_NyTri, [dim_nTri, dim_nTri2],  region%output%var_mesh_NyTri, &
                             long_name='Triangle neighbour function NyTri')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_Nx, [dim_nV, dim_nconmax2],  region%output%var_mesh_Nx, &
                             long_name='First order neighbours function Nx')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_Ny, [dim_nV, dim_nconmax2],  region%output%var_mesh_Ny, &
                             long_name='First order neighbours function Ny')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_Nxx, [dim_nV, dim_nconmax2],  region%output%var_mesh_Nxx, &
                             long_name='First order neighbours function Nxx')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_Nxy, [dim_nV, dim_nconmax2],  region%output%var_mesh_Nxy, &
                             long_name='First order neighbours function Nxy')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_Nyy, [dim_nV, dim_nconmax2],  region%output%var_mesh_Nyy, &
                             long_name='First order neighbours function Nyy')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_VAc, [dim_AcV, dim_nV2], region%output%var_mesh_VAc, &
                             long_name='Coordinates of Ac vertices')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_Aci, [dim_AcV, dim_AcV2], region%output%var_mesh_Aci, &
                             long_name='Arakawa A to C mesh map')
      CALL create_double_var( region%output%ncid, region%output%name_mesh_iAci, [dim_nV, dim_nconmax], region%output%var_mesh_iAci, &
                             long_name='Arakawa C to A mesh map')
    END IF
                           
    ! Climate and SMB components
    IF (C%save_SMB_components) THEN
      CALL create_double_var( region%output%ncid, region%output%name_Wind_WE, [dim_nV, dim_month, dim_t],  region%output%var_Wind_WE, &
                             long_name='Monthly 10 m west-east wind component', units='m/s')
      CALL create_double_var( region%output%ncid, region%output%name_Wind_SN, [dim_nV, dim_month, dim_t],  region%output%var_Wind_SN, &
                             long_name='Monthly 10 m south-north wind component', units='m/s')
      CALL create_double_var( region%output%ncid, region%output%name_Snowfall, [dim_nV, dim_month, dim_t],  region%output%var_Snowfall, &
                             long_name='Monthly total snowfall', units='m')
      CALL create_double_var( region%output%ncid, region%output%name_Melt, [dim_nV, dim_month, dim_t],  region%output%var_Melt, &
                             long_name='Monthly total melt', units='m')
      CALL create_double_var( region%output%ncid, region%output%name_Refreezing, [dim_nV, dim_t],  region%output%var_Refreezing, &
                             long_name='Yearly total refreezing', units='m')
      CALL create_double_var( region%output%ncid, region%output%name_Runoff, [dim_nV, dim_month, dim_t],  region%output%var_Runoff, &
                             long_name='Monthly total runoff', units='m')
    END IF
    
    ! Ice dynamics
    IF (C%save_ice_dynamics) THEN
      CALL create_double_var( region%output%ncid, region%output%name_dHs_dx, [dim_nV, dim_t],  region%output%var_dHs_dx, &
                             long_name='Surface slope in x direction', units='m/m')
      CALL create_double_var( region%output%ncid, region%output%name_dHs_dy, [dim_nV, dim_t],  region%output%var_dHs_dy, &
                             long_name='Surface slope in y direction', units='m/m')
      CALL create_double_var( region%output%ncid, region%output%name_U_SIA, [dim_nV, dim_t],  region%output%var_U_SIA, &
                             long_name='SIA ice x-velocity', units='m/yr')
      CALL create_double_var( region%output%ncid, region%output%name_V_SIA, [dim_nV, dim_t],  region%output%var_V_SIA, &
                             long_name='SIA ice y-velocity', units='m/yr')                        
      CALL create_double_var( region%output%ncid, region%output%name_Ux_ac, [dim_AcV, dim_t],  region%output%var_Ux_ac, &
                             long_name='Ux_ac', units='m')             
      CALL create_double_var( region%output%ncid, region%output%name_Uy_ac, [dim_AcV, dim_t],  region%output%var_Uy_ac, &
                             long_name='Uy_ac', units='m') 
      CALL create_double_var( region%output%ncid, region%output%name_U_3D, [dim_nV, dim_zeta, dim_t],  region%output%var_U_3D, &
                             long_name='U_3D', units='m/yr')
      CALL create_double_var( region%output%ncid, region%output%name_V_3D, [dim_nV, dim_zeta, dim_t],  region%output%var_V_3D, &
                             long_name='V_3D', units='m/yr')
      CALL create_double_var( region%output%ncid, region%output%name_W_3D, [dim_nV, dim_zeta, dim_t],  region%output%var_W_3D, &
                             long_name='W_3D', units='m/yr')
                             
    END IF

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( region%output%ncid))    
    
    ! Write general data and key mesh data:
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_V,               region%mesh%V             ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_Tri,             region%mesh%Tri           ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_nC,              region%mesh%nC            ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_C,               region%mesh%C             ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_niTri,           region%mesh%niTri         ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_iTri,            region%mesh%iTri          ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_edge_index,      region%mesh%edge_index    ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_Tricc,           region%mesh%Tricc         ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_TriC,            region%mesh%TriC          ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_Tri_edge_index,  region%mesh%Tri_edge_index))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_R,               region%mesh%R             ))
    
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_zeta,     C%zeta                                   ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_lat,      region%mesh%lat                          ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_lon,      region%mesh%lon                          ))
  
    ! Optional: write extra mesh data
    IF (C%save_mesh_data) THEN
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_A,               region%mesh%A             ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_Cw,              region%mesh%Cw            ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_TriA,            region%mesh%TriA          ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_NxTri,           region%mesh%NxTri         ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_NyTri,           region%mesh%NyTri         ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_Nx,              region%mesh%Nx            ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_Ny,              region%mesh%Ny            ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_Nxx,             region%mesh%Nxx           ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_Nxy,             region%mesh%Nxy           ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_Nyy,             region%mesh%Nyy           ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_VAc,             region%mesh%VAc           ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_Aci,             region%mesh%Aci           ))
    CALL handle_error(nf90_put_var( region%output%ncid, region%output%var_mesh_iAci,            region%mesh%iAci          ))
    END IF
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( region%output%ncid))
    
    ! Close the file
    CALL close_netcdf_file(region%output%ncid)
    
  END SUBROUTINE create_output_file
  
  SUBROUTINE read_PD_data_file(PD)
    ! Read the PD netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_PD_data_Fields), INTENT(INOUT) :: PD
    
    ! Open the netcdf file
    WRITE(0,*) '   Reading PD    data from file ', TRIM(PD%netcdf%filename), '...'
    CALL open_netcdf_file(PD%netcdf%filename, PD%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%var_x,      PD%x,           start=(/1   /)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%var_y,      PD%y,           start=(/1   /)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%var_Hi,     PD%Hi_cart,     start=(/1, 1/)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%var_Hb,     PD%Hb_cart,     start=(/1, 1/)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%var_Hs,     PD%Hs_cart,     start=(/1, 1/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD%netcdf%ncid)
    
  END SUBROUTINE read_PD_data_file 
  SUBROUTINE read_init_data_file_cart(init)
    ! Read the init netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_Fields), INTENT(INOUT) :: init
    
    ! Open the netcdf file
    WRITE(0,*) '   Reading init  data from file ', TRIM(init%netcdf%filename), '...'
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%var_x,      init%x,           start=(/1   /)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%var_y,      init%y,           start=(/1   /)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%var_Hi,     init%Hi_cart,     start=(/1, 1/)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%var_Hb,     init%Hb_cart,     start=(/1, 1/)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%var_Hs,     init%Hs_cart,     start=(/1, 1/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE read_init_data_file_cart
  SUBROUTINE read_init_data_file_mesh(init, mesh, netcdf, nt)
    ! Read the init netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    TYPE(type_mesh),                INTENT(INOUT) :: mesh
    TYPE(type_netcdf_file_output),  INTENT(INOUT) :: netcdf
    INTEGER,                        INTENT(IN)    :: nt
    
    ! Local variables:
    REAL(dp)                                      :: t_file
    
    ! Open the netcdf file
    WRITE(0,*) '   Reading init  data from file ', TRIM(init%netcdf%filename), '...'
    CALL open_netcdf_file(init%netcdf%filename, netcdf%ncid)
    
    ! Throw a warning when the current simulation doesn't start where the one we're initialising from ended    
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_time, t_file, start=(/nt/)))
    IF (t_file /= C%start_time_of_run) THEN
      WRITE(0,'(A,F8.2,A,F8.2,A)') '     WARNING: the init file ends at t = ', t_file, ', but the current simulation starts at t = ', C%start_time_of_run, '!'
    END IF
    
    ! Read mesh data
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_mesh_V,              mesh%V,                start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_mesh_Tri,            mesh%Tri,              start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_mesh_nC,             mesh%nC,               start=(/1   /)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_mesh_C,              mesh%C,                start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_mesh_niTri,          mesh%niTri,            start=(/1   /)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_mesh_iTri,           mesh%iTri,             start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_mesh_edge_index,     mesh%edge_index,       start=(/1   /)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_mesh_TriC,           mesh%TriC,             start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_mesh_Tricc,          mesh%Tricc,            start=(/1, 1/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_mesh_Tri_edge_index, mesh%Tri_edge_index,   start=(/1   /)))
    
    mesh%xmin = mesh%V(1,1)
    mesh%xmax = mesh%V(2,1)
    mesh%ymin = mesh%V(2,2)
    mesh%ymax = mesh%V(3,2)
    
    ! Read init data
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_Hi,                  init%Hi,               start=(/1,    nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_Hb,                  init%Hb,               start=(/1,    nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_Hs,                  init%Hs,               start=(/1,    nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_U_SSA,               init%U_SSA,            start=(/1,    nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_V_SSA,               init%V_SSA,            start=(/1,    nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_Ti,                  init%Ti,               start=(/1, 1, nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_Albedo,              init%Albedo,           start=(/1, 1, nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_FirnDepth,           init%FirnDepth,        start=(/1, 1, nt/)))
    CALL handle_error(nf90_get_var( netcdf%ncid, netcdf%var_MeltPreviousYear,    init%MeltPreviousYear, start=(/1,    nt/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(netcdf%ncid)
    
  END SUBROUTINE read_init_data_file_mesh
  SUBROUTINE read_PD_obs_data_file(PD_obs)
    ! Read the PD_obs0 netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: PD_obs
    
    ! Open the netcdf file
    WRITE(0,*) '   Reading PD observed climate data from file ', TRIM(PD_obs%netcdf%filename), '...'
    CALL open_netcdf_file(PD_obs%netcdf%filename, PD_obs%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%var_lon,     PD_obs%lon,     start=(/1      /)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%var_lat,     PD_obs%lat,     start=(/1      /)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%var_T2m,     PD_obs%T2m,     start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%var_Precip,  PD_obs%Precip,  start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%var_Hs,      PD_obs%Hs,      start=(/1, 1   /)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%var_Albedo,  PD_obs%Albedo,  start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%var_Wind_WE, PD_obs%Wind_WE, start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%var_Wind_SN, PD_obs%Wind_SN, start=(/1, 1, 1/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD_obs%netcdf%ncid)
    
  END SUBROUTINE read_PD_obs_data_file
  SUBROUTINE read_insolation_data_file(forcing, time) 
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
    REAL(dp),                INTENT(IN)    :: time
    
    ! Local variables
    INTEGER                                :: ti0, ti1, mi, li
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Q_temp0, Q_temp1
    
    IF (par%master) WRITE(0,*) ' Updating insolation forcing...'
    
    ! Initialise at zero
    forcing%ins_Q_TOA0 = 0._dp
    forcing%ins_Q_TOA1 = 0._dp
    
    ! Check if data for model time is available
    IF (time < forcing%ins_time(1)) THEN
      IF (par%master) WRITE(0,*) '  ERROR: insolation data only available between ', MINVAL(forcing%ins_time), ' y and ', MAXVAL(forcing%ins_time), ' y'
      STOP
    END IF
    
    ! Find time indices to be read
    IF (time <= forcing%ins_time(forcing%ins_nyears)) THEN
      IF (par%master) WRITE(0,*) ''
      ti1 = 1
      DO WHILE (forcing%ins_time(ti1) < time)
        ti1 = ti1 + 1
      END DO
      ti0 = ti1 - 1
      
      forcing%ins_t0 = forcing%ins_time(ti0)
      forcing%ins_t1 = forcing%ins_time(ti1)
    ELSE
      IF (par%master) WRITE(0,*) '  WARNING: using constant PD insolation for future projections!'
      IF (par%master) WRITE(0,*) ''
      ti0 = forcing%ins_nyears
      ti1 = forcing%ins_nyears
      
      forcing%ins_t0 = forcing%ins_time(ti0) - 1._dp
      forcing%ins_t1 = forcing%ins_time(ti1)
    END IF
    
    
    IF (par%master) THEN    
      ! Read data
      CALL open_netcdf_file(forcing%ins_netcdf%filename, forcing%ins_netcdf%ncid)
      
      ALLOCATE(Q_temp0(1, 12, forcing%ins_nlat))
      ALLOCATE(Q_temp1(1, 12, forcing%ins_nlat))
      
      CALL handle_error(nf90_get_var( forcing%ins_netcdf%ncid, forcing%ins_netcdf%var_Q_TOA, Q_temp0, &
                                      start=(/ti0, 1, 1/), count=(/ 1, 12, forcing%ins_nlat/), stride=(/1, 1, 1/) ))    
      CALL handle_error(nf90_get_var( forcing%ins_netcdf%ncid, forcing%ins_netcdf%var_Q_TOA, Q_temp1, &
                                      start=(/ti1, 1, 1/), count=(/ 1, 12, forcing%ins_nlat/), stride=(/1, 1, 1/) ))
      
      DO mi = 1, 12
      DO li = 1, forcing%ins_nlat
        forcing%ins_Q_TOA0(li,mi) = Q_temp0(1,mi,li)
        forcing%ins_Q_TOA1(li,mi) = Q_temp1(1,mi,li)
      END DO
      END DO
          
      DEALLOCATE(Q_temp0)
      DEALLOCATE(Q_temp1)
          
      ! Close the netcdf file
      CALL close_netcdf_file(forcing%ins_netcdf%ncid) 
    END IF 
    CALL sync
   
  END SUBROUTINE read_insolation_data_file
  SUBROUTINE read_insolation_data_file_time_lat(forcing) 
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
    
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%ins_netcdf%filename, forcing%ins_netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( forcing%ins_netcdf%ncid, forcing%ins_netcdf%var_time,    forcing%ins_time,    start=(/1      /)))
    CALL handle_error(nf90_get_var( forcing%ins_netcdf%ncid, forcing%ins_netcdf%var_lat,     forcing%ins_lat,     start=(/1      /)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%ins_netcdf%ncid)
    
  END SUBROUTINE read_insolation_data_file_time_lat

  SUBROUTINE inquire_PD_data_file(PD) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_PD_data_fields), INTENT(INOUT) :: PD
 
    ! Local variables:  
    INTEGER                               :: dim_x
    INTEGER                               :: dim_y
        
    ! Open the netcdf file
    CALL open_netcdf_file(PD%netcdf%filename, PD%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( PD%netcdf%ncid, PD%netcdf%name_x,      PD%nx,             dim_x)
    CALL inquire_dim( PD%netcdf%ncid, PD%netcdf%name_y,      PD%ny,             dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_x,      (/ dim_x /),         PD%netcdf%var_x)
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_y,      (/ dim_y /),         PD%netcdf%var_y)
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_Hi,     (/ dim_x, dim_y /),  PD%netcdf%var_Hi)
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_Hb,     (/ dim_x, dim_y /),  PD%netcdf%var_Hb)
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_Hs,     (/ dim_x, dim_y /),  PD%netcdf%var_Hs)
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD%netcdf%ncid)
    
  END SUBROUTINE inquire_PD_data_file
  SUBROUTINE inquire_init_data_file_cart(init) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields), INTENT(INOUT) :: init
 
    ! Local variables:  
    INTEGER                               :: dim_x
    INTEGER                               :: dim_y
        
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_x,      init%nx,             dim_x)
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_y,      init%ny,             dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_x,      (/ dim_x /),         init%netcdf%var_x)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_y,      (/ dim_y /),         init%netcdf%var_y)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_Hi,     (/ dim_x, dim_y /),  init%netcdf%var_Hi)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_Hb,     (/ dim_x, dim_y /),  init%netcdf%var_Hb)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_Hs,     (/ dim_x, dim_y /),  init%netcdf%var_Hs)
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE inquire_init_data_file_cart
  SUBROUTINE inquire_init_data_file_mesh(init, mesh, netcdf, nt)
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    TYPE(type_mesh),                INTENT(INOUT) :: mesh
    TYPE(type_netcdf_file_output),  INTENT(INOUT) :: netcdf
    INTEGER,                        INTENT(OUT)   :: nt
 
    ! Local variables:
    INTEGER                               :: int_dummy, dim_nV, dim_nV2, dim_nTri, dim_nTri2, dim_nconmax, dim_zeta, dim_month, dim_time
        
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( netcdf%ncid, netcdf%name_mesh_V,        mesh%nV,             dim_nV)
    CALL inquire_dim( netcdf%ncid, netcdf%name_mesh_V2,       int_dummy,           dim_nV2)
    CALL inquire_dim( netcdf%ncid, netcdf%name_mesh_Tri,      mesh%nTri,           dim_nTri)
    CALL inquire_dim( netcdf%ncid, netcdf%name_mesh_Tri2,     int_dummy,           dim_nTri2)
    CALL inquire_dim( netcdf%ncid, 'neighbour',               mesh%nconmax,        dim_nconmax)
    CALL inquire_dim( netcdf%ncid, netcdf%name_zeta,          int_dummy,           dim_zeta)
    CALL inquire_dim( netcdf%ncid, netcdf%name_month,         int_dummy,           dim_month)
    CALL inquire_dim( netcdf%ncid, netcdf%name_time,          nt,                  dim_time)
    
    mesh%nV_max   = mesh%nV
    mesh%nTri_max = mesh%nTri

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( netcdf%ncid, netcdf%name_mesh_V,              (/ dim_nV,   dim_nV2     /),      netcdf%var_mesh_V)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_mesh_Tri,            (/ dim_nTri, dim_nTri2   /),      netcdf%var_mesh_Tri)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_zeta,                (/ dim_zeta              /),      netcdf%var_zeta)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_month,               (/ dim_month             /),      netcdf%var_month)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_time,                (/ dim_time              /),      netcdf%var_time)
    
    CALL inquire_double_var( netcdf%ncid, netcdf%name_mesh_nC,             (/ dim_nV                /),      netcdf%var_mesh_nC)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_mesh_C,              (/ dim_nV,   dim_nconmax /),      netcdf%var_mesh_C)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_mesh_edge_index,     (/ dim_nV                /),      netcdf%var_mesh_edge_index)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_mesh_niTri,          (/ dim_nV                /),      netcdf%var_mesh_niTri)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_mesh_iTri,           (/ dim_nV,   dim_nconmax /),      netcdf%var_mesh_iTri)
    
    CALL inquire_double_var( netcdf%ncid, netcdf%name_mesh_Tricc,          (/ dim_nTri, dim_nV2     /),      netcdf%var_mesh_Tricc)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_mesh_TriC,           (/ dim_nTri, dim_nTri2   /),      netcdf%var_mesh_TriC)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_mesh_Tri_edge_index, (/ dim_nTri              /),      netcdf%var_mesh_Tri_edge_index)
        
    CALL inquire_double_var( netcdf%ncid, netcdf%name_Hi,                  (/ dim_nV,   dim_time    /),      netcdf%var_Hi)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_Hb,                  (/ dim_nV,   dim_time    /),      netcdf%var_Hb)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_Hs,                  (/ dim_nV,   dim_time    /),      netcdf%var_Hs)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_U_SSA,               (/ dim_nV,   dim_time    /),      netcdf%var_U_SSA)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_V_SSA,               (/ dim_nV,   dim_time    /),      netcdf%var_V_SSA)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_Ti,                  (/ dim_nV, dim_zeta,  dim_time /), netcdf%var_Ti)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_Albedo,              (/ dim_nV, dim_month, dim_time /), netcdf%var_Albedo)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_FirnDepth,           (/ dim_nV, dim_month, dim_time /), netcdf%var_FirnDepth)
    CALL inquire_double_var( netcdf%ncid, netcdf%name_MeltPreviousYear,    (/ dim_nV,   dim_time    /),      netcdf%var_MeltPreviousYear)
        
    ! Close the netcdf file
    CALL close_netcdf_file(netcdf%ncid)
    
  END SUBROUTINE inquire_init_data_file_mesh
  SUBROUTINE inquire_PD_obs_data_file(PD_obs) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: PD_obs
 
    ! Local variables:  
    INTEGER                               :: dim_lat
    INTEGER                               :: dim_lon
    INTEGER                               :: dim_month
    INTEGER                               :: int_dummy
        
    ! Open the netcdf file
    CALL open_netcdf_file(PD_obs%netcdf%filename, PD_obs%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_lat,      PD_obs%nlat,             dim_lat)
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_lon,      PD_obs%nlon,             dim_lon)
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_month,    int_dummy,               dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_lat,      (/ dim_lat /),                      PD_obs%netcdf%var_lat)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_lon,      (/ dim_lon /),                      PD_obs%netcdf%var_lon)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_T2m,      (/ dim_lon, dim_lat, dim_month /),  PD_obs%netcdf%var_T2m)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_Precip,   (/ dim_lon, dim_lat, dim_month /),  PD_obs%netcdf%var_Precip)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_Hs,       (/ dim_lon, dim_lat            /),  PD_obs%netcdf%var_Hs)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_Albedo,   (/ dim_lon, dim_lat, dim_month /),  PD_obs%netcdf%var_Albedo)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_Wind_WE,  (/ dim_lon, dim_lat, dim_month /),  PD_obs%netcdf%var_Wind_WE)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_Wind_SN,  (/ dim_lon, dim_lat, dim_month /),  PD_obs%netcdf%var_Wind_SN)
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD_obs%netcdf%ncid)
    
  END SUBROUTINE inquire_PD_obs_data_file 
  SUBROUTINE inquire_insolation_data_file(forcing)
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
 
    ! Local variables:  
    INTEGER                                :: dim_lat
    INTEGER                                :: dim_time
    INTEGER                                :: dim_month
    INTEGER                                :: int_dummy
            
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%ins_netcdf%filename, forcing%ins_netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( forcing%ins_netcdf%ncid, forcing%ins_netcdf%name_time,     forcing%ins_nyears,        dim_time)
    CALL inquire_dim( forcing%ins_netcdf%ncid, forcing%ins_netcdf%name_month,    int_dummy,                 dim_month)  
    CALL inquire_dim( forcing%ins_netcdf%ncid, forcing%ins_netcdf%name_lat,      forcing%ins_nlat,          dim_lat)
    
    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( forcing%ins_netcdf%ncid, forcing%ins_netcdf%name_time,     (/ dim_time /),                      forcing%ins_netcdf%var_time)
    CALL inquire_double_var( forcing%ins_netcdf%ncid, forcing%ins_netcdf%name_month,    (/ dim_month /),                     forcing%ins_netcdf%var_month)
    CALL inquire_double_var( forcing%ins_netcdf%ncid, forcing%ins_netcdf%name_lat,      (/ dim_lat /),                       forcing%ins_netcdf%var_lat)
    CALL inquire_double_var( forcing%ins_netcdf%ncid, forcing%ins_netcdf%name_Q_TOA,    (/ dim_time, dim_month, dim_lat /),  forcing%ins_netcdf%var_Q_TOA)
        
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%ins_netcdf%ncid)
    
  END SUBROUTINE inquire_insolation_data_file

  SUBROUTINE GetOutputFilename(region)
   
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

  SUBROUTINE open_netcdf_file(filename, ncid) 
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)  :: filename
    INTEGER,          INTENT(OUT) :: ncid
    
    ! Open netCDF file:
    CALL handle_error(nf90_open(filename, IOR(nf90_write,nf90_share), ncid))

  END SUBROUTINE open_netcdf_file
  SUBROUTINE close_netcdf_file(ncid)
    IMPLICIT NONE
  
    INTEGER, INTENT(INOUT) :: ncid

    ! Close netCDF file:
    CALL handle_error(nf90_close(ncid))
   
  END SUBROUTINE close_netcdf_file
  SUBROUTINE create_dim(ncid, dim_name, length, id_dim)
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
  SUBROUTINE create_double_var(ncid, var_name, id_dims, id_var, long_name, units, missing_value)
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
  SUBROUTINE inquire_dim(ncid, dim_name, dim_length, id_dim)
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
  SUBROUTINE inquire_double_var(ncid, var_name, id_dims, id_var)
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
  SUBROUTINE handle_error(stat, message)
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
