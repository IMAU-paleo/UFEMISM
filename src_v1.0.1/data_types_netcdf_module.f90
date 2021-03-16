MODULE data_types_netcdf_module
  ! Contains the TYPES for different NetCDf files read and written by UFEMISM.

  USE configuration_module,        ONLY: dp, C

  IMPLICIT NONE
    
  TYPE type_netcdf_output  
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! Index of time frame to be written to
    INTEGER :: ti
    
  ! ID's for dimensions:
  ! ===================
    
    ! Dimensions
    INTEGER :: id_dim_vi
    INTEGER :: id_dim_ti
    INTEGER :: id_dim_ci
    INTEGER :: id_dim_aci
    INTEGER :: id_dim_ciplusone
    INTEGER :: id_dim_two
    INTEGER :: id_dim_three
    INTEGER :: id_dim_four
    INTEGER :: id_dim_zeta
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
  
    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_aci                   = 'aci                  '
    CHARACTER(LEN=256) :: name_dim_ciplusone             = 'ciplusone            '
    CHARACTER(LEN=256) :: name_dim_two                   = 'two                  '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_four                  = 'four                 '
    CHARACTER(LEN=256) :: name_dim_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    
  ! Primary output - everything that's needed to restart a new run, and to plot/analyse data
  ! ========================================================================================
  
    INTEGER :: id_var_V
    INTEGER :: id_var_Tri
    INTEGER :: id_var_nC
    INTEGER :: id_var_C
    INTEGER :: id_var_niTri
    INTEGER :: id_var_iTri
    INTEGER :: id_var_edge_index
    INTEGER :: id_var_Tricc
    INTEGER :: id_var_TriC
    INTEGER :: id_var_Tri_edge_index
    INTEGER :: id_var_R
    INTEGER :: id_var_lat
    INTEGER :: id_var_lon
  
    INTEGER :: id_var_time
    INTEGER :: id_var_zeta
    INTEGER :: id_var_month
    
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hb
    INTEGER :: id_var_Hs
    INTEGER :: id_var_U_SIA
    INTEGER :: id_var_V_SIA
    INTEGER :: id_var_U_SSA
    INTEGER :: id_var_V_SSA
    INTEGER :: id_var_Ti
    INTEGER :: id_var_mask
    INTEGER :: id_var_T2m
    INTEGER :: id_var_Precip
    INTEGER :: id_var_Albedo
    INTEGER :: id_var_SMB
    INTEGER :: id_var_BMB
    INTEGER :: id_var_FirnDepth
    INTEGER :: id_var_MeltPreviousYear 
    INTEGER :: id_var_dHs_dx
    INTEGER :: id_var_dHs_dy 
    INTEGER :: id_var_D_SIA  
  
    CHARACTER(LEN=256) :: name_var_V                     = 'V                    '
    CHARACTER(LEN=256) :: name_var_Tri                   = 'Tri                  '    
    CHARACTER(LEN=256) :: name_var_nC                    = 'nC                   '
    CHARACTER(LEN=256) :: name_var_C                     = 'C                    '
    CHARACTER(LEN=256) :: name_var_niTri                 = 'niTri                '
    CHARACTER(LEN=256) :: name_var_iTri                  = 'iTri                 '
    CHARACTER(LEN=256) :: name_var_edge_index            = 'edge_index           '
    CHARACTER(LEN=256) :: name_var_Tricc                 = 'Tricc                '
    CHARACTER(LEN=256) :: name_var_TriC                  = 'TriC                 '
    CHARACTER(LEN=256) :: name_var_Tri_edge_index        = 'Tri_edge_index       '
    CHARACTER(LEN=256) :: name_var_R                     = 'R                    '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
    
    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   ' 
    CHARACTER(LEN=256) :: name_var_U_SIA                 = 'U_SIA                '
    CHARACTER(LEN=256) :: name_var_V_SIA                 = 'V_SIA                '
    CHARACTER(LEN=256) :: name_var_U_SSA                 = 'U_SSA                '
    CHARACTER(LEN=256) :: name_var_V_SSA                 = 'V_SSA                '
    CHARACTER(LEN=256) :: name_var_Ti                    = 'Ti                   '
    CHARACTER(LEN=256) :: name_var_mask                  = 'mask                 '
    CHARACTER(LEN=256) :: name_var_T2m                   = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_Precip                = 'Precip               '
    CHARACTER(LEN=256) :: name_var_Albedo                = 'Albedo               '
    CHARACTER(LEN=256) :: name_var_SMB                   = 'SMB                  '
    CHARACTER(LEN=256) :: name_var_BMB                   = 'BMB                  '
    CHARACTER(LEN=256) :: name_var_FirnDepth             = 'FirnDepth            '
    CHARACTER(LEN=256) :: name_var_MeltPreviousYear      = 'MeltPreviousYear     '
    CHARACTER(LEN=256) :: name_var_dHs_dx                = 'dHs_dx               '
    CHARACTER(LEN=256) :: name_var_dHs_dy                = 'dHs_dy               '
    CHARACTER(LEN=256) :: name_var_D_SIA                 = 'D_SIA                '
    
  ! Secondary output - mesh data
  ! ============================
    
    INTEGER :: id_var_A
    INTEGER :: id_var_Cw
    INTEGER :: id_var_TriA
    INTEGER :: id_var_NxTri
    INTEGER :: id_var_NyTri
    INTEGER :: id_var_Nx
    INTEGER :: id_var_Ny
    INTEGER :: id_var_Nxx
    INTEGER :: id_var_Nxy
    INTEGER :: id_var_Nyy
    INTEGER :: id_var_VAc
    INTEGER :: id_var_Aci
    INTEGER :: id_var_iAci
    
    CHARACTER(LEN=256) :: name_var_A                     = 'A                    '
    CHARACTER(LEN=256) :: name_var_Cw                    = 'Cw                   '
    CHARACTER(LEN=256) :: name_var_TriA                  = 'TriA                 '
    CHARACTER(LEN=256) :: name_var_NxTri                 = 'NxTri                '
    CHARACTER(LEN=256) :: name_var_NyTri                 = 'NyTri                '
    CHARACTER(LEN=256) :: name_var_Nx                    = 'Nx                   '
    CHARACTER(LEN=256) :: name_var_Ny                    = 'Ny                   '
    CHARACTER(LEN=256) :: name_var_Nxx                   = 'Nxx                  '
    CHARACTER(LEN=256) :: name_var_Nxy                   = 'Nxy                  '
    CHARACTER(LEN=256) :: name_var_Nyy                   = 'Nyy                  '
    CHARACTER(LEN=256) :: name_var_VAc                   = 'VAc                  '
    CHARACTER(LEN=256) :: name_var_Aci                   = 'Aci                  '
    CHARACTER(LEN=256) :: name_var_iAci                  = 'iAci                 '
    
  ! Secondary output - ice dynamics
  ! ===============================
  
    INTEGER :: id_var_dHs_dx_ac
    INTEGER :: id_var_dHs_dy_ac
    INTEGER :: id_var_D_SIA_ac
    INTEGER :: id_var_Ux_ac
    INTEGER :: id_var_Uy_ac
    INTEGER :: id_var_U_3D
    INTEGER :: id_var_V_3D
    INTEGER :: id_var_W_3D
    
    CHARACTER(LEN=256) :: name_var_dHs_dx_ac             = 'dHs_dx_ac            '
    CHARACTER(LEN=256) :: name_var_dHs_dy_ac             = 'dHs_dy_ac            '
    CHARACTER(LEN=256) :: name_var_D_SIA_ac              = 'D_SIA_ac             '
    CHARACTER(LEN=256) :: name_var_Ux_ac                 = 'Ux_ac                '
    CHARACTER(LEN=256) :: name_var_Uy_ac                 = 'Uy_ac                '
    CHARACTER(LEN=256) :: name_var_U_3D                  = 'U_3D                 '
    CHARACTER(LEN=256) :: name_var_V_3D                  = 'V_3D                 '
    CHARACTER(LEN=256) :: name_var_W_3D                  = 'W_3D                 '  
    
  ! Secondary output - climate and SMB components
  ! =============================================
  
    INTEGER :: id_var_Wind_WE
    INTEGER :: id_var_Wind_SN
    INTEGER :: id_var_Snowfall
    INTEGER :: id_var_Melt
    INTEGER :: id_var_Refreezing
    INTEGER :: id_var_Runoff
    
    CHARACTER(LEN=256) :: name_var_Wind_WE               = 'Wind_WE              '
    CHARACTER(LEN=256) :: name_var_Wind_SN               = 'Wind_SN              '
    CHARACTER(LEN=256) :: name_var_Snowfall              = 'Snowfall             '
    CHARACTER(LEN=256) :: name_var_Melt                  = 'Melt                 '
    CHARACTER(LEN=256) :: name_var_Refreezing            = 'Refreezing           '
    CHARACTER(LEN=256) :: name_var_Runoff                = 'Runoff               ' 
        
  END TYPE type_netcdf_output
    
  TYPE type_netcdf_PD_data
    ! For reading an input file describing a present-day model region, on a Cartesian grid
  
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    ! Dimensions
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_month
    
    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    
    ! Variables:
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hb
    INTEGER :: id_var_Hs
    
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   ' 
        
  END TYPE type_netcdf_PD_data
    
  TYPE type_netcdf_init_data
    ! For reading an input file describing the initial state of a model region, on a Cartesian grid
    
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    
    ! Variable names
    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    
    ! Key output
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hb
    INTEGER :: id_var_Hs
    
    ! Variable names
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   '
        
  END TYPE type_netcdf_init_data
    
  TYPE type_netcdf_climate_data
    ! For reading an input file containing either a GCM snapshot or a PD observations data set (e.g. ERA-40),
    ! describing the global climate with monthly fields on a lat/lon grid
  
    ! Integers describing open ports to different variables in an opened NetCDF file.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    ! Dimensions
    INTEGER :: id_dim_lon
    INTEGER :: id_dim_lat
    INTEGER :: id_dim_month
    
    CHARACTER(LEN=256) :: name_dim_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    
    ! Variables
    INTEGER :: id_var_lon
    INTEGER :: id_var_lat
    INTEGER :: id_var_Hs
    INTEGER :: id_var_T2m
    INTEGER :: id_var_Precip
    INTEGER :: id_var_Albedo
    INTEGER :: id_var_Wind_WE
    INTEGER :: id_var_Wind_SN
    
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   '
    CHARACTER(LEN=256) :: name_var_T2m                   = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_Precip                = 'Precip               '
    CHARACTER(LEN=256) :: name_var_Albedo                = 'Albedo               '
    CHARACTER(LEN=256) :: name_var_Wind_WE               = 'Wind_WE              '
    CHARACTER(LEN=256) :: name_var_Wind_SN               = 'Wind_SN              '
        
  END TYPE type_netcdf_climate_data
    
  TYPE type_netcdf_insolation
    ! For reading an input file containing an insolation history reconstruction (e.g. Lasker et al., 2004),
    ! describing top-of-the-atmosphere insolation for every month of the year at a latitudinal grid.
  
    ! Integers describing open ports to different variables in an opened NetCDF file.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    ! Dimensions
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
    INTEGER :: id_dim_lat
    
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '
    
    ! Variables
    INTEGER :: id_var_time
    INTEGER :: id_var_month
    INTEGER :: id_var_lat
    INTEGER :: id_var_Q_TOA
    
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_Q_TOA                 = 'Q_TOA                '
        
  END TYPE type_netcdf_insolation
  
CONTAINS

END MODULE data_types_netcdf_module
