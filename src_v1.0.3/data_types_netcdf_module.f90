MODULE data_types_netcdf_module
  ! Contains the TYPES for different NetCDf files read and written by UFEMISM.

  USE configuration_module,        ONLY: dp, C

  IMPLICIT NONE
    
  TYPE type_netcdf_restart
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! Index of time frame to be written to
    INTEGER :: ti
    
    ! Mesh data
    ! ==========
    
    INTEGER :: id_dim_vi
    INTEGER :: id_dim_ti
    INTEGER :: id_dim_ci
    INTEGER :: id_dim_aci
    INTEGER :: id_dim_ai
    INTEGER :: id_dim_tai
    INTEGER :: id_dim_ciplusone
    INTEGER :: id_dim_two
    INTEGER :: id_dim_three
    INTEGER :: id_dim_four
    INTEGER :: id_dim_vii_transect
  
    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_aci                   = 'aci                  '
    CHARACTER(LEN=256) :: name_dim_ciplusone             = 'ciplusone            '
    CHARACTER(LEN=256) :: name_dim_two                   = 'two                  '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_four                  = 'four                 '
    CHARACTER(LEN=256) :: name_dim_vii_transect          = 'vii                  '
    CHARACTER(LEN=256) :: name_dim_ai                    = 'ai                   '
    CHARACTER(LEN=256) :: name_dim_tai                   = 'tai                  '
  
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
    INTEGER :: id_var_VAc
    INTEGER :: id_var_Aci
    INTEGER :: id_var_iAci
    INTEGER :: id_var_VAaAc
    INTEGER :: id_var_TriAaAc
    INTEGER :: id_var_A
    INTEGER :: id_var_R
    INTEGER :: id_var_vi_transect
    INTEGER :: id_var_w_transect
  
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
    CHARACTER(LEN=256) :: name_var_VAc                   = 'VAc                  '
    CHARACTER(LEN=256) :: name_var_Aci                   = 'Aci                  '
    CHARACTER(LEN=256) :: name_var_iAci                  = 'iAci                 '
    CHARACTER(LEN=256) :: name_var_VAaAc                 = 'VAaAc                '
    CHARACTER(LEN=256) :: name_var_TriAaAc               = 'TriAaAc              '
    CHARACTER(LEN=256) :: name_var_A                     = 'A                    '
    CHARACTER(LEN=256) :: name_var_R                     = 'R                    '
    CHARACTER(LEN=256) :: name_var_vi_transect           = 'vi_transect          '
    CHARACTER(LEN=256) :: name_var_w_transect            = 'w_transect           '
    
    ! Grid data
    ! =========
    
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    
    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    
    ! Data dimensions
    ! ===============
    
    INTEGER :: id_dim_zeta
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
    
    CHARACTER(LEN=256) :: name_dim_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
  
    INTEGER :: id_var_time
    INTEGER :: id_var_zeta
    INTEGER :: id_var_month
    
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
    
    ! Variables
    ! =========
    
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hb
    INTEGER :: id_var_Hs
    INTEGER :: id_var_U_SIA
    INTEGER :: id_var_V_SIA
    INTEGER :: id_var_U_SSA
    INTEGER :: id_var_V_SSA
    INTEGER :: id_var_Ti
    INTEGER :: id_var_FirnDepth
    INTEGER :: id_var_MeltPreviousYear  
    
    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   ' 
    CHARACTER(LEN=256) :: name_var_U_SIA                 = 'U_SIA                '
    CHARACTER(LEN=256) :: name_var_V_SIA                 = 'V_SIA                '
    CHARACTER(LEN=256) :: name_var_U_SSA                 = 'U_SSA                '
    CHARACTER(LEN=256) :: name_var_V_SSA                 = 'V_SSA                '
    CHARACTER(LEN=256) :: name_var_Ti                    = 'Ti                   '
    CHARACTER(LEN=256) :: name_var_FirnDepth             = 'FirnDepth            '
    CHARACTER(LEN=256) :: name_var_MeltPreviousYear      = 'MeltPreviousYear     '
        
  END TYPE type_netcdf_restart
    
  TYPE type_netcdf_help_fields
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! Index of time frame to be written to
    INTEGER :: ti
    
    ! Mesh data
    ! ==========
    
    INTEGER :: id_dim_vi
    INTEGER :: id_dim_ti
    INTEGER :: id_dim_ci
    INTEGER :: id_dim_aci
    INTEGER :: id_dim_ai
    INTEGER :: id_dim_tai
    INTEGER :: id_dim_ciplusone
    INTEGER :: id_dim_two
    INTEGER :: id_dim_three
    INTEGER :: id_dim_four
    INTEGER :: id_dim_vii_transect
  
    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_aci                   = 'aci                  '
    CHARACTER(LEN=256) :: name_dim_ciplusone             = 'ciplusone            '
    CHARACTER(LEN=256) :: name_dim_two                   = 'two                  '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_four                  = 'four                 '
    CHARACTER(LEN=256) :: name_dim_vii_transect          = 'vii                  '
    CHARACTER(LEN=256) :: name_dim_ai                    = 'ai                   '
    CHARACTER(LEN=256) :: name_dim_tai                   = 'tai                  '
  
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
    INTEGER :: id_var_VAc
    INTEGER :: id_var_Aci
    INTEGER :: id_var_iAci
    INTEGER :: id_var_VAaAc
    INTEGER :: id_var_TriAaAc
    INTEGER :: id_var_A
    INTEGER :: id_var_R
    INTEGER :: id_var_vi_transect
    INTEGER :: id_var_w_transect
  
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
    CHARACTER(LEN=256) :: name_var_VAc                   = 'VAc                  '
    CHARACTER(LEN=256) :: name_var_Aci                   = 'Aci                  '
    CHARACTER(LEN=256) :: name_var_iAci                  = 'iAci                 '
    CHARACTER(LEN=256) :: name_var_VAaAc                 = 'VAaAc                '
    CHARACTER(LEN=256) :: name_var_TriAaAc               = 'TriAaAc              '
    CHARACTER(LEN=256) :: name_var_A                     = 'A                    '
    CHARACTER(LEN=256) :: name_var_R                     = 'R                    '
    CHARACTER(LEN=256) :: name_var_vi_transect           = 'vi_transect          '
    CHARACTER(LEN=256) :: name_var_w_transect            = 'w_transect           '
    
    ! Grid data
    ! =========
    
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    
    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    
    ! Data dimensions
    ! ===============
    
    INTEGER :: id_dim_zeta
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
    
    CHARACTER(LEN=256) :: name_dim_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
  
    INTEGER :: id_var_time
    INTEGER :: id_var_zeta
    INTEGER :: id_var_month
    
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
    
    ! Variables
    ! =========
    
    INTEGER :: id_help_field_01
    INTEGER :: id_help_field_02
    INTEGER :: id_help_field_03
    INTEGER :: id_help_field_04
    INTEGER :: id_help_field_05
    INTEGER :: id_help_field_06
    INTEGER :: id_help_field_07
    INTEGER :: id_help_field_08
    INTEGER :: id_help_field_09
    INTEGER :: id_help_field_10
    INTEGER :: id_help_field_11
    INTEGER :: id_help_field_12
    INTEGER :: id_help_field_13
    INTEGER :: id_help_field_14
    INTEGER :: id_help_field_15
    INTEGER :: id_help_field_16
    INTEGER :: id_help_field_17
    INTEGER :: id_help_field_18
    INTEGER :: id_help_field_19
    INTEGER :: id_help_field_20
    INTEGER :: id_help_field_21
    INTEGER :: id_help_field_22
    INTEGER :: id_help_field_23
    INTEGER :: id_help_field_24
    INTEGER :: id_help_field_25
    INTEGER :: id_help_field_26
    INTEGER :: id_help_field_27
    INTEGER :: id_help_field_28
    INTEGER :: id_help_field_29
    INTEGER :: id_help_field_30
    INTEGER :: id_help_field_31
    INTEGER :: id_help_field_32
    INTEGER :: id_help_field_33
    INTEGER :: id_help_field_34
    INTEGER :: id_help_field_35
    INTEGER :: id_help_field_36
    INTEGER :: id_help_field_37
    INTEGER :: id_help_field_38
    INTEGER :: id_help_field_39
    INTEGER :: id_help_field_40
    INTEGER :: id_help_field_41
    INTEGER :: id_help_field_42
    INTEGER :: id_help_field_43
    INTEGER :: id_help_field_44
    INTEGER :: id_help_field_45
    INTEGER :: id_help_field_46
    INTEGER :: id_help_field_47
    INTEGER :: id_help_field_48
    INTEGER :: id_help_field_49
    INTEGER :: id_help_field_50
        
  END TYPE type_netcdf_help_fields
    
  TYPE type_netcdf_debug
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! Mesh data
    ! ==========
    
    INTEGER :: id_dim_vi
    INTEGER :: id_dim_ti
    INTEGER :: id_dim_ci
    INTEGER :: id_dim_aci
    INTEGER :: id_dim_ai
    INTEGER :: id_dim_tai
    INTEGER :: id_dim_ciplusone
    INTEGER :: id_dim_two
    INTEGER :: id_dim_three
    INTEGER :: id_dim_four
    INTEGER :: id_dim_vii_transect
  
    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_aci                   = 'aci                  '
    CHARACTER(LEN=256) :: name_dim_ciplusone             = 'ciplusone            '
    CHARACTER(LEN=256) :: name_dim_two                   = 'two                  '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_four                  = 'four                 '
    CHARACTER(LEN=256) :: name_dim_vii_transect          = 'vii                  '
    CHARACTER(LEN=256) :: name_dim_ai                    = 'ai                   '
    CHARACTER(LEN=256) :: name_dim_tai                   = 'tai                  '
  
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
    INTEGER :: id_var_VAc
    INTEGER :: id_var_Aci
    INTEGER :: id_var_iAci
    INTEGER :: id_var_VAaAc
    INTEGER :: id_var_TriAaAc
    INTEGER :: id_var_A
    INTEGER :: id_var_R
    INTEGER :: id_var_vi_transect
    INTEGER :: id_var_w_transect
  
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
    CHARACTER(LEN=256) :: name_var_VAc                   = 'VAc                  '
    CHARACTER(LEN=256) :: name_var_Aci                   = 'Aci                  '
    CHARACTER(LEN=256) :: name_var_iAci                  = 'iAci                 '
    CHARACTER(LEN=256) :: name_var_VAaAc                 = 'VAaAc                '
    CHARACTER(LEN=256) :: name_var_TriAaAc               = 'TriAaAc              '
    CHARACTER(LEN=256) :: name_var_A                     = 'A                    '
    CHARACTER(LEN=256) :: name_var_R                     = 'R                    '
    CHARACTER(LEN=256) :: name_var_vi_transect           = 'vi_transect          '
    CHARACTER(LEN=256) :: name_var_w_transect            = 'w_transect           '
    
    ! Data dimensions
    ! ===============
    
    INTEGER :: id_dim_zeta
    INTEGER :: id_dim_month
    
    CHARACTER(LEN=256) :: name_dim_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
  
    INTEGER :: id_var_zeta
    INTEGER :: id_var_month
    
    CHARACTER(LEN=256) :: name_var_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
    
    ! Variables
    ! =========
    
    INTEGER :: id_var_int_2D_Aa_01
    INTEGER :: id_var_int_2D_Aa_02
    INTEGER :: id_var_int_2D_Aa_03
    INTEGER :: id_var_int_2D_Aa_04
    INTEGER :: id_var_int_2D_Aa_05
    INTEGER :: id_var_int_2D_Aa_06
    INTEGER :: id_var_int_2D_Aa_07
    INTEGER :: id_var_int_2D_Aa_08
    INTEGER :: id_var_int_2D_Aa_09
    INTEGER :: id_var_int_2D_Aa_10
    
    CHARACTER(LEN=256) :: name_var_int_2D_Aa_01        = 'int_2D_Aa_01      '
    CHARACTER(LEN=256) :: name_var_int_2D_Aa_02        = 'int_2D_Aa_02      '
    CHARACTER(LEN=256) :: name_var_int_2D_Aa_03        = 'int_2D_Aa_03      '
    CHARACTER(LEN=256) :: name_var_int_2D_Aa_04        = 'int_2D_Aa_04      '
    CHARACTER(LEN=256) :: name_var_int_2D_Aa_05        = 'int_2D_Aa_05      '
    CHARACTER(LEN=256) :: name_var_int_2D_Aa_06        = 'int_2D_Aa_06      '
    CHARACTER(LEN=256) :: name_var_int_2D_Aa_07        = 'int_2D_Aa_07      '
    CHARACTER(LEN=256) :: name_var_int_2D_Aa_08        = 'int_2D_Aa_08      '
    CHARACTER(LEN=256) :: name_var_int_2D_Aa_09        = 'int_2D_Aa_09      '
    CHARACTER(LEN=256) :: name_var_int_2D_Aa_10        = 'int_2D_Aa_10      '
    
    INTEGER :: id_var_int_2D_Ac_01
    INTEGER :: id_var_int_2D_Ac_02
    INTEGER :: id_var_int_2D_Ac_03
    INTEGER :: id_var_int_2D_Ac_04
    INTEGER :: id_var_int_2D_Ac_05
    INTEGER :: id_var_int_2D_Ac_06
    INTEGER :: id_var_int_2D_Ac_07
    INTEGER :: id_var_int_2D_Ac_08
    INTEGER :: id_var_int_2D_Ac_09
    INTEGER :: id_var_int_2D_Ac_10
    
    CHARACTER(LEN=256) :: name_var_int_2D_Ac_01        = 'int_2D_Ac_01      '
    CHARACTER(LEN=256) :: name_var_int_2D_Ac_02        = 'int_2D_Ac_02      '
    CHARACTER(LEN=256) :: name_var_int_2D_Ac_03        = 'int_2D_Ac_03      '
    CHARACTER(LEN=256) :: name_var_int_2D_Ac_04        = 'int_2D_Ac_04      '
    CHARACTER(LEN=256) :: name_var_int_2D_Ac_05        = 'int_2D_Ac_05      '
    CHARACTER(LEN=256) :: name_var_int_2D_Ac_06        = 'int_2D_Ac_06      '
    CHARACTER(LEN=256) :: name_var_int_2D_Ac_07        = 'int_2D_Ac_07      '
    CHARACTER(LEN=256) :: name_var_int_2D_Ac_08        = 'int_2D_Ac_08      '
    CHARACTER(LEN=256) :: name_var_int_2D_Ac_09        = 'int_2D_Ac_09      '
    CHARACTER(LEN=256) :: name_var_int_2D_Ac_10        = 'int_2D_Ac_10      '
    
    INTEGER :: id_var_int_2D_AaAc_01
    INTEGER :: id_var_int_2D_AaAc_02
    INTEGER :: id_var_int_2D_AaAc_03
    INTEGER :: id_var_int_2D_AaAc_04
    INTEGER :: id_var_int_2D_AaAc_05
    INTEGER :: id_var_int_2D_AaAc_06
    INTEGER :: id_var_int_2D_AaAc_07
    INTEGER :: id_var_int_2D_AaAc_08
    INTEGER :: id_var_int_2D_AaAc_09
    INTEGER :: id_var_int_2D_AaAc_10
    
    CHARACTER(LEN=256) :: name_var_int_2D_AaAc_01      = 'int_2D_AaAc_01      '
    CHARACTER(LEN=256) :: name_var_int_2D_AaAc_02      = 'int_2D_AaAc_02      '
    CHARACTER(LEN=256) :: name_var_int_2D_AaAc_03      = 'int_2D_AaAc_03      '
    CHARACTER(LEN=256) :: name_var_int_2D_AaAc_04      = 'int_2D_AaAc_04      '
    CHARACTER(LEN=256) :: name_var_int_2D_AaAc_05      = 'int_2D_AaAc_05      '
    CHARACTER(LEN=256) :: name_var_int_2D_AaAc_06      = 'int_2D_AaAc_06      '
    CHARACTER(LEN=256) :: name_var_int_2D_AaAc_07      = 'int_2D_AaAc_07      '
    CHARACTER(LEN=256) :: name_var_int_2D_AaAc_08      = 'int_2D_AaAc_08      '
    CHARACTER(LEN=256) :: name_var_int_2D_AaAc_09      = 'int_2D_AaAc_09      '
    CHARACTER(LEN=256) :: name_var_int_2D_AaAc_10      = 'int_2D_AaAc_10      '
    
    INTEGER :: id_var_dp_2D_Aa_01
    INTEGER :: id_var_dp_2D_Aa_02
    INTEGER :: id_var_dp_2D_Aa_03
    INTEGER :: id_var_dp_2D_Aa_04
    INTEGER :: id_var_dp_2D_Aa_05
    INTEGER :: id_var_dp_2D_Aa_06
    INTEGER :: id_var_dp_2D_Aa_07
    INTEGER :: id_var_dp_2D_Aa_08
    INTEGER :: id_var_dp_2D_Aa_09
    INTEGER :: id_var_dp_2D_Aa_10
    
    CHARACTER(LEN=256) :: name_var_dp_2D_Aa_01        = 'dp_2D_Aa_01      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Aa_02        = 'dp_2D_Aa_02      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Aa_03        = 'dp_2D_Aa_03      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Aa_04        = 'dp_2D_Aa_04      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Aa_05        = 'dp_2D_Aa_05      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Aa_06        = 'dp_2D_Aa_06      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Aa_07        = 'dp_2D_Aa_07      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Aa_08        = 'dp_2D_Aa_08      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Aa_09        = 'dp_2D_Aa_09      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Aa_10        = 'dp_2D_Aa_10      '
    
    INTEGER :: id_var_dp_2D_Ac_01
    INTEGER :: id_var_dp_2D_Ac_02
    INTEGER :: id_var_dp_2D_Ac_03
    INTEGER :: id_var_dp_2D_Ac_04
    INTEGER :: id_var_dp_2D_Ac_05
    INTEGER :: id_var_dp_2D_Ac_06
    INTEGER :: id_var_dp_2D_Ac_07
    INTEGER :: id_var_dp_2D_Ac_08
    INTEGER :: id_var_dp_2D_Ac_09
    INTEGER :: id_var_dp_2D_Ac_10
    
    CHARACTER(LEN=256) :: name_var_dp_2D_Ac_01        = 'dp_2D_Ac_01      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Ac_02        = 'dp_2D_Ac_02      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Ac_03        = 'dp_2D_Ac_03      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Ac_04        = 'dp_2D_Ac_04      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Ac_05        = 'dp_2D_Ac_05      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Ac_06        = 'dp_2D_Ac_06      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Ac_07        = 'dp_2D_Ac_07      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Ac_08        = 'dp_2D_Ac_08      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Ac_09        = 'dp_2D_Ac_09      '
    CHARACTER(LEN=256) :: name_var_dp_2D_Ac_10        = 'dp_2D_Ac_10      '
    
    INTEGER :: id_var_dp_2D_AaAc_01
    INTEGER :: id_var_dp_2D_AaAc_02
    INTEGER :: id_var_dp_2D_AaAc_03
    INTEGER :: id_var_dp_2D_AaAc_04
    INTEGER :: id_var_dp_2D_AaAc_05
    INTEGER :: id_var_dp_2D_AaAc_06
    INTEGER :: id_var_dp_2D_AaAc_07
    INTEGER :: id_var_dp_2D_AaAc_08
    INTEGER :: id_var_dp_2D_AaAc_09
    INTEGER :: id_var_dp_2D_AaAc_10
    
    CHARACTER(LEN=256) :: name_var_dp_2D_AaAc_01      = 'dp_2D_AaAc_01      '
    CHARACTER(LEN=256) :: name_var_dp_2D_AaAc_02      = 'dp_2D_AaAc_02      '
    CHARACTER(LEN=256) :: name_var_dp_2D_AaAc_03      = 'dp_2D_AaAc_03      '
    CHARACTER(LEN=256) :: name_var_dp_2D_AaAc_04      = 'dp_2D_AaAc_04      '
    CHARACTER(LEN=256) :: name_var_dp_2D_AaAc_05      = 'dp_2D_AaAc_05      '
    CHARACTER(LEN=256) :: name_var_dp_2D_AaAc_06      = 'dp_2D_AaAc_06      '
    CHARACTER(LEN=256) :: name_var_dp_2D_AaAc_07      = 'dp_2D_AaAc_07      '
    CHARACTER(LEN=256) :: name_var_dp_2D_AaAc_08      = 'dp_2D_AaAc_08      '
    CHARACTER(LEN=256) :: name_var_dp_2D_AaAc_09      = 'dp_2D_AaAc_09      '
    CHARACTER(LEN=256) :: name_var_dp_2D_AaAc_10      = 'dp_2D_AaAc_10      '
    
    INTEGER :: id_var_dp_3D_Aa_01
    INTEGER :: id_var_dp_3D_Aa_02
    INTEGER :: id_var_dp_3D_Aa_03
    INTEGER :: id_var_dp_3D_Aa_04
    INTEGER :: id_var_dp_3D_Aa_05
    INTEGER :: id_var_dp_3D_Aa_06
    INTEGER :: id_var_dp_3D_Aa_07
    INTEGER :: id_var_dp_3D_Aa_08
    INTEGER :: id_var_dp_3D_Aa_09
    INTEGER :: id_var_dp_3D_Aa_10
    
    CHARACTER(LEN=256) :: name_var_dp_3D_Aa_01         = 'dp_3D_Aa_01      '
    CHARACTER(LEN=256) :: name_var_dp_3D_Aa_02         = 'dp_3D_Aa_02      '
    CHARACTER(LEN=256) :: name_var_dp_3D_Aa_03         = 'dp_3D_Aa_03      '
    CHARACTER(LEN=256) :: name_var_dp_3D_Aa_04         = 'dp_3D_Aa_04      '
    CHARACTER(LEN=256) :: name_var_dp_3D_Aa_05         = 'dp_3D_Aa_05      '
    CHARACTER(LEN=256) :: name_var_dp_3D_Aa_06         = 'dp_3D_Aa_06      '
    CHARACTER(LEN=256) :: name_var_dp_3D_Aa_07         = 'dp_3D_Aa_07      '
    CHARACTER(LEN=256) :: name_var_dp_3D_Aa_08         = 'dp_3D_Aa_08      '
    CHARACTER(LEN=256) :: name_var_dp_3D_Aa_09         = 'dp_3D_Aa_09      '
    CHARACTER(LEN=256) :: name_var_dp_3D_Aa_10         = 'dp_3D_Aa_10      '
    
    INTEGER :: id_var_dp_2D_monthly_Aa_01
    INTEGER :: id_var_dp_2D_monthly_Aa_02
    INTEGER :: id_var_dp_2D_monthly_Aa_03
    INTEGER :: id_var_dp_2D_monthly_Aa_04
    INTEGER :: id_var_dp_2D_monthly_Aa_05
    INTEGER :: id_var_dp_2D_monthly_Aa_06
    INTEGER :: id_var_dp_2D_monthly_Aa_07
    INTEGER :: id_var_dp_2D_monthly_Aa_08
    INTEGER :: id_var_dp_2D_monthly_Aa_09
    INTEGER :: id_var_dp_2D_monthly_Aa_10
    
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_Aa_01 = 'dp_2D_monthly_Aa_01'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_Aa_02 = 'dp_2D_monthly_Aa_02'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_Aa_03 = 'dp_2D_monthly_Aa_03'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_Aa_04 = 'dp_2D_monthly_Aa_04'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_Aa_05 = 'dp_2D_monthly_Aa_05'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_Aa_06 = 'dp_2D_monthly_Aa_06'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_Aa_07 = 'dp_2D_monthly_Aa_07'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_Aa_08 = 'dp_2D_monthly_Aa_08'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_Aa_09 = 'dp_2D_monthly_Aa_09'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_Aa_10 = 'dp_2D_monthly_Aa_10'
        
  END TYPE type_netcdf_debug
    
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
