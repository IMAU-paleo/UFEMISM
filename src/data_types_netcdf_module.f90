MODULE data_types_netcdf_module
  ! Contains the TYPES for different NetCDf files read and written by UFEMISM.

! ===== USE modules =====
! =======================

  USE configuration_module,        ONLY: dp, C

  IMPLICIT NONE

! ===== Data types =====
! ======================

  ! == Grid
  ! ==========

  TYPE type_netcdf_grid
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

  END TYPE type_netcdf_grid

  ! == Restart
  ! ==========

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
    INTEGER :: id_dim_six
    INTEGER :: id_dim_vii_transect

    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_aci                   = 'aci                  '
    CHARACTER(LEN=256) :: name_dim_ciplusone             = 'ciplusone            '
    CHARACTER(LEN=256) :: name_dim_two                   = 'two                  '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_six                   = 'six                  '
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
    INTEGER :: id_var_SL
    INTEGER :: id_var_dHb
    INTEGER :: id_var_beta_sq
    INTEGER :: id_var_phi_fric
    INTEGER :: id_var_Ti
    INTEGER :: id_var_FirnDepth
    INTEGER :: id_var_MeltPreviousYear
    INTEGER :: id_var_C_abl_constant_inv
    INTEGER :: id_var_C_abl_Ts_inv
    INTEGER :: id_var_C_abl_Q_inv
    INTEGER :: id_var_C_refr_inv

    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   '
    CHARACTER(LEN=256) :: name_var_SL                    = 'SL                   '
    CHARACTER(LEN=256) :: name_var_dHb                   = 'dHb                  '
    CHARACTER(LEN=256) :: name_var_beta_sq               = 'beta_sq              '
    CHARACTER(LEN=256) :: name_var_phi_fric              = 'phi_fric             '
    CHARACTER(LEN=256) :: name_var_Ti                    = 'Ti                   '
    CHARACTER(LEN=256) :: name_var_FirnDepth             = 'FirnDepth            '
    CHARACTER(LEN=256) :: name_var_MeltPreviousYear      = 'MeltPreviousYear     '
    CHARACTER(LEN=256) :: name_var_C_abl_constant_inv    = 'C_abl_constant_inv   '
    CHARACTER(LEN=256) :: name_var_C_abl_Ts_inv          = 'C_abl_Ts_inv         '
    CHARACTER(LEN=256) :: name_var_C_abl_Q_inv           = 'C_abl_Q_inv          '
    CHARACTER(LEN=256) :: name_var_C_refr_inv            = 'C_refr_inv           '

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
    INTEGER :: id_dim_six
    INTEGER :: id_dim_vii_transect

    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_aci                   = 'aci                  '
    CHARACTER(LEN=256) :: name_dim_ciplusone             = 'ciplusone            '
    CHARACTER(LEN=256) :: name_dim_two                   = 'two                  '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_six                   = 'six                  '
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

  ! == Debugging/Profiling
  ! ======================

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
    INTEGER :: id_dim_six
    INTEGER :: id_dim_vii_transect

    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_aci                   = 'aci                  '
    CHARACTER(LEN=256) :: name_dim_ciplusone             = 'ciplusone            '
    CHARACTER(LEN=256) :: name_dim_two                   = 'two                  '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_six                   = 'six                  '
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

    INTEGER :: id_var_int_2D_a_01
    INTEGER :: id_var_int_2D_a_02
    INTEGER :: id_var_int_2D_a_03
    INTEGER :: id_var_int_2D_a_04
    INTEGER :: id_var_int_2D_a_05
    INTEGER :: id_var_int_2D_a_06
    INTEGER :: id_var_int_2D_a_07
    INTEGER :: id_var_int_2D_a_08
    INTEGER :: id_var_int_2D_a_09
    INTEGER :: id_var_int_2D_a_10

    CHARACTER(LEN=256) :: name_var_int_2D_a_01        = 'int_2D_a_01      '
    CHARACTER(LEN=256) :: name_var_int_2D_a_02        = 'int_2D_a_02      '
    CHARACTER(LEN=256) :: name_var_int_2D_a_03        = 'int_2D_a_03      '
    CHARACTER(LEN=256) :: name_var_int_2D_a_04        = 'int_2D_a_04      '
    CHARACTER(LEN=256) :: name_var_int_2D_a_05        = 'int_2D_a_05      '
    CHARACTER(LEN=256) :: name_var_int_2D_a_06        = 'int_2D_a_06      '
    CHARACTER(LEN=256) :: name_var_int_2D_a_07        = 'int_2D_a_07      '
    CHARACTER(LEN=256) :: name_var_int_2D_a_08        = 'int_2D_a_08      '
    CHARACTER(LEN=256) :: name_var_int_2D_a_09        = 'int_2D_a_09      '
    CHARACTER(LEN=256) :: name_var_int_2D_a_10        = 'int_2D_a_10      '

    INTEGER :: id_var_int_2D_b_01
    INTEGER :: id_var_int_2D_b_02
    INTEGER :: id_var_int_2D_b_03
    INTEGER :: id_var_int_2D_b_04
    INTEGER :: id_var_int_2D_b_05
    INTEGER :: id_var_int_2D_b_06
    INTEGER :: id_var_int_2D_b_07
    INTEGER :: id_var_int_2D_b_08
    INTEGER :: id_var_int_2D_b_09
    INTEGER :: id_var_int_2D_b_10

    CHARACTER(LEN=256) :: name_var_int_2D_b_01        = 'int_2D_b_01      '
    CHARACTER(LEN=256) :: name_var_int_2D_b_02        = 'int_2D_b_02      '
    CHARACTER(LEN=256) :: name_var_int_2D_b_03        = 'int_2D_b_03      '
    CHARACTER(LEN=256) :: name_var_int_2D_b_04        = 'int_2D_b_04      '
    CHARACTER(LEN=256) :: name_var_int_2D_b_05        = 'int_2D_b_05      '
    CHARACTER(LEN=256) :: name_var_int_2D_b_06        = 'int_2D_b_06      '
    CHARACTER(LEN=256) :: name_var_int_2D_b_07        = 'int_2D_b_07      '
    CHARACTER(LEN=256) :: name_var_int_2D_b_08        = 'int_2D_b_08      '
    CHARACTER(LEN=256) :: name_var_int_2D_b_09        = 'int_2D_b_09      '
    CHARACTER(LEN=256) :: name_var_int_2D_b_10        = 'int_2D_b_10      '

    INTEGER :: id_var_int_2D_c_01
    INTEGER :: id_var_int_2D_c_02
    INTEGER :: id_var_int_2D_c_03
    INTEGER :: id_var_int_2D_c_04
    INTEGER :: id_var_int_2D_c_05
    INTEGER :: id_var_int_2D_c_06
    INTEGER :: id_var_int_2D_c_07
    INTEGER :: id_var_int_2D_c_08
    INTEGER :: id_var_int_2D_c_09
    INTEGER :: id_var_int_2D_c_10

    CHARACTER(LEN=256) :: name_var_int_2D_c_01        = 'int_2D_c_01      '
    CHARACTER(LEN=256) :: name_var_int_2D_c_02        = 'int_2D_c_02      '
    CHARACTER(LEN=256) :: name_var_int_2D_c_03        = 'int_2D_c_03      '
    CHARACTER(LEN=256) :: name_var_int_2D_c_04        = 'int_2D_c_04      '
    CHARACTER(LEN=256) :: name_var_int_2D_c_05        = 'int_2D_c_05      '
    CHARACTER(LEN=256) :: name_var_int_2D_c_06        = 'int_2D_c_06      '
    CHARACTER(LEN=256) :: name_var_int_2D_c_07        = 'int_2D_c_07      '
    CHARACTER(LEN=256) :: name_var_int_2D_c_08        = 'int_2D_c_08      '
    CHARACTER(LEN=256) :: name_var_int_2D_c_09        = 'int_2D_c_09      '
    CHARACTER(LEN=256) :: name_var_int_2D_c_10        = 'int_2D_c_10      '

    INTEGER :: id_var_int_2D_ac_01
    INTEGER :: id_var_int_2D_ac_02
    INTEGER :: id_var_int_2D_ac_03
    INTEGER :: id_var_int_2D_ac_04
    INTEGER :: id_var_int_2D_ac_05
    INTEGER :: id_var_int_2D_ac_06
    INTEGER :: id_var_int_2D_ac_07
    INTEGER :: id_var_int_2D_ac_08
    INTEGER :: id_var_int_2D_ac_09
    INTEGER :: id_var_int_2D_ac_10

    CHARACTER(LEN=256) :: name_var_int_2D_ac_01        = 'int_2D_ac_01      '
    CHARACTER(LEN=256) :: name_var_int_2D_ac_02        = 'int_2D_ac_02      '
    CHARACTER(LEN=256) :: name_var_int_2D_ac_03        = 'int_2D_ac_03      '
    CHARACTER(LEN=256) :: name_var_int_2D_ac_04        = 'int_2D_ac_04      '
    CHARACTER(LEN=256) :: name_var_int_2D_ac_05        = 'int_2D_ac_05      '
    CHARACTER(LEN=256) :: name_var_int_2D_ac_06        = 'int_2D_ac_06      '
    CHARACTER(LEN=256) :: name_var_int_2D_ac_07        = 'int_2D_ac_07      '
    CHARACTER(LEN=256) :: name_var_int_2D_ac_08        = 'int_2D_ac_08      '
    CHARACTER(LEN=256) :: name_var_int_2D_ac_09        = 'int_2D_ac_09      '
    CHARACTER(LEN=256) :: name_var_int_2D_ac_10        = 'int_2D_ac_10      '

    INTEGER :: id_var_dp_2D_a_01
    INTEGER :: id_var_dp_2D_a_02
    INTEGER :: id_var_dp_2D_a_03
    INTEGER :: id_var_dp_2D_a_04
    INTEGER :: id_var_dp_2D_a_05
    INTEGER :: id_var_dp_2D_a_06
    INTEGER :: id_var_dp_2D_a_07
    INTEGER :: id_var_dp_2D_a_08
    INTEGER :: id_var_dp_2D_a_09
    INTEGER :: id_var_dp_2D_a_10

    CHARACTER(LEN=256) :: name_var_dp_2D_a_01        = 'dp_2D_a_01      '
    CHARACTER(LEN=256) :: name_var_dp_2D_a_02        = 'dp_2D_a_02      '
    CHARACTER(LEN=256) :: name_var_dp_2D_a_03        = 'dp_2D_a_03      '
    CHARACTER(LEN=256) :: name_var_dp_2D_a_04        = 'dp_2D_a_04      '
    CHARACTER(LEN=256) :: name_var_dp_2D_a_05        = 'dp_2D_a_05      '
    CHARACTER(LEN=256) :: name_var_dp_2D_a_06        = 'dp_2D_a_06      '
    CHARACTER(LEN=256) :: name_var_dp_2D_a_07        = 'dp_2D_a_07      '
    CHARACTER(LEN=256) :: name_var_dp_2D_a_08        = 'dp_2D_a_08      '
    CHARACTER(LEN=256) :: name_var_dp_2D_a_09        = 'dp_2D_a_09      '
    CHARACTER(LEN=256) :: name_var_dp_2D_a_10        = 'dp_2D_a_10      '

    INTEGER :: id_var_dp_2D_b_01
    INTEGER :: id_var_dp_2D_b_02
    INTEGER :: id_var_dp_2D_b_03
    INTEGER :: id_var_dp_2D_b_04
    INTEGER :: id_var_dp_2D_b_05
    INTEGER :: id_var_dp_2D_b_06
    INTEGER :: id_var_dp_2D_b_07
    INTEGER :: id_var_dp_2D_b_08
    INTEGER :: id_var_dp_2D_b_09
    INTEGER :: id_var_dp_2D_b_10

    CHARACTER(LEN=256) :: name_var_dp_2D_b_01        = 'dp_2D_b_01      '
    CHARACTER(LEN=256) :: name_var_dp_2D_b_02        = 'dp_2D_b_02      '
    CHARACTER(LEN=256) :: name_var_dp_2D_b_03        = 'dp_2D_b_03      '
    CHARACTER(LEN=256) :: name_var_dp_2D_b_04        = 'dp_2D_b_04      '
    CHARACTER(LEN=256) :: name_var_dp_2D_b_05        = 'dp_2D_b_05      '
    CHARACTER(LEN=256) :: name_var_dp_2D_b_06        = 'dp_2D_b_06      '
    CHARACTER(LEN=256) :: name_var_dp_2D_b_07        = 'dp_2D_b_07      '
    CHARACTER(LEN=256) :: name_var_dp_2D_b_08        = 'dp_2D_b_08      '
    CHARACTER(LEN=256) :: name_var_dp_2D_b_09        = 'dp_2D_b_09      '
    CHARACTER(LEN=256) :: name_var_dp_2D_b_10        = 'dp_2D_b_10      '

    INTEGER :: id_var_dp_2D_c_01
    INTEGER :: id_var_dp_2D_c_02
    INTEGER :: id_var_dp_2D_c_03
    INTEGER :: id_var_dp_2D_c_04
    INTEGER :: id_var_dp_2D_c_05
    INTEGER :: id_var_dp_2D_c_06
    INTEGER :: id_var_dp_2D_c_07
    INTEGER :: id_var_dp_2D_c_08
    INTEGER :: id_var_dp_2D_c_09
    INTEGER :: id_var_dp_2D_c_10

    CHARACTER(LEN=256) :: name_var_dp_2D_c_01        = 'dp_2D_c_01      '
    CHARACTER(LEN=256) :: name_var_dp_2D_c_02        = 'dp_2D_c_02      '
    CHARACTER(LEN=256) :: name_var_dp_2D_c_03        = 'dp_2D_c_03      '
    CHARACTER(LEN=256) :: name_var_dp_2D_c_04        = 'dp_2D_c_04      '
    CHARACTER(LEN=256) :: name_var_dp_2D_c_05        = 'dp_2D_c_05      '
    CHARACTER(LEN=256) :: name_var_dp_2D_c_06        = 'dp_2D_c_06      '
    CHARACTER(LEN=256) :: name_var_dp_2D_c_07        = 'dp_2D_c_07      '
    CHARACTER(LEN=256) :: name_var_dp_2D_c_08        = 'dp_2D_c_08      '
    CHARACTER(LEN=256) :: name_var_dp_2D_c_09        = 'dp_2D_c_09      '
    CHARACTER(LEN=256) :: name_var_dp_2D_c_10        = 'dp_2D_c_10      '

    INTEGER :: id_var_dp_2D_ac_01
    INTEGER :: id_var_dp_2D_ac_02
    INTEGER :: id_var_dp_2D_ac_03
    INTEGER :: id_var_dp_2D_ac_04
    INTEGER :: id_var_dp_2D_ac_05
    INTEGER :: id_var_dp_2D_ac_06
    INTEGER :: id_var_dp_2D_ac_07
    INTEGER :: id_var_dp_2D_ac_08
    INTEGER :: id_var_dp_2D_ac_09
    INTEGER :: id_var_dp_2D_ac_10

    CHARACTER(LEN=256) :: name_var_dp_2D_ac_01        = 'dp_2D_ac_01      '
    CHARACTER(LEN=256) :: name_var_dp_2D_ac_02        = 'dp_2D_ac_02      '
    CHARACTER(LEN=256) :: name_var_dp_2D_ac_03        = 'dp_2D_ac_03      '
    CHARACTER(LEN=256) :: name_var_dp_2D_ac_04        = 'dp_2D_ac_04      '
    CHARACTER(LEN=256) :: name_var_dp_2D_ac_05        = 'dp_2D_ac_05      '
    CHARACTER(LEN=256) :: name_var_dp_2D_ac_06        = 'dp_2D_ac_06      '
    CHARACTER(LEN=256) :: name_var_dp_2D_ac_07        = 'dp_2D_ac_07      '
    CHARACTER(LEN=256) :: name_var_dp_2D_ac_08        = 'dp_2D_ac_08      '
    CHARACTER(LEN=256) :: name_var_dp_2D_ac_09        = 'dp_2D_ac_09      '
    CHARACTER(LEN=256) :: name_var_dp_2D_ac_10        = 'dp_2D_ac_10      '

    INTEGER :: id_var_dp_3D_a_01
    INTEGER :: id_var_dp_3D_a_02
    INTEGER :: id_var_dp_3D_a_03
    INTEGER :: id_var_dp_3D_a_04
    INTEGER :: id_var_dp_3D_a_05
    INTEGER :: id_var_dp_3D_a_06
    INTEGER :: id_var_dp_3D_a_07
    INTEGER :: id_var_dp_3D_a_08
    INTEGER :: id_var_dp_3D_a_09
    INTEGER :: id_var_dp_3D_a_10

    CHARACTER(LEN=256) :: name_var_dp_3D_a_01        = 'dp_3D_a_01      '
    CHARACTER(LEN=256) :: name_var_dp_3D_a_02        = 'dp_3D_a_02      '
    CHARACTER(LEN=256) :: name_var_dp_3D_a_03        = 'dp_3D_a_03      '
    CHARACTER(LEN=256) :: name_var_dp_3D_a_04        = 'dp_3D_a_04      '
    CHARACTER(LEN=256) :: name_var_dp_3D_a_05        = 'dp_3D_a_05      '
    CHARACTER(LEN=256) :: name_var_dp_3D_a_06        = 'dp_3D_a_06      '
    CHARACTER(LEN=256) :: name_var_dp_3D_a_07        = 'dp_3D_a_07      '
    CHARACTER(LEN=256) :: name_var_dp_3D_a_08        = 'dp_3D_a_08      '
    CHARACTER(LEN=256) :: name_var_dp_3D_a_09        = 'dp_3D_a_09      '
    CHARACTER(LEN=256) :: name_var_dp_3D_a_10        = 'dp_3D_a_10      '

    INTEGER :: id_var_dp_2D_monthly_a_01
    INTEGER :: id_var_dp_2D_monthly_a_02
    INTEGER :: id_var_dp_2D_monthly_a_03
    INTEGER :: id_var_dp_2D_monthly_a_04
    INTEGER :: id_var_dp_2D_monthly_a_05
    INTEGER :: id_var_dp_2D_monthly_a_06
    INTEGER :: id_var_dp_2D_monthly_a_07
    INTEGER :: id_var_dp_2D_monthly_a_08
    INTEGER :: id_var_dp_2D_monthly_a_09
    INTEGER :: id_var_dp_2D_monthly_a_10

    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_a_01        = 'dp_2D_monthly_a_01      '
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_a_02        = 'dp_2D_monthly_a_02      '
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_a_03        = 'dp_2D_monthly_a_03      '
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_a_04        = 'dp_2D_monthly_a_04      '
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_a_05        = 'dp_2D_monthly_a_05      '
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_a_06        = 'dp_2D_monthly_a_06      '
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_a_07        = 'dp_2D_monthly_a_07      '
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_a_08        = 'dp_2D_monthly_a_08      '
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_a_09        = 'dp_2D_monthly_a_09      '
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_a_10        = 'dp_2D_monthly_a_10      '

  END TYPE type_netcdf_debug

  TYPE type_netcdf_resource_tracker
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! Index of time frame to be written to
    INTEGER :: ti

    ! Dimensions
    ! ==========

    INTEGER :: id_dim_time
    INTEGER :: id_dim_name_length

    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_name_length           = 'name_length          '

    INTEGER :: id_var_time

    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '

    ! Variables
    ! =========

    ! Total model resource use
    INTEGER :: id_var_tcomp_tot
    INTEGER :: id_var_mem_tot

    ! Per-subroutine resource use
    INTEGER, DIMENSION(:), ALLOCATABLE :: id_var_names
    INTEGER, DIMENSION(:), ALLOCATABLE :: id_var_tcomp
    INTEGER, DIMENSION(:), ALLOCATABLE :: id_var_mem

  END TYPE type_netcdf_resource_tracker

  ! == Reference geometries
  ! =======================

  TYPE type_netcdf_reference_geometry
    ! For reading an input file describing a reference ice-sheet geometry on a Cartesian grid

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! ID's for variables:
    ! ===================

    ! Dimensions
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

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

  END TYPE type_netcdf_reference_geometry

  ! == ???
  ! ======

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

  TYPE type_netcdf_ICE5G_data
    ! For reading one of the ICE5G time frames

    ! Integers describing open ports to different variables in an opened NetCDF file.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! ID's for variables:
    ! ===================

    ! Dimensions
    INTEGER :: id_dim_lon
    INTEGER :: id_dim_lat

    CHARACTER(LEN=256) :: name_dim_lon                   = 'long                 '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '

    ! Variables
    INTEGER :: id_var_lon
    INTEGER :: id_var_lat
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hb
    INTEGER :: id_var_mask_ice

    CHARACTER(LEN=256) :: name_var_lon                   = 'long                 '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_Hi                    = 'sftgit               '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'orog                 '
    CHARACTER(LEN=256) :: name_var_mask_ice              = 'sftgif               '

  END TYPE type_netcdf_ICE5G_data

  ! == Insolation
  ! =============

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

  ! == Geothermal heat flux
  ! =======================

  TYPE type_netcdf_geothermal_heat_flux
    ! For reading an input file containing geothermal heat flux (e.g. Shapiro and Ritzwoller, 2004),
    ! describing geothermal heat flux at a lon-lat grid.

    ! Integers describing open ports to different variables in an opened NetCDF file.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! ID's for variables:
    ! ===================

    ! Dimensions
    INTEGER :: id_dim_lon
    INTEGER :: id_dim_lat

    CHARACTER(LEN=256) :: name_dim_lon                   = 'Longitude            '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'Latitude             '

    ! Variables
    INTEGER :: id_var_lon
    INTEGER :: id_var_lat
    INTEGER :: id_var_ghf

    CHARACTER(LEN=256) :: name_var_lon                   = 'Longitude            '
    CHARACTER(LEN=256) :: name_var_lat                   = 'Latitude             '
    CHARACTER(LEN=256) :: name_var_ghf                   = 'hflux                '

  END TYPE type_netcdf_geothermal_heat_flux

  ! == Climate
  ! ==========

  TYPE type_netcdf_direct_climate_forcing_global
    ! For reading an input file containing climate data,
    ! describing 2-m air temperature and precipitation, on a global lon/lat-grid.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! ID's for variables:
    ! ===================

    ! Dimensions
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
    INTEGER :: id_dim_lon
    INTEGER :: id_dim_lat

    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    CHARACTER(LEN=256) :: name_dim_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '

    ! Variables
    INTEGER :: id_var_time
    INTEGER :: id_var_month
    INTEGER :: id_var_lon
    INTEGER :: id_var_lat

    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '

    INTEGER :: id_var_T2m
    INTEGER :: id_var_Precip

    CHARACTER(LEN=256) :: name_var_T2m                   = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_Precip                = 'Precip               '

  END TYPE type_netcdf_direct_climate_forcing_global

  TYPE type_netcdf_direct_climate_forcing_regional
    ! For reading an input file containing climate data,
    ! describing 2-m air temperature and precipitation, on a regional x/y-grid.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! ID's for variables:
    ! ===================

    ! Dimensions
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    CHARACTER(LEN=256) :: name_dim_x                     = 'NX                   '
    CHARACTER(LEN=256) :: name_dim_y                     = 'NY                   '

    ! Variables
    INTEGER :: id_var_time
    INTEGER :: id_var_month
    INTEGER :: id_var_x
    INTEGER :: id_var_y

    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

    INTEGER :: id_var_T2m
    INTEGER :: id_var_Precip

    CHARACTER(LEN=256) :: name_var_T2m                   = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_Precip                = 'Precip               '

  END TYPE type_netcdf_direct_climate_forcing_regional

  TYPE type_netcdf_direct_SMB_forcing_global
    ! For reading an input file containing SMB, on a global lon/lat-grid.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! ID's for variables:
    ! ===================

    ! Dimensions
    INTEGER :: id_dim_time
    INTEGER :: id_dim_lon
    INTEGER :: id_dim_lat

    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '

    ! Variables
    INTEGER :: id_var_time
    INTEGER :: id_var_lon
    INTEGER :: id_var_lat

    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '

    INTEGER :: id_var_T2m_year
    INTEGER :: id_var_SMB_year

    CHARACTER(LEN=256) :: name_var_T2m_year              = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_SMB_year              = 'SMB                  '

  END TYPE type_netcdf_direct_SMB_forcing_global

  TYPE type_netcdf_direct_SMB_forcing_regional
    ! For reading an input file containing SMB, on a regional x/y-grid.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! ID's for variables:
    ! ===================

    ! Dimensions
    INTEGER :: id_dim_time
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_x                     = 'NX                   '
    CHARACTER(LEN=256) :: name_dim_y                     = 'NY                   '

    ! Variables
    INTEGER :: id_var_time
    INTEGER :: id_var_x
    INTEGER :: id_var_y

    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

    INTEGER :: id_var_T2m_year
    INTEGER :: id_var_SMB_year

    CHARACTER(LEN=256) :: name_var_T2m_year              = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_SMB_year              = 'SMB                  '

  END TYPE type_netcdf_direct_SMB_forcing_regional

  TYPE type_netcdf_ocean_data
    ! For reading an input file containing either a GCM ocean snapshot or a PD observations data set (e.g. WOA18),
    ! describing the global ocean with yearly fields on a lat/lon/depth grid

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! ID's for variables:
    ! ===================

    ! Dimensions
    INTEGER :: id_dim_lon
    INTEGER :: id_dim_lat
    INTEGER :: id_dim_z_ocean

    CHARACTER(LEN=256) :: name_dim_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_dim_z_ocean               = 'depth                '

    ! Variables
    INTEGER :: id_var_lon
    INTEGER :: id_var_lat
    INTEGER :: id_var_z_ocean
    INTEGER :: id_var_T_ocean
    INTEGER :: id_var_S_ocean

    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_z_ocean               = 'depth                '

  END TYPE type_netcdf_ocean_data

  TYPE type_netcdf_extrapolated_ocean_data
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! Index of time frame to be written to
    INTEGER :: ti

    ! Dimensions
    ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_z_ocean

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_z_ocean               = 'z_ocean              '

    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_z_ocean

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_z_ocean               = 'z_ocean              '

    ! Variables
    ! =========

    INTEGER :: id_var_T_ocean
    INTEGER :: id_var_S_ocean

    CHARACTER(LEN=256) :: name_var_T_ocean               = 'T_ocean              '
    CHARACTER(LEN=256) :: name_var_S_ocean               = 'S_ocean              '

  END TYPE type_netcdf_extrapolated_ocean_data

  ! == SELEN
  ! ========

  TYPE type_netcdf_SELEN_global_topo
    ! A NETCDF file containing global topography data for SELEN on an irregular global mesh

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! Dimensions
    ! ==========

    INTEGER :: id_dim_vi
    INTEGER :: id_dim_ti
    INTEGER :: id_dim_ci
    INTEGER :: id_dim_three

    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '

    INTEGER :: id_var_V
    INTEGER :: id_var_Tri
    INTEGER :: id_var_nC
    INTEGER :: id_var_C
    INTEGER :: id_var_niTri
    INTEGER :: id_var_iTri

    CHARACTER(LEN=256) :: name_var_V                     = 'V                    '
    CHARACTER(LEN=256) :: name_var_Tri                   = 'Tri                  '
    CHARACTER(LEN=256) :: name_var_nC                    = 'nC                   '
    CHARACTER(LEN=256) :: name_var_C                     = 'C                    '
    CHARACTER(LEN=256) :: name_var_niTri                 = 'niTri                '
    CHARACTER(LEN=256) :: name_var_iTri                  = 'iTri                 '

    ! Variables
    ! =========

    INTEGER :: id_var_lat
    INTEGER :: id_var_lon
    INTEGER :: id_var_Hb
    INTEGER :: id_var_ianc

    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_ianc                  = 'ianc                 '

  END TYPE type_netcdf_SELEN_global_topo

  TYPE type_netcdf_SELEN_output
    ! A NETCDF file containing output of SELEN (ice loading, bed topography and geoid perturbation)
    ! on the global SELEN mesh

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! Index of time frame to be written to
    INTEGER :: ti

    ! Dimensions
    ! ==========

    INTEGER :: id_dim_vi
    INTEGER :: id_dim_ti
    INTEGER :: id_dim_ci
    INTEGER :: id_dim_three
    INTEGER :: id_dim_time
    INTEGER :: id_dim_ki

    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_ki                    = 'ki                   '

    INTEGER :: id_var_V
    INTEGER :: id_var_Tri
    INTEGER :: id_var_nC
    INTEGER :: id_var_C
    INTEGER :: id_var_niTri
    INTEGER :: id_var_iTri
    INTEGER :: id_var_time
    INTEGER :: id_var_ki

    CHARACTER(LEN=256) :: name_var_V                     = 'V                    '
    CHARACTER(LEN=256) :: name_var_Tri                   = 'Tri                  '
    CHARACTER(LEN=256) :: name_var_nC                    = 'nC                   '
    CHARACTER(LEN=256) :: name_var_C                     = 'C                    '
    CHARACTER(LEN=256) :: name_var_niTri                 = 'niTri                '
    CHARACTER(LEN=256) :: name_var_iTri                  = 'iTri                 '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_ki                    = 'ki                   '

    ! Variables
    ! =========

    INTEGER :: id_var_lat
    INTEGER :: id_var_lon
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hi_rel
    INTEGER :: id_var_N
    INTEGER :: id_var_U
    INTEGER :: id_var_ocean_function
    INTEGER :: id_var_load_history

    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hi_rel                = 'Hi_rel               '
    CHARACTER(LEN=256) :: name_var_N                     = 'N                    '
    CHARACTER(LEN=256) :: name_var_U                     = 'U                    '
    CHARACTER(LEN=256) :: name_var_ocean_function        = 'ocean_function       '

    CHARACTER(LEN=256) :: name_var_load_history          = 'load_history         '

  END TYPE type_netcdf_SELEN_output

  TYPE type_netcdf_BIV_bed_roughness
    ! A NetCDF file containing bed roughness resulting from a basal inversion routine

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    ! Grid
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

    ! Mesh
    INTEGER :: id_dim_vi
    INTEGER :: id_dim_ti
    INTEGER :: id_dim_ci
    INTEGER :: id_dim_aci
    INTEGER :: id_dim_ai
    INTEGER :: id_dim_tai
    INTEGER :: id_dim_ciplusone
    INTEGER :: id_dim_two
    INTEGER :: id_dim_three
    INTEGER :: id_dim_six
    INTEGER :: id_dim_vii_transect

    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_aci                   = 'aci                  '
    CHARACTER(LEN=256) :: name_dim_ciplusone             = 'ciplusone            '
    CHARACTER(LEN=256) :: name_dim_two                   = 'two                  '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_six                   = 'six                  '
    CHARACTER(LEN=256) :: name_dim_vii_transect          = 'vii                  '
    CHARACTER(LEN=256) :: name_dim_ai                    = 'ai                   '
    CHARACTER(LEN=256) :: name_dim_tai                   = 'tai                  '

  ! Variables
  ! =========

    ! Grid
    INTEGER :: id_var_x
    INTEGER :: id_var_y

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

    ! Mesh
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

    ! Field variables
    INTEGER :: id_var_phi_fric
    INTEGER :: id_var_alpha_sq
    INTEGER :: id_var_beta_sq

    CHARACTER(LEN=256) :: name_var_phi_fric              = 'phi_fric             '
    CHARACTER(LEN=256) :: name_var_alpha_sq              = 'alpha_sq             '
    CHARACTER(LEN=256) :: name_var_beta_sq               = 'beta_sq              '

  END TYPE type_netcdf_BIV_bed_roughness

  TYPE type_netcdf_BIV_target_velocity
    ! A NetCDF file containing surface [u,v]-velocity fields to be used as
    ! the target in a velocity-based basal inversion routine

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

    ! Field variables
    INTEGER :: id_var_u_surf
    INTEGER :: id_var_v_surf

    CHARACTER(LEN=256) :: name_var_u_surf                = 'u_surf               '
    CHARACTER(LEN=256) :: name_var_v_surf                = 'v_surf               '

  END TYPE type_netcdf_BIV_target_velocity

  TYPE type_netcdf_ISMIP_style_forcing
    ! NetCDF files containing the different ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing fields

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_time

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_time

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '

    ! Field variables
    INTEGER :: id_var_aSMB
    INTEGER :: id_var_dSMBdz
    INTEGER :: id_var_aST
    INTEGER :: id_var_dSTdz

    CHARACTER(LEN=256) :: name_var_aSMB                  = 'aSMB                 '
    CHARACTER(LEN=256) :: name_var_dSMBdz                = 'dSMBdz               '
    CHARACTER(LEN=256) :: name_var_aST                   = 'aST                  '
    CHARACTER(LEN=256) :: name_var_dSTdz                 = 'dSTdz                '

  END TYPE type_netcdf_ISMIP_style_forcing

  TYPE type_netcdf_ISMIP_style_baseline
    ! NetCDF file containing the baseline climate and orography for the ISMIP-style climate forcing

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

    ! Field variables
    INTEGER :: id_var_Hs
    INTEGER :: id_var_SMB
    INTEGER :: id_var_ST

    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   '
    CHARACTER(LEN=256) :: name_var_SMB                   = 'SMB                  '
    CHARACTER(LEN=256) :: name_var_ST                    = 'ST                   '

  END TYPE type_netcdf_ISMIP_style_baseline

  TYPE type_netcdf_prescribed_retreat_mask
    ! NetCDF files containing a prescribed retreat mask

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_time

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_time

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '

    CHARACTER(LEN=256) :: time_units  ! ISMIP uses "days since XXX", paleo stuff just uses "years"

    ! Field variables
    INTEGER :: id_var_ice_fraction

  END TYPE type_netcdf_prescribed_retreat_mask

  TYPE type_netcdf_prescribed_retreat_mask_refice
    ! NetCDF files containing the reference ice thickness for a prescribed retreat mask

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

    ! Field variables
    INTEGER :: id_var_Hi

    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '

  END TYPE type_netcdf_prescribed_retreat_mask_refice

CONTAINS

END MODULE data_types_netcdf_module