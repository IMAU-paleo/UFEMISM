MODULE data_types_netcdf_module_new
  ! Contains the TYPES for different NetCDf files read and written by UFEMISM.

  USE configuration_module,        ONLY: dp, C

  IMPLICIT NONE

  ! ===== Grids and mesh =====
  ! ==========================

  TYPE type_netcdf_xy_grid
    ! Default names for dimensions and variables

    ! NetCDF file identifier
    INTEGER            :: ncid

  ! Dimensions
  ! ==========

    ! Names
    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

    ! Identifiers
    INTEGER            :: id_dim_x
    INTEGER            :: id_dim_y

  ! Variables
  ! =========

    ! Names
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

    ! Identifiers
    INTEGER            :: id_var_x
    INTEGER            :: id_var_y

  END TYPE type_netcdf_xy_grid

  TYPE type_netcdf_lonlat_grid
    ! Default names for dimensions and variables

    ! NetCDF file identifier
    INTEGER            :: ncid

  ! Dimensions
  ! ==========

    ! Names
    CHARACTER(LEN=256) :: name_dim_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '

    ! Identifiers
    INTEGER            :: id_dim_lon
    INTEGER            :: id_dim_lat

  ! Variables
  ! =========

    ! Names
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '

    ! Identifiers
    INTEGER            :: id_var_lon
    INTEGER            :: id_var_lat

  END TYPE type_netcdf_lonlat_grid

  TYPE type_netcdf_mesh
    ! Default names for dimensions and variables

    ! NetCDF file identifier
    INTEGER            :: ncid

  ! Dimensions
  ! ==========

    ! Names
    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_aci                   = 'aci                  '
    CHARACTER(LEN=256) :: name_dim_two                   = 'two                  '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_six                   = 'six                  '

    ! Identifiers
    INTEGER :: id_dim_vi
    INTEGER :: id_dim_ti
    INTEGER :: id_dim_ci
    INTEGER :: id_dim_aci
    INTEGER :: id_dim_two
    INTEGER :: id_dim_three
    INTEGER :: id_dim_six

  ! Variables
  ! =========

    ! Names
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
    CHARACTER(LEN=256) :: name_var_A                     = 'A                    '
    CHARACTER(LEN=256) :: name_var_R                     = 'R                    '

    ! Identifiers
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
    INTEGER :: id_var_A
    INTEGER :: id_var_R

  END TYPE type_netcdf_mesh

CONTAINS

END MODULE data_types_netcdf_module_new