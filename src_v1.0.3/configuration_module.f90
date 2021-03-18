MODULE configuration_module
  
  ! The way it's done right now:
  ! Each config variable has two versions: one with the "_config" extension, which is
  ! an actual variable in this module only, and one without the extension, which is
  ! a field in the "C" type. THe "_config" variables are used to create a NAMELIST,
  ! which makes reading an external config file really easy - anything in the file that
  ! matches a variable in the namelist overwrites the default value. After that's done,
  ! the fields in the "C" type are replaced with the values of the "_config" variables,
  ! which now have either the default values, or those specified in the external config
  ! file.
  !
  ! While this is certainly very convenient when running the model, it does make adding
  ! new config parameters a bit tedious - you have to add the "_config" variable, add it
  ! as a field in the "C" type, add it to the namelist, and let the "C" type field be
  ! overwritten in the end.
  !
  ! Some day I'll figure out a more elegant solution for this...
  
  USE mpi

  IMPLICIT NONE
  
  INTEGER, PARAMETER                :: dp  = KIND(1.0D0)  ! Kind of double precision numbers. Reals should be declared as: REAL(dp) :: example

  ! ===================================================================================
  ! "_config  variables, which will be collected into a NAMELIST, and possibly replaced
  ! by the values in the external config file. Remember the "_config" extension!
  ! ===================================================================================

  ! Time steps and range
  ! ====================
  
  REAL(dp) :: start_time_of_run_config   = 0.0_dp       ! Start time (in years) of the simulations
  REAL(dp) :: end_time_of_run_config     = 50000.0_dp   ! End   time (in years) of the simulations
  REAL(dp) :: dt_coupling_config         = 100._dp      ! Interval of coupling (in years) between the four ice-sheets  
  REAL(dp) :: dt_max_config              = 10.0_dp      ! Maximum time step (in years) of the ice model
  REAL(dp) :: dt_thermo_config           = 10.0_dp      ! Time step (in years) for updating thermodynamics
  REAL(dp) :: dt_climate_config          = 10._dp       ! Time step (in years) for updating the climate
  REAL(dp) :: dt_SMB_config              = 10._dp       ! Time step (in years) for updating the SMB
  REAL(dp) :: dt_BMB_config              = 10._dp       ! Time step (in years) for updating the BMB
  REAL(dp) :: dt_output_config           = 5000.0_dp    ! Time step (in years) for writing output
  REAL(dp) :: dt_mesh_min_config         = 50._dp       ! Minimum amount of time (in years) between mesh updates
  
  ! Which ice sheets do we simulate?
  ! ================================
  
  LOGICAL :: do_NAM_config               = .FALSE.      ! North America
  LOGICAL :: do_EAS_config               = .FALSE.      ! Eurasia
  LOGICAL :: do_GRL_config               = .FALSE.      ! Greenland
  LOGICAL :: do_ANT_config               = .TRUE.       ! Antarctica
  
  ! Benchmark experiments
  ! =====================
  
  LOGICAL            :: do_benchmark_experiment_config          = .TRUE.
  CHARACTER(LEN=256) :: choice_benchmark_experiment_config      = 'EISMINT_I'
  REAL(dp)           :: halfar_solution_H0_config               = 5000._dp
  REAL(dp)           :: halfar_solution_R0_config               = 300000._dp
  REAL(dp)           :: bueler_solution_lambda_config           = 5._dp

  ! Mesh generation parameters
  ! ==========================
  
  INTEGER  :: nconmax_config                 = 32           ! Maximum number of vertex connections
  REAL(dp) :: alpha_min_config               = 0.4363_dp    ! Minimum internal angle of triangles (25 degrees)
  REAL(dp) :: dz_max_ice_config              = 20000._dp    ! Maximum allowed 2nd order surface deviation over ice 
  REAL(dp) :: res_max_config                 = 800._dp      ! Maximum allowed resolution                            [km]
  REAL(dp) :: res_max_margin_config          = 40._dp       ! Maximum allowed resolution over land-based ice margin [km]
  REAL(dp) :: res_max_gl_config              = 40._dp       !                                 grounding line        [km]
  REAL(dp) :: res_max_cf_config              = 40._dp       !                                 calving front         [km]
  REAL(dp) :: res_max_mountain_config        = 40._dp       !                                 mountains             [km]
  REAL(dp) :: res_max_coast_config           = 40._dp       !                                 coastline             [km]
  REAL(dp) :: mesh_fitness_threshold_config  = 0.95_dp      ! Minimum allowed mesh fitness (fraction of triangles that are not Bad) before mesh updating
  
  ! Resolutions of the different square grids
  ! =========================================
  
  REAL(dp) :: dx_grid_output_config          = 40000._dp    ! Resolution of the square grid used for writing output                       [m]
  REAL(dp) :: dx_grid_GIA_config             = 100000._dp   ! Resolution of the square grid used for GIA modelling (ELRA or SELEN)        [m]
  REAL(dp) :: dx_grid_smooth_config          = 50000._dp    ! Resolution of the square grid used for data smoothing in the climate matrix [m]
  
  ! High-resolution Points Of Interest (POIs)
  ! =========================================
  
  INTEGER                    :: nPOI_NAM_config            = 0       ! Number of POIs
  INTEGER                    :: nPOI_EAS_config            = 0
  INTEGER                    :: nPOI_GRL_config            = 0
  INTEGER                    :: nPOI_ANT_config            = 0
  REAL(dp), DIMENSION(200)   :: POI_NAM_coordinates_config = 0._dp   ! [lat,lon] coordinates op POIs (degrees)
  REAL(dp), DIMENSION(200)   :: POI_EAS_coordinates_config = 0._dp
  REAL(dp), DIMENSION(200)   :: POI_GRL_coordinates_config = 0._dp
  REAL(dp), DIMENSION(200)   :: POI_ANT_coordinates_config = 0._dp
  REAL(dp), DIMENSION(100)   :: POI_NAM_resolutions_config = 0._dp   ! Required resolution at POIs (km)
  REAL(dp), DIMENSION(100)   :: POI_EAS_resolutions_config = 0._dp
  REAL(dp), DIMENSION(100)   :: POI_GRL_resolutions_config = 0._dp
  REAL(dp), DIMENSION(100)   :: POI_ANT_resolutions_config = 0._dp
  
  ! Whether or not to let UFEMISM dynamically create its own output folder.
  ! This works fine locally, on LISA its better to use a fixed folder name.
  ! =======================================================================
  
  LOGICAL            :: create_new_output_dir_config = .TRUE.
  CHARACTER(LEN=256) :: output_dir_config            = 'results_UFEMISM'
  LOGICAL            :: do_write_debug_data_config   = .FALSE.

  ! The scaled vertical coordinate zeta, used mainly in thermodynamics
  ! ==================================================================
  
  INTEGER                        :: nz_config   = 15
  REAL(dp), DIMENSION(210), SAVE :: zeta_config = &
   (/0.00_dp, 0.10_dp, 0.20_dp, 0.30_dp, 0.40_dp, 0.50_dp, 0.60_dp, 0.70_dp, 0.80_dp, 0.90_dp, 0.925_dp, 0.95_dp, 0.975_dp, 0.99_dp, 1.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp  /)

  ! Input data file paths
  ! =====================
  
  ! Initial model state
  CHARACTER(LEN=256) :: filename_init_NAM_config = '/Users/berends/Documents/Models/UFEMISM/data/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_init_EAS_config = '/Users/berends/Documents/Models/UFEMISM/data/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_init_GRL_config = '/Users/berends/Documents/Models/UFEMISM/data/Greenland_BedMachine_5km_noEllesmere.nc'
  CHARACTER(LEN=256) :: filename_init_ANT_config = '/Users/berends/Documents/Models/UFEMISM/data/Antarctica_Bedmap2_Rignot_5km.nc'
  
  ! Whether the initial model state is read from a Cartesian grid file (e.g. Bedmachine, ETOPO1, etc.) or
  ! a mesh file from a previous UFEMISM run.
  CHARACTER(LEN=4)   :: filetype_init_NAM_config = 'cart'   ! Can be either "cart" or "mesh"
  CHARACTER(LEN=4)   :: filetype_init_EAS_config = 'cart'
  CHARACTER(LEN=4)   :: filetype_init_GRL_config = 'cart'
  CHARACTER(LEN=4)   :: filetype_init_ANT_config = 'cart'
  
  ! PD reference data (NetCDF)
  CHARACTER(LEN=256) :: filename_PD_NAM_config   = '/Users/berends/Documents/Models/UFEMISM/data/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_PD_EAS_config   = '/Users/berends/Documents/Models/UFEMISM/data/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_PD_GRL_config   = '/Users/berends/Documents/Models/UFEMISM/data/Greenland_BedMachine_5km_noEllesmere.nc'
  CHARACTER(LEN=256) :: filename_PD_ANT_config   = '/Users/berends/Documents/Models/UFEMISM/data/Antarctica_Bedmap2_Rignot_5km.nc'
  
  ! Present-day observed climate (ERA40) (NetCDF)
  CHARACTER(LEN=256) :: filename_PD_obs_climate_config = '/Users/berends/Documents/Datasets/ERA40/ERA40_climate_global.nc'  
  ! Insolation forcing (Laskar et al., 2004) (TXT)
  CHARACTER(LEN=256) :: filename_insolation_config = '/Users/berends/Documents/Datasets/Insolation_laskar/Insolation_Laskar_etal_2004.nc'

  ! Ice dynamics and thermodynamics
  ! ===============================
  REAL(dp)           :: m_enh_sia_config                        = 5.0_dp                 ! Ice flow enhancement factor in the SIA
  REAL(dp)           :: m_enh_ssa_config                        = 0.7_dp                 ! Ice flow enhancement factor in the SSA
  CHARACTER(LEN=256) :: choice_sliding_law_config               = 'Coulomb_regularised'  ! Choice of sliding law (currently only "Coulomb_regularised" is implemented)
  REAL(dp)           :: C_sliding_config                        = 1.0E7_dp               ! Factor   in Weertman sliding law
  REAL(dp)           :: m_sliding_config                        = 1._dp/3._dp            ! Exponent in Weertman sliding law
  LOGICAL            :: use_analytical_GL_flux_config           = .FALSE.                ! Whether or not the analytical grounding line flux solution is used
  REAL(dp)           :: geothermal_heat_flux_config             = 1.72E06_dp             ! Geothermal Heat flux [J m^-2 yr^-1] Sclater et al. (1980)
  
  ! Some parameters for numerically solving the SSA
  REAL(dp)           :: SSA_RN_tol_config                       = 1E-5_dp                ! Successive solutions of the effective viscosity iteration must not differ by more than this amount (norm(N-Nprev) / norm(N) < RN_tol)
  REAL(dp)           :: SSA_norm_dUV_tol_config                 = 0.5_dp                 ! Successive solutions of UV in the effective viscosity iteration must not differ by more than this amount (on average)
  INTEGER            :: SSA_max_outer_loops_config              = 50                     ! Maximum number of effective viscosity iterations
  REAL(dp)           :: SSA_max_grad_N_config                   = 315.5694_dp            ! = 1e-5 * sec_per_year according to Frank Pattyn. If set to zero, the solver becomes identical to the old ANICE solver.
  REAL(dp)           :: SSA_max_residual_UV_config              = 2.5_dp                 ! The maximum residual in U and V in the SOR solver of the linearised SSA must drop below this value [m yr^-1]
  REAL(dp)           :: SSA_SOR_omega_config                    = 1.2_dp                 ! The over-relaxation parameter in the SOR solver of the linearised SSA
  INTEGER            :: SSA_max_inner_loops_config              = 10000                  ! Maximum number of iterations in the SOR solver of the linearised SSA
  
  ! SMB melt tuning
  ! ===============
  
  REAL(dp)           :: C_abl_constant_config     = -49._dp
  REAL(dp)           :: C_abl_Ts_config           = 10._dp
  REAL(dp)           :: C_abl_Q_config            = 0.0227_dp
  REAL(dp)           :: C_refr_config             = 0.051_dp
  
  ! Which data fields will be written to the help_fields output file
  ! ================================================================
  
  CHARACTER(LEN=256) :: help_field_01_config    = 'lat'
  CHARACTER(LEN=256) :: help_field_02_config    = 'lon'
  CHARACTER(LEN=256) :: help_field_03_config    = 'U_surf'
  CHARACTER(LEN=256) :: help_field_04_config    = 'V_surf'
  CHARACTER(LEN=256) :: help_field_05_config    = 'U_base'
  CHARACTER(LEN=256) :: help_field_06_config    = 'V_base'
  CHARACTER(LEN=256) :: help_field_07_config    = 'U_vav'
  CHARACTER(LEN=256) :: help_field_08_config    = 'V_vav'
  CHARACTER(LEN=256) :: help_field_09_config    = 'U_3D'
  CHARACTER(LEN=256) :: help_field_10_config    = 'V_3D'
  CHARACTER(LEN=256) :: help_field_11_config    = 'W_3D'
  CHARACTER(LEN=256) :: help_field_12_config    = 'mask'
  CHARACTER(LEN=256) :: help_field_13_config    = 'dHs_dx'
  CHARACTER(LEN=256) :: help_field_14_config    = 'dHs_dy'
  CHARACTER(LEN=256) :: help_field_15_config    = 'T2m_year'
  CHARACTER(LEN=256) :: help_field_16_config    = 'Precip_year'
  CHARACTER(LEN=256) :: help_field_17_config    = 'Albedo_year'
  CHARACTER(LEN=256) :: help_field_18_config    = 'SMB_year'
  CHARACTER(LEN=256) :: help_field_19_config    = 'BMB'
  CHARACTER(LEN=256) :: help_field_20_config    = 'T2m'
  CHARACTER(LEN=256) :: help_field_21_config    = 'Precip'
  CHARACTER(LEN=256) :: help_field_22_config    = 'Albedo'
  CHARACTER(LEN=256) :: help_field_23_config    = 'SMB'
  CHARACTER(LEN=256) :: help_field_24_config    = 'none'
  CHARACTER(LEN=256) :: help_field_25_config    = 'none'
  CHARACTER(LEN=256) :: help_field_26_config    = 'none'
  CHARACTER(LEN=256) :: help_field_27_config    = 'none'
  CHARACTER(LEN=256) :: help_field_28_config    = 'none'
  CHARACTER(LEN=256) :: help_field_29_config    = 'none'
  CHARACTER(LEN=256) :: help_field_30_config    = 'none'
  CHARACTER(LEN=256) :: help_field_31_config    = 'none'
  CHARACTER(LEN=256) :: help_field_32_config    = 'none'
  CHARACTER(LEN=256) :: help_field_33_config    = 'none'
  CHARACTER(LEN=256) :: help_field_34_config    = 'none'
  CHARACTER(LEN=256) :: help_field_35_config    = 'none'
  CHARACTER(LEN=256) :: help_field_36_config    = 'none'
  CHARACTER(LEN=256) :: help_field_37_config    = 'none'
  CHARACTER(LEN=256) :: help_field_38_config    = 'none'
  CHARACTER(LEN=256) :: help_field_39_config    = 'none'
  CHARACTER(LEN=256) :: help_field_40_config    = 'none'
  CHARACTER(LEN=256) :: help_field_41_config    = 'none'
  CHARACTER(LEN=256) :: help_field_42_config    = 'none'
  CHARACTER(LEN=256) :: help_field_43_config    = 'none'
  CHARACTER(LEN=256) :: help_field_44_config    = 'none'
  CHARACTER(LEN=256) :: help_field_45_config    = 'none'
  CHARACTER(LEN=256) :: help_field_46_config    = 'none'
  CHARACTER(LEN=256) :: help_field_47_config    = 'none'
  CHARACTER(LEN=256) :: help_field_48_config    = 'none'
  CHARACTER(LEN=256) :: help_field_49_config    = 'none'
  CHARACTER(LEN=256) :: help_field_50_config    = 'none'



  ! ==========================================================================
  ! The "C" type, which contains all the config parameters as fields.
  ! These will all be overwritten with the values of the "_config" variables,
  ! which are either the default values specified above, are the values
  ! specified from the external config file.
  ! ==========================================================================
  
  TYPE constants_type

    ! Time steps and range
    ! =====================
   
    REAL(dp)                            :: start_time_of_run
    REAL(dp)                            :: end_time_of_run
    REAL(dp)                            :: dt_coupling
    REAL(dp)                            :: dt_max
    REAL(dp)                            :: dt_thermo
    REAL(dp)                            :: dt_climate
    REAL(dp)                            :: dt_SMB
    REAL(dp)                            :: dt_BMB
    REAL(dp)                            :: dt_output
    REAL(dp)                            :: dt_mesh_min
    
    ! Which ice sheets do we simulate?
    ! ================================
    
    LOGICAL                             :: do_NAM
    LOGICAL                             :: do_EAS
    LOGICAL                             :: do_GRL
    LOGICAL                             :: do_ANT
    
    ! Benchmark experiments
    ! =====================
    
    LOGICAL                             :: do_benchmark_experiment
    CHARACTER(LEN=256)                  :: choice_benchmark_experiment
    REAL(dp)                            :: halfar_solution_H0
    REAL(dp)                            :: halfar_solution_R0
    REAL(dp)                            :: bueler_solution_lambda

    ! Mesh generation parameters
    ! ==========================

    INTEGER                             :: nconmax
    REAL(dp)                            :: alpha_min
    REAL(dp)                            :: dz_max_ice
    REAL(dp)                            :: res_max
    REAL(dp)                            :: res_min
    REAL(dp)                            :: res_max_margin
    REAL(dp)                            :: res_max_gl
    REAL(dp)                            :: res_max_cf
    REAL(dp)                            :: res_max_mountain
    REAL(dp)                            :: res_max_coast
    REAL(dp)                            :: mesh_fitness_threshold
  
    ! Resolutions of the different square grids
    ! =========================================
    
    REAL(dp)                            :: dx_grid_output
    REAL(dp)                            :: dx_grid_GIA
    REAL(dp)                            :: dx_grid_smooth
    
    ! High-resolution Points Of Interest (POIs) and transects
    ! =======================================================
    
    INTEGER                             :: nPOI_NAM
    INTEGER                             :: nPOI_EAS
    INTEGER                             :: nPOI_GRL
    INTEGER                             :: nPOI_ANT
    REAL(dp), DIMENSION(200)            :: POI_NAM_coordinates
    REAL(dp), DIMENSION(200)            :: POI_EAS_coordinates
    REAL(dp), DIMENSION(200)            :: POI_GRL_coordinates
    REAL(dp), DIMENSION(200)            :: POI_ANT_coordinates
    REAL(dp), DIMENSION(100)            :: POI_NAM_resolutions
    REAL(dp), DIMENSION(100)            :: POI_EAS_resolutions
    REAL(dp), DIMENSION(100)            :: POI_GRL_resolutions
    REAL(dp), DIMENSION(100)            :: POI_ANT_resolutions

    ! Whether or not to let UFEMISM dynamically create its own output folder
   ! =======================================================================
   
    LOGICAL                  :: create_new_output_dir
    CHARACTER(LEN=256)       :: output_dir
    LOGICAL                  :: do_write_debug_data

    ! Scaled vertical coordinate zeta  
    ! ===============================
     
    INTEGER                             :: nz       ! Number of grid points in vertical direction for thermodynamics in ice sheet
    REAL(dp), DIMENSION(:), ALLOCATABLE :: zeta

    ! Input data file paths
    ! =====================
    
    CHARACTER(LEN=4)         :: filetype_init_NAM
    CHARACTER(LEN=4)         :: filetype_init_EAS
    CHARACTER(LEN=4)         :: filetype_init_GRL
    CHARACTER(LEN=4)         :: filetype_init_ANT
                            
    CHARACTER(LEN=256)       :: filename_init_NAM
    CHARACTER(LEN=256)       :: filename_init_EAS
    CHARACTER(LEN=256)       :: filename_init_GRL
    CHARACTER(LEN=256)       :: filename_init_ANT
    
    CHARACTER(LEN=256)       :: filename_PD_NAM
    CHARACTER(LEN=256)       :: filename_PD_EAS
    CHARACTER(LEN=256)       :: filename_PD_GRL
    CHARACTER(LEN=256)       :: filename_PD_ANT
    
    CHARACTER(LEN=256)       :: filename_PD_obs_climate
    
    CHARACTER(LEN=256)       :: filename_insolation

    ! Ice dynamics & thermodynamics
    ! =============================
    
    REAL(dp)                 :: m_enh_sia
    REAL(dp)                 :: m_enh_ssa
    CHARACTER(LEN=256)       :: choice_sliding_law
    REAL(dp)                 :: C_sliding
    REAL(dp)                 :: m_sliding
    LOGICAL                  :: use_analytical_GL_flux
    REAL(dp)                 :: geothermal_heat_flux
  
    ! Some parameters for numerically solving the SSA
    REAL(dp)                 :: SSA_RN_tol
    REAL(dp)                 :: SSA_norm_dUV_tol
    INTEGER                  :: SSA_max_outer_loops
    REAL(dp)                 :: SSA_max_grad_N
    REAL(dp)                 :: SSA_max_residual_UV
    REAL(dp)                 :: SSA_SOR_omega
    INTEGER                  :: SSA_max_inner_loops
    
    ! SMB melt tuning
    ! ===============
    
    REAL(dp)                 :: C_abl_constant
    REAL(dp)                 :: C_abl_Ts
    REAL(dp)                 :: C_abl_Q
    REAL(dp)                 :: C_refr
  
    ! Which data fields will be written to the help_fields output file
    ! ================================================================
    
    CHARACTER(LEN=256)       :: help_field_01  
    CHARACTER(LEN=256)       :: help_field_02  
    CHARACTER(LEN=256)       :: help_field_03  
    CHARACTER(LEN=256)       :: help_field_04  
    CHARACTER(LEN=256)       :: help_field_05  
    CHARACTER(LEN=256)       :: help_field_06  
    CHARACTER(LEN=256)       :: help_field_07  
    CHARACTER(LEN=256)       :: help_field_08  
    CHARACTER(LEN=256)       :: help_field_09
    CHARACTER(LEN=256)       :: help_field_10  
    CHARACTER(LEN=256)       :: help_field_11  
    CHARACTER(LEN=256)       :: help_field_12  
    CHARACTER(LEN=256)       :: help_field_13  
    CHARACTER(LEN=256)       :: help_field_14  
    CHARACTER(LEN=256)       :: help_field_15  
    CHARACTER(LEN=256)       :: help_field_16  
    CHARACTER(LEN=256)       :: help_field_17  
    CHARACTER(LEN=256)       :: help_field_18  
    CHARACTER(LEN=256)       :: help_field_19
    CHARACTER(LEN=256)       :: help_field_20  
    CHARACTER(LEN=256)       :: help_field_21  
    CHARACTER(LEN=256)       :: help_field_22  
    CHARACTER(LEN=256)       :: help_field_23  
    CHARACTER(LEN=256)       :: help_field_24  
    CHARACTER(LEN=256)       :: help_field_25  
    CHARACTER(LEN=256)       :: help_field_26  
    CHARACTER(LEN=256)       :: help_field_27  
    CHARACTER(LEN=256)       :: help_field_28  
    CHARACTER(LEN=256)       :: help_field_29
    CHARACTER(LEN=256)       :: help_field_30  
    CHARACTER(LEN=256)       :: help_field_31  
    CHARACTER(LEN=256)       :: help_field_32  
    CHARACTER(LEN=256)       :: help_field_33  
    CHARACTER(LEN=256)       :: help_field_34  
    CHARACTER(LEN=256)       :: help_field_35  
    CHARACTER(LEN=256)       :: help_field_36  
    CHARACTER(LEN=256)       :: help_field_37  
    CHARACTER(LEN=256)       :: help_field_38  
    CHARACTER(LEN=256)       :: help_field_39
    CHARACTER(LEN=256)       :: help_field_40  
    CHARACTER(LEN=256)       :: help_field_41  
    CHARACTER(LEN=256)       :: help_field_42  
    CHARACTER(LEN=256)       :: help_field_43  
    CHARACTER(LEN=256)       :: help_field_44  
    CHARACTER(LEN=256)       :: help_field_45  
    CHARACTER(LEN=256)       :: help_field_46  
    CHARACTER(LEN=256)       :: help_field_47  
    CHARACTER(LEN=256)       :: help_field_48  
    CHARACTER(LEN=256)       :: help_field_49
    CHARACTER(LEN=256)       :: help_field_50
    
    ! Values to be filled into the total mask (used only for diagnostic output)
    ! ==========================================================================
    
    INTEGER                  :: type_land
    INTEGER                  :: type_ocean
    INTEGER                  :: type_lake
    INTEGER                  :: type_sheet
    INTEGER                  :: type_shelf
    INTEGER                  :: type_coast
    INTEGER                  :: type_margin
    INTEGER                  :: type_groundingline
    INTEGER                  :: type_calvingfront
    
   ! Parameters of the polar stereographic projections of the four model regions
   ! (These have to match the values used to create the input files!)
   ! ===========================================================================     
                          
    REAL(dp)                 :: lambda_M_NAM                           
    REAL(dp)                 :: lambda_M_EAS                           
    REAL(dp)                 :: lambda_M_GRL                           
    REAL(dp)                 :: lambda_M_ANT
    REAL(dp)                 :: phi_M_NAM
    REAL(dp)                 :: phi_M_EAS
    REAL(dp)                 :: phi_M_GRL
    REAL(dp)                 :: phi_M_ANT
    REAL(dp)                 :: alpha_stereo_NAM
    REAL(dp)                 :: alpha_stereo_EAS
    REAL(dp)                 :: alpha_stereo_GRL
    REAL(dp)                 :: alpha_stereo_ANT

  END TYPE constants_type


  ! ===============================================
  ! "C" is an instance of the "constants_type" type
  ! ===============================================
  
  TYPE(constants_type), SAVE :: C


CONTAINS
  SUBROUTINE read_main_config_file( config_filename)
    ! Use a NAMELIST containing all the "_config" variables to read
    ! an external config file, and overwrite the default values of
    ! the specified variables with the values from the file.
    
    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256),INTENT(IN) :: config_filename
    
    INTEGER, PARAMETER :: config_unit = 28 ! Unit number which is used for the configuration file.
    INTEGER            :: ios
    
    ! The NAMELIST that's used to read the external config file.
    NAMELIST /CONFIG/start_time_of_run_config,                   &                
                     end_time_of_run_config,                     &
                     dt_coupling_config,                         &
                     dt_max_config,                              &
                     dt_thermo_config,                           &
                     dt_climate_config,                          &
                     dt_SMB_config,                              &
                     dt_BMB_config,                              &
                     dt_output_config,                           &
                     dt_mesh_min_config,                         &
                     do_NAM_config,                              &
                     do_EAS_config,                              &
                     do_GRL_config,                              &
                     do_ANT_config,                              &
                     do_benchmark_experiment_config,             &
                     choice_benchmark_experiment_config,         &
                     halfar_solution_H0_config,                  &
                     halfar_solution_R0_config,                  &
                     bueler_solution_lambda_config,              &
                     nconmax_config,                             &
                     alpha_min_config,                           &
                     dz_max_ice_config,                          &
                     res_max_config,                             &
                     res_max_margin_config,                      &
                     res_max_gl_config,                          &
                     res_max_cf_config,                          &
                     res_max_mountain_config,                    &
                     res_max_coast_config,                       &
                     mesh_fitness_threshold_config,              &
                     dx_grid_output_config,                      &
                     dx_grid_GIA_config,                         &
                     dx_grid_smooth_config,                      &
                     nPOI_NAM_config,                            &
                     nPOI_EAS_config,                            &
                     nPOI_GRL_config,                            &
                     nPOI_ANT_config,                            &
                     POI_NAM_coordinates_config,                 &
                     POI_EAS_coordinates_config,                 &
                     POI_GRL_coordinates_config,                 &
                     POI_ANT_coordinates_config,                 &
                     POI_NAM_resolutions_config,                 &
                     POI_EAS_resolutions_config,                 &
                     POI_GRL_resolutions_config,                 &
                     POI_ANT_resolutions_config,                 &
                     create_new_output_dir_config,               &
                     output_dir_config,                          &
                     do_write_debug_data_config,                 &
                     nz_config,                                  &
                     zeta_config,                                &
                     filename_init_NAM_config,                   &
                     filename_init_EAS_config,                   &
                     filename_init_GRL_config,                   &
                     filename_init_ANT_config,                   &
                     filetype_init_NAM_config,                   &
                     filetype_init_EAS_config,                   &
                     filetype_init_GRL_config,                   &
                     filetype_init_ANT_config,                   &
                     filename_PD_NAM_config,                     &
                     filename_PD_EAS_config,                     &
                     filename_PD_GRL_config,                     &
                     filename_PD_ANT_config,                     &
                     filename_PD_obs_climate_config,             &
                     filename_insolation_config,                 &
                     m_enh_sia_config,                           &
                     m_enh_ssa_config,                           &
                     choice_sliding_law_config,                  &
                     C_sliding_config,                           &
                     m_sliding_config,                           &
                     use_analytical_GL_flux_config,              &
                     geothermal_heat_flux_config,                &
                     SSA_RN_tol_config,                          &
                     SSA_norm_dUV_tol_config,                    &
                     SSA_max_outer_loops_config,                 &
                     SSA_max_grad_N_config,                      &
                     SSA_max_residual_UV_config,                 &
                     SSA_SOR_omega_config,                       &
                     SSA_max_inner_loops_config,                 &
                     C_abl_constant_config,                      &
                     C_abl_Ts_config,                            &
                     C_abl_Q_config,                             &
                     C_refr_config,                              &
                     help_field_01_config,                       &
                     help_field_02_config,                       &
                     help_field_03_config,                       &
                     help_field_04_config,                       &
                     help_field_05_config,                       &
                     help_field_06_config,                       &
                     help_field_07_config,                       &
                     help_field_08_config,                       &
                     help_field_09_config,                       &
                     help_field_10_config,                       &
                     help_field_11_config,                       &
                     help_field_12_config,                       &
                     help_field_13_config,                       &
                     help_field_14_config,                       &
                     help_field_15_config,                       &
                     help_field_16_config,                       &
                     help_field_17_config,                       &
                     help_field_18_config,                       &
                     help_field_19_config,                       &
                     help_field_20_config,                       &
                     help_field_21_config,                       &
                     help_field_22_config,                       &
                     help_field_23_config,                       &
                     help_field_24_config,                       &
                     help_field_25_config,                       &
                     help_field_26_config,                       &
                     help_field_27_config,                       &
                     help_field_28_config,                       &
                     help_field_29_config,                       &
                     help_field_30_config,                       &
                     help_field_31_config,                       &
                     help_field_32_config,                       &
                     help_field_33_config,                       &
                     help_field_34_config,                       &
                     help_field_35_config,                       &
                     help_field_36_config,                       &
                     help_field_37_config,                       &
                     help_field_38_config,                       &
                     help_field_39_config,                       &
                     help_field_40_config,                       &
                     help_field_41_config,                       &
                     help_field_42_config,                       &
                     help_field_43_config,                       &
                     help_field_44_config,                       &
                     help_field_45_config,                       &
                     help_field_46_config,                       &
                     help_field_47_config,                       &
                     help_field_48_config,                       &
                     help_field_49_config,                       &
                     help_field_50_config
                      
     IF (config_filename == '') RETURN
      
     OPEN(UNIT=config_unit, FILE=TRIM(config_filename), STATUS='OLD', ACTION='READ', iostat=ios)
     IF(ios /= 0) THEN
       WRITE(UNIT=*, FMT='(/3A/)') ' ERROR: Could not open the configuration file: ', TRIM(config_filename)
       STOP
     END IF

     ! In the following statement the entire configuration file is read, using the namelist (NML=CONFIG)
     READ(UNIT=config_unit, NML=CONFIG, IOSTAT=ios)
     CLOSE(UNIT=config_unit)

     IF(ios /= 0) THEN
       WRITE(UNIT=*, FMT='(/3A)') ' ERROR while reading configuration file: ', TRIM(config_filename)
       STOP
     END IF     

  END SUBROUTINE read_main_config_file

  SUBROUTINE initialize_main_constants
    ! Overwrite the values in the fields of the "C" type with the values
    ! of the "_config" variables, some which by now have had their default
    ! values overwritten by the values specified in the external config file.
    
    IMPLICIT NONE
    
    INTEGER :: cerr, ierr
    
   ! Time steps and range
   !=====================
   
    C%start_time_of_run                   = start_time_of_run_config
    C%end_time_of_run                     = end_time_of_run_config
    C%dt_coupling                         = dt_coupling_config
    C%dt_max                              = dt_max_config
    C%dt_thermo                           = dt_thermo_config
    C%dt_climate                          = dt_climate_config
    C%dt_SMB                              = dt_SMB_config
    C%dt_BMB                              = dt_BMB_config
    C%dt_output                           = dt_output_config
    C%dt_mesh_min                         = dt_mesh_min_config
    
    ! Which ice sheets do we simulate?
    ! ================================
    
    C%do_NAM                              = do_NAM_config
    C%do_EAS                              = do_EAS_config
    C%do_GRL                              = do_GRL_config
    C%do_ANT                              = do_ANT_config
    
    ! Benchmark experiments
    ! =====================
    
    C%do_benchmark_experiment             = do_benchmark_experiment_config
    C%choice_benchmark_experiment         = choice_benchmark_experiment_config
    
    ! A quick check to see if the specified benchmark experiment
    ! has actually been implemented
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        C%do_NAM = .FALSE.
        C%do_EAS = .FALSE.
        C%do_GRL = .FALSE.
        C%do_ANT = .TRUE.
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialize_main_constants!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    C%halfar_solution_H0                  = halfar_solution_H0_config
    C%halfar_solution_R0                  = halfar_solution_R0_config
    C%bueler_solution_lambda              = bueler_solution_lambda_config

    ! Mesh generation parameters
    ! ==========================
    
    C%nconmax                             = nconmax_config
    C%alpha_min                           = alpha_min_config
    C%dz_max_ice                          = dz_max_ice_config
    C%res_max                             = res_max_config
    C%res_max_margin                      = res_max_margin_config
    C%res_max_gl                          = res_max_gl_config
    C%res_max_cf                          = res_max_cf_config
    C%res_max_mountain                    = res_max_mountain_config
    C%res_max_coast                       = res_max_coast_config
    C%mesh_fitness_threshold              = mesh_fitness_threshold_config
    
    ! The smallest allowed resolution
    C%res_min = MIN( MIN( MIN( MIN( C%res_max_margin, C%res_max_gl), C%res_max_cf), C%res_max_mountain), C%res_max_coast)
  
    ! Resolutions of the different square grids
    ! =========================================
    
    C%dx_grid_output                      = dx_grid_output_config
    C%dx_grid_GIA                         = dx_grid_GIA_config
    C%dx_grid_smooth                      = dx_grid_smooth_config
    
    ! High-resolution Points Of Interest (POIs)
    ! =========================================
    
    C%nPOI_NAM                            = nPOI_NAM_config    
    C%nPOI_EAS                            = nPOI_EAS_config    
    C%nPOI_GRL                            = nPOI_GRL_config    
    C%nPOI_ANT                            = nPOI_ANT_config
    C%POI_NAM_coordinates                 = POI_NAM_coordinates_config
    C%POI_EAS_coordinates                 = POI_EAS_coordinates_config
    C%POI_GRL_coordinates                 = POI_GRL_coordinates_config
    C%POI_ANT_coordinates                 = POI_ANT_coordinates_config
    C%POI_NAM_resolutions                 = POI_NAM_resolutions_config
    C%POI_EAS_resolutions                 = POI_EAS_resolutions_config
    C%POI_GRL_resolutions                 = POI_GRL_resolutions_config
    C%POI_ANT_resolutions                 = POI_ANT_resolutions_config
    
    ! Whether or not to let UFEMISM dynamically create its own output folder
   ! =======================================================================
   
    C%create_new_output_dir               = create_new_output_dir_config
    C%output_dir                          = output_dir_config
    C%do_write_debug_data                 = do_write_debug_data_config

    ! Scaled vertical coordinate zeta  
    ! ===============================
    
    C%nz     = nz_config
    ALLOCATE( C%zeta( C%nz))
    C%zeta   = zeta_config( 1:C%nz)

    ! Input data file paths
    ! =====================
    
    C%filetype_init_NAM                   = filetype_init_NAM_config
    C%filetype_init_EAS                   = filetype_init_EAS_config
    C%filetype_init_GRL                   = filetype_init_GRL_config
    C%filetype_init_ANT                   = filetype_init_ANT_config
    
    IF (.NOT. (C%filetype_init_NAM == 'cart' .OR. C%filetype_init_NAM == 'mesh')) THEN
      WRITE(0,*) 'Config ERROR - filetype_init can only be "cart" or "mesh"!'
      STOP
    END IF
    IF (.NOT. (C%filetype_init_EAS == 'cart' .OR. C%filetype_init_EAS == 'mesh')) THEN
      WRITE(0,*) 'Config ERROR - filetype_init can only be "cart" or "mesh"!'
      STOP
    END IF
    IF (.NOT. (C%filetype_init_GRL == 'cart' .OR. C%filetype_init_GRL == 'mesh')) THEN
      WRITE(0,*) 'Config ERROR - filetype_init can only be "cart" or "mesh"!'
      STOP
    END IF
    IF (.NOT. (C%filetype_init_ANT == 'cart' .OR. C%filetype_init_ANT == 'mesh')) THEN
      WRITE(0,*) 'Config ERROR - filetype_init can only be "cart" or "mesh"!'
      STOP
    END IF
    
    C%filename_init_NAM                   = filename_init_NAM_config
    C%filename_init_EAS                   = filename_init_EAS_config
    C%filename_init_GRL                   = filename_init_GRL_config
    C%filename_init_ANT                   = filename_init_ANT_config
    
    C%filename_PD_NAM                     = filename_PD_NAM_config
    C%filename_PD_EAS                     = filename_PD_EAS_config
    C%filename_PD_GRL                     = filename_PD_GRL_config
    C%filename_PD_ANT                     = filename_PD_ANT_config
    
    C%filename_PD_obs_climate             = filename_PD_obs_climate_config
    
    C%filename_insolation                 = filename_insolation_config

    ! Ice dynamics & thermodynamics
    ! =============================
    
    C%m_enh_sia                           = m_enh_sia_config
    C%m_enh_ssa                           = m_enh_ssa_config    
    C%choice_sliding_law                  = choice_sliding_law_config
    C%C_sliding                           = C_sliding_config
    C%m_sliding                           = m_sliding_config
    C%use_analytical_GL_flux              = use_analytical_GL_flux_config
    C%geothermal_heat_flux                = geothermal_heat_flux_config
  
    ! Some parameters for numerically solving the SSA
    C%SSA_RN_tol                          = SSA_RN_tol_config
    C%SSA_norm_dUV_tol                    = SSA_norm_dUV_tol_config
    C%SSA_max_outer_loops                 = SSA_max_outer_loops_config
    C%SSA_max_grad_N                      = SSA_max_grad_N_config
    C%SSA_max_residual_UV                 = SSA_max_residual_UV_config
    C%SSA_SOR_omega                       = SSA_SOR_omega_config
    C%SSA_max_inner_loops                 = SSA_max_inner_loops_config
    
    ! SMB melt tuning
    ! ===============
    
    C%C_abl_constant                      = C_abl_constant_config
    C%C_abl_Ts                            = C_abl_Ts_config
    C%C_abl_Q                             = C_abl_Q_config
    C%C_refr                              = C_refr_config
  
    ! Which data fields will be written to the help_fields output file
    ! ================================================================
    
    C%help_field_01                       = help_field_01_config
    C%help_field_02                       = help_field_02_config
    C%help_field_03                       = help_field_03_config
    C%help_field_04                       = help_field_04_config
    C%help_field_05                       = help_field_05_config
    C%help_field_06                       = help_field_06_config
    C%help_field_07                       = help_field_07_config
    C%help_field_08                       = help_field_08_config
    C%help_field_09                       = help_field_09_config
    C%help_field_10                       = help_field_10_config
    C%help_field_11                       = help_field_11_config
    C%help_field_12                       = help_field_12_config
    C%help_field_13                       = help_field_13_config
    C%help_field_14                       = help_field_14_config
    C%help_field_15                       = help_field_15_config
    C%help_field_16                       = help_field_16_config
    C%help_field_17                       = help_field_17_config
    C%help_field_18                       = help_field_18_config
    C%help_field_19                       = help_field_19_config
    C%help_field_20                       = help_field_20_config
    C%help_field_21                       = help_field_21_config
    C%help_field_22                       = help_field_22_config
    C%help_field_23                       = help_field_23_config
    C%help_field_24                       = help_field_24_config
    C%help_field_25                       = help_field_25_config
    C%help_field_26                       = help_field_26_config
    C%help_field_27                       = help_field_27_config
    C%help_field_28                       = help_field_28_config
    C%help_field_29                       = help_field_29_config
    C%help_field_30                       = help_field_30_config
    C%help_field_31                       = help_field_31_config
    C%help_field_32                       = help_field_32_config
    C%help_field_33                       = help_field_33_config
    C%help_field_34                       = help_field_34_config
    C%help_field_35                       = help_field_35_config
    C%help_field_36                       = help_field_36_config
    C%help_field_37                       = help_field_37_config
    C%help_field_38                       = help_field_38_config
    C%help_field_39                       = help_field_39_config
    C%help_field_40                       = help_field_40_config
    C%help_field_41                       = help_field_41_config
    C%help_field_42                       = help_field_42_config
    C%help_field_43                       = help_field_43_config
    C%help_field_44                       = help_field_44_config
    C%help_field_45                       = help_field_45_config
    C%help_field_46                       = help_field_46_config
    C%help_field_47                       = help_field_47_config
    C%help_field_48                       = help_field_48_config
    C%help_field_49                       = help_field_49_config
    C%help_field_50                       = help_field_50_config
    
    ! Values to be filled into the total mask (used only for diagnostic output)
    ! ==========================================================================

    C%type_land           = 0
    C%type_ocean          = 1
    C%type_lake           = 2
    C%type_sheet          = 3
    C%type_shelf          = 4
    C%type_coast          = 5
    C%type_margin         = 6
    C%type_groundingline  = 7
    C%type_calvingfront   = 8
    
   ! Parameters of the polar stereographic projections of the four model regions
   ! (These have to match the values used to create the input files!)
   ! ===========================================================================  
  
    C%lambda_M_NAM     = 265._dp
    C%lambda_M_EAS     = 40._dp
    C%lambda_M_GRL     = 320._dp
    C%lambda_M_ANT     = 0._dp
    C%phi_M_NAM        = 62._dp
    C%phi_M_EAS        = 70._dp
    C%phi_M_GRL        = 72._dp
    C%phi_M_ANT        = -90._dp
    C%alpha_stereo_NAM = 165.0923_dp
    C%alpha_stereo_EAS = 165.04_dp
    C%alpha_stereo_GRL = 164.85_dp
    C%alpha_stereo_ANT = 165.0263_dp

  END SUBROUTINE initialize_main_constants

  SUBROUTINE create_output_dir

    IMPLICIT NONE

    CHARACTER(20)              :: output_folder_name

    INTEGER,    DIMENSION(8)   :: values
    LOGICAL                    :: ex

    CALL date_and_time(VALUES=values)

    ! Get proper year (assume we're still in the 21st century...)
    output_folder_name(1:10) = 'results_20'
    SELECT CASE( FLOOR(REAL(values(1))/10._dp)-200)
     CASE(0)
     output_folder_name(11:11) = '0'
     CASE(1)
     output_folder_name(11:11) = '1'
     CASE(2)
     output_folder_name(11:11) = '2'
     CASE(3)
     output_folder_name(11:11) = '3'
     CASE(4)
     output_folder_name(11:11) = '4'
     CASE(5)
     output_folder_name(11:11) = '5'
     CASE(6)
     output_folder_name(11:11) = '6'
     CASE(7)
     output_folder_name(11:11) = '7'
     CASE(8)
     output_folder_name(11:11) = '8'
     CASE(9)
     output_folder_name(11:11) = '9'
     CASE DEFAULT
     WRITE(0,*) 'make_output_folder: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( MOD(values(1),10))
     CASE(0)
     output_folder_name(12:12) = '0'
     CASE(1)
     output_folder_name(12:12) = '1'
     CASE(2)
     output_folder_name(12:12) = '2'
     CASE(3)
     output_folder_name(12:12) = '3'
     CASE(4)
     output_folder_name(12:12) = '4'
     CASE(5)
     output_folder_name(12:12) = '5'
     CASE(6)
     output_folder_name(12:12) = '6'
     CASE(7)
     output_folder_name(12:12) = '7'
     CASE(8)
     output_folder_name(12:12) = '8'
     CASE(9)
     output_folder_name(12:12) = '9'
     CASE DEFAULT
     WRITE(0,*) 'make_output_folder: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( values(2))
     CASE(1)
     output_folder_name(13:14) = '01'
     CASE(2)
     output_folder_name(13:14) = '02'
     CASE(3)
     output_folder_name(13:14) = '03'
     CASE(4)
     output_folder_name(13:14) = '04'
     CASE(5)
     output_folder_name(13:14) = '05'
     CASE(6)
     output_folder_name(13:14) = '06'
     CASE(7)
     output_folder_name(13:14) = '07'
     CASE(8)
     output_folder_name(13:14) = '08'
     CASE(9)
     output_folder_name(13:14) = '09'
     CASE(10)
     output_folder_name(13:14) = '10'
     CASE(11)
     output_folder_name(13:14) = '11'
     CASE(12)
     output_folder_name(13:14) = '12'
     CASE DEFAULT
     WRITE(0,*) 'make_output_folder: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( FLOOR(REAL(values(3))/10._dp))
     CASE(0)
     output_folder_name(15:15) = '0'
     CASE(1)
     output_folder_name(15:15) = '1'
     CASE(2)
     output_folder_name(15:15) = '2'
     CASE(3)
     output_folder_name(15:15) = '3'
     CASE DEFAULT
     WRITE(0,*) 'make_output_folder: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( MOD(values(3),10))
     CASE(0)
     output_folder_name(16:16) = '0'
     CASE(1)
     output_folder_name(16:16) = '1'
     CASE(2)
     output_folder_name(16:16) = '2'
     CASE(3)
     output_folder_name(16:16) = '3'
     CASE(4)
     output_folder_name(16:16) = '4'
     CASE(5)
     output_folder_name(16:16) = '5'
     CASE(6)
     output_folder_name(16:16) = '6'
     CASE(7)
     output_folder_name(16:16) = '7'
     CASE(8)
     output_folder_name(16:16) = '8'
     CASE(9)
     output_folder_name(16:16) = '9'
     CASE DEFAULT
     WRITE(0,*) 'make_output_folder: ERROR retrieving date and time!'
    END SELECT

    output_folder_name(17:20) = '_001'

    INQUIRE( FILE=TRIM(output_folder_name)//'/.', EXIST=ex )

    DO WHILE (ex)

     IF      (output_folder_name(20:20) == '0') THEN
      output_folder_name(20:20) = '1'
     ELSE IF (output_folder_name(20:20) == '1') THEN
      output_folder_name(20:20) = '2'
     ELSE IF (output_folder_name(20:20) == '2') THEN
      output_folder_name(20:20) = '3'
     ELSE IF (output_folder_name(20:20) == '3') THEN
      output_folder_name(20:20) = '4'
     ELSE IF (output_folder_name(20:20) == '4') THEN
      output_folder_name(20:20) = '5'
     ELSE IF (output_folder_name(20:20) == '5') THEN
      output_folder_name(20:20) = '6'
     ELSE IF (output_folder_name(20:20) == '6') THEN
      output_folder_name(20:20) = '7'
     ELSE IF (output_folder_name(20:20) == '7') THEN
      output_folder_name(20:20) = '8'
     ELSE IF (output_folder_name(20:20) == '8') THEN
      output_folder_name(20:20) = '9'
     ELSE IF (output_folder_name(20:20) == '9') THEN
      output_folder_name(20:20) = '0'

      IF      (output_folder_name(19:19) == '0') THEN
       output_folder_name(19:19) = '1'
      ELSE IF (output_folder_name(19:19) == '1') THEN
       output_folder_name(19:19) = '2'
      ELSE IF (output_folder_name(19:19) == '2') THEN
       output_folder_name(19:19) = '3'
      ELSE IF (output_folder_name(19:19) == '3') THEN
       output_folder_name(19:19) = '4'
      ELSE IF (output_folder_name(19:19) == '4') THEN
       output_folder_name(19:19) = '5'
      ELSE IF (output_folder_name(19:19) == '5') THEN
       output_folder_name(19:19) = '6'
      ELSE IF (output_folder_name(19:19) == '6') THEN
       output_folder_name(19:19) = '7'
      ELSE IF (output_folder_name(19:19) == '7') THEN
       output_folder_name(19:19) = '8'
      ELSE IF (output_folder_name(19:19) == '8') THEN
       output_folder_name(19:19) = '9'
      ELSE IF (output_folder_name(19:19) == '9') THEN
       output_folder_name(19:19) = '0'

       IF      (output_folder_name(18:18) == '0') THEN
        output_folder_name(18:18) = '1'
       ELSE IF (output_folder_name(18:18) == '1') THEN
        output_folder_name(18:18) = '2'
       ELSE IF (output_folder_name(18:18) == '2') THEN
        output_folder_name(18:18) = '3'
       ELSE IF (output_folder_name(18:18) == '3') THEN
        output_folder_name(18:18) = '4'
       ELSE IF (output_folder_name(18:18) == '4') THEN
        output_folder_name(18:18) = '5'
       ELSE IF (output_folder_name(18:18) == '5') THEN
        output_folder_name(18:18) = '6'
       ELSE IF (output_folder_name(18:18) == '6') THEN
        output_folder_name(18:18) = '7'
       ELSE IF (output_folder_name(18:18) == '7') THEN
        output_folder_name(18:18) = '8'
       ELSE IF (output_folder_name(18:18) == '8') THEN
        output_folder_name(18:18) = '9'
       ELSE IF (output_folder_name(18:18) == '9') THEN
        output_folder_name(18:18) = '0'
       END IF

      END IF

     END IF

     INQUIRE( FILE=TRIM(output_folder_name)//'/.', EXIST=ex )

    END DO

    IF (C%create_new_output_dir) THEN
      C%output_dir = TRIM(output_folder_name)
    END IF
    
    C%output_dir = TRIM(C%output_dir(1:255)) // '/'

    CALL system('mkdir ' // TRIM(C%output_dir))
    WRITE(0,*) 'Output directory: ', TRIM(C%output_dir)
    WRITE(0,*) ''

  END SUBROUTINE create_output_dir

END MODULE configuration_module
