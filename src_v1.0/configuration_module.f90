MODULE configuration_module
  IMPLICIT NONE

  ! PRECISION
  ! =========
  ! The kind of real numbers used by default throughout the program.
  ! Reals should be declared as:
  !
  ! REAL(dp) :: example
  !  or
  ! REAL(KIND=dp) :: example
  !
  ! dp must be a PARAMETER:
  INTEGER, PARAMETER                :: dp  = KIND(1.0D0)  ! Kind of double precision numbers.

  ! CONFIG VARIABLES:
  !==================
  ! Variables which are eventually set by the read_main_config_file subroutine (if in the CONFIG file):

  ! DIAGNOSTICA
  ! ===========
  LOGICAL :: save_workspace_data_config = .TRUE.
  
  LOGICAL            :: create_new_output_dir_config = .TRUE.
  CHARACTER(LEN=256) :: output_dir_config = 'results_UFEMISM'
  
  ! Output data to be saved
  LOGICAL :: save_mesh_data_config         = .TRUE.
  LOGICAL :: save_ice_dynamics_config      = .TRUE.
  LOGICAL :: save_SMB_components_config    = .TRUE.
  
  ! benchmark experiments
  LOGICAL          :: do_eismint_experiment_config            = .FALSE.
  CHARACTER(LEN=1) :: choice_eismint_experiment_config        = '0'
  REAL(dp)         :: eismint_xmin_config                     = -750000._dp
  REAL(dp)         :: eismint_xmax_config                     =  750000._dp
  REAL(dp)         :: eismint_ymin_config                     = -750000._dp
  REAL(dp)         :: eismint_ymax_config                     =  750000._dp
  LOGICAL          :: do_halfar_solution_config               = .FALSE.
  REAL(dp)         :: halfar_solution_H0_config               = 5000._dp
  REAL(dp)         :: halfar_solution_R0_config               = 300000._dp
  LOGICAL          :: do_bueler_solution_config               = .FALSE.
  REAL(dp)         :: bueler_solution_lambda_config           = 5._dp
  LOGICAL          :: do_mismip3d_experiment_config           = .FALSE.
  
  ! Which ice sheets do we simulate
  CHARACTER(LEN=4) :: which_icesheets_config = 'NEGA'

  ! MESH GENERATION
  ! ===============
  INTEGER  :: nconmax_config          = 32           ! Maximum number of vertex connections
  LOGICAL  :: do_rectangular_config   = .FALSE.      ! If we want ot use a semi-rectangular mesh
  REAL(dp) :: res_rectangular_config  = 40._dp       ! Horizontal resolution for the semi-rectangular mesh (km)
  REAL(dp) :: mesh_alpha_min_config   = 0.4363_dp    ! Minimum internal angle of triangles (25 degrees)
  REAL(dp) :: dz_max_ice_config       = 20000._dp    ! Maximum allowed 2nd order surface deviation over ice  
  REAL(dp) :: res_max_config          = 800._dp      ! Maximum allowed resolution                            (km)
  REAL(dp) :: res_max_margin_config   = 40._dp       ! Maximum allowed resolution over land-based ice margin (km)
  REAL(dp) :: res_max_gl_config       = 40._dp       !                                 grounding line        (km)
  REAL(dp) :: res_max_cf_config       = 40._dp       !                                 calving front         (km)
  REAL(dp) :: res_max_mountain_config = 40._dp       !                                 mountains             (km)
  REAL(dp) :: res_max_coast_config    = 40._dp       !                                 coastline             (km)
  REAL(dp) :: mesh_fitness_threshold_config = 0.95_dp ! Minimum allowed mesh fitness (fraction of triangles that are not Bad) before mesh updating
  
  ! High-resolution Points Of Interest
  INTEGER                    :: mesh_nPOI_NAM_config            = 0       ! Number of POIs
  INTEGER                    :: mesh_nPOI_EAS_config            = 0
  INTEGER                    :: mesh_nPOI_GRL_config            = 0
  INTEGER                    :: mesh_nPOI_ANT_config            = 0
  REAL(dp), DIMENSION(200)   :: mesh_POI_NAM_coordinates_config = 0._dp   ! [lat,lon] coordinates op POIs (degrees)
  REAL(dp), DIMENSION(200)   :: mesh_POI_EAS_coordinates_config = 0._dp
  REAL(dp), DIMENSION(200)   :: mesh_POI_GRL_coordinates_config = 0._dp
  REAL(dp), DIMENSION(200)   :: mesh_POI_ANT_coordinates_config = 0._dp
  REAL(dp), DIMENSION(100)   :: mesh_POI_NAM_resolutions_config = 0._dp   ! Required resolution at POIs (km)
  REAL(dp), DIMENSION(100)   :: mesh_POI_EAS_resolutions_config = 0._dp
  REAL(dp), DIMENSION(100)   :: mesh_POI_GRL_resolutions_config = 0._dp
  REAL(dp), DIMENSION(100)   :: mesh_POI_ANT_resolutions_config = 0._dp

  INTEGER                :: velocity_scaling_nvals_config = 4
  REAL(dp),DIMENSION(10) :: velocity_scaling_vel_config = &
   (/  -4._dp,   3._dp,   5._dp,  8._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp/)
  REAL(dp),DIMENSION(10) :: velocity_scaling_res_config = &
   (/1000._dp, 200._dp, 100._dp, 20._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp/)
   
  REAL(dp), DIMENSION(7) :: nested_mesh_resolutions_config = &
   (/200._dp, 100._dp, 50._dp, 25._dp, 12._dp, 6._dp, 3._dp/)

  ! Relative grid spacing in vertical direction in ice sheet. If k is a counter through the vertical layers, k=1
  ! at the surface corresponding with zeta=0, and k=NZ at the bottom of the ice sheet corresponding with zeta=1:
  ! This CONFIG variable zeta_config is declared as a large array, because fortran does not allow a CONFIG/NAMELIST
  ! variable which is ALLOCATABLE, only the C%NZ first elements of this array will be used and have to be specified
  ! in the CONFIG file.
  INTEGER                        :: NZ_config   = 15
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

  ! TIME STEPS AND LIMITS
  ! =====================
  REAL(dp) :: start_time_of_run_config   = 0.0_dp       ! Start time (in years) of the simulations
  REAL(dp) :: end_time_of_run_config     = 1000.0_dp    ! End   time (in years) of the simulations
  REAL(dp) :: dt_coupling_config         = 100._dp      ! Interval of coupling (in years) between the four ice-sheets  
  REAL(dp) :: dt_max_config              = 10.0_dp      ! Maximum time step (in years) of the ice model
  REAL(dp) :: dt_thermo_config           = 10.0_dp      ! Time step (in years) for updating thermodynamics
  REAL(dp) :: dt_climate_config          = 10._dp       ! Time step (in years) for updating the climate
  REAL(dp) :: dt_SMB_config              = 10._dp       ! Time step (in years) for updating the SMB
  REAL(dp) :: dt_BMB_config              = 10._dp       ! Time step (in years) for updating the BMB
  REAL(dp) :: dt_output_config           = 100.0_dp     ! Time step (in years) for writing output
  REAL(dp) :: dt_mesh_min_config         = 50._dp       ! Minimum amount of time (in years) between mesh updates

  ! Different data file paths
  ! =========================
  
  CHARACTER(LEN=4)   :: filetype_init_NAM_config = 'cart'
  CHARACTER(LEN=4)   :: filetype_init_EAS_config = 'cart'
  CHARACTER(LEN=4)   :: filetype_init_GRL_config = 'cart'
  CHARACTER(LEN=4)   :: filetype_init_ANT_config = 'cart'
  
  ! Initial data (NetCDF)
  CHARACTER(LEN=256) :: filename_init_NAM_config = '/Users/berends/Documents/Models/UFEMISM/data/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_init_EAS_config = '/Users/berends/Documents/Models/UFEMISM/data/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_init_GRL_config = '/Users/berends/Documents/Models/UFEMISM/data/Greenland_BedMachine_5km_noEllesmere.nc'
  CHARACTER(LEN=256) :: filename_init_ANT_config = '/Users/berends/Documents/Models/UFEMISM/data/Antarctica_Bedmap2_Rignot_5km.nc'
  ! PD reference data (NetCDF)
  CHARACTER(LEN=256) :: filename_PD_NAM_config   = '/Users/berends/Documents/Models/UFEMISM/data/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_PD_EAS_config   = '/Users/berends/Documents/Models/UFEMISM/data/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_PD_GRL_config   = '/Users/berends/Documents/Models/UFEMISM/data/Greenland_BedMachine_5km_noEllesmere.nc'
  CHARACTER(LEN=256) :: filename_PD_ANT_config   = '/Users/berends/Documents/Models/UFEMISM/data/Antarctica_Bedmap2_Rignot_5km.nc'  
  ! Present-day observed climate (ERA40) (NetCDF)
  CHARACTER(LEN=256) :: filename_PD_obs_climate_config = '/Users/berends/Documents/Datasets/ERA40/ERA40_climate_global.nc'  
  ! Insolation forcing (Laskar et al., 2004) (TXT)
  CHARACTER(LEN=256) :: filename_insolation_config = '/Users/berends/Documents/Datasets/Insolation_laskar/Insolation_Laskar_etal_2004.nc'

  ! ==========
  LOGICAL            :: create_new_output_dir
  CHARACTER(LEN=256) :: output_dir

  ! ICE DYNAMICS & THERMODYNAMICS
  ! =============================
  REAL(dp)           :: m_enh_sia_config                   = 5.0_dp
  REAL(dp)           :: m_enh_ssa_config                   = 0.7_dp
  CHARACTER(LEN=256) :: choice_sliding_law_config          = 'Coulomb'
  REAL(dp)           :: C_sliding_config                   = 1.0E7_dp     ! Factor   in Weertman sliding law
  REAL(dp)           :: m_sliding_config                   = 1._dp/3._dp  ! Exponent in Weertman sliding law
  LOGICAL            :: use_thermodynamics_config          = .TRUE.
  REAL(dp)           :: geothermal_heat_flux_config        = 1.72E06_dp   ! Geothermal Heat flux [J m^-2 yr^-1] Sclater et al. (1980)
  
  ! SMB melt tuning
  ! ===============
  REAL(dp)           :: C_abl_constant_config     = -49._dp
  REAL(dp)           :: C_abl_Ts_config           = 10._dp
  REAL(dp)           :: C_abl_Q_config            = 0.0227_dp
  REAL(dp)           :: C_refr_config             = 0.051_dp


  ! TYPE DEFINITIONS
  !=================

  ! This TYPE contains all the information once the CONFIG file is read never will change during the run of the program
  TYPE constants_type

    ! DIAGNOSTICA
    ! ===========
    LOGICAL                             :: save_workspace_data
    
    ! Output data to be saved
    LOGICAL                             :: save_mesh_data
    LOGICAL                             :: save_ice_dynamics
    LOGICAL                             :: save_SMB_components
    
    ! Which ice sheets we want to simulate
    CHARACTER(LEN=4)                    :: which_icesheets
    
    ! benchmark experiments
    LOGICAL                             :: do_eismint_experiment
    CHARACTER(LEN=1)                    :: choice_eismint_experiment
    REAL(dp)                            :: eismint_xmin
    REAL(dp)                            :: eismint_xmax
    REAL(dp)                            :: eismint_ymin
    REAL(dp)                            :: eismint_ymax
    LOGICAL                             :: do_halfar_solution
    REAL(dp)                            :: halfar_solution_H0
    REAL(dp)                            :: halfar_solution_R0 
    LOGICAL                             :: do_bueler_solution
    REAL(dp)                            :: bueler_solution_lambda
    LOGICAL                             :: do_mismip3d_experiment

    ! MESH GENERATION
    ! ===============

    INTEGER                             :: nconmax
    LOGICAL                             :: do_rectangular
    REAL(dp)                            :: res_rectangular
    REAL(dp)                            :: mesh_alpha_min
    REAL(dp)                            :: dz_max_ice
    REAL(dp)                            :: res_max
    REAL(dp)                            :: res_min
    REAL(dp)                            :: res_max_margin
    REAL(dp)                            :: res_max_gl
    REAL(dp)                            :: res_max_cf
    REAL(dp)                            :: res_max_mountain
    REAL(dp)                            :: res_max_coast
    INTEGER                             :: velocity_scaling_nvals
    REAL(dp), DIMENSION(10)             :: velocity_scaling_vel
    REAL(dp), DIMENSION(10)             :: velocity_scaling_res
    REAL(dp), DIMENSION(7)              :: nested_mesh_resolutions
    REAL(dp)                            :: mesh_fitness_threshold    
    INTEGER                             :: mesh_nPOI_NAM
    INTEGER                             :: mesh_nPOI_EAS
    INTEGER                             :: mesh_nPOI_GRL
    INTEGER                             :: mesh_nPOI_ANT
    REAL(dp), DIMENSION(200)            :: mesh_POI_NAM_coordinates
    REAL(dp), DIMENSION(200)            :: mesh_POI_EAS_coordinates
    REAL(dp), DIMENSION(200)            :: mesh_POI_GRL_coordinates
    REAL(dp), DIMENSION(200)            :: mesh_POI_ANT_coordinates
    REAL(dp), DIMENSION(100)            :: mesh_POI_NAM_resolutions
    REAL(dp), DIMENSION(100)            :: mesh_POI_EAS_resolutions
    REAL(dp), DIMENSION(100)            :: mesh_POI_GRL_resolutions
    REAL(dp), DIMENSION(100)            :: mesh_POI_ANT_resolutions

    ! Scaled vertical coordinate zeta   
    INTEGER                             :: NZ       ! Number of grid points in vertical direction for thermodynamics in ice sheet
    REAL(dp), DIMENSION(:), ALLOCATABLE :: zeta


                            ! TIME CONDITIONS
                            !===============================
    REAL(dp)                 :: start_time_of_run  ! Start time of run (in years)
    REAL(dp)                 :: end_time_of_run    ! End   time of run (in years)
    REAL(dp)                 :: dt_coupling        ! Interval of coupling (in years) between the four ice-sheets
    REAL(dp)                 :: dt_max             ! Maximum time step (in years) of the ice model
    REAL(dp)                 :: dt_thermo          ! Time step (in years) for updating the thermodynamics
    REAL(dp)                 :: dt_climate         ! Time step (in years) for updating the climate
    REAL(dp)                 :: dt_SMB             ! Time step (in years) for updating the SMB
    REAL(dp)                 :: dt_BMB             ! Time step (in years) for updating the BMB
    REAL(dp)                 :: dt_output          ! Time step (in years) for writing output
    REAL(dp)                 :: dt_mesh_min        ! Minimum amount of time (in years) between mesh updates

                            ! TYPES OF LAND (used in the mask)
                            !===============================
    INTEGER                  :: type_land
    INTEGER                  :: type_ocean
    INTEGER                  :: type_lake
    INTEGER                  :: type_sheet
    INTEGER                  :: type_shelf
    INTEGER                  :: type_coast
    INTEGER                  :: type_margin
    INTEGER                  :: type_groundingline
    INTEGER                  :: type_calvingfront


                            ! Source data file paths
                            !=================================================
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

                            ! The output folder (will be created by the model itself)
                            !========================================================
    LOGICAL                  :: create_new_output_dir
    CHARACTER(LEN=256)       :: output_dir

                            ! Ice dynamics & thermodynamics
                            ! =============================
    REAL(dp)                 :: m_enh_sia
    REAL(dp)                 :: m_enh_ssa
    CHARACTER(LEN=256)       :: choice_sliding_law
    REAL(dp)                 :: C_sliding
    REAL(dp)                 :: m_sliding
    LOGICAL                  :: use_thermodynamics
    REAL(dp)                 :: geothermal_heat_flux
    
                            ! SMB melt tuning
                            ! ===============
    REAL(dp)                 :: C_abl_constant
    REAL(dp)                 :: C_abl_Ts
    REAL(dp)                 :: C_abl_Q
    REAL(dp)                 :: C_refr

                            ! MATHEMATICAL AND PHYSICAL CONSTANTS
                            !====================================
    REAL(dp)                 :: earth_radius
    REAL(dp)                 :: pi
    REAL(dp)                 :: deg2rad             ! Conversion factor between radians and degrees
    REAL(dp)                 :: rad2deg             ! Conversion factor between degrees and radians
    REAL(dp)                 :: sec_per_year        ! seconds per year
    
                            ! Parameters of the polar stereographic projections of the four model regions
                            ! (used for mapping between global grids and regional meshes)
                            ! ===========================================================================                            
    REAL(dp)                 :: lambda_M_NAM                           
    REAL(dp)                 :: lambda_M_EAS                           
    REAL(dp)                 :: lambda_M_GRL                           
    REAL(dp)                 :: lambda_M_ANT
    REAL(dp)                 :: phi_M_NAM
    REAL(dp)                 :: phi_M_EAS
    REAL(dp)                 :: phi_M_GRL
    REAL(dp)                 :: phi_M_ANT
    REAL(dp)                 :: alpha_NAM
    REAL(dp)                 :: alpha_EAS
    REAL(dp)                 :: alpha_GRL
    REAL(dp)                 :: alpha_ANT

  END TYPE constants_type


  ! C is the 'struct' containing all the Constants from the CONFIG file and/or the defaults
  TYPE(constants_type), SAVE :: C



CONTAINS
  SUBROUTINE read_main_config_file(config_filename)
    ! This subroutine reads some of the variables defined in the grid subroutine
    ! from a configuration file. The name of the configuration file should be specified
    ! on the command line. If no name is specified on the command line, then the default
    ! values are used.
    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256),INTENT(IN) :: config_filename
    
    INTEGER, PARAMETER :: config_unit = 28 ! Unit number which is used for the configuration file.
    INTEGER            :: ios

     ! List of items in the configuration file:
     NAMELIST /CONFIG/save_mesh_data_config                     , &
                      save_workspace_data_config                , &
                      save_ice_dynamics_config                  , &
                      save_SMB_components_config                , &
                      create_new_output_dir_config              , &
                      output_dir_config                         , &
                      which_icesheets_config                    , &
                      do_eismint_experiment_config              , &
                      choice_eismint_experiment_config          , &
                      eismint_xmin_config                       , &
                      eismint_xmax_config                       , &
                      eismint_ymin_config                       , &
                      eismint_ymax_config                       , &
                      do_halfar_solution_config                 , &
                      halfar_solution_H0_config                 , &
                      halfar_solution_R0_config                 , &
                      do_bueler_solution_config                 , &
                      bueler_solution_lambda_config             , &
                      do_mismip3d_experiment_config             , &
                      nconmax_config                            , &
                      do_rectangular_config                     , &
                      res_rectangular_config                    , &
                      mesh_alpha_min_config                     , &
                      dz_max_ice_config                         , &
                      res_max_config                            , &
                      res_max_margin_config                     , &
                      res_max_gl_config                         , &
                      res_max_cf_config                         , &
                      res_max_mountain_config                   , &
                      res_max_coast_config                      , &
                      mesh_fitness_threshold_config             , &
                      mesh_nPOI_NAM_config                      , &
                      mesh_nPOI_EAS_config                      , &
                      mesh_nPOI_GRL_config                      , &
                      mesh_nPOI_ANT_config                      , &
                      mesh_POI_NAM_coordinates_config           , &
                      mesh_POI_EAS_coordinates_config           , &
                      mesh_POI_GRL_coordinates_config           , &
                      mesh_POI_ANT_coordinates_config           , &
                      mesh_POI_NAM_resolutions_config           , &
                      mesh_POI_EAS_resolutions_config           , &
                      mesh_POI_GRL_resolutions_config           , &
                      mesh_POI_ANT_resolutions_config           , &
                      velocity_scaling_nvals_config             , &
                      velocity_scaling_vel_config               , &
                      velocity_scaling_res_config               , &
                      nested_mesh_resolutions_config            , &
                      NZ_config                                 , &
                      zeta_config                               , &
                      start_time_of_run_config                  , &
                      end_time_of_run_config                    , &
                      dt_coupling_config                        , &
                      dt_max_config                             , &
                      dt_thermo_config                          , &
                      dt_climate_config                         , &
                      dt_SMB_config                             , &
                      dt_BMB_config                             , &
                      dt_output_config                          , &
                      dt_mesh_min_config                        , &
                      filetype_init_NAM_config                  , &
                      filetype_init_EAS_config                  , &
                      filetype_init_GRL_config                  , &
                      filetype_init_ANT_config                  , &
                      filename_init_NAM_config                  , &
                      filename_init_EAS_config                  , &
                      filename_init_GRL_config                  , &
                      filename_init_ANT_config                  , &
                      filename_PD_NAM_config                    , &
                      filename_PD_EAS_config                    , &
                      filename_PD_GRL_config                    , &
                      filename_PD_ANT_config                    , &
                      filename_PD_obs_climate_config            , &
                      filename_insolation_config                , &
                      m_enh_sia_config                          , &
                      m_enh_ssa_config                          , &
                      choice_sliding_law_config                 , &
                      C_sliding_config                          , &
                      m_sliding_config                          , &
                      use_thermodynamics_config                 , &
                      geothermal_heat_flux_config               , &
                      C_abl_constant_config                     , &
                      C_abl_Ts_config                           , &
                      C_abl_Q_config                            , &
                      C_refr_config
                      
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

  SUBROUTINE initialize_main_constants()
    ! This routine puts all the constants which will never change during the run after the CONFIG file
    ! has been read, into a special constant 'struct'
    IMPLICIT NONE

    ! Number of grid points in vertical direction for thermodynamics in ice sheet, if k is a counter
    ! through the vertical layers, k=1 at the surface and k=NZ at the bottom of the ice sheet:
    C%NZ     = NZ_config
    ALLOCATE(C%zeta(C%NZ))
    C%zeta   = zeta_config(1:C%NZ) ! Fortran does not allow a CONFIG/NAMELIST variable to be ALLOCATABLE, therefore this way

    ! Diagnostica
    C%save_workspace_data                 = save_workspace_data_config
    
    C%create_new_output_dir               = create_new_output_dir_config
    C%output_dir                          = output_dir_config
    
    ! EISMINT experiments
    C%do_eismint_experiment               = do_eismint_experiment_config
    C%choice_eismint_experiment           = choice_eismint_experiment_config
    C%eismint_xmin                        = eismint_xmin_config
    C%eismint_xmax                        = eismint_xmax_config
    C%eismint_ymin                        = eismint_ymin_config
    C%eismint_ymax                        = eismint_ymax_config
    C%do_halfar_solution                  = do_halfar_solution_config
    C%halfar_solution_H0                  = halfar_solution_H0_config
    C%halfar_solution_R0                  = halfar_solution_R0_config
    C%do_bueler_solution                  = do_bueler_solution_config
    C%bueler_solution_lambda              = bueler_solution_lambda_config
    C%do_mismip3d_experiment              = do_mismip3d_experiment_config
    
    ! Output data to be saved
    C%save_mesh_data                      = save_mesh_data_config
    C%save_ice_dynamics                   = save_ice_dynamics_config
    C%save_SMB_components                 = save_SMB_components_config
    
    ! Which ice sheets we want to simulate
    C%which_icesheets                     = which_icesheets_config
    IF (do_eismint_experiment_config)  C%which_icesheets = 'FFFA'
    IF (do_mismip3d_experiment_config) C%which_icesheets = 'FFFA'

    ! Mesh generation
    C%nconmax                             = nconmax_config
    C%do_rectangular                      = do_rectangular_config
    C%res_rectangular                     = res_rectangular_config
    C%mesh_alpha_min                      = mesh_alpha_min_config
    C%dz_max_ice                          = dz_max_ice_config
    C%res_max                             = res_max_config
    C%res_max_margin                      = res_max_margin_config
    C%res_max_gl                          = res_max_gl_config
    C%res_max_cf                          = res_max_cf_config
    C%res_max_mountain                    = res_max_mountain_config
    C%res_max_coast                       = res_max_coast_config
    
    ! The smallest allowed resolution
    C%res_min = MIN( MIN( MIN( MIN( C%res_max_margin, C%res_max_gl), C%res_max_cf), C%res_max_mountain), C%res_max_coast)
    
    C%mesh_fitness_threshold              = mesh_fitness_threshold_config
    C%velocity_scaling_nvals              = velocity_scaling_nvals_config
    C%velocity_scaling_vel                = velocity_scaling_vel_config
    C%velocity_scaling_res                = velocity_scaling_res_config
    C%nested_mesh_resolutions             = nested_mesh_resolutions_config
    
    C%mesh_nPOI_NAM                       = mesh_nPOI_NAM_config    
    C%mesh_nPOI_EAS                       = mesh_nPOI_EAS_config    
    C%mesh_nPOI_GRL                       = mesh_nPOI_GRL_config    
    C%mesh_nPOI_ANT                       = mesh_nPOI_ANT_config
    C%mesh_POI_NAM_coordinates            = mesh_POI_NAM_coordinates_config
    C%mesh_POI_EAS_coordinates            = mesh_POI_EAS_coordinates_config
    C%mesh_POI_GRL_coordinates            = mesh_POI_GRL_coordinates_config
    C%mesh_POI_ANT_coordinates            = mesh_POI_ANT_coordinates_config
    C%mesh_POI_NAM_resolutions            = mesh_POI_NAM_resolutions_config
    C%mesh_POI_EAS_resolutions            = mesh_POI_EAS_resolutions_config
    C%mesh_POI_GRL_resolutions            = mesh_POI_GRL_resolutions_config
    C%mesh_POI_ANT_resolutions            = mesh_POI_ANT_resolutions_config

    C%start_time_of_run                   = start_time_of_run_config            ! Start time of run; see the CONFIG file or the default
    C%end_time_of_run                     = end_time_of_run_config              ! End time of run; see the CONFIG file or the default
    C%dt_coupling                         = dt_coupling_config                  ! Interval of coupling (in years) between the four ice-sheets
    C%dt_max                              = dt_max_config                       ! Maximum time step (in years) of the ice model
    C%dt_thermo                           = dt_thermo_config                    ! Time step (in years) used for the thermodynamics
    C%dt_climate                          = dt_climate_config                   ! Time step (in years) used for the climate
    C%dt_SMB                              = dt_SMB_config                       ! Time step (in years) used for the SMB
    C%dt_BMB                              = dt_BMB_config                       ! Time step (in years) used for the BMB
    C%dt_output                           = dt_output_config                    ! Time step (in years) for writing output files
    C%dt_mesh_min                         = dt_mesh_min_config                  ! Minimum amount of time (in years) between mesh updates

  ! We use a mask array to store the type of land in each grid point. The following types of land are allowed:
  ! ==========================================================================================================
    C%type_land                           = 0
    C%type_ocean                          = 1
    C%type_lake                           = 2
    C%type_sheet                          = 3
    C%type_shelf                          = 4
    C%type_coast                          = 5
    C%type_margin                         = 6
    C%type_groundingline                  = 7
    C%type_calvingfront                   = 8

    C%m_enh_sia                           = m_enh_sia_config
    C%m_enh_ssa                           = m_enh_ssa_config
    
    C%choice_sliding_law                  = choice_sliding_law_config
    C%C_sliding                           = C_sliding_config
    C%m_sliding                           = m_sliding_config
    IF (.NOT. (C%choice_sliding_law == 'Coulomb' .OR. C%choice_sliding_law == 'Weertman')) THEN
      WRITE(0,*) ' ERROR: choice_sliding_law can only be "Coulomb" or "Weertman"!'
      STOP
    END IF
    
    C%use_thermodynamics                  = use_thermodynamics_config
    C%geothermal_heat_flux                = geothermal_heat_flux_config
    
    C%C_abl_constant                      = C_abl_constant_config
    C%C_abl_Ts                            = C_abl_Ts_config
    C%C_abl_Q                             = C_abl_Q_config
    C%C_refr                              = C_refr_config
    
    C%earth_radius                        = 6.371221e6_dp
    C%pi                                  = 2._dp*ACOS(0._dp)                ! Just pi=3.14159... exactly
    C%deg2rad                             = C%pi/180._dp                     ! Conversion factor between radians and degrees
    C%rad2deg                             = 180._dp/C%pi                     ! Conversion factor between degrees and radians
    C%sec_per_year                        = 31556943.36_dp                   ! = 365.2424 * 24 * 3600
    
    
  ! Parameters of the polar stereographic projections of the four model regions
  ! (used for mapping between global grids and regional meshes)
  ! =========================================================================== 
    C%lambda_M_NAM    = 265._dp
    C%lambda_M_EAS    = 40._dp
    C%lambda_M_GRL    = 320._dp
    C%lambda_M_ANT    = 0._dp
    C%phi_M_NAM       = 62._dp
    C%phi_M_EAS       = 70._dp
    C%phi_M_GRL       = 72._dp
    C%phi_M_ANT       = -90._dp
    C%alpha_NAM       = 165.0923_dp
    C%alpha_EAS       = 165.04_dp
    C%alpha_GRL       = 164.85_dp
    C%alpha_ANT       = 165.0263_dp

  ! Source data file paths
  ! ======================
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

  END SUBROUTINE initialize_main_constants

  SUBROUTINE create_output_dir()

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
