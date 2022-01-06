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
  
  REAL(dp)            :: start_time_of_run_config                    = 0.0_dp                           ! Start time (in years) of the simulations
  REAL(dp)            :: end_time_of_run_config                      = 50000.0_dp                       ! End   time (in years) of the simulations
  REAL(dp)            :: dt_coupling_config                          = 100._dp                          ! Interval of coupling (in years) between the four ice-sheets  
  REAL(dp)            :: dt_max_config                               = 10.0_dp                          ! Maximum time step (in years) of the ice model
  REAL(dp)            :: dt_thermo_config                            = 10.0_dp                          ! Time step (in years) for updating thermodynamics
  REAL(dp)            :: dt_climate_config                           = 10._dp                           ! Time step (in years) for updating the climate
  REAL(dp)            :: dt_ocean_config                             = 10._dp                           ! Time step (in years) for updating the ocean
  REAL(dp)            :: dt_SMB_config                               = 10._dp                           ! Time step (in years) for updating the SMB
  REAL(dp)            :: dt_BMB_config                               = 10._dp                           ! Time step (in years) for updating the BMB
  REAL(dp)            :: dt_output_config                            = 5000.0_dp                        ! Time step (in years) for writing output
  REAL(dp)            :: dt_mesh_min_config                          = 50._dp                           ! Minimum amount of time (in years) between mesh updates
  REAL(dp)            :: dt_bedrock_ELRA_config                      = 100._dp                          ! Time step (in years) for updating the bedrock deformation rate with the ELRA model
  
  ! Debugging
  ! =========
  
  LOGICAL             :: do_write_debug_data_config                  = .FALSE.                          ! Whether or not the debug NetCDF file should be created and written to
  LOGICAL             :: do_check_for_NaN_config                     = .FALSE.                          ! Whether or not fields should be checked for NaN values           
  LOGICAL             :: do_write_memory_tracker_config              = .FALSE.                          ! Whether or not the memory use tracker should be used                
  
  ! Which ice sheets do we simulate?
  ! ================================
  
  LOGICAL             :: do_NAM_config                               = .FALSE.                          ! North America
  LOGICAL             :: do_EAS_config                               = .FALSE.                          ! Eurasia
  LOGICAL             :: do_GRL_config                               = .FALSE.                          ! Greenland
  LOGICAL             :: do_ANT_config                               = .TRUE.                           ! Antarctica                  
  
  ! Domain size for the four regions
  ! ================================
  
  REAL(dp)            :: xmin_NAM_config                             = -3600000._dp                     ! Western  boundary     of the North America domain [m]
  REAL(dp)            :: xmax_NAM_config                             =  3600000._dp                     ! Eastern  boundary     of the North America domain [m]
  REAL(dp)            :: ymin_NAM_config                             = -2400000._dp                     ! Southern boundary     of the North America domain [m]
  REAL(dp)            :: ymax_NAM_config                             =  2400000._dp                     ! Northern boundary     of the North America domain [m]
  
  REAL(dp)            :: xmin_EAS_config                             = -3400000._dp                     ! Western  boundary     of the Eurasia domain [m]
  REAL(dp)            :: xmax_EAS_config                             =  3400000._dp                     ! Eastern  boundary     of the Eurasia domain [m]
  REAL(dp)            :: ymin_EAS_config                             = -2080000._dp                     ! Southern boundary     of the Eurasia domain [m]
  REAL(dp)            :: ymax_EAS_config                             =  2080000._dp                     ! Northern boundary     of the Eurasia domain [m]
  
  REAL(dp)            :: xmin_GRL_config                             =  -830000._dp                     ! Western  boundary     of the Greenland domain [m]
  REAL(dp)            :: xmax_GRL_config                             =   830000._dp                     ! Eastern  boundary     of the Greenland domain [m]
  REAL(dp)            :: ymin_GRL_config                             = -1430000._dp                     ! Southern boundary     of the Greenland domain [m]
  REAL(dp)            :: ymax_GRL_config                             =  1430000._dp                     ! Northern boundary     of the Greenland domain [m]
  
  REAL(dp)            :: xmin_ANT_config                             = -3300000._dp                     ! Western  boundary     of the Antarctica domain [m]
  REAL(dp)            :: xmax_ANT_config                             =  3300000._dp                     ! Eastern  boundary     of the Antarctica domain [m]
  REAL(dp)            :: ymin_ANT_config                             = -3300000._dp                     ! Southern boundary     of the Antarctica domain [m]
  REAL(dp)            :: ymax_ANT_config                             =  3300000._dp                     ! Northern boundary     of the Antarctica domain [m]
  
  ! Whether or not the simulation is a restart of a previous simulation
  ! ===================================================================
  
  LOGICAL             :: is_restart_config                           = .FALSE.
  REAL(dp)            :: time_to_restart_from_config                 = 0._dp                            ! Can be different from C%start_time_of_run, though this will issue a warning
  
  ! Initial model state when restarting from a previous run
  CHARACTER(LEN=256)  :: filename_restart_NAM_config                 = 'filename_restart_NAM_placeholder'
  CHARACTER(LEN=256)  :: filename_restart_EAS_config                 = 'filename_restart_EAS_placeholder'
  CHARACTER(LEN=256)  :: filename_restart_GRL_config                 = 'filename_restart_GRL_placeholder'
  CHARACTER(LEN=256)  :: filename_restart_ANT_config                 = 'filename_restart_ANT_placeholder'
  
  ! Benchmark experiments
  ! =====================
  
  LOGICAL             :: do_benchmark_experiment_config              = .TRUE.
  CHARACTER(LEN=256)  :: choice_benchmark_experiment_config          = 'EISMINT_I'
  REAL(dp)            :: SSA_icestream_m_config                      = 1                                ! Values tested by Schoof are 1, 10, and 20
  REAL(dp)            :: ISMIP_HOM_L_config                          = 160000.0                         ! Domain size of the ISMIP-HOM benchmarks
  CHARACTER(LEN=256)  :: ISMIP_HOM_E_Arolla_filename_config          = 'arolla100.dat'                  ! Path to the Haut Glacier d'Arolla input file
  LOGICAL             :: MISMIPplus_do_tune_A_for_GL_config          = .FALSE.                          ! Whether or not the flow factor A should be tuned for the GL position
  REAL(dp)            :: MISMIPplus_xGL_target_config                = 450000._dp                       ! Mid-channel GL position to tune the flow factor A for
  REAL(dp)            :: MISMIPplus_A_flow_initial_config            = 2.0E-17_dp                       ! Initial flow factor before tuning (or throughout the run when tuning is not used)
  CHARACTER(LEN=256)  :: MISMIPplus_scenario_config                  = ''                               ! Choose between the five MISMIP+  scenarios from Cornford   et al. (2020): ice0, ice1ra, ice1rr, ice2ra, ice2rr
  CHARACTER(LEN=256)  :: MISOMIP1_scenario_config                    = ''                               ! Choose between the four MISOMIP+ scenarios from Asay-Davis et al. (2016): IceOcean1ra, IceOcean1rr, IceOcean2ra, IceOcean2rr

  ! Mesh generation parameters
  ! ==========================
  
  INTEGER             :: nconmax_config                              = 16                               ! Maximum number of vertex connections
  REAL(dp)            :: alpha_min_config                            = 0.4363_dp                        ! Minimum internal angle of triangles (25 degrees)
  REAL(dp)            :: dz_max_ice_config                           = 20000._dp                        ! Maximum allowed 2nd order surface deviation over ice 
  REAL(dp)            :: res_max_config                              = 800._dp                          ! Maximum allowed resolution                            [km]
  REAL(dp)            :: res_max_margin_config                       = 40._dp                           ! Maximum allowed resolution over land-based ice margin [km]
  REAL(dp)            :: res_max_gl_config                           = 40._dp                           !                                 grounding line        [km]
  REAL(dp)            :: res_max_cf_config                           = 40._dp                           !                                 calving front         [km]
  REAL(dp)            :: res_max_mountain_config                     = 40._dp                           !                                 mountains             [km]
  REAL(dp)            :: res_max_coast_config                        = 40._dp                           !                                 coastline             [km]
  REAL(dp)            :: mesh_fitness_threshold_config               = 0.95_dp                          ! Minimum allowed mesh fitness (fraction of triangles that are not Bad) before mesh updating
  
  ! Resolutions of the different square grids
  ! =========================================
  
  REAL(dp)            :: dx_grid_output_config                       = 40000._dp                        ! Resolution of the square grid used for writing output                       [m]
  REAL(dp)            :: dx_grid_GIA_config                          = 100000._dp                       ! Resolution of the square grid used for GIA modelling (ELRA or SELEN)        [m]
  REAL(dp)            :: dx_grid_smooth_config                       = 50000._dp                        ! Resolution of the square grid used for data smoothing in the climate matrix [m]
  
  ! High-resolution Points Of Interest (POIs)
  ! =========================================
  
  INTEGER                    :: nPOI_NAM_config                      = 0                                ! Number of POIs
  INTEGER                    :: nPOI_EAS_config                      = 0
  INTEGER                    :: nPOI_GRL_config                      = 0
  INTEGER                    :: nPOI_ANT_config                      = 0
  REAL(dp), DIMENSION(200)   :: POI_NAM_coordinates_config           = 0._dp                            ! [lat,lon] coordinates op POIs (degrees)
  REAL(dp), DIMENSION(200)   :: POI_EAS_coordinates_config           = 0._dp
  REAL(dp), DIMENSION(200)   :: POI_GRL_coordinates_config           = 0._dp
  REAL(dp), DIMENSION(200)   :: POI_ANT_coordinates_config           = 0._dp
  REAL(dp), DIMENSION(100)   :: POI_NAM_resolutions_config           = 0._dp                            ! Required resolution at POIs (km)
  REAL(dp), DIMENSION(100)   :: POI_EAS_resolutions_config           = 0._dp
  REAL(dp), DIMENSION(100)   :: POI_GRL_resolutions_config           = 0._dp
  REAL(dp), DIMENSION(100)   :: POI_ANT_resolutions_config           = 0._dp
  
  ! Whether or not to let UFEMISM dynamically create its own output folder.
  ! This works fine locally, on LISA its better to use a fixed folder name.
  ! =======================================================================
  
  LOGICAL             :: create_new_output_dir_config                = .TRUE.
  CHARACTER(LEN=256)  :: output_dir_config                           = 'results_UFEMISM'

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

  ! Reference geometries (initial, present-day, and GIA equilibrium)
  ! ================================================================
  
  ! Initial geometry
  CHARACTER(LEN=256)  :: choice_refgeo_init_NAM_config               = 'realistic'                      ! Choice of initial geometry for North America; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_init_EAS_config               = 'realistic'                      ! Choice of initial geometry for Eurasia      ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_init_GRL_config               = 'realistic'                      ! Choice of initial geometry for Greenland    ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_init_ANT_config               = 'realistic'                      ! Choice of initial geometry for Antarctica   ; can be "idealised", "realistic", or "restart"
  REAL(dp)            :: time_to_restart_from_NAM_config             = 0._dp                            ! Can be different from C%start_time_of_run, though this will issue a warning
  REAL(dp)            :: time_to_restart_from_EAS_config             = 0._dp
  REAL(dp)            :: time_to_restart_from_GRL_config             = 0._dp
  REAL(dp)            :: time_to_restart_from_ANT_config             = 0._dp
  CHARACTER(LEN=256)  :: choice_refgeo_init_idealised_config         = 'flatearth'                      ! Choice of idealised initial geometry; see "generate_idealised_geometry" in reference_fields_module for options
  REAL(dp)            :: dx_refgeo_init_idealised_config             = 5000._dp                         ! Resolution of square grid used for idealised initial geometry
  CHARACTER(LEN=256)  :: filename_refgeo_init_NAM_config             = '/Users/berends/Documents/Datasets/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_init_EAS_config             = '/Users/berends/Documents/Datasets/ETOPO1/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_init_GRL_config             = '/Users/berends/Documents/Datasets/Bedmachine_Greenland/Greenland_BedMachine_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_init_ANT_config             = '/Users/berends/Documents/Datasets/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'
  
  ! Present-day geometry
  CHARACTER(LEN=256)  :: choice_refgeo_PD_NAM_config                 = 'realistic'                      ! Choice of present-day geometry for North America; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_PD_EAS_config                 = 'realistic'                      ! Choice of present-day geometry for Eurasia      ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_PD_GRL_config                 = 'realistic'                      ! Choice of present-day geometry for Greenland    ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_PD_ANT_config                 = 'realistic'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_PD_idealised_config           = 'flatearth'                      ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
  REAL(dp)            :: dx_refgeo_PD_idealised_config               = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry
  CHARACTER(LEN=256)  :: filename_refgeo_PD_NAM_config               = '/Users/berends/Documents/Datasets/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_PD_EAS_config               = '/Users/berends/Documents/Datasets/ETOPO1/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_PD_GRL_config               = '/Users/berends/Documents/Datasets/Bedmachine_Greenland/Greenland_BedMachine_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_PD_ANT_config               = '/Users/berends/Documents/Datasets/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'
  
  ! GIA equilibrium geometry
  CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_NAM_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for North America; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_EAS_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for Eurasia      ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_GRL_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for Greenland    ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_ANT_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for Antarctica   ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_idealised_config        = 'flatearth'                      ! Choice of idealised GIA equilibrium geometry; see "generate_idealised_geometry" in reference_fields_module for options
  REAL(dp)            :: dx_refgeo_GIAeq_idealised_config            = 5000._dp                         ! Resolution of square grid used for idealised GIA equilibrium geometry
  CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_NAM_config            = '/Users/berends/Documents/Datasets/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_EAS_config            = '/Users/berends/Documents/Datasets/ETOPO1/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_GRL_config            = '/Users/berends/Documents/Datasets/Bedmachine_Greenland/Greenland_BedMachine_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_ANT_config            = '/Users/berends/Documents/Datasets/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'

  LOGICAL             :: remove_Lake_Vostok_config                   = .TRUE.

  ! Input data file paths
  ! =====================
   
  ! Insolation forcing (NetCDF) (Laskar et al., 2004)
  CHARACTER(LEN=256)  :: filename_insolation_config                  = '/Datasets/Insolation_laskar/Insolation_Laskar_etal_2004.nc'
  
  ! CO2 record (ASCII text file, so the number of rows needs to be specified)
  CHARACTER(LEN=256)  :: filename_CO2_record_config                  = 'Datasets/CO2/EPICA_CO2_Bereiter_2015_100yr.dat'
  INTEGER             :: CO2_record_length_config                    = 8001
  
  ! d18O record (ASCII text file, so the number of rows needs to be specified)
  CHARACTER(LEN=256)  :: filename_d18O_record_config                 = 'Datasets/d18O/Ahn2017_d18O.dat'
  INTEGER             :: d18O_record_length_config                   = 2051
  
  ! Geothermal heat flux
  CHARACTER(LEN=256)  :: choice_geothermal_heat_flux_config          = 'spatial'                        ! Choice of geothermal heat flux; can be 'constant' or 'spatial'
  REAL(dp)            :: constant_geothermal_heat_flux_config        = 1.72E06_dp                       ! Geothermal Heat flux [J m^-2 yr^-1] Sclater et al. (1980)
  CHARACTER(LEN=256)  :: filename_geothermal_heat_flux_config        = '/Users/berends/Documents/Datasets/GHF/geothermal_heatflux_ShapiroRitzwoller2004_global_1x1_deg.nc'

  ! Ice dynamics - velocity
  ! =======================
  
  CHARACTER(LEN=256)  :: choice_ice_dynamics_config                  = 'DIVA'                           ! Choice of ice-dynamica approximation: "none" (= fixed geometry), "SIA", "SSA", "SIA/SSA", "DIVA"
  REAL(dp)            :: n_flow_config                               = 3.0_dp                           ! Exponent in Glen's flow law
  REAL(dp)            :: m_enh_sheet_config                          = 1.0_dp                           ! Ice flow enhancement factor for grounded ice
  REAL(dp)            :: m_enh_shelf_config                          = 1.0_dp                           ! Ice flow enhancement factor for floating ice
  CHARACTER(LEN=256)  :: choice_ice_margin_config                    = 'infinite_slab'                  ! Choice of ice margin boundary conditions: "BC", "infinite_slab"
  LOGICAL             :: include_SSADIVA_crossterms_config           = .TRUE.                           ! Whether or not to include the "cross-terms" of the SSA/DIVA
  LOGICAL             :: do_GL_subgrid_friction_config               = .TRUE.                           ! Whether or not to scale basal friction with the sub-grid grounded fraction (needed to get proper GL migration; only turn this off for showing the effect on the MISMIP_mod results!)
  LOGICAL             :: do_smooth_geometry_config                   = .FALSE.                          ! Whether or not to smooth the model geometry (bedrock + initial ice thickness)
  REAL(dp)            :: r_smooth_geometry_config                    = 0.5_dp                           ! Geometry smoothing radius (in number of grid cells)
  
  ! Some parameters for numerically solving the SSA/DIVA
  REAL(dp)            :: DIVA_visc_it_norm_dUV_tol_config            = 1E-2_dp                          ! Successive solutions of UV in the effective viscosity iteration must not differ by more than this amount (on average)
  INTEGER             :: DIVA_visc_it_nit_config                     = 50                               ! Maximum number of effective viscosity iterations
  REAL(dp)            :: DIVA_visc_it_relax_config                   = 0.7_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
  REAL(dp)            :: DIVA_beta_max_config                        = 1E20_dp                          ! The DIVA is not solved (i.e. velocities are assumed to be zero) for beta values larger than this
  REAL(dp)            :: DIVA_err_lim_config                         = 1E-5_dp                          ! The DIVA is not refined (i.e. velocities are no longer updated with SOR) wherever successive velocity iterations change the velocity by less than this amount
  REAL(dp)            :: DIVA_vel_max_config                         = 5000._dp                         ! DIVA velocities are limited to this value (u,v evaluated separately)
  REAL(dp)            :: DIVA_vel_min_config                         = 1E-5_dp                          ! DIVA velocities below this value are set to zero (u,v evaluated separately)
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_west_config              = 'infinite'                       ! Boundary conditions for the ice velocity field at the domain boundary in the DIVA
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_east_config              = 'infinite'                       ! Allowed choices: "infinite", "periodic", "zero"
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_south_config             = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_north_config             = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_west_config              = 'infinite'                       ! Boundary conditions for the ice velocity field at the domain boundary in the DIVA
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_east_config              = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_south_config             = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_north_config             = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_choice_matrix_solver_config            = 'PETSc'                          ! Choice of matrix solver for the ice velocity equations: "SOR", "PETSc"
  INTEGER             :: DIVA_SOR_nit_config                         = 10000                            ! DIVA SOR   solver - maximum number of iterations
  REAL(dp)            :: DIVA_SOR_tol_config                         = 2.5_dp                           ! DIVA SOR   solver - stop criterion, absolute difference
  REAL(dp)            :: DIVA_SOR_omega_config                       = 1.3_dp                           ! DIVA SOR   solver - over-relaxation parameter
  REAL(dp)            :: DIVA_PETSc_rtol_config                      = 0.01_dp                          ! DIVA PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
  REAL(dp)            :: DIVA_PETSc_abstol_config                    = 0.01_dp                          ! DIVA PETSc solver - stop criterion, absolute difference

  ! Ice dynamics - time integration
  ! ===============================
  
  CHARACTER(LEN=256)  :: choice_timestepping_config                  = 'pc'                             ! Choice of timestepping method: "direct", "pc" (NOTE: 'direct' does not work with DIVA ice dynamcis!)
  CHARACTER(LEN=256)  :: choice_ice_integration_method_config        = 'explicit'                       ! Choice of ice thickness integration scheme: "none" (i.e. unchanging geometry), "explicit", "semi-implicit"
  CHARACTER(LEN=256)  :: dHi_choice_matrix_solver_config             = 'SOR'                            ! Choice of matrix solver for the semi-implicit ice thickness equation: "SOR", "PETSc"
  INTEGER             :: dHi_SOR_nit_config                          = 3000                             ! dHi SOR   solver - maximum number of iterations
  REAL(dp)            :: dHi_SOR_tol_config                          = 2.5_dp                           ! dHi SOR   solver - stop criterion, absolute difference
  REAL(dp)            :: dHi_SOR_omega_config                        = 1.3_dp                           ! dHi SOR   solver - over-relaxation parameter
  REAL(dp)            :: dHi_PETSc_rtol_config                       = 0.001_dp                         ! dHi PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
  REAL(dp)            :: dHi_PETSc_abstol_config                     = 0.001_dp                         ! dHi PETSc solver - stop criterion, absolute difference
  
  ! Predictor-corrector ice-thickness update
  REAL(dp)            :: pc_epsilon_config                           = 3._dp                            ! Target truncation error in dHi_dt [m/yr] (epsilon in Robinson et al., 2020, Eq. 33)
  REAL(dp)            :: pc_k_I_config                               = 0.2_dp                           ! Exponent k_I in  Robinson et al., 2020, Eq. 33
  REAL(dp)            :: pc_k_p_config                               = 0.2_dp                           ! Exponent k_p in  Robinson et al., 2020, Eq. 33
  REAL(dp)            :: pc_eta_min_config                           = 1E-8_dp                          ! Normalisation term in estimation of the truncation error (Robinson et al., Eq. 32)
  INTEGER             :: pc_max_timestep_iterations_config           = 5                                ! Maximum number of iterations of each time step
  REAL(dp)            :: pc_redo_tol_config                          = 10._dp                           ! Maximum allowed truncation error (any higher and the timestep is decreased)
  REAL(dp)            :: dt_min_config                               = 0.01_dp                          ! Smallest allowed time step [yr]
  
  ! Ice thickness boundary conditions
  CHARACTER(LEN=256)  :: ice_thickness_west_BC_config                = 'zero'                           ! Choice of boundary conditions for ice thickness at the domain boundary: "infinite", "periodic", "zero", "ISMIP_HOM_F"
  CHARACTER(LEN=256)  :: ice_thickness_east_BC_config                = 'zero' 
  CHARACTER(LEN=256)  :: ice_thickness_south_BC_config               = 'zero'
  CHARACTER(LEN=256)  :: ice_thickness_north_BC_config               = 'zero'
  CHARACTER(LEN=256)  :: choice_mask_noice_NAM_config                = 'NAM_remove_GRL'                 ! Choice of mask_noice configuration
  CHARACTER(LEN=256)  :: choice_mask_noice_EAS_config                = 'EAS_remove_GRL'
  CHARACTER(LEN=256)  :: choice_mask_noice_GRL_config                = 'GRL_remove_Ellesmere'
  CHARACTER(LEN=256)  :: choice_mask_noice_ANT_config                = 'none'                           ! For Antarctica, additional choices are included for certain idealised-geometry experiments: "MISMIP_mod", "MISMIP+"

  ! Ice dynamics - basal conditions and sliding
  ! ===========================================
  
  ! Sliding laws
  CHARACTER(LEN=256)  :: choice_sliding_law_config                   = 'Coulomb_regularised'            ! Choice of sliding law: "no_sliding", "idealised", "Coulomb", "Coulomb_regularised", "Weertman", "Tsai2015", "Schoof2005", "Zoet-Iverson"
  CHARACTER(LEN=256)  :: choice_idealised_sliding_law_config         = ''                               ! "ISMIP_HOM_C", "ISMIP_HOM_D", "ISMIP_HOM_E", "ISMIP_HOM_F"
  REAL(dp)            :: slid_delta_v_config                         = 1.0E-3_dp                        ! Normalisation parameter to prevent errors when velocity is zero
  REAL(dp)            :: slid_Weertman_m_config                      = 3._dp                            ! Exponent in Weertman sliding law
  REAL(dp)            :: slid_Coulomb_reg_q_plastic_config           = 0.3_dp                           ! Scaling exponent   in regularised Coulomb sliding law
  REAL(dp)            :: slid_Coulomb_reg_u_threshold_config         = 100._dp                          ! Threshold velocity in regularised Coulomb sliding law
  REAL(dp)            :: slid_ZI_ut_config                           = 200._dp                          ! (uniform) transition velocity used in the Zoet-Iverson sliding law [m/yr]
  REAL(dp)            :: slid_ZI_p_config                            = 5._dp                            ! Velocity exponent             used in the Zoet-Iverson sliding law
  
  ! Basal hydrology
  CHARACTER(LEN=256)  :: choice_basal_hydrology_config               = 'Martin2011'                     ! Choice of basal conditions: "saturated", "Martin2011"
  REAL(dp)            :: Martin2011_hydro_Hb_min_config              = 0._dp                            ! Martin et al. (2011) basal hydrology model: low-end  Hb  value of bedrock-dependent pore-water pressure
  REAL(dp)            :: Martin2011_hydro_Hb_max_config              = 1000._dp                         ! Martin et al. (2011) basal hydrology model: high-end Hb  value of bedrock-dependent pore-water pressure
  
  ! Basal roughness / friction
  CHARACTER(LEN=256)  :: choice_basal_roughness_config               = 'parameterised'                  ! "uniform", "parameterised", "prescribed"
  REAL(dp)            :: slid_Weertman_beta_sq_uniform_config        = 1.0E4_dp                         ! Uniform value for beta_sq  in Weertman sliding law
  REAL(dp)            :: slid_Coulomb_phi_fric_uniform_config        = 15._dp                           ! Uniform value for phi_fric in (regularised) Coulomb sliding law
  REAL(dp)            :: slid_Tsai2015_alpha_sq_uniform_config       = 0.5_dp                           ! Uniform value for alpha_sq in the Tsai2015 sliding law
  REAL(dp)            :: slid_Tsai2015_beta_sq_uniform_config        = 1.0E4_dp                         ! Uniform value for beta_sq  in the Tsai2015 sliding law
  REAL(dp)            :: slid_Schoof2005_alpha_sq_uniform_config     = 0.5_dp                           ! Uniform value for alpha_sq in the Schoof2005 sliding law
  REAL(dp)            :: slid_Schoof2005_beta_sq_uniform_config      = 1.0E4_dp                         ! Uniform value for beta_sq  in the Schoof2005 sliding law
  REAL(dp)            :: slid_ZI_phi_fric_uniform_config             = 15._dp                           ! Uniform value for phi_fric in the Zoet-Iverson sliding law
  CHARACTER(LEN=256)  :: choice_param_basal_roughness_config         = 'Martin2011'                     ! "Martin2011", "SSA_icestream", "MISMIP+", "BIVMIP_A", "BIVMIP_B", "BIVMIP_C"
  REAL(dp)            :: Martin2011till_phi_Hb_min_config            = -1000._dp                        ! Martin et al. (2011) bed roughness model: low-end  Hb  value of bedrock-dependent till friction angle
  REAL(dp)            :: Martin2011till_phi_Hb_max_config            = 0._dp                            ! Martin et al. (2011) bed roughness model: high-end Hb  value of bedrock-dependent till friction angle
  REAL(dp)            :: Martin2011till_phi_min_config               = 5._dp                            ! Martin et al. (2011) bed roughness model: low-end  phi value of bedrock-dependent till friction angle
  REAL(dp)            :: Martin2011till_phi_max_config               = 20._dp                           ! Martin et al. (2011) bed roughness model: high-end phi value of bedrock-dependent till friction angle
  CHARACTER(LEN=256)  :: basal_roughness_filename_config             = ''                               ! NetCDF file containing a basal roughness field for the chosen sliding law
  
  ! Ice dynamics - calving
  ! ======================
  
  CHARACTER(LEN=256)  :: choice_calving_law_config                   = 'threshold_thickness'            ! Choice of calving law: "none", "threshold_thickness"
  REAL(dp)            :: calving_threshold_thickness_config          = 200._dp                          ! Threshold ice thickness in the "threshold_thickness" calving law (200m taken from ANICE)
  LOGICAL             :: do_remove_shelves_config                    = .FALSE.                          ! If set to TRUE, all floating ice is always instantly removed (used in the ABUMIP-ABUK experiment)
  LOGICAL             :: remove_shelves_larger_than_PD_config        = .FALSE.                          ! If set to TRUE, all floating ice beyond the present-day calving front is removed (used for some Antarctic spin-ups)
  LOGICAL             :: continental_shelf_calving_config            = .FALSE.                          ! If set to TRUE, all ice beyond the continental shelf edge (set by a maximum depth) is removed
  REAL(dp)            :: continental_shelf_min_height_config         = -2000._dp                        ! Maximum depth of the continental shelf
  
  ! Thermodynamics and rheology
  ! ===========================
  
  CHARACTER(LEN=256)  :: choice_initial_ice_temperature_config       = 'Robin'                          ! Choice of initial ice temperature profile: "uniform", "linear", "Robin", "restart"
  REAL(dp)            :: uniform_ice_temperature_config              = 270._dp                          ! Uniform ice temperature (applied when choice_initial_ice_temperature_config = "uniform")
  CHARACTER(LEN=256)  :: choice_thermo_model_config                  = '3D_heat_equation'               ! Choice of thermodynamical model: "none", "3D_heat_equation"
  CHARACTER(LEN=256)  :: choice_ice_rheology_config                  = 'Huybrechts1992'                 ! Choice of ice rheology model: "uniform", "Huybrechts1992", "MISMIP_mod"
  REAL(dp)            :: uniform_flow_factor_config                  = 1E-16_dp                         ! Uniform ice flow factor (applied when choice_ice_rheology_model_config = "uniform")
  CHARACTER(LEN=256)  :: choice_ice_heat_capacity_config             = 'Pounder1965'                    ! Choice of ice heat capacity model: "uniform", "Pounder1965"
  REAL(dp)            :: uniform_ice_heat_capacity_config            = 2009._dp                         ! Uniform ice heat capacity (applied when choice_ice_heat_capacity_config = "uniform")
  CHARACTER(LEN=256)  :: choice_ice_thermal_conductivity_config      = 'Ritz1987'                       ! Choice of ice heat capacity model: "uniform", "Ritz1987"
  REAL(dp)            :: uniform_ice_thermal_conductivity_config     = 6.626958E7_dp                    ! Uniform ice thermal conductivity (applied when choice_ice_thermal_conductivity_config = "uniform")
  
  ! Sea level and GIA
  ! =================
  
  LOGICAL             :: do_ocean_floodfill_config                   = .TRUE.                           ! Use a flood-fill to determine the ocean mask, so that (pro-/sub-glacial) lakes dont exist
  CHARACTER(LEN=256)  :: choice_sealevel_model_config                = 'eustatic'                       ! Can be "fixed", "prescribed", "eustatic", or "SELEN"
  REAL(dp)            :: fixed_sealevel_config                       = 0._dp
  CHARACTER(LEN=256)  :: filename_sealevel_record_config             = 'name_of_file.dat'
  INTEGER             :: sealevel_record_length_config               = 1
  
  CHARACTER(LEN=256)  :: choice_GIA_model_config                     = 'ELRA'                           ! Can be "none", "ELRA", or "SELEN"
  REAL(dp)            :: ELRA_lithosphere_flex_rigidity_config       = 1.0E+25_dp                       ! Lithospheric flexural rigidity [kg m^2 s^-2]
  REAL(dp)            :: ELRA_bedrock_relaxation_time_config         = 3000.0_dp                        ! Relaxation time for bedrock adjustment [yr]
  REAL(dp)            :: ELRA_mantle_density_config                  = 3300.0_dp                        ! Mantle density [kg m^-3]
  
  ! Climate matrix
  ! ==============
  
  ! Present-day observed climate (ERA40) (NetCDF)
  CHARACTER(LEN=256)  :: filename_PD_obs_climate_config              = 'Datasets/ERA40/ERA40_climate_global.nc'
  
  ! GCM snapshots
  CHARACTER(LEN=256)  :: choice_climate_matrix_config                = 'PI_LGM'                         ! 'PI_LGM' uses 2 snapshots
  CHARACTER(LEN=256)  :: filename_GCM_snapshot_PI_config             = 'Datasets/GCM_snapshots/Singarayer_Valdes_2010_PI_Control.nc'
  CHARACTER(LEN=256)  :: filename_GCM_snapshot_LGM_config            = 'Datasets/GCM_snapshots/Singarayer_Valdes_2010_LGM.nc'
  
  ! Ice5G ice-sheet geometry
  CHARACTER(LEN=256)  :: filename_ICE5G_PD_config                    = 'Datasets/ICE5G/ice5g_v1.2_00.0k_1deg.nc'
  CHARACTER(LEN=256)  :: filename_ICE5G_LGM_config                   = 'Datasets/ICE5G/ice5g_v1.2_21.0k_1deg.nc'
  
  REAL(dp)            :: constant_lapserate_config                   = 0.008_dp                         ! Constant atmospheric lapse rate [K m^-1]
  
  ! Forcing
  ! =======
  
  ! The choice of forcing:
  ! 'd18O_inverse_dT_glob' : Use the inverse routine with the specified d18O record to calculate a global temperature offset (e.g. de Boer et al., 2013)
  ! 'CO2_direct'           : Use the specified CO2 record to force the climate matrix (e.g. Berends et al., 2018)
  ! 'd18O_inverse_CO2'     : Use the inverse routine with the specified d18O record to calculate CO2 and then force the climate matrix (e.g. Berends et al., 2019)
  CHARACTER(LEN=256)  :: choice_forcing_method_config                = 'CO2_direct'
  
  REAL(dp)            :: dT_deepwater_averaging_window_config        = 3000                             ! Time window (in yr) over which global mean temperature anomaly is averaged to find the deep-water temperature anomaly
  REAL(dp)            :: dT_deepwater_dT_surf_ratio_config           = 0.25_dp                          ! Ratio between global mean surface temperature change and deep-water temperature change
  REAL(dp)            :: d18O_dT_deepwater_ratio_config              = -0.28_dp                         ! Ratio between deep-water temperature change and benthic d18O change
  
  REAL(dp)            :: dT_glob_inverse_averaging_window_config     = 2000._dp                         ! Time window (in yr) over which global mean temperature anomaly is averaged before changing it with the inverse routine
  REAL(dp)            :: inverse_d18O_to_dT_glob_scaling_config      = 20._dp                           ! Scaling factor between modelled d18O anomaly and prescribed temperature anomaly change (value from de Boer et al., 2013)
  REAL(dp)            :: CO2_inverse_averaging_window_config         = 2000._dp                         ! Time window (in yr) over which CO2                             is averaged before changing it with the inverse routine
  REAL(dp)            :: inverse_d18O_to_CO2_scaling_config          = 68._dp                           ! Scaling factor between modelled d18O anomaly and modelled CO2 change (value from Berends et al., 2019)
  REAL(dp)            :: inverse_d18O_to_CO2_initial_CO2_config      = 280._dp                          ! CO2 value at the start of the simulation when using the inverse method to calculate CO2
  
  ! SMB tuning
  ! ==========
  
  REAL(dp)            :: C_abl_constant_NAM_config                   = -49._dp                          ! 34._dp    (commented values are old ANICE defaults, but since refreezing was not calculated right
  REAL(dp)            :: C_abl_constant_EAS_config                   = -49._dp                          !            and this has since been fixed, these values will still not give the same results as
  REAL(dp)            :: C_abl_constant_GRL_config                   = -49._dp                          !            they used to in ANICE.)
  REAL(dp)            :: C_abl_constant_ANT_config                   = -49._dp
  REAL(dp)            :: C_abl_Ts_NAM_config                         = 10._dp                           ! 10._dp
  REAL(dp)            :: C_abl_Ts_EAS_config                         = 10._dp
  REAL(dp)            :: C_abl_Ts_GRL_config                         = 10._dp
  REAL(dp)            :: C_abl_Ts_ANT_config                         = 10._dp
  REAL(dp)            :: C_abl_Q_NAM_config                          = 0.0227_dp                        ! 0.513_dp
  REAL(dp)            :: C_abl_Q_EAS_config                          = 0.0227_dp
  REAL(dp)            :: C_abl_Q_GRL_config                          = 0.0227_dp
  REAL(dp)            :: C_abl_Q_ANT_config                          = 0.0227_dp
  REAL(dp)            :: C_refr_NAM_config                           = 0.051_dp                         ! 0.012_dp
  REAL(dp)            :: C_refr_EAS_config                           = 0.051_dp 
  REAL(dp)            :: C_refr_GRL_config                           = 0.051_dp 
  REAL(dp)            :: C_refr_ANT_config                           = 0.051_dp
  
  ! Sub-shelf melt parameterisation
  ! ===============================
  
  ! Ocean temperature (used for both thermodynamics and basal melt)
  CHARACTER(LEN=256)  :: choice_ocean_temperature_model_config       = 'scaled'                         ! Can be "fixed" (use PD value) or "scaled" (scale between "PD", "warm", and "cold" values based on forcing (prescribed or inverse-modelled))
  REAL(dp)            :: ocean_temperature_PD_config                 = 271.46_dp                        ! present day temperature of the ocean beneath the shelves [K; -1.7 Celsius]
  REAL(dp)            :: ocean_temperature_cold_config               = 268.16_dp                        ! cold period temperature of the ocean beneath the shelves [K; -5.0 Celcius]
  REAL(dp)            :: ocean_temperature_warm_config               = 275.16_dp                        ! warm period temperature of the ocean beneath the shelves [K;  2.0 Celcius]
  
  ! Mean ocean temperatures and deep ocean & exposed shelf melt rates were tuned by Bas de Boer in 2011.
  REAL(dp)            :: T_ocean_mean_PD_NAM_config                  = -1.7_dp                          ! Present day temperature of the ocean beneath the shelves [Celcius]
  REAL(dp)            :: T_ocean_mean_PD_EAS_config                  = -1.7_dp
  REAL(dp)            :: T_ocean_mean_PD_GRL_config                  =  2.0_dp
  REAL(dp)            :: T_ocean_mean_PD_ANT_config                  = -1.7_dp
  REAL(dp)            :: T_ocean_mean_cold_NAM_config                = -5.0_dp                          ! Cold period temperature of the ocean beneath the shelves [Celcius]
  REAL(dp)            :: T_ocean_mean_cold_EAS_config                = -5.0_dp
  REAL(dp)            :: T_ocean_mean_cold_GRL_config                =  0.0_dp
  REAL(dp)            :: T_ocean_mean_cold_ANT_config                = -5.0_dp
  REAL(dp)            :: T_ocean_mean_warm_NAM_config                =  2.0_dp                          ! Warm period temperature of the ocean beneath the shelves [Celcius]
  REAL(dp)            :: T_ocean_mean_warm_EAS_config                =  2.0_dp
  REAL(dp)            :: T_ocean_mean_warm_GRL_config                =  4.0_dp
  REAL(dp)            :: T_ocean_mean_warm_ANT_config                =  2.0_dp
            
  REAL(dp)            :: BMB_deepocean_PD_NAM_config                 =  -5._dp                          ! Present-day sub-shelf melt rate for deep-ocean areas [m/year]
  REAL(dp)            :: BMB_deepocean_PD_EAS_config                 =  -5._dp
  REAL(dp)            :: BMB_deepocean_PD_GRL_config                 =  -5._dp
  REAL(dp)            :: BMB_deepocean_PD_ANT_config                 =  -5._dp
  REAL(dp)            :: BMB_deepocean_cold_NAM_config               =  -2._dp                          ! Cold period sub-shelf melt rate for deep-ocean areas [m/year]
  REAL(dp)            :: BMB_deepocean_cold_EAS_config               =  -2._dp
  REAL(dp)            :: BMB_deepocean_cold_GRL_config               =  -2._dp
  REAL(dp)            :: BMB_deepocean_cold_ANT_config               =  -2._dp
  REAL(dp)            :: BMB_deepocean_warm_NAM_config               = -10._dp                          ! Warm period sub-shelf melt rate for deep-ocean areas [m/year]    
  REAL(dp)            :: BMB_deepocean_warm_EAS_config               = -10._dp
  REAL(dp)            :: BMB_deepocean_warm_GRL_config               = -10._dp
  REAL(dp)            :: BMB_deepocean_warm_ANT_config               = -10._dp

  REAL(dp)            :: BMB_shelf_exposed_PD_NAM_config             =  -3._dp                          ! Present-day sub-shelf melt rate for exposed areas    [m/year]
  REAL(dp)            :: BMB_shelf_exposed_PD_EAS_config             =  -3._dp
  REAL(dp)            :: BMB_shelf_exposed_PD_GRL_config             =  -3._dp
  REAL(dp)            :: BMB_shelf_exposed_PD_ANT_config             =  -3._dp
  REAL(dp)            :: BMB_shelf_exposed_cold_NAM_config           =  -0._dp                          ! Cold period sub-shelf melt rate for exposed areas    [m/year]
  REAL(dp)            :: BMB_shelf_exposed_cold_EAS_config           =  -0._dp
  REAL(dp)            :: BMB_shelf_exposed_cold_GRL_config           =  -0._dp
  REAL(dp)            :: BMB_shelf_exposed_cold_ANT_config           =  -0._dp
  REAL(dp)            :: BMB_shelf_exposed_warm_NAM_config           =  -6._dp                          ! Warm period sub-shelf melt rate for exposed areas    [m/year]
  REAL(dp)            :: BMB_shelf_exposed_warm_EAS_config           =  -6._dp
  REAL(dp)            :: BMB_shelf_exposed_warm_GRL_config           =  -6._dp
  REAL(dp)            :: BMB_shelf_exposed_warm_ANT_config           =  -6._dp
    
  REAL(dp)            :: subshelf_melt_factor_NAM_config             = 0.005_dp                         ! Overall tuning factor for sub-shelf melt rate
  REAL(dp)            :: subshelf_melt_factor_EAS_config             = 0.005_dp
  REAL(dp)            :: subshelf_melt_factor_GRL_config             = 0.005_dp
  REAL(dp)            :: subshelf_melt_factor_ANT_config             = 0.005_dp
  
  REAL(dp)            :: deep_ocean_threshold_depth_NAM_config       = 1200._dp                         ! Threshold water depth for "deep ocean" (as opposed to continental shelf);
  REAL(dp)            :: deep_ocean_threshold_depth_EAS_config       = 800._dp                          ! this mostly prevents ice shelves from growing beyond the continental shelf
  REAL(dp)            :: deep_ocean_threshold_depth_GRL_config       = 800._dp                          ! Different depths for different regions is a bit ad hoc, but in reality
  REAL(dp)            :: deep_ocean_threshold_depth_ANT_config       = 1800._dp                         ! the different surface ocean temperatures probably result in the same effect...
  
  ! Which data fields will be written to the help_fields output file
  ! ================================================================
  
  CHARACTER(LEN=256)  :: help_field_01_config                        = 'lat'
  CHARACTER(LEN=256)  :: help_field_02_config                        = 'lon'
  CHARACTER(LEN=256)  :: help_field_03_config                        = 'resolution'
  CHARACTER(LEN=256)  :: help_field_04_config                        = 'U_surf'
  CHARACTER(LEN=256)  :: help_field_05_config                        = 'V_surf'
  CHARACTER(LEN=256)  :: help_field_06_config                        = 'U_base'
  CHARACTER(LEN=256)  :: help_field_07_config                        = 'V_base'
  CHARACTER(LEN=256)  :: help_field_08_config                        = 'U_vav'
  CHARACTER(LEN=256)  :: help_field_09_config                        = 'V_vav'
  CHARACTER(LEN=256)  :: help_field_10_config                        = 'U_3D'
  CHARACTER(LEN=256)  :: help_field_11_config                        = 'V_3D'
  CHARACTER(LEN=256)  :: help_field_12_config                        = 'W_3D'
  CHARACTER(LEN=256)  :: help_field_13_config                        = 'mask'
  CHARACTER(LEN=256)  :: help_field_14_config                        = 'dHs_dx'
  CHARACTER(LEN=256)  :: help_field_15_config                        = 'dHs_dy'
  CHARACTER(LEN=256)  :: help_field_16_config                        = 'T2m_year'
  CHARACTER(LEN=256)  :: help_field_17_config                        = 'Precip_year'
  CHARACTER(LEN=256)  :: help_field_18_config                        = 'Albedo_year'
  CHARACTER(LEN=256)  :: help_field_19_config                        = 'SMB_year'
  CHARACTER(LEN=256)  :: help_field_20_config                        = 'BMB'
  CHARACTER(LEN=256)  :: help_field_21_config                        = 'T2m'
  CHARACTER(LEN=256)  :: help_field_22_config                        = 'Precip'
  CHARACTER(LEN=256)  :: help_field_23_config                        = 'Albedo'
  CHARACTER(LEN=256)  :: help_field_24_config                        = 'SMB'
  CHARACTER(LEN=256)  :: help_field_25_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_26_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_27_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_28_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_29_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_30_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_31_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_32_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_33_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_34_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_35_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_36_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_37_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_38_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_39_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_40_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_41_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_42_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_43_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_44_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_45_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_46_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_47_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_48_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_49_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_50_config                        = 'none'



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
    REAL(dp)                            :: dt_ocean
    REAL(dp)                            :: dt_SMB
    REAL(dp)                            :: dt_BMB
    REAL(dp)                            :: dt_output
    REAL(dp)                            :: dt_mesh_min
    REAL(dp)                            :: dt_bedrock_ELRA
  
    ! Debugging
    ! =========
    
    LOGICAL                             :: do_write_debug_data
    LOGICAL                             :: do_check_for_NaN
    LOGICAL                             :: do_write_memory_tracker
    
    ! Which ice sheets do we simulate?
    ! ================================
    
    LOGICAL                             :: do_NAM
    LOGICAL                             :: do_EAS
    LOGICAL                             :: do_GRL
    LOGICAL                             :: do_ANT  

    ! Domain size for the four regions
    ! ================================

    REAL(dp)                            :: xmin_NAM
    REAL(dp)                            :: xmax_NAM
    REAL(dp)                            :: ymin_NAM
    REAL(dp)                            :: ymax_NAM

    REAL(dp)                            :: xmin_EAS
    REAL(dp)                            :: xmax_EAS
    REAL(dp)                            :: ymin_EAS
    REAL(dp)                            :: ymax_EAS

    REAL(dp)                            :: xmin_GRL
    REAL(dp)                            :: xmax_GRL
    REAL(dp)                            :: ymin_GRL
    REAL(dp)                            :: ymax_GRL

    REAL(dp)                            :: xmin_ANT
    REAL(dp)                            :: xmax_ANT
    REAL(dp)                            :: ymin_ANT
    REAL(dp)                            :: ymax_ANT
  
    ! Whether or not the simulation is a restart of a previous simulation
    ! ===================================================================
    
    LOGICAL                             :: is_restart
    REAL(dp)                            :: time_to_restart_from
    
    ! Initial model state when restarting from a previous run
    CHARACTER(LEN=256)                  :: filename_restart_NAM
    CHARACTER(LEN=256)                  :: filename_restart_EAS
    CHARACTER(LEN=256)                  :: filename_restart_GRL
    CHARACTER(LEN=256)                  :: filename_restart_ANT
    
    ! Benchmark experiments
    ! =====================
    
    LOGICAL                             :: do_benchmark_experiment
    CHARACTER(LEN=256)                  :: choice_benchmark_experiment
    REAL(dp)                            :: SSA_icestream_m
    REAL(dp)                            :: ISMIP_HOM_L
    CHARACTER(LEN=256)                  :: ISMIP_HOM_E_Arolla_filename
    LOGICAL                             :: MISMIPplus_do_tune_A_for_GL
    REAL(dp)                            :: MISMIPplus_xGL_target
    REAL(dp)                            :: MISMIPplus_A_flow_initial
    CHARACTER(LEN=256)                  :: MISMIPplus_scenario
    CHARACTER(LEN=256)                  :: MISOMIP1_scenario
    
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
   
    LOGICAL                             :: create_new_output_dir
    CHARACTER(LEN=256)                  :: output_dir

    ! Scaled vertical coordinate zeta  
    ! ===============================
     
    INTEGER                             :: nz       ! Number of grid points in vertical direction for thermodynamics in ice sheet
    REAL(dp), DIMENSION(:), ALLOCATABLE :: zeta

    ! Reference geometries (initial, present-day, and GIA equilibrium)
    ! ================================================================
    
    ! Initial geometry
    CHARACTER(LEN=256)                  :: choice_refgeo_init_NAM
    CHARACTER(LEN=256)                  :: choice_refgeo_init_EAS
    CHARACTER(LEN=256)                  :: choice_refgeo_init_GRL
    CHARACTER(LEN=256)                  :: choice_refgeo_init_ANT
    REAL(dp)                            :: time_to_restart_from_NAM
    REAL(dp)                            :: time_to_restart_from_EAS
    REAL(dp)                            :: time_to_restart_from_GRL
    REAL(dp)                            :: time_to_restart_from_ANT
    CHARACTER(LEN=256)                  :: choice_refgeo_init_idealised
    REAL(dp)                            :: dx_refgeo_init_idealised
    CHARACTER(LEN=256)                  :: filename_refgeo_init_NAM
    CHARACTER(LEN=256)                  :: filename_refgeo_init_EAS
    CHARACTER(LEN=256)                  :: filename_refgeo_init_GRL
    CHARACTER(LEN=256)                  :: filename_refgeo_init_ANT

    ! Present-day geometry
    CHARACTER(LEN=256)                  :: choice_refgeo_PD_NAM
    CHARACTER(LEN=256)                  :: choice_refgeo_PD_EAS
    CHARACTER(LEN=256)                  :: choice_refgeo_PD_GRL
    CHARACTER(LEN=256)                  :: choice_refgeo_PD_ANT
    CHARACTER(LEN=256)                  :: choice_refgeo_PD_idealised
    REAL(dp)                            :: dx_refgeo_PD_idealised
    CHARACTER(LEN=256)                  :: filename_refgeo_PD_NAM
    CHARACTER(LEN=256)                  :: filename_refgeo_PD_EAS
    CHARACTER(LEN=256)                  :: filename_refgeo_PD_GRL
    CHARACTER(LEN=256)                  :: filename_refgeo_PD_ANT

    ! GIA equilibrium geometry
    CHARACTER(LEN=256)                  :: choice_refgeo_GIAeq_NAM
    CHARACTER(LEN=256)                  :: choice_refgeo_GIAeq_EAS
    CHARACTER(LEN=256)                  :: choice_refgeo_GIAeq_GRL
    CHARACTER(LEN=256)                  :: choice_refgeo_GIAeq_ANT
    CHARACTER(LEN=256)                  :: choice_refgeo_GIAeq_idealised
    REAL(dp)                            :: dx_refgeo_GIAeq_idealised
    CHARACTER(LEN=256)                  :: filename_refgeo_GIAeq_NAM
    CHARACTER(LEN=256)                  :: filename_refgeo_GIAeq_EAS
    CHARACTER(LEN=256)                  :: filename_refgeo_GIAeq_GRL
    CHARACTER(LEN=256)                  :: filename_refgeo_GIAeq_ANT
    
    LOGICAL                             :: remove_Lake_Vostok

    ! Input data file paths
    ! =====================
    
    CHARACTER(LEN=256)                  :: filename_insolation
    
    CHARACTER(LEN=256)                  :: filename_CO2_record
    INTEGER                             :: CO2_record_length
    CHARACTER(LEN=256)                  :: filename_d18O_record
    INTEGER                             :: d18O_record_length
    
    ! Geothermal heat flux
    CHARACTER(LEN=256)                  :: choice_geothermal_heat_flux
    REAL(dp)                            :: constant_geothermal_heat_flux
    CHARACTER(LEN=256)                  :: filename_geothermal_heat_flux

    ! Ice dynamics - velocity
    ! =======================
    
    CHARACTER(LEN=256)                  :: choice_ice_dynamics
    REAL(dp)                            :: n_flow
    REAL(dp)                            :: m_enh_sheet
    REAL(dp)                            :: m_enh_shelf
    CHARACTER(LEN=256)                  :: choice_ice_margin
    LOGICAL                             :: include_SSADIVA_crossterms
    LOGICAL                             :: do_GL_subgrid_friction
    LOGICAL                             :: do_smooth_geometry
    REAL(dp)                            :: r_smooth_geometry
    
    ! Some parameters for numerically solving the SSA/DIVA
    REAL(dp)                            :: DIVA_visc_it_norm_dUV_tol
    INTEGER                             :: DIVA_visc_it_nit
    REAL(dp)                            :: DIVA_visc_it_relax
    REAL(dp)                            :: DIVA_beta_max
    REAL(dp)                            :: DIVA_err_lim
    REAL(dp)                            :: DIVA_vel_max
    REAL(dp)                            :: DIVA_vel_min
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_west
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_east
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_south
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_north
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_west
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_east
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_south
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_north
    CHARACTER(LEN=256)                  :: DIVA_choice_matrix_solver
    INTEGER                             :: DIVA_SOR_nit
    REAL(dp)                            :: DIVA_SOR_tol
    REAL(dp)                            :: DIVA_SOR_omega
    REAL(dp)                            :: DIVA_PETSc_rtol
    REAL(dp)                            :: DIVA_PETSc_abstol
  
    ! Ice dynamics - time integration
    ! ===============================
    
    CHARACTER(LEN=256)                  :: choice_timestepping
    CHARACTER(LEN=256)                  :: choice_ice_integration_method
    CHARACTER(LEN=256)                  :: dHi_choice_matrix_solver
    INTEGER                             :: dHi_SOR_nit
    REAL(dp)                            :: dHi_SOR_tol
    REAL(dp)                            :: dHi_SOR_omega
    REAL(dp)                            :: dHi_PETSc_rtol
    REAL(dp)                            :: dHi_PETSc_abstol
    
    ! Predictor-corrector ice-thickness update
    REAL(dp)                            :: pc_epsilon
    REAL(dp)                            :: pc_k_I
    REAL(dp)                            :: pc_k_p
    REAL(dp)                            :: pc_eta_min
    INTEGER                             :: pc_max_timestep_iterations
    REAL(dp)                            :: pc_redo_tol
    REAL(dp)                            :: dt_min
  
    ! Ice thickness boundary conditions
    CHARACTER(LEN=256)                  :: ice_thickness_west_BC
    CHARACTER(LEN=256)                  :: ice_thickness_east_BC
    CHARACTER(LEN=256)                  :: ice_thickness_south_BC
    CHARACTER(LEN=256)                  :: ice_thickness_north_BC
    CHARACTER(LEN=256)                  :: choice_mask_noice_NAM
    CHARACTER(LEN=256)                  :: choice_mask_noice_EAS
    CHARACTER(LEN=256)                  :: choice_mask_noice_GRL
    CHARACTER(LEN=256)                  :: choice_mask_noice_ANT
    
    ! Ice dynamics - basal conditions and sliding
    ! ===========================================

    ! Sliding laws
    CHARACTER(LEN=256)                  :: choice_sliding_law
    CHARACTER(LEN=256)                  :: choice_idealised_sliding_law
    REAL(dp)                            :: slid_delta_v
    REAL(dp)                            :: slid_Weertman_m
    REAL(dp)                            :: slid_Coulomb_reg_q_plastic
    REAL(dp)                            :: slid_Coulomb_reg_u_threshold
    REAL(dp)                            :: slid_ZI_ut
    REAL(dp)                            :: slid_ZI_p

    ! Basal hydrology
    CHARACTER(LEN=256)                  :: choice_basal_hydrology
    REAL(dp)                            :: Martin2011_hydro_Hb_min
    REAL(dp)                            :: Martin2011_hydro_Hb_max

    ! Basal roughness / friction
    CHARACTER(LEN=256)                  :: choice_basal_roughness
    REAL(dp)                            :: slid_Weertman_beta_sq_uniform
    REAL(dp)                            :: slid_Coulomb_phi_fric_uniform
    REAL(dp)                            :: slid_Tsai2015_alpha_sq_uniform
    REAL(dp)                            :: slid_Tsai2015_beta_sq_uniform
    REAL(dp)                            :: slid_Schoof2005_alpha_sq_uniform
    REAL(dp)                            :: slid_Schoof2005_beta_sq_uniform
    REAL(dp)                            :: slid_ZI_phi_fric_uniform
    CHARACTER(LEN=256)                  :: choice_param_basal_roughness
    REAL(dp)                            :: Martin2011till_phi_Hb_min
    REAL(dp)                            :: Martin2011till_phi_Hb_max
    REAL(dp)                            :: Martin2011till_phi_min
    REAL(dp)                            :: Martin2011till_phi_max
    CHARACTER(LEN=256)                  :: basal_roughness_filename
    
    ! Ice dynamics - calving
    ! ======================
    
    CHARACTER(LEN=256)                  :: choice_calving_law
    REAL(dp)                            :: calving_threshold_thickness
    LOGICAL                             :: do_remove_shelves
    LOGICAL                             :: remove_shelves_larger_than_PD
    LOGICAL                             :: continental_shelf_calving
    REAL(dp)                            :: continental_shelf_min_height

    ! Thermodynamics and rheology
    ! ===========================
    
    CHARACTER(LEN=256)                  :: choice_initial_ice_temperature
    REAL(dp)                            :: uniform_ice_temperature
    CHARACTER(LEN=256)                  :: choice_thermo_model
    CHARACTER(LEN=256)                  :: choice_ice_rheology
    REAL(dp)                            :: uniform_flow_factor
    CHARACTER(LEN=256)                  :: choice_ice_heat_capacity
    REAL(dp)                            :: uniform_ice_heat_capacity
    CHARACTER(LEN=256)                  :: choice_ice_thermal_conductivity
    REAL(dp)                            :: uniform_ice_thermal_conductivity
  
    ! Sea level and GIA
    ! =================
    
    LOGICAL                             :: do_ocean_floodfill
    CHARACTER(LEN=256)                  :: choice_sealevel_model
    REAL(dp)                            :: fixed_sealevel
    CHARACTER(LEN=256)                  :: filename_sealevel_record
    INTEGER                             :: sealevel_record_length
  
    CHARACTER(LEN=256)                  :: choice_GIA_model
    REAL(dp)                            :: ELRA_lithosphere_flex_rigidity
    REAL(dp)                            :: ELRA_bedrock_relaxation_time
    REAL(dp)                            :: ELRA_mantle_density
    
    ! Climate matrix
    ! ==============
    
    CHARACTER(LEN=256)                  :: filename_PD_obs_climate
    CHARACTER(LEN=256)                  :: choice_climate_matrix
    CHARACTER(LEN=256)                  :: filename_GCM_snapshot_PI
    CHARACTER(LEN=256)                  :: filename_GCM_snapshot_LGM
    CHARACTER(LEN=256)                  :: filename_ICE5G_PD
    CHARACTER(LEN=256)                  :: filename_ICE5G_LGM
    
    CHARACTER(LEN=256)                  :: choice_ocean_temperature_model
    REAL(dp)                            :: ocean_temperature_PD
    REAL(dp)                            :: ocean_temperature_cold
    REAL(dp)                            :: ocean_temperature_warm
    
    REAL(dp)                            :: constant_lapserate
    
    ! Forcing
    ! =======
    
    CHARACTER(LEN=256)                  :: choice_forcing_method
    
    REAL(dp)                            :: dT_deepwater_averaging_window
    REAL(dp)                            :: dT_deepwater_dT_surf_ratio
    REAL(dp)                            :: d18O_dT_deepwater_ratio
    
    REAL(dp)                            :: dT_glob_inverse_averaging_window
    REAL(dp)                            :: inverse_d18O_to_dT_glob_scaling
    REAL(dp)                            :: CO2_inverse_averaging_window
    REAL(dp)                            :: inverse_d18O_to_CO2_scaling
    REAL(dp)                            :: inverse_d18O_to_CO2_initial_CO2
    
    ! SMB melt tuning
    ! ===============
    
    REAL(dp)                            :: C_abl_constant_NAM
    REAL(dp)                            :: C_abl_constant_EAS
    REAL(dp)                            :: C_abl_constant_GRL
    REAL(dp)                            :: C_abl_constant_ANT
    REAL(dp)                            :: C_abl_Ts_NAM
    REAL(dp)                            :: C_abl_Ts_EAS
    REAL(dp)                            :: C_abl_Ts_GRL
    REAL(dp)                            :: C_abl_Ts_ANT
    REAL(dp)                            :: C_abl_Q_NAM
    REAL(dp)                            :: C_abl_Q_EAS
    REAL(dp)                            :: C_abl_Q_GRL
    REAL(dp)                            :: C_abl_Q_ANT
    REAL(dp)                            :: C_refr_NAM
    REAL(dp)                            :: C_refr_EAS
    REAL(dp)                            :: C_refr_GRL
    REAL(dp)                            :: C_refr_ANT
  
    ! Sub-shelf melt parameterisation
    ! ===============================
    
    REAL(dp)                            :: T_ocean_mean_PD_NAM
    REAL(dp)                            :: T_ocean_mean_PD_EAS
    REAL(dp)                            :: T_ocean_mean_PD_GRL
    REAL(dp)                            :: T_ocean_mean_PD_ANT
    REAL(dp)                            :: T_ocean_mean_cold_NAM
    REAL(dp)                            :: T_ocean_mean_cold_EAS
    REAL(dp)                            :: T_ocean_mean_cold_GRL
    REAL(dp)                            :: T_ocean_mean_cold_ANT
    REAL(dp)                            :: T_ocean_mean_warm_NAM
    REAL(dp)                            :: T_ocean_mean_warm_EAS
    REAL(dp)                            :: T_ocean_mean_warm_GRL
    REAL(dp)                            :: T_ocean_mean_warm_ANT
              
    REAL(dp)                            :: BMB_deepocean_PD_NAM
    REAL(dp)                            :: BMB_deepocean_PD_EAS
    REAL(dp)                            :: BMB_deepocean_PD_GRL
    REAL(dp)                            :: BMB_deepocean_PD_ANT
    REAL(dp)                            :: BMB_deepocean_cold_NAM
    REAL(dp)                            :: BMB_deepocean_cold_EAS
    REAL(dp)                            :: BMB_deepocean_cold_GRL
    REAL(dp)                            :: BMB_deepocean_cold_ANT
    REAL(dp)                            :: BMB_deepocean_warm_NAM
    REAL(dp)                            :: BMB_deepocean_warm_EAS
    REAL(dp)                            :: BMB_deepocean_warm_GRL
    REAL(dp)                            :: BMB_deepocean_warm_ANT
  
    REAL(dp)                            :: BMB_shelf_exposed_PD_NAM
    REAL(dp)                            :: BMB_shelf_exposed_PD_EAS
    REAL(dp)                            :: BMB_shelf_exposed_PD_GRL
    REAL(dp)                            :: BMB_shelf_exposed_PD_ANT
    REAL(dp)                            :: BMB_shelf_exposed_cold_NAM
    REAL(dp)                            :: BMB_shelf_exposed_cold_EAS
    REAL(dp)                            :: BMB_shelf_exposed_cold_GRL
    REAL(dp)                            :: BMB_shelf_exposed_cold_ANT
    REAL(dp)                            :: BMB_shelf_exposed_warm_NAM
    REAL(dp)                            :: BMB_shelf_exposed_warm_EAS
    REAL(dp)                            :: BMB_shelf_exposed_warm_GRL
    REAL(dp)                            :: BMB_shelf_exposed_warm_ANT
      
    REAL(dp)                            :: subshelf_melt_factor_NAM
    REAL(dp)                            :: subshelf_melt_factor_EAS
    REAL(dp)                            :: subshelf_melt_factor_GRL
    REAL(dp)                            :: subshelf_melt_factor_ANT
    
    REAL(dp)                            :: deep_ocean_threshold_depth_NAM
    REAL(dp)                            :: deep_ocean_threshold_depth_EAS
    REAL(dp)                            :: deep_ocean_threshold_depth_GRL
    REAL(dp)                            :: deep_ocean_threshold_depth_ANT
  
    ! Which data fields will be written to the help_fields output file
    ! ================================================================
    
    CHARACTER(LEN=256)                  :: help_field_01  
    CHARACTER(LEN=256)                  :: help_field_02  
    CHARACTER(LEN=256)                  :: help_field_03  
    CHARACTER(LEN=256)                  :: help_field_04  
    CHARACTER(LEN=256)                  :: help_field_05  
    CHARACTER(LEN=256)                  :: help_field_06  
    CHARACTER(LEN=256)                  :: help_field_07  
    CHARACTER(LEN=256)                  :: help_field_08  
    CHARACTER(LEN=256)                  :: help_field_09
    CHARACTER(LEN=256)                  :: help_field_10  
    CHARACTER(LEN=256)                  :: help_field_11  
    CHARACTER(LEN=256)                  :: help_field_12  
    CHARACTER(LEN=256)                  :: help_field_13  
    CHARACTER(LEN=256)                  :: help_field_14  
    CHARACTER(LEN=256)                  :: help_field_15  
    CHARACTER(LEN=256)                  :: help_field_16  
    CHARACTER(LEN=256)                  :: help_field_17  
    CHARACTER(LEN=256)                  :: help_field_18  
    CHARACTER(LEN=256)                  :: help_field_19
    CHARACTER(LEN=256)                  :: help_field_20  
    CHARACTER(LEN=256)                  :: help_field_21  
    CHARACTER(LEN=256)                  :: help_field_22  
    CHARACTER(LEN=256)                  :: help_field_23  
    CHARACTER(LEN=256)                  :: help_field_24  
    CHARACTER(LEN=256)                  :: help_field_25  
    CHARACTER(LEN=256)                  :: help_field_26  
    CHARACTER(LEN=256)                  :: help_field_27  
    CHARACTER(LEN=256)                  :: help_field_28  
    CHARACTER(LEN=256)                  :: help_field_29
    CHARACTER(LEN=256)                  :: help_field_30  
    CHARACTER(LEN=256)                  :: help_field_31  
    CHARACTER(LEN=256)                  :: help_field_32  
    CHARACTER(LEN=256)                  :: help_field_33  
    CHARACTER(LEN=256)                  :: help_field_34  
    CHARACTER(LEN=256)                  :: help_field_35  
    CHARACTER(LEN=256)                  :: help_field_36  
    CHARACTER(LEN=256)                  :: help_field_37  
    CHARACTER(LEN=256)                  :: help_field_38  
    CHARACTER(LEN=256)                  :: help_field_39
    CHARACTER(LEN=256)                  :: help_field_40  
    CHARACTER(LEN=256)                  :: help_field_41  
    CHARACTER(LEN=256)                  :: help_field_42  
    CHARACTER(LEN=256)                  :: help_field_43  
    CHARACTER(LEN=256)                  :: help_field_44  
    CHARACTER(LEN=256)                  :: help_field_45  
    CHARACTER(LEN=256)                  :: help_field_46  
    CHARACTER(LEN=256)                  :: help_field_47  
    CHARACTER(LEN=256)                  :: help_field_48  
    CHARACTER(LEN=256)                  :: help_field_49
    CHARACTER(LEN=256)                  :: help_field_50
    
    ! Values to be filled into the total mask (used only for diagnostic output)
    ! ==========================================================================
    
    INTEGER                             :: type_land
    INTEGER                             :: type_ocean
    INTEGER                             :: type_lake
    INTEGER                             :: type_sheet
    INTEGER                             :: type_shelf
    INTEGER                             :: type_coast
    INTEGER                             :: type_margin
    INTEGER                             :: type_groundingline
    INTEGER                             :: type_calvingfront
    
   ! Parameters of the polar stereographic projections of the four model regions
   ! (These have to match the values used to create the input files!)
   ! ===========================================================================     
                          
    REAL(dp)                            :: lambda_M_NAM                           
    REAL(dp)                            :: lambda_M_EAS                           
    REAL(dp)                            :: lambda_M_GRL                           
    REAL(dp)                            :: lambda_M_ANT
    REAL(dp)                            :: phi_M_NAM
    REAL(dp)                            :: phi_M_EAS
    REAL(dp)                            :: phi_M_GRL
    REAL(dp)                            :: phi_M_ANT
    REAL(dp)                            :: alpha_stereo_NAM
    REAL(dp)                            :: alpha_stereo_EAS
    REAL(dp)                            :: alpha_stereo_GRL
    REAL(dp)                            :: alpha_stereo_ANT

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
    NAMELIST /CONFIG/start_time_of_run_config,                        &
                     end_time_of_run_config,                          &
                     dt_coupling_config,                              &
                     dt_max_config,                                   &
                     dt_thermo_config,                                &
                     dt_climate_config,                               &
                     dt_ocean_config,                                 &
                     dt_SMB_config,                                   &
                     dt_BMB_config,                                   &
                     dt_output_config,                                &
                     dt_mesh_min_config,                              &
                     dt_bedrock_ELRA_config,                          &
                     do_write_debug_data_config,                      &
                     do_check_for_NaN_config,                         &
                     do_write_memory_tracker_config,                  &
                     do_NAM_config,                                   &
                     do_EAS_config,                                   &
                     do_GRL_config,                                   &
                     do_ANT_config,                                   & 
                     xmin_NAM_config,                                 &
                     xmax_NAM_config,                                 &
                     ymin_NAM_config,                                 &
                     ymax_NAM_config,                                 &
                     xmin_EAS_config,                                 &
                     xmax_EAS_config,                                 &
                     ymin_EAS_config,                                 &
                     ymax_EAS_config,                                 &
                     xmin_GRL_config,                                 &
                     xmax_GRL_config,                                 &
                     ymin_GRL_config,                                 &
                     ymax_GRL_config,                                 &
                     xmin_ANT_config,                                 &
                     xmax_ANT_config,                                 &
                     ymin_ANT_config,                                 &
                     ymax_ANT_config,                                 &
                     is_restart_config,                               &
                     time_to_restart_from_config,                     &
                     filename_restart_NAM_config,                     &
                     filename_restart_EAS_config,                     &
                     filename_restart_GRL_config,                     &
                     filename_restart_ANT_config,                     &
                     do_benchmark_experiment_config,                  &
                     choice_benchmark_experiment_config,              &
                     SSA_icestream_m_config,                          &
                     ISMIP_HOM_L_config,                              &
                     ISMIP_HOM_E_Arolla_filename_config,              &
                     MISMIPplus_do_tune_A_for_GL_config,              &
                     MISMIPplus_xGL_target_config,                    &
                     MISMIPplus_A_flow_initial_config,                &
                     MISMIPplus_scenario_config,                      &
                     MISOMIP1_scenario_config,                        &
                     nconmax_config,                                  &
                     alpha_min_config,                                &
                     dz_max_ice_config,                               &
                     res_max_config,                                  &
                     res_max_margin_config,                           &
                     res_max_gl_config,                               &
                     res_max_cf_config,                               &
                     res_max_mountain_config,                         &
                     res_max_coast_config,                            &
                     mesh_fitness_threshold_config,                   &
                     dx_grid_output_config,                           &
                     dx_grid_GIA_config,                              &
                     dx_grid_smooth_config,                           &
                     nPOI_NAM_config,                                 &
                     nPOI_EAS_config,                                 &
                     nPOI_GRL_config,                                 &
                     nPOI_ANT_config,                                 &
                     POI_NAM_coordinates_config,                      &
                     POI_EAS_coordinates_config,                      &
                     POI_GRL_coordinates_config,                      &
                     POI_ANT_coordinates_config,                      &
                     POI_NAM_resolutions_config,                      &
                     POI_EAS_resolutions_config,                      &
                     POI_GRL_resolutions_config,                      &
                     POI_ANT_resolutions_config,                      &
                     create_new_output_dir_config,                    &
                     output_dir_config,                               &
                     nz_config,                                       &
                     zeta_config,                                     &
                     choice_refgeo_init_NAM_config,                   &
                     choice_refgeo_init_EAS_config,                   &
                     choice_refgeo_init_GRL_config,                   &
                     choice_refgeo_init_ANT_config,                   &
                     time_to_restart_from_NAM_config,                 &
                     time_to_restart_from_EAS_config,                 &
                     time_to_restart_from_GRL_config,                 &
                     time_to_restart_from_ANT_config,                 &
                     choice_refgeo_init_idealised_config,             &
                     dx_refgeo_init_idealised_config,                 &
                     filename_refgeo_init_NAM_config,                 &
                     filename_refgeo_init_EAS_config,                 &
                     filename_refgeo_init_GRL_config,                 &
                     filename_refgeo_init_ANT_config,                 &
                     choice_refgeo_PD_NAM_config,                     &
                     choice_refgeo_PD_EAS_config,                     &
                     choice_refgeo_PD_GRL_config,                     &
                     choice_refgeo_PD_ANT_config,                     &
                     choice_refgeo_PD_idealised_config,               &
                     dx_refgeo_PD_idealised_config,                   &
                     filename_refgeo_PD_NAM_config,                   &
                     filename_refgeo_PD_EAS_config,                   &
                     filename_refgeo_PD_GRL_config,                   &
                     filename_refgeo_PD_ANT_config,                   &
                     choice_refgeo_GIAeq_NAM_config,                  &
                     choice_refgeo_GIAeq_EAS_config,                  &
                     choice_refgeo_GIAeq_GRL_config,                  &
                     choice_refgeo_GIAeq_ANT_config,                  &
                     choice_refgeo_GIAeq_idealised_config,            &
                     dx_refgeo_GIAeq_idealised_config,                &
                     filename_refgeo_GIAeq_NAM_config,                &
                     filename_refgeo_GIAeq_EAS_config,                &
                     filename_refgeo_GIAeq_GRL_config,                &
                     filename_refgeo_GIAeq_ANT_config,                &
                     remove_Lake_Vostok_config,                       &
                     filename_insolation_config,                      &
                     filename_CO2_record_config,                      &
                     CO2_record_length_config,                        &
                     filename_d18O_record_config,                     &
                     d18O_record_length_config,                       &
                     choice_geothermal_heat_flux_config,              &
                     constant_geothermal_heat_flux_config,            &
                     filename_geothermal_heat_flux_config,            &
                     choice_ice_dynamics_config,                      &
                     n_flow_config,                                   &
                     m_enh_sheet_config,                              &
                     m_enh_shelf_config,                              &
                     choice_ice_margin_config,                        &
                     include_SSADIVA_crossterms_config,               &
                     do_GL_subgrid_friction_config,                   &
                     do_smooth_geometry_config,                       &
                     r_smooth_geometry_config,                        &
                     DIVA_visc_it_norm_dUV_tol_config,                &
                     DIVA_visc_it_nit_config,                         &
                     DIVA_visc_it_relax_config,                       &
                     DIVA_beta_max_config,                            &
                     DIVA_err_lim_config,                             &
                     DIVA_vel_max_config,                             &
                     DIVA_vel_min_config,                             &
                     DIVA_boundary_BC_u_west_config,                  &
                     DIVA_boundary_BC_u_east_config,                  &
                     DIVA_boundary_BC_u_south_config,                 &
                     DIVA_boundary_BC_u_north_config,                 &
                     DIVA_boundary_BC_v_west_config,                  &
                     DIVA_boundary_BC_v_east_config,                  &
                     DIVA_boundary_BC_v_south_config,                 &
                     DIVA_boundary_BC_v_north_config,                 &
                     DIVA_choice_matrix_solver_config,                &
                     DIVA_SOR_nit_config,                             &
                     DIVA_SOR_tol_config,                             &
                     DIVA_SOR_omega_config,                           &
                     DIVA_PETSc_rtol_config,                          &
                     DIVA_PETSc_abstol_config,                        &
                     choice_timestepping_config,                      &
                     choice_ice_integration_method_config,            &
                     dHi_choice_matrix_solver_config,                 &
                     dHi_SOR_nit_config,                              &
                     dHi_SOR_tol_config,                              &
                     dHi_SOR_omega_config,                            &
                     dHi_PETSc_rtol_config,                           &
                     dHi_PETSc_abstol_config,                         &
                     pc_epsilon_config,                               &
                     pc_k_I_config,                                   &
                     pc_k_p_config,                                   &
                     pc_eta_min_config,                               &
                     pc_max_timestep_iterations_config,               &
                     pc_redo_tol_config,                              &
                     dt_min_config,                                   &
                     ice_thickness_west_BC_config,                    &
                     ice_thickness_east_BC_config,                    &
                     ice_thickness_south_BC_config,                   &
                     ice_thickness_north_BC_config,                   &
                     choice_mask_noice_NAM_config,                    &
                     choice_mask_noice_EAS_config,                    &
                     choice_mask_noice_GRL_config,                    &
                     choice_mask_noice_ANT_config,                    &
                     choice_sliding_law_config,                       &
                     choice_idealised_sliding_law_config,             &
                     slid_delta_v_config,                             &
                     slid_Weertman_m_config,                          &
                     slid_Coulomb_reg_q_plastic_config,               &
                     slid_Coulomb_reg_u_threshold_config,             &
                     slid_ZI_ut_config,                               &
                     slid_ZI_p_config,                                &
                     choice_basal_hydrology_config,                   &
                     Martin2011_hydro_Hb_min_config,                  &
                     Martin2011_hydro_Hb_max_config,                  &
                     choice_basal_roughness_config,                   &
                     slid_Weertman_beta_sq_uniform_config,            &
                     slid_Coulomb_phi_fric_uniform_config,            &
                     slid_Tsai2015_alpha_sq_uniform_config,           &
                     slid_Tsai2015_beta_sq_uniform_config,            &
                     slid_Schoof2005_alpha_sq_uniform_config,         &
                     slid_Schoof2005_beta_sq_uniform_config,          &
                     slid_ZI_phi_fric_uniform_config,                 &
                     choice_param_basal_roughness_config,             &
                     Martin2011till_phi_Hb_min_config,                &
                     Martin2011till_phi_Hb_max_config,                &
                     Martin2011till_phi_min_config,                   &
                     Martin2011till_phi_max_config,                   &
                     basal_roughness_filename_config,                 &
                     choice_calving_law_config,                       &
                     calving_threshold_thickness_config,              &
                     do_remove_shelves_config,                        &
                     remove_shelves_larger_than_PD_config,            &
                     continental_shelf_calving_config,                &
                     continental_shelf_min_height_config,             &
                     choice_initial_ice_temperature_config,           &
                     uniform_ice_temperature_config,                  &
                     choice_thermo_model_config,                      &
                     choice_ice_rheology_config,                      &
                     uniform_flow_factor_config,                      &
                     choice_ice_heat_capacity_config,                 &
                     uniform_ice_heat_capacity_config,                &
                     choice_ice_thermal_conductivity_config,          &
                     uniform_ice_thermal_conductivity_config,         &
                     do_ocean_floodfill_config,                       &
                     choice_sealevel_model_config,                    &
                     fixed_sealevel_config,                           &
                     filename_sealevel_record_config,                 &
                     sealevel_record_length_config,                   &
                     choice_GIA_model_config,                         &
                     ELRA_lithosphere_flex_rigidity_config,           &
                     ELRA_bedrock_relaxation_time_config,             &
                     ELRA_mantle_density_config,                      &
                     filename_PD_obs_climate_config,                  &
                     choice_climate_matrix_config,                    &
                     filename_GCM_snapshot_PI_config,                 &
                     filename_GCM_snapshot_LGM_config,                &
                     filename_ICE5G_PD_config,                        &
                     filename_ICE5G_LGM_config,                       &
                     choice_ocean_temperature_model_config,           &
                     ocean_temperature_PD_config,                     &
                     ocean_temperature_cold_config,                   &
                     ocean_temperature_warm_config,                   &
                     constant_lapserate_config,                       &
                     choice_forcing_method_config,                    &
                     dT_deepwater_averaging_window_config,            &
                     dT_deepwater_dT_surf_ratio_config,               &
                     d18O_dT_deepwater_ratio_config,                  &
                     dT_glob_inverse_averaging_window_config,         &
                     inverse_d18O_to_dT_glob_scaling_config,          &
                     CO2_inverse_averaging_window_config,             &
                     inverse_d18O_to_CO2_scaling_config,              &
                     inverse_d18O_to_CO2_initial_CO2_config,          &
                     C_abl_constant_NAM_config,                       &
                     C_abl_constant_EAS_config,                       &
                     C_abl_constant_GRL_config,                       &
                     C_abl_constant_ANT_config,                       &
                     C_abl_Ts_NAM_config,                             &
                     C_abl_Ts_EAS_config,                             &
                     C_abl_Ts_GRL_config,                             &
                     C_abl_Ts_ANT_config,                             &
                     C_abl_Q_NAM_config,                              &
                     C_abl_Q_EAS_config,                              &
                     C_abl_Q_GRL_config,                              &
                     C_abl_Q_ANT_config,                              &
                     C_refr_NAM_config,                               &
                     C_refr_EAS_config,                               &
                     C_refr_GRL_config,                               &
                     C_refr_ANT_config,                               &
                     T_ocean_mean_PD_NAM_config,                      &
                     T_ocean_mean_PD_EAS_config,                      &
                     T_ocean_mean_PD_GRL_config,                      &
                     T_ocean_mean_PD_ANT_config,                      &
                     T_ocean_mean_cold_NAM_config,                    &
                     T_ocean_mean_cold_EAS_config,                    &
                     T_ocean_mean_cold_GRL_config,                    &
                     T_ocean_mean_cold_ANT_config,                    &
                     T_ocean_mean_warm_NAM_config,                    &
                     T_ocean_mean_warm_EAS_config,                    &
                     T_ocean_mean_warm_GRL_config,                    &
                     T_ocean_mean_warm_ANT_config,                    &
                     BMB_deepocean_PD_NAM_config,                     &
                     BMB_deepocean_PD_EAS_config,                     &
                     BMB_deepocean_PD_GRL_config,                     &
                     BMB_deepocean_PD_ANT_config,                     &
                     BMB_deepocean_cold_NAM_config,                   &
                     BMB_deepocean_cold_EAS_config,                   &
                     BMB_deepocean_cold_GRL_config,                   &
                     BMB_deepocean_cold_ANT_config,                   &
                     BMB_deepocean_warm_NAM_config,                   &
                     BMB_deepocean_warm_EAS_config,                   &
                     BMB_deepocean_warm_GRL_config,                   &
                     BMB_deepocean_warm_ANT_config,                   &
                     BMB_shelf_exposed_PD_NAM_config,                 &
                     BMB_shelf_exposed_PD_EAS_config,                 &
                     BMB_shelf_exposed_PD_GRL_config,                 &
                     BMB_shelf_exposed_PD_ANT_config,                 &
                     BMB_shelf_exposed_cold_NAM_config,               &
                     BMB_shelf_exposed_cold_EAS_config,               &
                     BMB_shelf_exposed_cold_GRL_config,               &
                     BMB_shelf_exposed_cold_ANT_config,               &
                     BMB_shelf_exposed_warm_NAM_config,               &
                     BMB_shelf_exposed_warm_EAS_config,               &
                     BMB_shelf_exposed_warm_GRL_config,               &
                     BMB_shelf_exposed_warm_ANT_config,               &
                     subshelf_melt_factor_NAM_config,                 &
                     subshelf_melt_factor_EAS_config,                 &
                     subshelf_melt_factor_GRL_config,                 &
                     subshelf_melt_factor_ANT_config,                 &
                     deep_ocean_threshold_depth_NAM_config,           &
                     deep_ocean_threshold_depth_EAS_config,           &
                     deep_ocean_threshold_depth_GRL_config,           &
                     deep_ocean_threshold_depth_ANT_config,           &
                     help_field_01_config,                            &
                     help_field_02_config,                            &
                     help_field_03_config,                            &
                     help_field_04_config,                            &
                     help_field_05_config,                            &
                     help_field_06_config,                            &
                     help_field_07_config,                            &
                     help_field_08_config,                            &
                     help_field_09_config,                            &
                     help_field_10_config,                            &
                     help_field_11_config,                            &
                     help_field_12_config,                            &
                     help_field_13_config,                            &
                     help_field_14_config,                            &
                     help_field_15_config,                            &
                     help_field_16_config,                            &
                     help_field_17_config,                            &
                     help_field_18_config,                            &
                     help_field_19_config,                            &
                     help_field_20_config,                            &
                     help_field_21_config,                            &
                     help_field_22_config,                            &
                     help_field_23_config,                            &
                     help_field_24_config,                            &
                     help_field_25_config,                            &
                     help_field_26_config,                            &
                     help_field_27_config,                            &
                     help_field_28_config,                            &
                     help_field_29_config,                            &
                     help_field_30_config,                            &
                     help_field_31_config,                            &
                     help_field_32_config,                            &
                     help_field_33_config,                            &
                     help_field_34_config,                            &
                     help_field_35_config,                            &
                     help_field_36_config,                            &
                     help_field_37_config,                            &
                     help_field_38_config,                            &
                     help_field_39_config,                            &
                     help_field_40_config,                            &
                     help_field_41_config,                            &
                     help_field_42_config,                            &
                     help_field_43_config,                            &
                     help_field_44_config,                            &
                     help_field_45_config,                            &
                     help_field_46_config,                            &
                     help_field_47_config,                            &
                     help_field_48_config,                            &
                     help_field_49_config,                            &
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
    
    ! Time steps and range
    ! ====================
   
    C%start_time_of_run                        = start_time_of_run_config
    C%end_time_of_run                          = end_time_of_run_config
    C%dt_coupling                              = dt_coupling_config
    C%dt_max                                   = dt_max_config
    C%dt_thermo                                = dt_thermo_config
    C%dt_climate                               = dt_climate_config
    C%dt_ocean                                 = dt_ocean_config
    C%dt_SMB                                   = dt_SMB_config
    C%dt_BMB                                   = dt_BMB_config
    C%dt_output                                = dt_output_config
    C%dt_mesh_min                              = dt_mesh_min_config
    C%dt_bedrock_ELRA                          = dt_bedrock_ELRA_config
  
    ! Debugging
    ! =========
    
    C%do_write_debug_data                      = do_write_debug_data_config
    C%do_check_for_NaN                         = do_check_for_NaN_config
    C%do_write_memory_tracker                  = do_write_memory_tracker_config
    
    ! Which ice sheets do we simulate?
    ! ================================
    
    C%do_NAM                                   = do_NAM_config
    C%do_EAS                                   = do_EAS_config
    C%do_GRL                                   = do_GRL_config
    C%do_ANT                                   = do_ANT_config 

    ! Domain size for the four regions
    ! ================================

    C%xmin_NAM                                 = xmin_NAM_config
    C%xmax_NAM                                 = xmax_NAM_config
    C%ymin_NAM                                 = ymin_NAM_config
    C%ymax_NAM                                 = ymax_NAM_config

    C%xmin_EAS                                 = xmin_EAS_config
    C%xmax_EAS                                 = xmax_EAS_config
    C%ymin_EAS                                 = ymin_EAS_config
    C%ymax_EAS                                 = ymax_EAS_config

    C%xmin_GRL                                 = xmin_GRL_config
    C%xmax_GRL                                 = xmax_GRL_config
    C%ymin_GRL                                 = ymin_GRL_config
    C%ymax_GRL                                 = ymax_GRL_config

    C%xmin_ANT                                 = xmin_ANT_config
    C%xmax_ANT                                 = xmax_ANT_config
    C%ymin_ANT                                 = ymin_ANT_config
    C%ymax_ANT                                 = ymax_ANT_config
  
    ! Whether or not the simulation is a restart of a previous simulation
    ! ===================================================================
    
    C%is_restart                               = is_restart_config
    C%time_to_restart_from                     = time_to_restart_from_config
    
    C%filename_restart_NAM                     = filename_restart_NAM_config
    C%filename_restart_EAS                     = filename_restart_EAS_config
    C%filename_restart_GRL                     = filename_restart_GRL_config
    C%filename_restart_ANT                     = filename_restart_ANT_config
    
    ! Benchmark experiments
    ! =====================
    
    C%do_benchmark_experiment                  = do_benchmark_experiment_config
    C%choice_benchmark_experiment              = choice_benchmark_experiment_config
    C%SSA_icestream_m                          = SSA_icestream_m_config
    C%ISMIP_HOM_L                              = ISMIP_HOM_L_config
    C%ISMIP_HOM_E_Arolla_filename              = ISMIP_HOM_E_Arolla_filename_config
    C%MISMIPplus_do_tune_A_for_GL              = MISMIPplus_do_tune_A_for_GL_config
    C%MISMIPplus_xGL_target                    = MISMIPplus_xGL_target_config
    C%MISMIPplus_A_flow_initial                = MISMIPplus_A_flow_initial_config
    C%MISMIPplus_scenario                      = MISMIPplus_scenario_config
    C%MISOMIP1_scenario                        = MISOMIP1_scenario_config

    ! Mesh generation parameters
    ! ==========================
    
    C%nconmax                                  = nconmax_config
    C%alpha_min                                = alpha_min_config
    C%dz_max_ice                               = dz_max_ice_config
    C%res_max                                  = res_max_config
    C%res_max_margin                           = res_max_margin_config
    C%res_max_gl                               = res_max_gl_config
    C%res_max_cf                               = res_max_cf_config
    C%res_max_mountain                         = res_max_mountain_config
    C%res_max_coast                            = res_max_coast_config
    C%mesh_fitness_threshold                   = mesh_fitness_threshold_config
    
    ! The smallest allowed resolution
    C%res_min = MIN( MIN( MIN( MIN( C%res_max_margin, C%res_max_gl), C%res_max_cf), C%res_max_mountain), C%res_max_coast)
  
    ! Resolutions of the different square grids
    ! =========================================
    
    C%dx_grid_output                           = dx_grid_output_config
    C%dx_grid_GIA                              = dx_grid_GIA_config
    C%dx_grid_smooth                           = dx_grid_smooth_config
    
    ! High-resolution Points Of Interest (POIs)
    ! =========================================
    
    C%nPOI_NAM                                 = nPOI_NAM_config    
    C%nPOI_EAS                                 = nPOI_EAS_config    
    C%nPOI_GRL                                 = nPOI_GRL_config    
    C%nPOI_ANT                                 = nPOI_ANT_config
    C%POI_NAM_coordinates                      = POI_NAM_coordinates_config
    C%POI_EAS_coordinates                      = POI_EAS_coordinates_config
    C%POI_GRL_coordinates                      = POI_GRL_coordinates_config
    C%POI_ANT_coordinates                      = POI_ANT_coordinates_config
    C%POI_NAM_resolutions                      = POI_NAM_resolutions_config
    C%POI_EAS_resolutions                      = POI_EAS_resolutions_config
    C%POI_GRL_resolutions                      = POI_GRL_resolutions_config
    C%POI_ANT_resolutions                      = POI_ANT_resolutions_config
    
    ! Whether or not to let UFEMISM dynamically create its own output folder
   ! =======================================================================
   
    C%create_new_output_dir                    = create_new_output_dir_config
    C%output_dir                               = output_dir_config

    ! Scaled vertical coordinate zeta  
    ! ===============================
    
    C%nz     = nz_config
    ALLOCATE( C%zeta( C%nz))
    C%zeta   = zeta_config( 1:C%nz)

    ! Reference geometries (initial, present-day, and GIA equilibrium)
    ! ================================================================
    
    ! Initial geometry
    C%choice_refgeo_init_NAM                   = choice_refgeo_init_NAM_config
    C%choice_refgeo_init_EAS                   = choice_refgeo_init_EAS_config
    C%choice_refgeo_init_GRL                   = choice_refgeo_init_GRL_config
    C%choice_refgeo_init_ANT                   = choice_refgeo_init_ANT_config
    C%time_to_restart_from_NAM                 = time_to_restart_from_NAM_config
    C%time_to_restart_from_EAS                 = time_to_restart_from_EAS_config
    C%time_to_restart_from_GRL                 = time_to_restart_from_GRL_config
    C%time_to_restart_from_ANT                 = time_to_restart_from_ANT_config
    C%choice_refgeo_init_idealised             = choice_refgeo_init_idealised_config
    C%dx_refgeo_init_idealised                 = dx_refgeo_init_idealised_config
    C%filename_refgeo_init_NAM                 = filename_refgeo_init_NAM_config
    C%filename_refgeo_init_EAS                 = filename_refgeo_init_EAS_config
    C%filename_refgeo_init_GRL                 = filename_refgeo_init_GRL_config
    C%filename_refgeo_init_ANT                 = filename_refgeo_init_ANT_config

    ! Present-day geometry
    C%choice_refgeo_PD_NAM                     = choice_refgeo_PD_NAM_config
    C%choice_refgeo_PD_EAS                     = choice_refgeo_PD_EAS_config
    C%choice_refgeo_PD_GRL                     = choice_refgeo_PD_GRL_config
    C%choice_refgeo_PD_ANT                     = choice_refgeo_PD_ANT_config
    C%choice_refgeo_PD_idealised               = choice_refgeo_PD_idealised_config
    C%dx_refgeo_PD_idealised                   = dx_refgeo_PD_idealised_config
    C%filename_refgeo_PD_NAM                   = filename_refgeo_PD_NAM_config
    C%filename_refgeo_PD_EAS                   = filename_refgeo_PD_EAS_config
    C%filename_refgeo_PD_GRL                   = filename_refgeo_PD_GRL_config
    C%filename_refgeo_PD_ANT                   = filename_refgeo_PD_ANT_config

    ! GIA equilibrium geometry
    C%choice_refgeo_GIAeq_NAM                  = choice_refgeo_GIAeq_NAM_config
    C%choice_refgeo_GIAeq_EAS                  = choice_refgeo_GIAeq_EAS_config
    C%choice_refgeo_GIAeq_GRL                  = choice_refgeo_GIAeq_GRL_config
    C%choice_refgeo_GIAeq_ANT                  = choice_refgeo_GIAeq_ANT_config
    C%choice_refgeo_GIAeq_idealised            = choice_refgeo_GIAeq_idealised_config
    C%dx_refgeo_GIAeq_idealised                = dx_refgeo_GIAeq_idealised_config
    C%filename_refgeo_GIAeq_NAM                = filename_refgeo_GIAeq_NAM_config
    C%filename_refgeo_GIAeq_EAS                = filename_refgeo_GIAeq_EAS_config
    C%filename_refgeo_GIAeq_GRL                = filename_refgeo_GIAeq_GRL_config
    C%filename_refgeo_GIAeq_ANT                = filename_refgeo_GIAeq_ANT_config
    
    C%remove_Lake_Vostok                       = remove_Lake_Vostok_config

    ! Input data file paths
    ! =====================
    
    C%filename_insolation                      = filename_insolation_config
    
    C%filename_CO2_record                      = filename_CO2_record_config
    C%CO2_record_length                        = CO2_record_length_config
    C%filename_d18O_record                     = filename_d18O_record_config
    C%d18O_record_length                       = d18O_record_length_config
    
    ! Geothermal heat flux
    C%choice_geothermal_heat_flux              = choice_geothermal_heat_flux_config
    C%constant_geothermal_heat_flux            = constant_geothermal_heat_flux_config
    C%filename_geothermal_heat_flux            = filename_geothermal_heat_flux_config

    ! Ice dynamics - velocity
    ! =======================
    
    C%choice_ice_dynamics                      = choice_ice_dynamics_config
    C%n_flow                                   = n_flow_config
    C%m_enh_sheet                              = m_enh_sheet_config
    C%m_enh_shelf                              = m_enh_shelf_config
    C%choice_ice_margin                        = choice_ice_margin_config
    C%include_SSADIVA_crossterms               = include_SSADIVA_crossterms_config
    C%do_GL_subgrid_friction                   = do_GL_subgrid_friction_config
    C%do_smooth_geometry                       = do_smooth_geometry_config
    C%r_smooth_geometry                        = r_smooth_geometry_config
    
    ! Some parameters for numerically solving the SSA/DIVA
    C%DIVA_visc_it_norm_dUV_tol                = DIVA_visc_it_norm_dUV_tol_config
    C%DIVA_visc_it_nit                         = DIVA_visc_it_nit_config
    C%DIVA_visc_it_relax                       = DIVA_visc_it_relax_config
    C%DIVA_beta_max                            = DIVA_beta_max_config
    C%DIVA_err_lim                             = DIVA_err_lim_config
    C%DIVA_vel_max                             = DIVA_vel_max_config
    C%DIVA_vel_min                             = DIVA_vel_min_config
    C%DIVA_boundary_BC_u_west                  = DIVA_boundary_BC_u_west_config
    C%DIVA_boundary_BC_u_east                  = DIVA_boundary_BC_u_east_config
    C%DIVA_boundary_BC_u_south                 = DIVA_boundary_BC_u_south_config
    C%DIVA_boundary_BC_u_north                 = DIVA_boundary_BC_u_north_config
    C%DIVA_boundary_BC_v_west                  = DIVA_boundary_BC_v_west_config
    C%DIVA_boundary_BC_v_east                  = DIVA_boundary_BC_v_east_config
    C%DIVA_boundary_BC_v_south                 = DIVA_boundary_BC_v_south_config
    C%DIVA_boundary_BC_v_north                 = DIVA_boundary_BC_v_north_config
    C%DIVA_choice_matrix_solver                = DIVA_choice_matrix_solver_config
    C%DIVA_SOR_nit                             = DIVA_SOR_nit_config
    C%DIVA_SOR_tol                             = DIVA_SOR_tol_config
    C%DIVA_SOR_omega                           = DIVA_SOR_omega_config
    C%DIVA_PETSc_rtol                          = DIVA_PETSc_rtol_config
    C%DIVA_PETSc_abstol                        = DIVA_PETSc_abstol_config
  
    ! Ice dynamics - time integration
    ! ===============================
    
    C%choice_timestepping                      = choice_timestepping_config
    C%choice_ice_integration_method            = choice_ice_integration_method_config
    C%dHi_choice_matrix_solver                 = dHi_choice_matrix_solver_config
    C%dHi_SOR_nit                              = dHi_SOR_nit_config
    C%dHi_SOR_tol                              = dHi_SOR_tol_config
    C%dHi_SOR_omega                            = dHi_SOR_omega_config
    C%dHi_PETSc_rtol                           = dHi_PETSc_rtol_config
    C%dHi_PETSc_abstol                         = dHi_PETSc_abstol_config
    
    ! Predictor-corrector ice-thickness update
    C%pc_epsilon                               = pc_epsilon_config
    C%pc_k_I                                   = pc_k_I_config
    C%pc_k_p                                   = pc_k_p_config
    C%pc_eta_min                               = pc_eta_min_config
    C%pc_max_timestep_iterations               = pc_max_timestep_iterations_config
    C%pc_redo_tol                              = pc_redo_tol_config
    C%dt_min                                   = dt_min_config
  
    ! Ice thickness boundary conditions
    C%ice_thickness_west_BC                    = ice_thickness_west_BC_config
    C%ice_thickness_east_BC                    = ice_thickness_east_BC_config
    C%ice_thickness_south_BC                   = ice_thickness_south_BC_config
    C%ice_thickness_north_BC                   = ice_thickness_north_BC_config
    C%choice_mask_noice_NAM                    = choice_mask_noice_NAM_config
    C%choice_mask_noice_EAS                    = choice_mask_noice_EAS_config
    C%choice_mask_noice_GRL                    = choice_mask_noice_GRL_config
    C%choice_mask_noice_ANT                    = choice_mask_noice_ANT_config

    ! Ice dynamics - basal conditions and sliding
    ! ===========================================

    ! Sliding laws
    C%choice_sliding_law                       = choice_sliding_law_config
    C%choice_idealised_sliding_law             = choice_idealised_sliding_law_config
    C%slid_delta_v                             = slid_delta_v_config
    C%slid_Weertman_m                          = slid_Weertman_m_config
    C%slid_Coulomb_reg_q_plastic               = slid_Coulomb_reg_q_plastic_config
    C%slid_Coulomb_reg_u_threshold             = slid_Coulomb_reg_u_threshold_config
    C%slid_ZI_ut                               = slid_ZI_ut_config
    C%slid_ZI_p                                = slid_ZI_p_config

    ! Basal hydrology
    C%choice_basal_hydrology                   = choice_basal_hydrology_config
    C%Martin2011_hydro_Hb_min                  = Martin2011_hydro_Hb_min_config
    C%Martin2011_hydro_Hb_max                  = Martin2011_hydro_Hb_max_config

    ! Basal roughness / friction
    C%choice_basal_roughness                   = choice_basal_roughness_config
    C%slid_Weertman_beta_sq_uniform            = slid_Weertman_beta_sq_uniform_config
    C%slid_Coulomb_phi_fric_uniform            = slid_Coulomb_phi_fric_uniform_config
    C%slid_Tsai2015_alpha_sq_uniform           = slid_Tsai2015_alpha_sq_uniform_config
    C%slid_Tsai2015_beta_sq_uniform            = slid_Tsai2015_beta_sq_uniform_config
    C%slid_Schoof2005_alpha_sq_uniform         = slid_Schoof2005_alpha_sq_uniform_config
    C%slid_Schoof2005_beta_sq_uniform          = slid_Schoof2005_beta_sq_uniform_config
    C%slid_ZI_phi_fric_uniform                 = slid_ZI_phi_fric_uniform_config
    C%choice_param_basal_roughness             = choice_param_basal_roughness_config
    C%Martin2011till_phi_Hb_min                = Martin2011till_phi_Hb_min_config
    C%Martin2011till_phi_Hb_max                = Martin2011till_phi_Hb_max_config
    C%Martin2011till_phi_min                   = Martin2011till_phi_min_config
    C%Martin2011till_phi_max                   = Martin2011till_phi_max_config
    C%basal_roughness_filename                 = basal_roughness_filename_config
  
    ! Ice dynamics - calving
    ! ======================
    
    C%choice_calving_law                       = choice_calving_law_config
    C%calving_threshold_thickness              = calving_threshold_thickness_config
    C%do_remove_shelves                        = do_remove_shelves_config
    C%remove_shelves_larger_than_PD            = remove_shelves_larger_than_PD_config
    C%continental_shelf_calving                = continental_shelf_calving_config
    C%continental_shelf_min_height             = continental_shelf_min_height_config
  
    ! Thermodynamics and rheology
    ! ===========================
    
    C%choice_initial_ice_temperature           = choice_initial_ice_temperature_config
    C%uniform_ice_temperature                  = uniform_ice_temperature_config
    C%choice_thermo_model                      = choice_thermo_model_config
    C%choice_ice_rheology                      = choice_ice_rheology_config
    C%uniform_flow_factor                      = uniform_flow_factor_config
    C%choice_ice_heat_capacity                 = choice_ice_heat_capacity_config
    C%uniform_ice_heat_capacity                = uniform_ice_heat_capacity_config
    C%choice_ice_thermal_conductivity          = choice_ice_thermal_conductivity_config
    C%uniform_ice_thermal_conductivity         = uniform_ice_thermal_conductivity_config
  
    ! Sea level and GIA
    ! =================
    
    C%do_ocean_floodfill                       = do_ocean_floodfill_config
    C%choice_sealevel_model                    = choice_sealevel_model_config
    C%fixed_sealevel                           = fixed_sealevel_config
    C%filename_sealevel_record                 = filename_sealevel_record_config
    C%sealevel_record_length                   = sealevel_record_length_config
  
    C%choice_GIA_model                         = choice_GIA_model_config
    C%ELRA_lithosphere_flex_rigidity           = ELRA_lithosphere_flex_rigidity_config
    C%ELRA_bedrock_relaxation_time             = ELRA_bedrock_relaxation_time_config
    C%ELRA_mantle_density                      = ELRA_mantle_density_config
    
    ! Climate matrix
    ! ==============
    
    C%filename_PD_obs_climate                  = filename_PD_obs_climate_config
    C%choice_climate_matrix                    = choice_climate_matrix_config
    C%filename_GCM_snapshot_PI                 = filename_GCM_snapshot_PI_config
    C%filename_GCM_snapshot_LGM                = filename_GCM_snapshot_LGM_config
    C%filename_ICE5G_PD                        = filename_ICE5G_PD_config
    C%filename_ICE5G_LGM                       = filename_ICE5G_LGM_config
    
    C%choice_ocean_temperature_model           = choice_ocean_temperature_model_config
    C%ocean_temperature_PD                     = ocean_temperature_PD_config
    C%ocean_temperature_cold                   = ocean_temperature_cold_config
    C%ocean_temperature_warm                   = ocean_temperature_warm_config
    
    C%constant_lapserate                       = constant_lapserate_config
    
    ! Forcing
    ! =======
    
    C%choice_forcing_method                    = choice_forcing_method_config
    
    C%dT_deepwater_averaging_window            = dT_deepwater_averaging_window_config
    C%dT_deepwater_dT_surf_ratio               = dT_deepwater_dT_surf_ratio_config
    C%d18O_dT_deepwater_ratio                  = d18O_dT_deepwater_ratio_config
    
    C%dT_glob_inverse_averaging_window         = dT_glob_inverse_averaging_window_config
    C%inverse_d18O_to_dT_glob_scaling          = inverse_d18O_to_dT_glob_scaling_config
    C%CO2_inverse_averaging_window             = CO2_inverse_averaging_window_config
    C%inverse_d18O_to_CO2_scaling              = inverse_d18O_to_CO2_scaling_config
    C%inverse_d18O_to_CO2_initial_CO2          = inverse_d18O_to_CO2_initial_CO2_config
    
    ! SMB melt tuning
    ! ===============
    
    C%C_abl_constant_NAM                       = C_abl_constant_NAM_config
    C%C_abl_constant_EAS                       = C_abl_constant_EAS_config
    C%C_abl_constant_GRL                       = C_abl_constant_GRL_config
    C%C_abl_constant_ANT                       = C_abl_constant_ANT_config
    C%C_abl_Ts_NAM                             = C_abl_Ts_NAM_config
    C%C_abl_Ts_EAS                             = C_abl_Ts_EAS_config
    C%C_abl_Ts_GRL                             = C_abl_Ts_GRL_config
    C%C_abl_Ts_ANT                             = C_abl_Ts_ANT_config
    C%C_abl_Q_NAM                              = C_abl_Q_NAM_config
    C%C_abl_Q_EAS                              = C_abl_Q_EAS_config
    C%C_abl_Q_GRL                              = C_abl_Q_GRL_config
    C%C_abl_Q_ANT                              = C_abl_Q_ANT_config
    C%C_refr_NAM                               = C_refr_NAM_config
    C%C_refr_EAS                               = C_refr_EAS_config
    C%C_refr_GRL                               = C_refr_GRL_config
    C%C_refr_ANT                               = C_refr_ANT_config
    
    ! Sub-shelf melt parameterisation
    ! ===============================
    
    C%T_ocean_mean_PD_NAM                      = T_ocean_mean_PD_NAM_config
    C%T_ocean_mean_PD_EAS                      = T_ocean_mean_PD_EAS_config
    C%T_ocean_mean_PD_GRL                      = T_ocean_mean_PD_GRL_config
    C%T_ocean_mean_PD_ANT                      = T_ocean_mean_PD_ANT_config
    C%T_ocean_mean_cold_NAM                    = T_ocean_mean_cold_NAM_config
    C%T_ocean_mean_cold_EAS                    = T_ocean_mean_cold_EAS_config
    C%T_ocean_mean_cold_GRL                    = T_ocean_mean_cold_GRL_config
    C%T_ocean_mean_cold_ANT                    = T_ocean_mean_cold_ANT_config
    C%T_ocean_mean_warm_NAM                    = T_ocean_mean_warm_NAM_config
    C%T_ocean_mean_warm_EAS                    = T_ocean_mean_warm_EAS_config
    C%T_ocean_mean_warm_GRL                    = T_ocean_mean_warm_GRL_config
    C%T_ocean_mean_warm_ANT                    = T_ocean_mean_warm_ANT_config
    
    C%BMB_deepocean_PD_NAM                     = BMB_deepocean_PD_NAM_config
    C%BMB_deepocean_PD_EAS                     = BMB_deepocean_PD_EAS_config
    C%BMB_deepocean_PD_GRL                     = BMB_deepocean_PD_GRL_config
    C%BMB_deepocean_PD_ANT                     = BMB_deepocean_PD_ANT_config
    C%BMB_deepocean_cold_NAM                   = BMB_deepocean_cold_NAM_config
    C%BMB_deepocean_cold_EAS                   = BMB_deepocean_cold_EAS_config
    C%BMB_deepocean_cold_GRL                   = BMB_deepocean_cold_GRL_config
    C%BMB_deepocean_cold_ANT                   = BMB_deepocean_cold_ANT_config
    C%BMB_deepocean_warm_NAM                   = BMB_deepocean_warm_NAM_config
    C%BMB_deepocean_warm_EAS                   = BMB_deepocean_warm_EAS_config
    C%BMB_deepocean_warm_GRL                   = BMB_deepocean_warm_GRL_config
    C%BMB_deepocean_warm_ANT                   = BMB_deepocean_warm_ANT_config
    
    C%BMB_shelf_exposed_PD_NAM                 = BMB_shelf_exposed_PD_NAM_config
    C%BMB_shelf_exposed_PD_EAS                 = BMB_shelf_exposed_PD_EAS_config
    C%BMB_shelf_exposed_PD_GRL                 = BMB_shelf_exposed_PD_GRL_config
    C%BMB_shelf_exposed_PD_ANT                 = BMB_shelf_exposed_PD_ANT_config
    C%BMB_shelf_exposed_cold_NAM               = BMB_shelf_exposed_cold_NAM_config
    C%BMB_shelf_exposed_cold_EAS               = BMB_shelf_exposed_cold_EAS_config
    C%BMB_shelf_exposed_cold_GRL               = BMB_shelf_exposed_cold_GRL_config
    C%BMB_shelf_exposed_warm_NAM               = BMB_shelf_exposed_warm_NAM_config
    C%BMB_shelf_exposed_cold_ANT               = BMB_shelf_exposed_cold_ANT_config
    C%BMB_shelf_exposed_warm_EAS               = BMB_shelf_exposed_warm_EAS_config
    C%BMB_shelf_exposed_warm_GRL               = BMB_shelf_exposed_warm_GRL_config
    C%BMB_shelf_exposed_warm_ANT               = BMB_shelf_exposed_warm_ANT_config
    
    C%subshelf_melt_factor_NAM                 = subshelf_melt_factor_NAM_config
    C%subshelf_melt_factor_EAS                 = subshelf_melt_factor_EAS_config
    C%subshelf_melt_factor_GRL                 = subshelf_melt_factor_GRL_config
    C%subshelf_melt_factor_ANT                 = subshelf_melt_factor_ANT_config
    
    C%deep_ocean_threshold_depth_NAM           = deep_ocean_threshold_depth_NAM_config
    C%deep_ocean_threshold_depth_EAS           = deep_ocean_threshold_depth_EAS_config
    C%deep_ocean_threshold_depth_GRL           = deep_ocean_threshold_depth_GRL_config
    C%deep_ocean_threshold_depth_ANT           = deep_ocean_threshold_depth_ANT_config
  
    ! Which data fields will be written to the help_fields output file
    ! ================================================================
    
    C%help_field_01                            = help_field_01_config
    C%help_field_02                            = help_field_02_config
    C%help_field_03                            = help_field_03_config
    C%help_field_04                            = help_field_04_config
    C%help_field_05                            = help_field_05_config
    C%help_field_06                            = help_field_06_config
    C%help_field_07                            = help_field_07_config
    C%help_field_08                            = help_field_08_config
    C%help_field_09                            = help_field_09_config
    C%help_field_10                            = help_field_10_config
    C%help_field_11                            = help_field_11_config
    C%help_field_12                            = help_field_12_config
    C%help_field_13                            = help_field_13_config
    C%help_field_14                            = help_field_14_config
    C%help_field_15                            = help_field_15_config
    C%help_field_16                            = help_field_16_config
    C%help_field_17                            = help_field_17_config
    C%help_field_18                            = help_field_18_config
    C%help_field_19                            = help_field_19_config
    C%help_field_20                            = help_field_20_config
    C%help_field_21                            = help_field_21_config
    C%help_field_22                            = help_field_22_config
    C%help_field_23                            = help_field_23_config
    C%help_field_24                            = help_field_24_config
    C%help_field_25                            = help_field_25_config
    C%help_field_26                            = help_field_26_config
    C%help_field_27                            = help_field_27_config
    C%help_field_28                            = help_field_28_config
    C%help_field_29                            = help_field_29_config
    C%help_field_30                            = help_field_30_config
    C%help_field_31                            = help_field_31_config
    C%help_field_32                            = help_field_32_config
    C%help_field_33                            = help_field_33_config
    C%help_field_34                            = help_field_34_config
    C%help_field_35                            = help_field_35_config
    C%help_field_36                            = help_field_36_config
    C%help_field_37                            = help_field_37_config
    C%help_field_38                            = help_field_38_config
    C%help_field_39                            = help_field_39_config
    C%help_field_40                            = help_field_40_config
    C%help_field_41                            = help_field_41_config
    C%help_field_42                            = help_field_42_config
    C%help_field_43                            = help_field_43_config
    C%help_field_44                            = help_field_44_config
    C%help_field_45                            = help_field_45_config
    C%help_field_46                            = help_field_46_config
    C%help_field_47                            = help_field_47_config
    C%help_field_48                            = help_field_48_config
    C%help_field_49                            = help_field_49_config
    C%help_field_50                            = help_field_50_config
    
    ! Values to be filled into the total mask (used only for diagnostic output)
    ! ==========================================================================

    C%type_land                                = 0
    C%type_ocean                               = 1
    C%type_lake                                = 2
    C%type_sheet                               = 3
    C%type_shelf                               = 4
    C%type_coast                               = 5
    C%type_margin                              = 6
    C%type_groundingline                       = 7
    C%type_calvingfront                        = 8
    
   ! Parameters of the polar stereographic projections of the four model regions
   ! (These have to match the values used to create the input files!)
   ! ===========================================================================  
  
    C%lambda_M_NAM                             = 265._dp
    C%lambda_M_EAS                             = 40._dp
    C%lambda_M_GRL                             = 320._dp
    C%lambda_M_ANT                             = 0._dp
    C%phi_M_NAM                                = 62._dp
    C%phi_M_EAS                                = 70._dp
    C%phi_M_GRL                                = 72._dp
    C%phi_M_ANT                                = -90._dp
    C%alpha_stereo_NAM                         = 165.0923_dp
    C%alpha_stereo_EAS                         = 165.04_dp
    C%alpha_stereo_GRL                         = 164.85_dp
    C%alpha_stereo_ANT                         = 165.0263_dp

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
