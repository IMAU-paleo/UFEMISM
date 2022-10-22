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

! ===== USE modules =====
! =======================

  USE mpi

  IMPLICIT NONE

! ===== Precision =====
! =====================

  INTEGER, PARAMETER  :: dp  = KIND(1.0D0)  ! Kind of double precision numbers. Reals should be declared as: REAL(dp) :: example

! ===== Error messaging / debugging / profiling system =====
! ==========================================================

  CHARACTER(LEN=1024) :: routine_path
  INTEGER             :: n_MPI_windows
  REAL(dp)            :: mem_use_tot, mem_use_tot_max

  TYPE subroutine_resource_tracker
    ! Track the resource use (computation time, memory) of a single subroutine
    CHARACTER(LEN = 2048) :: routine_path
    REAL(dp)              :: tstart, tcomp
    INTEGER               :: n_MPI_windows_init, n_MPI_windows_final
    REAL(dp)              :: mem_use, mem_use_max
  END TYPE subroutine_resource_tracker

  TYPE( subroutine_resource_tracker), DIMENSION(:), ALLOCATABLE :: resource_tracker

! ===== Configuration variables =====
! ===================================
  ! The "_config  variables, which will be collected into a NAMELIST, and possibly replaced
  ! by the values in the external config file. Remember the "_config" extension!

  ! == Time steps and range
  ! =======================

    ! General
    REAL(dp)            :: start_time_of_run_config                    = 0.0_dp                           ! Start time (in years) of the simulations
    REAL(dp)            :: end_time_of_run_config                      = 50000.0_dp                       ! End   time (in years) of the simulations
    REAL(dp)            :: dt_mesh_min_config                          = 50._dp                           ! Minimum amount of time (in years) between mesh updates
    REAL(dp)            :: dt_coupling_config                          = 100._dp                          ! Interval of coupling (in years) between the four ice-sheets

    ! Ice dynamics
    REAL(dp)            :: dt_max_config                               = 10.0_dp                          ! Maximum time step (in years) of the ice model
    REAL(dp)            :: dt_min_config                               = 0.01_dp                          ! Smallest allowed time step [yr]
    REAL(dp)            :: dt_startup_phase_config                     = 10._dp                           ! Length of time window (in years) after start_time when dt = dt_min, to ensure smooth restarts
    REAL(dp)            :: dt_cooldown_phase_config                    = 10._dp                           ! Length of time window (in years) before end_time when dt = dt_min, to ensure smooth restarts

    ! Sub-models
    REAL(dp)            :: dt_thermo_config                            = 10.0_dp                          ! Time step (in years) for updating thermodynamics
    REAL(dp)            :: dt_climate_config                           = 10._dp                           ! Time step (in years) for updating the climate
    REAL(dp)            :: dt_ocean_config                             = 10._dp                           ! Time step (in years) for updating the ocean
    REAL(dp)            :: dt_SMB_config                               = 10._dp                           ! Time step (in years) for updating the SMB
    REAL(dp)            :: dt_BMB_config                               = 10._dp                           ! Time step (in years) for updating the BMB
    REAL(dp)            :: dt_bedrock_ELRA_config                      = 100._dp                          ! Time step (in years) for updating the bedrock deformation rate with the ELRA model
    REAL(dp)            :: dt_SELEN_config                             = 1000._dp                         ! Time step (in years) for calling SELEN

    ! Inversions
    REAL(dp)            :: dt_slid_inv_config                          = 10._dp                           ! Time step (in years) for calling the iterative inversion of basal roughness
    REAL(dp)            :: dt_SMB_inv_config                           = 50._dp                           ! Time step (in years) for calling the iterative inversion of the IMAU-ITM SMB parameters

    ! Output
    REAL(dp)            :: dt_output_config                            = 5000.0_dp                        ! Time step (in years) for writing output

  ! == Which ice sheets do we simulate?
  ! ===================================

    LOGICAL             :: do_NAM_config                               = .FALSE.                          ! North America
    LOGICAL             :: do_EAS_config                               = .FALSE.                          ! Eurasia
    LOGICAL             :: do_GRL_config                               = .FALSE.                          ! Greenland
    LOGICAL             :: do_ANT_config                               = .TRUE.                           ! Antarctica

  ! == Benchmark experiments
  ! ========================

    ! SSA_icestream (see Schoof 2006, and also Bueler and Brown 2009)
    REAL(dp)            :: SSA_icestream_m_config                      = 1                                ! Values tested by Schoof are 1, 10, and 20

    ! ISMIP-HOM (see Pattyn et al. 2008)
    REAL(dp)            :: ISMIP_HOM_L_config                          = 160000.0                         ! Domain size of the ISMIP-HOM benchmarks
    CHARACTER(LEN=256)  :: ISMIP_HOM_E_Arolla_filename_config          = 'arolla100.dat'                  ! Path to the Haut Glacier d'Arolla input file

    ! MISMIP+ (see Asay-Davis et al., 2016)
    LOGICAL             :: MISMIPplus_do_tune_A_for_GL_config          = .FALSE.                          ! Whether or not the flow factor A should be tuned for the GL position
    REAL(dp)            :: MISMIPplus_xGL_target_config                = 450000._dp                       ! Mid-channel GL position to tune the flow factor A for
    REAL(dp)            :: MISMIPplus_A_flow_initial_config            = 2.0E-17_dp                       ! Initial flow factor before tuning (or throughout the run when tuning is not used)
    CHARACTER(LEN=256)  :: MISMIPplus_scenario_config                  = ''                               ! Choose between the five MISMIP+  scenarios from Cornford   et al. (2020): ice0, ice1ra, ice1rr, ice2ra, ice2rr

    ! MISOMIP1 (see Asay-Davis et al., 2016)
    CHARACTER(LEN=256)  :: MISOMIP1_scenario_config                    = ''                               ! Choose between the four MISOMIP+ scenarios from Asay-Davis et al. (2016): IceOcean1ra, IceOcean1rr, IceOcean2ra, IceOcean2rr

  ! == Whether or not to let UFEMISM dynamically create its own output folder
  ! =========================================================================

    LOGICAL             :: create_procedural_output_dir_config         = .TRUE.                           ! Automatically create an output directory with a procedural name (e.g. results_20210720_001/)
    CHARACTER(LEN=256)  :: fixed_output_dir_config                     = 'results_UFEMISM'                ! If not, create a directory with this name instead (stops the program if this directory already exists)
    CHARACTER(LEN=256)  :: fixed_output_dir_suffix_config              = ''                               ! Suffix to put after the fixed output directory name, useful when doing ensemble runs with the template+variation set-up
    LOGICAL             :: do_write_regional_scalar_output_config      = .TRUE.
    LOGICAL             :: do_write_global_scalar_output_config        = .TRUE.

  ! == Debugging
  ! ============

    LOGICAL             :: do_write_debug_data_config                  = .FALSE.                          ! Whether or not the debug NetCDF file should be created and written to
    LOGICAL             :: do_write_grid_data_config                   = .FALSE.                          ! Whether or not the grid NetCDF files should be created and written to
    LOGICAL             :: do_check_for_NaN_config                     = .FALSE.                          ! Whether or not fields should be checked for NaN values
    LOGICAL             :: do_time_display_config                      = .FALSE.                          ! Print current model time to screen

  ! == The four model regions
  ! =========================

    ! North America
    REAL(dp)            :: lambda_M_NAM_config                         = 265._dp                          ! Longitude of the pole of the stereographic projection for the North America domain [degrees east]
    REAL(dp)            :: phi_M_NAM_config                            = 62._dp                           ! Latitude  of the pole of the stereographic projection for the North America domain [degrees north]
    REAL(dp)            :: beta_stereo_NAM_config                      = 71._dp                           ! Standard parallel     of the stereographic projection for the North America domain [degrees]
    REAL(dp)            :: xmin_NAM_config                             = -3600000._dp                     ! Western  boundary        of the North America domain [m]
    REAL(dp)            :: xmax_NAM_config                             =  3600000._dp                     ! Eastern  boundary     of the North America domain [m]
    REAL(dp)            :: ymin_NAM_config                             = -2400000._dp                     ! Southern boundary     of the North America domain [m]
    REAL(dp)            :: ymax_NAM_config                             =  2400000._dp                     ! Northern boundary     of the North America domain [m]

    ! Eurasia
    REAL(dp)            :: lambda_M_EAS_config                         = 40._dp                           ! Longitude of the pole of the stereographic projection for the Eurasia domain [degrees east]
    REAL(dp)            :: phi_M_EAS_config                            = 70._dp                           ! Latitude  of the pole of the stereographic projection for the Eurasia domain [degrees north]
    REAL(dp)            :: beta_stereo_EAS_config                      = 71._dp                           ! Standard parallel     of the stereographic projection for the Eurasia domain [degrees]
    REAL(dp)            :: xmin_EAS_config                             = -3400000._dp                     ! Western  boundary     of the Eurasia domain [m]
    REAL(dp)            :: xmax_EAS_config                             =  3400000._dp                     ! Eastern  boundary     of the Eurasia domain [m]
    REAL(dp)            :: ymin_EAS_config                             = -2080000._dp                     ! Southern boundary     of the Eurasia domain [m]
    REAL(dp)            :: ymax_EAS_config                             =  2080000._dp                     ! Northern boundary     of the Eurasia domain [m]

    ! Greenland
    REAL(dp)            :: lambda_M_GRL_config                         = -45._dp                          ! Longitude of the pole of the stereographic projection for the Greenland domain [degrees east]
    REAL(dp)            :: phi_M_GRL_config                            = 90._dp                           ! Latitude  of the pole of the stereographic projection for the Greenland domain [degrees north]
    REAL(dp)            :: beta_stereo_GRL_config                      = 70._dp                           ! Standard parallel     of the stereographic projection for the Greenland domain [degrees]
    REAL(dp)            :: xmin_GRL_config                             =  -720000._dp                     ! Western  boundary     of the Greenland domain [m]
    REAL(dp)            :: xmax_GRL_config                             =   960000._dp                     ! Eastern  boundary     of the Greenland domain [m]
    REAL(dp)            :: ymin_GRL_config                             = -3450000._dp                     ! Southern boundary     of the Greenland domain [m]
    REAL(dp)            :: ymax_GRL_config                             =  -570000._dp                     ! Northern boundary     of the Greenland domain [m]

    ! Antarctica
    REAL(dp)            :: lambda_M_ANT_config                         = 0._dp                            ! Longitude of the pole of the stereographic projection for the Antarctica domain [degrees east]
    REAL(dp)            :: phi_M_ANT_config                            = -90._dp                          ! Latitude  of the pole of the stereographic projection for the Antarctica domain [degrees north]
    REAL(dp)            :: beta_stereo_ANT_config                      = 71._dp                           ! Standard parallel     of the stereographic projection for the Antarctica domain [degrees]
    REAL(dp)            :: xmin_ANT_config                             = -3300000._dp                     ! Western  boundary     of the Antarctica domain [m]
    REAL(dp)            :: xmax_ANT_config                             =  3300000._dp                     ! Eastern  boundary     of the Antarctica domain [m]
    REAL(dp)            :: ymin_ANT_config                             = -3300000._dp                     ! Southern boundary     of the Antarctica domain [m]
    REAL(dp)            :: ymax_ANT_config                             =  3300000._dp                     ! Northern boundary     of the Antarctica domain [m]

  ! == Mesh
  ! =======

    ! Generation
    LOGICAL             :: use_submesh_config                          = .FALSE.                          ! Generate the mesh in parallel using multi-process method

    ! Parameters
    INTEGER             :: nconmax_config                              = 32                               ! Maximum number of vertex connections
    REAL(dp)            :: alpha_min_config                            = 0.55_dp                          ! Minimum internal angle of triangles (0.4363 = 25 degrees)
    REAL(dp)            :: dz_max_ice_config                           = 20000._dp                        ! Maximum allowed 2nd order surface deviation over ice
    REAL(dp)            :: res_max_config                              = 800._dp                          ! Maximum allowed resolution                            [km]
    REAL(dp)            :: res_max_ice_config                          = 40._dp                           ! Maximum allowed resolution over ice                   [km]
    REAL(dp)            :: res_max_margin_config                       = 40._dp                           ! Maximum allowed resolution over land-based ice margin [km]
    REAL(dp)            :: res_max_gl_config                           = 40._dp                           !                                 grounding line        [km]
    REAL(dp)            :: res_max_cf_config                           = 40._dp                           !                                 calving front         [km]
    REAL(dp)            :: res_max_mountain_config                     = 40._dp                           !                                 mountains             [km]
    REAL(dp)            :: res_max_coast_config                        = 40._dp                           !                                 coastline             [km]
    REAL(dp)            :: mesh_fitness_threshold_config               = 0.95_dp                          ! Minimum allowed mesh fitness (fraction of triangles that are not Bad) before mesh updating

    ! Forced mesh update
    LOGICAL             :: do_force_mesh_update_config                 = .FALSE.                          ! Force a mesh update once after do_force_mesh_update_after years
    REAL(dp)            :: do_force_mesh_update_after_config           = 10._dp                           ! Minimum years after which a force update occurs (model needs some time before first update)

  ! == Resolutions of the different square grids
  ! ============================================

    REAL(dp)            :: dx_grid_output_config                       = 40000._dp                        ! Resolution of the square grid used for writing output                [m]
    REAL(dp)            :: dx_grid_GIA_config                          = 100000._dp                       ! Resolution of the square grid used for GIA modelling (ELRA or SELEN) [m]
    REAL(dp)            :: dx_grid_smooth_config                       = 50000._dp                        ! Resolution of the square grid used for data smoothing                [m]

  ! == High-resolution Points Of Interest (POIs)
  ! ============================================

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

  ! ISMIP-style output
  ! ==================

    LOGICAL             :: do_write_ISMIP_output_config                = .FALSE.                          ! Whether or not to create a set of ISMIP output files
    CHARACTER(LEN=256)  :: ISMIP_output_group_code_config              = 'IMAU'                           ! Code for the group      name in the ISMIP output file names
    CHARACTER(LEN=256)  :: ISMIP_output_model_code_config              = 'UFEMISM'                        ! Code for the model      name in the ISMIP output file names
    CHARACTER(LEN=256)  :: ISMIP_output_experiment_code_config         = 'test'                           ! Code for the experiment name in the ISMIP output file names
    CHARACTER(LEN=256)  :: ISMIP_output_basetime_config                = 'YYYY-MM-DD'                     ! Basetime for the ISMIP output files (e.g. '1900-01-01')

  ! == The scaled vertical coordinate zeta, used mainly in thermodynamics
  ! =====================================================================

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

  ! == Reference geometries (initial, present-day, and GIA equilibrium)
  ! ===================================================================

    ! Some pre-processing stuff for reference ice geometry
    REAL(dp)            :: refgeo_Hi_min_config                        = 2.0_dp                           ! Remove ice thinner than this value in the reference ice geometry. Particularly useful for BedMachine Greenland, which somehow covers the entire tundra with half a meter of ice...
    LOGICAL             :: remove_Lake_Vostok_config                   = .TRUE.                           ! Remove Lake Vostok when running Antarctic simulations

    ! == Initial geometry
    ! ===================

    CHARACTER(LEN=256)  :: choice_refgeo_init_NAM_config               = 'realistic'                      ! Choice of initial geometry for North America; can be "idealised", "realistic", or "restart"
    CHARACTER(LEN=256)  :: choice_refgeo_init_EAS_config               = 'realistic'                      ! Choice of initial geometry for Eurasia      ; can be "idealised", "realistic", or "restart"
    CHARACTER(LEN=256)  :: choice_refgeo_init_GRL_config               = 'realistic'                      ! Choice of initial geometry for Greenland    ; can be "idealised", "realistic", or "restart"
    CHARACTER(LEN=256)  :: choice_refgeo_init_ANT_config               = 'realistic'                      ! Choice of initial geometry for Antarctica   ; can be "idealised", "realistic", or "restart"
    ! Idealised settings
    CHARACTER(LEN=256)  :: choice_refgeo_init_idealised_config         = 'flatearth'                      ! Choice of idealised initial geometry; see "generate_idealised_geometry" in reference_fields_module for options
    REAL(dp)            :: dx_refgeo_init_idealised_config             = 5000._dp                         ! Resolution of square grid used for idealised initial geometry
    ! Realistic settings
    CHARACTER(LEN=256)  :: filename_refgeo_init_NAM_config             = 'data/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_init_EAS_config             = 'data/ETOPO1/Eurasia_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_init_GRL_config             = 'data/Bedmachine_Greenland/BedMachine_Greenland_v4_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_init_ANT_config             = 'data/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'

    ! == Present-day geometry
    ! =======================

    CHARACTER(LEN=256)  :: choice_refgeo_PD_NAM_config                 = 'realistic'                      ! Choice of present-day geometry for North America; can be "idealised", or "realistic"
    CHARACTER(LEN=256)  :: choice_refgeo_PD_EAS_config                 = 'realistic'                      ! Choice of present-day geometry for Eurasia      ; can be "idealised", or "realistic"
    CHARACTER(LEN=256)  :: choice_refgeo_PD_GRL_config                 = 'realistic'                      ! Choice of present-day geometry for Greenland    ; can be "idealised", or "realistic"
    CHARACTER(LEN=256)  :: choice_refgeo_PD_ANT_config                 = 'realistic'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "realistic"
    ! Idealised settings
    CHARACTER(LEN=256)  :: choice_refgeo_PD_idealised_config           = 'flatearth'                      ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    REAL(dp)            :: dx_refgeo_PD_idealised_config               = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry
    ! Realistic settings
    CHARACTER(LEN=256)  :: filename_refgeo_PD_NAM_config               = 'data/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_PD_EAS_config               = 'data/ETOPO1/Eurasia_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_PD_GRL_config               = 'data/Bedmachine_Greenland/BedMachine_Greenland_v4_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_PD_ANT_config               = 'data/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'

    ! == GIA equilibrium geometry
    ! ===========================

    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_NAM_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for North America; can be "idealised", or "realistic"
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_EAS_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for Eurasia      ; can be "idealised", or "realistic"
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_GRL_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for Greenland    ; can be "idealised", or "realistic"
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_ANT_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for Antarctica   ; can be "idealised", or "realistic"
    ! Idealised settings
    CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_idealised_config        = 'flatearth'                      ! Choice of idealised GIA equilibrium geometry; see "generate_idealised_geometry" in reference_fields_module for options
    REAL(dp)            :: dx_refgeo_GIAeq_idealised_config            = 5000._dp                         ! Resolution of square grid used for idealised GIA equilibrium geometry
    ! Realistic settings
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_NAM_config            = 'data/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_EAS_config            = 'data/ETOPO1/Eurasia_ETOPO1_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_GRL_config            = 'data/Bedmachine_Greenland/BedMachine_Greenland_v4_5km.nc'
    CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_ANT_config            = 'data/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'

  ! == Whether or not the simulation is a restart of a previous simulation
  ! ======================================================================

    LOGICAL             :: is_restart_config                           = .FALSE.

    ! Time in the restart file from which the data will be extracted
    REAL(dp)            :: time_to_restart_from_NAM_config             = 0._dp                            ! Can be different from C%start_time_of_run, be careful though
    REAL(dp)            :: time_to_restart_from_EAS_config             = 0._dp
    REAL(dp)            :: time_to_restart_from_GRL_config             = 0._dp
    REAL(dp)            :: time_to_restart_from_ANT_config             = 0._dp

    ! Initial model state when restarting from a previous run
    CHARACTER(LEN=256)  :: filename_restart_NAM_config                 = 'filename_restart_NAM_placeholder'
    CHARACTER(LEN=256)  :: filename_restart_EAS_config                 = 'filename_restart_EAS_placeholder'
    CHARACTER(LEN=256)  :: filename_restart_GRL_config                 = 'filename_restart_GRL_placeholder'
    CHARACTER(LEN=256)  :: filename_restart_ANT_config                 = 'filename_restart_ANT_placeholder'

  ! == Global forcing (insolation, CO2, d18O, geothermal heat flux)
  ! ===============================================================

    ! Possible choice_forcing_method options:
    ! 'none'                 : No global forcing used at all; climate or SMB are fully parameterised or directly prescribed
    ! 'd18O_inverse_dT_glob' : Use the inverse routine with the specified d18O record to calculate a global temperature offset (e.g. de Boer et al., 2013)
    ! 'CO2_direct'           : Use the specified CO2 record to force the climate matrix (e.g. Berends et al., 2018)
    ! 'd18O_inverse_CO2'     : Use the inverse routine with the specified d18O record to calculate CO2 and then force the climate matrix (e.g. Berends et al., 2019)
    CHARACTER(LEN=256)  :: choice_forcing_method_config                = 'CO2_direct'

    ! Insolation forcing (NetCDF)
    CHARACTER(LEN=256)  :: choice_insolation_forcing_config            = 'realistic'                      ! Choice of insolation forcing: "none", "static", "realistic"
    REAL(dp)            :: static_insolation_time_config               = 0._dp                            ! Keep insolation values fixed to this time when choice_insolation_forcing = 'static'
    CHARACTER(LEN=256)  :: filename_insolation_config                  = 'data/Insolation/Laskar_etal_2004_insolation.nc'


    ! CO2 record (ASCII text file, so the number of rows needs to be specified)
    CHARACTER(LEN=256)  :: filename_CO2_record_config                  = 'data/CO2/EPICA_CO2_Bereiter_2015_100yr.dat'
    INTEGER             :: CO2_record_length_config                    = 8001

    ! d18O record (ASCII text file, so the number of rows needs to be specified)
    CHARACTER(LEN=256)  :: filename_d18O_record_config                 = 'data/d18O/Ahn2017_d18O.dat'
    INTEGER             :: d18O_record_length_config                   = 2051

    ! Geothermal heat flux
    CHARACTER(LEN=256)  :: choice_geothermal_heat_flux_config          = 'spatial'                        ! Choice of geothermal heat flux; can be 'constant' or 'spatial'
    REAL(dp)            :: constant_geothermal_heat_flux_config        = 1.72E06_dp                       ! Geothermal Heat flux [J m^-2 yr^-1] Sclater et al. (1980)
    CHARACTER(LEN=256)  :: filename_geothermal_heat_flux_config        = 'data/GHF/geothermal_heatflux_ShapiroRitzwoller2004_global_1x1_deg.nc'

    ! Parameters for calculating modelled benthic d18O
    LOGICAL             :: do_calculate_benthic_d18O_config            = .FALSE.                          ! Whether or not to calculate modelled benthic d18O (set to .FALSE. for e.g. idealised-geometry experiments, future projections)
    REAL(dp)            :: dT_deepwater_averaging_window_config        = 3000                             ! Time window (in yr) over which global mean temperature anomaly is averaged to find the deep-water temperature anomaly
    REAL(dp)            :: dT_deepwater_dT_surf_ratio_config           = 0.25_dp                          ! Ratio between global mean surface temperature change and deep-water temperature change
    REAL(dp)            :: d18O_dT_deepwater_ratio_config              = -0.28_dp                         ! Ratio between deep-water temperature change and benthic d18O change

    ! Parameters for the inverse routine
    REAL(dp)            :: dT_glob_inverse_averaging_window_config     = 2000._dp                         ! Time window (in yr) over which global mean temperature anomaly is averaged before changing it with the inverse routine
    REAL(dp)            :: inverse_d18O_to_dT_glob_scaling_config      = 20._dp                           ! Scaling factor between modelled d18O anomaly and prescribed temperature anomaly change (value from de Boer et al., 2013)
    REAL(dp)            :: CO2_inverse_averaging_window_config         = 2000._dp                         ! Time window (in yr) over which CO2                             is averaged before changing it with the inverse routine
    REAL(dp)            :: inverse_d18O_to_CO2_scaling_config          = 68._dp                           ! Scaling factor between modelled d18O anomaly and modelled CO2 change (value from Berends et al., 2019)
    REAL(dp)            :: inverse_d18O_to_CO2_initial_CO2_config      = 280._dp                          ! CO2 value at the start of the simulation when using the inverse method to calculate CO2

  ! == Ice dynamics - velocity
  ! ==========================

    CHARACTER(LEN=256)  :: choice_ice_dynamics_config                  = 'DIVA'                           ! Choice of ice-dynamical approximation: "none" (= fixed geometry), "SIA", "SSA", "SIA/SSA", "DIVA"
    REAL(dp)            :: n_flow_config                               = 3.0_dp                           ! Exponent in Glen's flow law
    REAL(dp)            :: m_enh_sheet_config                          = 1.0_dp                           ! Ice flow enhancement factor for grounded ice
    REAL(dp)            :: m_enh_shelf_config                          = 1.0_dp                           ! Ice flow enhancement factor for floating ice
    CHARACTER(LEN=256)  :: choice_ice_margin_config                    = 'infinite_slab'                  ! Choice of ice margin boundary conditions: "BC", "infinite_slab"
    LOGICAL             :: do_hybrid_Bernales2017_config               = .TRUE.                           ! Apply SStA + reduced SIA hybrid scheme when using the SIA/SSA method
    REAL(dp)            :: vel_ref_Bernales2017_config                 = 30._dp                           ! Reference "onset" velocity for an ice stream (point of half SIA reduction)
    LOGICAL             :: include_SSADIVA_crossterms_config           = .TRUE.                           ! Whether or not to include the "cross-terms" of the SSA/DIVA
    LOGICAL             :: do_GL_subgrid_friction_config               = .TRUE.                           ! Whether or not to scale basal friction with the sub-grid grounded fraction (needed to get proper GL migration; only turn this off for showing the effect on the MISMIP_mod results!)
    LOGICAL             :: do_smooth_geometry_config                   = .FALSE.                          ! Whether or not to smooth the model geometry (bedrock + initial ice thickness)
    REAL(dp)            :: r_smooth_geometry_config                    = 0.5_dp                           ! Geometry smoothing radius (in number of grid cells)

    ! Some parameters for numerically solving the SSA/DIVA
    REAL(dp)            :: DIVA_visc_it_norm_dUV_tol_config            = 1E-2_dp                          ! Successive solutions of UV in the effective viscosity iteration must not differ by more than this amount (on average)
    INTEGER             :: DIVA_visc_it_nit_config                     = 50                               ! Maximum number of effective viscosity iterations
    REAL(dp)            :: DIVA_visc_it_relax_config                   = 0.4_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
    REAL(dp)            :: DIVA_epsilon_sq_0_config                    = 1E-15_dp                         ! Normalisation term so that zero velocity gives non-zero viscosity
    REAL(dp)            :: DIVA_visc_eff_min_config                    = 1E3_dp                           ! Minimum value for effective viscosity
    REAL(dp)            :: DIVA_beta_max_config                        = 1E20_dp                          ! Maximum value for basal friction coefficient
    REAL(dp)            :: DIVA_vel_max_config                         = 5000._dp                         ! DIVA velocities are limited to this value
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
    REAL(dp)            :: DIVA_PETSc_abstol_config                    = 2.5_dp                           ! DIVA PETSc solver - stop criterion, absolute difference

    ! Velocity wind-up
    REAL(dp)            :: windup_total_years_config                   = 0._dp                            ! Go back in time and run the velocity solver for this amount of years before actual run starts

  ! == Ice dynamics - time integration
  ! ==================================

    CHARACTER(LEN=256)  :: choice_timestepping_config                  = 'pc'                             ! Choice of timestepping method: "direct", "pc" (NOTE: 'direct' does not work with DIVA ice dynamcis!)
    CHARACTER(LEN=256)  :: choice_ice_integration_method_config        = 'explicit'                       ! Choice of ice thickness integration scheme: "none" (i.e. unchanging geometry), "explicit", "semi-implicit"
    CHARACTER(LEN=256)  :: dHi_choice_matrix_solver_config             = 'SOR'                            ! Choice of matrix solver for the semi-implicit ice thickness equation: "SOR", "PETSc" [NOT USED]
    INTEGER             :: dHi_SOR_nit_config                          = 3000                             ! dHi SOR   solver - maximum number of iterations [NOT USED]
    REAL(dp)            :: dHi_SOR_tol_config                          = 2.5_dp                           ! dHi SOR   solver - stop criterion, absolute difference [NOT USED]
    REAL(dp)            :: dHi_SOR_omega_config                        = 1.3_dp                           ! dHi SOR   solver - over-relaxation parameter [NOT USED]
    REAL(dp)            :: dHi_PETSc_rtol_config                       = 0.001_dp                         ! dHi PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached) [NOT USED]
    REAL(dp)            :: dHi_PETSc_abstol_config                     = 0.001_dp                         ! dHi PETSc solver - stop criterion, absolute difference [NOT USED]

    ! Predictor-corrector ice-thickness update
    REAL(dp)            :: pc_epsilon_config                           = 3._dp                            ! Target truncation error in dHi_dt [m/yr] (epsilon in Robinson et al., 2020, Eq. 33)
    REAL(dp)            :: pc_k_I_config                               = 0.2_dp                           ! Exponent k_I in  Robinson et al., 2020, Eq. 33
    REAL(dp)            :: pc_k_p_config                               = 0.2_dp                           ! Exponent k_p in  Robinson et al., 2020, Eq. 33
    REAL(dp)            :: pc_eta_min_config                           = 1E-8_dp                          ! Normalisation term in estimation of the truncation error (Robinson et al., Eq. 32)
    INTEGER             :: pc_max_timestep_iterations_config           = 5                                ! Maximum number of iterations of each time step [NOT USED]
    REAL(dp)            :: pc_redo_tol_config                          = 10._dp                           ! Maximum allowed truncation error (any higher and the timestep is decreased) [NOT USED]

    ! Ice thickness boundary conditions
    CHARACTER(LEN=256)  :: ice_thickness_west_BC_config                = 'zero'                           ! Choice of boundary conditions for ice thickness at the domain boundary: "infinite", "periodic", "zero", "ISMIP_HOM_F"
    CHARACTER(LEN=256)  :: ice_thickness_east_BC_config                = 'zero'
    CHARACTER(LEN=256)  :: ice_thickness_south_BC_config               = 'zero'
    CHARACTER(LEN=256)  :: ice_thickness_north_BC_config               = 'zero'
    CHARACTER(LEN=256)  :: choice_mask_noice_NAM_config                = 'NAM_remove_GRL'                 ! Choice of mask_noice configuration
    CHARACTER(LEN=256)  :: choice_mask_noice_EAS_config                = 'EAS_remove_GRL'
    CHARACTER(LEN=256)  :: choice_mask_noice_GRL_config                = 'GRL_remove_Ellesmere'
    CHARACTER(LEN=256)  :: choice_mask_noice_ANT_config                = 'none'                           ! For Antarctica, additional choices are included for certain idealised-geometry experiments: "MISMIP_mod", "MISMIP+"

    ! Partially fixed geometry, useful for initialisation and inversion runs
    ! A value of 1 means fully fixed, 0 means fully free, values in between
    ! create a "delayed" ice thickness evolution.
    REAL(dp)            :: fixed_sheet_geometry_config                 = 0._dp                            ! Keep geometry of grounded ice fixed
    REAL(dp)            :: fixed_shelf_geometry_config                 = 0._dp                            ! Keep geometry of floating ice fixed
    REAL(dp)            :: fixed_grounding_line_g_config               = 0._dp                            ! Keep ice thickness at the grounded side of grounding line fixed
    REAL(dp)            :: fixed_grounding_line_f_config               = 0._dp                            ! Keep ice thickness at the floating side of grounding line fixed

    ! Memory of first dHi_dt of simulation
    INTEGER             :: dHi_dt_window_size_config                   = 1000                             ! Number of previous time steps used to compute a running average of dHi_dt

  ! == Ice dynamics - basal conditions and sliding
  ! ==============================================

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
    CHARACTER(LEN=256)  :: choice_basal_roughness_config               = 'parameterised'                  ! "uniform", "parameterised", "prescribed", or "restart"
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
    CHARACTER(LEN=256)  :: basal_roughness_restart_type_config         = 'average'                        ! Values from previous run: "last" (last output) or "average" (running average)
    LOGICAL             :: do_basal_roughness_remap_adjustment_config  = .TRUE.                           ! If TRUE, adjust bed roughness based on previous dH_dt history after a mesh update/remap

    ! Bed roughness inversion
    LOGICAL             :: do_slid_inv_config                          = .FALSE.                          ! Whether or not to perform an iterative inversion of bed roughness
    CHARACTER(LEN=256)  :: choice_slid_inv_method_config               = 'Bernales2017'                   ! Choice of iterative inversion method: "Bernales2017", "Berends2022"
    REAL(dp)            :: slid_inv_t_start_config                     = -9.9E9_dp                        ! Minimum model time when the inversion is allowed
    REAL(dp)            :: slid_inv_t_end_config                       = +9.9E9_dp                        ! Maximum model time when the inversion is allowed
    REAL(dp)            :: slid_inv_phi_min_config                     = 2._dp                            ! Minimum value of phi_fric allowed during inversion
    REAL(dp)            :: slid_inv_phi_max_config                     = 30._dp                           ! Maximum value of phi_fric allowed during inversion
    CHARACTER(LEN=256)  :: slid_inv_filename_output_config             = 'bed_roughness_inv.nc'           ! NetCDF file where the final inverted basal roughness will be saved
    INTEGER             :: slid_inv_window_size_config                 = 1000                             ! Number of previous time steps used to compute a running average of inverted values

    LOGICAL             :: do_slid_inv_Bernales2017_smooth_config      = .FALSE.                          ! If set to TRUE, inverted basal roughness is smoothed
    LOGICAL             :: do_slid_inv_Bernales2017_extrap_config      = .FALSE.                          ! If set to TRUE, inverted basal roughness is extrapolated over ice-free regions
    REAL(dp)            :: slid_inv_Bernales2017_hi_scale_config       = 10000._dp                        ! Scaling constant for inversion procedure [m]
    REAL(dp)            :: slid_inv_Bernales2017_smooth_r_config       = 500._dp                          ! Smoothing radius for inversion procedure [m]
    REAL(dp)            :: slid_inv_Bernales2017_smooth_w_config       = .01_dp                           ! Weight given to the smoothed roughness (1 = full smoothing applied)
    REAL(dp)            :: slid_inv_Bernales2017_tol_diff_config       = 100._dp                          ! Minimum ice thickness difference [m] that triggers inversion
    REAL(dp)            :: slid_inv_Bernales2017_tol_frac_config       = 1.0_dp                           ! Minimum ratio between ice thickness difference and reference value that triggers inversion

    REAL(dp)            :: slid_inv_Berends2022_tauc_config            = 10._dp                           ! Timescale       in the Berends2022 geometry/velocity-based basal inversion method [yr]
    REAL(dp)            :: slid_inv_Berends2022_H0_config              = 100._dp                          ! First  thickness scale in the Berends2022 geometry/velocity-based basal inversion method [m]
    REAL(dp)            :: slid_inv_Berends2022_u0_config              = 250._dp                          ! First  velocity  scale in the Berends2022 geometry/velocity-based basal inversion method [m/yr]
    REAL(dp)            :: slid_inv_Berends2022_Hi_scale_config        = 300._dp                          ! Second thickness scale in the Berends2022 geometry/velocity-based basal inversion method [m]
    REAL(dp)            :: slid_inv_Berends2022_u_scale_config         = 3000._dp                         ! Second velocity  scale in the Berends2022 geometry/velocity-based basal inversion method [m/yr]
    CHARACTER(LEN=256)  :: slid_inv_target_velocity_filename_config    = ''                               ! NetCDF file where the target velocities are read in the Berends2022 geometry/velocity-based basal inversion methods

  ! == Ice dynamics - calving
  ! =========================

    CHARACTER(LEN=256)  :: choice_calving_law_config                   = 'threshold_thickness'            ! Choice of calving law: "none", "threshold_thickness"
    REAL(dp)            :: calving_threshold_thickness_shelf_config    = 100._dp                          ! Threshold ice thickness for ice shelf calving front in the "threshold_thickness" calving law
    REAL(dp)            :: calving_threshold_thickness_sheet_config    = 0._dp                            ! Threshold ice thickness for ice sheet calving front in the "threshold_thickness" calving law
    INTEGER             :: max_calving_rounds_config                   = 20                               ! Maximum number of calving loops during chain reaction
    LOGICAL             :: do_remove_shelves_config                    = .FALSE.                          ! If set to TRUE, all floating ice is always instantly removed (used in the ABUMIP-ABUK experiment)
    LOGICAL             :: remove_shelves_larger_than_PD_config        = .FALSE.                          ! If set to TRUE, all floating ice beyond the present-day calving front is removed (used for some Antarctic spin-ups)
    LOGICAL             :: continental_shelf_calving_config            = .FALSE.                          ! If set to TRUE, all ice beyond the continental shelf edge (set by a maximum depth) is removed
    REAL(dp)            :: continental_shelf_min_height_config         = -2000._dp                        ! Maximum depth of the continental shelf
    REAL(dp)            :: minimum_ice_thickness_config                = 0._dp                            ! If ice anywhere is thinner than this, remove it

  ! == Thermodynamics and rheology
  ! ==============================

    CHARACTER(LEN=256)  :: choice_initial_ice_temperature_config       = 'Robin'                          ! Choice of initial ice temperature profile: "uniform", "linear", "Robin", "restart"
    REAL(dp)            :: uniform_ice_temperature_config              = 270._dp                          ! Uniform ice temperature (applied when choice_initial_ice_temperature_config = "uniform")
    CHARACTER(LEN=256)  :: choice_thermo_model_config                  = '3D_heat_equation'               ! Choice of thermodynamical model: "none", "3D_heat_equation"
    CHARACTER(LEN=256)  :: choice_ice_rheology_config                  = 'Huybrechts1992'                 ! Choice of ice rheology model: "uniform", "Huybrechts1992", "MISMIP_mod"
    REAL(dp)            :: uniform_flow_factor_config                  = 1E-16_dp                         ! Uniform ice flow factor (applied when choice_ice_rheology_model_config = "uniform")
    CHARACTER(LEN=256)  :: choice_ice_heat_capacity_config             = 'Pounder1965'                    ! Choice of ice heat capacity model: "uniform", "Pounder1965"
    REAL(dp)            :: uniform_ice_heat_capacity_config            = 2009._dp                         ! Uniform ice heat capacity (applied when choice_ice_heat_capacity_config = "uniform")
    CHARACTER(LEN=256)  :: choice_ice_thermal_conductivity_config      = 'Ritz1987'                       ! Choice of ice heat capacity model: "uniform", "Ritz1987"
    REAL(dp)            :: uniform_ice_thermal_conductivity_config     = 6.626958E7_dp                    ! Uniform ice thermal conductivity (applied when choice_ice_thermal_conductivity_config = "uniform")

  ! == Climate
  ! ==========

    CHARACTER(LEN=256)  :: choice_climate_model_config                 = 'matrix_warm_cold'               ! Choice of climate model: "none", "idealised", "PD_obs", "PD_dTglob", "matrix_warm_cold", "direct_global", "direct_regional"
    CHARACTER(LEN=256)  :: choice_idealised_climate_config             = 'EISMINT1_A'

    ! NetCDF files containing direct global/regional climate forcing
    CHARACTER(LEN=256)  :: filename_direct_global_climate_config       = ''
    CHARACTER(LEN=256)  :: filename_direct_regional_climate_NAM_config = ''
    CHARACTER(LEN=256)  :: filename_direct_regional_climate_EAS_config = ''
    CHARACTER(LEN=256)  :: filename_direct_regional_climate_GRL_config = ''
    CHARACTER(LEN=256)  :: filename_direct_regional_climate_ANT_config = ''

    ! NetCDF file containing the present-day observed climate (e.g. ERA40)
    CHARACTER(LEN=256)  :: filename_PD_obs_climate_config              = 'data/ERA40/ERA40_climate_global.nc'

    ! GCM snapshots in the matrix_warm_cold option
    CHARACTER(LEN=256)  :: filename_climate_snapshot_PI_config         = 'data/GCM_snapshots/Singarayer_Valdes_2010_PI_Control.nc'
    CHARACTER(LEN=256)  :: filename_climate_snapshot_warm_config       = 'data/GCM_snapshots/Singarayer_Valdes_2010_PI_Control.nc'
    CHARACTER(LEN=256)  :: filename_climate_snapshot_cold_config       = 'data/GCM_snapshots/Singarayer_Valdes_2010_LGM.nc'

    REAL(dp)            :: constant_lapserate_config                   = 0.008_dp                         ! Constant atmospheric lapse rate [K m^-1]

    ! Orbit time and CO2 concentration of the warm and cold snapshots
    REAL(dp)            :: matrix_high_CO2_level_config                = 280._dp                          ! CO2 level  pertaining to the warm climate (PI  level default)
    REAL(dp)            :: matrix_low_CO2_level_config                 = 190._dp                          ! CO2 level  pertaining to the cold climate (LGM level default)
    REAL(dp)            :: matrix_warm_orbit_time_config               = 0._dp                            ! Orbit time pertaining to the warm climate (PI default)
    REAL(dp)            :: matrix_cold_orbit_time_config               = -21000._dp                       ! Orbit time pertaining to the cold climate (LGM default)

    ! Whether or not to apply a bias correction to the GCM snapshots
    LOGICAL             :: climate_matrix_biascorrect_warm_config      = .TRUE.                           ! Whether or not to apply a bias correction (modelled vs observed PI climate) to the "warm" GCM snapshot
    LOGICAL             :: climate_matrix_biascorrect_cold_config      = .TRUE.                           ! Whether or not to apply a bias correction (modelled vs observed PI climate) to the "cold" GCM snapshot

    ! Glacial index
    LOGICAL             :: switch_glacial_index_config                 = .FALSE.                          ! If a glacial index is used, warm/cold weights will depend only on CO2

    ! Iterative adjustment of precipitation and temperature
    LOGICAL             :: do_clim_inv_config                          = .FALSE.                          ! Whether or not to adjust climate fields in tricky areas
    REAL(dp)            :: clim_inv_t_start_config                     = -9.9E9_dp                        ! Minimum model time when the inversion is allowed
    REAL(dp)            :: clim_inv_t_end_config                       = +9.9E9_dp                        ! Maximum model time when the inversion is allowed
    REAL(dp)            :: clim_inv_Hb_min_config                      = 1500._dp                         ! Minimum bedrock elevation where the inversion operates
    REAL(dp)            :: clim_inv_Hi_max_config                      = 100._dp                          ! Maximum ice thickness where the inversion operates
    INTEGER             :: clim_inv_window_size_config                 = 500                              ! Number of previous time steps used to compute a running average of inverted values

  ! == Ocean
  ! ========

    CHARACTER(LEN=256)  :: choice_ocean_model_config                   = 'matrix_warm_cold'               ! Choice of ocean model: "none", "idealised", "uniform_warm_cold", "PD_obs", "matrix_warm_cold"
    CHARACTER(LEN=256)  :: choice_idealised_ocean_config               = 'MISMIP+_warm'                   ! Choice of idealised ocean: 'MISMIP+_warm', 'MISMIP+_cold', 'MISOMIP1', 'Reese2018_ANT'

    ! NetCDF file containing the present-day observed ocean (WOA18) (NetCDF)
    CHARACTER(LEN=256)  :: filename_PD_obs_ocean_config                = 'data/WOA/woa18_decav_ts00_04_remapcon_r360x180_NaN.nc'
    CHARACTER(LEN=256)  :: name_ocean_temperature_obs_config           = 't_an'                           ! E.g. objectively analysed mean (t_an) or statistical mean (t_mn)
    CHARACTER(LEN=256)  :: name_ocean_salinity_obs_config              = 's_an'                           ! E.g. objectively analysed mean (s_an) or statistical mean (s_mn)
    LOGICAL             :: use_inverted_ocean_config                   = .FALSE.                          ! Whether to combine the PD ocean data with inverted values

    ! GCM snapshots in the matrix_warm_cold option
    CHARACTER(LEN=256)  :: filename_GCM_ocean_snapshot_PI_config       = 'data/COSMOS_ocean_examples/COSMOS_PI_oceanTS_prep.nc'
    CHARACTER(LEN=256)  :: filename_GCM_ocean_snapshot_warm_config     = 'data/COSMOS_ocean_examples/COSMOS_PI_oceanTS_prep.nc'
    CHARACTER(LEN=256)  :: filename_GCM_ocean_snapshot_cold_config     = 'data/COSMOS_ocean_examples/COSMOS_LGM_oceanTS_prep.nc'
    CHARACTER(LEN=256)  :: name_ocean_temperature_GCM_config           = 't_ocean'
    CHARACTER(LEN=256)  :: name_ocean_salinity_GCM_config              = 's_ocean'

    ! Uniform ocean temperature values used when choice_ocean_model = "uniform_warm_cold"
    REAL(dp)            :: ocean_temperature_PD_config                 = 271.46_dp                        ! present day temperature of the ocean beneath the shelves [K; -1.7 Celsius]
    REAL(dp)            :: ocean_temperature_cold_config               = 268.16_dp                        ! cold period temperature of the ocean beneath the shelves [K; -5.0 Celcius]
    REAL(dp)            :: ocean_temperature_warm_config               = 275.16_dp                        ! warm period temperature of the ocean beneath the shelves [K;  2.0 Celcius]

    ! Parameters used when choice_ocean_model = "matrix_warm_cold"
    CHARACTER(LEN=256)  :: choice_ocean_vertical_grid_config           = 'regular'                        ! Choice of vertical grid to be used for ocean data
    REAL(dp)            :: ocean_vertical_grid_max_depth_config        = 1500._dp                         ! Maximum depth           to be used for ocean data
    REAL(dp)            :: ocean_regular_grid_dz_config                = 150._dp                          ! Vertical grid spacing   to be used for ocean data when choice_ocean_vertical_grid_config = 'regular'
    CHARACTER(LEN=256)  :: ocean_extrap_dir_config                     = 'data/extrapolated_ocean_files'  ! Directory where extrapolated ocean files are stored
    REAL(dp)            :: ocean_extrap_res_config                     = 5000._dp                         ! High resolution at which the ocean data extrapolation should be performed
    REAL(dp)            :: ocean_extrap_Gauss_sigma_config             = 8000._dp                         ! 1-sigma of the Gaussian smoothing operation used to extrapolate the ocean data
    CHARACTER(LEN=256)  :: ocean_extrap_hires_geo_filename_NAM_config  = 'data/ETOPO1/NorthAmerica_ETOPO1_5km.nc'                     ! Path to a NetCDF file containing
    CHARACTER(LEN=256)  :: ocean_extrap_hires_geo_filename_EAS_config  = 'data/ETOPO1/Eurasia_ETOPO1_5km.nc'                          ! (present-day) geometry at high
    CHARACTER(LEN=256)  :: ocean_extrap_hires_geo_filename_GRL_config  = 'data/Bedmachine_Greenland/BedMachine_Greenland_v4_5km.nc'      ! resolution, used for ocean
    CHARACTER(LEN=256)  :: ocean_extrap_hires_geo_filename_ANT_config  = 'data/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc' ! data extrapolation
    REAL(dp)            :: ocean_w_tot_hist_averaging_window_config    = 1500._dp                         ! Time window (in yr) over which the weighing fields for sea-water temperature at maximum depth are averaged

  ! == Surface mass balance
  ! =======================

    CHARACTER(LEN=256)  :: choice_SMB_model_config                     = 'IMAU-ITM'                       ! Choice of SMB model: "uniform", "idealised", "IMAU-ITM", "direct_global", "direct_regional"
    CHARACTER(LEN=256)  :: choice_idealised_SMB_config                 = 'EISMINT1_A'
    REAL(dp)            :: SMB_uniform_config                          = 0._dp                            ! Uniform SMB, applied when choice_SMB_model = "uniform" [mie/yr]

    ! NetCDF file containing direct global/regional climate forcing
    CHARACTER(LEN=256)  :: filename_direct_global_SMB_config           = ''
    CHARACTER(LEN=256)  :: filename_direct_regional_SMB_NAM_config     = ''
    CHARACTER(LEN=256)  :: filename_direct_regional_SMB_EAS_config     = ''
    CHARACTER(LEN=256)  :: filename_direct_regional_SMB_GRL_config     = ''
    CHARACTER(LEN=256)  :: filename_direct_regional_SMB_ANT_config     = ''

    ! Firn layer
    CHARACTER(LEN=256)  :: SMB_IMAUITM_choice_init_firn_NAM_config     = 'uniform'                        ! How to initialise the firn layer in the IMAU-ITM SMB model: "uniform", "restart"
    CHARACTER(LEN=256)  :: SMB_IMAUITM_choice_init_firn_EAS_config     = 'uniform'
    CHARACTER(LEN=256)  :: SMB_IMAUITM_choice_init_firn_GRL_config     = 'uniform'
    CHARACTER(LEN=256)  :: SMB_IMAUITM_choice_init_firn_ANT_config     = 'uniform'
    REAL(dp)            :: SMB_IMAUITM_initial_firn_thickness_config   = 1._dp                            ! Initial firn thickness of the IMAU-ITEM SMB model [m] (used when SMB_IMAUITM_choice_init_firn = "uniform")

    ! IMAU-ITM SMB model parameters                                                                       ! Commented values are Tijn` calibration against RACMO
    REAL(dp)            :: SMB_IMAUITM_C_abl_constant_NAM_config       = 0._dp                            ! ??.? : Homogeneous-reduction factor during melt computation
    REAL(dp)            :: SMB_IMAUITM_C_abl_constant_EAS_config       = 0._dp
    REAL(dp)            :: SMB_IMAUITM_C_abl_constant_GRL_config       = 0._dp
    REAL(dp)            :: SMB_IMAUITM_C_abl_constant_ANT_config       = 0._dp
    REAL(dp)            :: SMB_IMAUITM_C_abl_Ts_NAM_config             = 10._dp                           ! 10.0 : Temperature-based melt factor
    REAL(dp)            :: SMB_IMAUITM_C_abl_Ts_EAS_config             = 10._dp
    REAL(dp)            :: SMB_IMAUITM_C_abl_Ts_GRL_config             = 10._dp
    REAL(dp)            :: SMB_IMAUITM_C_abl_Ts_ANT_config             = 10._dp
    REAL(dp)            :: SMB_IMAUITM_C_abl_Q_NAM_config              = 0.0227_dp                        ! 0.0227 : Insolation-based melt factor
    REAL(dp)            :: SMB_IMAUITM_C_abl_Q_EAS_config              = 0.0227_dp
    REAL(dp)            :: SMB_IMAUITM_C_abl_Q_GRL_config              = 0.0227_dp
    REAL(dp)            :: SMB_IMAUITM_C_abl_Q_ANT_config              = 0.0227_dp
    REAL(dp)            :: SMB_IMAUITM_C_refr_NAM_config               = 0.051_dp                         ! 0.051 : Temperature-meltwater-based refreezing factor
    REAL(dp)            :: SMB_IMAUITM_C_refr_EAS_config               = 0.051_dp
    REAL(dp)            :: SMB_IMAUITM_C_refr_GRL_config               = 0.051_dp
    REAL(dp)            :: SMB_IMAUITM_C_refr_ANT_config               = 0.051_dp

    ! IMAU-ITM SMB model inversion
    LOGICAL             :: do_SMB_IMAUITM_inversion_config             = .FALSE.                          ! If set to TRUE, basal roughness is iteratively adjusted to match initial ice thickness
    CHARACTER(LEN=256)  :: SMB_IMAUITM_inv_choice_init_C_config        = 'uniform'                        ! How to initialise the C parameters in the IMAU-ITM SMB inversion: "uniform", "restart"
    REAL(dp)            :: SMB_IMAUITM_inv_scale_config                = 10000._dp                        ! Scaling constant for inversion procedure [m]
    REAL(dp)            :: SMB_IMAUITM_inv_C_abl_constant_min_config   = -50._dp                          ! Minimum value of C_abl_constant allowed during inversion
    REAL(dp)            :: SMB_IMAUITM_inv_C_abl_constant_max_config   = 0._dp                            ! Maximum value of C_abl_constant allowed during inversion
    REAL(dp)            :: SMB_IMAUITM_inv_C_abl_Ts_min_config         = 0._dp                            ! Minimum value of C_abl_Ts       allowed during inversion
    REAL(dp)            :: SMB_IMAUITM_inv_C_abl_Ts_max_config         = 50._dp                           ! Maximum value of C_abl_Ts       allowed during inversion
    REAL(dp)            :: SMB_IMAUITM_inv_C_abl_Q_min_config          = 0._dp                            ! Minimum value of C_abl_Q        allowed during inversion
    REAL(dp)            :: SMB_IMAUITM_inv_C_abl_Q_max_config          = 1.0_dp                           ! Maximum value of C_abl_Q        allowed during inversion
    REAL(dp)            :: SMB_IMAUITM_inv_C_refr_min_config           = 0._dp                            ! Minimum value of C_refr         allowed during inversion
    REAL(dp)            :: SMB_IMAUITM_inv_C_refr_max_config           = 0.1_dp                           ! Maximum value of C_refr         allowed during inversion

  ! ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
  ! ==============================================================

    CHARACTER(LEN=256)  :: ISMIP_forcing_filename_baseline_config      = ''                              ! NetCDF file containing the baseline climate
    CHARACTER(LEN=256)  :: ISMIP_forcing_foldername_aSMB_config        = ''                              ! Folder containing the single-year NetCDF files of the SMB anomaly
    CHARACTER(LEN=256)  :: ISMIP_forcing_basefilename_aSMB_config      = ''                              ! Filename without the year (e.g. if the actual file is "aSMB_MARv3.12-yearly-CESM2-ssp585-1950.nc",   then this variable should be "aSMB_MARv3.12-yearly-CESM2-ssp585-"
    CHARACTER(LEN=256)  :: ISMIP_forcing_foldername_dSMBdz_config      = ''                              ! Folder containing the single-year NetCDF files of the SMB lapse rate
    CHARACTER(LEN=256)  :: ISMIP_forcing_basefilename_dSMBdz_config    = ''                              ! Filename without the year (e.g. if the actual file is "dSMBdz_MARv3.12-yearly-CESM2-ssp585-1950.nc", then this variable should be "dSMBdz_MARv3.12-yearly-CESM2-ssp585-"
    CHARACTER(LEN=256)  :: ISMIP_forcing_foldername_aST_config         = ''                              ! Folder containing the single-year NetCDF files of the temperature anomaly
    CHARACTER(LEN=256)  :: ISMIP_forcing_basefilename_aST_config       = ''                              ! Filename without the year (e.g. if the actual file is "aST_MARv3.12-yearly-CESM2-ssp585-1950.nc",    then this variable should be "aST_MARv3.12-yearly-CESM2-ssp585-"
    CHARACTER(LEN=256)  :: ISMIP_forcing_foldername_dSTdz_config       = ''                              ! Folder containing the single-year NetCDF files of the temperature lapse rate
    CHARACTER(LEN=256)  :: ISMIP_forcing_basefilename_dSTdz_config     = ''                              ! Filename without the year (e.g. if the actual file is "dSTdz_MARv3.12-yearly-CESM2-ssp585-1950.nc",  then this variable should be "dSTdz_MARv3.12-yearly-CESM2-ssp585-"

  ! == Basal mass balance
  ! =====================

    CHARACTER(LEN=256)  :: choice_BMB_shelf_model_config               = 'Favier2019_quad'                ! Choice of shelf BMB: "uniform", "idealised", "ANICE_legacy", "Favier2019_lin", "Favier2019_quad", "Favier2019_Mplus", "Lazeroms2018_plume", "PICO", "PICOP", 'inversion'
    CHARACTER(LEN=256)  :: choice_idealised_BMB_shelf_config           = 'MISMIP+'
    CHARACTER(LEN=256)  :: choice_BMB_sheet_model_config               = 'uniform'                        ! Choice of sheet BMB: "uniform"
    REAL(dp)            :: BMB_shelf_uniform_config                    = 0._dp                            ! Uniform shelf BMB, applied when choice_BMB_shelf_model = "uniform" [mie/yr]
    REAL(dp)            :: BMB_sheet_uniform_config                    = 0._dp                            ! Uniform sheet BMB, applied when choice_BMB_sheet_model = "uniform" [mie/yr]
    CHARACTER(LEN=256)  :: choice_BMB_subgrid_config                   = 'FCMP'                           ! Choice of sub-grid BMB scheme: "FCMP", "PMP", "NMP" (following Leguy et al., 2021)
    REAL(dp)            :: BMB_max_config                              = 20._dp                           ! Maximum amount of allowed basal mass balance [mie/yr] (+ is refreezing)
    REAL(dp)            :: BMB_min_config                              = -200._dp                         ! Minimum amount of allowed basal mass balance [mie/yr] (- is melting)

    LOGICAL             :: do_ocean_temperature_inversion_config       = .FALSE.                          ! Whether or not to invert for ocean temperatures in the Bernales202X shelf BMB model
    REAL(dp)            :: ocean_temperature_inv_t_start_config        = -9.9E9_dp                        ! Minimum model time when the inversion is allowed
    REAL(dp)            :: ocean_temperature_inv_t_end_config          = +9.9E9_dp                        ! Maximum model time when the inversion is allowed
    INTEGER             :: T_base_window_size_config                   = 200                              ! Number of previous time steps used to compute a running average of inverted T_ocean_base

    LOGICAL             :: BMB_inv_use_restart_field_config            = .FALSE.                          ! Whether or not to use BMB_shelf field from a the restart file
    REAL(dp)            :: BMB_inv_scale_shelf_config                  = 200._dp                          ! Scaling constant for inversion procedure over shelves [m]
    REAL(dp)            :: BMB_inv_scale_ocean_config                  = 100._dp                          ! Scaling constant for inversion procedure over open ocean [m]

    CHARACTER(LEN=256)  :: choice_basin_scheme_NAM_config              = 'none'                           ! Choice of basin ID scheme; can be 'none' or 'file'
    CHARACTER(LEN=256)  :: choice_basin_scheme_EAS_config              = 'none'
    CHARACTER(LEN=256)  :: choice_basin_scheme_GRL_config              = 'file'
    CHARACTER(LEN=256)  :: choice_basin_scheme_ANT_config              = 'none'
    CHARACTER(LEN=256)  :: filename_basins_NAM_config                  = ''                               ! Path to a text file containing polygons of drainage basins
    CHARACTER(LEN=256)  :: filename_basins_EAS_config                  = ''
    CHARACTER(LEN=256)  :: filename_basins_GRL_config                  = 'data/drainage_basins/grndrainagesystems_ekholm.txt'
    CHARACTER(LEN=256)  :: filename_basins_ANT_config                  = 'data/drainage_basins/ant_full_drainagesystem_polygons.txt'
    LOGICAL             :: do_merge_basins_ANT_config                  = .TRUE.                           ! Whether or not to merge some of the Antarctic basins
    LOGICAL             :: do_merge_basins_GRL_config                  = .TRUE.                           ! Whether or not to merge some of the Greenland basins

    CHARACTER(LEN=256)       ::  choice_BMB_shelf_amplification_config        = 'uniform'                 ! Choice of method to determine BMB amplification factors: "uniform", "basin"
    INTEGER                  ::  basin_BMB_amplification_n_ANT_config         = 17                        ! Number of basins used for ANT
    REAL(dp), DIMENSION(17)  ::  basin_BMB_amplification_factor_ANT_config    = &                         ! BMB amplification factor for each basin for ANT
    (/ 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &
       1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp /)
    INTEGER                  ::  basin_BMB_amplification_n_GRL_config         = 8                         ! Number of basins used for GRL
    REAL(dp), DIMENSION(8)   ::  basin_BMB_amplification_factor_GRL_config    = &                         ! BMB amplification factor for each basin for GRL
    (/ 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp /)

    ! Parameters for the three simple melt parameterisations from Favier et al. (2019)
    REAL(dp)            :: BMB_Favier2019_lin_GammaT_config            = 3.3314E-5_dp   ! 2.03E-5_dp      ! Heat exchange velocity [m s^-1]
    REAL(dp)            :: BMB_Favier2019_quad_GammaT_config           = 111.6E-5_dp    ! 99.32E-5_dp     ! Commented values are from Favier et al. (2019), Table 3
    REAL(dp)            :: BMB_Favier2019_Mplus_GammaT_config          = 108.6E-5_dp    ! 132.9E-5_dp     ! Actual value are re-tuned for IMAU-ICE, following the same approach (see Asay-Davis et al., 2016, ISOMIP+)

    ! Parameters for the Lazeroms et al. (2018) plume-parameterisation BMB model
    REAL(dp)            :: BMB_Lazeroms2018_GammaT_config              = 3.7506E-04_dp  ! 1.1E-3_dp       ! Thermal exchange velocity; tuned following ISOMIP+ protocol (Asay-Davis et al., 2016, Sect. 3.2.1), commented value from Lazeroms et al. (2018)
    CHARACTER(LEN=256)  :: BMB_Lazeroms2018_find_GL_scheme_config      = 'along_ice_flow'                 ! How to determine the GL origin of a plume: "GL_average", "along_ice_flow"

    ! Parameters for the PICO BMB model
    INTEGER             :: BMB_PICO_nboxes_config                      = 5                                ! Number of sub-shelf ocean boxes used by PICO
    REAL(dp)            :: BMB_PICO_GammaTstar_config                  = 3.6131E-05_dp  ! 2.0E-5_dp       ! Effective turbulent temperature exchange velocity [m s^-1]; tuned following ISOMIP+ protocol (Asay-Davis et al., 2016, Sect. 3.2.1), commented value from Reese et al. (2018)

    ! Parameters for the ANICE_legacy sub-shelf melt model
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

  ! == Englacial isotope tracing
  ! ============================

    CHARACTER(LEN=256)  :: choice_ice_isotopes_model_config             = 'ANICE_legacy'                  ! Choice of englacial isotopes model: "none", "uniform", "ANICE_legacy"
    REAL(dp)            :: uniform_ice_d18O_config                      = 0._dp                           ! Uniform englacial d18O-value (used when choice_ice_isotope_model_config = "uniform")

  ! == Sea level and GIA
  ! ====================

    LOGICAL             :: do_ocean_floodfill_config                   = .TRUE.                           ! Use a flood-fill to determine the ocean mask, so that (pro-/sub-glacial) lakes dont exist
    CHARACTER(LEN=256)  :: choice_sealevel_model_config                = 'eustatic'                       ! Can be "fixed", "prescribed", "eustatic", or "SELEN"
    REAL(dp)            :: fixed_sealevel_config                       = 0._dp                            ! Sea level for "fixed" method
    REAL(dp)            :: initial_guess_sealevel_config               = 0._dp                            ! Initial sea-level guess value for "eustatic" and "SELEN" methods
    CHARACTER(LEN=256)  :: filename_sealevel_record_config             = 'name_of_file.dat'
    INTEGER             :: sealevel_record_length_config               = 1

    CHARACTER(LEN=256)  :: choice_GIA_model_config                     = 'ELRA'                           ! Can be "none", "ELRA", or "SELEN"
    REAL(dp)            :: ELRA_lithosphere_flex_rigidity_config       = 1.0E+25_dp                       ! Lithospheric flexural rigidity [kg m^2 s^-2]
    REAL(dp)            :: ELRA_bedrock_relaxation_time_config         = 3000.0_dp                        ! Relaxation time for bedrock adjustment [yr]
    REAL(dp)            :: ELRA_mantle_density_config                  = 3300.0_dp                        ! Mantle density [kg m^-3]

  ! == SELEN
  ! ========

    LOGICAL             :: SELEN_run_at_t_start_config                  = .FALSE.                         ! Whether or not to run SELEN in the first coupling loop (needed for some benchmark experiments)
    INTEGER             :: SELEN_n_TDOF_iterations_config               = 1                               ! Number of Time-Dependent Ocean Function iterations
    INTEGER             :: SELEN_n_recursion_iterations_config          = 1                               ! Number of recursion iterations
    LOGICAL             :: SELEN_use_rotational_feedback_config         = .FALSE.                         ! If TRUE, rotational feedback is included
    INTEGER             :: SELEN_n_harmonics_config                     = 128                             ! Maximum number of harmonic degrees
    LOGICAL             :: SELEN_display_progress_config                = .FALSE.                         ! Whether or not to display the progress of the big loops to the screen (doesn't work on Cartesius!)

    CHARACTER(LEN=256)  :: SELEN_dir_config                             = 'data/SELEN'                    ! Directory where SELEN initial files and spherical harmonics are stored
    CHARACTER(LEN=256)  :: SELEN_global_topo_filename_config            = 'SELEN_global_topography.nc'    ! Filename for the SELEN global topography file (located in SELEN_dir)
    CHARACTER(LEN=256)  :: SELEN_TABOO_init_filename_config             = 'SELEN_TABOO_initial_file.dat'  ! Filename for the TABOO initial file           (idem                )
    CHARACTER(LEN=256)  :: SELEN_LMJ_VALUES_filename_config             = 'SELEN_lmj_values.bin'          ! Filename for the LJ and MJ values file        (idem                )

    INTEGER                  :: SELEN_irreg_time_n_config               = 15                              ! Number of entries in the irregular moving time window
    REAL(dp), DIMENSION(50)  :: SELEN_irreg_time_window_config          = &                               ! Values of entries in the irregular moving time window
   (/20._dp, 20._dp, 20._dp, 5._dp, 5._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp  /)

    REAL(dp)            :: SELEN_lith_thickness_config                  = 100._dp                         ! Thickness of the elastic lithosphere [km]
    INTEGER             :: SELEN_visc_n_config                          = 3                               ! Number      of viscous asthenosphere layers
    REAL(dp), DIMENSION(3) :: SELEN_visc_prof_config                    = (/ 3._dp, 0.6_dp, 0.3_dp /)     ! Viscosities of viscous asthenosphere layers [?]

    ! Settings for the TABOO Earth deformation model
    INTEGER             :: SELEN_TABOO_CDE_config                       = 0                               ! code of the model (see taboo for explanation)
    INTEGER             :: SELEN_TABOO_TLOVE_config                     = 1                               ! Tidal love numbers yes/no
    INTEGER             :: SELEN_TABOO_DEG1_config                      = 1                               ! Tidal love numbers degree
    REAL(dp)            :: SELEN_TABOO_RCMB_config                      = 3480._dp                        ! Radius of CMB (km)

  ! == Which data fields will be written to the help_fields output file
  ! ===================================================================

    CHARACTER(LEN=256)  :: help_field_01_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_02_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_03_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_04_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_05_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_06_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_07_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_08_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_09_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_10_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_11_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_12_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_13_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_14_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_15_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_16_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_17_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_18_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_19_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_20_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_21_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_22_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_23_config                        = 'none'
    CHARACTER(LEN=256)  :: help_field_24_config                        = 'none'
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

! ===== The C type =====
! ======================
  ! The "C" type, which contains all the config parameters as fields.
  ! These will all be overwritten with the values of the "_config" variables,
  ! which are either the default values specified above, are the values
  ! specified from the external config file.

  TYPE constants_type

    ! Time steps and range
    ! =====================

    REAL(dp)                            :: start_time_of_run
    REAL(dp)                            :: end_time_of_run
    REAL(dp)                            :: dt_coupling
    REAL(dp)                            :: dt_max
    REAL(dp)                            :: dt_min
    REAL(dp)                            :: dt_startup_phase
    REAL(dp)                            :: dt_cooldown_phase
    REAL(dp)                            :: dt_thermo
    REAL(dp)                            :: dt_climate
    REAL(dp)                            :: dt_ocean
    REAL(dp)                            :: dt_SMB
    REAL(dp)                            :: dt_BMB
    REAL(dp)                            :: dt_output
    REAL(dp)                            :: dt_mesh_min
    REAL(dp)                            :: dt_bedrock_ELRA
    REAL(dp)                            :: dt_SELEN
    REAL(dp)                            :: dt_slid_inv
    REAL(dp)                            :: dt_SMB_inv

    ! Which ice sheets do we simulate?
    ! ================================

    LOGICAL                             :: do_NAM
    LOGICAL                             :: do_EAS
    LOGICAL                             :: do_GRL
    LOGICAL                             :: do_ANT

    ! Benchmark experiments
    ! =====================

    REAL(dp)                            :: SSA_icestream_m
    REAL(dp)                            :: ISMIP_HOM_L
    CHARACTER(LEN=256)                  :: ISMIP_HOM_E_Arolla_filename
    LOGICAL                             :: MISMIPplus_do_tune_A_for_GL
    REAL(dp)                            :: MISMIPplus_xGL_target
    REAL(dp)                            :: MISMIPplus_A_flow_initial
    CHARACTER(LEN=256)                  :: MISMIPplus_scenario
    CHARACTER(LEN=256)                  :: MISOMIP1_scenario

    ! Whether or not to let UFEMISM dynamically create its own output folder
    ! =======================================================================

    LOGICAL                             :: create_procedural_output_dir
    CHARACTER(LEN=256)                  :: fixed_output_dir
    CHARACTER(LEN=256)                  :: fixed_output_dir_suffix
    LOGICAL                             :: do_write_regional_scalar_output
    LOGICAL                             :: do_write_global_scalar_output

    ! Debugging
    ! =========

    LOGICAL                             :: do_write_debug_data
    LOGICAL                             :: do_write_grid_data
    LOGICAL                             :: do_check_for_NaN
    LOGICAL                             :: do_time_display

    ! == The four model regions
    ! =========================

    ! North America
    REAL(dp)                            :: lambda_M_NAM
    REAL(dp)                            :: phi_M_NAM
    REAL(dp)                            :: beta_stereo_NAM
    REAL(dp)                            :: xmin_NAM
    REAL(dp)                            :: xmax_NAM
    REAL(dp)                            :: ymin_NAM
    REAL(dp)                            :: ymax_NAM

    ! Eurasia
    REAL(dp)                            :: lambda_M_EAS
    REAL(dp)                            :: phi_M_EAS
    REAL(dp)                            :: beta_stereo_EAS
    REAL(dp)                            :: xmin_EAS
    REAL(dp)                            :: xmax_EAS
    REAL(dp)                            :: ymin_EAS
    REAL(dp)                            :: ymax_EAS

    ! Greenland
    REAL(dp)                            :: lambda_M_GRL
    REAL(dp)                            :: phi_M_GRL
    REAL(dp)                            :: beta_stereo_GRL
    REAL(dp)                            :: xmin_GRL
    REAL(dp)                            :: xmax_GRL
    REAL(dp)                            :: ymin_GRL
    REAL(dp)                            :: ymax_GRL

    ! Antarctica
    REAL(dp)                            :: lambda_M_ANT
    REAL(dp)                            :: phi_M_ANT
    REAL(dp)                            :: beta_stereo_ANT
    REAL(dp)                            :: xmin_ANT
    REAL(dp)                            :: xmax_ANT
    REAL(dp)                            :: ymin_ANT
    REAL(dp)                            :: ymax_ANT

    ! Mesh generation parameters
    ! ==========================

    LOGICAL                             :: use_submesh
    INTEGER                             :: nconmax
    REAL(dp)                            :: alpha_min
    REAL(dp)                            :: dz_max_ice
    REAL(dp)                            :: res_max
    REAL(dp)                            :: res_min
    REAL(dp)                            :: res_max_ice
    REAL(dp)                            :: res_max_margin
    REAL(dp)                            :: res_max_gl
    REAL(dp)                            :: res_max_cf
    REAL(dp)                            :: res_max_mountain
    REAL(dp)                            :: res_max_coast
    REAL(dp)                            :: mesh_fitness_threshold
    LOGICAL                             :: do_force_mesh_update
    REAL(dp)                            :: do_force_mesh_update_after

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

    ! ISMIP-style output
    ! ==================

    LOGICAL                             :: do_write_ISMIP_output
    CHARACTER(LEN=256)                  :: ISMIP_output_group_code
    CHARACTER(LEN=256)                  :: ISMIP_output_model_code
    CHARACTER(LEN=256)                  :: ISMIP_output_experiment_code
    CHARACTER(LEN=256)                  :: ISMIP_output_basetime

    ! Scaled vertical coordinate zeta
    ! ===============================

    INTEGER                             :: nz       ! Number of grid points in vertical direction for thermodynamics in ice sheet
    REAL(dp), DIMENSION(:), ALLOCATABLE :: zeta

    ! Reference geometries (initial, present-day, and GIA equilibrium)
    ! ================================================================

    ! Some pre-processing stuff for reference ice geometry
    REAL(dp)                            :: refgeo_Hi_min
    LOGICAL                             :: remove_Lake_Vostok

    ! Initial geometry
    CHARACTER(LEN=256)                  :: choice_refgeo_init_NAM
    CHARACTER(LEN=256)                  :: choice_refgeo_init_EAS
    CHARACTER(LEN=256)                  :: choice_refgeo_init_GRL
    CHARACTER(LEN=256)                  :: choice_refgeo_init_ANT
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

    ! Whether or not the simulation is a restart of a previous simulation
    ! ===================================================================

    LOGICAL                             :: is_restart

    ! Time in the restart file from which the data will be extracted
    REAL(dp)                            :: time_to_restart_from_NAM
    REAL(dp)                            :: time_to_restart_from_EAS
    REAL(dp)                            :: time_to_restart_from_GRL
    REAL(dp)                            :: time_to_restart_from_ANT

    ! Initial model state when restarting from a previous run
    CHARACTER(LEN=256)                  :: filename_restart_NAM
    CHARACTER(LEN=256)                  :: filename_restart_EAS
    CHARACTER(LEN=256)                  :: filename_restart_GRL
    CHARACTER(LEN=256)                  :: filename_restart_ANT

    ! Global forcing (insolation, CO2, d18O, geothermal heat flux)
    ! ============================================================

    CHARACTER(LEN=256)                  :: choice_forcing_method

    ! Insolation forcing (NetCDF)
    CHARACTER(LEN=256)                  :: choice_insolation_forcing
    REAL(dp)                            :: static_insolation_time
    CHARACTER(LEN=256)                  :: filename_insolation

    ! CO2 record (ASCII text file, so the number of rows needs to be specified)
    CHARACTER(LEN=256)                  :: filename_CO2_record
    INTEGER                             :: CO2_record_length

    ! d18O record (ASCII text file, so the number of rows needs to be specified)
    CHARACTER(LEN=256)                  :: filename_d18O_record
    INTEGER                             :: d18O_record_length

    ! Geothermal heat flux
    CHARACTER(LEN=256)                  :: choice_geothermal_heat_flux
    REAL(dp)                            :: constant_geothermal_heat_flux
    CHARACTER(LEN=256)                  :: filename_geothermal_heat_flux

    ! Parameters for calculating modelled benthic d18O
    LOGICAL                             :: do_calculate_benthic_d18O
    REAL(dp)                            :: dT_deepwater_averaging_window
    REAL(dp)                            :: dT_deepwater_dT_surf_ratio
    REAL(dp)                            :: d18O_dT_deepwater_ratio

    ! Parameters for the inverse routine
    REAL(dp)                            :: dT_glob_inverse_averaging_window
    REAL(dp)                            :: inverse_d18O_to_dT_glob_scaling
    REAL(dp)                            :: CO2_inverse_averaging_window
    REAL(dp)                            :: inverse_d18O_to_CO2_scaling
    REAL(dp)                            :: inverse_d18O_to_CO2_initial_CO2

    ! Ice dynamics - velocity
    ! =======================

    CHARACTER(LEN=256)                  :: choice_ice_dynamics
    REAL(dp)                            :: n_flow
    REAL(dp)                            :: m_enh_sheet
    REAL(dp)                            :: m_enh_shelf
    CHARACTER(LEN=256)                  :: choice_ice_margin
    LOGICAL                             :: do_hybrid_Bernales2017
    REAL(dp)                            :: vel_ref_Bernales2017
    LOGICAL                             :: include_SSADIVA_crossterms
    LOGICAL                             :: do_GL_subgrid_friction
    LOGICAL                             :: do_smooth_geometry
    REAL(dp)                            :: r_smooth_geometry

    ! Some parameters for numerically solving the SSA/DIVA
    REAL(dp)                            :: DIVA_visc_it_norm_dUV_tol
    INTEGER                             :: DIVA_visc_it_nit
    REAL(dp)                            :: DIVA_visc_it_relax
    REAL(dp)                            :: DIVA_epsilon_sq_0
    REAL(dp)                            :: DIVA_visc_eff_min
    REAL(dp)                            :: DIVA_beta_max
    REAL(dp)                            :: DIVA_vel_max
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

    ! Velocity wind-up before a restart
    REAL(dp)                            :: windup_total_years

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

    ! Ice thickness boundary conditions
    CHARACTER(LEN=256)                  :: ice_thickness_west_BC
    CHARACTER(LEN=256)                  :: ice_thickness_east_BC
    CHARACTER(LEN=256)                  :: ice_thickness_south_BC
    CHARACTER(LEN=256)                  :: ice_thickness_north_BC
    CHARACTER(LEN=256)                  :: choice_mask_noice_NAM
    CHARACTER(LEN=256)                  :: choice_mask_noice_EAS
    CHARACTER(LEN=256)                  :: choice_mask_noice_GRL
    CHARACTER(LEN=256)                  :: choice_mask_noice_ANT

    ! Partially fixed geometry, useful for initialisation and inversion runs
    REAL(dp)                            :: fixed_sheet_geometry
    REAL(dp)                            :: fixed_shelf_geometry
    REAL(dp)                            :: fixed_grounding_line_g
    REAL(dp)                            :: fixed_grounding_line_f

    ! Memory of previous run during a restart
    INTEGER                             :: dHi_dt_window_size

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
    CHARACTER(LEN=256)                  :: basal_roughness_restart_type
    LOGICAL                             :: do_basal_roughness_remap_adjustment

    ! Bed roughness inversion
    LOGICAL                             :: do_slid_inv
    CHARACTER(LEN=256)                  :: choice_slid_inv_method
    REAL(dp)                            :: slid_inv_t_start
    REAL(dp)                            :: slid_inv_t_end
    REAL(dp)                            :: slid_inv_phi_min
    REAL(dp)                            :: slid_inv_phi_max
    CHARACTER(LEN=256)                  :: slid_inv_filename_output
    INTEGER                             :: slid_inv_window_size

    LOGICAL                             :: do_slid_inv_Bernales2017_smooth
    LOGICAL                             :: do_slid_inv_Bernales2017_extrap
    REAL(dp)                            :: slid_inv_Bernales2017_hi_scale
    REAL(dp)                            :: slid_inv_Bernales2017_smooth_r
    REAL(dp)                            :: slid_inv_Bernales2017_smooth_w
    REAL(dp)                            :: slid_inv_Bernales2017_tol_diff
    REAL(dp)                            :: slid_inv_Bernales2017_tol_frac

    REAL(dp)                            :: slid_inv_Berends2022_tauc
    REAL(dp)                            :: slid_inv_Berends2022_H0
    REAL(dp)                            :: slid_inv_Berends2022_u0
    REAL(dp)                            :: slid_inv_Berends2022_Hi_scale
    REAL(dp)                            :: slid_inv_Berends2022_u_scale
    CHARACTER(LEN=256)                  :: slid_inv_target_velocity_filename

    ! Ice dynamics - calving
    ! ======================

    CHARACTER(LEN=256)                  :: choice_calving_law
    REAL(dp)                            :: calving_threshold_thickness_shelf
    REAL(dp)                            :: calving_threshold_thickness_sheet
    INTEGER                             :: max_calving_rounds
    LOGICAL                             :: do_remove_shelves
    LOGICAL                             :: remove_shelves_larger_than_PD
    LOGICAL                             :: continental_shelf_calving
    REAL(dp)                            :: continental_shelf_min_height
    REAL(dp)                            :: minimum_ice_thickness

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

    ! Climate
    ! =======

    CHARACTER(LEN=256)                  :: choice_climate_model
    CHARACTER(LEN=256)                  :: choice_idealised_climate

    ! NetCDF files containing direct global/regional climate forcing
    CHARACTER(LEN=256)                  :: filename_direct_global_climate
    CHARACTER(LEN=256)                  :: filename_direct_regional_climate_NAM
    CHARACTER(LEN=256)                  :: filename_direct_regional_climate_EAS
    CHARACTER(LEN=256)                  :: filename_direct_regional_climate_GRL
    CHARACTER(LEN=256)                  :: filename_direct_regional_climate_ANT

    ! NetCDF file containing the present-day observed climate (e.g. ERA40)
    CHARACTER(LEN=256)                  :: filename_PD_obs_climate

    ! GCM snapshots in the matrix_warm_cold option
    CHARACTER(LEN=256)                  :: filename_climate_snapshot_PI
    CHARACTER(LEN=256)                  :: filename_climate_snapshot_warm
    CHARACTER(LEN=256)                  :: filename_climate_snapshot_cold

    REAL(dp)                            :: constant_lapserate

    ! Orbit time and CO2 concentration of the warm and cold snapshots
    REAL(dp)                            :: matrix_high_CO2_level
    REAL(dp)                            :: matrix_low_CO2_level
    REAL(dp)                            :: matrix_warm_orbit_time
    REAL(dp)                            :: matrix_cold_orbit_time

    ! Whether or not to apply a bias correction to the GCM snapshots
    LOGICAL                             :: climate_matrix_biascorrect_warm
    LOGICAL                             :: climate_matrix_biascorrect_cold

    ! Glacial index
    LOGICAL                             :: switch_glacial_index

    ! Iterative adjustment of precipitation and temperature
    LOGICAL                             :: do_clim_inv
    REAL(dp)                            :: clim_inv_t_start
    REAL(dp)                            :: clim_inv_t_end
    REAL(dp)                            :: clim_inv_Hb_min
    REAL(dp)                            :: clim_inv_Hi_max
    INTEGER                             :: clim_inv_window_size

    ! Ocean
    ! =====

    CHARACTER(LEN=256)                  :: choice_ocean_model
    CHARACTER(LEN=256)                  :: choice_idealised_ocean

    ! NetCDF file containing the present-day observed ocean (WOA18) (NetCDF)
    CHARACTER(LEN=256)                  :: filename_PD_obs_ocean
    CHARACTER(LEN=256)                  :: name_ocean_temperature_obs
    CHARACTER(LEN=256)                  :: name_ocean_salinity_obs
    LOGICAL                             :: use_inverted_ocean

    ! GCM snapshots in the matrix_warm_cold option
    CHARACTER(LEN=256)                  :: filename_GCM_ocean_snapshot_PI
    CHARACTER(LEN=256)                  :: filename_GCM_ocean_snapshot_warm
    CHARACTER(LEN=256)                  :: filename_GCM_ocean_snapshot_cold
    CHARACTER(LEN=256)                  :: name_ocean_temperature_GCM
    CHARACTER(LEN=256)                  :: name_ocean_salinity_GCM

    ! Uniform ocean temperature values used when choice_ocean_model = "uniform_warm_cold"
    REAL(dp)                            :: ocean_temperature_PD
    REAL(dp)                            :: ocean_temperature_cold
    REAL(dp)                            :: ocean_temperature_warm

    ! Parameters used when choice_ocean_model = "matrix_warm_cold"
    CHARACTER(LEN=256)                  :: choice_ocean_vertical_grid
    REAL(dp)                            :: ocean_vertical_grid_max_depth
    REAL(dp)                            :: ocean_regular_grid_dz
    INTEGER                             :: nz_ocean ! NOTE: nz_ocean and z_ocean cannot be set through the config file, but are filled in by the "initialise_ocean_vertical_grid" in the ocean_module!
    REAL(dp), DIMENSION(:), ALLOCATABLE :: z_ocean
    CHARACTER(LEN=256)                  :: ocean_extrap_dir
    REAL(dp)                            :: ocean_extrap_res
    REAL(dp)                            :: ocean_extrap_Gauss_sigma
    CHARACTER(LEN=256)                  :: ocean_extrap_hires_geo_filename_NAM
    CHARACTER(LEN=256)                  :: ocean_extrap_hires_geo_filename_EAS
    CHARACTER(LEN=256)                  :: ocean_extrap_hires_geo_filename_GRL
    CHARACTER(LEN=256)                  :: ocean_extrap_hires_geo_filename_ANT
    REAL(dp)                            :: ocean_w_tot_hist_averaging_window

    ! Surface mass balance
    ! ====================

    CHARACTER(LEN=256)                  :: choice_SMB_model
    CHARACTER(LEN=256)                  :: choice_idealised_SMB
    REAL(dp)                            :: SMB_uniform

    ! NetCDF file containing direct global/regional SMB forcing
    CHARACTER(LEN=256)                  :: filename_direct_global_SMB
    CHARACTER(LEN=256)                  :: filename_direct_regional_SMB_NAM
    CHARACTER(LEN=256)                  :: filename_direct_regional_SMB_EAS
    CHARACTER(LEN=256)                  :: filename_direct_regional_SMB_GRL
    CHARACTER(LEN=256)                  :: filename_direct_regional_SMB_ANT

    ! Tuning parameters for the IMAU-ITM SMB model
    CHARACTER(LEN=256)                  :: SMB_IMAUITM_choice_init_firn_NAM
    CHARACTER(LEN=256)                  :: SMB_IMAUITM_choice_init_firn_EAS
    CHARACTER(LEN=256)                  :: SMB_IMAUITM_choice_init_firn_GRL
    CHARACTER(LEN=256)                  :: SMB_IMAUITM_choice_init_firn_ANT
    REAL(dp)                            :: SMB_IMAUITM_initial_firn_thickness
    REAL(dp)                            :: SMB_IMAUITM_C_abl_constant_NAM
    REAL(dp)                            :: SMB_IMAUITM_C_abl_constant_EAS
    REAL(dp)                            :: SMB_IMAUITM_C_abl_constant_GRL
    REAL(dp)                            :: SMB_IMAUITM_C_abl_constant_ANT
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Ts_NAM
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Ts_EAS
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Ts_GRL
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Ts_ANT
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Q_NAM
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Q_EAS
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Q_GRL
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Q_ANT
    REAL(dp)                            :: SMB_IMAUITM_C_refr_NAM
    REAL(dp)                            :: SMB_IMAUITM_C_refr_EAS
    REAL(dp)                            :: SMB_IMAUITM_C_refr_GRL
    REAL(dp)                            :: SMB_IMAUITM_C_refr_ANT

    ! IMAU-ITM SMB model inversion
    LOGICAL                             :: do_SMB_IMAUITM_inversion
    CHARACTER(LEN=256)                  :: SMB_IMAUITM_inv_choice_init_C
    REAL(dp)                            :: SMB_IMAUITM_inv_scale
    REAL(dp)                            :: SMB_IMAUITM_inv_C_abl_constant_min
    REAL(dp)                            :: SMB_IMAUITM_inv_C_abl_constant_max
    REAL(dp)                            :: SMB_IMAUITM_inv_C_abl_Ts_min
    REAL(dp)                            :: SMB_IMAUITM_inv_C_abl_Ts_max
    REAL(dp)                            :: SMB_IMAUITM_inv_C_abl_Q_min
    REAL(dp)                            :: SMB_IMAUITM_inv_C_abl_Q_max
    REAL(dp)                            :: SMB_IMAUITM_inv_C_refr_min
    REAL(dp)                            :: SMB_IMAUITM_inv_C_refr_max

    ! ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
    ! ==============================================================

    CHARACTER(LEN=256)                  :: ISMIP_forcing_filename_baseline
    CHARACTER(LEN=256)                  :: ISMIP_forcing_foldername_aSMB
    CHARACTER(LEN=256)                  :: ISMIP_forcing_basefilename_aSMB
    CHARACTER(LEN=256)                  :: ISMIP_forcing_foldername_dSMBdz
    CHARACTER(LEN=256)                  :: ISMIP_forcing_basefilename_dSMBdz
    CHARACTER(LEN=256)                  :: ISMIP_forcing_foldername_aST
    CHARACTER(LEN=256)                  :: ISMIP_forcing_basefilename_aST
    CHARACTER(LEN=256)                  :: ISMIP_forcing_foldername_dSTdz
    CHARACTER(LEN=256)                  :: ISMIP_forcing_basefilename_dSTdz

    ! Basal mass balance - sub-shelf melt
    ! ===================================

    CHARACTER(LEN=256)                  :: choice_BMB_shelf_model
    CHARACTER(LEN=256)                  :: choice_idealised_BMB_shelf
    CHARACTER(LEN=256)                  :: choice_BMB_sheet_model
    REAL(dp)                            :: BMB_shelf_uniform
    REAL(dp)                            :: BMB_sheet_uniform
    CHARACTER(LEN=256)                  :: choice_BMB_subgrid
    REAL(dp)                            :: BMB_max
    REAL(dp)                            :: BMB_min

    LOGICAL                             :: do_ocean_temperature_inversion
    REAL(dp)                            :: ocean_temperature_inv_t_start
    REAL(dp)                            :: ocean_temperature_inv_t_end
    INTEGER                             :: T_base_window_size

    LOGICAL                             :: BMB_inv_use_restart_field
    REAL(dp)                            :: BMB_inv_scale_shelf
    REAL(dp)                            :: BMB_inv_scale_ocean

    CHARACTER(LEN=256)                  :: choice_basin_scheme_NAM
    CHARACTER(LEN=256)                  :: choice_basin_scheme_EAS
    CHARACTER(LEN=256)                  :: choice_basin_scheme_GRL
    CHARACTER(LEN=256)                  :: choice_basin_scheme_ANT
    CHARACTER(LEN=256)                  :: filename_basins_NAM
    CHARACTER(LEN=256)                  :: filename_basins_EAS
    CHARACTER(LEN=256)                  :: filename_basins_GRL
    CHARACTER(LEN=256)                  :: filename_basins_ANT
    LOGICAL                             :: do_merge_basins_ANT
    LOGICAL                             :: do_merge_basins_GRL

    CHARACTER(LEN=256)                  :: choice_BMB_shelf_amplification
    INTEGER                             :: basin_BMB_amplification_n_ANT
    REAL(dp), DIMENSION(:), ALLOCATABLE :: basin_BMB_amplification_factor_ANT
    INTEGER                             :: basin_BMB_amplification_n_GRL
    REAL(dp), DIMENSION(:), ALLOCATABLE :: basin_BMB_amplification_factor_GRL

    ! Parameters for the three simple melt parameterisations from Favier et al. (2019)
    REAL(dp)                            :: BMB_Favier2019_lin_GammaT
    REAL(dp)                            :: BMB_Favier2019_quad_GammaT
    REAL(dp)                            :: BMB_Favier2019_Mplus_GammaT

    ! Parameters for the Lazeroms et al. (2018) plume-parameterisation BMB model
    REAL(dp)                            :: BMB_Lazeroms2018_GammaT
    CHARACTER(LEN=256)                  :: BMB_Lazeroms2018_find_GL_scheme

    ! Parameters for the PICO BMB model
    INTEGER                             :: BMB_PICO_nboxes
    REAL(dp)                            :: BMB_PICO_GammaTstar

    ! Parameters for the ANICE_legacy sub-shelf melt model
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

    ! Englacial isotope tracing
    ! ========================

    CHARACTER(LEN=256)                  :: choice_ice_isotopes_model
    REAL(dp)                            :: uniform_ice_d18O

    ! Sea level and GIA
    ! =================

    LOGICAL                             :: do_ocean_floodfill
    CHARACTER(LEN=256)                  :: choice_sealevel_model
    REAL(dp)                            :: fixed_sealevel
    REAL(dp)                            :: initial_guess_sealevel
    CHARACTER(LEN=256)                  :: filename_sealevel_record
    INTEGER                             :: sealevel_record_length

    CHARACTER(LEN=256)                  :: choice_GIA_model
    REAL(dp)                            :: ELRA_lithosphere_flex_rigidity
    REAL(dp)                            :: ELRA_bedrock_relaxation_time
    REAL(dp)                            :: ELRA_mantle_density

    ! SELEN
    ! =====

    LOGICAL                             :: SELEN_run_at_t_start
    INTEGER                             :: SELEN_n_TDOF_iterations
    INTEGER                             :: SELEN_n_recursion_iterations
    LOGICAL                             :: SELEN_use_rotational_feedback
    INTEGER                             :: SELEN_n_harmonics
    LOGICAL                             :: SELEN_display_progress

    CHARACTER(LEN=256)                  :: SELEN_dir
    CHARACTER(LEN=256)                  :: SELEN_global_topo_filename
    CHARACTER(LEN=256)                  :: SELEN_TABOO_init_filename
    CHARACTER(LEN=256)                  :: SELEN_LMJ_VALUES_filename

    INTEGER                             :: SELEN_irreg_time_n
    REAL(dp), DIMENSION(:), ALLOCATABLE :: SELEN_irreg_time_window

    REAL(dp)                            :: SELEN_lith_thickness
    INTEGER                             :: SELEN_visc_n
    REAL(dp), DIMENSION(:), ALLOCATABLE :: SELEN_visc_prof

    INTEGER                             :: SELEN_TABOO_CDE
    INTEGER                             :: SELEN_TABOO_TLOVE
    INTEGER                             :: SELEN_TABOO_DEG1
    REAL(dp)                            :: SELEN_TABOO_RCMB

    ! Some derived values
    INTEGER                             :: SELEN_i1, SELEN_i2           ! Parallelisation of loops over global grid pixels
    INTEGER                             :: SELEN_j1, SELEN_j2           ! Parallelisation of loops over harmonic degrees
    REAL(dp)                            :: SELEN_alfa
    INTEGER                             :: SELEN_jmax
    INTEGER                             :: SELEN_reg_time_n

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

    ! The output directory
    ! ====================

    CHARACTER(LEN=256)                  :: output_dir

  END TYPE constants_type

  ! ===============================================
  ! "C" is an instance of the "constants_type" type
  ! ===============================================

  TYPE(constants_type), SAVE :: C

! ===== The TABOO type =====
! ==========================
  ! Since some of the TABOO routines have variables named C (thanks, Giorgio...),
  ! we cannot use the regular config structure there. Collect the required config
  ! parameters into a smaller separate structure called C_TABOO

  TYPE constants_type_TABOO

    INTEGER                             :: IMODE            ! SELEN integration mode
    INTEGER                             :: NV               ! Number of viscoelastic layers
    REAL(dp), DIMENSION(:), ALLOCATABLE :: VSC              ! Viscosity profile
    INTEGER                             :: CDE              ! Code of the model (see taboo for explanation)
    INTEGER                             :: TLOVE            ! Tidal love numbers yes/no
    INTEGER                             :: DEG1             ! Tidal love numbers degree
    REAL(dp)                            :: LTH              ! Lithospheric thickness [km]
    REAL(dp)                            :: RCMB             ! Radius of CMB (km)

  END TYPE constants_type_TABOO

  ! ===========================================================
  ! "C_TABOO" is an instance of the "constants_type_TABOO" type
  ! ===========================================================

  TYPE(constants_type_TABOO), SAVE :: C_TABOO

CONTAINS

! ===== Model configuration =====
! ===============================

  SUBROUTINE initialise_model_configuration( version_number)
    ! Initialise the C (configuration) structure from one or two external config text files,
    ! set up the output directory (either procedurally from the current date, or directly
    ! from the config-specified folder name), and copy the config file(s) there.

    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: version_number

    ! Local variables:
    INTEGER                                            :: ierr, cerr, process_rank, number_of_processes, p
    LOGICAL                                            :: master
    CHARACTER(LEN=256)                                 :: config_filename, template_filename, variation_filename, config_mode
    INTEGER                                            :: i,n
    CHARACTER(LEN=20)                                  :: output_dir_procedural
    LOGICAL                                            :: ex

    ! Get rank of current process and total number of processes
    ! (needed because the configuration_module cannot access the par structure)
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, process_rank, ierr)
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, number_of_processes, ierr)
    master = (process_rank == 0)

    ! ===== Set up the config structure =====
    ! =======================================

    ! The name(s) of the config file(s) are provided as input arguments when calling the UFEMISM_program
    ! executable. After calling MPI_INIT, only the master process "sees" these arguments, so they need to be
    ! broadcast to the other processes.

    IF (master) THEN

      config_filename       = ''
      template_filename     = ''
      variation_filename    = ''
      config_mode           = ''

      IF     (iargc() == 0) THEN

        WRITE(0,*) ' ERROR: IMAU-ICE v', TRIM(version_number), ' needs at least one config file to run!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

      ELSEIF (iargc() == 1) THEN

        ! Run the model with a single config file
        CALL getarg( 1, config_filename)
        config_mode = 'single_config'

        WRITE(0,*) ''
        WRITE(0,*) ' Simulation settings from configuration file: ', TRIM(config_filename)
        WRITE(0,*) ''

      ELSEIF (iargc() == 2) THEN

        ! Run the model with two config files (template+variation)
        CALL getarg( 1, template_filename )
        CALL getarg( 2, variation_filename)
        config_mode = 'template+variation'

        WRITE(0,*) ''
        WRITE(0,*) ' Simulation settings from configuration files: ', TRIM(template_filename), ' & ', TRIM(variation_filename)
        WRITE(0,*) ''

      ELSE

        WRITE(0,*) ' ERROR: IMAU-ICE v', TRIM(version_number), ' can take either one or two config files to run!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

      END IF

    END IF ! IF (master) THEN

    CALL MPI_BCAST( config_filename,    256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( template_filename,  256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( variation_filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( config_mode,        256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Let each of the processors read the config file in turns so there's no access conflicts
    IF (config_mode == 'single_config') THEN
      ! Read only a single config file

      DO p = 0, number_of_processes-1
        IF (p == process_rank) THEN

          ! Read the external file, use a Fortran NAMELIST to overwrite the default
          ! values of the XXX_config variables
          CALL read_config_file( config_filename)

          ! Copy values from the XXX_config variables to the C structure
          CALL copy_variables_to_struct

        END IF
        CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
      END DO

    ELSEIF (config_mode == 'template+variation') THEN
      ! Read two config file consecutively: one "template" and one "variation"

      DO p = 0, number_of_processes-1
        IF (p == process_rank) THEN

          ! Read the external file, use a Fortran NAMELIST to overwrite the default
          ! values of the XXX_config variables

          ! First the template, then the variation
          CALL read_config_file( template_filename)
          CALL read_config_file( variation_filename)

          ! Copy values from the XXX_config variables to the C structure
          CALL copy_variables_to_struct

        END IF
        CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
      END DO

    ELSE ! IF (config_mode == 'single_config') THEN

      IF (master) WRITE(0,*) 'initialise_model_configuration - ERROR: unknown config_mode "', TRIM(config_mode), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    END IF ! IF (config_mode == 'single_config') THEN

    ! ===== Set up the output directory =====
    ! =======================================

    ! First get the name of the output directory (either procedural, or provided in the config file)

    DO n = 1, 256
      C%output_dir(n:n) = ' '
    END DO

    IF (C%create_procedural_output_dir) THEN
      ! Automatically create an output directory with a procedural name (e.g. results_20210720_001/)

      IF (master) THEN
        CALL get_procedural_output_dir_name( output_dir_procedural)
        C%output_dir(1:21) = TRIM(output_dir_procedural) // '/'
      END IF
      CALL MPI_BCAST( C%output_dir, 256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ELSE
      ! Use the provided name (return an error if this directory already exists)

      C%output_dir = TRIM(C%fixed_output_dir) // TRIM(C%fixed_output_dir_suffix) // '/'

      INQUIRE( FILE = TRIM(C%output_dir)//'/.', EXIST=ex)
      IF (ex) THEN
        WRITE(0,*) ' ERROR: fixed_output_dir_config ', TRIM(C%output_dir), ' already exists!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

    END IF

    ! Create the directory
    IF (master) THEN
      CALL system('mkdir ' // TRIM(C%output_dir))
      WRITE(0,*) ''
      WRITE(0,*) ' Output directory: ', TRIM(C%output_dir)
      WRITE(0,*) ''
    END IF
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)

    ! Copy the config file to the output directory
    IF (master) THEN
      IF     (config_mode == 'single_config') THEN
        CALL system('cp ' // config_filename    // ' ' // TRIM(C%output_dir))
      ELSEIF (config_mode == 'template+variation') THEN
        CALL system('cp ' // template_filename  // ' ' // TRIM(C%output_dir))
        CALL system('cp ' // variation_filename // ' ' // TRIM(C%output_dir))
      ELSE
        WRITE(0,*) ' initialise_model_configuration - ERROR: unknown config_mode "', TRIM(config_mode), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF ! IF (config_mode == 'single_config') THEN
    END IF ! IF (master) THEN
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)

    ! Set up the resource tracker
    ! ===========================

    ! Allocate space to track up to 2,000 subroutines per region
    ! That should be enough for a while...
    n = 0
    IF (C%do_NAM) n = n + 2000
    IF (C%do_EAS) n = n + 2000
    IF (C%do_GRL) n = n + 2000
    IF (C%do_ANT) n = n + 2000
    ALLOCATE( resource_tracker( n))

    ! Initialise values
    DO i = 1, n
      resource_tracker( i)%routine_path = 'subroutine_placeholder'
      resource_tracker( i)%tstart       = 0._dp
      resource_tracker( i)%tcomp        = 0._dp
      resource_tracker( i)%mem_use      = 0._dp
      resource_tracker( i)%mem_use_max  = 0._dp
    END DO

    ! Initialise the total model memory use tracker
    mem_use_tot     = 0._dp
    mem_use_tot_max = 0._dp

  END SUBROUTINE initialise_model_configuration

  SUBROUTINE read_config_file( config_filename)
    ! Use a NAMELIST containing all the "_config" variables to read
    ! an external config file, and overwrite the default values of
    ! the specified variables with the values from the file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),INTENT(IN) :: config_filename

    ! Local variables:
    CHARACTER(LEN=256)            :: namelist_filename
    INTEGER, PARAMETER            :: config_unit   = 1337
    INTEGER, PARAMETER            :: namelist_unit = 1338
    INTEGER                       :: ios, ierr, cerr

    ! The NAMELIST that's used to read the external config file.
    NAMELIST /CONFIG/start_time_of_run_config,                        &
                     end_time_of_run_config,                          &
                     dt_coupling_config,                              &
                     dt_max_config,                                   &
                     dt_min_config,                                   &
                     dt_startup_phase_config,                         &
                     dt_cooldown_phase_config,                        &
                     dt_thermo_config,                                &
                     dt_climate_config,                               &
                     dt_ocean_config,                                 &
                     dt_SMB_config,                                   &
                     dt_BMB_config,                                   &
                     dt_output_config,                                &
                     dt_mesh_min_config,                              &
                     dt_bedrock_ELRA_config,                          &
                     dt_SELEN_config,                                 &
                     dt_slid_inv_config,                              &
                     dt_SMB_inv_config,                               &
                     do_NAM_config,                                   &
                     do_EAS_config,                                   &
                     do_GRL_config,                                   &
                     do_ANT_config,                                   &
                     SSA_icestream_m_config,                          &
                     ISMIP_HOM_L_config,                              &
                     ISMIP_HOM_E_Arolla_filename_config,              &
                     MISMIPplus_do_tune_A_for_GL_config,              &
                     MISMIPplus_xGL_target_config,                    &
                     MISMIPplus_A_flow_initial_config,                &
                     MISMIPplus_scenario_config,                      &
                     MISOMIP1_scenario_config,                        &
                     create_procedural_output_dir_config,             &
                     fixed_output_dir_config,                         &
                     fixed_output_dir_suffix_config,                  &
                     do_write_regional_scalar_output_config,          &
                     do_write_global_scalar_output_config,            &
                     do_write_debug_data_config,                      &
                     do_write_grid_data_config,                       &
                     do_check_for_NaN_config,                         &
                     do_time_display_config,                          &
                     lambda_M_NAM_config,                             &
                     phi_M_NAM_config,                                &
                     beta_stereo_NAM_config,                          &
                     xmin_NAM_config,                                 &
                     xmax_NAM_config,                                 &
                     ymin_NAM_config,                                 &
                     ymax_NAM_config,                                 &
                     lambda_M_EAS_config,                             &
                     phi_M_EAS_config,                                &
                     beta_stereo_EAS_config,                          &
                     xmin_EAS_config,                                 &
                     xmax_EAS_config,                                 &
                     ymin_EAS_config,                                 &
                     ymax_EAS_config,                                 &
                     lambda_M_GRL_config,                             &
                     phi_M_GRL_config,                                &
                     beta_stereo_GRL_config,                          &
                     xmin_GRL_config,                                 &
                     xmax_GRL_config,                                 &
                     ymin_GRL_config,                                 &
                     ymax_GRL_config,                                 &
                     lambda_M_ANT_config,                             &
                     phi_M_ANT_config,                                &
                     beta_stereo_ANT_config,                          &
                     xmin_ANT_config,                                 &
                     xmax_ANT_config,                                 &
                     ymin_ANT_config,                                 &
                     ymax_ANT_config,                                 &
                     nz_config,                                       &
                     zeta_config,                                     &
                     refgeo_Hi_min_config,                            &
                     remove_Lake_Vostok_config,                       &
                     choice_refgeo_init_NAM_config,                   &
                     choice_refgeo_init_EAS_config,                   &
                     choice_refgeo_init_GRL_config,                   &
                     choice_refgeo_init_ANT_config,                   &
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
                     is_restart_config,                               &
                     time_to_restart_from_NAM_config,                 &
                     time_to_restart_from_EAS_config,                 &
                     time_to_restart_from_GRL_config,                 &
                     time_to_restart_from_ANT_config,                 &
                     filename_restart_NAM_config,                     &
                     filename_restart_EAS_config,                     &
                     filename_restart_GRL_config,                     &
                     filename_restart_ANT_config,                     &
                     choice_forcing_method_config,                    &
                     choice_insolation_forcing_config,                &
                     static_insolation_time_config,                   &
                     filename_insolation_config,                      &
                     filename_CO2_record_config,                      &
                     CO2_record_length_config,                        &
                     filename_d18O_record_config,                     &
                     d18O_record_length_config,                       &
                     choice_geothermal_heat_flux_config,              &
                     constant_geothermal_heat_flux_config,            &
                     filename_geothermal_heat_flux_config,            &
                     do_calculate_benthic_d18O_config,                &
                     dT_deepwater_averaging_window_config,            &
                     dT_deepwater_dT_surf_ratio_config,               &
                     d18O_dT_deepwater_ratio_config,                  &
                     dT_glob_inverse_averaging_window_config,         &
                     inverse_d18O_to_dT_glob_scaling_config,          &
                     CO2_inverse_averaging_window_config,             &
                     inverse_d18O_to_CO2_scaling_config,              &
                     inverse_d18O_to_CO2_initial_CO2_config,          &
                     choice_ice_dynamics_config,                      &
                     n_flow_config,                                   &
                     m_enh_sheet_config,                              &
                     m_enh_shelf_config,                              &
                     choice_ice_margin_config,                        &
                     do_hybrid_Bernales2017_config,                   &
                     vel_ref_Bernales2017_config,                     &
                     include_SSADIVA_crossterms_config,               &
                     do_GL_subgrid_friction_config,                   &
                     do_smooth_geometry_config,                       &
                     r_smooth_geometry_config,                        &
                     DIVA_visc_it_norm_dUV_tol_config,                &
                     DIVA_visc_it_nit_config,                         &
                     DIVA_visc_it_relax_config,                       &
                     DIVA_epsilon_sq_0_config,                        &
                     DIVA_visc_eff_min_config,                        &
                     DIVA_beta_max_config,                            &
                     DIVA_vel_max_config,                             &
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
                     windup_total_years_config,                       &
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
                     ice_thickness_west_BC_config,                    &
                     ice_thickness_east_BC_config,                    &
                     ice_thickness_south_BC_config,                   &
                     ice_thickness_north_BC_config,                   &
                     choice_mask_noice_NAM_config,                    &
                     choice_mask_noice_EAS_config,                    &
                     choice_mask_noice_GRL_config,                    &
                     choice_mask_noice_ANT_config,                    &
                     fixed_sheet_geometry_config,                     &
                     fixed_shelf_geometry_config,                     &
                     fixed_grounding_line_g_config,                   &
                     fixed_grounding_line_f_config,                   &
                     dHi_dt_window_size_config,                       &
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
                     basal_roughness_restart_type_config,             &
                     do_basal_roughness_remap_adjustment_config,      &
                     do_slid_inv_config,                              &
                     choice_slid_inv_method_config,                   &
                     slid_inv_t_start_config,                         &
                     slid_inv_t_end_config,                           &
                     slid_inv_phi_min_config,                         &
                     slid_inv_phi_max_config,                         &
                     slid_inv_filename_output_config,                 &
                     slid_inv_window_size_config,                     &
                     do_slid_inv_Bernales2017_smooth_config,          &
                     do_slid_inv_Bernales2017_extrap_config,          &
                     slid_inv_Bernales2017_hi_scale_config,           &
                     slid_inv_Bernales2017_smooth_r_config,           &
                     slid_inv_Bernales2017_smooth_w_config,           &
                     slid_inv_Bernales2017_tol_diff_config,           &
                     slid_inv_Bernales2017_tol_frac_config,           &
                     slid_inv_Berends2022_tauc_config,                &
                     slid_inv_Berends2022_H0_config,                  &
                     slid_inv_Berends2022_u0_config,                  &
                     slid_inv_Berends2022_Hi_scale_config,            &
                     slid_inv_Berends2022_u_scale_config,             &
                     slid_inv_target_velocity_filename_config,        &
                     choice_calving_law_config,                       &
                     calving_threshold_thickness_shelf_config,        &
                     calving_threshold_thickness_sheet_config,        &
                     max_calving_rounds_config,                       &
                     do_remove_shelves_config,                        &
                     remove_shelves_larger_than_PD_config,            &
                     continental_shelf_calving_config,                &
                     continental_shelf_min_height_config,             &
                     minimum_ice_thickness_config,                    &
                     use_submesh_config,                              &
                     nconmax_config,                                  &
                     alpha_min_config,                                &
                     dz_max_ice_config,                               &
                     res_max_config,                                  &
                     res_max_ice_config,                              &
                     res_max_margin_config,                           &
                     res_max_gl_config,                               &
                     res_max_cf_config,                               &
                     res_max_mountain_config,                         &
                     res_max_coast_config,                            &
                     mesh_fitness_threshold_config,                   &
                     do_force_mesh_update_config,                     &
                     do_force_mesh_update_after_config,               &
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
                     do_write_ISMIP_output_config,                    &
                     ISMIP_output_group_code_config,                  &
                     ISMIP_output_model_code_config,                  &
                     ISMIP_output_experiment_code_config,             &
                     ISMIP_output_basetime_config,                    &
                     choice_initial_ice_temperature_config,           &
                     uniform_ice_temperature_config,                  &
                     choice_thermo_model_config,                      &
                     choice_ice_rheology_config,                      &
                     uniform_flow_factor_config,                      &
                     choice_ice_heat_capacity_config,                 &
                     uniform_ice_heat_capacity_config,                &
                     choice_ice_thermal_conductivity_config,          &
                     uniform_ice_thermal_conductivity_config,         &
                     choice_climate_model_config,                     &
                     choice_idealised_climate_config,                 &
                     filename_direct_global_climate_config,           &
                     filename_direct_regional_climate_NAM_config,     &
                     filename_direct_regional_climate_EAS_config,     &
                     filename_direct_regional_climate_GRL_config,     &
                     filename_direct_regional_climate_ANT_config,     &
                     filename_PD_obs_climate_config,                  &
                     filename_climate_snapshot_PI_config,             &
                     filename_climate_snapshot_warm_config,           &
                     filename_climate_snapshot_cold_config,           &
                     constant_lapserate_config,                       &
                     matrix_high_CO2_level_config,                    &
                     matrix_low_CO2_level_config,                     &
                     matrix_warm_orbit_time_config,                   &
                     matrix_cold_orbit_time_config,                   &
                     climate_matrix_biascorrect_warm_config,          &
                     climate_matrix_biascorrect_cold_config,          &
                     switch_glacial_index_config,                     &
                     do_clim_inv_config,                              &
                     clim_inv_t_start_config,                         &
                     clim_inv_t_end_config,                           &
                     clim_inv_Hb_min_config,                          &
                     clim_inv_Hi_max_config,                          &
                     clim_inv_window_size_config,                     &
                     choice_ocean_model_config,                       &
                     choice_idealised_ocean_config,                   &
                     filename_PD_obs_ocean_config,                    &
                     name_ocean_temperature_obs_config,               &
                     name_ocean_salinity_obs_config,                  &
                     use_inverted_ocean_config,                       &
                     filename_GCM_ocean_snapshot_PI_config,           &
                     filename_GCM_ocean_snapshot_warm_config,         &
                     filename_GCM_ocean_snapshot_cold_config,         &
                     name_ocean_temperature_GCM_config,               &
                     name_ocean_salinity_GCM_config,                  &
                     ocean_temperature_PD_config,                     &
                     ocean_temperature_cold_config,                   &
                     ocean_temperature_warm_config,                   &
                     choice_ocean_vertical_grid_config,               &
                     ocean_vertical_grid_max_depth_config,            &
                     ocean_regular_grid_dz_config,                    &
                     ocean_extrap_dir_config,                         &
                     ocean_extrap_res_config,                         &
                     ocean_extrap_Gauss_sigma_config,                 &
                     ocean_extrap_hires_geo_filename_NAM_config,      &
                     ocean_extrap_hires_geo_filename_EAS_config,      &
                     ocean_extrap_hires_geo_filename_GRL_config,      &
                     ocean_extrap_hires_geo_filename_ANT_config,      &
                     ocean_w_tot_hist_averaging_window_config,        &
                     choice_SMB_model_config,                         &
                     choice_idealised_SMB_config,                     &
                     SMB_uniform_config,                              &
                     filename_direct_global_SMB_config,               &
                     filename_direct_regional_SMB_NAM_config,         &
                     filename_direct_regional_SMB_EAS_config,         &
                     filename_direct_regional_SMB_GRL_config,         &
                     filename_direct_regional_SMB_ANT_config,         &
                     SMB_IMAUITM_choice_init_firn_NAM_config,         &
                     SMB_IMAUITM_choice_init_firn_EAS_config,         &
                     SMB_IMAUITM_choice_init_firn_GRL_config,         &
                     SMB_IMAUITM_choice_init_firn_ANT_config,         &
                     SMB_IMAUITM_initial_firn_thickness_config,       &
                     SMB_IMAUITM_C_abl_constant_NAM_config,           &
                     SMB_IMAUITM_C_abl_constant_EAS_config,           &
                     SMB_IMAUITM_C_abl_constant_GRL_config,           &
                     SMB_IMAUITM_C_abl_constant_ANT_config,           &
                     SMB_IMAUITM_C_abl_Ts_NAM_config,                 &
                     SMB_IMAUITM_C_abl_Ts_EAS_config,                 &
                     SMB_IMAUITM_C_abl_Ts_GRL_config,                 &
                     SMB_IMAUITM_C_abl_Ts_ANT_config,                 &
                     SMB_IMAUITM_C_abl_Q_NAM_config,                  &
                     SMB_IMAUITM_C_abl_Q_EAS_config,                  &
                     SMB_IMAUITM_C_abl_Q_GRL_config,                  &
                     SMB_IMAUITM_C_abl_Q_ANT_config,                  &
                     SMB_IMAUITM_C_refr_NAM_config,                   &
                     SMB_IMAUITM_C_refr_EAS_config,                   &
                     SMB_IMAUITM_C_refr_GRL_config,                   &
                     SMB_IMAUITM_C_refr_ANT_config,                   &
                     do_SMB_IMAUITM_inversion_config,                 &
                     SMB_IMAUITM_inv_choice_init_C_config,            &
                     SMB_IMAUITM_inv_scale_config,                    &
                     SMB_IMAUITM_inv_C_abl_constant_min_config,       &
                     SMB_IMAUITM_inv_C_abl_constant_max_config,       &
                     SMB_IMAUITM_inv_C_abl_Ts_min_config,             &
                     SMB_IMAUITM_inv_C_abl_Ts_max_config,             &
                     SMB_IMAUITM_inv_C_abl_Q_min_config,              &
                     SMB_IMAUITM_inv_C_abl_Q_max_config,              &
                     SMB_IMAUITM_inv_C_refr_min_config,               &
                     SMB_IMAUITM_inv_C_refr_max_config,               &
                     ISMIP_forcing_filename_baseline_config,          &
                     ISMIP_forcing_foldername_aSMB_config,            &
                     ISMIP_forcing_basefilename_aSMB_config,          &
                     ISMIP_forcing_foldername_dSMBdz_config,          &
                     ISMIP_forcing_basefilename_dSMBdz_config,        &
                     ISMIP_forcing_foldername_aST_config,             &
                     ISMIP_forcing_basefilename_aST_config,           &
                     ISMIP_forcing_foldername_dSTdz_config,           &
                     ISMIP_forcing_basefilename_dSTdz_config,         &
                     choice_BMB_shelf_model_config,                   &
                     choice_idealised_BMB_shelf_config,               &
                     choice_BMB_sheet_model_config,                   &
                     BMB_shelf_uniform_config,                        &
                     BMB_sheet_uniform_config,                        &
                     choice_BMB_subgrid_config,                       &
                     BMB_max_config,                                  &
                     BMB_min_config,                                  &
                     do_ocean_temperature_inversion_config,           &
                     ocean_temperature_inv_t_start_config,            &
                     ocean_temperature_inv_t_end_config,              &
                     T_base_window_size_config,                       &
                     BMB_inv_use_restart_field_config,                &
                     BMB_inv_scale_shelf_config,                      &
                     BMB_inv_scale_ocean_config,                      &
                     choice_basin_scheme_NAM_config,                  &
                     choice_basin_scheme_EAS_config,                  &
                     choice_basin_scheme_GRL_config,                  &
                     choice_basin_scheme_ANT_config,                  &
                     filename_basins_NAM_config,                      &
                     filename_basins_EAS_config,                      &
                     filename_basins_GRL_config,                      &
                     filename_basins_ANT_config,                      &
                     choice_BMB_shelf_amplification_config,           &
                     basin_BMB_amplification_n_ANT_config,            &
                     basin_BMB_amplification_factor_ANT_config,       &
                     basin_BMB_amplification_n_GRL_config,            &
                     basin_BMB_amplification_factor_GRL_config,       &
                     do_merge_basins_ANT_config,                      &
                     do_merge_basins_GRL_config,                      &
                     BMB_Favier2019_lin_GammaT_config,                &
                     BMB_Favier2019_quad_GammaT_config,               &
                     BMB_Favier2019_Mplus_GammaT_config,              &
                     BMB_Lazeroms2018_GammaT_config,                  &
                     BMB_Lazeroms2018_find_GL_scheme_config,          &
                     BMB_PICO_nboxes_config,                          &
                     BMB_PICO_GammaTstar_config,                      &
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
                     choice_ice_isotopes_model_config,                &
                     uniform_ice_d18O_config,                         &
                     do_ocean_floodfill_config,                       &
                     choice_sealevel_model_config,                    &
                     fixed_sealevel_config,                           &
                     initial_guess_sealevel_config,                   &
                     filename_sealevel_record_config,                 &
                     sealevel_record_length_config,                   &
                     choice_GIA_model_config,                         &
                     ELRA_lithosphere_flex_rigidity_config,           &
                     ELRA_bedrock_relaxation_time_config,             &
                     ELRA_mantle_density_config,                      &
                     SELEN_run_at_t_start_config,                     &
                     SELEN_n_TDOF_iterations_config,                  &
                     SELEN_n_recursion_iterations_config,             &
                     SELEN_use_rotational_feedback_config,            &
                     SELEN_n_harmonics_config,                        &
                     SELEN_display_progress_config,                   &
                     SELEN_dir_config,                                &
                     SELEN_global_topo_filename_config,               &
                     SELEN_TABOO_init_filename_config,                &
                     SELEN_LMJ_VALUES_filename_config,                &
                     SELEN_irreg_time_n_config,                       &
                     SELEN_irreg_time_window_config,                  &
                     SELEN_lith_thickness_config,                     &
                     SELEN_visc_n_config,                             &
                     SELEN_visc_prof_config,                          &
                     SELEN_TABOO_CDE_config,                          &
                     SELEN_TABOO_TLOVE_config,                        &
                     SELEN_TABOO_DEG1_config,                         &
                     SELEN_TABOO_RCMB_config,                         &
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

    ! Write the CONFIG namelist to a temporary file
    namelist_filename = 'config_namelist_temp.txt'
    OPEN(  UNIT = namelist_unit, FILE = TRIM( namelist_filename))
    WRITE( UNIT = namelist_unit, NML  = CONFIG)
    CLOSE( UNIT = namelist_unit)

    ! Check the config file for validity
    CALL check_config_file_validity( config_filename, namelist_filename)

    ! Delete the temporary CONFIG namelist file
    CALL system('rm -f ' // TRIM( namelist_filename))

    ! Open the config file
    OPEN(  UNIT = config_unit, FILE = TRIM( config_filename), STATUS = 'OLD', ACTION = 'READ', IOSTAT = ios)
    IF (ios /= 0) THEN
      WRITE(0,'(A,A,A)') colour_string('ERROR: config file "' // TRIM( config_filename),'red') // '" not found!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Read the config file using the CONFIG namelist
    READ(  UNIT = config_unit, NML = CONFIG, IOSTAT = ios)
    IF (ios /= 0) THEN
      WRITE(0,'(A,A,A)') colour_string('ERROR: error while reading config file "' // TRIM( config_filename),'red') // '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Close the config file
    CLOSE( UNIT = config_unit)

  END SUBROUTINE read_config_file

  SUBROUTINE check_config_file_validity( config_filename, namelist_filename)
    ! Check if the provided config file is valid
    !
    ! Do this by reading one line at a time of the config file, determining the name of the variable
    ! declared in that line, and checking if that variable also exists in the namelist file
    !
    ! Assumes that the CONFIG namelist has already been written to the specified file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),INTENT(IN) :: config_filename, namelist_filename

    ! Local variables:
    INTEGER, PARAMETER            :: config_unit   = 1337
    INTEGER, PARAMETER            :: namelist_unit = 1338
    INTEGER                       :: ios, ierr, cerr
    LOGICAL                       :: found_end_of_file_config, found_end_of_file_namelist
    CHARACTER(256)                :: single_line_config      , single_line_namelist
    INTEGER                       :: line_counter_config     , line_counter_namelist
    LOGICAL                       :: found_match, found_mismatch

    ! Open the config and namelist files
    OPEN( UNIT = config_unit, FILE = config_filename, IOSTAT = ios)
    IF (ios /= 0) THEN
      WRITE(0,'(A)') colour_string('ERROR','red') // ': config file "' // TRIM( config_filename) // '" not found!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Read one line at a time of the config file, determine the name of the variable
    ! declared in that line, and check if that variable also exists in the namelist file

    found_end_of_file_config = .FALSE.
    line_counter_config      = 0
    found_mismatch           = .FALSE.

    DO WHILE (.NOT. found_end_of_file_config)

      line_counter_config = line_counter_config + 1

      ! Read a single line from the config file
      READ( UNIT = config_unit, FMT = '(A)', IOSTAT = ios) single_line_config

      ! If we've reached the end of the file before finding the terminating forward slash, this config file is not valid.
      IF (ios < 0) THEN
        WRITE(0,'(A)') colour_string('ERROR','red') // ': config file "' // TRIM( config_filename) // '" is not terminated with a forward slash!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      ! Remove all leading spaces
      CALL remove_leading_spaces( single_line_config)

      ! The variable name is the part of the string left of the first (, =, or space.
      single_line_config = single_line_config( 1:SCAN( single_line_config, '( =')-1)

      ! Get config variable in all caps for case-insensitive comparison
      CALL capitalise_string( single_line_config)

      ! The forward slash at the end terminates the config file
      IF (single_line_config == '/') THEN
        found_end_of_file_config = .TRUE.
      END IF

      ! Disregard empty lines, commented lines, and the header line
      IF (single_line_config == '' .OR. single_line_config == '&CONFIG' .OR. single_line_config( 1:1) == '!') THEN
        CYCLE
      END IF

      ! Open the namelist file
      OPEN( UNIT = namelist_unit, FILE = namelist_filename)
      IF (ios /= 0) THEN
        WRITE(0,'(A)') colour_string('ERROR','red') // ': namelist file "' // TRIM( namelist_filename) // '" not found!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      ! Read all variables from the namelist file and check if any of them match the current config variable

      found_end_of_file_namelist = .FALSE.
      line_counter_namelist      = 0
      found_match                = .FALSE.

      DO WHILE ((.NOT. found_end_of_file_namelist) .AND. (.NOT. found_match))

        line_counter_namelist = line_counter_namelist + 1

        ! Read a single line from the namelist file
        READ( UNIT = namelist_unit, FMT = '(A)', IOSTAT = ios) single_line_namelist

        ! If we've reached the end of the file before finding the terminating forward slash, this namelist file is not valid.
        IF (ios < 0) THEN
          WRITE(0,'(A)') colour_string('ERROR','red') // ': namelist file "' // TRIM( namelist_filename) // '" is not terminated with a forward slash!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF

        ! Remove all leading spaces
        CALL remove_leading_spaces( single_line_namelist)

        ! The variable name is the part of the string left of the first (, =, or space.
        single_line_namelist = single_line_namelist( 1:SCAN( single_line_namelist, '( =')-1)

        ! Get namelist variable in all caps for case-insensitive comparison
        CALL capitalise_string( single_line_namelist)

        ! The forward slash at the end terminates the config file
        IF (single_line_namelist == '/') THEN
          found_end_of_file_namelist = .TRUE.
        END IF

        ! Disregard empty lines, commented lines, and the header line
        IF (single_line_namelist == '' .OR. single_line_namelist == '&CONFIG' .OR. single_line_namelist( 1:1) == '!') THEN
          CYCLE
        END IF

        ! Check if this namelist variable matches the config variable
        IF (single_line_namelist == single_line_config) THEN
          found_match = .TRUE.
        END IF

      END DO ! DO WHILE ((.NOT. found_end_of_file_namelist) .AND. (.NOT. found_match))

      ! If no matching variable was found in the namelist file, print an error
      IF (.NOT. found_match) THEN
        WRITE(0,'(A,I4)') colour_string('ERROR','red') // ': invalid config variable "' // TRIM( single_line_config) // &
          '" in file "' // TRIM( config_filename) // '", line ', line_counter_config
        found_mismatch = .TRUE.
      END IF

      ! Close the namelist file
      CLOSE( UNIT = namelist_unit)

    END DO ! DO WHILE (.NOT. found_end_of_file_config)

    ! Close the config file
    CLOSE( UNIT = config_unit)

    ! If an invalid config variable was found, crash.
    IF (found_mismatch) CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

  END SUBROUTINE check_config_file_validity

  SUBROUTINE copy_variables_to_struct
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
    C%dt_min                                   = dt_min_config
    C%dt_startup_phase                         = dt_startup_phase_config
    C%dt_cooldown_phase                        = dt_cooldown_phase_config
    C%dt_thermo                                = dt_thermo_config
    C%dt_climate                               = dt_climate_config
    C%dt_ocean                                 = dt_ocean_config
    C%dt_SMB                                   = dt_SMB_config
    C%dt_BMB                                   = dt_BMB_config
    C%dt_output                                = dt_output_config
    C%dt_mesh_min                              = dt_mesh_min_config
    C%dt_bedrock_ELRA                          = dt_bedrock_ELRA_config
    C%dt_SELEN                                 = dt_SELEN_config
    C%dt_slid_inv                              = dt_slid_inv_config
    C%dt_SMB_inv                               = dt_SMB_inv_config

    ! Which ice sheets do we simulate?
    ! ================================

    C%do_NAM                                   = do_NAM_config
    C%do_EAS                                   = do_EAS_config
    C%do_GRL                                   = do_GRL_config
    C%do_ANT                                   = do_ANT_config

    ! Benchmark experiments
    ! =====================

    C%SSA_icestream_m                          = SSA_icestream_m_config
    C%ISMIP_HOM_L                              = ISMIP_HOM_L_config
    C%ISMIP_HOM_E_Arolla_filename              = ISMIP_HOM_E_Arolla_filename_config
    C%MISMIPplus_do_tune_A_for_GL              = MISMIPplus_do_tune_A_for_GL_config
    C%MISMIPplus_xGL_target                    = MISMIPplus_xGL_target_config
    C%MISMIPplus_A_flow_initial                = MISMIPplus_A_flow_initial_config
    C%MISMIPplus_scenario                      = MISMIPplus_scenario_config
    C%MISOMIP1_scenario                        = MISOMIP1_scenario_config

    ! Whether or not to let UFEMISM dynamically create its own output folder
    ! =======================================================================

    C%create_procedural_output_dir             = create_procedural_output_dir_config
    C%fixed_output_dir                         = fixed_output_dir_config
    C%fixed_output_dir_suffix                  = fixed_output_dir_suffix_config
    C%do_write_regional_scalar_output          = do_write_regional_scalar_output_config
    C%do_write_global_scalar_output            = do_write_global_scalar_output_config

    ! Debugging
    ! =========

    C%do_write_debug_data                      = do_write_debug_data_config
    C%do_write_grid_data                       = do_write_grid_data_config
    C%do_check_for_NaN                         = do_check_for_NaN_config
    C%do_time_display                          = do_time_display_config

    ! == The four model regions
    ! =========================

    ! North America
    C%lambda_M_NAM                             = lambda_M_NAM_config
    C%phi_M_NAM                                = phi_M_NAM_config
    C%beta_stereo_NAM                          = beta_stereo_NAM_config
    C%xmin_NAM                                 = xmin_NAM_config
    C%xmax_NAM                                 = xmax_NAM_config
    C%ymin_NAM                                 = ymin_NAM_config
    C%ymax_NAM                                 = ymax_NAM_config

    ! Eurasia
    C%lambda_M_EAS                             = lambda_M_EAS_config
    C%phi_M_EAS                                = phi_M_EAS_config
    C%beta_stereo_EAS                          = beta_stereo_EAS_config
    C%xmin_EAS                                 = xmin_EAS_config
    C%xmax_EAS                                 = xmax_EAS_config
    C%ymin_EAS                                 = ymin_EAS_config
    C%ymax_EAS                                 = ymax_EAS_config

    ! Greenland
    C%lambda_M_GRL                             = lambda_M_GRL_config
    C%phi_M_GRL                                = phi_M_GRL_config
    C%beta_stereo_GRL                          = beta_stereo_GRL_config
    C%xmin_GRL                                 = xmin_GRL_config
    C%xmax_GRL                                 = xmax_GRL_config
    C%ymin_GRL                                 = ymin_GRL_config
    C%ymax_GRL                                 = ymax_GRL_config

    ! Antarctica
    C%lambda_M_ANT                             = lambda_M_ANT_config
    C%phi_M_ANT                                = phi_M_ANT_config
    C%beta_stereo_ANT                          = beta_stereo_ANT_config
    C%xmin_ANT                                 = xmin_ANT_config
    C%xmax_ANT                                 = xmax_ANT_config
    C%ymin_ANT                                 = ymin_ANT_config
    C%ymax_ANT                                 = ymax_ANT_config

    ! Mesh generation parameters
    ! ==========================

    C%use_submesh                              = use_submesh_config
    C%nconmax                                  = nconmax_config
    C%alpha_min                                = alpha_min_config
    C%dz_max_ice                               = dz_max_ice_config
    C%res_max                                  = res_max_config
    C%res_max_ice                              = res_max_ice_config
    C%res_max_margin                           = res_max_margin_config
    C%res_max_gl                               = res_max_gl_config
    C%res_max_cf                               = res_max_cf_config
    C%res_max_mountain                         = res_max_mountain_config
    C%res_max_coast                            = res_max_coast_config
    C%mesh_fitness_threshold                   = mesh_fitness_threshold_config
    C%do_force_mesh_update                     = do_force_mesh_update_config
    C%do_force_mesh_update_after               = do_force_mesh_update_after_config

    ! The smallest allowed resolution
    C%res_min = MIN( MIN( MIN( MIN( C%res_max_margin, C%res_max_gl), C%res_max_cf), C%res_max_mountain), C%res_max_coast, C%res_max_ice)

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

    ! ISMIP-style output
    ! ==================

    C%do_write_ISMIP_output                    = do_write_ISMIP_output_config
    C%ISMIP_output_group_code                  = ISMIP_output_group_code_config
    C%ISMIP_output_model_code                  = ISMIP_output_model_code_config
    C%ISMIP_output_experiment_code             = ISMIP_output_experiment_code_config
    C%ISMIP_output_basetime                    = ISMIP_output_basetime_config

    ! Scaled vertical coordinate zeta
    ! ===============================

    C%nz     = nz_config
    ALLOCATE( C%zeta( C%nz))
    C%zeta   = zeta_config( 1:C%nz)

    ! Reference geometries (initial, present-day, and GIA equilibrium)
    ! ================================================================

    ! Some pre-processing stuff for reference ice geometry
    C%refgeo_Hi_min                            = refgeo_Hi_min_config
    C%remove_Lake_Vostok                       = remove_Lake_Vostok_config

    ! Initial geometry
    C%choice_refgeo_init_NAM                   = choice_refgeo_init_NAM_config
    C%choice_refgeo_init_EAS                   = choice_refgeo_init_EAS_config
    C%choice_refgeo_init_GRL                   = choice_refgeo_init_GRL_config
    C%choice_refgeo_init_ANT                   = choice_refgeo_init_ANT_config
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

    ! Whether or not the simulation is a restart of a previous simulation
    ! ===================================================================

    C%is_restart                               = is_restart_config

    C%time_to_restart_from_NAM                 = time_to_restart_from_NAM_config
    C%time_to_restart_from_EAS                 = time_to_restart_from_EAS_config
    C%time_to_restart_from_GRL                 = time_to_restart_from_GRL_config
    C%time_to_restart_from_ANT                 = time_to_restart_from_ANT_config

    C%filename_restart_NAM                     = filename_restart_NAM_config
    C%filename_restart_EAS                     = filename_restart_EAS_config
    C%filename_restart_GRL                     = filename_restart_GRL_config
    C%filename_restart_ANT                     = filename_restart_ANT_config

    ! Global forcing (insolation, CO2, d18O, geothermal heat flux)
    ! ============================================================

    C%choice_forcing_method                    = choice_forcing_method_config

    ! Insolation forcing (NetCDF)
    C%choice_insolation_forcing                = choice_insolation_forcing_config
    C%static_insolation_time                   = static_insolation_time_config
    C%filename_insolation                      = filename_insolation_config

    ! CO2 record (ASCII text file, so the number of rows needs to be specified)
    C%filename_CO2_record                      = filename_CO2_record_config
    C%CO2_record_length                        = CO2_record_length_config

    ! d18O record (ASCII text file, so the number of rows needs to be specified)
    C%filename_d18O_record                     = filename_d18O_record_config
    C%d18O_record_length                       = d18O_record_length_config

    ! Geothermal heat flux
    C%choice_geothermal_heat_flux              = choice_geothermal_heat_flux_config
    C%constant_geothermal_heat_flux            = constant_geothermal_heat_flux_config
    C%filename_geothermal_heat_flux            = filename_geothermal_heat_flux_config

    ! Parameters for calculating modelled benthic d18O
    C%do_calculate_benthic_d18O                = do_calculate_benthic_d18O_config
    C%dT_deepwater_averaging_window            = dT_deepwater_averaging_window_config
    C%dT_deepwater_dT_surf_ratio               = dT_deepwater_dT_surf_ratio_config
    C%d18O_dT_deepwater_ratio                  = d18O_dT_deepwater_ratio_config

    ! Parameters for the inverse routine
    C%dT_glob_inverse_averaging_window         = dT_glob_inverse_averaging_window_config
    C%inverse_d18O_to_dT_glob_scaling          = inverse_d18O_to_dT_glob_scaling_config
    C%CO2_inverse_averaging_window             = CO2_inverse_averaging_window_config
    C%inverse_d18O_to_CO2_scaling              = inverse_d18O_to_CO2_scaling_config
    C%inverse_d18O_to_CO2_initial_CO2          = inverse_d18O_to_CO2_initial_CO2_config

    ! Ice dynamics - velocity
    ! =======================

    C%choice_ice_dynamics                      = choice_ice_dynamics_config
    C%n_flow                                   = n_flow_config
    C%m_enh_sheet                              = m_enh_sheet_config
    C%m_enh_shelf                              = m_enh_shelf_config
    C%choice_ice_margin                        = choice_ice_margin_config
    C%do_hybrid_Bernales2017                   = do_hybrid_Bernales2017_config
    C%vel_ref_Bernales2017                     = vel_ref_Bernales2017_config
    C%include_SSADIVA_crossterms               = include_SSADIVA_crossterms_config
    C%do_GL_subgrid_friction                   = do_GL_subgrid_friction_config
    C%do_smooth_geometry                       = do_smooth_geometry_config
    C%r_smooth_geometry                        = r_smooth_geometry_config

    ! Some parameters for numerically solving the SSA/DIVA
    C%DIVA_visc_it_norm_dUV_tol                = DIVA_visc_it_norm_dUV_tol_config
    C%DIVA_visc_it_nit                         = DIVA_visc_it_nit_config
    C%DIVA_visc_it_relax                       = DIVA_visc_it_relax_config
    C%DIVA_epsilon_sq_0                        = DIVA_epsilon_sq_0_config
    C%DIVA_visc_eff_min                        = DIVA_visc_eff_min_config
    C%DIVA_beta_max                            = DIVA_beta_max_config
    C%DIVA_vel_max                             = DIVA_vel_max_config
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

    ! Velocity wind-up before a restart
    C%windup_total_years                       = windup_total_years_config

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

    ! Ice thickness boundary conditions
    C%ice_thickness_west_BC                    = ice_thickness_west_BC_config
    C%ice_thickness_east_BC                    = ice_thickness_east_BC_config
    C%ice_thickness_south_BC                   = ice_thickness_south_BC_config
    C%ice_thickness_north_BC                   = ice_thickness_north_BC_config
    C%choice_mask_noice_NAM                    = choice_mask_noice_NAM_config
    C%choice_mask_noice_EAS                    = choice_mask_noice_EAS_config
    C%choice_mask_noice_GRL                    = choice_mask_noice_GRL_config
    C%choice_mask_noice_ANT                    = choice_mask_noice_ANT_config

    ! Partially fixed geometry, useful for initialisation and inversion runs
    C%fixed_sheet_geometry                     = fixed_sheet_geometry_config
    C%fixed_shelf_geometry                     = fixed_shelf_geometry_config
    C%fixed_grounding_line_g                   = fixed_grounding_line_g_config
    C%fixed_grounding_line_f                   = fixed_grounding_line_f_config

    ! Memory of previous run during a restart
    C%dHi_dt_window_size                       = dHi_dt_window_size_config

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
    C%basal_roughness_restart_type             = basal_roughness_restart_type_config
    C%do_basal_roughness_remap_adjustment      = do_basal_roughness_remap_adjustment_config

    ! Bed roughness inversion
    C%do_slid_inv                              = do_slid_inv_config
    C%choice_slid_inv_method                   = choice_slid_inv_method_config
    C%slid_inv_t_start                         = slid_inv_t_start_config
    C%slid_inv_t_end                           = slid_inv_t_end_config
    C%slid_inv_phi_min                         = slid_inv_phi_min_config
    C%slid_inv_phi_max                         = slid_inv_phi_max_config
    C%slid_inv_filename_output                 = slid_inv_filename_output_config
    C%slid_inv_window_size                     = slid_inv_window_size_config
    C%do_slid_inv_Bernales2017_smooth          = do_slid_inv_Bernales2017_smooth_config
    C%do_slid_inv_Bernales2017_extrap          = do_slid_inv_Bernales2017_extrap_config
    C%slid_inv_Bernales2017_hi_scale           = slid_inv_Bernales2017_hi_scale_config
    C%slid_inv_Bernales2017_smooth_r           = slid_inv_Bernales2017_smooth_r_config
    C%slid_inv_Bernales2017_smooth_w           = slid_inv_Bernales2017_smooth_w_config
    C%slid_inv_Bernales2017_tol_diff           = slid_inv_Bernales2017_tol_diff_config
    C%slid_inv_Bernales2017_tol_frac           = slid_inv_Bernales2017_tol_frac_config
    C%slid_inv_Berends2022_tauc                = slid_inv_Berends2022_tauc_config
    C%slid_inv_Berends2022_H0                  = slid_inv_Berends2022_H0_config
    C%slid_inv_Berends2022_u0                  = slid_inv_Berends2022_u0_config
    C%slid_inv_Berends2022_Hi_scale            = slid_inv_Berends2022_Hi_scale_config
    C%slid_inv_Berends2022_u_scale             = slid_inv_Berends2022_u_scale_config
    C%slid_inv_target_velocity_filename        = slid_inv_target_velocity_filename_config

    ! Ice dynamics - calving
    ! ======================

    C%choice_calving_law                       = choice_calving_law_config
    C%calving_threshold_thickness_shelf        = calving_threshold_thickness_shelf_config
    C%calving_threshold_thickness_sheet        = calving_threshold_thickness_sheet_config
    C%max_calving_rounds                       = max_calving_rounds_config
    C%do_remove_shelves                        = do_remove_shelves_config
    C%remove_shelves_larger_than_PD            = remove_shelves_larger_than_PD_config
    C%continental_shelf_calving                = continental_shelf_calving_config
    C%continental_shelf_min_height             = continental_shelf_min_height_config
    C%minimum_ice_thickness                    = minimum_ice_thickness_config

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

    ! Climate
    ! =======

    C%choice_climate_model                     = choice_climate_model_config
    C%choice_idealised_climate                 = choice_idealised_climate_config

    ! NetCDF files containing direct global/regional climate forcing
    C%filename_direct_global_climate           = filename_direct_global_climate_config
    C%filename_direct_regional_climate_NAM     = filename_direct_regional_climate_NAM_config
    C%filename_direct_regional_climate_EAS     = filename_direct_regional_climate_EAS_config
    C%filename_direct_regional_climate_GRL     = filename_direct_regional_climate_GRL_config
    C%filename_direct_regional_climate_ANT     = filename_direct_regional_climate_ANT_config

    ! NetCDF file containing the present-day observed climate (e.g. ERA40)
    C%filename_PD_obs_climate                  = filename_PD_obs_climate_config

    ! GCM snapshots in the matrix_warm_cold option
    C%filename_climate_snapshot_PI             = filename_climate_snapshot_PI_config
    C%filename_climate_snapshot_warm           = filename_climate_snapshot_warm_config
    C%filename_climate_snapshot_cold           = filename_climate_snapshot_cold_config

    C%constant_lapserate                       = constant_lapserate_config

    ! Orbit time and CO2 concentration of the warm and cold snapshots
    C%matrix_high_CO2_level                    = matrix_high_CO2_level_config
    C%matrix_low_CO2_level                     = matrix_low_CO2_level_config
    C%matrix_warm_orbit_time                   = matrix_warm_orbit_time_config
    C%matrix_cold_orbit_time                   = matrix_cold_orbit_time_config

    ! Whether or not to apply a bias correction to the GCM snapshots
    C%climate_matrix_biascorrect_warm          = climate_matrix_biascorrect_warm_config
    C%climate_matrix_biascorrect_cold          = climate_matrix_biascorrect_cold_config

    ! Glacial index
    C%switch_glacial_index                     = switch_glacial_index_config

    ! Iterative adjustment of precipitation and temperature
    C%do_clim_inv                              = do_clim_inv_config
    C%clim_inv_t_start                         = clim_inv_t_start_config
    C%clim_inv_t_end                           = clim_inv_t_end_config
    C%clim_inv_Hb_min                          = clim_inv_Hb_min_config
    C%clim_inv_Hi_max                          = clim_inv_Hi_max_config
    C%clim_inv_window_size                     = clim_inv_window_size_config

    ! Ocean
    ! =====

    C%choice_ocean_model                       = choice_ocean_model_config
    C%choice_idealised_ocean                   = choice_idealised_ocean_config

    ! NetCDF file containing the present-day observed ocean (WOA18) (NetCDF)
    C%filename_PD_obs_ocean                    = filename_PD_obs_ocean_config
    C%name_ocean_temperature_obs               = name_ocean_temperature_obs_config
    C%name_ocean_salinity_obs                  = name_ocean_salinity_obs_config
    C%use_inverted_ocean                       = use_inverted_ocean_config

    ! GCM snapshots in the matrix_warm_cold option
    C%filename_GCM_ocean_snapshot_PI           = filename_GCM_ocean_snapshot_PI_config
    C%filename_GCM_ocean_snapshot_warm         = filename_GCM_ocean_snapshot_warm_config
    C%filename_GCM_ocean_snapshot_cold         = filename_GCM_ocean_snapshot_cold_config
    C%name_ocean_temperature_GCM               = name_ocean_temperature_GCM_config
    C%name_ocean_salinity_GCM                  = name_ocean_salinity_GCM_config

    ! Uniform ocean temperature values used when choice_ocean_model = "uniform_warm_cold"
    C%ocean_temperature_PD                     = ocean_temperature_PD_config
    C%ocean_temperature_cold                   = ocean_temperature_cold_config
    C%ocean_temperature_warm                   = ocean_temperature_warm_config

    ! Parameters used when choice_ocean_model = "matrix_warm_cold"
    C%choice_ocean_vertical_grid               = choice_ocean_vertical_grid_config
    C%ocean_vertical_grid_max_depth            = ocean_vertical_grid_max_depth_config
    C%ocean_regular_grid_dz                    = ocean_regular_grid_dz_config
    C%ocean_extrap_dir                         = ocean_extrap_dir_config
    C%ocean_extrap_res                         = ocean_extrap_res_config
    C%ocean_extrap_Gauss_sigma                 = ocean_extrap_Gauss_sigma_config
    C%ocean_extrap_hires_geo_filename_NAM      = ocean_extrap_hires_geo_filename_NAM_config
    C%ocean_extrap_hires_geo_filename_EAS      = ocean_extrap_hires_geo_filename_EAS_config
    C%ocean_extrap_hires_geo_filename_GRL      = ocean_extrap_hires_geo_filename_GRL_config
    C%ocean_extrap_hires_geo_filename_ANT      = ocean_extrap_hires_geo_filename_ANT_config
    C%ocean_w_tot_hist_averaging_window        = ocean_w_tot_hist_averaging_window_config

    ! Surface mass balance
    ! ====================

    C%choice_SMB_model                         = choice_SMB_model_config
    C%choice_idealised_SMB                     = choice_idealised_SMB_config
    C%SMB_uniform                              = SMB_uniform_config

    ! NetCDF file containing direct global/regional SMB forcing
    C%filename_direct_global_SMB               = filename_direct_global_SMB_config
    C%filename_direct_regional_SMB_NAM         = filename_direct_regional_SMB_NAM_config
    C%filename_direct_regional_SMB_EAS         = filename_direct_regional_SMB_EAS_config
    C%filename_direct_regional_SMB_GRL         = filename_direct_regional_SMB_GRL_config
    C%filename_direct_regional_SMB_ANT         = filename_direct_regional_SMB_ANT_config

    ! Tuning parameters for the IMAU-ITM SMB model
    C%SMB_IMAUITM_choice_init_firn_NAM         = SMB_IMAUITM_choice_init_firn_NAM_config
    C%SMB_IMAUITM_choice_init_firn_EAS         = SMB_IMAUITM_choice_init_firn_EAS_config
    C%SMB_IMAUITM_choice_init_firn_GRL         = SMB_IMAUITM_choice_init_firn_GRL_config
    C%SMB_IMAUITM_choice_init_firn_ANT         = SMB_IMAUITM_choice_init_firn_ANT_config
    C%SMB_IMAUITM_initial_firn_thickness       = SMB_IMAUITM_initial_firn_thickness_config
    C%SMB_IMAUITM_C_abl_constant_NAM           = SMB_IMAUITM_C_abl_constant_NAM_config
    C%SMB_IMAUITM_C_abl_constant_EAS           = SMB_IMAUITM_C_abl_constant_EAS_config
    C%SMB_IMAUITM_C_abl_constant_GRL           = SMB_IMAUITM_C_abl_constant_GRL_config
    C%SMB_IMAUITM_C_abl_constant_ANT           = SMB_IMAUITM_C_abl_constant_ANT_config
    C%SMB_IMAUITM_C_abl_Ts_NAM                 = SMB_IMAUITM_C_abl_Ts_NAM_config
    C%SMB_IMAUITM_C_abl_Ts_EAS                 = SMB_IMAUITM_C_abl_Ts_EAS_config
    C%SMB_IMAUITM_C_abl_Ts_GRL                 = SMB_IMAUITM_C_abl_Ts_GRL_config
    C%SMB_IMAUITM_C_abl_Ts_ANT                 = SMB_IMAUITM_C_abl_Ts_ANT_config
    C%SMB_IMAUITM_C_abl_Q_NAM                  = SMB_IMAUITM_C_abl_Q_NAM_config
    C%SMB_IMAUITM_C_abl_Q_EAS                  = SMB_IMAUITM_C_abl_Q_EAS_config
    C%SMB_IMAUITM_C_abl_Q_GRL                  = SMB_IMAUITM_C_abl_Q_GRL_config
    C%SMB_IMAUITM_C_abl_Q_ANT                  = SMB_IMAUITM_C_abl_Q_ANT_config
    C%SMB_IMAUITM_C_refr_NAM                   = SMB_IMAUITM_C_refr_NAM_config
    C%SMB_IMAUITM_C_refr_EAS                   = SMB_IMAUITM_C_refr_EAS_config
    C%SMB_IMAUITM_C_refr_GRL                   = SMB_IMAUITM_C_refr_GRL_config
    C%SMB_IMAUITM_C_refr_ANT                   = SMB_IMAUITM_C_refr_ANT_config

    ! IMAU-ITM SMB model inversion
    C%do_SMB_IMAUITM_inversion                 = do_SMB_IMAUITM_inversion_config
    C%SMB_IMAUITM_inv_choice_init_C            = SMB_IMAUITM_inv_choice_init_C_config
    C%SMB_IMAUITM_inv_scale                    = SMB_IMAUITM_inv_scale_config
    C%SMB_IMAUITM_inv_C_abl_constant_min       = SMB_IMAUITM_inv_C_abl_constant_min_config
    C%SMB_IMAUITM_inv_C_abl_constant_max       = SMB_IMAUITM_inv_C_abl_constant_max_config
    C%SMB_IMAUITM_inv_C_abl_Ts_min             = SMB_IMAUITM_inv_C_abl_Ts_min_config
    C%SMB_IMAUITM_inv_C_abl_Ts_max             = SMB_IMAUITM_inv_C_abl_Ts_max_config
    C%SMB_IMAUITM_inv_C_abl_Q_min              = SMB_IMAUITM_inv_C_abl_Q_min_config
    C%SMB_IMAUITM_inv_C_abl_Q_max              = SMB_IMAUITM_inv_C_abl_Q_max_config
    C%SMB_IMAUITM_inv_C_refr_min               = SMB_IMAUITM_inv_C_refr_min_config
    C%SMB_IMAUITM_inv_C_refr_max               = SMB_IMAUITM_inv_C_refr_max_config

    ! ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
    ! ==============================================================

    C%ISMIP_forcing_filename_baseline          = ISMIP_forcing_filename_baseline_config
    C%ISMIP_forcing_foldername_aSMB            = ISMIP_forcing_foldername_aSMB_config
    C%ISMIP_forcing_basefilename_aSMB          = ISMIP_forcing_basefilename_aSMB_config
    C%ISMIP_forcing_foldername_dSMBdz          = ISMIP_forcing_foldername_dSMBdz_config
    C%ISMIP_forcing_basefilename_dSMBdz        = ISMIP_forcing_basefilename_dSMBdz_config
    C%ISMIP_forcing_foldername_aST             = ISMIP_forcing_foldername_aST_config
    C%ISMIP_forcing_basefilename_aST           = ISMIP_forcing_basefilename_aST_config
    C%ISMIP_forcing_foldername_dSTdz           = ISMIP_forcing_foldername_dSTdz_config
    C%ISMIP_forcing_basefilename_dSTdz         = ISMIP_forcing_basefilename_dSTdz_config

    ! Basal mass balance - sub-shelf melt
    ! ===================================

    C%choice_BMB_shelf_model                   = choice_BMB_shelf_model_config
    C%choice_idealised_BMB_shelf               = choice_idealised_BMB_shelf_config
    C%choice_BMB_sheet_model                   = choice_BMB_sheet_model_config
    C%BMB_shelf_uniform                        = BMB_shelf_uniform_config
    C%BMB_sheet_uniform                        = BMB_sheet_uniform_config
    C%choice_BMB_subgrid                       = choice_BMB_subgrid_config
    C%BMB_max                                  = BMB_max_config
    C%BMB_min                                  = BMB_min_config

    C%do_ocean_temperature_inversion           = do_ocean_temperature_inversion_config
    C%ocean_temperature_inv_t_start            = ocean_temperature_inv_t_start_config
    C%ocean_temperature_inv_t_end              = ocean_temperature_inv_t_end_config
    C%T_base_window_size                       = T_base_window_size_config

    C%BMB_inv_use_restart_field                = BMB_inv_use_restart_field_config
    C%BMB_inv_scale_shelf                      = BMB_inv_scale_shelf_config
    C%BMB_inv_scale_ocean                      = BMB_inv_scale_ocean_config

    C%choice_basin_scheme_NAM                  = choice_basin_scheme_NAM_config
    C%choice_basin_scheme_EAS                  = choice_basin_scheme_EAS_config
    C%choice_basin_scheme_GRL                  = choice_basin_scheme_GRL_config
    C%choice_basin_scheme_ANT                  = choice_basin_scheme_ANT_config
    C%filename_basins_NAM                      = filename_basins_NAM_config
    C%filename_basins_EAS                      = filename_basins_EAS_config
    C%filename_basins_GRL                      = filename_basins_GRL_config
    C%filename_basins_ANT                      = filename_basins_ANT_config
    C%do_merge_basins_ANT                      = do_merge_basins_ANT_config
    C%do_merge_basins_GRL                      = do_merge_basins_GRL_config

    C%choice_BMB_shelf_amplification           = choice_BMB_shelf_amplification_config
    C%basin_BMB_amplification_n_ANT            = basin_BMB_amplification_n_ANT_config
    ALLOCATE( C%basin_BMB_amplification_factor_ANT( C%basin_BMB_amplification_n_ANT))
    C%basin_BMB_amplification_factor_ANT       = basin_BMB_amplification_factor_ANT_config( 1:C%basin_BMB_amplification_n_ANT)
    C%basin_BMB_amplification_n_GRL            = basin_BMB_amplification_n_GRL_config
    ALLOCATE( C%basin_BMB_amplification_factor_GRL( C%basin_BMB_amplification_n_GRL))
    C%basin_BMB_amplification_factor_GRL       = basin_BMB_amplification_factor_GRL_config( 1:C%basin_BMB_amplification_n_GRL)

    ! Parameters for the three simple melt parameterisations from Favier et al. (2019)
    C%BMB_Favier2019_lin_GammaT                = BMB_Favier2019_lin_GammaT_config
    C%BMB_Favier2019_quad_GammaT               = BMB_Favier2019_quad_GammaT_config
    C%BMB_Favier2019_Mplus_GammaT              = BMB_Favier2019_Mplus_GammaT_config

    ! Parameters for the Lazeroms et al. (2018) plume-parameterisation BMB model
    C%BMB_Lazeroms2018_GammaT                  = BMB_Lazeroms2018_GammaT_config
    C%BMB_Lazeroms2018_find_GL_scheme          = BMB_Lazeroms2018_find_GL_scheme_config

    ! Parameters for the PICO BMB model
    C%BMB_PICO_nboxes                          = BMB_PICO_nboxes_config
    C%BMB_PICO_GammaTstar                      = BMB_PICO_GammaTstar_config

    ! Parameters for the ANICE_legacy sub-shelf melt model
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

    ! Englacial isotope tracing
    ! ========================

    C%choice_ice_isotopes_model                = choice_ice_isotopes_model_config
    C%uniform_ice_d18O                         = uniform_ice_d18O_config

    ! Sea level and GIA
    ! =================

    C%do_ocean_floodfill                       = do_ocean_floodfill_config
    C%choice_sealevel_model                    = choice_sealevel_model_config
    C%fixed_sealevel                           = fixed_sealevel_config
    C%initial_guess_sealevel                   = initial_guess_sealevel_config
    C%filename_sealevel_record                 = filename_sealevel_record_config
    C%sealevel_record_length                   = sealevel_record_length_config

    C%choice_GIA_model                         = choice_GIA_model_config
    C%ELRA_lithosphere_flex_rigidity           = ELRA_lithosphere_flex_rigidity_config
    C%ELRA_bedrock_relaxation_time             = ELRA_bedrock_relaxation_time_config
    C%ELRA_mantle_density                      = ELRA_mantle_density_config

    ! SELEN
    ! =====

    C%SELEN_run_at_t_start                     = SELEN_run_at_t_start_config
    C%SELEN_n_TDOF_iterations                  = SELEN_n_TDOF_iterations_config
    C%SELEN_n_recursion_iterations             = SELEN_n_recursion_iterations_config
    C%SELEN_use_rotational_feedback            = SELEN_use_rotational_feedback_config
    C%SELEN_n_harmonics                        = SELEN_n_harmonics_config
    C%SELEN_display_progress                   = SELEN_display_progress_config

    C%SELEN_dir                                = SELEN_dir_config
    C%SELEN_global_topo_filename               = SELEN_global_topo_filename_config
    C%SELEN_TABOO_init_filename                = SELEN_TABOO_init_filename_config
    C%SELEN_LMJ_VALUES_filename                = SELEN_LMJ_VALUES_filename_config

    C%SELEN_irreg_time_n                       = SELEN_irreg_time_n_config
    ALLOCATE( C%SELEN_irreg_time_window( C%SELEN_irreg_time_n))
    C%SELEN_irreg_time_window                  = SELEN_irreg_time_window_config( 1:C%SELEN_irreg_time_n)

    C%SELEN_lith_thickness                     = SELEN_lith_thickness_config
    C%SELEN_visc_n                             = SELEN_visc_n_config
    ALLOCATE( C%SELEN_visc_prof( C%SELEN_visc_n))
    C%SELEN_visc_prof      = SELEN_visc_prof_config( 1:C%SELEN_visc_n)

    C%SELEN_TABOO_CDE                          = SELEN_TABOO_CDE_config
    C%SELEN_TABOO_TLOVE                        = SELEN_TABOO_TLOVE_config
    C%SELEN_TABOO_DEG1                         = SELEN_TABOO_DEG1_config
    C%SELEN_TABOO_RCMB                         = SELEN_TABOO_RCMB_config

    ! Fill in some derived values
    C%SELEN_jmax       = (C%SELEN_n_harmonics + 1) * (C%SELEN_n_harmonics + 2) / 2
    C%SELEN_reg_time_n = CEILING(MAX(1._dp, SUM(C%SELEN_irreg_time_window( 1:C%SELEN_irreg_time_n))) * 1000._dp / C%dt_SELEN)

    CALL initialize_TABOO_config

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

  END SUBROUTINE copy_variables_to_struct

  SUBROUTINE initialize_TABOO_config

    ! Taboo settings
    C_TABOO%IMODE            = 1
    C_TABOO%NV               = C%SELEN_visc_n
    ALLOCATE( C_TABOO%VSC( C%SELEN_visc_n))
    C_TABOO%VSC              = C%SELEN_visc_prof
    C_TABOO%CDE              = C%SELEN_TABOO_CDE
    C_TABOO%TLOVE            = C%SELEN_TABOO_TLOVE
    C_TABOO%DEG1             = C%SELEN_TABOO_DEG1
    C_TABOO%LTH              = C%SELEN_lith_thickness
    C_TABOO%RCMB             = C%SELEN_TABOO_RCMB

  END SUBROUTINE initialize_TABOO_config

  SUBROUTINE get_procedural_output_dir_name( output_dir)
    ! Generate a procedural output directory for the current date (e.g. results_20210721_001)
    ! Keep increasing the counter at the end until a directory is available.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(20),                       INTENT(INOUT) :: output_dir

    ! Local variables:
    INTEGER,  DIMENSION(8)                             :: values
    LOGICAL                                            :: ex

    CALL date_and_time(VALUES=values)

    ! Get proper year (assume we're still in the 21st century...)
    output_dir(1:10) = 'results_20'
    SELECT CASE( FLOOR(REAL(values(1))/10._dp)-200)
    CASE(0)
      output_dir(11:11) = '0'
    CASE(1)
      output_dir(11:11) = '1'
    CASE(2)
      output_dir(11:11) = '2'
    CASE(3)
      output_dir(11:11) = '3'
    CASE(4)
      output_dir(11:11) = '4'
    CASE(5)
      output_dir(11:11) = '5'
    CASE(6)
      output_dir(11:11) = '6'
    CASE(7)
      output_dir(11:11) = '7'
    CASE(8)
      output_dir(11:11) = '8'
    CASE(9)
      output_dir(11:11) = '9'
    CASE DEFAULT
      WRITE(0,*) 'get_procedural_output_dir: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( MOD(values(1),10))
    CASE(0)
      output_dir(12:12) = '0'
    CASE(1)
      output_dir(12:12) = '1'
    CASE(2)
      output_dir(12:12) = '2'
    CASE(3)
      output_dir(12:12) = '3'
    CASE(4)
      output_dir(12:12) = '4'
    CASE(5)
      output_dir(12:12) = '5'
    CASE(6)
      output_dir(12:12) = '6'
    CASE(7)
      output_dir(12:12) = '7'
    CASE(8)
      output_dir(12:12) = '8'
    CASE(9)
      output_dir(12:12) = '9'
    CASE DEFAULT
      WRITE(0,*) 'get_procedural_output_dir: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( values(2))
    CASE(1)
      output_dir(13:14) = '01'
    CASE(2)
      output_dir(13:14) = '02'
    CASE(3)
      output_dir(13:14) = '03'
    CASE(4)
      output_dir(13:14) = '04'
    CASE(5)
      output_dir(13:14) = '05'
    CASE(6)
      output_dir(13:14) = '06'
    CASE(7)
      output_dir(13:14) = '07'
    CASE(8)
      output_dir(13:14) = '08'
    CASE(9)
      output_dir(13:14) = '09'
    CASE(10)
      output_dir(13:14) = '10'
    CASE(11)
      output_dir(13:14) = '11'
    CASE(12)
      output_dir(13:14) = '12'
    CASE DEFAULT
      WRITE(0,*) 'get_procedural_output_dir: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( FLOOR(REAL(values(3))/10._dp))
    CASE(0)
      output_dir(15:15) = '0'
    CASE(1)
      output_dir(15:15) = '1'
    CASE(2)
      output_dir(15:15) = '2'
    CASE(3)
      output_dir(15:15) = '3'
    CASE DEFAULT
      WRITE(0,*) 'get_procedural_output_dir: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( MOD(values(3),10))
    CASE(0)
      output_dir(16:16) = '0'
    CASE(1)
      output_dir(16:16) = '1'
    CASE(2)
      output_dir(16:16) = '2'
    CASE(3)
      output_dir(16:16) = '3'
    CASE(4)
      output_dir(16:16) = '4'
    CASE(5)
      output_dir(16:16) = '5'
    CASE(6)
      output_dir(16:16) = '6'
    CASE(7)
      output_dir(16:16) = '7'
    CASE(8)
      output_dir(16:16) = '8'
    CASE(9)
      output_dir(16:16) = '9'
    CASE DEFAULT
      WRITE(0,*) 'get_procedural_output_dir: ERROR retrieving date and time!'
    END SELECT

    output_dir(17:20) = '_001'

    INQUIRE( FILE = TRIM(output_dir)//'/.', EXIST=ex )

    DO WHILE (ex)

     IF      (output_dir(20:20) == '0') THEN
       output_dir(20:20) = '1'
     ELSE IF (output_dir(20:20) == '1') THEN
       output_dir(20:20) = '2'
     ELSE IF (output_dir(20:20) == '2') THEN
       output_dir(20:20) = '3'
     ELSE IF (output_dir(20:20) == '3') THEN
       output_dir(20:20) = '4'
     ELSE IF (output_dir(20:20) == '4') THEN
       output_dir(20:20) = '5'
     ELSE IF (output_dir(20:20) == '5') THEN
       output_dir(20:20) = '6'
     ELSE IF (output_dir(20:20) == '6') THEN
       output_dir(20:20) = '7'
     ELSE IF (output_dir(20:20) == '7') THEN
       output_dir(20:20) = '8'
     ELSE IF (output_dir(20:20) == '8') THEN
       output_dir(20:20) = '9'
     ELSE IF (output_dir(20:20) == '9') THEN
       output_dir(20:20) = '0'

       IF      (output_dir(19:19) == '0') THEN
         output_dir(19:19) = '1'
       ELSE IF (output_dir(19:19) == '1') THEN
         output_dir(19:19) = '2'
       ELSE IF (output_dir(19:19) == '2') THEN
         output_dir(19:19) = '3'
       ELSE IF (output_dir(19:19) == '3') THEN
         output_dir(19:19) = '4'
       ELSE IF (output_dir(19:19) == '4') THEN
         output_dir(19:19) = '5'
       ELSE IF (output_dir(19:19) == '5') THEN
         output_dir(19:19) = '6'
       ELSE IF (output_dir(19:19) == '6') THEN
         output_dir(19:19) = '7'
       ELSE IF (output_dir(19:19) == '7') THEN
         output_dir(19:19) = '8'
       ELSE IF (output_dir(19:19) == '8') THEN
         output_dir(19:19) = '9'
       ELSE IF (output_dir(19:19) == '9') THEN
         output_dir(19:19) = '0'

         IF      (output_dir(18:18) == '0') THEN
           output_dir(18:18) = '1'
         ELSE IF (output_dir(18:18) == '1') THEN
           output_dir(18:18) = '2'
         ELSE IF (output_dir(18:18) == '2') THEN
           output_dir(18:18) = '3'
         ELSE IF (output_dir(18:18) == '3') THEN
           output_dir(18:18) = '4'
         ELSE IF (output_dir(18:18) == '4') THEN
           output_dir(18:18) = '5'
         ELSE IF (output_dir(18:18) == '5') THEN
           output_dir(18:18) = '6'
         ELSE IF (output_dir(18:18) == '6') THEN
           output_dir(18:18) = '7'
         ELSE IF (output_dir(18:18) == '7') THEN
           output_dir(18:18) = '8'
         ELSE IF (output_dir(18:18) == '8') THEN
           output_dir(18:18) = '9'
         ELSE IF (output_dir(18:18) == '9') THEN
           output_dir(18:18) = '0'
         END IF

       END IF

     END IF

     INQUIRE( FILE=TRIM(output_dir)//'/.', EXIST=ex )

    END DO

  END SUBROUTINE get_procedural_output_dir_name

  SUBROUTINE write_total_model_time_to_screen( tstart, tstop)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: tstart, tstop

    ! Local variables
    REAL(dp)                                           :: dt
    INTEGER                                            :: nr, ns, nm, nh, nd

    dt = tstop - tstart

    ns = CEILING(dt)

    nr = MOD(ns, 60*60*24)
    nd = (ns - nr) / (60*60*24)
    ns = ns - (nd*60*60*24)

    nr = MOD(ns, 60*60)
    nh = (ns - nr) / (60*60)
    ns = ns - (nh*60*60)

    nr = MOD(ns, 60)
    nm = (ns - nr) / (60)
    ns = ns - (nm*60)

    WRITE(0,'(A)') ''
    WRITE(0,'(A)') ' ================================================================================'
    WRITE(0,'(A,I2,A,I2,A,I2,A,I2,A)') ' ===== Simulation finished in ', nd, ' days, ', nh, ' hours, ', nm, ' minutes and ', ns, ' seconds! ====='
    WRITE(0,'(A)') ' ================================================================================'
    WRITE(0,'(A)') ''

  END SUBROUTINE write_total_model_time_to_screen

! ===== Extended error messaging / debugging system =====
! =======================================================

  SUBROUTINE init_routine( routine_name, do_track_resource_use)
    ! Initialise an IMAU-ICE subroutine; update the routine path

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: routine_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: do_track_resource_use

    ! Local variables:
    INTEGER                                            :: len_path_tot, len_path_used, len_name
    INTEGER                                            :: ierr, cerr
    INTEGER                                            :: i
    LOGICAL                                            :: do_track_resource_use_loc

    ! Check if routine_name has enough memory
    len_path_tot  = LEN(      routine_path)
    len_path_used = LEN_TRIM( routine_path)
    len_name      = LEN_TRIM( routine_name)

    IF (len_path_used + 1 + len_name > len_path_tot) THEN
      WRITE(0,*) 'init_routine - ERROR: routine_path = "', TRIM( routine_path), '", no more space to append routine_name = "', TRIM( routine_name), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Append this routine to the routine path
    routine_path = TRIM( routine_path) // '/' // TRIM( routine_name)

    ! Check if resource use for this subroutine should be tracked
    ! (especially for the NetCDF routines we don't want to do this, as there are
    ! a great many of them and the resource tracker output file will become annoyingly big)

    IF (PRESENT( do_track_resource_use)) THEN
      do_track_resource_use_loc = do_track_resource_use
    ELSE
      do_track_resource_use_loc = .TRUE.
    END IF

    IF (do_track_resource_use_loc) THEN

      ! Initialise the computation time tracker
      CALL find_subroutine_in_resource_tracker( i)
      resource_tracker( i)%tstart = MPI_WTIME()

      ! Check maximum MPI window at the start of the routine
      resource_tracker( i)%n_MPI_windows_init = n_MPI_windows

    ELSE

      routine_path = TRIM( routine_path) // '_NOTRACK'

    END IF


  END SUBROUTINE init_routine

  SUBROUTINE finalise_routine( routine_name, n_extra_windows_expected)
    ! Finalise; remove the current routine name from the routine path

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: routine_name
    INTEGER,  INTENT(IN), OPTIONAL                     :: n_extra_windows_expected

    ! Local variables:
    LOGICAL                                            :: do_track_resource_use
    INTEGER                                            :: len_path_tot, i, ii
    INTEGER                                            :: ierr, cerr
    REAL(dp)                                           :: dt
    INTEGER                                            :: n_extra_windows_expected_loc, n_extra_windows_found

    ! Check if resource use should be tracked for this subroutine
    i = INDEX( routine_path, '_NOTRACK')
    IF ( i == 0) THEN
      do_track_resource_use = .TRUE.
    ELSE
      do_track_resource_use = .FALSE.
    END IF

    IF (do_track_resource_use) THEN
      ! Resource use for this subroutine should be tracked

      ! Add computation time to the resource tracker
      CALL find_subroutine_in_resource_tracker( i)
      dt = MPI_WTIME() - resource_tracker( i)%tstart
      resource_tracker( i)%tcomp = resource_tracker( i)%tcomp + dt

      ! Check maximum MPI window at the end of the routine
      resource_tracker( i)%n_MPI_windows_final = n_MPI_windows

      ! If it is larger than expected, warn that there might be a memory leak
      n_extra_windows_expected_loc = 0
      IF (PRESENT( n_extra_windows_expected)) n_extra_windows_expected_loc = n_extra_windows_expected
      n_extra_windows_found = resource_tracker( i)%n_MPI_windows_final - resource_tracker( i)%n_MPI_windows_init

      ii = INDEX( routine_path, 'UFEMISM_program/initialise_')
      IF (ii == 0 .AND. n_extra_windows_found > n_extra_windows_expected_loc) THEN
        ! This subroutine has more memory allocated at the start than at the beginning.
        CALL warning('more memory was allocated and not freed than expected; possible memory leak! (expected {int_01} extra windows, found {int_02})', &
          int_01 = n_extra_windows_expected_loc, int_02 = n_extra_windows_found)
      END IF

      ! Find where in the string exactly the current routine name is located
      i = INDEX( routine_path, routine_name)

      IF (i == 0) THEN
        WRITE(0,*) 'finalise_routine - ERROR: routine_name = "', TRIM( routine_name), '" not found in routine_path = "', TRIM( routine_path), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      ! Remove the current routine name from the routine path
      len_path_tot = LEN( routine_path)
      routine_path( i-1:len_path_tot) = ' '

    ELSE ! IF (do_track_resource_use) THEN
      ! Resource use for this subroutine should not be tracked

      ! Find where in the string exactly the current routine name is located
      i = INDEX( routine_path, TRIM( routine_name) // '_NOTRACK')

      IF (i == 0) THEN
        WRITE(0,*) 'finalise_routine - ERROR: routine_name = "', TRIM( routine_name), '" not found in routine_path = "', TRIM( routine_path), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      ! Remove the current routine name from the routine path
      len_path_tot = LEN( routine_path)
      routine_path( i-1:len_path_tot) = ' '

    END IF ! IF (do_track_resource_use) THEN

  END SUBROUTINE finalise_routine

  SUBROUTINE find_subroutine_in_resource_tracker( i)
    ! Find the current subroutine in the resource tracker. If it's not there yet, add it.

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                             INTENT(OUT)   :: i

    ! Local variables:
    INTEGER                                            :: n

    n = SIZE( resource_tracker)

    DO i = 1, n
      IF     (resource_tracker( i)%routine_path == routine_path) THEN
        ! The current subroutine is listed at this position in the resource tracker
        RETURN
      ELSEIF (resource_tracker( i)%routine_path == 'subroutine_placeholder') THEN
        ! We've checked all listed subroutines and haven't found the current one; add it
        resource_tracker( i)%routine_path = routine_path
        RETURN
      END IF
    END DO

    ! If we've reached this point, then the resource tracker is overflowing
    CALL crash('Resource tracker overflows! Allocate more memory for it in initialise_model_configuration.')

  END SUBROUTINE find_subroutine_in_resource_tracker

  SUBROUTINE reset_resource_tracker
    ! Reset the computation times and maximum memory use for all subroutines in the resource tracker

    IMPLICIT NONE

    ! Local variables:
    INTEGER                                            :: i,n

    mem_use_tot_max = 0._dp

    n = SIZE( resource_tracker)

    DO i = 1, n
      resource_tracker( i)%tstart      = 0._dp
      resource_tracker( i)%tcomp       = 0._dp
      resource_tracker( i)%mem_use_max = 0._dp
    END DO

  END SUBROUTINE reset_resource_tracker

  SUBROUTINE crash( err_msg, int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
                              dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    ! Crash the model, write the error message to the screen

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: err_msg
    INTEGER,  INTENT(IN), OPTIONAL                     :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    REAL(dp), INTENT(IN), OPTIONAL                     ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    ! Local variables:
    CHARACTER(LEN=1024)                                :: err_msg_loc
    INTEGER                                            :: ierr, cerr, pari, parn, nc
    CHARACTER(LEN=9)                                   :: fmt
    CHARACTER(LEN=:), ALLOCATABLE                      :: process_str

    ! Get local, edit-able copy of error message string
    err_msg_loc = err_msg

    ! Get rank of current process and total number of processes
    ! (needed because the configuration_module cannot access the par structure)
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, pari, ierr)
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, parn, ierr)

    ! Set the process string (e.g. "05/16")
    IF     (parn < 10) THEN
      nc = 1
    ELSEIF (parn < 100) THEN
      nc = 2
    ELSEIF (parn < 1000) THEN
      nc = 3
    ELSEIF (parn < 10000) THEN
      nc = 4
    ELSE
      nc = 5
    END IF

    WRITE( fmt,'(A,I1,A,I1,A)') '(I', nc, ',A,I', nc, ')'
    ALLOCATE(CHARACTER(2*nc+1) :: process_str)
    WRITE( process_str,fmt) pari, '/', parn

    ! Insert numbers into string if needed
    IF (PRESENT( int_01)) CALL insert_val_into_string_int( err_msg_loc, '{int_01}', int_01)
    IF (PRESENT( int_02)) CALL insert_val_into_string_int( err_msg_loc, '{int_02}', int_02)
    IF (PRESENT( int_03)) CALL insert_val_into_string_int( err_msg_loc, '{int_03}', int_03)
    IF (PRESENT( int_04)) CALL insert_val_into_string_int( err_msg_loc, '{int_04}', int_04)
    IF (PRESENT( int_05)) CALL insert_val_into_string_int( err_msg_loc, '{int_05}', int_05)
    IF (PRESENT( int_06)) CALL insert_val_into_string_int( err_msg_loc, '{int_06}', int_06)
    IF (PRESENT( int_07)) CALL insert_val_into_string_int( err_msg_loc, '{int_07}', int_07)
    IF (PRESENT( int_08)) CALL insert_val_into_string_int( err_msg_loc, '{int_08}', int_08)
    IF (PRESENT( int_09)) CALL insert_val_into_string_int( err_msg_loc, '{int_09}', int_09)
    IF (PRESENT( int_10)) CALL insert_val_into_string_int( err_msg_loc, '{int_10}', int_10)

    IF (PRESENT( dp_01 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_01}' , dp_01 )
    IF (PRESENT( dp_02 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_02}' , dp_02 )
    IF (PRESENT( dp_03 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_03}' , dp_03 )
    IF (PRESENT( dp_04 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_04}' , dp_04 )
    IF (PRESENT( dp_05 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_05}' , dp_05 )
    IF (PRESENT( dp_06 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_06}' , dp_06 )
    IF (PRESENT( dp_07 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_07}' , dp_07 )
    IF (PRESENT( dp_08 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_08}' , dp_08 )
    IF (PRESENT( dp_09 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_09}' , dp_09 )
    IF (PRESENT( dp_10 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_10}' , dp_10 )

    ! Write the error to the screen
    WRITE(0,'(A,A,A,A,A,A)') colour_string(' ERROR: ' // TRIM( err_msg_loc),'red') // ' in ' // colour_string( TRIM(routine_path),'light blue') // &
      ' on process ', colour_string( process_str,'light blue'), ' (0 = master)'

    ! Stop the program
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

  END SUBROUTINE crash

  SUBROUTINE warning( err_msg, int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
                                dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    ! Write the warning message to the screen, but don't crash the model

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: err_msg
    INTEGER,  INTENT(IN), OPTIONAL                     :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    REAL(dp), INTENT(IN), OPTIONAL                     ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    ! Local variables:
    CHARACTER(LEN=1024)                                :: err_msg_loc
    INTEGER                                            :: ierr, pari, parn, nc
    CHARACTER(LEN=9)                                   :: fmt
    CHARACTER(LEN=:), ALLOCATABLE                      :: process_str

    ! Get local, edit-able copy of error message string
    err_msg_loc = err_msg

    ! Get rank of current process and total number of processes
    ! (needed because the configuration_module cannot access the par structure)
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, pari, ierr)
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, parn, ierr)

    ! Set the process string (e.g. "05/16")
    IF     (parn < 10) THEN
      nc = 1
    ELSEIF (parn < 100) THEN
      nc = 2
    ELSEIF (parn < 1000) THEN
      nc = 3
    ELSEIF (parn < 10000) THEN
      nc = 4
    ELSE
      nc = 5
    END IF

    WRITE( fmt,'(A,I1,A,I1,A)') '(I', nc, ',A,I', nc, ')'
    ALLOCATE(CHARACTER(2*nc+1) :: process_str)
    WRITE( process_str,fmt) pari, '/', parn

    ! Insert numbers into string if needed
    IF (PRESENT( int_01)) CALL insert_val_into_string_int( err_msg_loc, '{int_01}', int_01)
    IF (PRESENT( int_02)) CALL insert_val_into_string_int( err_msg_loc, '{int_02}', int_02)
    IF (PRESENT( int_03)) CALL insert_val_into_string_int( err_msg_loc, '{int_03}', int_03)
    IF (PRESENT( int_04)) CALL insert_val_into_string_int( err_msg_loc, '{int_04}', int_04)
    IF (PRESENT( int_05)) CALL insert_val_into_string_int( err_msg_loc, '{int_05}', int_05)
    IF (PRESENT( int_06)) CALL insert_val_into_string_int( err_msg_loc, '{int_06}', int_06)
    IF (PRESENT( int_07)) CALL insert_val_into_string_int( err_msg_loc, '{int_07}', int_07)
    IF (PRESENT( int_08)) CALL insert_val_into_string_int( err_msg_loc, '{int_08}', int_08)
    IF (PRESENT( int_09)) CALL insert_val_into_string_int( err_msg_loc, '{int_09}', int_09)
    IF (PRESENT( int_10)) CALL insert_val_into_string_int( err_msg_loc, '{int_10}', int_10)

    IF (PRESENT( dp_01 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_01}' , dp_01 )
    IF (PRESENT( dp_02 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_02}' , dp_02 )
    IF (PRESENT( dp_03 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_03}' , dp_03 )
    IF (PRESENT( dp_04 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_04}' , dp_04 )
    IF (PRESENT( dp_05 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_05}' , dp_05 )
    IF (PRESENT( dp_06 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_06}' , dp_06 )
    IF (PRESENT( dp_07 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_07}' , dp_07 )
    IF (PRESENT( dp_08 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_08}' , dp_08 )
    IF (PRESENT( dp_09 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_09}' , dp_09 )
    IF (PRESENT( dp_10 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_10}' , dp_10 )

    ! Write the error to the screen
    WRITE(0,'(A,A,A,A,A,A)') colour_string(' WARNING: ' // TRIM( err_msg_loc),'yellow') // ' in ' // colour_string( TRIM(routine_path),'light blue') // &
      ' on process ', colour_string( process_str,'light blue'), ' (0 = master)'

    ! Clean up after yourself
    DEALLOCATE( process_str)

  END SUBROUTINE warning

  FUNCTION colour_string( str, col) RESULT( str_col)
    ! Add colour to a string for writing to the terminal

    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: str, col

    ! Result variables:
    CHARACTER(LEN=:), ALLOCATABLE                      :: str_col

    ALLOCATE(CHARACTER(LEN(str)+9) :: str_col)       ! The +9 is just enough to store the color characters

    ! The 91m gives red, 0m sets the default back
    ! Available colors: 90:gray, 91:red, 92:green, 93:yellow, 94:blue, 95:pink, 96:light blue
    IF     (col == 'gray') THEN
      str_col = achar(27)//'[90m'//str//achar(27)//'[0m'
    ELSEIF (col == 'red') THEN
      str_col = achar(27)//'[91m'//str//achar(27)//'[0m'
    ELSEIF (col == 'green') THEN
      str_col = achar(27)//'[92m'//str//achar(27)//'[0m'
    ELSEIF (col == 'yellow') THEN
      str_col = achar(27)//'[93m'//str//achar(27)//'[0m'
    ELSEIF (col == 'blue') THEN
      str_col = achar(27)//'[94m'//str//achar(27)//'[0m'
    ELSEIF (col == 'pink') THEN
      str_col = achar(27)//'[95m'//str//achar(27)//'[0m'
    ELSEIF (col == 'light blue') THEN
      str_col = achar(27)//'[96m'//str//achar(27)//'[0m'
    ELSE
      WRITE(0,*) ''
    END IF

  END FUNCTION colour_string

  SUBROUTINE insert_val_into_string_int( str, marker, val)
    ! Replace marker in str with val (where val is an integer)
    !
    ! Example: str    = 'Johnny has {int_01} apples.'
    !          marker = '{int_01}'
    !          val    = 5
    !
    ! This returns: str = 'Johnny has 5 apples'

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(INOUT) :: str
    CHARACTER(LEN=*),                    INTENT(IN)    :: marker
    INTEGER,                             INTENT(IN)    :: val

    ! Local variables:
    INTEGER                                            :: ci
    INTEGER                                            :: nc
    CHARACTER(LEN=4)                                   :: fmt
    CHARACTER(LEN=:), ALLOCATABLE                      :: val_str
    INTEGER                                            :: len_str, len_marker

    ! Find position ci in str where i_str occurs
    ci = INDEX( str, marker)

    ! Safety
    IF (ci == 0) CALL crash('insert_val_into_string_int: couldnt find marker "' // TRIM( marker) // '" in string "' // TRIM( str) // '"!')

    ! Write val to a string
    IF     (ABS( val) < 10) THEN
      nc = 1
    ELSEIF (ABS( val) < 100) THEN
      nc = 2
    ELSEIF (ABS( val) < 1000) THEN
      nc = 3
    ELSEIF (ABS( val) < 10000) THEN
      nc = 4
    ELSEIF (ABS( val) < 100000) THEN
      nc = 5
    ELSEIF (ABS( val) < 1000000) THEN
      nc = 6
    ELSEIF (ABS( val) < 10000000) THEN
      nc = 7
    ELSEIF (ABS( val) < 100000000) THEN
      nc = 8
    ELSE
      nc = 9
    END IF
    ! Add room for a minus sign if needed
    IF (val < 0) nc = nc + 1

    WRITE( fmt,'(A,I1,A)') '(I', nc, ')'
    ALLOCATE(CHARACTER(nc) :: val_str)
    WRITE( val_str,fmt) val

    ! Find total string length right now
    len_str    = LEN( str)
    len_marker = LEN( marker)

    ! Insert the integer string into the string
    str = str(1:ci-1) // val_str // str(ci+len_marker:len_str)

    ! Clean up after yourself
    DEALLOCATE( val_str)

  END SUBROUTINE insert_val_into_string_int

  SUBROUTINE insert_val_into_string_dp( str, marker, val)
    ! Replace marker in str with val (where val is a double-precision number)
    !
    ! Example: str    = 'Johnny weighs {dp_01} kg.'
    !          marker = '{dp_01}'
    !          val    = 57.098
    !
    ! This returns: str = 'Johnny weighs 57.098 kg'

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(INOUT) :: str
    CHARACTER(LEN=*),                    INTENT(IN)    :: marker
    REAL(dp),                            INTENT(IN)    :: val

    ! Local variables:
    INTEGER                                            :: ci
    CHARACTER(LEN=11)                                  :: val_str
    INTEGER                                            :: len_str, len_marker

    ! Find position ci in str where i_str occurs
    ci = INDEX( str, marker)

    ! Safety
    IF (ci == 0) CALL crash('insert_val_into_string_dp: couldnt find marker "' // TRIM( marker) // '" in string "' // TRIM( str) // '"!')

    ! Write val to a string
    WRITE( val_str,'(E11.5)') val

    ! Find total string length right now
    len_str    = LEN( str)
    len_marker = LEN( marker)

    ! Insert the integer string into the string
    str = str(1:ci-1) // val_str // str(ci+len_marker:len_str)

  END SUBROUTINE insert_val_into_string_dp

  SUBROUTINE capitalise_string( str)
    ! Changes a string which contains lower case letters to a string with upper case letters
    ! Useful for case-insensitive string comparison

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                   INTENT(INOUT) :: str

    ! Local variables:
    INTEGER                                             :: i, index_cap

    CHARACTER(26), PARAMETER                            :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(26), PARAMETER                            :: low = 'abcdefghijklmnopqrstuvwxyz'

    DO i = 1, LEN_TRIM( str)
      index_cap = INDEX( low, str( i:i))
      IF (index_cap > 0) str( i:i) = cap( index_cap:index_cap)
    END DO

  END SUBROUTINE capitalise_string

  SUBROUTINE remove_leading_spaces( str)
    ! Changes a string which contains lower case letters to a string with upper case letters
    ! Useful for case-insensitive string comparison

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                   INTENT(INOUT) :: str

    ! Local variables:
    INTEGER                                             :: lstr

    DO WHILE (str( 1:1) == ' ' .AND. LEN_TRIM( str) > 0)
      lstr = LEN_TRIM( str)
      str( 1:lstr-1) = str( 2:lstr)
      str( lstr:lstr) = ' '
    END DO

  END SUBROUTINE remove_leading_spaces

END MODULE configuration_module