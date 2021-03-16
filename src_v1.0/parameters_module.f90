MODULE parameters_module
  USE configuration_module, ONLY : dp, C
  IMPLICIT NONE

  ! This module contains often used parameterizations and often used physical parameters.

  ! Physical parameters:
  REAL(dp), PARAMETER :: CC                   = 8.7E-04_dp    ! Clausius Clapeyron gradient [K m^-1]
  REAL(dp), PARAMETER :: T0                   = 273.16_dp     ! Triple point of water [K]
  REAL(dp), SAVE      :: m_enh_sheet          = 4.0_dp        ! Flow enhancement parameter for the sheet [-]
  REAL(dp), PARAMETER :: m_enh_shelf          = 1.0_dp        ! Flow enhancement parameter for the shelf [-]
  REAL(dp), PARAMETER :: n_flow               = 3.0_dp        ! Parameter used in flow law [-]
  REAL(dp), PARAMETER :: m_flow               = 3.0_dp        ! Parameter used in Weertman type sliding law []
  REAL(dp), PARAMETER :: q_plastic            = 0.30_dp       ! Parameter used for basal stress (inverse of m_flow)
  REAL(dp), PARAMETER :: u_threshold          = 100._dp       ! scaling of tau_yield to get the correct unit (function of q_plastic)
  REAL(dp), SAVE      :: Fghf                 = 1.72E06_dp    ! Geothermal Heat flux [J m^-2 yr^-1] Sclater et al. (1980), see Huybrechts (4.44)
  REAL(dp), PARAMETER :: grav                 = 9.81_dp       ! Acceleration of gravity [m s^-2]
  REAL(dp), PARAMETER :: De                   = 1.0E+25_dp    ! Lithospheric flexural rigidity [kg m^2 s^-2]
  REAL(dp), PARAMETER :: tau                  = 3000.0_dp     ! Relaxation time for bedrock adjustment [yr]
  REAL(dp), PARAMETER :: SMT                  = 271.15_dp     ! Seawater temperature [K]
  REAL(dp), PARAMETER :: ESOA                 = 3.62E+14_dp   ! Earth Surface Ocean Area [m^2]
  REAL(dp), PARAMETER :: k0                   = 0.9728_dp     ! Scale factor for the stereographic projection [-], calculate_k0() is more precise
  REAL(dp), PARAMETER :: kr                   = 1.041E+08_dp  ! Conductivity of rock [J m^-1 K^-1 yr^-1] mean of Turcotte and Schubert (1982), see Huybrechts (4.44)
  REAL(dp), PARAMETER :: rock_heat_capacity   = 1000.0_dp     ! Heat capacity of rock [J kg^-1 K^-1]
  REAL(dp), PARAMETER :: A_sliding            = 1.8E-10_dp    ! Sliding coefficient inversely proportional to the bed roughness [m^8 yr^-1 N^-3]
  REAL(dp), PARAMETER :: ice_density          =  910.0_dp     ! Ice density [kg m^-3]
  REAL(dp), PARAMETER :: seawater_density     = 1028.0_dp     ! Seawater density [kg m^-3]
  REAL(dp), PARAMETER :: rock_density         = 3000.0_dp     ! Rock density [kg m^-3]
  REAL(dp), PARAMETER :: mantle_density       = 3300.0_dp     ! Mantle density [kg m^-3]
  REAL(dp), PARAMETER :: Atmos_lapserate_low  = 0.005102_dp   ! Athmospheric lapse rate [K m^-1] if z < 1500 m
  REAL(dp), PARAMETER :: Atmos_lapserate_high = 0.014285_dp   ! Atmospheric lapse rate  [K m^-1] if z > 1500 m
  REAL(dp), PARAMETER :: Atmos_lapserate_NH   = 0.008_dp      ! A constant atmospheric lapse rate  [K m^-1] - calculated with RACMO (MH: 0.0073743)
  REAL(dp), PARAMETER :: R_earth              = 6.371221E6_dp ! Earth Radius [m] ! also eismint
  REAL(dp), PARAMETER :: L_fusion             = 3.335E+5_dp   ! Latent heat of fusion [J kg-1]
  REAL(dp), PARAMETER :: epsilon_sq_0         = 1.0E-30_dp    ! to prevent a zero viscocity (C_uv)
  REAL(dp), PARAMETER :: delta_v              = 0.01_dp       ! to prevent a zero basal drag (beta_base), same as B&B 09 (PISM)
  REAL(dp), PARAMETER :: ocean_area           = 3.611E14_dp   ! World ocean area [m^2]

  REAL(dp), PARAMETER :: Docean               = 4000._dp      ! average depth of the ocean (meters)
  REAL(dp), PARAMETER :: gamma_dw             = -0.28_dp      ! slope of d18O to deep-water temperature (permille / degC)

  REAL(dp), PARAMETER :: hsliso_ant_PD        = -60.373_dp    ! present day ice volume used for d18O of Antarctica (mseq)
  REAL(dp), PARAMETER :: hsliso_grl_PD        =  -7.089_dp    ! present day ice volume used for d18O of Greenland (mseq)
  
CONTAINS

END MODULE parameters_module
