MODULE data_types_module
  ! Contains all the different types for storing data. Put all together in a separate module so that
  ! all subroutines can use all types without interdependency conflicts, and also to make the
  ! modules with the actual physics code more readable.
  ! If only Types could be collapsed in BBEdit...

#include <petsc/finclude/petscksp.h>

! ===== USE modules =====
! =======================

  USE mpi
  USE petscksp
  USE configuration_module,        ONLY: dp, C
  USE data_types_netcdf_module

  IMPLICIT NONE

! ===== Data types =====
! ======================

  ! == Ice model
  ! ============

  TYPE type_sparse_matrix_CSR_dp
    ! Compressed Sparse Row (CSR) format matrix

    INTEGER,                    POINTER     :: m,n                         ! A = [m-by-n]
    INTEGER,                    POINTER     :: nnz_max                     ! Maximum number of non-zero entries in A (determines how much memory is allocated)
    INTEGER,                    POINTER     :: nnz                         ! Number         of non-zero entries in A (determines how much memory is allocated)
    INTEGER,  DIMENSION(:    ), POINTER     :: ptr
    INTEGER,  DIMENSION(:    ), POINTER     :: index
    REAL(dp), DIMENSION(:    ), POINTER     :: val
    INTEGER :: wm, wn, wnnz_max, wnnz, wptr, windex, wval

  END TYPE type_sparse_matrix_CSR_dp

  TYPE type_ice_model
    ! The ice dynamics sub-model data structure.

    ! Basic data - ice thickness, bedrock & surface elevation, sea level (geoid elevation), englacial temperature, and ice velocities
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_a                        ! Ice thickness [m]
    REAL(dp), DIMENSION(:    ), POINTER     :: Hb_a                        ! Bedrock elevation [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:    ), POINTER     :: Hs_a                        ! Surface elevation [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:    ), POINTER     :: SL_a                        ! Sea level (geoid elevation) [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:    ), POINTER     :: TAF_a                       ! Thickness above flotation [m]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Ti_a                        ! Englacial temperature [K]
    INTEGER :: wHi_a, wHb_a, wHs_a, wSL_a, wTAF_a, wTi_a

    ! Ice velocities
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_3D_a                      ! 3-D ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_3D_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_3D_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_3D_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: w_3D_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: w_3D_b
    INTEGER :: wu_3D_a, wv_3D_a, wu_3D_b, wv_3D_b, ww_3D_a, ww_3D_b

    REAL(dp), DIMENSION(:    ), POINTER     :: u_vav_a                     ! Vertically averaged ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:    ), POINTER     :: v_vav_a
    REAL(dp), DIMENSION(:    ), POINTER     :: u_vav_b
    REAL(dp), DIMENSION(:    ), POINTER     :: v_vav_b
    REAL(dp), DIMENSION(:    ), POINTER     :: uabs_vav_a
    REAL(dp), DIMENSION(:    ), POINTER     :: uabs_vav_b
    INTEGER :: wu_vav_a, wv_vav_a, wu_vav_b, wv_vav_b, wuabs_vav_a, wuabs_vav_b

    REAL(dp), DIMENSION(:    ), POINTER     :: u_surf_a                    ! Ice velocity at the surface [m yr^-1]
    REAL(dp), DIMENSION(:    ), POINTER     :: v_surf_a
    REAL(dp), DIMENSION(:    ), POINTER     :: u_surf_b
    REAL(dp), DIMENSION(:    ), POINTER     :: v_surf_b
    REAL(dp), DIMENSION(:    ), POINTER     :: w_surf_a
    REAL(dp), DIMENSION(:    ), POINTER     :: uabs_surf_a
    REAL(dp), DIMENSION(:    ), POINTER     :: uabs_surf_b
    INTEGER :: wu_surf_a, wv_surf_a, wu_surf_b, wv_surf_b, ww_surf_a, wuabs_surf_a, wuabs_surf_b

    REAL(dp), DIMENSION(:    ), POINTER     :: u_base_a                    ! Ice velocity at the base [m yr^-1]
    REAL(dp), DIMENSION(:    ), POINTER     :: v_base_a
    REAL(dp), DIMENSION(:    ), POINTER     :: u_base_b
    REAL(dp), DIMENSION(:    ), POINTER     :: v_base_b
    REAL(dp), DIMENSION(:    ), POINTER     :: w_base_a
    REAL(dp), DIMENSION(:    ), POINTER     :: uabs_base_a
    REAL(dp), DIMENSION(:    ), POINTER     :: uabs_base_b
    INTEGER :: wu_base_a, wv_base_a, wu_base_b, wv_base_b, ww_base_a, wuabs_base_a, wuabs_base_b

    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_3D_SIA_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_3D_SIA_b
    REAL(dp), DIMENSION(:    ), POINTER     :: u_base_SSA_b
    REAL(dp), DIMENSION(:    ), POINTER     :: v_base_SSA_b
    INTEGER :: wu_3D_SIA_b, wv_3D_SIA_b, wu_base_SSA_b, wv_base_SSA_b

    ! Different masks
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_land_a
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_ocean_a
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_lake_a
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_ice_a
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_sheet_a
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_shelf_a
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_coast_a
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_margin_a
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_gl_a
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_cf_a
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_a
    REAL(dp), DIMENSION(:    ), POINTER     :: f_grnd_a
    REAL(dp), DIMENSION(:    ), POINTER     :: f_grnd_b
    INTEGER,  DIMENSION(:    ), POINTER     :: basin_ID                    ! The drainage basin to which each grid cell belongs
    INTEGER,                    POINTER     :: nbasins                     ! Total number of basins defined for this region
    INTEGER :: wmask_land_a, wmask_ocean_a, wmask_lake_a, wmask_ice_a, wmask_sheet_a, wmask_shelf_a
    INTEGER :: wmask_coast_a, wmask_margin_a, wmask_gl_a, wmask_cf_a, wmask_a, wf_grnd_a, wf_grnd_b
    INTEGER :: wbasin_ID, wnbasins

    ! Ice physical properties
    REAL(dp), DIMENSION(:,:  ), POINTER     :: A_flow_3D_a                 ! Flow parameter [Pa^-3 y^-1]
    REAL(dp), DIMENSION(:    ), POINTER     :: A_flow_vav_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Ti_pmp_a                    ! The pressure melting point temperature [K]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Cpi_a                       ! Specific heat capacity of ice [J kg^-1 K^-1].
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Ki_a                        ! Conductivity of ice [J m^-1 K^-1 yr^-1].
    INTEGER :: wA_flow_3D_a, wA_flow_vav_a, wTi_pmp_a, wCpi_a, wKi_a

    ! Zeta derivatives
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dzeta_dt_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dzeta_dx_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dzeta_dy_a
    REAL(dp), DIMENSION(:    ), POINTER     :: dzeta_dz_a
    INTEGER :: wdzeta_dt_a, wdzeta_dx_a, wdzeta_dy_a, wdzeta_dz_a

    ! Ice dynamics - basal hydrology
    REAL(dp), DIMENSION(:    ), POINTER     :: overburden_pressure_a       ! Overburden pressure ( = H * rho_i * g) [Pa]
    REAL(dp), DIMENSION(:    ), POINTER     :: pore_water_pressure_a       ! Pore water pressure (determined by basal hydrology model) [Pa]
    REAL(dp), DIMENSION(:    ), POINTER     :: Neff_a                      ! Effective pressure ( = overburden pressure - pore water pressure) [Pa]
    INTEGER :: woverburden_pressure_a, wpore_water_pressure_a, wNeff_a

    ! Ice dynamics - basal roughness / friction
    REAL(dp), DIMENSION(:    ), POINTER     :: phi_fric_a                  ! Till friction angle (degrees)
    REAL(dp), DIMENSION(:    ), POINTER     :: phi_fric_inv_a              ! Inverted till friction angle (degrees)
    REAL(dp), DIMENSION(:    ), POINTER     :: tauc_a                      ! Till yield stress tauc   (used when choice_sliding_law = "Coloumb" or "Coulomb_regularised")
    REAL(dp), DIMENSION(:    ), POINTER     :: alpha_sq_a                  ! Coulomb-law friction coefficient [unitless]         (used when choice_sliding_law =             "Tsai2015", or "Schoof2005")
    REAL(dp), DIMENSION(:    ), POINTER     :: beta_sq_a                   ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3] (used when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005")
    REAL(dp), DIMENSION(:    ), POINTER     :: beta_sq_inv_a               ! Inverted power-law friction coefficient [Pa m^−1/3 yr^1/3]
    INTEGER :: wphi_fric_a, wphi_fric_inv_a, wtauc_a, walpha_sq_a, wbeta_sq_a, wbeta_sq_inv_a

    ! Ice dynamics - basal inversion
    REAL(dp), DIMENSION(:    ), POINTER     :: BIV_uabs_surf_target  ! Target surface velocity for the basal inversion [m/yr]
    INTEGER :: wBIV_uabs_surf_target

    ! Ice dynamics - physical terms in the SSA/DIVA
    REAL(dp), DIMENSION(:    ), POINTER     :: taudx_b                     ! x-component of the driving stress
    REAL(dp), DIMENSION(:    ), POINTER     :: taudy_b                     ! x-component of the driving stress
    REAL(dp), DIMENSION(:    ), POINTER     :: du_dx_a                     ! Vertically averaged   xx strain rate
    REAL(dp), DIMENSION(:    ), POINTER     :: du_dy_a                     ! Vertically averaged   xy strain rate
    REAL(dp), DIMENSION(:    ), POINTER     :: dv_dx_a                     ! Vertically averaged   yy strain rate
    REAL(dp), DIMENSION(:    ), POINTER     :: dv_dy_a                     ! Vertically averaged   yy strain rate
    REAL(dp), DIMENSION(:,:  ), POINTER     :: du_dz_3D_b                  ! 3-D                   xz strain rate
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dv_dz_3D_b                  ! 3-D                   yz strain rate
    REAL(dp), DIMENSION(:,:  ), POINTER     :: visc_eff_3D_a               ! 3-D                   effective viscosity
    REAL(dp), DIMENSION(:    ), POINTER     :: visc_eff_int_a              ! Vertically integrated effective viscosity
    REAL(dp), DIMENSION(:    ), POINTER     :: N_a                         ! Product term N = eta * H
    REAL(dp), DIMENSION(:    ), POINTER     :: beta_a                      ! Sliding term beta (as in, [basal shear stress] = [beta] * [basal velocity])
    REAL(dp), DIMENSION(:    ), POINTER     :: beta_eff_a                  ! Beta_eff, appearing in the DIVA
    REAL(dp), DIMENSION(:    ), POINTER     :: beta_eff_b
    REAL(dp), DIMENSION(:    ), POINTER     :: taubx_b                     ! x-component of the basal shear stress
    REAL(dp), DIMENSION(:    ), POINTER     :: tauby_b                     ! y-component of the basal shear stress
    REAL(dp), DIMENSION(:    ), POINTER     :: F2_a                        ! F2, appearing in the DIVA
    REAL(dp), DIMENSION(:    ), POINTER     :: u_prev_b
    REAL(dp), DIMENSION(:    ), POINTER     :: v_prev_b
    INTEGER :: wtaudx_b, wtaudy_b
    INTEGER :: wdu_dx_a, wdu_dy_a, wdv_dx_a, wdv_dy_a, wdu_dz_3D_b, wdv_dz_3D_b, wvisc_eff_3D_a, wvisc_eff_int_a, wN_a
    INTEGER :: wbeta_a, wbeta_eff_a, wbeta_eff_b, wtaubx_b, wtauby_b, wF2_a
    INTEGER :: wu_prev_b, wv_prev_b

    ! Ice dynamics - some administrative stuff to make solving the SSA/DIVA more efficient
    INTEGER,  DIMENSION(:    ), POINTER     :: ti2n_u, ti2n_v
    INTEGER,  DIMENSION(:,:  ), POINTER     :: n2ti_uv
    TYPE(type_sparse_matrix_CSR_dp)         :: M_SSADIVA
    INTEGER                                 :: DIVA_SOR_nit
    REAL(dp)                                :: DIVA_SOR_tol
    REAL(dp)                                :: DIVA_SOR_omega
    REAL(dp)                                :: DIVA_PETSc_rtol
    REAL(dp)                                :: DIVA_PETSc_abstol
    INTEGER :: wti2n_u, wti2n_v, wn2ti_uv

    ! Ice dynamics - ice thickness calculation
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dVi_in
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dVi_out
    REAL(dp), DIMENSION(:    ), POINTER     :: dHi_dt_a
    REAL(dp), DIMENSION(:    ), POINTER     :: dHs_dt_a
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_tplusdt_a
    INTEGER :: wdVi_in, wdVi_out, wdHi_dt_a, wdHs_dt_a, wHi_tplusdt_a

    ! Ice dynamics - calving
    REAL(dp), DIMENSION(:    ), POINTER     :: float_margin_frac_a         ! Ice-covered fraction for calving front pixels
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_eff_cf_a                 ! Effective ice thickness at calving front pixels (= Hi of thinnest non-calving-front neighbour)
    INTEGER :: wfloat_margin_frac_a, wHi_eff_cf_a

    ! Ice dynamics - predictor/corrector ice thickness update
    REAL(dp),                   POINTER     :: pc_zeta
    REAL(dp), DIMENSION(:    ), POINTER     :: pc_tau
    REAL(dp), DIMENSION(:    ), POINTER     :: pc_fcb
    REAL(dp),                   POINTER     :: pc_eta
    REAL(dp),                   POINTER     :: pc_eta_prev
    REAL(dp),                   POINTER     :: pc_beta1
    REAL(dp),                   POINTER     :: pc_beta2
    REAL(dp),                   POINTER     :: pc_beta3
    REAL(dp),                   POINTER     :: pc_beta4
    REAL(dp), DIMENSION(:    ), POINTER     :: pc_f1
    REAL(dp), DIMENSION(:    ), POINTER     :: pc_f2
    REAL(dp), DIMENSION(:    ), POINTER     :: pc_f3
    REAL(dp), DIMENSION(:    ), POINTER     :: pc_f4
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_old
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_pred
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_corr
    INTEGER :: wpc_zeta, wpc_tau, wpc_fcb, wpc_eta, wpc_eta_prev, wpc_beta1, wpc_beta2, wpc_beta3, wpc_beta4
    INTEGER :: wpc_f1, wpc_f2, wpc_f3, wpc_f4, wHi_old, wHi_pred, wHi_corr

    ! Thermodynamics
    INTEGER,  DIMENSION(:    ), POINTER     :: mask_ice_a_prev        ! Ice mask from previous time step
    REAL(dp), DIMENSION(:,:  ), POINTER     :: internal_heating_a     ! Internal heating due to deformation
    REAL(dp), DIMENSION(:    ), POINTER     :: frictional_heating_a   ! Friction heating due to basal sliding
    REAL(dp), DIMENSION(:    ), POINTER     :: GHF_a                  ! Geothermal heat flux
    INTEGER :: wmask_ice_a_prev, winternal_heating_a, wfrictional_heating_a, wGHF_a

    ! Isotope content
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_a_prev
    REAL(dp), DIMENSION(:    ), POINTER     :: IsoRef
    REAL(dp), DIMENSION(:    ), POINTER     :: IsoSurf
    REAL(dp), DIMENSION(:    ), POINTER     :: MB_iso
    REAL(dp), DIMENSION(:    ), POINTER     :: IsoIce
    REAL(dp), DIMENSION(:    ), POINTER     :: IsoIce_new
    INTEGER :: wHi_a_prev, wIsoRef, wIsoSurf, wMB_iso, wIsoIce, wIsoIce_new

    ! ELRA GIA model
    INTEGER,                    POINTER     :: flex_prof_rad
    REAL(dp), DIMENSION(:,:  ), POINTER     :: flex_prof_grid
    REAL(dp), DIMENSION(:    ), POINTER     :: surface_load_PD_mesh
    REAL(dp), DIMENSION(:    ), POINTER     :: surface_load_mesh
    REAL(dp), DIMENSION(:    ), POINTER     :: surface_load_rel_mesh
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surface_load_rel_grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surface_load_rel_ext_grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_eq_grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_dt_grid
    REAL(dp), DIMENSION(:    ), POINTER     :: dHb_a
    REAL(dp), DIMENSION(:    ), POINTER     :: dHb_dt_a
    REAL(dp), DIMENSION(:    ), POINTER     :: dSL_dt_a
    INTEGER :: wflex_prof_rad, wflex_prof_grid
    INTEGER :: wsurface_load_PD_mesh, wsurface_load_mesh, wsurface_load_rel_mesh, wsurface_load_rel_grid, wsurface_load_rel_ext_grid
    INTEGER :: wdHb_eq_grid, wdHb_grid, wdHb_dt_grid
    INTEGER :: wdHb_a, wdHb_dt_a, wdSL_dt_a

    ! Mesh adaptation data
    REAL(dp), DIMENSION(:    ), POINTER     :: surf_curv
    REAL(dp), DIMENSION(:    ), POINTER     :: log_velocity
    INTEGER :: wsurf_curv, wlog_velocity

    ! Useful extra stuff
    REAL(dp), DIMENSION(:    ), POINTER     :: dHi_a   ! Ice thickness difference w.r.t. PD
    REAL(dp), DIMENSION(:    ), POINTER     :: dHs_a   ! Ice elevation difference w.r.t. PD
    INTEGER :: wdHi_a, wdHs_a

  END TYPE type_ice_model

  ! == Mesh
  ! =======

  TYPE type_mesh
    ! The unstructured triangular mesh.

    ! Basic meta properties
    ! =====================

    CHARACTER(LEN=3)                        :: region_name                   ! NAM, EAS, GRL, ANT
    REAL(dp),                   POINTER     :: lambda_M                      ! Oblique stereographic projection parameters
    REAL(dp),                   POINTER     :: phi_M
    REAL(dp),                   POINTER     :: beta_stereo
    REAL(dp),                   POINTER     :: xmin                          ! X and Y range of the square covered by the mesh
    REAL(dp),                   POINTER     :: xmax
    REAL(dp),                   POINTER     :: ymin
    REAL(dp),                   POINTER     :: ymax
    REAL(dp),                   POINTER     :: tol_dist                      ! Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    INTEGER,                    POINTER     :: nV_mem                        ! Size of allocated memory for vertices
    INTEGER,                    POINTER     :: nTri_mem                      ! Size of allocated memory for triangles
    INTEGER,                    POINTER     :: nC_mem                        ! Maximum allowed number of connections per vertex
    INTEGER,                    POINTER     :: nV                            ! Number of vertices
    INTEGER,                    POINTER     :: nTri                          ! Number of triangles
    INTEGER,                    POINTER     :: perturb_dir                   ! Perturbation direction (0 = anticlockwise, 1 = clockwise)
    REAL(dp),                   POINTER     :: alpha_min                     ! Sharpest inner angle allowed by Rupperts algorithm
    REAL(dp),                   POINTER     :: dz_max_ice                    ! Maximum allowed vertical error in Rupperts algorithm over ice
    REAL(dp),                   POINTER     :: res_max                       ! Maximum resolution anywhere
    REAL(dp),                   POINTER     :: res_max_margin                ! Maximum resolution over the ice margin
    REAL(dp),                   POINTER     :: res_max_gl                    ! Maximum resolution over the grounding line
    REAL(dp),                   POINTER     :: res_max_cf                    ! Maximum resolution over the calving front
    REAL(dp),                   POINTER     :: res_max_mountain              ! Maximum resolution over ice-free mountains  (important for getting the inception right)
    REAL(dp),                   POINTER     :: res_max_coast                 ! Maximum resolution over ice-free coast line (to make plots look nicer)
    REAL(dp),                   POINTER     :: res_min                       ! Minimum resolution anywhere ( = MINVAL([res_max_margin, res_max_gl, res_max_cf]))
    REAL(dp),                   POINTER     :: resolution_min                ! Finest   resolution of the mesh ( = MINVAL(R), where R = distance to nearest neighbour)
    REAL(dp),                   POINTER     :: resolution_max                ! Coarsest resolution of the mesh
    INTEGER :: wlambda_M, wphi_M, wbeta_stereo, wxmin, wxmax, wymin, wymax, wtol_dist, wnV_mem, wnTri_mem, wnC_mem, wnV, wnTri, wperturb_dir
    INTEGER :: walpha_min, wdz_max_ice, wres_max, wres_max_margin, wres_max_gl, wres_max_cf, wres_max_mountain, wres_max_coast, wres_min, wresolution_min, wresolution_max

    ! Primary mesh data (needed for mesh creation & refinement)
    ! =========================================================

    ! Vertex data
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V                             ! The X and Y coordinates of all the vertices
    INTEGER,  DIMENSION(:    ), POINTER     :: nC                            ! The number of other vertices this vertex is connected to
    INTEGER,  DIMENSION(:,:  ), POINTER     :: C                             ! The list   of other vertices this vertex is connected to (ordered counter-clockwise, from edge to edge for edge vertices)
    INTEGER,  DIMENSION(:    ), POINTER     :: niTri                         ! The number of triangles this vertex is a part of
    INTEGER,  DIMENSION(:,:  ), POINTER     :: iTri                          ! The list   of triangles this vertex is a part of (ordered counter-clockwise)
    INTEGER,  DIMENSION(:    ), POINTER     :: edge_index                    ! Each vertex's Edge Index; 0 = free, 1 = north, 2 = northeast, etc.
    INTEGER,  DIMENSION(:    ), POINTER     :: mesh_old_ti_in                ! When creating a new mesh: for every new mesh vertex, the old mesh triangle that contains it
    INTEGER :: wV, wnC, wC, wniTri, wiTri, wedge_index, wmesh_old_ti_in

    ! Triangle data
    INTEGER,  DIMENSION(:,:  ), POINTER     :: Tri                           ! The triangle array: Tri(ti) = [vi1, vi2, vi3] (vertices ordered counter-clockwise)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: TriCC                         ! The X,Y-coordinates of each triangle's circumcenter
    INTEGER,  DIMENSION(:,:  ), POINTER     :: TriC                          ! The (up to) three neighbour triangles (order across from 1st, 2nd and 3d vertex, respectively)
    INTEGER :: wTri, wTriCC, wTriC

    ! Refinement lists
    INTEGER,  DIMENSION(:,:  ), POINTER     :: Triflip                       ! List of triangles to flip, used in updating Delaunay triangulation
    INTEGER,  DIMENSION(:    ), POINTER     :: RefMap                        ! Map   of which triangles      have been marked for refining
    INTEGER,  DIMENSION(:    ), POINTER     :: RefStack                      ! Stack of       triangles that have been marked for refining
    INTEGER,                    POINTER     :: RefStackN
    INTEGER :: wTriflip, wRefMap, wRefStack, wRefStackN

    ! Maps+stacks for FloodFill-style searching
    INTEGER,  DIMENSION(:    ), POINTER     :: VMap, VStack1, VStack2
    INTEGER,  DIMENSION(:    ), POINTER     :: TriMap, TriStack1, TriStack2
    INTEGER :: wVMap, wVStack1, wVStack2, wTriMap, wTriStack1, wTriStack2
    INTEGER :: VStackN1, VStackN2
    INTEGER :: TriStackN1, TriStackN2

    ! Points-of-Interest
    INTEGER,                    POINTER     :: nPOI                          ! Number of Points of Interest (POI) in this mesh
    REAL(dp), DIMENSION(:,:  ), POINTER     :: POI_coordinates               ! Lat-lon coordinates of a POI
    REAL(dp), DIMENSION(:,:  ), POINTER     :: POI_XY_coordinates            ! X-Y     coordinates of a POI
    REAL(dp), DIMENSION(:    ), POINTER     :: POI_resolutions               ! Resolution          of a POI (km)
    INTEGER,  DIMENSION(:,:  ), POINTER     :: POI_vi                        ! The three vertices surrounding a POI
    REAL(dp), DIMENSION(:,:  ), POINTER     :: POI_w                         ! Their relative weights in trilinear interpolation
    INTEGER :: wnPOI, wPOI_coordinates, wPOI_XY_coordinates, wPOI_resolutions, wPOI_vi, wPOI_w

    ! Secondary mesh data
    ! ===================

    ! Derived geometry data
    REAL(dp), DIMENSION(:    ), POINTER     :: A                             ! The area             of each vertex's Voronoi cell
    REAL(dp), DIMENSION(:,:  ), POINTER     :: VorGC                         ! The geometric centre of each vertex's Voronoi cell
    REAL(dp), DIMENSION(:    ), POINTER     :: R                             ! The resolution (defined as distance to nearest neighbour)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Cw                            ! The width of all vertex connections (= length of the shared Voronoi cell edge)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: TriGC                         ! The X,Y-coordinates of each triangle's geometric centre
    INTEGER,  DIMENSION(:    ), POINTER     :: Tri_edge_index                ! Same as for vertices, only NE, SE, etc. aren't used
    REAL(dp), DIMENSION(:    ), POINTER     :: TriA                          ! The area of each triangle
    INTEGER :: wA, wVorGC, wR, wCw, wTriGC, wTri_edge_index, wTriA

    ! Staggered (Arakawa C) mesh
    INTEGER,                    POINTER     :: nAc                           ! The number of Ac vertices
    REAL(dp), DIMENSION(:,:  ), POINTER     :: VAc                           ! x,y coordinates of the Ac vertices
    INTEGER,  DIMENSION(:,:  ), POINTER     :: Aci                           ! Mapping array from the Aa to the Ac mesh
    INTEGER,  DIMENSION(:,:  ), POINTER     :: iAci                          ! Mapping array from the Ac to the Aa mesh
    INTEGER,  DIMENSION(:    ), POINTER     :: edge_index_Ac
    INTEGER :: wnAc, wVAc, wAci, wiAci, wedge_index_Ac

    ! Matrix operators: mapping
    TYPE(tMat)                              :: M_map_a_b                     ! Operation: map     from the a-grid to the b-grid
    TYPE(tMat)                              :: M_map_a_c                     ! Operation: map     from the a-grid to the c-grid
    TYPE(tMat)                              :: M_map_b_a                     ! Operation: map     from the b-grid to the a-grid
    TYPE(tMat)                              :: M_map_b_c                     ! Operation: map     from the b-grid to the c-grid
    TYPE(tMat)                              :: M_map_c_a                     ! Operation: map     from the c-grid to the a-grid
    TYPE(tMat)                              :: M_map_c_b                     ! Operation: map     from the c-grid to the b-grid

    ! Matrix operators: d/dx
    TYPE(tMat)                              :: M_ddx_a_a                     ! Operation: d/dx    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddx_a_b                     ! Operation: d/dx    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddx_a_c                     ! Operation: d/dx    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddx_b_a                     ! Operation: d/dx    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddx_b_b                     ! Operation: d/dx    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddx_b_c                     ! Operation: d/dx    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddx_c_a                     ! Operation: d/dx    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddx_c_b                     ! Operation: d/dx    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddx_c_c                     ! Operation: d/dx    from the a-grid to the a-grid

    ! Matrix operators: d/dy
    TYPE(tMat)                              :: M_ddy_a_a                     ! Operation: d/dy    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddy_a_b                     ! Operation: d/dy    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddy_a_c                     ! Operation: d/dy    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddy_b_a                     ! Operation: d/dy    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddy_b_b                     ! Operation: d/dy    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddy_b_c                     ! Operation: d/dy    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddy_c_a                     ! Operation: d/dy    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddy_c_b                     ! Operation: d/dy    from the a-grid to the a-grid
    TYPE(tMat)                              :: M_ddy_c_c                     ! Operation: d/dy    from the a-grid to the a-grid

    ! 2nd-order accurate matrix operators on the b-grid
    TYPE(tMat)                              :: M2_ddx_b_b                    ! Operation: d/dx    from the b-grid to the b-grid
    TYPE(tMat)                              :: M2_ddy_b_b                    ! Operation: d/dy    from the b-grid to the b-grid
    TYPE(tMat)                              :: M2_d2dx2_b_b                  ! Operation: d2/dx2  from the b-grid to the b-grid
    TYPE(tMat)                              :: M2_d2dxdy_b_b                 ! Operation: d2/dxdy from the b-grid to the b-grid
    TYPE(tMat)                              :: M2_d2dy2_b_b                  ! Operation: d2/dy2  from the b-grid to the b-grid

    TYPE(type_sparse_matrix_CSR_dp)         :: M2_ddx_b_b_CSR                ! Operation: d/dx    from the b-grid to the b-grid
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_ddy_b_b_CSR                ! Operation: d/dy    from the b-grid to the b-grid
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_d2dx2_b_b_CSR              ! Operation: d2/dx2  from the b-grid to the b-grid
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_d2dxdy_b_b_CSR             ! Operation: d2/dxdy from the b-grid to the b-grid
    TYPE(type_sparse_matrix_CSR_dp)         :: M2_d2dy2_b_b_CSR              ! Operation: d2/dy2  from the b-grid to the b-grid

    ! Matrix operator for applying Neumann boundary conditions to triangles at the domain border
    TYPE(tMat)                              :: M_Neumann_BC_b
    TYPE(type_sparse_matrix_CSR_dp)         :: M_Neumann_BC_b_CSR

    ! Lat/lon coordinates
    REAL(dp), DIMENSION(:    ), POINTER     :: lat
    REAL(dp), DIMENSION(:    ), POINTER     :: lon
    INTEGER :: wlat, wlon

    ! Transect
    INTEGER,                    POINTER     :: nV_transect                   ! Number of vertex pairs for the transect
    INTEGER,  DIMENSION(:,:  ), POINTER     :: vi_transect                   ! List   of vertex pairs for the transect
    REAL(dp), DIMENSION(:,:  ), POINTER     :: w_transect                    ! Interpolation weights for the vertex pairs
    INTEGER :: wnV_transect, wvi_transect, ww_transect

    ! Parallelisation
    INTEGER                                 :: vi1, vi2                      ! Vertices
    INTEGER                                 :: ti1, ti2                      ! Triangles
    INTEGER                                 :: ci1, ci2                      ! Edges

  END TYPE type_mesh

  ! == Regular x/y-grid
  ! ===================

  TYPE type_grid
    ! A regular square grid covering a model region

    ! Basic grid data
    INTEGER,                    POINTER     :: nx, ny, n
    REAL(dp),                   POINTER     :: dx
    REAL(dp), DIMENSION(:    ), POINTER     :: x, y
    REAL(dp),                   POINTER     :: xmin, xmax, ymin, ymax
    INTEGER :: wnx, wny, wn, wdx, wx, wy, wxmin, wxmax, wymin, wymax

    ! Parallelisation by domain decomposition
    INTEGER                                 :: i1, i2, j1, j2

    ! Sparse matrices representing the remapping operations between a mesh and a grid
    TYPE(tMat)                              :: M_map_grid2mesh              ! Remapping from a grid to a mesh using second-order conservative remapping
    TYPE(tMat)                              :: M_map_mesh2grid              ! Remapping from a mesh to a grid using second-order conservative remapping
    REAL(dp),                   POINTER     :: tol_dist
    INTEGER :: wtol_dist

    ! Conversion tables for grid-form vs. vector-form data
    INTEGER,  DIMENSION(:,:  ), POINTER     :: ij2n, n2ij
    INTEGER :: wij2n, wn2ij

    ! Lat-lon coordinates
    REAL(dp), DIMENSION(:,:  ), POINTER     :: lat, lon
    INTEGER :: wlat, wlon

    ! Projection parameters for the grid
    REAL(dp),                   POINTER     :: lambda_M
    REAL(dp),                   POINTER     :: phi_M
    REAL(dp),                   POINTER     :: beta_stereo
    INTEGER :: wlambda_M, wphi_M, wbeta_stereo

  END TYPE type_grid

  ! == Regular lon/lat-grid
  ! =======================

  TYPE type_grid_lonlat
    ! A regular square lon/lat-grid

    ! Basic grid data
    INTEGER,                    POINTER     :: nlon, nlat
    REAL(dp),                   POINTER     :: dlon, dlat
    REAL(dp), DIMENSION(:    ), POINTER     :: lon, lat
    REAL(dp),                   POINTER     :: lonmin, lonmax, latmin, latmax
    INTEGER :: wnlon, wnlat, wn, wdlon, wdlat, wlon, wlat, wlonmin, wlonmax, wlatmin, wlatmax

    ! Parallelisation by domain decomposition
    INTEGER                                 :: i1, i2, j1, j2

    ! Projection parameters for the grid
    REAL(dp),                   POINTER     :: lambda_M
    REAL(dp),                   POINTER     :: phi_M
    REAL(dp),                   POINTER     :: beta_stereo
    INTEGER :: wlambda_M, wphi_M, wbeta_stereo

  END TYPE type_grid_lonlat

  ! == Reference geometries
  ! =======================

  TYPE type_reference_geometry
    ! Data structure containing a reference ice-sheet geometry (either schematic or read from an external file).

    ! NetCDF file containing the data
    TYPE(type_netcdf_reference_geometry)    :: netcdf

    ! Raw data as read from a NetCDF file
    TYPE(type_grid)                         :: grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_grid
    INTEGER :: wHi_grid, wHb_grid, wHs_grid

    ! Derived data on the grid (surface curvature and masks, needed for mesh creation)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surf_curv
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_land
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ocean
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_sheet
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_shelf
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_margin
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_gl
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_cf
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_coast
    INTEGER :: wsurf_curv, wmask_land, wmask_ocean, wmask_ice, wmask_sheet, wmask_shelf, wmask_margin, wmask_gl, wmask_cf, wmask_coast

    ! Data on the model mesh
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi
    REAL(dp), DIMENSION(:    ), POINTER     :: Hb
    REAL(dp), DIMENSION(:    ), POINTER     :: Hs
    INTEGER :: wHi, wHb, wHs

  END TYPE type_reference_geometry

  ! == Mapping
  ! ==========

  TYPE type_single_row_mapping_matrices
    ! Results from integrating a single set of lines

    INTEGER                                 :: n_max
    INTEGER                                 :: n
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: index_left
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: LI_xdy, LI_mxydx, LI_xydy

  END TYPE type_single_row_mapping_matrices

  TYPE type_remapping_mesh_mesh
    ! Sparse matrices representing the remapping operations between two meshes for different remapping methods

    INTEGER                                 :: int_dummy
    TYPE(tMat)                              :: M_trilin                     ! Remapping using trilinear interpolation
    TYPE(tMat)                              :: M_nearest_neighbour          ! Remapping using nearest-neighbour interpolation
    TYPE(tMat)                              :: M_cons_1st_order             ! Remapping using first-order conservative remapping
    TYPE(tMat)                              :: M_cons_2nd_order             ! Remapping using second-order conservative remapping

  END TYPE type_remapping_mesh_mesh

  TYPE type_remapping_lonlat2mesh
    ! Indices and weights for mapping data from a global lon/lat-grid to the model mesh using bilinear interpolation

    INTEGER,  DIMENSION(:    ), POINTER     :: ilat1, ilat2, ilon1, ilon2
    REAL(dp), DIMENSION(:    ), POINTER     :: wlat1, wlat2, wlon1, wlon2
    INTEGER :: wilat1, wilat2, wilon1, wilon2, wwlat1, wwlat2, wwlon1, wwlon2

  END TYPE type_remapping_lonlat2mesh

  ! == Forcing
  ! ==========

  TYPE type_forcing_data
    ! Data structure containing model forcing data - CO2 record, d18O record, (global) insolation record

    ! Data for the inverse routine
    REAL(dp),                   POINTER     :: dT_glob                                   ! Modelled global mean annual surface temperature anomaly w.r.t. PD
    REAL(dp), DIMENSION(:    ), POINTER     :: dT_glob_history                           ! Time window (length set in config) listing previous dT_glob values
    INTEGER,                    POINTER     :: ndT_glob_history                          ! Number of entries (= length of time window / dt_coupling)
    REAL(dp),                   POINTER     :: dT_deepwater                              ! Modelled deep-water temperature anomaly (= window averaged dT_glob * scaling factor), scaling factor set in config
    INTEGER :: wdT_glob, wdT_glob_history, wndT_glob_history, wdT_deepwater

    REAL(dp),                   POINTER     :: d18O_NAM, d18O_EAS, d18O_GRL, d18O_ANT    ! Modelled benthic d18O contributions from ice volume in the model regions
    REAL(dp),                   POINTER     :: d18O_from_ice_volume_mod                  ! Modelled benthic d18O contribution from global ice volume
    REAL(dp),                   POINTER     :: d18O_from_temperature_mod                 ! Modelled benthic d18O contribution from deep-water temperature change
    REAL(dp),                   POINTER     :: d18O_mod                                  ! Modelled benthic d18O
    INTEGER :: wd18O_NAM, wd18O_EAS, wd18O_GRL, wd18O_ANT, wd18O_from_ice_volume_mod, wd18O_from_temperature_mod, wd18O_mod

    REAL(dp),                   POINTER     :: dT_glob_inverse                           ! Global mean annual surface temperature anomaly resulting from the inverse method
    REAL(dp), DIMENSION(:    ), POINTER     :: dT_glob_inverse_history                   ! Time window (length set in config) listing previous dT_glob values
    INTEGER,                    POINTER     :: ndT_glob_inverse_history                  ! Number of entries (= length of time window / dt_coupling)
    REAL(dp),                   POINTER     :: CO2_inverse                               ! CO2 resulting from the inverse method
    REAL(dp), DIMENSION(:    ), POINTER     :: CO2_inverse_history                       ! Time window (length set in config) listing previous CO2_inverse values
    INTEGER,                    POINTER     :: nCO2_inverse_history                      ! Number of entries (= length of time window / dt_coupling)
    REAL(dp),                   POINTER     :: CO2_mod                                   ! Either equal to CO2_obs or to CO2_inverse, for easier writing to output.
    INTEGER :: wdT_glob_inverse, wdT_glob_inverse_history, wndT_glob_inverse_history, wCO2_inverse, wCO2_inverse_history, wnCO2_inverse_history, wCO2_mod

    ! External forcing: CO2 record
    REAL(dp), DIMENSION(:    ), POINTER     :: CO2_time
    REAL(dp), DIMENSION(:    ), POINTER     :: CO2_record
    REAL(dp),                   POINTER     :: CO2_obs
    INTEGER :: wCO2_time, wCO2_record, wCO2_obs

    ! External forcing: d18O record
    REAL(dp), DIMENSION(:    ), POINTER     :: d18O_time
    REAL(dp), DIMENSION(:    ), POINTER     :: d18O_record
    REAL(dp),                   POINTER     :: d18O_obs
    REAL(dp),                   POINTER     :: d18O_obs_PD
    INTEGER :: wd18O_time, wd18O_record, wd18O_obs, wd18O_obs_PD

    ! External forcing: insolation
    TYPE(type_netcdf_insolation)            :: netcdf_ins
    INTEGER,                    POINTER     :: ins_nyears
    INTEGER,                    POINTER     :: ins_nlat
    REAL(dp), DIMENSION(:    ), POINTER     :: ins_time
    REAL(dp), DIMENSION(:    ), POINTER     :: ins_lat
    REAL(dp),                   POINTER     :: ins_t0, ins_t1
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ins_Q_TOA0, ins_Q_TOA1
    INTEGER :: wins_nyears, wins_nlat, wins_time, wins_lat, wins_t0, wins_t1, wins_Q_TOA0, wins_Q_TOA1

    ! External forcing: geothermal heat flux
    TYPE(type_netcdf_geothermal_heat_flux)  :: netcdf_ghf
    TYPE(type_grid_lonlat)                  :: grid_ghf
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ghf_ghf
    INTEGER :: wghf_ghf

    ! External forcing: sea level record
    REAL(dp), DIMENSION(:    ), POINTER     :: sealevel_time
    REAL(dp), DIMENSION(:    ), POINTER     :: sealevel_record
    REAL(dp),                   POINTER     :: sealevel_obs
    INTEGER :: wsealevel_time, wsealevel_record, wsealevel_obs

  END TYPE type_forcing_data

  ! == Global climate
  !==================

  TYPE type_climate_snapshot_global
    ! Global climate snapshot, either from present-day observations (e.g. ERA40) or from a GCM snapshot.

    CHARACTER(LEN=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.

    ! NetCDF file containing the data
    TYPE(type_netcdf_climate_data)          :: netcdf

    ! Grid
    INTEGER,                    POINTER     :: nlat, nlon
    REAL(dp), DIMENSION(:    ), POINTER     :: lat
    REAL(dp), DIMENSION(:    ), POINTER     :: lon
    INTEGER :: wnlat, wnlon, wlat, wlon

    ! General forcing info (not relevant for PD observations)
    REAL(dp),                   POINTER     :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp),                   POINTER     :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   POINTER     :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp),                   POINTER     :: orbit_obl
    REAL(dp),                   POINTER     :: orbit_pre
    INTEGER :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre

    ! Actual GCM data
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs                            ! Orography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T2m                           ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Precip                        ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_WE                       ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_SN                       ! Monthly mean south_north wind speed (m/s)
    INTEGER :: wHs, wT2m, wPrecip, wWind_WE, wWind_SN

    ! Paralelisation
    INTEGER                                 :: i1, i2                        ! Grid domain (:,i1:i2) of each process

  END TYPE type_climate_snapshot_global

  TYPE type_direct_climate_forcing_global
    ! Direct global climate forcing from a provided NetCDF file

    ! NetCDF file containing the data
    TYPE(type_netcdf_direct_climate_forcing_global) :: netcdf

    ! Grid
    INTEGER,                    POINTER     :: nlat, nlon
    REAL(dp), DIMENSION(:    ), POINTER     :: lat
    REAL(dp), DIMENSION(:    ), POINTER     :: lon
    INTEGER :: wnlat, wnlon, wlat, wlon

    ! Time dimension
    INTEGER,                    POINTER     :: nyears
    REAL(dp), DIMENSION(:    ), POINTER     :: time
    REAL(dp),                   POINTER     :: t0, t1
    INTEGER :: wnyears, wtime, wt0, wt1

    ! Actual data (two timeframes enveloping the model time)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs0,      Hs1                 ! Orography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T2m0,     T2m1                ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Precip0,  Precip1             ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_WE0, Wind_WE1            ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_SN0, Wind_SN1            ! Monthly mean south_north wind speed (m/s)
    INTEGER :: wHs0, wT2m0, wPrecip0, wWind_WE0, wWind_SN0
    INTEGER :: wHs1, wT2m1, wPrecip1, wWind_WE1, wWind_SN1

  END TYPE type_direct_climate_forcing_global

  TYPE type_direct_SMB_forcing_global
    ! Direct global SMB forcing from a provided NetCDF file

    ! NetCDF file containing the data
    TYPE(type_netcdf_direct_SMB_forcing_global) :: netcdf

    ! Grid
    INTEGER,                    POINTER     :: nlat, nlon
    REAL(dp), DIMENSION(:    ), POINTER     :: lat
    REAL(dp), DIMENSION(:    ), POINTER     :: lon
    INTEGER :: wnlat, wnlon, wlat, wlon

    ! Time dimension
    INTEGER,                    POINTER     :: nyears
    REAL(dp), DIMENSION(:    ), POINTER     :: time
    REAL(dp),                   POINTER     :: t0, t1
    INTEGER :: wnyears, wtime, wt0, wt1

    ! Actual data (two timeframes enveloping the model time)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: T2m_year0, T2m_year1          ! Monthly mean 2m air temperature [K] (needed for thermodynamics)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SMB_year0, SMB_year1          ! Yearly total SMB (m.i.e.)
    INTEGER :: wT2m_year0, wT2m_year1, wSMB_year0, wSMB_year1

  END TYPE type_direct_SMB_forcing_global

  TYPE type_climate_matrix_global
    ! The climate matrix data structure. Contains all the different global GCM snapshots.

    ! The present-day observed climate (e.g. ERA40)
    TYPE(type_climate_snapshot_global)       :: PD_obs

    ! The GCM climate snapshots
    TYPE(type_climate_snapshot_global)       :: GCM_PI
    TYPE(type_climate_snapshot_global)       :: GCM_warm
    TYPE(type_climate_snapshot_global)       :: GCM_cold

    ! ! Direct climate forcing
    TYPE(type_direct_climate_forcing_global) :: direct
    TYPE(type_direct_SMB_forcing_global)     :: SMB_direct

  END TYPE type_climate_matrix_global

    TYPE type_ocean_snapshot_global
    ! Global ocean snapshot, either from present-day observations (e.g. WOA18) or from a GCM ocean snapshot.

    CHARACTER(LEN=256)                      :: name                          ! 'WOA', 'COSMOS_LGM', etc.

    ! NetCDF file containing the data
    TYPE(type_netcdf_ocean_data)            :: netcdf

    ! Grid
    INTEGER,                    POINTER     :: nlat, nlon
    REAL(dp), DIMENSION(:    ), POINTER     :: lat
    REAL(dp), DIMENSION(:    ), POINTER     :: lon
    INTEGER :: wnlat, wnlon, wlat, wlon

    ! General forcing info (not relevant for PD observations)
    REAL(dp),                   POINTER     :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp),                   POINTER     :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   POINTER     :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp),                   POINTER     :: orbit_obl
    REAL(dp),                   POINTER     :: orbit_pre
    INTEGER :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre

    ! Global 3-D ocean data on the original vertical grid
    REAL(dp),                   POINTER     :: T_ocean_mean                  ! Regional mean ocean temperature (used for basal melt when no ocean temperature data is provided)
    REAL(dp), DIMENSION(:    ), POINTER     :: z_ocean_raw                   ! Vertical coordinate of the 3-D ocean fields [m below sea surface]
    INTEGER,                    POINTER     :: nz_ocean_raw                  ! Number of vertical layers in the 3-D ocean fields
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T_ocean_raw                   ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: S_ocean_raw                   ! 3-D annual mean ocean salinity    [PSU]
    INTEGER :: wT_ocean_mean, wz_ocean_raw, wnz_ocean_raw, wT_ocean_raw, wS_ocean_raw

    ! Global 3-D ocean data on the ice-model vertical grid
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T_ocean                       ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: S_ocean                       ! 3-D annual mean ocean salinity    [PSU]
    INTEGER :: wT_ocean, wS_ocean

    ! Paralelisation
    INTEGER                                 :: i1, i2                        ! Grid domain (:,i1:i2) of each process

  END TYPE type_ocean_snapshot_global

  TYPE type_ocean_matrix_global
    ! The ocean matrix data structure. Contains the PD observations and all the different global GCM snapshots.

    ! The present-day observed ocean (e.g. WOA18)
    TYPE(type_ocean_snapshot_global)        :: PD_obs

    ! The GCM ocean snapshots for climate and ocean.
    TYPE(type_ocean_snapshot_global)        :: GCM_PI
    TYPE(type_ocean_snapshot_global)        :: GCM_warm
    TYPE(type_ocean_snapshot_global)        :: GCM_cold

  END TYPE type_ocean_matrix_global

  ! == Regional climate
  !====================

  TYPE type_climate_snapshot_regional
    ! Global climate snapshot, either from present-day observations (e.g. ERA40) or from a GCM snapshot, projected onto the regional model mesh.

    CHARACTER(LEN=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.

    ! General forcing info (not relevant for PD observations)
    REAL(dp),                   POINTER     :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp),                   POINTER     :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   POINTER     :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp),                   POINTER     :: orbit_obl
    REAL(dp),                   POINTER     :: orbit_pre
    REAL(dp),                   POINTER     :: sealevel
    INTEGER :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre, wsealevel

    ! Actual observed / GCM data (read from external NetCDF file)
    REAL(dp), DIMENSION(:    ), POINTER     :: Hs                            ! Orography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: T2m                           ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Precip                        ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Wind_WE                       ! Monthly mean west-east   wind speed (m/s)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Wind_SN                       ! Monthly mean south-north wind speed (m/s)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Wind_LR                       ! Monthly mean wind speed in the x-direction (m/s)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Wind_DU                       ! Monthly mean wind speed in the y-direction (m/s)
    INTEGER :: wHs, wT2m, wPrecip, wHs_ref, wWind_WE, wWind_SN, wWind_LR, wWind_DU ! MPI windows to all these memory spaces

    ! Spatially variable lapse rate for GCM snapshots (see Berends et al., 2018)
    REAL(dp), DIMENSION(:    ), POINTER     :: lambda
    INTEGER :: wlambda

    ! Bias-corrected GCM data
    REAL(dp), DIMENSION(:,:  ), POINTER     :: T2m_corr                      ! Bias-corrected monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Precip_corr                   ! Bias-corrected monthly mean precipitation (m)
    INTEGER :: wT2m_corr, wPrecip_corr

    ! Reference absorbed insolation (for GCM snapshots), or insolation at model time for the applied climate
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Q_TOA                         ! Monthly mean insolation at the top of the atmosphere (W/m2) (taken from the prescribed insolation solution at orbit_time)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Albedo                        ! Monthly mean surface albedo (calculated using our own SMB scheme for consistency)
    REAL(dp), DIMENSION(:    ), POINTER     :: I_abs                         ! Total yearly absorbed insolation, used in the climate matrix for interpolation
    INTEGER :: wQ_TOA, wAlbedo, wI_abs

  END TYPE type_climate_snapshot_regional

  TYPE type_direct_climate_forcing_regional
    ! Direct regional climate forcing from a provided NetCDF file

    ! NetCDF file containing the data
    TYPE(type_netcdf_direct_climate_forcing_regional) :: netcdf

    ! Grid
    TYPE(type_grid)                         :: grid

    ! Time dimension
    INTEGER,                    POINTER     :: nyears
    REAL(dp), DIMENSION(:    ), POINTER     :: time
    REAL(dp),                   POINTER     :: t0, t1
    INTEGER :: wnyears, wtime, wt0, wt1

    ! Actual data (two timeframes enveloping the model time) on the source grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs0_raw,      Hs1_raw         ! Orography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T2m0_raw,     T2m1_raw        ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Precip0_raw,  Precip1_raw     ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_WE0_raw, Wind_WE1_raw    ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_SN0_raw, Wind_SN1_raw    ! Monthly mean south_north wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_LR0_raw, Wind_LR1_raw    ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_DU0_raw, Wind_DU1_raw    ! Monthly mean south_north wind speed (m/s)
    INTEGER :: wHs0_raw, wT2m0_raw, wPrecip0_raw, wWind_WE0_raw, wWind_SN0_raw, wWind_LR0_raw, wWind_DU0_raw
    INTEGER :: wHs1_raw, wT2m1_raw, wPrecip1_raw, wWind_WE1_raw, wWind_SN1_raw, wWind_LR1_raw, wWind_DU1_raw

    ! Actual data (two timeframes enveloping the model time) on the model mesh
    REAL(dp), DIMENSION(:    ), POINTER     :: Hs0,      Hs1                 ! Orography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: T2m0,     T2m1                ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Precip0,  Precip1             ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Wind_WE0, Wind_WE1            ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Wind_SN0, Wind_SN1            ! Monthly mean south_north wind speed (m/s)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Wind_LR0, Wind_LR1            ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Wind_DU0, Wind_DU1            ! Monthly mean south_north wind speed (m/s)
    INTEGER :: wHs0, wT2m0, wPrecip0, wWind_WE0, wWind_SN0, wWind_LR0, wWind_DU0
    INTEGER :: wHs1, wT2m1, wPrecip1, wWind_WE1, wWind_SN1, wWind_LR1, wWind_DU1

  END TYPE type_direct_climate_forcing_regional

  TYPE type_direct_SMB_forcing_regional
    ! Direct regional SMB forcing from a provided NetCDF file

    ! NetCDF file containing the data
    TYPE(type_netcdf_direct_SMB_forcing_regional) :: netcdf

    ! Grid
    TYPE(type_grid)                         :: grid

    ! Time dimension
    INTEGER,                    POINTER     :: nyears
    REAL(dp), DIMENSION(:    ), POINTER     :: time
    REAL(dp),                   POINTER     :: t0, t1
    INTEGER :: wnyears, wtime, wt0, wt1

    ! Actual data (two timeframes enveloping the model time) on the source grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: T2m_year0_raw, T2m_year1_raw  ! Monthly mean 2m air temperature [K] (needed for thermodynamics)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SMB_year0_raw, SMB_year1_raw  ! Yearly total SMB (m.i.e.)
    INTEGER :: wT2m_year0_raw, wT2m_year1_raw, wSMB_year0_raw, wSMB_year1_raw

    ! Actual data (two timeframes enveloping the model time) on the model mesh
    REAL(dp), DIMENSION(:    ), POINTER     :: T2m_year0, T2m_year1          ! Monthly mean 2m air temperature [K] (needed for thermodynamics)
    REAL(dp), DIMENSION(:    ), POINTER     :: SMB_year0, SMB_year1          ! Yearly total SMB (m.i.e.)
    INTEGER :: wT2m_year0, wT2m_year1, wSMB_year0, wSMB_year1

  END TYPE type_direct_SMB_forcing_regional

  TYPE type_ISMIP_style_forcing
    ! Data fields for the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing

    ! The grid used by the raw data files
    TYPE(type_grid)                         :: grid_raw

    ! The baseline elevastion, SMB, and surface temperature
    REAL(dp), DIMENSION(:    ), POINTER     :: Hs_ref
    REAL(dp), DIMENSION(:    ), POINTER     :: SMB_ref
    REAL(dp), DIMENSION(:    ), POINTER     :: ST_ref
    INTEGER :: wHs_ref, wSMB_ref, wST_ref

    ! Timestamps of the two timeframes
    REAL(dp),                   POINTER     :: t0, t1
    INTEGER :: wt0, wt1

    ! The two timeframes enveloping the model time (on the model mesh)
    REAL(dp), DIMENSION(:    ), POINTER     :: aSMB0  , aSMB1
    REAL(dp), DIMENSION(:    ), POINTER     :: dSMBdz0, dSMBdz1
    REAL(dp), DIMENSION(:    ), POINTER     :: aST0   , aST1
    REAL(dp), DIMENSION(:    ), POINTER     :: dSTdz0 , dSTdz1
    INTEGER :: waSMB0, waSMB1, wdSMBdz0, wdSMBdz1, waST0, waST1, wdSTdz0, wdSTdz1

    ! The time-interpolated values of aSMB, dSMBdz, ST, and dSTdz
    REAL(dp), DIMENSION(:    ), POINTER     :: aSMB
    REAL(dp), DIMENSION(:    ), POINTER     :: dSMBdz
    REAL(dp), DIMENSION(:    ), POINTER     :: aST
    REAL(dp), DIMENSION(:    ), POINTER     :: dSTdz
    INTEGER :: waSMB, wdSMBdz, waST, wdSTdz

    ! The applied values of SMB and ST (i.e. after applying the anomaly and elevation correction)
    REAL(dp), DIMENSION(:    ), POINTER     :: SMB
    REAL(dp), DIMENSION(:    ), POINTER     :: ST
    INTEGER :: wSMB, wST

  END TYPE type_ISMIP_style_forcing

  TYPE type_climate_matrix_regional
    ! All the relevant climate data fields (PD observations, GCM snapshots, and final, applied climate) on the model region grid

    ! The present-day observed climate (e.g. ERA40)
    TYPE(type_climate_snapshot_regional)    :: PD_obs                        ! PD observations (e.g. ERA40)

    ! The GCM climate snapshots
    TYPE(type_climate_snapshot_regional)    :: GCM_PI                        ! Pre-industrial  (e.g. HadCM3, Singarayer & Valdes, 2010), for bias correction
    TYPE(type_climate_snapshot_regional)    :: GCM_warm                      ! Warm            (e.g. HadCM3, Singarayer & Valdes, 2010)
    TYPE(type_climate_snapshot_regional)    :: GCM_cold                      ! Cold            (e.g. HadCM3, Singarayer & Valdes, 2010)
    TYPE(type_climate_snapshot_regional)    :: applied                       ! Final applied climate

    ! GCM bias
    REAL(dp), DIMENSION(:,:  ), POINTER     :: GCM_bias_T2m                  ! GCM temperature   bias (= [modelled PI temperature  ] - [observed PI temperature  ])
    REAL(dp), DIMENSION(:,:  ), POINTER     :: GCM_bias_Precip               ! GCM precipitation bias (= [modelled PI precipitation] / [observed PI precipitation])
    INTEGER :: wGCM_bias_T2m, wGCM_bias_Precip

    ! Direct climate/SMB forcing
    TYPE(type_direct_climate_forcing_regional) :: direct
    TYPE(type_direct_SMB_forcing_regional)     :: SMB_direct

    ! Data fields for the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
    TYPE(type_ISMIP_style_forcing)          :: ISMIP_style

  END TYPE type_climate_matrix_regional

  TYPE type_ocean_snapshot_regional
    ! Global ocean snapshot, either from present-day observations (e.g. WOA18) or from a GCM snapshot, projected onto the regional model mesh.

    CHARACTER(LEN=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.

    ! General forcing info (not relevant for PD observations)
    REAL(dp),                   POINTER     :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp),                   POINTER     :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   POINTER     :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp),                   POINTER     :: orbit_obl
    REAL(dp),                   POINTER     :: orbit_pre
    REAL(dp),                   POINTER     :: sealevel
    INTEGER :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre, wsealevel

    ! Ocean data
    REAL(dp),                   POINTER     :: T_ocean_mean                  ! Regional mean ocean temperature (used for basal melt when no ocean temperature data is provided)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: T_ocean                       ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: S_ocean                       ! 3-D annual mean ocean salinity    [PSU]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: T_ocean_ext                   ! 3-D annual mean ocean temperature, extrapolated beneath ice shelves [K]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: S_ocean_ext                   ! 3-D annual mean ocean salinity   , extrapolated beneath ice shelves [PSU]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: T_ocean_corr_ext              ! Bias-corrected 3-D annual mean ocean temperature, extrapolated beneath ice shelves [K]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: S_ocean_corr_ext              ! Bias-corrected 3-D annual mean ocean salinity,    extrapolated beneath ice shelves [PSU]
    CHARACTER(LEN=256)                      :: hires_ocean_foldername        ! Name of folder containing the extrapolated ocean data (not bias-corrected)
    INTEGER :: wT_ocean_mean, wT_ocean, wS_ocean, wT_ocean_ext, wS_ocean_ext, wT_ocean_corr_ext, wS_ocean_corr_ext

    ! History of the weighing fields
    REAL(dp), DIMENSION(:,:  ), POINTER     :: w_tot_history
    INTEGER,                    POINTER     :: nw_tot_history
    INTEGER :: ww_tot_history, wnw_tot_history

  END TYPE type_ocean_snapshot_regional

  TYPE type_ocean_matrix_regional
    ! All the relevant ocean data fields (PD observations, GCM snapshots, and final, applied ocean) on the model region grid

    TYPE(type_ocean_snapshot_regional)      :: PD_obs                        ! PD observations (e.g. WOA18)
    TYPE(type_ocean_snapshot_regional)      :: GCM_PI                        ! Pre-industrial reference snapshot, for bias correction
    TYPE(type_ocean_snapshot_regional)      :: GCM_warm                      ! Warm snapshot
    TYPE(type_ocean_snapshot_regional)      :: GCM_cold                      ! Cold snapshot
    TYPE(type_ocean_snapshot_regional)      :: applied                       ! Final applied ocean

  END TYPE type_ocean_matrix_regional

  TYPE type_highres_ocean_data
    ! High-resolution (extrapolated) regional ocean data

    ! NetCDF files
    TYPE( type_netcdf_reference_geometry)      :: netcdf_geo
    TYPE( type_netcdf_extrapolated_ocean_data) :: netcdf

    ! Grid
    TYPE( type_grid)                        :: grid

    ! Geometry data
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi                            ! Ice thickness [m]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb                            ! Bedrock elevation [m]
    INTEGER :: wHi, wHb

    ! Ice basins
    INTEGER,  DIMENSION(:,:  ), POINTER     :: basin_ID                      ! The drainage basin to which each grid cell belongs
    INTEGER,                    POINTER     :: nbasins                       ! Total number of basins defined for this region
    INTEGER :: wbasin_ID, wnbasins

    ! Raw and extrapolated ocean data
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T_ocean                       ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: S_ocean                       ! 3-D annual mean ocean salinity    [PSU]
    INTEGER ::wT_ocean, wS_ocean

  END TYPE type_highres_ocean_data

  ! == Mass balance
  ! ===============

  TYPE type_SMB_model
    ! The different SMB components, calculated from the prescribed climate

    ! Tuning parameters (different for each region, set from config)
    REAL(dp),                   POINTER     :: C_abl_constant
    REAL(dp),                   POINTER     :: C_abl_Ts
    REAL(dp),                   POINTER     :: C_abl_Q
    REAL(dp),                   POINTER     :: C_refr
    INTEGER :: wC_abl_constant, wC_abl_Ts, wC_abl_Q, wC_refr

    ! Inverted tuning parameters (different for each region, 1-D set from config)
    REAL(dp), DIMENSION(:    ), POINTER     :: C_abl_constant_inv
    REAL(dp), DIMENSION(:    ), POINTER     :: C_abl_Ts_inv
    REAL(dp), DIMENSION(:    ), POINTER     :: C_abl_Q_inv
    REAL(dp), DIMENSION(:    ), POINTER     :: C_refr_inv
    INTEGER :: wC_abl_constant_inv, wC_abl_Ts_inv, wC_abl_Q_inv, wC_refr_inv

    ! Data fields
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Q_TOA                         ! The prescribed monthly insolation, from an external forcing file
    REAL(dp), DIMENSION(:    ), POINTER     :: AlbedoSurf                    ! Surface albedo underneath the snow layer (water, rock or ice)
    REAL(dp), DIMENSION(:    ), POINTER     :: MeltPreviousYear              ! Total melt that occurred during the previous year (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: FirnDepth                     ! Depth of the firn layer (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Rainfall                      ! Monthly rainfall (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Snowfall                      ! Monthly snowfall (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: AddedFirn                     ! Monthly added firn (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Melt                          ! Monthly melt (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Refreezing                    ! Monthly refreezing (m)
    REAL(dp), DIMENSION(:    ), POINTER     :: Refreezing_year               ! Yearly  refreezing (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Runoff                        ! Monthly runoff (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Albedo                        ! Monthly albedo
    REAL(dp), DIMENSION(:    ), POINTER     :: Albedo_year                   ! Yearly albedo
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SMB                           ! Monthly SMB (m)
    REAL(dp), DIMENSION(:    ), POINTER     :: SMB_year                      ! Yearly  SMB (m)
    INTEGER :: wQ_TOA, wAlbedoSurf, wMeltPreviousYear, wFirnDepth, wRainfall, wSnowfall
    INTEGER :: wAddedFirn, wMelt, wRefreezing, wRefreezing_year, wRunoff, wAlbedo, wAlbedo_year, wSMB, wSMB_year

  END TYPE type_SMB_model

  TYPE type_BMB_model
    ! The different BMB components

    ! General data fields
    !====================

    REAL(dp), DIMENSION(:    ), POINTER     :: BMB                           ! The basal mass balance (same as SMB: negative means ice loss, positive means ice gain!) [m/yr]
    REAL(dp), DIMENSION(:    ), POINTER     :: BMB_sheet                     ! The basal mass balance underneath the land-based ice sheet [m/yr]
    REAL(dp), DIMENSION(:    ), POINTER     :: BMB_shelf                     ! The basal mass balance underneath the floating   ice shelf [m/yr]

    ! The ANICE_legacy BMB model
    ! ==========================

    ! Tuning parameters (different for each region, set from config)
    REAL(dp),                   POINTER     :: T_ocean_mean_PD
    REAL(dp),                   POINTER     :: T_ocean_mean_cold
    REAL(dp),                   POINTER     :: T_ocean_mean_warm
    REAL(dp),                   POINTER     :: BMB_deepocean_PD
    REAL(dp),                   POINTER     :: BMB_deepocean_cold
    REAL(dp),                   POINTER     :: BMB_deepocean_warm
    REAL(dp),                   POINTER     :: BMB_shelf_exposed_PD
    REAL(dp),                   POINTER     :: BMB_shelf_exposed_cold
    REAL(dp),                   POINTER     :: BMB_shelf_exposed_warm
    REAL(dp),                   POINTER     :: subshelf_melt_factor
    REAL(dp),                   POINTER     :: deep_ocean_threshold_depth
    INTEGER :: wT_ocean_mean_PD, wT_ocean_mean_cold, wT_ocean_mean_warm
    INTEGER :: wBMB_deepocean_PD, wBMB_deepocean_cold, wBMB_deepocean_warm
    INTEGER :: wBMB_shelf_exposed_PD, wBMB_shelf_exposed_cold, wBMB_shelf_exposed_warm
    INTEGER :: wsubshelf_melt_factor, wdeep_ocean_threshold_depth

    ! The linear/quadratic models from Favier et al. (2019)
    ! =====================================================

    REAL(dp), DIMENSION(:    ), POINTER     :: T_ocean_base                  ! Ocean temperature    at the ice shelf base
    REAL(dp), DIMENSION(:    ), POINTER     :: T_ocean_freeze_base           ! Ocean freezing point at the ice shelf base (depends on pressure and salinity)
    INTEGER :: wT_ocean_base, wT_ocean_freeze_base

    ! The Lazeroms (2018) plume model
    ! ===============================

    ! NOTE: also uses T_ocean_base!

    INTEGER,  DIMENSION(:,:  ), POINTER     :: search_directions             ! The 16 search directions
    REAL(dp), DIMENSION(:    ), POINTER     :: eff_plume_source_depth        ! Effective plume source depth (average source depth over all valid plume paths)
    REAL(dp), DIMENSION(:    ), POINTER     :: eff_basal_slope               ! Effective basal slope        (average slope        over all valid plume paths)
    INTEGER :: wsearch_directions, weff_plume_source_depth, weff_basal_slope

    ! The PICO model
    ! ==============

    ! NOTE: also uses search_directions!

    REAL(dp), DIMENSION(:    ), POINTER     :: PICO_d_GL                     ! Distance to grounding line [m]
    REAL(dp), DIMENSION(:    ), POINTER     :: PICO_d_IF                     ! Distance to ice front      [m]
    REAL(dp), DIMENSION(:    ), POINTER     :: PICO_r                        ! Relative distance to grounding line [0-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_A                        ! Area covered by each ocean box in each basin [n_basins x n_boxes]
    INTEGER,  DIMENSION(:    ), POINTER     :: PICO_n_D                      ! Number of ocean boxes for each ice basin
    INTEGER,  DIMENSION(:    ), POINTER     :: PICO_k                        ! PICO ocean box number to which the shelf grid cells belong
    INTEGER :: wPICO_d_GL, wPICO_d_IF, wPICO_r, wPICO_A, wPICO_n_D, wPICO_k

    REAL(dp), DIMENSION(:    ), POINTER     :: PICO_T                        ! 2-D     ambient temperature [K]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_Tk                       ! Average ambient temperature within each basin-box
    REAL(dp), DIMENSION(:    ), POINTER     :: PICO_S                        ! 2-D     ambient salinity    [PSU]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_Sk                       ! Average ambient salinity    within each basin-box
    REAL(dp), DIMENSION(:    ), POINTER     :: PICO_p                        ! 2-D     basal pressure      [Pa]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_pk                       ! Average basal pressure      within each basin-box
    REAL(dp), DIMENSION(:    ), POINTER     :: PICO_m                        ! 2-D     melt rate           [m/yr]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_mk                       ! Average melt rate           within each basin-box
    INTEGER :: wPICO_T, wPICO_Tk, wPICO_S, wPICO_Sk, wPICO_p, wPICO_pk, wPICO_m, wPICO_mk

    ! Additional data fields
    !=======================

    REAL(dp), DIMENSION(:    ), POINTER     :: sub_angle                     ! "subtended angle"      for the sub-shelf melt parameterisation
    REAL(dp), DIMENSION(:    ), POINTER     :: dist_open                     ! distance to open ocean for the sub-shelf melt parameterisation
    INTEGER :: wBMB, wBMB_sheet, wBMB_shelf, wsub_angle, wdist_open

  END TYPE type_BMB_model

  ! == SELEN
  ! ========

  TYPE type_SELEN_mesh
    ! The global SELEN mesh, in UFEMISM mesh format

    INTEGER,                    POINTER     :: nV                            ! Number of vertices
    INTEGER,                    POINTER     :: nTri                          ! Number of triangles
    INTEGER,                    POINTER     :: nC_mem                        ! Maximum allowed number of connections per vertex
    INTEGER :: wnV, wnTri, wnC_mem

    ! Actual mesh data
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V                             ! The X and Y coordinates of all the vertices
    INTEGER,  DIMENSION(:,:  ), POINTER     :: Tri                           ! The triangle array: Tri(ti) = [vi1, vi2, vi3] (vertices ordered counter-clockwise)
    INTEGER,  DIMENSION(:    ), POINTER     :: nC                            ! The number of other vertices this vertex is connected to
    INTEGER,  DIMENSION(:,:  ), POINTER     :: C                             ! The list   of other vertices this vertex is connected to (ordered counter-clockwise, from edge to edge for edge vertices)
    INTEGER,  DIMENSION(:    ), POINTER     :: niTri                         ! The number of triangles this vertex is a part of
    INTEGER,  DIMENSION(:,:  ), POINTER     :: iTri                          ! The list   of triangles this vertex is a part of (ordered counter-clockwise)
    INTEGER :: wV, wTri, wnC, wC, wniTri, wiTri

    ! Lat/lon
    REAL(dp), DIMENSION(:    ), POINTER     :: lat                           ! Latitude  (degrees north)
    REAL(dp), DIMENSION(:    ), POINTER     :: lon                           ! Longitude (degrees east)
    INTEGER :: wlat, wlon

    ! "anchor points" (used in some places to reduce spherical harmonics calculations)
    INTEGER,                    POINTER     :: nanc                          ! Number of anchor points
    REAL(dp), DIMENSION(:    ), POINTER     :: ancplist_lat                  ! Latitude of each anchor point
    INTEGER,  DIMENSION(:,:  ), POINTER     :: ancplist_isl                  ! Range of global grid pixels for each anchor point
    INTEGER,  DIMENSION(:    ), POINTER     :: ianc                          ! Anchor point index of each global grid pixel
    INTEGER :: wnanc, wancplist_lat, wancplist_isl, wianc

  END TYPE type_SELEN_mesh

  TYPE type_SELEN_global
    ! Global SELEN data fields, spherical harmonics, etc.

    ! Timers to trigger solving the SLE in the main program
    REAL(dp)                                :: t0_SLE, t1_SLE

    ! NetCDF global topography input
    TYPE(type_netcdf_SELEN_global_topo)     :: netcdf_topo

    ! NetCDF output
    TYPE(type_netcdf_SELEN_output)          :: output

    ! The irregular global mesh
    TYPE( type_SELEN_mesh)                  :: mesh

    ! Global reference fields
    REAL(dp), DIMENSION(:    ), POINTER     :: topo_ref                      ! Reference topography
    REAL(dp), DIMENSION(:    ), POINTER     :: load_ref                      ! Reference load
    INTEGER,  DIMENSION(:    ), POINTER     :: of_ref                        ! Reference ocean function
    INTEGER,  DIMENSION(:,:  ), POINTER     :: icemodel_region               ! To which ice model region a SELEN pixel belongs (if any), and their ice model subgrid index
    INTEGER,                    POINTER     :: nel_icemodel                  ! Number of pixels that belong to ice model regions
    INTEGER,  DIMENSION(:,:  ), POINTER     :: isl_icemodel                  ! List   of pixels that belong to ice model regions
    INTEGER :: wtopo_ref, wload_ref, wof_ref, wicemodel_region, wnel_icemodel, wisl_icemodel

    ! Global end results (ice loading, geoid, bed deformation, geoid rate, bed deformation rate)
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_glob
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_rel_glob
    REAL(dp), DIMENSION(:    ), POINTER     :: N_glob
    REAL(dp), DIMENSION(:    ), POINTER     :: U_glob
    INTEGER,  DIMENSION(:    ), POINTER     :: of_glob
    INTEGER :: wHi_glob, wHi_rel_glob, wN_glob, wU_glob, wof_glob

    ! MJ and LJ values
    INTEGER,  DIMENSION(:    ), POINTER     ::  MJ_VAL,  LJ_VAL
    INTEGER                                 :: wMJ_VAL, wLJ_VAL

    ! Reference data for SELEN in SH
    REAL(dp),   DIMENSION(:,:), POINTER     :: ALF
    COMPLEX*16, DIMENSION(:  ), POINTER     :: dLONG_TABLE
    COMPLEX*16, DIMENSION(:,:), POINTER     :: LONG_TABLE
    INTEGER :: wALF, wdLONG_TABLE, wLONG_TABLE

    ! Ice loading history on the global grid on the irregular moving time window
    REAL(dp),   DIMENSION(:    ), POINTER   :: dice_loading_history_irreg_glob
    REAL(dp),   DIMENSION(:,:  ), POINTER   :: ice_loading_history_irreg_glob
    INTEGER :: wdice_loading_history_irreg_glob, wice_loading_history_irreg_glob

    ! Memory of SELEN
    COMPLEX*16, DIMENSION(:    ), POINTER   :: MEM_S
    COMPLEX*16, DIMENSION(:    ), POINTER   :: MEM_U
    INTEGER                                 :: wMEM_S, wMEM_U

    ! SELEN internal data
    INTEGER,    DIMENSION(:    ), POINTER   :: Dm                  ! INTEGER,    DIMENSION( C_SLE%JMAX                )     :: DM
    REAL(dp),   DIMENSION(:    ), POINTER   :: X                   ! REAL(dp),   DIMENSION( C_SLE%NP,  0:C_SLE%TCONV  )     :: X
    INTEGER :: wDm, wX

    COMPLEX*16, DIMENSION(:    ), POINTER     :: ZE                  ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: ZE                  ! Eustatic Z array for water loading
    COMPLEX*16, DIMENSION(:    ), POINTER     :: SE                  ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: SE                  ! Eustatic S array for ice loading
    COMPLEX*16, DIMENSION(:    ), POINTER     :: AAAA                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: AAAA
    COMPLEX*16, DIMENSION(:    ), POINTER     :: AAAA_MOD            ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: AAAA_MOD
    COMPLEX*16, DIMENSION(:    ), POINTER     :: BBBB                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: BBBB
    COMPLEX*16, DIMENSION(:    ), POINTER     :: BBBB_MOD            ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: BBBB_MOD
    COMPLEX*16, DIMENSION(:    ), POINTER     :: HHHH                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: HHHH
    COMPLEX*16, DIMENSION(:    ), POINTER     :: KKKK                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: KKKK
    COMPLEX*16, DIMENSION(:    ), POINTER     :: TTTT                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: IIII                ! Ice loading in SH
    COMPLEX*16, DIMENSION(:    ), POINTER     :: IIII                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: TTTT                ! Relative ice loading in SH
    INTEGER :: wZE, wSE, wAAAA, wAAAA_MOD, wBBBB, wBBBB_MOD, wHHHH, wKKKK, wTTTT, wIIII

    ! New arrays for the rotational feedback
    COMPLEX*16, DIMENSION(:    ), POINTER     :: load_ice_and_water  ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: load_ice_and_water  ! Ice and water loading in SH for rotational feedback
    COMPLEX*16, DIMENSION(:    ), POINTER     :: ZROT                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: ZROT
    COMPLEX*16, DIMENSION(:    ), POINTER     :: ZROT_MOD            ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: ZROT_MOD
    COMPLEX*16, DIMENSION(:    ), POINTER     :: SROT                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: SROT
    COMPLEX*16, DIMENSION(:    ), POINTER     :: SROT_MOD            ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: SROT_MOD
    COMPLEX*16, DIMENSION(:    ), POINTER     :: ZROTVV              ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: ZROTVV
    COMPLEX*16, DIMENSION(:    ), POINTER     :: SROTVV              ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: SROTVV
    INTEGER :: wload_ice_and_water, wZROT, wZROT_MOD, wSROT, wSROT_MOD, wZROTVV, wSROTVV

    ! Results
    COMPLEX*16, DIMENSION(:    ), POINTER     :: S                   ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: S                   ! relative sea level
    COMPLEX*16, DIMENSION(:    ), POINTER     :: U                   ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: U                   ! Solid Earth deformation (N = S+U)
    COMPLEX*16, DIMENSION(:    ), POINTER     :: Z                   ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV,0:C%selen_max_iteration)     :: Z
    INTEGER :: wS, wU, wZ

    ! New arrays for tdof
    REAL(dp),   DIMENSION(:    ), POINTER     :: newalf              ! REAL(dp),   DIMENSION( C_SLE%JMAX,  C_SLE%NANCH  )     :: newalf              ! New values of alf
    REAL(dp),   DIMENSION(:    ), POINTER     :: slc                 ! REAL(dp),   DIMENSION( C_SLE%NP,  0:C_SLE%TCONV  )     :: slc                 ! sea level forward in time for all elements
    INTEGER,    DIMENSION(:    ), POINTER     :: WET                 ! INTEGER,    DIMENSION( C_SLE%NP,  0:C_SLE%TCONV  )     :: WET                 ! new ocean function (0 or 1) forward in time
    REAL(dp),   DIMENSION(:    ), POINTER     :: newtopo             ! REAL(dp),   DIMENSION( C_SLE%NP,  0:C_SLE%TCONV  )     :: newtopo             ! new topograhpy forward in time
    INTEGER :: wnewalf, wslc, wWET, wnewtopo

    INTEGER,    DIMENSION(:    ), POINTER     :: fj                  ! INTEGER,    DIMENSION( C_SLE%NP,  0:C_SLE%TCONV  )     :: fj                  ! floating function for ice thickness (float or not)
    INTEGER,    DIMENSION(:    ), POINTER     :: xcf                 ! INTEGER,    DIMENSION( C_SLE%NP,  0:C_SLE%TCONV+1)     :: xcf
    INTEGER :: wfj, wxcf

    INTEGER,    DIMENSION(:    ), POINTER     :: newet2              ! INTEGER,    DIMENSION( C_SLE%NP                  )     :: newet2              ! new ocean function after iteration
    REAL(dp),   DIMENSION(:    ), POINTER     :: newtopo2            ! REAL(dp),   DIMENSION( C_SLE%NP                  )     :: newtopo2            ! new topograhpy after iteration
    INTEGER :: wnewet2, wnewtopo2

    REAL(dp),   DIMENSION(:    ), POINTER     :: S_global, U_global  ! REAL(dp),   DIMENSION( C_SLE%NP                  )     :: S_global, U_global  ! Global fields of S and U for output
    INTEGER :: wS_global, wU_global

    COMPLEX*16, DIMENSION(:    ), POINTER     :: MSint               ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: MSint               ! Interpolated memory of relative sea level
    COMPLEX*16, DIMENSION(:    ), POINTER     :: MUint               ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: MUint               ! Interpolationd memory of solid earth deformation (N is geoid)
    INTEGER :: wMSint, wMUint

    COMPLEX*16, DIMENSION(:    ), POINTER     :: MSread              ! COMPLEX*16, DIMENSION( C_SLE%JMAX                )     :: MSread              ! Initial memory of relative sea level
    INTEGER :: wMSread

    ! New variables for spline interpolation
    REAL(dp),   DIMENSION(:    ), POINTER     :: int_time            ! REAL(dp),   DIMENSION(            0:C_SLE%TCONV  )     :: int_time            ! time points (depending on DELTA) for interpolation
    COMPLEX*16, DIMENSION(:    ), POINTER     :: sp1, spn, preSint   ! COMPLEX*16, DIMENSION( C_SLE%JMAX)                     :: sp1, spn, preSint   ! spline in- and output for S
    COMPLEX*16, DIMENSION(:    ), POINTER     :: up1, upn, preUint   ! COMPLEX*16, DIMENSION( C_SLE%JMAX)                     :: up1, upn, preUint   ! spline in- and output for U
    COMPLEX*16, DIMENSION(:    ), POINTER     :: s2,u2               ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: s2,u2               ! spline output
    INTEGER :: wint_time, wsp1, wspn, wpreSint, wup1, wupn, wpreUint, ws2, wu2

    COMPLEX*16, DIMENSION(:    ), POINTER     :: varpreoc            ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: varpreoc            ! pre determined SH coefficients ocean function
    COMPLEX*16, DIMENSION(:    ), POINTER     :: varoc               ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: varoc               ! variable SH coefficient of ocean function
    COMPLEX*16, DIMENSION(:    ), POINTER     :: varoc_inv           ! COMPLEX*16, DIMENSION(            0:C_SLE%TCONV  )     :: varoc_inv           ! 1 / varoc(1,:) - the first degree
    INTEGER :: wvarpreoc, wvaroc, wvaroc_inv

  END TYPE type_SELEN_global

  TYPE type_SELEN_regional
    ! SELEN data for this model region

    ! Mapping between the icemodel square grid and the SELEN global grid
    REAL(dp),                   POINTER :: rad                    ! Radius of the circular disc describing a square grid element (used for mapping)
    INTEGER,  DIMENSION(:,:  ), POINTER :: map_nisl               ! For each ice model square grid cell, the number of SELEN global grid cells to which it contributes
    INTEGER,  DIMENSION(:,:,:), POINTER :: map_isl                ! For each ice model square grid cell, the list   of SELEN global grid cells to which it contributes
    REAL(dp), DIMENSION(:,:,:), POINTER :: map_w                  ! The weights of those contributions
    INTEGER,                    POINTER :: nr                     ! The number of SELEN global grid pixels that lie inside this ice model region (accounting for region overlap)
    INTEGER                             :: ir1, ir2               ! Parallelisation process domains
    INTEGER,  DIMENSION(:    ), POINTER :: map_isl_region2glob    ! For each SELEN global grid pixel inside this ice model region, the corresponding global grid pixel index
    INTEGER :: wrad, wmap_nisl, wmap_isl, wmap_w, wnslr, wmap_isl_region2glob

    ! Reference fields
    REAL(dp), DIMENSION(:,:  ), POINTER :: load_ref               ! Reference loading field    (on the regional square grid)
    REAL(dp), DIMENSION(:,:  ), POINTER :: topo_ref               ! Reference topography field (on the regional square grid)
    INTEGER :: wload_ref, wtopo_ref

    ! Path to directory containing spherical harmonics binary files for this region
    CHARACTER(LEN=256)                  :: sh_foldername

    ! Ice loading history (on the regional square grid, both on the regular and irregular moving time window)
    REAL(dp), DIMENSION(:,:,:), POINTER :: ice_loading_history_reg_sq
    REAL(dp), DIMENSION(:,:,:), POINTER :: ice_loading_history_irreg_sq
    INTEGER :: wice_loading_history_reg_sq, wice_loading_history_irreg_sq

    ! SELEN results
    REAL(dp), DIMENSION(:,:  ), POINTER :: dHb_t_grid_GIA
    REAL(dp), DIMENSION(:,:  ), POINTER :: dHb_tplusdt_grid_GIA
    REAL(dp), DIMENSION(:,:  ), POINTER :: SL_t_grid_GIA
    REAL(dp), DIMENSION(:,:  ), POINTER :: SL_tplusdt_grid_GIA
    REAL(dp), DIMENSION(:    ), POINTER :: dHb_t
    REAL(dp), DIMENSION(:    ), POINTER :: dHb_tplusdt
    REAL(dp), DIMENSION(:    ), POINTER :: SL_t
    REAL(dp), DIMENSION(:    ), POINTER :: SL_tplusdt
    INTEGER :: wdHb_t_grid_GIA, wdHb_tplusdt_grid_GIA, wSL_t_grid_GIA, wSL_tplusdt_grid_GIA
    INTEGER :: wdHb_t,          wdHb_tplusdt,          wSL_t,          wSL_tplusdt

  END TYPE type_SELEN_regional

  ! == Restart
  ! ==========

  TYPE type_restart_data
    ! Restart data and NetCDF file

    ! NetCDF file
    TYPE(type_netcdf_restart)               :: netcdf

    ! Dimensions
    INTEGER,                    POINTER     :: nz, nt
    REAL(dp), DIMENSION(:    ), POINTER     :: zeta, time
    INTEGER :: wnz, wnt, wzeta, wtime

    ! Ice dynamics
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi
    REAL(dp), DIMENSION(:    ), POINTER     :: Hb
    REAL(dp), DIMENSION(:    ), POINTER     :: Hs
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Ti
    INTEGER :: wHi, wHb, wHs, wTi

    ! Bed roughness
    REAL(dp), DIMENSION(:    ), POINTER     :: beta_sq
    REAL(dp), DIMENSION(:    ), POINTER     :: phi_fric
    INTEGER :: wbeta_sq, wphi_fric

    ! GIA
    REAL(dp), DIMENSION(:    ), POINTER     :: SL
    REAL(dp), DIMENSION(:    ), POINTER     :: dHb
    INTEGER :: wSL, wdHb

    ! SMB
    REAL(dp), DIMENSION(:,:  ), POINTER     :: FirnDepth
    REAL(dp), DIMENSION(:    ), POINTER     :: MeltPreviousYear
    REAL(dp), DIMENSION(:    ), POINTER     :: C_abl_constant_inv
    REAL(dp), DIMENSION(:    ), POINTER     :: C_abl_Ts_inv
    REAL(dp), DIMENSION(:    ), POINTER     :: C_abl_Q_inv
    REAL(dp), DIMENSION(:    ), POINTER     :: C_refr_inv
    INTEGER :: wFirnDepth, wMeltPreviousYear, wC_abl_constant_inv, wC_abl_Ts_inv, wC_abl_Q_inv, wC_refr_inv

    ! Isotopes
    REAL(dp), DIMENSION(:,:  ), POINTER     :: IsoIce
    INTEGER :: wIsoIce

  END TYPE type_restart_data

  ! == Model region
  ! ===============

  TYPE type_model_region
    ! Contains all the different data structures, organised by sub-model (ice, climate)

    ! Metadata
    CHARACTER(LEN=3)                        :: name                                           ! NAM, EAS, GRL, ANT
    CHARACTER(LEN=256)                      :: long_name                                      ! North America, Eurasia, Greenland, Antarctica

    ! The current time (step) of this particular region.
    REAL(dp), POINTER                       :: time
    REAL(dp), POINTER                       :: dt
    REAL(dp), POINTER                       :: dt_prev
    INTEGER :: wtime, wdt, wdt_prev

    ! Timers and switches for determining which modules need to be called at what points in time during the simulation
    REAL(dp), POINTER                       :: dt_crit_SIA
    REAL(dp), POINTER                       :: dt_crit_SSA
    REAL(dp), POINTER                       :: dt_crit_ice, dt_crit_ice_prev
    REAL(dp), POINTER                       :: t_last_mesh,    t_next_mesh
    REAL(dp), POINTER                       :: t_last_SIA,     t_next_SIA
    REAL(dp), POINTER                       :: t_last_SSA,     t_next_SSA
    REAL(dp), POINTER                       :: t_last_DIVA,    t_next_DIVA
    REAL(dp), POINTER                       :: t_last_thermo,  t_next_thermo
    REAL(dp), POINTER                       :: t_last_output,  t_next_output
    REAL(dp), POINTER                       :: t_last_climate, t_next_climate
    REAL(dp), POINTER                       :: t_last_ocean,   t_next_ocean
    REAL(dp), POINTER                       :: t_last_SMB,     t_next_SMB
    REAL(dp), POINTER                       :: t_last_BMB,     t_next_BMB
    REAL(dp), POINTER                       :: t_last_ELRA,    t_next_ELRA
    REAL(dp), POINTER                       :: t_last_basal,   t_next_basal
    REAL(dp), POINTER                       :: t_last_SMB_inv, t_next_SMB_inv
    LOGICAL,  POINTER                       :: do_mesh
    LOGICAL,  POINTER                       :: do_SIA
    LOGICAL,  POINTER                       :: do_SSA
    LOGICAL,  POINTER                       :: do_DIVA
    LOGICAL,  POINTER                       :: do_thermo
    LOGICAL,  POINTER                       :: do_climate
    LOGICAL,  POINTER                       :: do_ocean
    LOGICAL,  POINTER                       :: do_SMB
    LOGICAL,  POINTER                       :: do_BMB
    LOGICAL,  POINTER                       :: do_output
    LOGICAL,  POINTER                       :: do_ELRA
    LOGICAL,  POINTER                       :: do_basal
    LOGICAL,  POINTER                       :: do_SMB_inv
    INTEGER :: wdt_crit_SIA, wdt_crit_SSA, wdt_crit_ice, wdt_crit_ice_prev
    INTEGER :: wt_last_mesh, wt_last_SIA, wt_last_SSA, wt_last_DIVA, wt_last_thermo, wt_last_output, wt_last_climate, wt_last_ocean, wt_last_SMB, wt_last_BMB, wt_last_ELRA, wt_last_basal, wt_last_SMB_inv
    INTEGER :: wt_next_mesh, wt_next_SIA, wt_next_SSA, wt_next_DIVA, wt_next_thermo, wt_next_output, wt_next_climate, wt_next_ocean, wt_next_SMB, wt_next_BMB, wt_next_ELRA, wt_next_basal, wt_next_SMB_inv
    INTEGER ::     wdo_mesh,     wdo_SIA,     wdo_SSA,     wdo_DIVA,     wdo_thermo,     wdo_output,     wdo_climate,     wdo_ocean,     wdo_SMB,     wdo_BMB,     wdo_ELRA,     wdo_basal,     wdo_SMB_inv

    ! The region's ice sheet's volume and volume above flotation (in mSLE, so the second one is the ice sheets GMSL contribution)
    REAL(dp), POINTER                       :: ice_area
    REAL(dp), POINTER                       :: ice_volume
    REAL(dp), POINTER                       :: ice_volume_PD
    REAL(dp), POINTER                       :: ice_volume_above_flotation
    REAL(dp), POINTER                       :: ice_volume_above_flotation_PD
    REAL(dp), POINTER                       :: GMSL_contribution
    INTEGER :: wice_area, wice_volume, wice_volume_PD, wice_volume_above_flotation, wice_volume_above_flotation_PD, wGMSL_contribution

    ! Regionally integrated mass balance components
    REAL(dp), POINTER                       :: int_T2m
    REAL(dp), POINTER                       :: int_snowfall
    REAL(dp), POINTER                       :: int_rainfall
    REAL(dp), POINTER                       :: int_melt
    REAL(dp), POINTER                       :: int_refreezing
    REAL(dp), POINTER                       :: int_runoff
    REAL(dp), POINTER                       :: int_SMB
    REAL(dp), POINTER                       :: int_BMB
    REAL(dp), POINTER                       :: int_MB
    INTEGER :: wint_T2m, wint_snowfall, wint_rainfall, wint_melt, wint_refreezing, wint_runoff, wint_SMB, wint_BMB, wint_MB

    ! Variables related to the englacial isotope content
    REAL(dp), POINTER                       :: mean_isotope_content
    REAL(dp), POINTER                       :: mean_isotope_content_PD
    REAL(dp), POINTER                       :: d18O_contribution
    REAL(dp), POINTER                       :: d18O_contribution_PD
    INTEGER :: wmean_isotope_content, wmean_isotope_content_PD, wd18O_contribution, wd18O_contribution_PD

    ! Reference geometries
    TYPE(type_reference_geometry)           :: refgeo_init                               ! Initial         ice-sheet geometry
    TYPE(type_reference_geometry)           :: refgeo_PD                                 ! Present-day     ice-sheet geometry
    TYPE(type_reference_geometry)           :: refgeo_GIAeq                              ! GIA equilibrium ice-sheet geometry

    ! Mask where ice is not allowed to form (so Greenland is not included in NAM and EAS, and Ellesmere is not included in GRL)
    INTEGER,  DIMENSION(:), POINTER         :: mask_noice
    INTEGER                                 :: wmask_noice

    ! Sub-models
    TYPE(type_mesh)                         :: mesh                                      ! The finite element mesh for this model region
    TYPE(type_mesh)                         :: mesh_new                                  ! The new mesh after updating (so that the old one can be kept until data has been mapped)
    TYPE(type_ice_model)                    :: ice                                       ! All the ice model data for this model region
    TYPE(type_climate_matrix_regional)      :: climate_matrix                            ! All the climate data for this model region (new version)
    TYPE(type_ocean_matrix_regional)        :: ocean_matrix                              ! All the ocean data for this model region
    TYPE(type_SMB_model)                    :: SMB                                       ! The different SMB components for this model region
    TYPE(type_BMB_model)                    :: BMB                                       ! The different BMB components for this model region
    TYPE(type_SELEN_regional)               :: SELEN            ! SELEN input and output data for this model region

    ! Output netcdf files
    LOGICAL                                 :: output_file_exists
    TYPE(type_netcdf_restart)               :: restart_mesh
    TYPE(type_netcdf_restart)               :: restart_grid
    TYPE(type_netcdf_help_fields)           :: help_fields_mesh
    TYPE(type_netcdf_help_fields)           :: help_fields_grid

    ! Different square grids
    TYPE(type_grid)                         :: grid_output                               ! For the "_grid" output files
    TYPE(type_grid)                         :: grid_GIA                                  ! For either the ELRA model or SELEN
    TYPE(type_grid)                         :: grid_smooth                               ! For smoothing data fields (used in the climate matrix)

    ! Computation times
    REAL(dp), POINTER                       :: tcomp_total
    REAL(dp), POINTER                       :: tcomp_ice
    REAL(dp), POINTER                       :: tcomp_thermo
    REAL(dp), POINTER                       :: tcomp_climate
    REAL(dp), POINTER                       :: tcomp_GIA
    REAL(dp), POINTER                       :: tcomp_mesh
    INTEGER :: wtcomp_total, wtcomp_ice, wtcomp_thermo, wtcomp_climate, wtcomp_GIA, wtcomp_mesh

    ! Restart data
    TYPE(type_restart_data)                 :: restart

  END TYPE type_model_region

  ! == Debug
  ! ========

  TYPE type_debug_fields
    ! Dummy variables for debugging

    ! NetCDF debug file
    TYPE(type_netcdf_debug)                 :: netcdf

    ! NetCDF debug file name
    CHARACTER(LEN=256)                      :: filename

    ! Data
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_a_01
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_a_02
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_a_03
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_a_04
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_a_05
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_a_06
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_a_07
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_a_08
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_a_09
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_a_10
    INTEGER :: wint_2D_a_01, wint_2D_a_02, wint_2D_a_03, wint_2D_a_04, wint_2D_a_05
    INTEGER :: wint_2D_a_06, wint_2D_a_07, wint_2D_a_08, wint_2D_a_09, wint_2D_a_10

    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_b_01
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_b_02
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_b_03
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_b_04
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_b_05
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_b_06
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_b_07
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_b_08
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_b_09
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_b_10
    INTEGER :: wint_2D_b_01, wint_2D_b_02, wint_2D_b_03, wint_2D_b_04, wint_2D_b_05
    INTEGER :: wint_2D_b_06, wint_2D_b_07, wint_2D_b_08, wint_2D_b_09, wint_2D_b_10

    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_c_01
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_c_02
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_c_03
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_c_04
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_c_05
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_c_06
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_c_07
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_c_08
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_c_09
    INTEGER,  DIMENSION(:    ), POINTER     :: int_2D_c_10
    INTEGER :: wint_2D_c_01, wint_2D_c_02, wint_2D_c_03, wint_2D_c_04, wint_2D_c_05
    INTEGER :: wint_2D_c_06, wint_2D_c_07, wint_2D_c_08, wint_2D_c_09, wint_2D_c_10

    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_a_01
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_a_02
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_a_03
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_a_04
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_a_05
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_a_06
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_a_07
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_a_08
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_a_09
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_a_10
    INTEGER :: wdp_2D_a_01, wdp_2D_a_02, wdp_2D_a_03, wdp_2D_a_04, wdp_2D_a_05
    INTEGER :: wdp_2D_a_06, wdp_2D_a_07, wdp_2D_a_08, wdp_2D_a_09, wdp_2D_a_10

    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_b_01
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_b_02
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_b_03
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_b_04
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_b_05
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_b_06
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_b_07
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_b_08
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_b_09
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_b_10
    INTEGER :: wdp_2D_b_01, wdp_2D_b_02, wdp_2D_b_03, wdp_2D_b_04, wdp_2D_b_05
    INTEGER :: wdp_2D_b_06, wdp_2D_b_07, wdp_2D_b_08, wdp_2D_b_09, wdp_2D_b_10

    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_c_01
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_c_02
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_c_03
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_c_04
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_c_05
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_c_06
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_c_07
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_c_08
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_c_09
    REAL(dp), DIMENSION(:    ), POINTER     :: dp_2D_c_10
    INTEGER :: wdp_2D_c_01, wdp_2D_c_02, wdp_2D_c_03, wdp_2D_c_04, wdp_2D_c_05
    INTEGER :: wdp_2D_c_06, wdp_2D_c_07, wdp_2D_c_08, wdp_2D_c_09, wdp_2D_c_10

    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_3D_a_01
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_3D_a_02
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_3D_a_03
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_3D_a_04
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_3D_a_05
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_3D_a_06
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_3D_a_07
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_3D_a_08
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_3D_a_09
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_3D_a_10
    INTEGER :: wdp_3D_a_01, wdp_3D_a_02, wdp_3D_a_03, wdp_3D_a_04, wdp_3D_a_05
    INTEGER :: wdp_3D_a_06, wdp_3D_a_07, wdp_3D_a_08, wdp_3D_a_09, wdp_3D_a_10

    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_monthly_a_01
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_monthly_a_02
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_monthly_a_03
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_monthly_a_04
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_monthly_a_05
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_monthly_a_06
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_monthly_a_07
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_monthly_a_08
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_monthly_a_09
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_monthly_a_10
    INTEGER :: wdp_2D_monthly_a_01, wdp_2D_monthly_a_02, wdp_2D_monthly_a_03, wdp_2D_monthly_a_04, wdp_2D_monthly_a_05
    INTEGER :: wdp_2D_monthly_a_06, wdp_2D_monthly_a_07, wdp_2D_monthly_a_08, wdp_2D_monthly_a_09, wdp_2D_monthly_a_10

  END TYPE type_debug_fields

  TYPE type_memory_use_tracker

    ! Memory use history
    INTEGER(KIND=MPI_ADDRESS_KIND)      :: total                     ! Total amount of allocated shared memory (in bytes)
    INTEGER                             :: n                         ! Number of entries
    INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(:), ALLOCATABLE :: h   ! Memory use history over the past coupling interval

  END TYPE type_memory_use_tracker

CONTAINS

END MODULE data_types_module