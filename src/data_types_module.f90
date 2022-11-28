module data_types_module
  ! Contains all the different types for storing data. Put all together in a separate module so that
  ! all subroutines can use all types without interdependency conflicts, and also to make the
  ! modules with the actual physics code more readable.
  ! If only Types could be collapsed in BBEdit...

#include <petsc/finclude/petscksp.h>
  USE petscksp
  USE mpi
  USE configuration_module,        ONLY: dp, C
  USE data_types_netcdf_module,    ONLY: type_netcdf_climate_data, type_netcdf_reference_geometry, &
                                         type_netcdf_insolation, type_netcdf_restart, type_netcdf_help_fields, &
                                         type_netcdf_debug, type_netcdf_ICE5G_data, type_netcdf_geothermal_heat_flux, &
                                         type_netcdf_direct_climate_forcing_global, type_netcdf_direct_SMB_forcing_global, &
                                         type_netcdf_direct_climate_forcing_regional, type_netcdf_direct_SMB_forcing_regional, &
                                         type_netcdf_ocean_data, type_netcdf_extrapolated_ocean_data, &
                                         type_netcdf_resource_tracker, type_netcdf_scalars_global, type_netcdf_scalars_regional

  IMPLICIT NONE

  TYPE type_sparse_matrix_CSR_dp
    ! Compressed Sparse Row (CSR) format matrix

    INTEGER                                 :: m,n                         ! A = [m-by-n]
    INTEGER                                 :: nnz_max                     ! Maximum number of non-zero entries in A (determines how much memory is allocated)
    INTEGER                                 :: nnz                         ! Number         of non-zero entries in A locally (determines how much memory is allocated)
    integer                                 :: nnz_tot                     ! globally nnz entries
    INTEGER,  DIMENSION(:    ), allocatable :: ptr
    integer                                 :: ptr_i1, ptr_i2              ! set by allocate, which rows are on this processor, defaults to normal partition_list
    logical                                 :: balanced = .true.           ! SSADIVA is not 'balanced' like petsc wants it to be
    INTEGER,  DIMENSION(:    ), allocatable :: index
    REAL(dp), DIMENSION(:    ), allocatable :: val

  END TYPE type_sparse_matrix_CSR_dp

  TYPE type_ice_model
    ! The ice dynamics sub-model data structure.

    ! Basic data - ice thickness, bedrock & surface elevation, sea level (geoid elevation), englacial temperature, and ice velocities
    REAL(dp), DIMENSION(:    ), allocatable :: Hi_a                        ! Ice thickness [m]
    REAL(dp), DIMENSION(:    ), allocatable :: Hb_a                        ! Bedrock elevation [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:    ), allocatable :: Hs_a                        ! Surface elevation [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:    ), allocatable :: SL_a                        ! Sea level (geoid elevation) [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:    ), allocatable :: TAF_a                       ! Thickness above flotation [m]
    REAL(dp), DIMENSION(:,:  ), allocatable :: Ti_a                        ! Englacial temperature [K]

    ! Ice velocities
    REAL(dp), DIMENSION(:,:  ), allocatable :: u_3D_a                      ! 3-D ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), allocatable :: v_3D_a
    REAL(dp), DIMENSION(:,:  ), allocatable :: u_3D_b
    REAL(dp), DIMENSION(:,:  ), allocatable :: v_3D_b
    REAL(dp), DIMENSION(:,:  ), allocatable :: w_3D_a
    REAL(dp), DIMENSION(:,:  ), allocatable :: w_3D_b

    REAL(dp), DIMENSION(:    ), allocatable :: u_vav_a                     ! Vertically averaged ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:    ), allocatable :: v_vav_a
    REAL(dp), DIMENSION(:    ), allocatable :: u_vav_b
    REAL(dp), DIMENSION(:    ), allocatable :: v_vav_b
    REAL(dp), DIMENSION(:    ), allocatable :: uabs_vav_a
    REAL(dp), DIMENSION(:    ), allocatable :: uabs_vav_b

    REAL(dp), DIMENSION(:    ), allocatable :: u_surf_a                    ! Ice velocity at the surface [m yr^-1]
    REAL(dp), DIMENSION(:    ), allocatable :: v_surf_a
    REAL(dp), DIMENSION(:    ), allocatable :: u_surf_b
    REAL(dp), DIMENSION(:    ), allocatable :: v_surf_b
    REAL(dp), DIMENSION(:    ), allocatable :: uabs_surf_a
    REAL(dp), DIMENSION(:    ), allocatable :: uabs_surf_b

    REAL(dp), DIMENSION(:    ), allocatable :: u_base_a                    ! Ice velocity at the base [m yr^-1]
    REAL(dp), DIMENSION(:    ), allocatable :: v_base_a
    REAL(dp), DIMENSION(:    ), allocatable :: u_base_b
    REAL(dp), DIMENSION(:    ), allocatable :: v_base_b
    REAL(dp), DIMENSION(:    ), allocatable :: uabs_base_a
    REAL(dp), DIMENSION(:    ), allocatable :: uabs_base_b

    REAL(dp), DIMENSION(:,:  ), allocatable :: u_3D_SIA_b
    REAL(dp), DIMENSION(:,:  ), allocatable :: v_3D_SIA_b
    REAL(dp), DIMENSION(:    ), allocatable :: u_base_SSA_b
    REAL(dp), DIMENSION(:    ), allocatable :: v_base_SSA_b

    ! Different masks
    INTEGER,  DIMENSION(:    ), allocatable :: mask_land_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_ocean_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_lake_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_ice_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_sheet_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_shelf_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_coast_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_margin_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_gl_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_glf_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_cf_a
    INTEGER,  DIMENSION(:    ), allocatable :: mask_a
    REAL(dp), DIMENSION(:    ), allocatable :: f_grnd_a
    REAL(dp), DIMENSION(:    ), allocatable :: f_grnd_b
    INTEGER,  DIMENSION(:    ), allocatable :: basin_ID                    ! The drainage basin to which each grid cell belongs
    INTEGER                                 :: nbasins                     ! Total number of basins defined for this region

    ! Ice physical properties
    REAL(dp), DIMENSION(:,:  ), allocatable :: A_flow_3D_a                 ! Flow parameter [Pa^-3 y^-1]
    REAL(dp), DIMENSION(:    ), allocatable :: A_flow_vav_a
    REAL(dp), DIMENSION(:,:  ), allocatable :: Ti_pmp_a                    ! The pressure melting point temperature [K]
    REAL(dp), DIMENSION(:,:  ), allocatable :: Cpi_a                       ! Specific heat capacity of ice [J kg^-1 K^-1].
    REAL(dp), DIMENSION(:,:  ), allocatable :: Ki_a                        ! Conductivity of ice [J m^-1 K^-1 yr^-1].

    ! Zeta derivatives
    REAL(dp), DIMENSION(:,:  ), allocatable :: dzeta_dt_a
    REAL(dp), DIMENSION(:,:  ), allocatable :: dzeta_dx_a
    REAL(dp), DIMENSION(:,:  ), allocatable :: dzeta_dy_a
    REAL(dp), DIMENSION(:    ), allocatable :: dzeta_dz_a

    ! Ice dynamics - basal hydroallocatable
    REAL(dp), DIMENSION(:    ), allocatable :: overburden_pressure_a       ! Overburden pressure ( = H * rho_i * g) [Pa]
    REAL(dp), DIMENSION(:    ), allocatable :: pore_water_pressure_a       ! Pore water pressure (determined by basal hydrology model) [Pa]
    REAL(dp), DIMENSION(:    ), allocatable :: Neff_a                      ! Effective pressure ( = overburden pressure - pore water pressure) [Pa]

    ! Ice dynamics - basal roughness / friction
    REAL(dp), DIMENSION(:    ), allocatable :: phi_fric_a                  ! Till friction angle (degrees)
    REAL(dp), DIMENSION(:    ), allocatable :: phi_fric_inv_a              ! Inverted till friction angle (degrees)
    REAL(dp), DIMENSION(:    ), allocatable :: tauc_a                      ! Till yield stress tauc   (used when choice_sliding_law = "Coloumb" or "Coulomb_regularised")
    REAL(dp), DIMENSION(:    ), allocatable :: alpha_sq_a                  ! Coulomb-law friction coefficient [unitless]         (used when choice_sliding_law =             "Tsai2015", or "Schoof2005")
    REAL(dp), DIMENSION(:    ), allocatable :: beta_sq_a                   ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3] (used when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005")

    ! Ice dynamics - physical terms in the SSA/DIVA
    REAL(dp), DIMENSION(:    ), allocatable :: taudx_b                     ! x-component of the driving stress
    REAL(dp), DIMENSION(:    ), allocatable :: taudy_b                     ! x-component of the driving stress
    REAL(dp), DIMENSION(:    ), allocatable :: du_dx_a                     ! Vertically averaged   xx strain rate
    REAL(dp), DIMENSION(:    ), allocatable :: du_dy_a                     ! Vertically averaged   xy strain rate
    REAL(dp), DIMENSION(:    ), allocatable :: dv_dx_a                     ! Vertically averaged   yy strain rate
    REAL(dp), DIMENSION(:    ), allocatable :: dv_dy_a                     ! Vertically averaged   yy strain rate
    REAL(dp), DIMENSION(:,:  ), allocatable :: du_dz_3D_b                  ! 3-D                   xz strain rate
    REAL(dp), DIMENSION(:,:  ), allocatable :: dv_dz_3D_b                  ! 3-D                   yz strain rate
    REAL(dp), DIMENSION(:,:  ), allocatable :: visc_eff_3D_a               ! 3-D                   effective viscosity
    REAL(dp), DIMENSION(:    ), allocatable :: visc_eff_int_a              ! Vertically integrated effective viscosity
    REAL(dp), DIMENSION(:    ), allocatable :: N_a                         ! Product term N = eta * H
    REAL(dp), DIMENSION(:    ), allocatable :: beta_a                      ! Sliding term beta (as in, [basal shear stress] = [beta] * [basal velocity])
    REAL(dp), DIMENSION(:    ), allocatable :: beta_eff_a                  ! Beta_eff, appearing in the DIVA
    REAL(dp), DIMENSION(:    ), allocatable :: beta_eff_b
    REAL(dp), DIMENSION(:    ), allocatable :: taubx_b                     ! x-component of the basal shear stress
    REAL(dp), DIMENSION(:    ), allocatable :: tauby_b                     ! y-component of the basal shear stress
    REAL(dp), DIMENSION(:    ), allocatable :: F2_a                        ! F2, appearing in the DIVA
    REAL(dp), DIMENSION(:    ), allocatable :: u_prev_b
    REAL(dp), DIMENSION(:    ), allocatable :: v_prev_b

    ! Ice dynamics - some administrative stuff to make solving the SSA/DIVA more efficient
    INTEGER,  DIMENSION(:    ), allocatable :: ti2n_u, ti2n_v
    INTEGER,  DIMENSION(:,:  ), allocatable :: n2ti_uv
    TYPE(type_sparse_matrix_CSR_dp)         :: M_SSADIVA

    ! Ice dynamics - ice thickness calculation
    REAL(dp), DIMENSION(:,:  ), allocatable :: dVi_in
    REAL(dp), DIMENSION(:    ), allocatable :: dHs_dt_a
    REAL(dp), DIMENSION(:    ), allocatable :: Hi_tplusdt_a
    REAL(dp), DIMENSION(:    ), allocatable :: dHi_dt_a

    ! Ice dynamics - calving
    REAL(dp), DIMENSION(:    ), allocatable :: float_margin_frac_a         ! Ice-covered fraction for calving front pixels
    REAL(dp), DIMENSION(:    ), allocatable :: Hi_eff_cf_a                 ! Effective ice thickness at calving front pixels (= Hi of thinnest non-calving-front neighbour)

    ! Ice dynamics - predictor/corrector ice thickness update
    REAL(dp)                                :: pc_zeta
    REAL(dp), DIMENSION(:    ), allocatable :: pc_tau
    REAL(dp), DIMENSION(:    ), allocatable :: pc_fcb
    REAL(dp)                                :: pc_eta = 1.
    REAL(dp)                                :: pc_eta_prev
    REAL(dp)                                :: pc_beta1
    REAL(dp)                                :: pc_beta2
    REAL(dp)                                :: pc_beta3
    REAL(dp)                                :: pc_beta4
    REAL(dp), DIMENSION(:    ), allocatable :: pc_f1
    REAL(dp), DIMENSION(:    ), allocatable :: pc_f2
    REAL(dp), DIMENSION(:    ), allocatable :: pc_f3
    REAL(dp), DIMENSION(:    ), allocatable :: pc_f4
    REAL(dp), DIMENSION(:    ), allocatable :: Hi_old
    REAL(dp), DIMENSION(:    ), allocatable :: Hi_pred
    REAL(dp), DIMENSION(:    ), allocatable :: Hi_corr

    ! Thermodynamics
    INTEGER,  DIMENSION(:    ), allocatable :: mask_ice_a_prev        ! Ice mask from previous time step
    REAL(dp), DIMENSION(:,:  ), allocatable :: internal_heating_a     ! Internal heating due to deformation
    REAL(dp), DIMENSION(:    ), allocatable :: frictional_heating_a   ! Friction heating due to basal sliding
    REAL(dp), DIMENSION(:    ), allocatable :: GHF_a                  ! Geothermal heat flux

    ! Isotope content
    REAL(dp), DIMENSION(:    ), allocatable :: Hi_a_prev
    REAL(dp), DIMENSION(:    ), allocatable :: IsoRef
    REAL(dp), DIMENSION(:    ), allocatable :: IsoSurf
    REAL(dp), DIMENSION(:    ), allocatable :: MB_iso
    REAL(dp), DIMENSION(:    ), allocatable :: IsoIce
    REAL(dp), DIMENSION(:    ), allocatable :: IsoIce_new

    ! ELRA GIA model
    INTEGER                                 :: flex_prof_rad
    REAL(dp), DIMENSION(:,:  ), allocatable :: flex_prof_grid
    REAL(dp), DIMENSION(:    ), allocatable :: surface_load_PD_mesh
    REAL(dp), DIMENSION(:    ), allocatable :: surface_load_mesh
    REAL(dp), DIMENSION(:    ), allocatable :: surface_load_rel_mesh
    REAL(dp), DIMENSION(:,:  ), allocatable :: surface_load_rel_grid
    REAL(dp), DIMENSION(:,:  ), allocatable :: surface_load_rel_ext_grid
    REAL(dp), DIMENSION(:,:  ), allocatable :: dHb_eq_grid
    REAL(dp), DIMENSION(:,:  ), allocatable :: dHb_grid
    REAL(dp), DIMENSION(:,:  ), allocatable :: dHb_dt_grid
    REAL(dp), DIMENSION(:    ), allocatable :: dHb_a
    REAL(dp), DIMENSION(:    ), allocatable :: dHb_dt_a
    REAL(dp), DIMENSION(:    ), allocatable :: dSL_dt_a

    ! Mesh adaptation data
    REAL(dp), DIMENSION(:    ), allocatable :: surf_curv
    REAL(dp), DIMENSION(:    ), allocatable :: log_velocity

    ! Useful extra stuff
    REAL(dp), DIMENSION(:    ), allocatable :: dHi_a   ! Ice thickness difference w.r.t. PD
    REAL(dp), DIMENSION(:    ), allocatable :: dHs_a   ! Ice elevation difference w.r.t. PD

  END TYPE type_ice_model

  TYPE type_mesh
    ! The unstructured triangular mesh.

    ! Basic meta properties
    ! =====================

    CHARACTER(LEN=3)                        :: region_name                   ! NAM, EAS, GRL, ANT
    REAL(dp)                                :: lambda_M                      ! Oblique stereographic projection parameters
    REAL(dp)                                :: phi_M
    REAL(dp)                                :: alpha_stereo
    REAL(dp)                                :: xmin                          ! X and Y range of the square covered by the mesh
    REAL(dp)                                :: xmax
    REAL(dp)                                :: ymin
    REAL(dp)                                :: ymax
    REAL(dp)                                :: tol_dist                      ! Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    INTEGER                                 :: nV_mem                        ! Size of allocated memory for vertices
    INTEGER                                 :: nTri_mem                      ! Size of allocated memory for triangles
    INTEGER                                 :: nC_mem                        ! Maximum allowed number of connections per vertex
    INTEGER                                 :: nV                            ! Number of vertices
    INTEGER                                 :: nTri                          ! Number of triangles
    INTEGER                                 :: perturb_dir                   ! Perturbation direction (0 = anticlockwise, 1 = clockwise)
    REAL(dp)                                :: alpha_min = 0                 ! Sharpest inner angle allowed by Rupperts algorithm (default because it was assumed by shared alloc)
    REAL(dp)                                :: dz_max_ice                    ! Maximum allowed vertical error in Rupperts algorithm over ice
    REAL(dp)                                :: res_max                       ! Maximum resolution anywhere
    REAL(dp)                                :: res_max_margin                ! Maximum resolution over the ice margin
    REAL(dp)                                :: res_max_gl                    ! Maximum resolution over the grounding line
    REAL(dp)                                :: res_max_cf                    ! Maximum resolution over the calving front
    REAL(dp)                                :: res_max_mountain              ! Maximum resolution over ice-free mountains  (important for getting the inception right)
    REAL(dp)                                :: res_max_coast                 ! Maximum resolution over ice-free coast line (to make plots look nicer)
    REAL(dp)                                :: res_min                       ! Minimum resolution anywhere ( = MINVAL([res_max_margin, res_max_gl, res_max_cf]))
    REAL(dp)                                :: resolution_min                ! Finest   resolution of the mesh ( = MINVAL(R), where R = distance to nearest neighbour)
    REAL(dp)                                :: resolution_max                ! Coarsest resolution of the mesh

    ! Primary mesh data (needed for mesh creation & refinement)
    ! =========================================================

    ! Vertex data
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: V                             ! The X and Y coordinates of all the vertices
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: nC                            ! The number of other vertices this vertex is connected to
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: C                             ! The list   of other vertices this vertex is connected to (ordered counter-clockwise, from edge to edge for edge vertices)
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: niTri                         ! The number of triangles this vertex is a part of
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: iTri                          ! The list   of triangles this vertex is a part of (ordered counter-clockwise)
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: edge_index                    ! Each vertex's Edge Index; 0 = free, 1 = north, 2 = northeast, etc.
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: mesh_old_ti_in                ! When creating a new mesh: for every new mesh vertex, the old mesh triangle that contains it

    ! Triangle data
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: Tri                           ! The triangle array: Tri(ti) = [vi1, vi2, vi3] (vertices ordered counter-clockwise)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: TriCC                         ! The X,Y-coordinates of each triangle's circumcenter
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: TriC                          ! The (up to) three neighbour triangles (order across from 1st, 2nd and 3d vertex, respectively)

    ! Refinement lists
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: Triflip                       ! List of triangles to flip, used in updating Delaunay triangulation
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: RefMap                        ! Map   of which triangles      have been marked for refining
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: RefStack                      ! Stack of       triangles that have been marked for refining
    INTEGER                                 :: RefStackN

    ! Maps+stacks for FloodFill-ALLOCATABLEhing
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: VMap, VStack1, VStack2
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: TriMap, TriStack1, TriStack2
    INTEGER                                 :: VStackN1, VStackN2
    INTEGER                                 :: TriStackN1, TriStackN2

    ! Points-of-Interest
    INTEGER,                    ALLOCATABLE :: nPOI                          ! Number of Points of Interest (POI) in this mesh
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: POI_coordinates               ! Lat-lon coordinates of a POI
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: POI_XY_coordinates            ! X-Y     coordinates of a POI
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: POI_resolutions               ! Resolution          of a POI (km)
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: POI_vi                        ! The three vertices surrounding a POI
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: POI_w                         ! Their relative weights in trilinear interpolation

    ! Secondary mesh data
    ! ===================

    ! Derived geometry data
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: A                             ! The area             of each vertex's Voronoi cell
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: VorGC                         ! The geometric centre of each vertex's Voronoi cell
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: R                             ! The resolution (defined as distance to nearest neighbour)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Cw                            ! The width of all vertex connections (= length of the shared Voronoi cell edge)
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: TriGC                         ! The X,Y-coordinates of each triangle's geometric centre
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: Tri_edge_index                ! Same as for vertices, only NE, SE, etc. aren't used
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: TriA                          ! The area of each triangle

    ! Staggered (Arakawa C) meshALLOCATABLE
    INTEGER,                    ALLOCATABLE :: nAc                           ! The number of Ac vertices
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: VAc                           ! x,y coordinates of the Ac vertices
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: Aci                           ! Mapping array from the Aa to the Ac mesh
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: iAci                          ! Mapping array from the Ac to the Aa mesh
    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: edge_index_Ac

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
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: lat
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: lon

    ! Transect
    INTEGER                                 :: nV_transect                   ! Number of vertex pairs for the transect
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE :: vi_transect                   ! List   of vertex pairs for the transect
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: w_transect                    ! Interpolation weights for the vertex pairs

    ! Parallelisation
    INTEGER                                 :: vi1, vi2                      ! Vertices
    INTEGER                                 :: ti1, ti2                      ! Triangles
    INTEGER                                 :: ci1, ci2                      ! Edges

  END TYPE type_mesh

  TYPE type_debug_fields
    ! Dummy variables for debugging

    ! NetCDF debug file
    TYPE(type_netcdf_debug)                 :: netcdf

    ! Data
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_a_01
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_a_02
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_a_03
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_a_04
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_a_05
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_a_06
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_a_07
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_a_08
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_a_09
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_a_10

    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_b_01
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_b_02
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_b_03
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_b_04
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_b_05
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_b_06
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_b_07
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_b_08
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_b_09
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_b_10

    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_c_01
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_c_02
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_c_03
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_c_04
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_c_05
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_c_06
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_c_07
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_c_08
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_c_09
    INTEGER,  DIMENSION(:    ), pointer     :: int_2D_c_10

    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_a_01
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_a_02
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_a_03
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_a_04
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_a_05
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_a_06
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_a_07
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_a_08
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_a_09
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_a_10

    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_b_01
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_b_02
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_b_03
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_b_04
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_b_05
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_b_06
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_b_07
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_b_08
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_b_09
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_b_10

    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_c_01
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_c_02
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_c_03
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_c_04
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_c_05
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_c_06
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_c_07
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_c_08
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_c_09
    REAL(dp), DIMENSION(:    ), pointer     :: dp_2D_c_10

    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_3D_a_01
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_3D_a_02
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_3D_a_03
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_3D_a_04
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_3D_a_05
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_3D_a_06
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_3D_a_07
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_3D_a_08
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_3D_a_09
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_3D_a_10

    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_2D_monthly_a_01
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_2D_monthly_a_02
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_2D_monthly_a_03
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_2D_monthly_a_04
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_2D_monthly_a_05
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_2D_monthly_a_06
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_2D_monthly_a_07
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_2D_monthly_a_08
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_2D_monthly_a_09
    REAL(dp), DIMENSION(:,:  ), pointer     :: dp_2D_monthly_a_10

  END TYPE type_debug_fields

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

  TYPE type_grid
    ! A regular square grid covering a model region

    ! Basic grid data
    INTEGER                                 :: nx, ny, n
    REAL(dp),                   allocatable :: dx
    REAL(dp), DIMENSION(:    ), allocatable :: x, y
    REAL(dp)                                :: xmin, xmax, ymin, ymax

    ! Parallelisation by domain decomposition
    INTEGER                                 :: i1, i2, j1, j2

    ! Sparse matrices representing the remapping operations between a mesh and a grid
    TYPE(tMat)                              :: M_map_grid2mesh              ! Remapping from a grid to a mesh using second-order conservative remapping
    TYPE(tMat)                              :: M_map_mesh2grid              ! Remapping from a mesh to a grid using second-order conservative remapping
    REAL(dp)                                :: tol_dist

    ! Conversion tables for grid-form vs. vector-form data
    INTEGER,  DIMENSION(:,:  ), allocatable :: ij2n, n2ij

    ! Lat-lon coordinates
    REAL(dp), DIMENSION(:,:  ), allocatable :: lat, lon

    ! Projection parameters for the grid
    REAL(dp)                                :: lambda_M
    REAL(dp)                                :: phi_M
    REAL(dp)                                :: alpha_stereo

  END TYPE type_grid

  TYPE type_remapping_latlon2mesh
    ! Indices and weights for mapping data from a global lat-lon grid to the model mesh using bilinear interpolation

    INTEGER,  DIMENSION(:    ), ALLOCATABLE :: ilat1, ilat2, ilon1, ilon2
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: wlat1, wlat2, wlon1, wlon2

  END TYPE type_remapping_latlon2mesh

  TYPE type_latlongrid
    ! A global lat-lon grid

    INTEGER                                 :: nlat, nlon
    REAL(dp), DIMENSION(:    ), allocatable :: lat, lon
    REAL(dp)                                :: dlat, dlon

    INTEGER                                 :: i1, i2 ! Parallelisation by domain decomposition

  END TYPE type_latlongrid

  TYPE type_SMB_model
    ! The different SMB components, calculated from the prescribed climate

    ! Tuning parameters (different for each region, set from config)
    REAL(dp)                                :: C_abl_constant
    REAL(dp)                                :: C_abl_Ts
    REAL(dp)                                :: C_abl_Q
    REAL(dp)                                :: C_refr

    ! Data fields
    REAL(dp), DIMENSION(:,:  ), allocatable :: Q_TOA                         ! The prescribed monthly insolation, from an external forcing file
    REAL(dp), DIMENSION(:    ), allocatable :: AlbedoSurf                    ! Surface albedo underneath the snow layer (water, rock or ice)
    REAL(dp), DIMENSION(:    ), allocatable :: MeltPreviousYear              ! Total melt that occurred during the previous year (m)
    REAL(dp), DIMENSION(:,:  ), allocatable :: FirnDepth                     ! Depth of the firn layer (m)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Rainfall                      ! Monthly rainfall (m)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Snowfall                      ! Monthly snowfall (m)
    REAL(dp), DIMENSION(:,:  ), allocatable :: AddedFirn                     ! Monthly added firn (m)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Melt                          ! Monthly melt (m)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Refreezing                    ! Monthly refreezing (m)
    REAL(dp), DIMENSION(:    ), allocatable :: Refreezing_year               ! Yearly  refreezing (m)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Runoff                        ! Monthly runoff (m)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Albedo                        ! Monthly albedo
    REAL(dp), DIMENSION(:    ), allocatable :: Albedo_year                   ! Yearly albedo
    REAL(dp), DIMENSION(:,:  ), allocatable :: SMB                           ! Monthly SMB (m)
    REAL(dp), DIMENSION(:    ), allocatable :: SMB_year                      ! Yearly  SMB (m)

  END TYPE type_SMB_model

  TYPE type_BMB_model
    ! The different BMB components

    ! General data fields
    !====================

    REAL(dp), DIMENSION(:    ), allocatable :: BMB                           ! The basal mass balance (same as SMB: negative means ice loss, positive means ice gain!) [m/yr]
    REAL(dp), DIMENSION(:    ), allocatable :: BMB_sheet                     ! The basal mass balance underneath the land-based ice sheet [m/yr]
    REAL(dp), DIMENSION(:    ), allocatable :: BMB_shelf                     ! The basal mass balance underneath the floating   ice shelf [m/yr]

    ! The ANICE_legacy BMB modelallocatable
    ! ==========================allocatable

    ! Tuning parameters (differeallocatable region, set from config)
    REAL(dp)                                :: T_ocean_mean_PD
    REAL(dp)                                :: T_ocean_mean_cold
    REAL(dp)                                :: T_ocean_mean_warm
    REAL(dp)                                :: BMB_deepocean_PD
    REAL(dp)                                :: BMB_deepocean_cold
    REAL(dp)                                :: BMB_deepocean_warm
    REAL(dp)                                :: BMB_shelf_exposed_PD
    REAL(dp)                                :: BMB_shelf_exposed_cold
    REAL(dp)                                :: BMB_shelf_exposed_warm
    REAL(dp)                                :: subshelf_melt_factor
    REAL(dp)                                :: deep_ocean_threshold_depth

    ! The linear/quadratic modelallocatableer et al. (2019)
    ! ==========================allocatable================

    REAL(dp), DIMENSION(:    ), allocatable :: T_ocean_base                  ! Ocean temperature    at the ice shelf base
    REAL(dp), DIMENSION(:    ), allocatable :: T_ocean_freeze_base           ! Ocean freezing point at the ice shelf base (depends on pressure and salinity)

    ! The Lazeroms (2018) plume allocatable
    ! ==========================allocatable

    ! NOTE: also uses T_ocean_baallocatable

    INTEGER,  DIMENSION(:,:  ), allocatable :: search_directions             ! The 16 search directions
    REAL(dp), DIMENSION(:    ), allocatable :: eff_plume_source_depth        ! Effective plume source depth (average source depth over all valid plume paths)
    REAL(dp), DIMENSION(:    ), allocatable :: eff_basal_slope               ! Effective basal slope        (average slope        over all valid plume paths)

    ! The PICO model
    ! ==============

    ! NOTE: also uses search_directions!

    REAL(dp), DIMENSION(:    ), allocatable :: PICO_d_GL                     ! Distance to grounding line [m]
    REAL(dp), DIMENSION(:    ), allocatable :: PICO_d_IF                     ! Distance to ice front      [m]
    REAL(dp), DIMENSION(:    ), allocatable :: PICO_r                        ! Relative distance to grounding line [0-1]
    REAL(dp), DIMENSION(:,:  ), allocatable :: PICO_A                        ! Area covered by each ocean box in each basin [n_basins x n_boxes]
    INTEGER,  DIMENSION(:    ), allocatable :: PICO_n_D                      ! Number of ocean boxes for each ice basin
    INTEGER,  DIMENSION(:    ), allocatable :: PICO_k                        ! PICO ocean box number to which the shelf grid cells belong

    REAL(dp), DIMENSION(:    ), allocatable :: PICO_T                        ! 2-D     ambient temperature [K]
    REAL(dp), DIMENSION(:,:  ), allocatable :: PICO_Tk                       ! Average ambient temperature within each basin-box
    REAL(dp), DIMENSION(:    ), allocatable :: PICO_S                        ! 2-D     ambient salinity    [PSU]
    REAL(dp), DIMENSION(:,:  ), allocatable :: PICO_Sk                       ! Average ambient salinity    within each basin-box
    REAL(dp), DIMENSION(:    ), allocatable :: PICO_p                        ! 2-D     basal pressure      [Pa]
    REAL(dp), DIMENSION(:,:  ), allocatable :: PICO_pk                       ! Average basal pressure      within each basin-box
    REAL(dp), DIMENSION(:    ), allocatable :: PICO_m                        ! 2-D     melt rate           [m/yr]
    REAL(dp), DIMENSION(:,:  ), allocatable :: PICO_mk                       ! Average melt rate           within each basin-box

    ! Additional data fields
    !=======================

    REAL(dp), DIMENSION(:    ), allocatable :: sub_angle                     ! "subtended angle"      for the sub-shelf melt parameterisation
    REAL(dp), DIMENSION(:    ), allocatable :: dist_open                     ! distance to open ocean for the sub-shelf melt parameterisation

  END TYPE type_BMB_model

  TYPE type_reference_geometry
    ! Data structure containing a reference ice-sheet geometry (either schematic or read from an external file).

    ! NetCDF file containing the data
    TYPE(type_netcdf_reference_geometry)    :: netcdf

    ! Raw data as read from a NetCDF file
    TYPE(type_grid)                         :: grid
    REAL(dp), DIMENSION(:,:  ), allocatable :: Hi_grid
    REAL(dp), DIMENSION(:,:  ), allocatable :: Hb_grid
    REAL(dp), DIMENSION(:,:  ), allocatable :: Hs_grid

    ! Derived data on the grid (allocatablevature and masks, needed for mesh creation)
    REAL(dp), DIMENSION(:,:  ), allocatable :: surf_curv
    INTEGER,  DIMENSION(:,:  ), allocatable :: mask_land
    INTEGER,  DIMENSION(:,:  ), allocatable :: mask_ocean
    INTEGER,  DIMENSION(:,:  ), allocatable :: mask_ice
    INTEGER,  DIMENSION(:,:  ), allocatable :: mask_sheet
    INTEGER,  DIMENSION(:,:  ), allocatable :: mask_shelf
    INTEGER,  DIMENSION(:,:  ), allocatable :: mask_margin
    INTEGER,  DIMENSION(:,:  ), allocatable :: mask_gl
    INTEGER,  DIMENSION(:,:  ), allocatable :: mask_cf
    INTEGER,  DIMENSION(:,:  ), allocatable :: mask_coast

    ! Data on the model mesh
    REAL(dp), DIMENSION(:    ), allocatable :: Hi
    REAL(dp), DIMENSION(:    ), allocatable :: Hb
    REAL(dp), DIMENSION(:    ), allocatable :: Hs

  END TYPE type_reference_geometry

  TYPE type_forcing_data
    ! Data structure containing model forcing data - CO2 record, d18O record, (global) insolation record

    ! Data for the inverse routine
    REAL(dp)                                :: dT_glob                                   ! Modelled global mean annual surface temperature anomaly w.r.t. PD
    REAL(dp), DIMENSION(:    ), allocatable :: dT_glob_history                           ! Time window (length set in config) listing previous dT_glob values
    INTEGER                                 :: ndT_glob_history                          ! Number of entries (= length of time window / dt_coupling)
    REAL(dp)                                :: dT_deepwater                              ! Modelled deep-water temperature anomaly (= window averaged dT_glob * scaling factor), scaling factor set in config

    REAL(dp)                                :: d18O_NAM, d18O_EAS, d18O_GRL, d18O_ANT    ! Modelled benthic d18O contributions from ice volume in the model regions
    REAL(dp)                                :: d18O_from_ice_volume_mod                  ! Modelled benthic d18O contribution from global ice volume
    REAL(dp)                                :: d18O_from_temperature_mod                 ! Modelled benthic d18O contribution from deep-water temperature change
    REAL(dp)                                :: d18O_mod                                  ! Modelled benthic d18O

    REAL(dp)                                :: dT_glob_inverse                           ! Global mean annual surface temperature anomaly resulting from the inverse method
    REAL(dp), DIMENSION(:    ), allocatable :: dT_glob_inverse_history                   ! Time window (length set in config) listing previous dT_glob values
    INTEGER                                 :: ndT_glob_inverse_history                  ! Number of entries (= length of time window / dt_coupling)
    REAL(dp)                                :: CO2_inverse                               ! CO2 resulting from the inverse method
    REAL(dp), DIMENSION(:    ), allocatable :: CO2_inverse_history                       ! Time window (length set in config) listing previous CO2_inverse values
    INTEGER                                 :: nCO2_inverse_history                      ! Number of entries (= length of time window / dt_coupling)
    REAL(dp)                                :: CO2_mod                                   ! Either equal to CO2_obs or to CO2_inverse, for easier writing to output.

    ! External forcing: CO2 record
    REAL(dp), DIMENSION(:    ), allocatable :: CO2_time
    REAL(dp), DIMENSION(:    ), allocatable :: CO2_record
    REAL(dp)                                :: CO2_obs

    ! External forcing: d18O record
    REAL(dp), DIMENSION(:    ), allocatable :: d18O_time
    REAL(dp), DIMENSION(:    ), allocatable :: d18O_record
    REAL(dp)                                :: d18O_obs
    REAL(dp)                                :: d18O_obs_PD

    ! External forcing: insolation
    TYPE(type_netcdf_insolation)            :: netcdf_ins
    INTEGER                                 :: ins_nyears
    INTEGER                                 :: ins_nlat
    REAL(dp), DIMENSION(:    ), allocatable :: ins_time
    REAL(dp), DIMENSION(:    ), allocatable :: ins_lat
    REAL(dp)                                :: ins_t0, ins_t1
    REAL(dp), DIMENSION(:,:  ), allocatable :: ins_Q_TOA0, ins_Q_TOA1

    ! External forcing: geothermal heat flux
    TYPE(type_netcdf_geothermal_heat_flux)  :: netcdf_ghf
    TYPE(type_latlongrid)                   :: grid_ghf
    REAL(dp), DIMENSION(:,:  ), allocatable :: ghf_ghf

    ! External forcing: sea level record
    real(dp), dimension(:    ), allocatable :: sealevel_time
    real(dp), dimension(:    ), allocatable :: sealevel_record
    real(dp)                                :: sealevel_obs

  END TYPE type_forcing_data

  TYPE type_memory_use_tracker

    ! Memory use history
    INTEGER(KIND=MPI_ADDRESS_KIND)      :: total                     ! Total amount of allocated shared memory (in bytes)
    INTEGER                             :: n                         ! Number of entries
    INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(:), ALLOCATABLE :: h   ! Memory use history over the past coupling interval

  END TYPE type_memory_use_tracker

  ! Global climate types
  !=====================

  TYPE type_climate_snapshot_global
    ! Global climate snapshot, either from present-day observations (e.g. ERA40) or from a GCM snapshot.

    CHARACTER(LEN=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.

    ! NetCDF file containing the data
    TYPE(type_netcdf_climate_data)          :: netcdf

    ! Grid
    INTEGER                                 :: nlat, nlon
    REAL(dp), DIMENSION(:    ), allocatable :: lat
    REAL(dp), DIMENSION(:    ), allocatable :: lon

    ! General forcing info (not relevant for PD observations)
    REAL(dp)                                :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp)                                :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp)                                :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp)                                :: orbit_obl
    REAL(dp)                                :: orbit_pre

    ! Actual GCM data
    REAL(dp), DIMENSION(:,:  ), allocatable :: Hs                            ! Orography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:,:), allocatable :: T2m                           ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:,:), allocatable :: Precip                        ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:,:), allocatable :: Wind_WE                       ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), allocatable :: Wind_SN                       ! Monthly mean south_north wind speed (m/s)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Mask_ice                      ! Ice mask: 1 ice, 0 no ice

    ! Paralelisation
    INTEGER                                 :: i1, i2                        ! Grid domain (:,i1:i2) of each process

  END TYPE type_climate_snapshot_global

  TYPE type_climate_matrix_global
    ! The climate matrix data structure. Contains all the different global climate snapshots.

    ! The present-day observed climate (e.g. ERA40)
    TYPE(type_climate_snapshot_global)       :: PD_obs

    ! The GCM climate snapshots
    TYPE(type_climate_snapshot_global)       :: GCM_PI
    TYPE(type_climate_snapshot_global)       :: GCM_warm
    TYPE(type_climate_snapshot_global)       :: GCM_cold

  END TYPE type_climate_matrix_global

  ! Regional climate types
  !=======================

  TYPE type_climate_snapshot_regional
    ! Global climate snapshot, either from present-day observations (e.g. ERA40) or from a GCM snapshot, projected onto the regional model mesh.

    CHARACTER(LEN=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.

    ! General forcing info (not relevant for PD observations)
    REAL(dp)                                :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp)                                :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp)                                :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp)                                :: orbit_obl
    REAL(dp)                                :: orbit_pre
    REAL(dp)                                :: sealevel

    ! Actual observed / GCM data (read from external NetCDF file)
    REAL(dp), DIMENSION(:    ), allocatable :: Hs                            ! Orography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:  ), allocatable :: T2m                           ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Precip                        ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Wind_WE                       ! Monthly mean west-east   wind speed (m/s)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Wind_SN                       ! Monthly mean south-north wind speed (m/s)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Wind_LR                       ! Monthly mean wind speed in the x-direction (m/s)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Wind_DU                       ! Monthly mean wind speed in the y-direction (m/s)
    REAL(dp), DIMENSION(:    ), allocatable :: Mask_ice                      ! Ice mask: 1 ice, 0 no ice

    ! Spatially variable lapse rallocatable snapshots (see Berends et al., 2018)
    REAL(dp), DIMENSION(:    ), allocatable :: lambda

    ! Bias-corrected GCM data   allocatable
    REAL(dp), DIMENSION(:,:  ), allocatable :: T2m_corr                      ! Bias-corrected monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Precip_corr                   ! Bias-corrected monthly mean precipitation (m)
    REAL(dp), DIMENSION(:    ), allocatable :: Hs_corr                       ! Bias-corrected surface elevation (m)

    ! Reference absorbed insolatallocatableM snapshots), or insolation at model time for the applied climate
    REAL(dp), DIMENSION(:,:  ), allocatable :: Q_TOA                         ! Monthly mean insolation at the top of the atmosphere (W/m2) (taken from the prescribed insolation solution at orbit_time)
    REAL(dp), DIMENSION(:,:  ), allocatable :: Albedo                        ! Monthly mean surface albedo (calculated using our own SMB scheme for consistency)
    REAL(dp), DIMENSION(:    ), allocatable :: I_abs                         ! Total yearly absorbed insolation, used in the climate matrix for interpolation

  END TYPE type_climate_snapshot_regional

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
    REAL(dp), DIMENSION(:,:  ), allocatable :: GCM_bias_T2m                  ! GCM temperature   bias (= [modelled PI temperature  ] - [observed PI temperature  ])
    REAL(dp), DIMENSION(:,:  ), allocatable :: GCM_bias_Precip               ! GCM precipitation bias (= [modelled PI precipitation] / [observed PI precipitation])

  END TYPE type_climate_matrix_regional

  TYPE type_ocean_snapshot_global
    ! Global ocean snapshot, either from present-day observations (e.g. WOA18) or from a GCM ocean snapshot.

    CHARACTER(LEN=256)                      :: name                          ! 'WOA', 'COSMOS_LGM', etc.

    ! NetCDF file containing the data
    TYPE(type_netcdf_ocean_data)            :: netcdf

    ! Grid
    INTEGER                                 :: nlat, nlon
    REAL(dp), DIMENSION(:    ), allocatable :: lat
    REAL(dp), DIMENSION(:    ), allocatable :: lon

    ! General forcing info (not relevant for PD observations)
    REAL(dp)                                :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp)                                :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp)                                :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp)                                :: orbit_obl
    REAL(dp)                                :: orbit_pre

    ! Global 3-D ocean data on the original vertical grid
    REAL(dp)                                :: T_ocean_mean                  ! Regional mean ocean temperature (used for basal melt when no ocean temperature data is provided)
    REAL(dp), DIMENSION(:    ), allocatable :: z_ocean_raw                   ! Vertical coordinate of the 3-D ocean fields [m below sea surface]
    INTEGER                                 :: nz_ocean_raw                  ! Number of vertical layers in the 3-D ocean fields
    REAL(dp), DIMENSION(:,:,:), allocatable :: T_ocean_raw                   ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:,:), allocatable :: S_ocean_raw                   ! 3-D annual mean ocean salinity    [PSU]

    ! Global 3-D ocean data on tallocatablel vertical grid
    REAL(dp), DIMENSION(:,:,:), allocatable :: T_ocean                       ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:,:), allocatable :: S_ocean                       ! 3-D annual mean ocean salinity    [PSU]

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

  type type_ocean_snapshot_regional
    ! Global ocean snapshot, either from present-day observations (e.g. WOA18) or from a GCM snapshot, projected onto the regional model mesh.

    character(len=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.

    ! General forcing info (not relevant for PD observations)
    real(dp)                                :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    real(dp)                                :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    real(dp)                                :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    real(dp)                                :: orbit_obl
    real(dp)                                :: orbit_pre
    real(dp)                                :: sealevel

    ! Ocean data
    real(dp)                                :: T_ocean_mean                  ! Regional mean ocean temperature (used for basal melt when no ocean temperature data is provided)
    real(dp), dimension(:,:  ), allocatable :: T_ocean                       ! 3-D annual mean ocean temperature [K]
    real(dp), dimension(:,:  ), allocatable :: S_ocean                       ! 3-D annual mean ocean salinity    [PSU]
    real(dp), dimension(:,:  ), allocatable :: T_ocean_ext                   ! 3-D annual mean ocean temperature, extrapolated beneath ice shelves [K]
    real(dp), dimension(:,:  ), allocatable :: S_ocean_ext                   ! 3-D annual mean ocean salinity   , extrapolated beneath ice shelves [PSU]
    real(dp), dimension(:,:  ), allocatable :: T_ocean_corr_ext              ! Bias-corrected 3-D annual mean ocean temperature, extrapolated beneath ice shelves [K]
    real(dp), dimension(:,:  ), allocatable :: S_ocean_corr_ext              ! Bias-corrected 3-D annual mean ocean salinity,    extrapolated beneath ice shelves [PSU]
    real(dp), dimension(:    ), allocatable :: T_ocean_inv                   ! Inverted annual mean ocean temperature at the base of ice shelves [K]
    character(len=256)                      :: hires_ocean_foldername        ! Name of folder containing the extrapolated ocean data (not bias-corrected)

    ! History of the weighing fiallocatable
    real(dp), dimension(:,:  ), allocatable :: w_tot_history
    integer                                 :: nw_tot_history

  end type type_ocean_snapshot_regional

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
    REAL(dp), DIMENSION(:,:  ), allocatable :: Hi                            ! Ice thickness [m]
    REAL(dp), DIMENSION(:,:  ), allocatable :: Hb                            ! Bedrock elevation [m]

    ! Ice basins
    INTEGER,  DIMENSION(:,:  ), allocatable :: basin_ID                      ! The drainage basin to which each grid cell belongs
    INTEGER                                 :: nbasins                       ! Total number of basins defined for this region

    ! Raw and extrapolated ocean data
    REAL(dp), DIMENSION(:,:,:), allocatable :: T_ocean                       ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:,:), allocatable :: S_ocean                       ! 3-D annual mean ocean salinity    [PSU]

  END TYPE type_highres_ocean_data

  ! Updated model region type including new climate matrix version (move it up at the end of the climate port)
  ! ==========================================================================================================

  TYPE type_model_region
    ! Contains all the different data structures, organised by sub-model (ice, climate)

    ! Metadata
    CHARACTER(LEN=3)                        :: name                                           ! NAM, EAS, GRL, ANT
    CHARACTER(LEN=256)                      :: long_name                                      ! North America, Eurasia, Greenland, Antarctica

    ! The current time (step) of this particular region.
    REAL(dp)                                :: time
    REAL(dp)                                :: dt
    REAL(dp)                                :: dt_prev

    ! Timers and switches for determining which modules need to be called at what points in time during the simulation
    REAL(dp)                                :: dt_crit_SIA
    REAL(dp)                                :: dt_crit_SSA
    REAL(dp)                                :: dt_crit_ice, dt_crit_ice_prev
    REAL(dp)                                :: t_last_mesh,    t_next_mesh
    REAL(dp)                                :: t_last_SIA,     t_next_SIA
    REAL(dp)                                :: t_last_SSA,     t_next_SSA
    REAL(dp)                                :: t_last_DIVA,    t_next_DIVA
    REAL(dp)                                :: t_last_thermo,  t_next_thermo
    REAL(dp)                                :: t_last_output,  t_next_output
    REAL(dp)                                :: t_last_climate, t_next_climate
    REAL(dp)                                :: t_last_ocean,   t_next_ocean
    REAL(dp)                                :: t_last_SMB,     t_next_SMB
    REAL(dp)                                :: t_last_BMB,     t_next_BMB
    REAL(dp)                                :: t_last_ELRA,    t_next_ELRA
    REAL(dp)                                :: t_last_slid_inv, t_next_slid_inv
    LOGICAL                                 :: do_mesh
    LOGICAL                                 :: do_SIA
    LOGICAL                                 :: do_SSA
    LOGICAL                                 :: do_DIVA
    LOGICAL                                 :: do_thermo
    LOGICAL                                 :: do_climate
    LOGICAL                                 :: do_ocean
    LOGICAL                                 :: do_SMB
    LOGICAL                                 :: do_BMB
    LOGICAL                                 :: do_output
    LOGICAL                                 :: do_ELRA
    LOGICAL                                 :: do_slid_inv

    ! The region's ice sheet's volume and volume above flotation (in mSLE, so the second one is the ice sheets GMSL contribution)
    REAL(dp)                                :: ice_area
    REAL(dp)                                :: ice_volume
    REAL(dp)                                :: ice_volume_PD
    REAL(dp)                                :: ice_volume_above_flotation
    REAL(dp)                                :: ice_volume_above_flotation_PD
    REAL(dp)                                :: GMSL_contribution

    ! Regionally integrated mass balance components
    REAL(dp)                                :: int_T2m
    REAL(dp)                                :: int_snowfall
    REAL(dp)                                :: int_rainfall
    REAL(dp)                                :: int_melt
    REAL(dp)                                :: int_refreezing
    REAL(dp)                                :: int_runoff
    REAL(dp)                                :: int_SMB
    REAL(dp)                                :: int_BMB
    REAL(dp)                                :: int_MB

    ! Variables related to the englacial isotope content
    REAL(dp)                                :: mean_isotope_content
    REAL(dp)                                :: mean_isotope_content_PD
    REAL(dp)                                :: d18O_contribution
    REAL(dp)                                :: d18O_contribution_PD

    ! Reference geometries
    TYPE(type_reference_geometry)           :: refgeo_init                               ! Initial         ice-sheet geometry
    TYPE(type_reference_geometry)           :: refgeo_PD                                 ! Present-day     ice-sheet geometry
    TYPE(type_reference_geometry)           :: refgeo_GIAeq                              ! GIA equilibrium ice-sheet geometry

    ! Mask where ice is not allowed to form (so Greenland is not included in NAM and EAS, and Ellesmere is not included in GRL)
    INTEGER,  DIMENSION(:), allocatable     :: mask_noice

    ! Sub-models
    TYPE(type_mesh)                         :: mesh                                      ! The finite element mesh for this model region
    TYPE(type_mesh)                         :: mesh_new                                  ! The new mesh after updating (so that the old one can be kept until data has been mapped)
    TYPE(type_ice_model)                    :: ice                                       ! All the ice model data for this model region
    TYPE(type_climate_matrix_regional)      :: climate_matrix                            ! All the climate data for this model region (new version)
    TYPE(type_ocean_matrix_regional)        :: ocean_matrix                              ! All the ocean data for this model region
    TYPE(type_SMB_model)                    :: SMB                                       ! The different SMB components for this model region
    TYPE(type_BMB_model)                    :: BMB                                       ! The different BMB components for this model region

    ! Output netcdf files
    LOGICAL                                 :: output_file_exists
    TYPE(type_netcdf_restart)               :: restart_mesh
    TYPE(type_netcdf_restart)               :: restart_grid
    TYPE(type_netcdf_help_fields)           :: help_fields_mesh
    TYPE(type_netcdf_help_fields)           :: help_fields_grid
    TYPE(type_netcdf_scalars_regional)      :: scalars

    ! Different square grids
    TYPE(type_grid)                         :: grid_output                               ! For the "_grid" output files
    TYPE(type_grid)                         :: grid_GIA                                  ! For either the ELRA model or SELEN
    TYPE(type_grid)                         :: grid_smooth                               ! For smoothing data fields (used in the climate matrix)

    ! Computation times
    REAL(dp)                                :: tcomp_total    = 0.
    REAL(dp)                                :: tcomp_ice      = 0.
    REAL(dp)                                :: tcomp_thermo   = 0.
    REAL(dp)                                :: tcomp_climate  = 0.
    REAL(dp)                                :: tcomp_GIA      = 0.
    REAL(dp)                                :: tcomp_mesh     = 0.

  END TYPE type_model_region

  TYPE type_restart_data
    ! Restart data and NetCDF file

    ! NetCDF file
    TYPE(type_netcdf_restart)               :: netcdf

    ! Grid
    TYPE(type_grid)                         :: grid       ! Needed for the mapping from grid to mesh
    INTEGER                                 :: nz, nt
    REAL(dp), DIMENSION(:    ), allocatable :: zeta, time

    ! Data

    ! Ice dynamics
    REAL(dp), DIMENSION(:,:  ), allocatable :: Hi
    REAL(dp), DIMENSION(:,:  ), allocatable :: Hb
    REAL(dp), DIMENSION(:,:  ), allocatable :: Hs
    REAL(dp), DIMENSION(:,:,:), allocatable :: Ti

    ! GIA
    REAL(dp), DIMENSION(:,:  ), allocatable :: SL
    REAL(dp), DIMENSION(:,:  ), allocatable :: dHb

    ! SMB
    REAL(dp), DIMENSION(:,:,:), allocatable :: FirnDepth
    REAL(dp), DIMENSION(:,:  ), allocatable :: MeltPreviousYear

    ! Isotopes
    REAL(dp), DIMENSION(:,:  ), allocatable :: IsoIce

  END TYPE type_restart_data

  type type_global_scalar_data
    ! Structure containing some global scalar values: sea level, CO2, d18O components, computation times, etc.

    ! Netcdf file
    type(type_netcdf_scalars_global)        :: netcdf

    ! Sea level
    real(dp)                                :: GMSL                                      ! Global mean sea level change
    real(dp)                                :: GMSL_NAM                                  ! Global mean sea level change (contribution from ice in North America)
    real(dp)                                :: GMSL_EAS                                  ! Global mean sea level change (contribution from ice in Eurasia)
    real(dp)                                :: GMSL_GRL                                  ! Global mean sea level change (contribution from ice in Greenland)
    real(dp)                                :: GMSL_ANT                                  ! Global mean sea level change (contribution from ice in Antarctica)

    ! CO2
    real(dp)                                :: CO2_obs                                   ! Observed atmospheric CO2
    real(dp)                                :: CO2_mod                                   ! Modelled atmospheric CO2

    ! d18O
    real(dp)                                :: d18O_obs                                  ! Observed benthic d18O
    real(dp)                                :: d18O_mod                                  ! Modelled benthic d18O
    real(dp)                                :: d18O_ice                                  ! Contribution to benthic d18O from global ice volume
    real(dp)                                :: d18O_Tdw                                  ! Contribution to benthic d18O from deep-water temperature
    real(dp)                                :: d18O_NAM                                  ! Contribution to benthic d18O from ice in North America
    real(dp)                                :: d18O_EAS                                  ! Contribution to benthic d18O from ice in Eurasia
    real(dp)                                :: d18O_GRL                                  ! Contribution to benthic d18O from ice in Greenland
    real(dp)                                :: d18O_ANT                                  ! Contribution to benthic d18O from ice in Antarctica

    ! Temperature
    real(dp)                                :: dT_glob                                   ! Global mean annual surface temperature change
    real(dp)                                :: dT_dw                                     ! Deep-water temperature change

    ! Computation times for all regions combined
    real(dp)                                :: tcomp_total    = 0.
    real(dp)                                :: tcomp_ice      = 0.
    real(dp)                                :: tcomp_thermo   = 0.
    real(dp)                                :: tcomp_climate  = 0.
    real(dp)                                :: tcomp_GIA      = 0.

  end type type_global_scalar_data

contains

end module data_types_module
