MODULE data_types_module
  ! Contains all the different types for storing data. Put all together in a separate module so that
  ! all subroutines can use all types without interdependency conflicts, and also to make the 
  ! modules with the actual physics code more readable.
  ! If only Types could be collapsed in BBEdit...

  USE configuration_module,        ONLY: dp, C

  IMPLICIT NONE
  
  TYPE type_ice_model
    ! The ice dynamics sub-model data structure.
    
    ! Basic data - ice thickness, bedrock height, surface height, mask, surface ice velocities and temperature
    REAL(dp), DIMENSION(:    ), POINTER :: Hi
    REAL(dp), DIMENSION(:    ), POINTER :: Hb
    REAL(dp), DIMENSION(:    ), POINTER :: Hs
    REAL(dp), DIMENSION(:    ), POINTER :: sealevel
    REAL(dp), DIMENSION(:    ), POINTER :: Us
    REAL(dp), DIMENSION(:    ), POINTER :: Vs
    REAL(dp), DIMENSION(:,:  ), POINTER :: Ti
    INTEGER                             :: wHi, wHb, wHs, wsealevel, wUs, wVs, wTi
    
    ! Different masks
    INTEGER,  DIMENSION(:    ), POINTER :: mask
    INTEGER,  DIMENSION(:    ), POINTER :: mask_land
    INTEGER,  DIMENSION(:    ), POINTER :: mask_ocean
    INTEGER,  DIMENSION(:    ), POINTER :: mask_lake
    INTEGER,  DIMENSION(:    ), POINTER :: mask_ice
    INTEGER,  DIMENSION(:    ), POINTER :: mask_sheet
    INTEGER,  DIMENSION(:    ), POINTER :: mask_shelf
    INTEGER,  DIMENSION(:    ), POINTER :: mask_coast
    INTEGER,  DIMENSION(:    ), POINTER :: mask_margin
    INTEGER,  DIMENSION(:    ), POINTER :: mask_groundingline
    INTEGER,  DIMENSION(:    ), POINTER :: mask_calvingfront
    INTEGER                             :: wmask, wmask_land, wmask_ocean, wmask_lake, wmask_ice, wmask_sheet, wmask_shelf
    INTEGER                             :: wmask_coast, wmask_margin, wmask_groundingline, wmask_calvingfront
    
    ! Ice physical properties
    REAL(dp), DIMENSION(:,:  ), POINTER :: Ti_pmp      ! The pressure melting point temperature [K]
    REAL(dp), DIMENSION(:,:  ), POINTER :: A_flow      ! Flow parameter [Pa^-3 y^-1]
    REAL(dp), DIMENSION(:    ), POINTER :: A_flow_mean ! Vertically mean of the flow parameter [Pa^-3 y^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER :: Cpi         ! Specific heat capacity of ice [J kg^-1 K^-1].
    REAL(dp), DIMENSION(:,:  ), POINTER :: Ki          ! Conductivity of ice [J m^-1 K^-1 yr^-1].
    INTEGER                             :: wTi_pmp, wA_flow, wA_flow_mean, wCpi, wKi
    
    ! Spatial derivatives and curvatures
    REAL(dp), DIMENSION(:    ), POINTER :: dHi_dx, dHi_dy, dHi_dt
    REAL(dp), DIMENSION(:    ), POINTER :: dHb_dx, dHb_dy, dHb_dt
    REAL(dp), DIMENSION(:    ), POINTER :: dHs_dx, dHs_dy, dHs_dt
    REAL(dp), DIMENSION(:    ), POINTER :: dHs_dx_shelf, dHs_dy_shelf
    INTEGER                             :: wdHi_dx, wdHi_dy, wdHi_dt, wdHb_dx, wdHb_dy, wdHb_dt, wdHs_dx, wdHs_dy, wdHs_dt, wdHs_dx_shelf, wdHs_dy_shelf
    
    ! Zeta derivatives
    REAL(dp), DIMENSION(:,:  ), POINTER :: dzeta_dt
    REAL(dp), DIMENSION(:,:  ), POINTER :: dzeta_dx
    REAL(dp), DIMENSION(:,:  ), POINTER :: dzeta_dy
    REAL(dp), DIMENSION(:    ), POINTER :: dzeta_dz
    INTEGER                             :: wdzeta_dt, wdzeta_dx, wdzeta_dy, wdzeta_dz
        
    ! Ice dynamics - SIA
    REAL(dp), DIMENSION(:    ), POINTER :: Hi_ac
    REAL(dp), DIMENSION(:    ), POINTER :: Hb_ac
    REAL(dp), DIMENSION(:    ), POINTER :: Hs_ac
    REAL(dp), DIMENSION(:    ), POINTER :: sealevel_ac
    INTEGER,  DIMENSION(:    ), POINTER :: mask_ac
    INTEGER,  DIMENSION(:    ), POINTER :: mask_ice_ac
    INTEGER,  DIMENSION(:    ), POINTER :: mask_sheet_ac
    INTEGER,  DIMENSION(:    ), POINTER :: mask_shelf_ac
    INTEGER,  DIMENSION(:    ), POINTER :: mask_ocean_ac
    INTEGER,  DIMENSION(:    ), POINTER :: mask_groundingline_ac
    INTEGER,  DIMENSION(:    ), POINTER :: mask_calvingfront_ac
    REAL(dp), DIMENSION(:,:  ), POINTER :: Ti_ac
    REAL(dp), DIMENSION(:,:  ), POINTER :: Ti_pmp_ac      ! The pressure melting point temperature [K]
    REAL(dp), DIMENSION(:,:  ), POINTER :: A_flow_ac      ! Flow parameter [Pa^-3 y^-1]
    REAL(dp), DIMENSION(:    ), POINTER :: A_flow_mean_ac ! Vertically mean of the flow parameter [Pa^-3 y^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER :: Cpi_ac         ! Specific heat capacity of ice [J kg^-1 K^-1].
    REAL(dp), DIMENSION(:,:  ), POINTER :: Ki_ac          ! Conductivity of ice [J m^-1 K^-1 yr^-1].
    REAL(dp), DIMENSION(:    ), POINTER :: dHs_dp_ac
    REAL(dp), DIMENSION(:    ), POINTER :: dHs_do_ac
    REAL(dp), DIMENSION(:    ), POINTER :: dHs_dx_ac
    REAL(dp), DIMENSION(:    ), POINTER :: dHs_dy_ac
    REAL(dp), DIMENSION(:,:  ), POINTER :: D_SIA_3D_ac
    REAL(dp), DIMENSION(:    ), POINTER :: D_SIA_ac
    REAL(dp), DIMENSION(:    ), POINTER :: Up_ac
    REAL(dp), DIMENSION(:    ), POINTER :: Uo_ac
    REAL(dp), DIMENSION(:    ), POINTER :: Ux_ac
    REAL(dp), DIMENSION(:    ), POINTER :: Uy_ac
    REAL(dp), DIMENSION(:    ), POINTER :: U_SIA
    REAL(dp), DIMENSION(:    ), POINTER :: V_SIA
    INTEGER                             :: wHi_ac, wHb_ac, wHs_ac, wsealevel_ac, wTi_ac, wTi_pmp_ac, wA_flow_ac, wA_flow_mean_ac, wCpi_ac, wKi_ac, wdHs_dp_ac, wdHs_do_ac
    INTEGER                             :: wmask_ac, wmask_ice_ac, wmask_sheet_ac, wmask_shelf_ac, wmask_ocean_ac, wmask_groundingline_ac, wmask_calvingfront_ac
    INTEGER                             :: wdHs_dx_ac, wdHs_dy_ac, wD_SIA_3D_ac, wD_SIA_ac, wUp_ac, wUo_ac, wUx_ac, wUy_ac , wU_SIA, wV_SIA
    
    ! Ice dynamics - SSA
    REAL(dp), DIMENSION(:    ), POINTER :: U_SSA
    REAL(dp), DIMENSION(:    ), POINTER :: V_SSA
    REAL(dp), DIMENSION(:    ), POINTER :: tau_yield
    REAL(dp), DIMENSION(:    ), POINTER :: dUdx, dUdy, dVdx, dVdy
    REAL(dp), DIMENSION(:    ), POINTER :: nu
    REAL(dp), DIMENSION(:    ), POINTER :: beta_base
    REAL(dp), DIMENSION(:    ), POINTER :: rhs_x
    REAL(dp), DIMENSION(:    ), POINTER :: rhs_y
    REAL(dp), DIMENSION(:    ), POINTER :: eu_i
    REAL(dp), DIMENSION(:    ), POINTER :: ev_i
    REAL(dp), DIMENSION(:    ), POINTER :: Rx
    REAL(dp), DIMENSION(:    ), POINTER :: Ry
    INTEGER                             :: wU_SSA, wV_SSA, wtau_yield, wdUdx, wdUdy, wdVdx, wdVdy, wnu, wbeta_base, wrhs_x, wrhs_y, weu_i, wev_i, wRx, wRy
    
    ! Ice dynamics - ice thickness calculation
    REAL(dp), DIMENSION(:    ), POINTER :: Ux_SSA_ac
    REAL(dp), DIMENSION(:    ), POINTER :: Uy_SSA_ac
    REAL(dp), DIMENSION(:    ), POINTER :: Up_SSA_ac
    REAL(dp), DIMENSION(:    ), POINTER :: Uo_SSA_ac
    REAL(dp), DIMENSION(:,:  ), POINTER :: dVi_in
    REAL(dp), DIMENSION(:,:  ), POINTER :: dVi_out
    INTEGER                             :: wdVi_in, wdVi_out, wUx_SSA_ac, wUy_SSA_ac, wUp_SSA_ac, wUo_SSA_ac
    
    ! Thermodynamics
    REAL(dp), DIMENSION(:    ), POINTER :: frictional_heating
    REAL(dp), DIMENSION(:    ), POINTER :: Fr
    REAL(dp), DIMENSION(:,:  ), POINTER :: U_3D
    REAL(dp), DIMENSION(:,:  ), POINTER :: V_3D
    REAL(dp), DIMENSION(:,:  ), POINTER :: W_3D
    INTEGER                             :: wfrictional_heating, wFr, wU_3D, wV_3D, wW_3D
    
    ! Mesh adaptation data
    REAL(dp), DIMENSION(:    ), POINTER :: surf_curv
    REAL(dp), DIMENSION(:    ), POINTER :: log_velocity
    INTEGER                             :: wsurf_curv, wlog_velocity
          
  END TYPE type_ice_model

  TYPE type_mesh
    ! The unstructured triangular mesh.

    ! Basic meta properties
    CHARACTER(LEN=3)                    :: region_name                   ! NAM, EAS, GRL, ANT
    REAL(dp)                            :: lambda_M                      ! Oblique stereographic projection parameters
    REAL(dp)                            :: phi_M
    REAL(dp)                            :: alpha_stereographic
    REAL(dp),                   POINTER :: xmin                          ! X and Y range of the square covered by the mesh
    REAL(dp),                   POINTER :: xmax
    REAL(dp),                   POINTER :: ymin
    REAL(dp),                   POINTER :: ymax
    REAL(dp),                   POINTER :: tol_dist                      ! Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    INTEGER,                    POINTER :: nconmax                       ! Maximum allowed number of connections per vertex - highest observed value is 11
    INTEGER,                    POINTER :: nV_max                        ! Size of allocated memory for vertices
    INTEGER,                    POINTER :: nTri_max                      ! Size of allocated memory for triangles
    INTEGER,                    POINTER :: nV                            ! Number of vertices
    INTEGER,                    POINTER :: nTri                          ! Number of triangles
    INTEGER,                    POINTER :: perturb_dir                   ! Perturbation direction (0 = anticlockwise, 1 = clockwise)
    REAL(dp),                   POINTER :: alpha_min                     ! Sharpest inner angle allowed by Rupperts algorithm
    REAL(dp),                   POINTER :: dz_max_ice                    ! Maximum allowed vertical error in Rupperts algorithm over ice
    REAL(dp),                   POINTER :: resolution_min                ! Finest   resolution of the mesh (= sqrt(min(A)/pi))
    REAL(dp),                   POINTER :: resolution_max                ! Coarsest resolution of the mesh (= sqrt(max(A)/pi))
    INTEGER :: wxmin, wxmax, wymin, wymax, wtol_dist, wnconmax, wnV_max, wnTri_max, wnV, wnTri, wperturb_dir, walpha_min, wdz_max_ice, wresolution_min, wresolution_max

    ! Actual mesh data
    REAL(dp), DIMENSION(:,:  ), POINTER :: V                             ! The X and Y coordinates of all the vertices
    REAL(dp), DIMENSION(:    ), POINTER :: A                             ! The area of each vertex's Voronoi cell
    REAL(dp), DIMENSION(:,:  ), POINTER :: VorGC                         ! The geometric center of each vertex's Voronoi cell
    REAL(dp), DIMENSION(:    ), POINTER :: R                             ! The resolution (defined as distance to nearest neighbour)
    INTEGER,  DIMENSION(:    ), POINTER :: nC                            ! The number of other vertices this vertex is connected to
    INTEGER,  DIMENSION(:,:  ), POINTER :: C                             ! The list   of other vertices this vertex is connected to (ordered counter-clockwise, from edge to edge for edge vertices)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Cw                            ! The width of those connections (=length of the shared Voronoi cell edge), needed for ice fluxes
    INTEGER,  DIMENSION(:    ), POINTER :: niTri                         ! The number of triangles this vertex is a part of
    INTEGER,  DIMENSION(:,:  ), POINTER :: iTri                          ! The list   of triangles this vertex is a part of (ordered counter-clockwise)
    INTEGER,  DIMENSION(:    ), POINTER :: edge_index                    ! Each vertex's Edge Index; 0 = free, 1 = north, 2 = northeast, etc.
    INTEGER,  DIMENSION(:    ), POINTER :: mesh_old_ti_in                ! When creating a new mesh: for every new mesh vertex, the old mesh triangle that contains it
    INTEGER                             :: wV, wA, wVorGC, wR, wnC, wC, wCw, wniTri, wiTri, wedge_index, wmesh_old_ti_in ! MPI windows to all these memory spaces

    INTEGER,  DIMENSION(:,:  ), POINTER :: Tri                           ! The triangle array: Tri(ti) = [vi1, vi2, vi3] (vertices ordered counter-clockwise)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Tricc                         ! The X,Y-coordinates of each triangle's circumcenter
    INTEGER,  DIMENSION(:,:  ), POINTER :: TriC                          ! The (up to) three neighbour triangles (order across from 1st, 2nd and 3d vertex, respectively)
    INTEGER,  DIMENSION(:    ), POINTER :: Tri_edge_index                ! Same as for vertices, only NE, SE, etc. aren't used
    REAL(dp), DIMENSION(:    ), POINTER :: TriA                          ! The area of each triangle
    INTEGER                             :: wTri, wTricc, wTriC, wTri_edge_index, wTriA ! MPI windows to all these memory spaces
    
    INTEGER,  DIMENSION(:,:  ), POINTER :: Triflip                       ! List of triangles to flip, used in updating Delaunay triangulation
    INTEGER,  DIMENSION(:    ), POINTER :: RefMap                        ! Map   of which triangles      have been marked for refining
    INTEGER,  DIMENSION(:    ), POINTER :: RefStack                      ! Stack of       triangles that have been marked for refining
    INTEGER,                    POINTER :: RefStackN
    INTEGER                             :: wTriflip, wRefMap, wRefStack, wRefStackN ! MPI windows to all these memory spaces

    INTEGER,  DIMENSION(:    ), POINTER :: VMap, VStack1, VStack2        ! Map and stacks for doing a FloodFill-style search on the vertices.  To be used as distributed shared memory.
    INTEGER,  DIMENSION(:    ), POINTER :: TriMap, TriStack1, TriStack2  ! Map and stacks for doing a FloodFill-style search on the triangles. To be used as distributed shared memory.
    INTEGER                             :: wVMap, wVStack1, wVStack2, wTriMap, wTriStack1, wTriStack2
    INTEGER                             :: VStackN1, VStackN2
    INTEGER                             :: TriStackN1, TriStackN2

    REAL(dp), DIMENSION(:,:  ), POINTER :: NxTri                         ! The first order neighbour function for each triangle
    REAL(dp), DIMENSION(:,:  ), POINTER :: NyTri                         ! See subroutine "Calculate_Neighbour_Functions" for explanation    
    REAL(dp), DIMENSION(:,:  ), POINTER :: Nx                            ! The first order neighbour function for each vertex
    REAL(dp), DIMENSION(:,:  ), POINTER :: Ny
    REAL(dp), DIMENSION(:,:  ), POINTER :: Nxx                           ! The second order neighbour functions for each vertex
    REAL(dp), DIMENSION(:,:  ), POINTER :: Nxy
    REAL(dp), DIMENSION(:,:  ), POINTER :: Nyy    
    REAL(dp), DIMENSION(:,:  ), POINTER :: Nxm                           ! The "masked" first order neighbour function for each vertex
    REAL(dp), DIMENSION(:,:  ), POINTER :: Nym
    REAL(dp), DIMENSION(:,:  ), POINTER :: Nxxm                          ! The "masked" second order neighbour functions for each vertex
    REAL(dp), DIMENSION(:,:  ), POINTER :: Nxym
    REAL(dp), DIMENSION(:,:  ), POINTER :: Nyym
    INTEGER                             :: wNxTri, wNyTri, wNx, wNy, wNxx, wNxy, wNyy, wNxm, wNym, wNxxm, wNxym, wNyym ! MPI windows to all these memory spaces
    
    INTEGER,                    POINTER :: nAC                           ! The number of Ac vertices
    REAL(dp), DIMENSION(:,:  ), POINTER :: VAc                           ! x,y coordinates of the Ac vertices
    INTEGER,  DIMENSION(:,:  ), POINTER :: Aci                           ! Mapping array from the Aa to the Ac mesh
    INTEGER,  DIMENSION(:,:  ), POINTER :: iAci                          ! Mapping array from the Ac to the Aa mesh
    REAL(dp), DIMENSION(:,:),   POINTER :: Nx_ac                         ! x          neighbour functions on the Ac mesh
    REAL(dp), DIMENSION(:,:),   POINTER :: Ny_ac                         ! y          neighbour functions on the Ac mesh
    REAL(dp), DIMENSION(:  ),   POINTER :: Np_ac                         ! Parallel   neighbour functions on the Ac mesh
    REAL(dp), DIMENSION(:,:),   POINTER :: No_ac                         ! Orthogonal neighbour functions on the Ac mesh
    INTEGER                             :: wnAC, wVAc, wAci, wiAci, wNx_ac, wNy_ac, wNp_ac, wNo_ac
    
    REAL(dp), DIMENSION(:    ), POINTER :: lat                           ! Lat, lon coordinates
    REAL(dp), DIMENSION(:    ), POINTER :: lon
    INTEGER                             :: wlat, wlon
    
    REAL(dp), DIMENSION(:    ), POINTER :: omega                         ! Relaxation parameter used in SOR solving of SSA
    INTEGER                             :: womega
    
    INTEGER,                    POINTER :: nPOI                          ! Number of Points of Interest (POI) in this mesh
    REAL(dp), DIMENSION(:,:  ), POINTER :: POI_coordinates               ! Lat-lon coordinates of a POI
    REAL(dp), DIMENSION(:,:  ), POINTER :: POI_XY_coordinates            ! X-Y     coordinates of a POI
    REAL(dp), DIMENSION(:    ), POINTER :: POI_resolutions               ! Resolution          of a POI (km)
    INTEGER,  DIMENSION(:,:  ), POINTER :: POI_vi                        ! The three vertices surrounding a POI
    REAL(dp), DIMENSION(:,:  ), POINTER :: POI_w                         ! Their relative weights in trilinear interpolation
    INTEGER                             :: wnPOI, wPOI_coordinates, wPOI_XY_coordinates, wPOI_resolutions, wPOI_vi, wPOI_w
    
    ! Parallelisation
    INTEGER                             :: v1, v2                        ! Vertex domain [v1,v2] of each process (equal number of vertices)
    INTEGER                             :: t1, t2                        ! Vertex domain [v1,v2] of each process (equal number of triangles)
    INTEGER                             :: c1, c2                        ! Arakawa C vertex domain [c1,c2] of each process

  END TYPE type_mesh
  
  TYPE type_mesh_config
    ! Config containing the different mesh resolutions, so they can be changed for the different nested meshes for the multigrid SSA solver
       
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
  
  END TYPE type_mesh_config
  
  TYPE type_remapping_trilin
    ! Indices and weights for mapping data between two meshes, using trilinear interpolation
    
    INTEGER,  DIMENSION(:,:  ), POINTER :: vi   ! Indices of vertices from mesh_src contributing to this mesh_dst vertex
    REAL(dp), DIMENSION(:,:  ), POINTER :: w    ! Weights of vertices from mesh_src contributing to this mesh_dst vertex
    
    INTEGER :: wvi, ww
    
  END TYPE type_remapping_trilin
  
  TYPE type_remapping_nearest_neighbour
    ! Indices and weights for mapping data between two meshes, using nearest-neighbour interpolation
    
    INTEGER,  DIMENSION(:    ), POINTER :: vi   ! Index of mesh_src vertex nearest to this mesh_dst vertex
    
    INTEGER :: wvi
    
  END TYPE type_remapping_nearest_neighbour
  
  TYPE type_remapping_conservative
    ! Indices and weights for mapping data between two meshes, using 1st or 2nd order conservative remapping
    
    INTEGER,  DIMENSION(:    ), POINTER :: nV   ! Number  of vertices from mesh_src contributing to this mesh_dst vertex
    INTEGER,  DIMENSION(:,:  ), POINTER :: vi   ! Indices of vertices from mesh_src contributing to this mesh_dst vertex
    REAL(dp), DIMENSION(:,:  ), POINTER :: w0   ! Weights of vertices from mesh_src contributing to this mesh_dst vertex
    REAL(dp), DIMENSION(:,:  ), POINTER :: w1x  ! Weights of vertices from mesh_src contributing to this mesh_dst vertex
    REAL(dp), DIMENSION(:,:  ), POINTER :: w1y  ! Weights of vertices from mesh_src contributing to this mesh_dst vertex
    INTEGER :: wnV, wvi, ww0, ww1x, ww1y
    
  END TYPE type_remapping_conservative
  
  TYPE type_remapping
    ! Mapping index and weight arrays for remapping data between two meshes, for different remapping methods
    
    TYPE(type_remapping_trilin)                 :: trilin
    TYPE(type_remapping_nearest_neighbour)      :: nearest_neighbour
    TYPE(type_remapping_conservative)           :: conservative
  
  END TYPE type_remapping
    
  TYPE type_netcdf_file_output
    ! Integers describing open ports to different variables in an opened NetCDF file.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! Index of time frame to be written to
    INTEGER :: ti
  
    ! Number of time frames
    CHARACTER(LEN=256) :: name_nt                    = 'nt                   '
    
    ! ID's for variables:
    ! ===================
    
    ! Dimension variables (x and y not used for output, but for reading Cartesian input data)
    INTEGER :: var_time
    INTEGER :: var_lat
    INTEGER :: var_lon
    INTEGER :: var_zeta
    INTEGER :: var_month
  
    CHARACTER(LEN=256) :: name_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_month                 = 'month                '
    CHARACTER(LEN=256) :: name_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_lon                   = 'lon                  '
    
    ! Mesh data - V and Tri always saved, rest only when specified through config
    INTEGER :: var_mesh_V
    INTEGER :: var_mesh_Tri
        
    INTEGER :: var_mesh_VAc
    INTEGER :: var_mesh_Aci
    INTEGER :: var_mesh_iAci

    INTEGER :: var_mesh_A
    INTEGER :: var_mesh_R
    INTEGER :: var_mesh_nC
    INTEGER :: var_mesh_C
    INTEGER :: var_mesh_Cw
    INTEGER :: var_mesh_niTri
    INTEGER :: var_mesh_iTri
    INTEGER :: var_mesh_edge_index
    INTEGER :: var_mesh_Tricc
    INTEGER :: var_mesh_TriC
    INTEGER :: var_mesh_Tri_edge_index
    INTEGER :: var_mesh_TriA
    INTEGER :: var_mesh_NxTri
    INTEGER :: var_mesh_NyTri
    INTEGER :: var_mesh_Nx
    INTEGER :: var_mesh_Ny
    INTEGER :: var_mesh_Nxx
    INTEGER :: var_mesh_Nxy
    INTEGER :: var_mesh_Nyy
  
    CHARACTER(LEN=256) :: name_mesh_nV               = 'nV                   '
    CHARACTER(LEN=256) :: name_mesh_nTri             = 'nTri                 '
    CHARACTER(LEN=256) :: name_mesh_V                = 'vertex               '
    CHARACTER(LEN=256) :: name_mesh_V2               = 'vertex2              '
    CHARACTER(LEN=256) :: name_mesh_Tri              = 'triangle             '
    CHARACTER(LEN=256) :: name_mesh_Tri2             = 'triangle2            '
    CHARACTER(LEN=256) :: name_mesh_AcV              = 'AcV                  '
    CHARACTER(LEN=256) :: name_mesh_VAc              = 'VAc                  '
    CHARACTER(LEN=256) :: name_mesh_Aci              = 'Aci                  '
    CHARACTER(LEN=256) :: name_mesh_iAci             = 'iAci                 '
    CHARACTER(LEN=256) :: name_mesh_A                = 'A                    '
    CHARACTER(LEN=256) :: name_mesh_R                = 'R                    '
    CHARACTER(LEN=256) :: name_mesh_nC               = 'nC                   '
    CHARACTER(LEN=256) :: name_mesh_C                = 'C                    '
    CHARACTER(LEN=256) :: name_mesh_Cw               = 'Cw                   '
    CHARACTER(LEN=256) :: name_mesh_niTri            = 'niTri                '
    CHARACTER(LEN=256) :: name_mesh_iTri             = 'iTri                 '
    CHARACTER(LEN=256) :: name_mesh_edge_index       = 'edge_index           '
    CHARACTER(LEN=256) :: name_mesh_Tricc            = 'Tricc                '
    CHARACTER(LEN=256) :: name_mesh_TriC             = 'TriC                 '
    CHARACTER(LEN=256) :: name_mesh_Tri_edge_index   = 'Tri_edge_index       '
    CHARACTER(LEN=256) :: name_mesh_TriA             = 'TriA                 '
    CHARACTER(LEN=256) :: name_mesh_NxTri            = 'NxTri                '
    CHARACTER(LEN=256) :: name_mesh_NyTri            = 'NyTri                '
    CHARACTER(LEN=256) :: name_mesh_Nx               = 'Nx                   '
    CHARACTER(LEN=256) :: name_mesh_Ny               = 'Ny                   '
    CHARACTER(LEN=256) :: name_mesh_Nxx              = 'Nxx                  '
    CHARACTER(LEN=256) :: name_mesh_Nxy              = 'Nxy                  '
    CHARACTER(LEN=256) :: name_mesh_Nyy              = 'Nyy                  '
    
    ! Key output
    INTEGER :: var_Hi
    INTEGER :: var_Hb
    INTEGER :: var_Hs
    INTEGER :: var_Us
    INTEGER :: var_Vs
    INTEGER :: var_U_SSA
    INTEGER :: var_V_SSA
    INTEGER :: var_Ti
    INTEGER :: var_mask
    INTEGER :: var_T2m
    INTEGER :: var_Precip
    INTEGER :: var_SMB
    INTEGER :: var_BMB
    INTEGER :: var_FirnDepth
    INTEGER :: var_MeltPreviousYear
    
    CHARACTER(LEN=256) :: name_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_Hs                    = 'Hs                   ' 
    CHARACTER(LEN=256) :: name_Us                    = 'Us                   '
    CHARACTER(LEN=256) :: name_Vs                    = 'Vs                   '
    CHARACTER(LEN=256) :: name_U_SSA                 = 'U_SSA                '
    CHARACTER(LEN=256) :: name_V_SSA                 = 'V_SSA                '
    CHARACTER(LEN=256) :: name_Ti                    = 'Ti                   '
    CHARACTER(LEN=256) :: name_mask                  = 'mask                 '
    CHARACTER(LEN=256) :: name_T2m                   = 'T2m                  '
    CHARACTER(LEN=256) :: name_Precip                = 'Precip               '
    CHARACTER(LEN=256) :: name_SMB                   = 'SMB                  '
    CHARACTER(LEN=256) :: name_BMB                   = 'BMB                  '
    CHARACTER(LEN=256) :: name_FirnDepth             = 'FirnDepth            '
    CHARACTER(LEN=256) :: name_MeltPreviousYear      = 'MeltPreviousYear     '
    
    ! Climate
    INTEGER :: var_Wind_WE
    INTEGER :: var_Wind_SN
    
    CHARACTER(LEN=256) :: name_Wind_WE               = 'Wind_WE              '
    CHARACTER(LEN=256) :: name_Wind_SN               = 'Wind_SN              '
       
    ! SMB
    INTEGER :: var_Q_TOA
    INTEGER :: var_Snowfall
    INTEGER :: var_Melt
    INTEGER :: var_Refreezing
    INTEGER :: var_Runoff
    INTEGER :: var_Albedo
    
    CHARACTER(LEN=256) :: name_Q_TOA                 = 'Q_TOA                '
    CHARACTER(LEN=256) :: name_Snowfall              = 'Snowfall             '
    CHARACTER(LEN=256) :: name_Melt                  = 'Melt                 '
    CHARACTER(LEN=256) :: name_Refreezing            = 'Refreezing           '
    CHARACTER(LEN=256) :: name_Runoff                = 'Runoff               '
    CHARACTER(LEN=256) :: name_Albedo                = 'Albedo               '
    
    ! Ice dynamics
    INTEGER :: var_dHs_dx
    INTEGER :: var_dHs_dy    
    INTEGER :: var_U_SIA
    INTEGER :: var_V_SIA  
    INTEGER :: var_Ux_ac
    INTEGER :: var_Uy_ac
    INTEGER :: var_U_3D
    INTEGER :: var_V_3D
    INTEGER :: var_W_3D
    
    CHARACTER(LEN=256) :: name_dHs_dx                = 'dHs_dx               '
    CHARACTER(LEN=256) :: name_dHs_dy                = 'dHs_dy               '
    CHARACTER(LEN=256) :: name_Ux_ac                 = 'Ux_ac                '
    CHARACTER(LEN=256) :: name_Uy_ac                 = 'Uy_ac                '
    CHARACTER(LEN=256) :: name_U_SIA                 = 'U_SIA                '
    CHARACTER(LEN=256) :: name_V_SIA                 = 'V_SIA                '
    CHARACTER(LEN=256) :: name_U_3D                  = 'U_3D                 '
    CHARACTER(LEN=256) :: name_V_3D                  = 'V_3D                 '
    CHARACTER(LEN=256) :: name_W_3D                  = 'W_3D                 '    
        
  END TYPE type_netcdf_file_output
    
  TYPE type_netcdf_file_PD_data
    ! Integers describing open ports to different variables in an opened NetCDF file.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    ! Dimension variables (x and y not used for output, but for reading Cartesian input data)
    INTEGER :: var_x
    INTEGER :: var_y
    INTEGER :: var_month
    
    ! Variable names
    CHARACTER(LEN=256) :: name_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_month                 = 'month                '
    
    ! Key output
    INTEGER :: var_Hi
    INTEGER :: var_Hb
    INTEGER :: var_Hs
    
    ! Variable names
    CHARACTER(LEN=256) :: name_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_Hs                    = 'Hs                   ' 
        
  END TYPE type_netcdf_file_PD_data
    
  TYPE type_netcdf_file_init_PD_data
    ! Integers describing open ports to different variables in an opened NetCDF file.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    ! Dimension variables (x and y not used for output, but for reading Cartesian input data)
    INTEGER :: var_x
    INTEGER :: var_y
    
    ! Variable names
    CHARACTER(LEN=256) :: name_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_y                     = 'y                    '
    
    ! Key output
    INTEGER :: var_Hi
    INTEGER :: var_Hb
    INTEGER :: var_Hs
    INTEGER :: var_firn_depth
    
    ! Variable names
    CHARACTER(LEN=256) :: name_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_Hs                    = 'Hs                   '
    CHARACTER(LEN=256) :: name_Firn_depth            = 'Firn_depth           '
        
  END TYPE type_netcdf_file_init_PD_data
    
  TYPE type_netcdf_file_climate_data
    ! Integers describing open ports to different variables in an opened NetCDF file.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    ! Dimension variables (x and y not used for output, but for reading Cartesian input data)
    INTEGER :: var_lon
    INTEGER :: var_lat
    INTEGER :: var_month
    
    ! Variable names
    CHARACTER(LEN=256) :: name_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_month                 = 'month                '
    
    ! Key output
    INTEGER :: var_T2m
    INTEGER :: var_Precip
    INTEGER :: var_Hs
    INTEGER :: var_Albedo
    INTEGER :: var_Wind_WE
    INTEGER :: var_Wind_SN
    
    ! Variable names
    CHARACTER(LEN=256) :: name_T2m                   = 'T2m                  '
    CHARACTER(LEN=256) :: name_Precip                = 'Precip               '
    CHARACTER(LEN=256) :: name_Hs                    = 'Hs                   '
    CHARACTER(LEN=256) :: name_Albedo                = 'Albedo               '
    CHARACTER(LEN=256) :: name_Wind_WE               = 'Wind_WE              '
    CHARACTER(LEN=256) :: name_Wind_SN               = 'Wind_SN              '
        
  END TYPE type_netcdf_file_climate_data
    
  TYPE type_netcdf_file_insolation
    ! Integers describing open ports to different variables in an opened NetCDF file.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    ! Dimension variables (x and y not used for output, but for reading Cartesian input data)
    INTEGER :: var_time
    INTEGER :: var_month
    INTEGER :: var_lat
    
    ! Variable names
    CHARACTER(LEN=256) :: name_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_month                 = 'month                '
    CHARACTER(LEN=256) :: name_lat                   = 'lat                  '
    
    ! Key output
    INTEGER :: var_Q_TOA
    
    ! Variable names
    CHARACTER(LEN=256) :: name_Q_TOA                 = 'Q_TOA                '
        
  END TYPE type_netcdf_file_insolation
  
  TYPE type_subclimate_global
    ! Global climate data, either from present-day observations (ERA40) or from a GCM snapshot.
    
    CHARACTER(LEN=256)                  :: name             ! 'ERA40', 'HadCM3_PI', etc.
    
    ! NetCDF file containing the data
    TYPE(type_netcdf_file_climate_data) :: netcdf
    
    ! Grid
    INTEGER,                    POINTER :: nlat, nlon
    REAL(dp), DIMENSION(:    ), POINTER :: lat
    REAL(dp), DIMENSION(:    ), POINTER :: lon
    INTEGER                             :: wnlat, wnlon, wlat, wlon
        
    ! General forcing info (not relevant for PD observations)
    REAL(dp),                   POINTER :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp),                   POINTER :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   POINTER :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp),                   POINTER :: orbit_obl
    REAL(dp),                   POINTER :: orbit_pre
    INTEGER                             :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre ! MPI windows to all these memory spaces
    
    ! Actual GCM data
    REAL(dp), DIMENSION(:,:,:), POINTER :: T2m                           ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:,:), POINTER :: Precip                        ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Hs                            ! Orography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:,:), POINTER :: Q_TOA                         ! Monthly mean insolation at the top of the atmosphere (W/m2)
    REAL(dp), DIMENSION(:,:,:), POINTER :: Albedo                        ! Monthly mean surface albedo (calculated using UFEMISM's own mass balance scheme for GCM snapshots)
    REAL(dp), DIMENSION(:,:,:), POINTER :: Wind_WE                       ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER :: Wind_SN                       ! Monthly mean south_north wind speed (m/s)
    INTEGER                             :: wT2m, wPrecip, wHs, wQ_TOA, wAlbedo, wWind_WE, wWind_SN ! MPI windows to all these memory spaces
    
    ! Paralelisation
    INTEGER                             :: i1, i2                      ! Grid domain (i1:i2,:) of each process
  
  END TYPE type_subclimate_global
  
  TYPE type_climate_matrix
    ! The climate matrix data structure. Contains all the different global GCM snapshots.
    
    ! The present-day observed climate (ERA40)
    TYPE(type_subclimate_global)        :: PD_obs
    
    ! The GCM snapshots.
    TYPE(type_subclimate_global)        :: PI
    TYPE(type_subclimate_global)        :: LGM  
  
  END TYPE type_climate_matrix
  
  TYPE type_subclimate_mesh
    ! Data from PD observations or a GCM snapshot, projected from their initial global grid onto the mesh.
        
    ! General forcing info (not relevant for PD observations)
    REAL(dp),                   POINTER :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp),                   POINTER :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   POINTER :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp),                   POINTER :: orbit_obl
    REAL(dp),                   POINTER :: orbit_pre
    INTEGER                             :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre ! MPI windows to all these memory spaces
    
    ! Actual GCM data
    REAL(dp), DIMENSION(:,:  ), POINTER :: T2m                           ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Precip                        ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:    ), POINTER :: Hs                            ! Orography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Q_TOA                         ! Monthly mean insolation at the top of the atmosphere (W/m2)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Albedo                        ! Monthly mean surface albedo (calculated using UFEMISM's own mass balance scheme for GCM snapshots)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Wind_WE                       ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Wind_SN                       ! Monthly mean south_north wind speed (m/s)
    INTEGER                             :: wT2m, wPrecip, wHs, wQ_TOA, wAlbedo, wWind_WE, wWind_SN ! MPI windows to all these memory spaces
  
  END TYPE type_subclimate_mesh
  
  TYPE type_climate_model
    ! All the relevant climate data fields (PD observations, GCM snapshots, and final, applied climate) on the mesh
    
    TYPE(type_subclimate_mesh)          :: PD_obs
    TYPE(type_subclimate_mesh)          :: GCM_PI
    TYPE(type_subclimate_mesh)          :: GCM_LGM
    TYPE(type_subclimate_mesh)          :: applied
          
  END TYPE type_climate_model
  
  TYPE type_SMB_model
    ! The different SMB components, calculated from the prescribed climate
    
    REAL(dp), DIMENSION(:,:  ), POINTER :: Q_TOA                         ! The prescribed monthly insolation, from an external forcing file
    REAL(dp), DIMENSION(:    ), POINTER :: AlbedoSurf                    ! Surface albedo underneath the snow layer (water, rock or ice)
    REAL(dp), DIMENSION(:    ), POINTER :: MeltPreviousYear              ! Total melt that occurred during the previous year (m)
    REAL(dp), DIMENSION(:,:  ), POINTER :: FirnDepth                     ! Depth of the firn layer (m)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Rainfall                      ! Monthly rainfall (m)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Snowfall                      ! Monthly snowfall (m)
    REAL(dp), DIMENSION(:,:  ), POINTER :: AddedFirn                     ! Monthly added firn (m)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Melt                          ! Monthly melt (m)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Refreezing                    ! Monthly refreezing (m)
    REAL(dp), DIMENSION(:    ), POINTER :: Refreezing_year               ! Yearly  refreezing (m)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Runoff                        ! Monthly runoff (m)
    REAL(dp), DIMENSION(:,:  ), POINTER :: Albedo                        ! Monthly albedo
    REAL(dp), DIMENSION(:,:  ), POINTER :: SMB                           ! Monthly SMB (m)
    REAL(dp), DIMENSION(:    ), POINTER :: SMB_year                      ! Yearly  SMB (m)
    INTEGER :: wQ_TOA, wAlbedoSurf, wMeltPreviousYear, wFirnDepth, wRainfall, wSnowfall
    INTEGER :: wAddedFirn, wMelt, wRefreezing, wRefreezing_year, wRunoff, wAlbedo, wSMB, wSMB_year
  
  END TYPE type_SMB_model
  
  TYPE type_BMB_model
    ! The different BMB components
    
    REAL(dp), DIMENSION(:    ), POINTER :: BMB                           ! The basal mass balance (same as SMB: negative means melt)
    REAL(dp), DIMENSION(:    ), POINTER :: BMB_sheet                     ! The basal mass balance underneath the land-based ice sheet
    REAL(dp), DIMENSION(:    ), POINTER :: BMB_shelf                     ! The basal mass balance underneath the floating   ice shelf
    REAL(dp), DIMENSION(:    ), POINTER :: sub_angle                     ! "subtended angle"      for the sub-shelf melt parameterisation
    REAL(dp), DIMENSION(:    ), POINTER :: dist_open                     ! distance to open ocean for the sub-shelf melt parameterisation
    REAL(dp), DIMENSION(:    ), POINTER :: ocean_intrusion               ! "ocean intrusion" coefficient for the new sub-shelf melt parameterisation
    INTEGER :: wBMB, wBMB_sheet, wBMB_shelf, wsub_angle, wdist_open, wocean_intrusion
  
  END TYPE type_BMB_model
  
  TYPE type_PD_data_fields
    ! Data structure containing data fields describing the present-day world.
    
    ! NetCDF file containing the data
    TYPE(type_netcdf_file_PD_data)      :: netcdf
    
    ! Grid
    INTEGER,                    POINTER :: nx, ny
    REAL(dp), DIMENSION(:    ), POINTER :: x
    REAL(dp), DIMENSION(:    ), POINTER :: y
    INTEGER                             :: wnx, wny, wx, wy
    
    ! On Cartesian grid
    REAL(dp), DIMENSION(:,:  ), POINTER :: Hi_cart
    REAL(dp), DIMENSION(:,:  ), POINTER :: Hb_cart
    REAL(dp), DIMENSION(:,:  ), POINTER :: Hs_cart
    INTEGER,  DIMENSION(:,:  ), POINTER :: mask_cart
    INTEGER                             :: wHi_cart, wHb_cart, wHs_cart, wmask_cart
    
    ! On the mesh
    REAL(dp), DIMENSION(:    ), POINTER :: Hi
    REAL(dp), DIMENSION(:    ), POINTER :: Hb
    REAL(dp), DIMENSION(:    ), POINTER :: Hs
    INTEGER,  DIMENSION(:    ), POINTER :: mask
    INTEGER                             :: wHi, wHb, wHs, wmask
    
    ! Paralelisation
    INTEGER                             :: i1, i2                      ! Grid domain (i1:i2,:) of each process
              
  END TYPE type_PD_data_fields
  
  TYPE type_init_data_fields
    ! Data structure containing data fields describing the present-day world.
    
    ! NetCDF file containing the data
    TYPE(type_netcdf_file_init_PD_data) :: netcdf
    
    INTEGER,                    POINTER           :: nx, ny
    REAL(dp), DIMENSION(:    ), POINTER :: x
    REAL(dp), DIMENSION(:    ), POINTER :: y
    INTEGER                             :: wnx, wny, wx, wy
    
    ! On Cartesian grid
    REAL(dp), DIMENSION(:,:  ), POINTER :: Hi_cart
    REAL(dp), DIMENSION(:,:  ), POINTER :: Hb_cart
    REAL(dp), DIMENSION(:,:  ), POINTER :: Hs_cart
    REAL(dp), DIMENSION(:,:  ), POINTER :: Firn_depth_cart
    REAL(dp), DIMENSION(:,:  ), POINTER :: surface_curvature_cart
    INTEGER,  DIMENSION(:,:  ), POINTER :: mask_cart
    INTEGER,  DIMENSION(:,:  ), POINTER :: mask_ice_cart
    INTEGER,  DIMENSION(:,:  ), POINTER :: mask_gl_cart
    INTEGER,  DIMENSION(:,:  ), POINTER :: mask_cf_cart
    INTEGER,  DIMENSION(:,:  ), POINTER :: mask_coast_cart
    INTEGER                             :: wHi_cart, wHb_cart, wHs_cart, wFirn_depth_cart, wsurface_curvature_cart
    INTEGER                             :: wmask_cart, wmask_ice_cart, wmask_gl_cart, wmask_cf_cart, wmask_coast_cart
    
    ! On the mesh
    REAL(dp), DIMENSION(:    ), POINTER :: Hi
    REAL(dp), DIMENSION(:    ), POINTER :: Hb
    REAL(dp), DIMENSION(:    ), POINTER :: Hs
    INTEGER,  DIMENSION(:    ), POINTER :: mask
    REAL(dp), DIMENSION(:    ), POINTER :: U_SSA
    REAL(dp), DIMENSION(:    ), POINTER :: V_SSA
    REAL(dp), DIMENSION(:,:  ), POINTER :: Ti
    REAL(dp), DIMENSION(:,:  ), POINTER :: Albedo
    REAL(dp), DIMENSION(:,:  ), POINTER :: FirnDepth
    REAL(dp), DIMENSION(:    ), POINTER :: MeltPreviousYear
    INTEGER                             :: wHi, wHb, wHs, wmask, wU_SSA, wV_SSA, wTi, wAlbedo, wFirnDepth, wMeltPReviousYear
    
    ! Paralelisation
    INTEGER                             :: i1, i2                      ! Grid domain (i1:i2,:) of each process
              
  END TYPE type_init_data_fields
  
  TYPE type_forcing_data
    ! Data structure containing model forcing data - CO2 record, d18O record, (global) insolation record
    
    ! Insolation forcing (Laskar et al., 2004)
    TYPE(type_netcdf_file_insolation)   :: ins_netcdf
    INTEGER                             :: ins_nyears
    INTEGER                             :: ins_nlat
    REAL(dp), DIMENSION(:    ), POINTER :: ins_time
    REAL(dp), DIMENSION(:    ), POINTER :: ins_lat
    REAL(dp),                   POINTER :: ins_t0, ins_t1
    REAL(dp), DIMENSION(:,:  ), POINTER :: ins_Q_TOA0, ins_Q_TOA1
    INTEGER                             :: wins_time, wins_lat, wins_t0, wins_t1, wins_Q_TOA0, wins_Q_TOA1
    
  END TYPE type_forcing_data
  
  TYPE type_model_region
    ! Contains all the different data structures, organised by sub-model (ice, climate)
    
    ! Metadata
    CHARACTER(LEN=3)                    :: name                  ! NAM, EAS, GRL, ANT
    CHARACTER(LEN=256)                  :: long_name             ! North America, Eurasia, Greenland, Antarctica
    
    ! Timers and switches for determining which modules need to be called at what points in time during the simulation
    REAL(dp)                            :: time                   ! The current time of this particular region.
    REAL(dp)                            :: dt, dt_prev            ! The (dynamic) time step    
    REAL(dp)                            :: dt_SIA                 ! The critical timestep following from the SIA ice velocities
    REAL(dp)                            :: dt_SSA                 ! The critical timestep following from the SSA ice velocities    
    REAL(dp)                            :: t0_SIA,     t1_SIA     ! Time of last and next update of SIA velocity
    REAL(dp)                            :: t0_SSA,     t1_SSA     !                                 SSA velocity
    REAL(dp)                            :: t0_thermo,  t1_Thermo  !                                 thermodynamics
    REAL(dp)                            :: t0_output,  t1_output  !                                 output file
    REAL(dp)                            :: t0_mesh,    t1_mesh    !                                 mesh
    REAL(dp)                            :: t0_climate, t1_climate !                                 climate
    REAL(dp)                            :: t0_SMB,     t1_SMB     !                                 SMB
    REAL(dp)                            :: t0_BMB,     t1_BMB     !                                 BMB
    
    LOGICAL                             :: do_solve_SIA
    LOGICAL                             :: do_solve_SSA
    LOGICAL                             :: do_thermodynamics
    LOGICAL                             :: do_climate
    LOGICAL                             :: do_SMB
    LOGICAL                             :: do_BMB
    LOGICAL                             :: do_write_output
    
    ! Keep track of whether or not a NetCDF output file has been created (since we usually update the mesh more often than the output)
    LOGICAL                             :: output_file_exists
    
    ! The region's ice sheet's volume and volume above flotation (in mSLE, so the second one is the ice sheets GMSL contribution)
    REAL(dp), POINTER                   :: ice_area
    REAL(dp), POINTER                   :: ice_volume
    REAL(dp), POINTER                   :: ice_volume_above_flotation
    REAL(dp), POINTER                   :: ice_volume_above_flotation_PD
    REAL(dp), POINTER                   :: GMSL_contribution
    INTEGER                             :: wice_area, wice_volume, wice_volume_above_flotation, wice_volume_above_flotation_PD, wGMSL_contribution
        
    ! Reference data fields
    TYPE(type_PD_data_fields)           :: PD               ! The present-day data fields for this model region, on a high-res Cartesian grid
    TYPE(type_init_data_fields)         :: init             ! The initial     data fields for this model region, on a high-res Cartesian grid
        
    ! Sub-models
    TYPE(type_mesh)                     :: mesh             ! The finite element mesh for this model region
    TYPE(type_mesh)                     :: mesh_new         ! The new mesh after updating (so that the old one can be kept until data has been mapped)
    TYPE(type_ice_model)                :: ice              ! All the ice model data for this model region
    TYPE(type_climate_model)            :: climate          ! All the climate data for this model region
    TYPE(type_SMB_model)                :: SMB              ! The different SMB components for this model region
    TYPE(type_BMB_model)                :: BMB              ! The different BMB components for this model region
    
    ! Output netcdf file
    TYPE(type_netcdf_file_output)       :: output
    
    ! Computation times
    REAL(dp)                            :: tcomp_total
    REAL(dp)                            :: tcomp_mesh
    REAL(dp)                            :: tcomp_SIA
    REAL(dp)                            :: tcomp_SSA
    REAL(dp)                            :: tcomp_thermo
    REAL(dp)                            :: tcomp_climate
    REAL(dp)                            :: tcomp_output
    
!    ! Low-resolution "nested" meshes used in the multigrid SSA solver
!    INTEGER                             :: n_nested_meshes  ! Not all of them will always be used    
!    TYPE(type_mesh)                     :: mesh_low_1       ! 200 km
!    TYPE(type_mesh)                     :: mesh_low_2       ! 100 km
!    TYPE(type_mesh)                     :: mesh_low_3       !  50 km
!    TYPE(type_mesh)                     :: mesh_low_4       !  25 km
!    TYPE(type_mesh)                     :: mesh_low_5       !  12 km
!    TYPE(type_mesh)                     :: mesh_low_6       !   6 km
!    TYPE(type_mesh)                     :: mesh_low_7       !   3 km
!    
!    ! Mapping arrays for mapping the model data and intermediate solutions between the nested meshes
!    TYPE(type_remapping)                :: map_mesh_12,         map_mesh_21
!    TYPE(type_remapping)                :: map_mesh_23,         map_mesh_32
!    TYPE(type_remapping)                :: map_mesh_34,         map_mesh_43
!    TYPE(type_remapping)                :: map_mesh_45,         map_mesh_54
!    TYPE(type_remapping)                :: map_mesh_56,         map_mesh_65
!    TYPE(type_remapping)                :: map_mesh_67,         map_mesh_76
!    TYPE(type_remapping)                :: map_mesh_nest_final, map_mesh_final_nest
!    
!    ! Low resolution ice model data used in the multigrid SSA solver
!    TYPE(type_ice_model)                :: ice_low_1
!    TYPE(type_ice_model)                :: ice_low_2
!    TYPE(type_ice_model)                :: ice_low_3
!    TYPE(type_ice_model)                :: ice_low_4
!    TYPE(type_ice_model)                :: ice_low_5
!    TYPE(type_ice_model)                :: ice_low_6
!    TYPE(type_ice_model)                :: ice_low_7
    
  END TYPE type_model_region
  
CONTAINS

END MODULE data_types_module
