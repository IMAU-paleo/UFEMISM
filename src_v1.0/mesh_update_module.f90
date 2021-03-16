MODULE mesh_update_module
  ! Routines for creating a new mesh based on forcing data on an old mesh.

  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync, allocate_shared_memory_int_0D, allocate_shared_memory_dp_0D, &
                                                    allocate_shared_memory_int_1D, allocate_shared_memory_dp_1D, &
                                                    allocate_shared_memory_int_2D, allocate_shared_memory_dp_2D, &
                                                    allocate_shared_memory_int_3D, allocate_shared_memory_dp_3D, &
                                                    allocate_shared_memory_bool_1D, deallocate_shared_memory, &
                                                    adapt_shared_memory_int_1D,    adapt_shared_memory_dp_1D, &
                                                    adapt_shared_memory_int_2D,    adapt_shared_memory_dp_2D, &
                                                    adapt_shared_memory_int_3D,    adapt_shared_memory_dp_3D, &
                                                    adapt_shared_memory_bool_1D                   
  USE data_types_module,           ONLY: type_mesh, type_mesh_config, type_model_region, type_ice_model, type_PD_data_fields
  USE mesh_help_functions_module,  ONLY: FindConnectionWidths, FindTriangleAreas, FindVoronoiCellAreas, GetLatLonCoordinates, DetermineMeshResolution, &
                                         CheckMesh, Cross, Cart_Bilinear_dp, MaxCartOverTriangle_dp, Mesh_Bilinear_dp, Mesh_Bilinear_int, IsInTriangle, &
                                         DetermineProcessXYDomain_regular, DetermineProcessXYDomain_balanced, NewTriangleContainsOldMask, WriteMeshToTextFile, &
                                         IsBoundarySegment, IsEncroachedUpon
  USE mesh_memory_module,          ONLY: AllocateMesh, AllocateSubmesh, AllocateMesh_extra, ExtendSubmesh, CropMeshMemory, &
                                         DeallocateSubmesh, MoveDataFromSubmeshToMesh
  USE mesh_Delaunay_module,        ONLY: SplitTriangle
  USE mesh_derivatives_module,     ONLY: GetNeighbourFunctions, GetMeshCurvatures, GetMeshCurvaturesVertex
  USE mesh_creation_module,        ONLY: InitialiseDummyMesh, PerturbDummyMesh, AlignAllSubmeshes, RefineSubmesh_geo_only, &
                                         MergeSubmeshesMeta, CreateFinalMeshFromMergedSubmesh
  USE mesh_ArakawaC_module,        ONLY: MakeAcMesh

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'

  CONTAINS
  
  ! === Check whether or not a triangle meets all the fitness criteria.
  ! If you want to change the rules for mesh creation, this is where to do it.
  SUBROUTINE IsGoodTriangle(mesh_new, ti, mesh_old, ice, PD, IsGood, C_mesh)
    ! Check if triangle ti of the mesh is Bad (and should be refined by Ruppert's Algorithm)
    ! A triangle is Bad if:
    !   - its smallest internal angle is too small
    !   - its 2nd order surface deviation (=max(curvature)*typical_length) is too large
    !   - its area exceeds the limits based on ice velocity, grounding line or calving front

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh_new
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh_old       ! The old mesh
    TYPE(type_ice_model),       INTENT(IN)        :: ice
    TYPE(type_PD_data_fields),  INTENT(IN)        :: PD
    TYPE(type_mesh_config),     INTENT(IN)        :: C_mesh

    INTEGER,                    INTENT(IN)        :: ti             ! The triangle we want to check
    LOGICAL,                    INTENT(OUT)       :: IsGood         ! The result

    REAL(dp)                                      :: max_curv    
    REAL(dp)                                      :: mean_mask_dp_p, mean_mask_dp_q, mean_mask_dp_r, mean_mask_dp_m    
    REAL(dp), DIMENSION(2)                        :: p,q,r,m
    INTEGER                                       :: vp,vq,vr
    REAL(dp), DIMENSION(2)                        :: pq,qr,rp
    LOGICAL                                       :: isso
    REAL(dp)                                      :: dmax, dz
    INTEGER                                       :: n, vi
    REAL(dp)                                      :: ap,aq,ar,alpha
    INTEGER                                       :: trinel  
    LOGICAL                                       :: contains_ice, contains_land, contains_ocean, contains_sheet, contains_shelf
    LOGICAL                                       :: contains_coast, contains_margin, contains_gl, contains_cf
    REAL(dp), PARAMETER                           :: Hb_lo = 500._dp
    REAL(dp), PARAMETER                           :: Hb_hi = 1500._dp  
    REAL(dp)                                      :: Hb_max, w_Hb, lr_lo, lr_hi, r_crit
    REAL(dp)                                      :: x_range_mesh, y_range_mesh

    IsGood = .TRUE.
    
    x_range_mesh = mesh_new%xmax - mesh_new%xmin
    y_range_mesh = mesh_new%ymax - mesh_new%ymin

    ! Triangle geometry (the basis of the original version of Rupperts Algorithm)
    ! ===========================================================================
    
    p = [MIN(mesh_new%xmax - x_range_mesh/1E6_dp, MAX(mesh_new%xmin + x_range_mesh/1E6_dp, mesh_new%V(mesh_new%Tri(ti,1),1) )), &
         MIN(mesh_new%ymax - y_range_mesh/1E6_dp, MAX(mesh_new%ymin + y_range_mesh/1E6_dp, mesh_new%V(mesh_new%Tri(ti,1),2) ))]
    q = [MIN(mesh_new%xmax - x_range_mesh/1E6_dp, MAX(mesh_new%xmin + x_range_mesh/1E6_dp, mesh_new%V(mesh_new%Tri(ti,2),1) )), &
         MIN(mesh_new%ymax - y_range_mesh/1E6_dp, MAX(mesh_new%ymin + y_range_mesh/1E6_dp, mesh_new%V(mesh_new%Tri(ti,2),2) ))]
    r = [MIN(mesh_new%xmax - x_range_mesh/1E6_dp, MAX(mesh_new%xmin + x_range_mesh/1E6_dp, mesh_new%V(mesh_new%Tri(ti,3),1) )), &
         MIN(mesh_new%ymax - y_range_mesh/1E6_dp, MAX(mesh_new%ymin + y_range_mesh/1E6_dp, mesh_new%V(mesh_new%Tri(ti,3),2) ))]
    m = (p+q+r)/3._dp

    vp = mesh_new%Tri(ti,1)
    vq = mesh_new%Tri(ti,2)
    vr = mesh_new%Tri(ti,3)

    ! Triangle legs
    pq = p-q
    qr = q-r
    rp = r-p

    ! Longest triangle leg
    dmax = MAXVAL([NORM2(pq), NORM2(qr), NORM2(rp)])

    ! Internal angles
    ap = ACOS(-(rp(1)*pq(1) + rp(2)*pq(2))/(NORM2(rp)*NORM2(pq)))
    aq = ACOS(-(pq(1)*qr(1) + pq(2)*qr(2))/(NORM2(pq)*NORM2(qr)))
    ar = ACOS(-(rp(1)*qr(1) + rp(2)*qr(2))/(NORM2(rp)*NORM2(qr)))

    ! Smallest internal angle
    alpha = MINVAL([ap, aq, ar])

    IF (alpha < C%mesh_alpha_min) THEN
      IsGood = .FALSE.
      RETURN
    END IF
    
    ! If its an edge triangle, check if the third vertex encroaches on the edge segment
    IF (IsBoundarySegment( mesh_new, vp, vq)) THEN
      CALL IsEncroachedUpon( mesh_new, vp, vq, isso)
      IF (isso) THEN
        IsGood = .FALSE.
        RETURN
      END IF
    ELSEIF (IsBoundarySegment( mesh_new, vq, vr)) THEN
      CALL IsEncroachedUpon( mesh_new, vq, vr, isso)
      IF (isso) THEN
        IsGood = .FALSE.
        RETURN
      END IF
    ELSEIF (IsBoundarySegment( mesh_new, vr, vp)) THEN
      CALL IsEncroachedUpon( mesh_new, vr, vp, isso)
      IF (isso) THEN
        IsGood = .FALSE.
        RETURN
      END IF
    END IF
    
    ! Coarsest allowed resolution
    ! ============================
    
    IF (dmax > C_mesh%res_max * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN
    END IF
    
    
    ! Finest allowed resolution
    ! =========================
    
    IF (dmax < C_mesh%res_min * 2._dp * 1000._dp) THEN
      IsGood = .TRUE.
      RETURN
    END IF
    
    ! Resolution at points of interest
    ! ================================
    
    DO n = 1, mesh_new%nPOI
      IF (IsInTriangle( p, q, r, mesh_new%POI_XY_coordinates(n,:)) .AND. dmax > mesh_new%POI_resolutions(n) * 1.5_dp * 1000._dp) THEN
        IsGood = .FALSE.
        RETURN
      END IF
    END DO

    ! Determine what's inside the triangle
    ! ====================================
    
    vi = 1
    CALL NewTriangleContainsOldMask( mesh_old, p, q, r, ice%mask_ice,           vi, contains_ice   )
    CALL NewTriangleContainsOldMask( mesh_old, p, q, r, ice%mask_land,          vi, contains_land  )
    CALL NewTriangleContainsOldMask( mesh_old, p, q, r, ice%mask_ocean,         vi, contains_ocean )
    CALL NewTriangleContainsOldMask( mesh_old, p, q, r, ice%mask_sheet,         vi, contains_sheet )
    CALL NewTriangleContainsOldMask( mesh_old, p, q, r, ice%mask_shelf,         vi, contains_shelf )
    CALL NewTriangleContainsOldMask( mesh_old, p, q, r, ice%mask_coast,         vi, contains_coast )
    CALL NewTriangleContainsOldMask( mesh_old, p, q, r, ice%mask_margin,        vi, contains_margin)
    CALL NewTriangleContainsOldMask( mesh_old, p, q, r, ice%mask_groundingline, vi, contains_gl    )
    CALL NewTriangleContainsOldMask( mesh_old, p, q, r, ice%mask_calvingfront,  vi, contains_cf    )
        
    ! Special area resolution - ice margin, grounding line, calving front
    ! ===================================================================
    
    ! Coastline
    IF (contains_coast .AND. dmax > C_mesh%res_max_coast * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN
    END IF
    
    ! Ice margin
    IF (contains_margin .AND. dmax > C_mesh%res_max_margin * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN
    END IF

    ! Grounding line
    IF (contains_gl .AND. dmax > C_mesh%res_max_gl * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN
    END IF

    ! Calving front
    IF (contains_cf .AND. dmax > C_mesh%res_max_cf * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN
    END IF

    ! Second-order surface deviation (curvature times size)
    ! =====================================================

    CALL Mesh_Bilinear_dp(mesh_old, ice%surf_curv, p, mesh_new%mesh_old_ti_in(vp), mean_mask_dp_p)
    CALL Mesh_Bilinear_dp(mesh_old, ice%surf_curv, q, mesh_new%mesh_old_ti_in(vq), mean_mask_dp_q)
    CALL Mesh_Bilinear_dp(mesh_old, ice%surf_curv, r, mesh_new%mesh_old_ti_in(vr), mean_mask_dp_r)
    CALL Mesh_Bilinear_dp(mesh_old, ice%surf_curv, m, mesh_new%mesh_old_ti_in(vp), mean_mask_dp_m)

    max_curv = MAXVAL([mean_mask_dp_p, mean_mask_dp_q, mean_mask_dp_r, mean_mask_dp_m])
    dz = 0.5_dp * max_curv * dmax**2

    IF (contains_ice .AND. dz > C_mesh%dz_max_ice) THEN
      IsGood = .FALSE.
     RETURN
    END IF
    
    ! Ice-free bed topography (higher res for mountains so inception is captured better)
    ! ==================================================================================
    
    CALL MaxCartOverTriangle_dp( p, q, r, PD%Hb_cart, PD%x, PD%y, PD%nx, PD%ny, Hb_max, trinel)
    IF (trinel==0) THEN
      CALL Cart_Bilinear_dp( PD%Hb_cart, PD%x, PD%y, PD%nx, PD%ny, (p+q+r)/3._dp, Hb_max)
    END IF
    
    lr_lo = LOG(C_mesh%res_max)
    lr_hi = LOG(C_mesh%res_max_mountain)
    
    w_Hb = MIN(1._dp, MAX(0._dp, (Hb_max - Hb_lo) / (Hb_hi - Hb_lo)))    
    r_crit = EXP( (w_Hb * lr_hi) + ((1._dp - w_Hb) * lr_lo))
    
    IF (contains_land .AND. dmax > r_crit * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN      
    END IF


  END SUBROUTINE IsGoodTriangle
  
  ! == Mesh creation routines ==
  SUBROUTINE CreateNewMesh(region)
    ! Called by all processes, but only the master can actually do stuff. Other processes
    ! only need to run the memory extension routines to make sure pointers are kept up to date.
  
    ! Input variables
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    ! Local variables
    INTEGER                                       :: vi
    REAL(dp)                                      :: d2dx2, d2dxdy, d2dy2
    
    INTEGER                                       :: orientation
    TYPE(type_mesh)                               :: submesh
    REAL(dp)                                      :: xmin, xmax, ymin, ymax
    TYPE(type_mesh_config)                        :: C_mesh_inc, C_mesh_final
    REAL(dp)                                      :: res_min
    INTEGER                                       :: rit
    
    CHARACTER(LEN=256) :: submesh_filename, str1, str2
    
    IF (par%master) WRITE(0,*) '  Creating a new mesh for region ', region%mesh%region_name, '...'
    
    ! Orientation of domain partitioning: east-west for GRL, north-south everywhere else
    IF (region%name == 'GRL') THEN
      orientation = 1
    ELSE
      orientation = 0
    END IF
    
    ! Make the ice margin 1 row of elements wider. This slightly increases the vertex count of
    ! the new mesh, but means the new mesh stays Fit much longer.
    CALL WidenHighResZones( region%ice, region%mesh, region%time)
    
    ! Calculate surface curvature
    DO vi = region%mesh%v1, region%mesh%v2
      CALL GetMeshCurvaturesVertex( region%mesh, region%ice%Hs, d2dx2, d2dxdy, d2dy2, vi)
      region%ice%surf_curv(vi) = MAX(-1E-6, MIN(1E-6, SQRT(d2dx2**2 + d2dy2**2 + d2dxdy**2)))
    END DO
    
    ! Determine the domain of this process' submesh (based on distribution of vertices in
    ! the previous mesh, which works very well for workload balancing)
    CALL DetermineProcessXYDomain_regular( region,      par%i, orientation, xmin, xmax, ymin, ymax)
!    CALL DetermineProcessXYDomain_balanced(region%mesh, par%i, orientation, xmin, xmax, ymin, ymax)
    
    ! Allocate memory, initialise a dummy mesh, refine it, and crop the memory
    CALL AllocateSubmesh(     submesh, region%name, 10)  
    CALL InitialiseDummyMesh( submesh, xmin, xmax, ymin, ymax)
    CALL PerturbDummyMesh(    submesh, 1 - region%mesh%perturb_dir)
        
    ! Determine the parameters for the final mesh
    C_mesh_final%dz_max_ice       = C%dz_max_ice
    C_mesh_final%res_min          = C%res_min
    C_mesh_final%res_max          = C%res_max
    C_mesh_final%res_max_margin   = C%res_max_margin
    C_mesh_final%res_max_gl       = C%res_max_gl
    C_mesh_final%res_max_cf       = C%res_max_cf
    C_mesh_final%res_max_mountain = C%res_max_mountain
    C_mesh_final%res_max_coast    = C%res_max_coast
    
    ! Start with everything at the maximum allowed resolution
    C_mesh_inc%dz_max_ice = C_mesh_final%dz_max_ice
    C_mesh_inc%res_max    = C_mesh_final%res_max    
    res_min               = C_mesh_final%res_max
    
    rit = 0
    DO WHILE (res_min > C_mesh_final%res_min)
    
      rit = rit + 1
!      IF (par%master) WRITE(0,'(A,I3,A,F7.2,A)') '    Resolution increment ', rit, ': res_min = ', res_min, ' km'
    
      ! Determine resolutions
      C_mesh_inc%res_min          = MAX( C_mesh_final%res_min,          res_min)
      C_mesh_inc%res_max_margin   = MAX( C_mesh_final%res_max_margin,   res_min)
      C_mesh_inc%res_max_gl       = MAX( C_mesh_final%res_max_gl,       res_min)
      C_mesh_inc%res_max_cf       = MAX( C_mesh_final%res_max_cf,       res_min)
      C_mesh_inc%res_max_mountain = MAX( C_mesh_final%res_max_mountain, res_min)
      C_mesh_inc%res_max_coast    = MAX( C_mesh_final%res_max_coast,    res_min)
      
      ! Refine the process submesh
      CALL RefineMesh( submesh, region%mesh, region%ice, region%PD, C_mesh_inc) 
      
      ! Align with neighbouring submeshes
      CALL AlignAllSubmeshes(      submesh, orientation)
      CALL RefineSubmesh_geo_only( submesh)
      
      ! DENK DROM
      WRITE(str1,'(I2)') par%i;   str1 = ADJUSTL(str1)
      WRITE(str2,'(I2)') rit;     str2 = ADJUSTL(str2)
      submesh_filename = 'submesh_p' // TRIM(str1) // '_ri' // TRIM(str2) // '.txt'
!      CALL WriteMeshToTextFile( submesh, submesh_filename)
      
      ! Increase resolution
      res_min = res_min / 4._dp
    
    END DO
    
    ! One final pass with the desired resolution settings 
    CALL RefineMesh( submesh, region%mesh, region%ice, region%PD, C_mesh_final) 
    
    ! Merge the process submeshes, create the final shared-memory mesh
    CALL MergeSubmeshesMeta( submesh, orientation)
    CALL CreateFinalMeshFromMergedSubmesh( submesh, region%mesh_new)

    IF (par%master) THEN
      WRITE(0,'(A)')                '   Finished creating final mesh.'
      WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh_new%nV
      WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh_new%nTri
      WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh_new%resolution_min/1000._dp, ' - ', region%mesh_new%resolution_max/1000._dp, ' km'
    END IF
    
  END SUBROUTINE CreateNewMesh
  
  SUBROUTINE WidenHighResZones( ice, mesh, time)
    ! Make the three lines of interest (ice margin, grounding line, calving front) wider by twice their
    ! config-specified resolution. This makes the high-resolution area in the new mesh slightly wider,
    ! and prevents issues where a mesh update results in slight relaxation of the grounding line,
    ! moving it out of the high-res area, triggering a mesh update, etc.
    ! Only done by master, no easy way to parallelise this.
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables
    INTEGER                                            :: vi, ci, vc
    REAL(dp)                                           :: D
    INTEGER, DIMENSION(mesh%nV)                        :: mask_widened
    REAL(dp), PARAMETER                                :: widening_factor = 1.2_dp
    
    ! Exception for the first mesh update in benchmark experiments. Since these start with
    ! zero ice thickness, at the time of the first mesh update the ice margin will lie in a
    ! very low resolution area. Widening the margin will result in a huge swath of the domain getting
    ! a high resolution after the update, which will prevent any further updates.
    IF (time < C%start_time_of_run + 100._dp) RETURN
            
    IF (par%master) THEN
    
  ! Ice margin
  ! ==========
    
      mask_widened = ice%mask_margin
      
      DO vi = 1, mesh%nV
        IF (ice%mask_margin(vi)==1) THEN
          DO ci = 1, mesh%nC(vi)
            vc = mesh%C(vi,ci)
            D = NORM2( mesh%V(vi,:) - mesh%V(vc,:))
            IF (D < C%res_max_margin * 1000._dp * widening_factor) THEN
              mask_widened(vc) = 1
            END IF
          END DO
        END IF
      END DO
      
      ice%mask_margin = mask_widened
    
  ! Grounding line
  ! ==============
    
      mask_widened = ice%mask_groundingline
      
      DO vi = 1, mesh%nV
        IF (ice%mask_groundingline(vi)==1) THEN
          DO ci = 1, mesh%nC(vi)
            vc = mesh%C(vi,ci)
            D = NORM2( mesh%V(vi,:) - mesh%V(vc,:))
            IF (D < C%res_max_gl * 1000._dp * widening_factor) THEN
              mask_widened(vc) = 1
            END IF
          END DO
        END IF
      END DO
      
      ice%mask_groundingline = mask_widened
    
  ! Calving front
  ! =============
    
      mask_widened = ice%mask_calvingfront
      
      DO vi = 1, mesh%nV
        IF (ice%mask_calvingfront(vi)==1) THEN
          DO ci = 1, mesh%nC(vi)
            vc = mesh%C(vi,ci)
            D = NORM2( mesh%V(vi,:) - mesh%V(vc,:))
            IF (D < C%res_max_cf * 1000._dp * widening_factor) THEN
              mask_widened(vc) = 1
            END IF
          END DO
        END IF
      END DO
      
      ice%mask_calvingfront = mask_widened
      
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE WidenHighResZones
  SUBROUTINE WidenHighResZones_old( ice, mesh, time)
    ! Make the three lines of interest (ice margin, grounding line, calving front) one triangle
    ! wider in the relevant masks. This makes the high-resolution area in the new mesh slightly wider,
    ! and prevents issues where a mesh update results in slight relaxation of the grounding line,
    ! moving it out of the high-res area, triggering a mesh update, etc.
    ! Only done by master, no easy way to parallelise this.
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables
    INTEGER, PARAMETER                                 :: number_of_widening_loops = 1
    INTEGER                                            :: vi, ci, vc
    INTEGER, DIMENSION(mesh%nV)                        :: mask_widened
    INTEGER                                            :: it
    
    ! Exception for the first mesh update in benchmark experiments. Since these start with
    ! zero ice thickness, at the time of the first mesh update the ice margin will lie in a
    ! very low resolution area. Widening the margin will result in a huge swath of the domain getting
    ! a high resolution after the update, which will prevent any further updates.
    IF (time < C%start_time_of_run + 100._dp) RETURN
            
    IF (par%master) THEN
    
      DO it = 1, number_of_widening_loops
    
  ! Ice margin
  ! ==========
    
      mask_widened = ice%mask_margin
      
      DO vi = 1, mesh%nV
        IF (ice%mask_margin(vi)==1) THEN
          DO ci = 1, mesh%nC(vi)
            vc = mesh%C(vi,ci)
            mask_widened(vc) = 1
          END DO
        END IF
      END DO
      
      ice%mask_margin = mask_widened
    
  ! Grounding line
  ! ==============
    
      mask_widened = ice%mask_groundingline
      
      DO vi = 1, mesh%nV
        IF (ice%mask_groundingline(vi)==1) THEN
          DO ci = 1, mesh%nC(vi)
            vc = mesh%C(vi,ci)
            mask_widened(vc) = 1
          END DO
        END IF
      END DO
      
      ice%mask_groundingline = mask_widened
    
  ! Calving front
  ! =============
    
      mask_widened = ice%mask_calvingfront
      
      DO vi = 1, mesh%nV
        IF (ice%mask_calvingfront(vi)==1) THEN
          DO ci = 1, mesh%nC(vi)
            vc = mesh%C(vi,ci)
            mask_widened(vc) = 1
          END DO
        END IF
      END DO
      
      ice%mask_calvingfront = mask_widened
      
      END DO
      
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE WidenHighResZones_old
  
  ! == Extended Ruppert's algorithm
  SUBROUTINE RefineMesh( mesh, mesh_old, ice, PD, C_mesh)
    ! Refine a mesh. Single-core, but multiple submeshes can be done in parallel on different cores.
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh_old
    TYPE(type_ice_model),       INTENT(IN)        :: ice
    TYPE(type_PD_data_fields),  INTENT(IN)        :: PD
    TYPE(type_mesh_config),     INTENT(IN)        :: C_mesh
    
    ! Local variables
    INTEGER                                       :: ti, i, ierr, status(MPI_STATUS_SIZE)
    LOGICAL                                       :: IsGood, FinishedRefining, DoExtendMemory
    REAL(dp), DIMENSION(2)                        :: p
    
    FinishedRefining = .FALSE.
    DoExtendMemory   = .FALSE.
    
    CALL ExtendSubmesh( mesh, mesh%nV + 1000)
    
    mesh%RefMap    = 0
    mesh%RefStack  = 0
    mesh%RefStackN = 0
    DO ti = 1, mesh%nTri
      mesh%RefMap(ti)               = 1
      mesh%RefStackN                = mesh%RefStackN + 1
      mesh%RefStack(mesh%RefStackN) = ti
    END DO
    
    mesh%mesh_old_ti_in = 1
    
    DO WHILE (.NOT. FinishedRefining)
    
      ! Refine the mesh until it's done, or until it's almost out of memory.
      ! ====================================================================
      
      DoExtendMemory = .FALSE.
          
      DO WHILE (mesh%RefStackN > 0)
    
        ! Check the last triangle list in the RefineStack. If it's
        ! Bad, refine it and add the affected triangles to the RefineStack.
      
        ti = mesh%RefStack( mesh%RefStackN)
        CALL IsGoodTriangle(mesh, ti, mesh_old, ice, PD, IsGood, C_mesh)
              
        IF (IsGood) THEN
          ! Remove this triangle from the stack
          mesh%RefMap(ti) = 0
          mesh%RefStack(mesh%RefStackN) = 0
          mesh%RefStackN = mesh%RefStackN - 1
        ELSE
          ! Spit this triangle, add the affected triangles to the stack
          p = mesh%Tricc(ti,:)
          CALL SplitTriangle( mesh, ti, p)
        END IF
        
        ! If we're reaching the memory limit, stop refining and extend the memory.
        IF (mesh%nV > mesh%nV_max - 10) THEN
          DoExtendMemory = .TRUE.
          EXIT
        END IF
        
      END DO ! DO WHILE (mesh%RefStackN > 0)
      
      ! Check if all processes finished refining. If so, exit.
      ! ======================================================
      
      FinishedRefining = .FALSE.
      IF (mesh%RefStackN == 0) FinishedRefining = .TRUE.      
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, FinishedRefining, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)      
      IF (FinishedRefining) EXIT  
      
      ! Check if any process needs to extend their memory.
      ! ==================================================
      
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, DoExtendMemory, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr) 
            
      ! By extending the memory to mesh%nV + 1000, we ensure that processes that 
      ! have already finished refining do not keep adding useless extra memory.
      
      IF (DoExtendMemory) THEN
        CALL ExtendSubmesh( mesh, mesh%nV + 1000)
      END IF
    
    END DO ! DO WHILE (.NOT. FinishedRefining)
  
  END SUBROUTINE RefineMesh
  
  ! == Determine if mesh updating is needed
  SUBROUTINE DetermineMeshFitness(mesh, ice, fitness)
    ! Determine how "fit" the current mesh is.
  
    ! Input variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_ice_model),       INTENT(INOUT)     :: ice
    
    ! Output variables
    REAL(dp),                   INTENT(OUT)       :: fitness
    
    ! Local variables
    INTEGER                                       :: ti, i, ierr, status(MPI_STATUS_SIZE)
    INTEGER                                       :: v1, v2, v3
    REAL(dp), DIMENSION(2)                        :: p,q,r,pq,qr,rp
    REAL(dp)                                      :: dmax
    REAL(dp), PARAMETER                           :: res_tol = 1.2_dp ! Resolution tolerance factor
    INTEGER                                       :: ncoast, nucoast, nmargin, numargin, ngl, nugl, ncf, nucf
    REAL(dp)                                      :: lcoast, lucoast, lmargin, lumargin, lgl, lugl, lcf, lucf
    REAL(dp)                                      :: fcoast, fmargin, fgl, fcf
        
    fitness = 1._dp
   
    ! Determine fraction of fit triangles
    
    ncoast   = 0 ! Number of       coastline triangles
    nmargin  = 0 ! Number of       margin triangles
    ngl      = 0 ! Mumber of       grounding line triangles
    ncf      = 0 ! Number of       calving front triangles
    nucoast  = 0 ! Number of unfit coastline triangles
    numargin = 0 ! Number of unfit margin triangles
    nugl     = 0 ! Mumber of unfit grounding line triangles
    nucf     = 0 ! Number of unfit calving front triangles
    
    lcoast   = 0._dp ! Total length of coastline
    lmargin  = 0._dp
    lgl      = 0._dp
    lcf      = 0._dp
    lucoast  = 0._dp ! Unfit length of coastline
    lumargin = 0._dp
    lugl     = 0._dp
    lucf     = 0._dp
    
    DO ti = mesh%t1, mesh%t2
     
      ! Triangle vertex indices
      v1 = mesh%Tri(ti,1)
      v2 = mesh%Tri(ti,2)
      v3 = mesh%Tri(ti,3)
    
      ! Triangle vertex coordinates
      p = mesh%V(v1,:)
      q = mesh%V(v2,:)
      r = mesh%V(v3,:)

      ! Triangle legs
      pq = p-q
      qr = q-r
      rp = r-p
      
      ! Longest triangle leg
      dmax = MAXVAL([SQRT(pq(1)**2+pq(2)**2), SQRT(qr(1)**2+qr(2)**2), SQRT(rp(1)**2+rp(2)**2)])
      
      IF (ice%mask_coast(v1)==1 .OR. ice%mask_coast(v2)==1 .OR. ice%mask_coast(v3)==1) THEN
        ncoast = ncoast + 1
        lcoast = lcoast + dmax
        IF (dmax > C%res_max_coast*2.0_dp*1000._dp*res_tol) THEN
          nucoast = nucoast + 1
          lucoast = lucoast + dmax
        END IF
      END IF
      
      IF (ice%mask_margin(v1)==1 .OR. ice%mask_margin(v2)==1 .OR. ice%mask_margin(v3)==1) THEN
        nmargin = nmargin + 1
        lmargin = lmargin + dmax
        IF (dmax > C%res_max_margin*2.0_dp*1000._dp*res_tol) THEN
          numargin = numargin + 1
          lumargin = lumargin + dmax
        END IF
      END IF
      IF (ice%mask_groundingline(v1)==1 .OR. ice%mask_groundingline(v2)==1 .OR. ice%mask_groundingline(v3)==1) THEN
        ngl = ngl + 1
        lgl = lgl + dmax
        IF (dmax > C%res_max_gl*2.0_dp*1000._dp*res_tol) THEN
          nugl = nugl + 1
          lugl = lugl + dmax
        END IF
      END IF
      IF (ice%mask_calvingfront(v1)==1 .OR. ice%mask_calvingfront(v2)==1 .OR. ice%mask_calvingfront(v3)==1) THEN
        ncf = ncf + 1
        lcf = lcf + dmax
        IF (dmax > C%res_max_cf*2.0_dp*1000._dp*res_tol) THEN
          nucf = nucf + 1
          lucf = lucf + dmax
        END IF
      END IF
      
    END DO ! DO ti = mesh%t1, mesh%t2
    CALL sync
    
    ! Gather mesh fitness data from all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, ncoast,   1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nmargin,  1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, ngl,      1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, ncf,      1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nucoast,  1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, numargin, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nugl,     1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nucf,     1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lcoast,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lmargin,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lgl,      1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lcf,      1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lucoast,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lumargin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lugl,     1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lucf,     1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    ! Calculate mesh fitness    
    fcoast  = 1._dp - lucoast  / lcoast
    fmargin = 1._dp - lumargin / lmargin
    fgl     = 1._dp - lugl     / lgl
    fcf     = 1._dp - lucf     / lcf
    
    IF (ncoast ==0) fcoast  = 1._dp
    IF (nmargin==0) fmargin = 1._dp
    IF (ngl    ==0) fgl     = 1._dp
    if (ncf    ==0) fcf     = 1._dp
    
    fitness = MIN( MIN( MIN( fcoast, fmargin), fgl), fcf)
    
    IF (par%master) THEN
      DO i = 1, par%n-1
        CALL MPI_SEND( fitness, 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, ierr)
      END DO      
    ELSE
      CALL MPI_RECV( fitness, 1, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
    END IF
    CALL sync
      
    !IF (par%master) WRITE(0,'(A,I3,A)') '   Mesh fitness: ', NINT(meshfitness * 100._dp), ' %'
    
        
  END SUBROUTINE DetermineMeshFitness

END MODULE mesh_update_module
