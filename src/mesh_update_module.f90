MODULE mesh_update_module

  ! Routines for creating a new mesh based on forcing data on an old mesh.

  ! Import basic functionality
  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, write_to_memory_log, &
                                             allocate_shared_int_0D,   allocate_shared_dp_0D, &
                                             allocate_shared_int_1D,   allocate_shared_dp_1D, &
                                             allocate_shared_int_2D,   allocate_shared_dp_2D, &
                                             allocate_shared_int_3D,   allocate_shared_dp_3D, &
                                             allocate_shared_bool_0D,  allocate_shared_bool_1D, &
                                             reallocate_shared_int_0D, reallocate_shared_dp_0D, &
                                             reallocate_shared_int_1D, reallocate_shared_dp_1D, &
                                             reallocate_shared_int_2D, reallocate_shared_dp_2D, &
                                             reallocate_shared_int_3D, reallocate_shared_dp_3D, &
                                             deallocate_shared
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  
  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh, type_model_region, type_ice_model, type_reference_geometry
  USE mesh_help_functions_module,      ONLY: cart_bilinear_dp, max_cart_over_triangle_dp, mesh_bilinear_dp, is_in_triangle, check_mesh, &
                                             partition_domain_regular, partition_domain_x_balanced, partition_domain_y_balanced, is_walltowall, &
                                             new_triangle_contains_old_mask, write_mesh_to_text_file, is_boundary_segment, is_encroached_upon
  USE mesh_memory_module,              ONLY: allocate_submesh_primary, extend_submesh_primary
  USE mesh_Delaunay_module,            ONLY: split_triangle
  USE mesh_creation_module,            ONLY: initialise_dummy_mesh, perturb_dummy_mesh, align_all_submeshes, refine_submesh_geo_only, debug_mesh_creation, &
                                             merge_all_submeshes, create_final_mesh_from_merged_submesh, Lloyds_algorithm_single_iteration_submesh
  USE mesh_operators_module,           ONLY: d2dx2_a_to_a_2D, d2dxdy_a_to_a_2D, d2dy2_a_to_a_2D

  IMPLICIT NONE

  CONTAINS
  
  ! === Check whether or not a triangle meets all the fitness criteria.
  ! If you want to change the rules for mesh creation, this is where to do it.
  SUBROUTINE is_good_triangle( mesh_new, ti, mesh_old, ice, refgeo_PD, IsGood)
    ! Check if triangle ti of the mesh is Bad (and should be refined by Ruppert's Algorithm)
    ! A triangle is Bad if:
    !   - its smallest internal angle is too small
    !   - its 2nd order surface deviation (=max(curvature)*typical_length) is too large
    !   - its area exceeds the limits based on ice velocity, grounding line or calving front
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh_new
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh_old       ! The old mesh
    TYPE(type_ice_model),       INTENT(IN)        :: ice
    TYPE(type_reference_geometry), INTENT(IN)     :: refgeo_PD

    INTEGER,                    INTENT(IN)        :: ti             ! The triangle we want to check
    LOGICAL,                    INTENT(OUT)       :: IsGood         ! The result

    REAL(dp)                                      :: max_curv    
    REAL(dp)                                      :: mean_mask_dp_p, mean_mask_dp_q, mean_mask_dp_r, mean_mask_dp_m    
    INTEGER                                       :: vp, vq, vr
    REAL(dp), DIMENSION(2)                        :: p, q, r, m, POI
    REAL(dp), DIMENSION(2)                        :: pq, qr, rp
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

    IF (alpha < mesh_new%alpha_min) THEN
      IsGood = .FALSE.
      RETURN
    END IF
    
    ! If its an edge triangle, check if the third vertex encroaches on the edge segment
    IF (is_boundary_segment( mesh_new, vp, vq)) THEN
      CALL is_encroached_upon( mesh_new, vp, vq, isso)
      IF (isso) THEN
        IsGood = .FALSE.
        RETURN
      END IF
    ELSEIF (is_boundary_segment( mesh_new, vq, vr)) THEN
      CALL is_encroached_upon( mesh_new, vq, vr, isso)
      IF (isso) THEN
        IsGood = .FALSE.
        RETURN
      END IF
    ELSEIF (is_boundary_segment( mesh_new, vr, vp)) THEN
      CALL is_encroached_upon( mesh_new, vr, vp, isso)
      IF (isso) THEN
        IsGood = .FALSE.
        RETURN
      END IF
    END IF
    
    ! Forbid "wall to wall" triangles (i.e. triangles with vertices lying on opposite domain boundaries)
    IF (is_walltowall( mesh_new, ti)) THEN
      IsGood = .FALSE.
      RETURN
    END IF
    
    ! Coarsest allowed resolution
    ! ============================
    
    IF (dmax > mesh_new%res_max * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN
    END IF
    
    
    ! Finest allowed resolution
    ! =========================
    
    IF (dmax < mesh_new%res_min * 2._dp * 1000._dp) THEN
      IsGood = .TRUE.
      RETURN
    END IF
    
    ! Resolution at points of interest
    ! ================================
    
    DO n = 1, mesh_new%nPOI
      POI = mesh_new%POI_XY_coordinates(n,:)
      IF (is_in_triangle( p, q, r, POI) .AND. dmax > mesh_new%POI_resolutions(n) * 1.5_dp * 1000._dp) THEN
        IsGood = .FALSE.
        RETURN
      END IF
    END DO

    ! Determine what's inside the triangle
    ! ====================================
    
    vi = 1
    CALL new_triangle_contains_old_mask( mesh_old, p, q, r, ice%mask_ice_a,    vi, contains_ice   )
    CALL new_triangle_contains_old_mask( mesh_old, p, q, r, ice%mask_land_a,   vi, contains_land  )
    CALL new_triangle_contains_old_mask( mesh_old, p, q, r, ice%mask_ocean_a,  vi, contains_ocean )
    CALL new_triangle_contains_old_mask( mesh_old, p, q, r, ice%mask_sheet_a,  vi, contains_sheet )
    CALL new_triangle_contains_old_mask( mesh_old, p, q, r, ice%mask_shelf_a,  vi, contains_shelf )
    CALL new_triangle_contains_old_mask( mesh_old, p, q, r, ice%mask_coast_a,  vi, contains_coast )
    CALL new_triangle_contains_old_mask( mesh_old, p, q, r, ice%mask_margin_a, vi, contains_margin)
    CALL new_triangle_contains_old_mask( mesh_old, p, q, r, ice%mask_gl_a,     vi, contains_gl    )
    CALL new_triangle_contains_old_mask( mesh_old, p, q, r, ice%mask_cf_a,     vi, contains_cf    )
        
    ! Special area resolution - ice margin, grounding line, calving front
    ! ===================================================================
    
    ! Coastline
    IF (contains_coast .AND. dmax > mesh_new%res_max_coast * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN
    END IF
    
    ! Ice margin
    IF (contains_margin .AND. dmax > mesh_new%res_max_margin * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN
    END IF

    ! Grounding line
    IF (contains_gl .AND. dmax > mesh_new%res_max_gl * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN
    END IF

    ! Calving front
    IF (contains_cf .AND. dmax > mesh_new%res_max_cf * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN
    END IF

    ! Second-order surface deviation (curvature times size)
    ! =====================================================

    CALL mesh_bilinear_dp( mesh_old, ice%surf_curv, p, mesh_new%mesh_old_ti_in( vp), mean_mask_dp_p)
    CALL mesh_bilinear_dp( mesh_old, ice%surf_curv, q, mesh_new%mesh_old_ti_in( vq), mean_mask_dp_q)
    CALL mesh_bilinear_dp( mesh_old, ice%surf_curv, r, mesh_new%mesh_old_ti_in( vr), mean_mask_dp_r)
    CALL mesh_bilinear_dp( mesh_old, ice%surf_curv, m, mesh_new%mesh_old_ti_in( vp), mean_mask_dp_m)

    max_curv = MAXVAL([mean_mask_dp_p, mean_mask_dp_q, mean_mask_dp_r, mean_mask_dp_m])
    dz = 0.5_dp * max_curv * dmax**2

    IF (contains_ice .AND. dz > mesh_new%dz_max_ice) THEN
      IsGood = .FALSE.
      RETURN
    END IF
    
    ! Ice-free bed topography (higher res for mountains so inception is captured better)
    ! ==================================================================================
    
    CALL max_cart_over_triangle_dp( p, q, r, refgeo_PD%Hb_grid, refgeo_PD%grid%x, refgeo_PD%grid%y, refgeo_PD%grid%nx, refgeo_PD%grid%ny, Hb_max, trinel)
    IF (trinel==0) THEN
      CALL cart_bilinear_dp( refgeo_PD%Hb_grid, refgeo_PD%grid%x, refgeo_PD%grid%y, refgeo_PD%grid%nx, refgeo_PD%grid%ny, (p+q+r)/3._dp, Hb_max)
    END IF
    
    lr_lo = LOG( mesh_new%res_max)
    lr_hi = LOG( mesh_new%res_max_mountain)
    
    w_Hb = MIN(1._dp, MAX(0._dp, (Hb_max - Hb_lo) / (Hb_hi - Hb_lo)))    
    r_crit = EXP( (w_Hb * lr_hi) + ((1._dp - w_Hb) * lr_lo))
    
    IF (contains_land .AND. dmax > r_crit * 2._dp * 1000._dp) THEN
      IsGood = .FALSE.
      RETURN      
    END IF

  END SUBROUTINE is_good_triangle
  
  ! == Mesh creation routines ==
  SUBROUTINE create_new_mesh( region)
    ! Called by all processes, but only the master can actually do stuff. Other processes
    ! only need to run the memory extension routines to make sure pointers are kept up to date.
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    ! Local variables
    CHARACTER(LEN=64), PARAMETER                  :: routine_name = 'create_new_mesh'
    INTEGER                                       :: n1, n2
    INTEGER                                       :: vi
    REAL(dp), DIMENSION(:    ), POINTER           ::  d2dx2,  d2dxdy,  d2dy2
    INTEGER                                       :: wd2dx2, wd2dxdy, wd2dy2
    INTEGER                                       :: orientation, it
    TYPE(type_mesh)                               :: submesh
    REAL(dp)                                      :: xmin, xmax, ymin, ymax
    REAL(dp)                                      :: res_min_inc
    CHARACTER(LEN=2)                              :: str_processid
    
    n1 = par%mem%n
    
    IF (par%master) WRITE(0,*) '  Creating a new mesh for region ', region%mesh%region_name, '...'
    
    ! Orientation of domain partitioning: east-west for GRL, north-south everywhere else
    IF (region%name == 'GRL') THEN
      orientation = 1
    ELSE
      orientation = 0
    END IF
    
    ! Make the ice margin 1 row of elements wider. This slightly increases the vertex count of
    ! the new mesh, but means the new mesh stays Fit much longer.
    CALL widen_high_res_zones( region%ice, region%mesh, region%time)
    
    ! Calculate surface curvature
    CALL allocate_shared_dp_1D( region%mesh%nV, d2dx2,  wd2dx2)
    CALL allocate_shared_dp_1D( region%mesh%nV, d2dxdy, wd2dxdy)
    CALL allocate_shared_dp_1D( region%mesh%nV, d2dy2,  wd2dy2)
    CALL d2dx2_a_to_a_2D(  region%mesh, region%ice%Hs_a, d2dx2 )
    CALL d2dxdy_a_to_a_2D( region%mesh, region%ice%Hs_a, d2dxdy)
    CALL d2dy2_a_to_a_2D(  region%mesh, region%ice%Hs_a, d2dy2 )
    DO vi = region%mesh%vi1, region%mesh%vi2
      region%ice%surf_curv( vi) = MAX(-1E-6, MIN(1E-6, SQRT(d2dx2( vi)**2 + d2dy2( vi)**2 + d2dxdy( vi)**2)))
    END DO
    CALL sync
    CALL deallocate_shared( wd2dx2)
    CALL deallocate_shared( wd2dxdy)
    CALL deallocate_shared( wd2dy2)
    
    ! Determine the domain of this process' submesh (based on distribution of vertices in
    ! the previous mesh, which works very well for workload balancing)
    IF (orientation == 0) THEN
!      CALL partition_domain_regular( region%mesh%xmin, region%mesh%xmax, par%i, par%n, xmin, xmax)
      CALL partition_domain_x_balanced( region%mesh, region%mesh%ymin, region%mesh%ymax, par%i, par%n, xmin, xmax)
      ymin = region%mesh%ymin
      ymax = region%mesh%ymax
    ELSE
!      CALL partition_domain_regular( region%mesh%ymin, region%mesh%ymax, par%i, par%n, ymin, ymax)
      CALL partition_domain_y_balanced( region%mesh, region%mesh%xmin, region%mesh%xmax, par%i, par%n, ymin, ymax)
      xmin = region%mesh%xmin
      xmax = region%mesh%xmax
    END IF
    
    ! Allocate memory, initialise a dummy mesh, refine it, and crop the memory
    CALL allocate_submesh_primary( submesh, region%name, 10, 20, C%nconmax)  
    CALL initialise_dummy_mesh(    submesh, xmin, xmax, ymin, ymax)
    CALL perturb_dummy_mesh(       submesh, 1 - region%mesh%perturb_dir)
    
    res_min_inc = C%res_max * 2._dp
    
    it = 0
    DO WHILE (res_min_inc > C%res_min)
      it = it + 1
      
      ! Increase resolution
      res_min_inc = res_min_inc / 2._dp
   
      ! Determine resolutions
      submesh%res_min          = MAX( C%res_min,          res_min_inc)
      submesh%res_max_margin   = MAX( C%res_max_margin,   res_min_inc)
      submesh%res_max_gl       = MAX( C%res_max_gl,       res_min_inc)
      submesh%res_max_cf       = MAX( C%res_max_cf,       res_min_inc)
      submesh%res_max_mountain = MAX( C%res_max_mountain, res_min_inc)
      submesh%res_max_coast    = MAX( C%res_max_coast,    res_min_inc)
      
      IF (debug_mesh_creation) WRITE(0,*) '  Process ', par%i, ' refining submesh to ', submesh%res_max_gl, ' km...'
      
      ! Refine the process submesh
      CALL refine_mesh( submesh, region%mesh, region%ice, region%refgeo_PD) 
      
      ! Align with neighbouring submeshes
      CALL align_all_submeshes( submesh, orientation)
      
      ! Split any new triangles (added during alignment) that are too sharp
      CALL refine_submesh_geo_only( submesh)
      
      ! Smooth the submesh using Lloyd' algorithm
      CALL Lloyds_algorithm_single_iteration_submesh( submesh)
      
      ! Write submesh to text file for debugging
      WRITE(str_processid,'(I2)') par%i;   str_processid = ADJUSTL(str_processid)
      IF (debug_mesh_creation) CALL write_mesh_to_text_file( submesh, 'submesh_proc_' // TRIM(str_processid) // '.txt')
      
      ! Check if everything went correctly
      CALL check_mesh( submesh)
    
    END DO
    
    ! Merge the process submeshes, create the final shared-memory mesh
    IF (debug_mesh_creation .AND. par%master) WRITE(0,*) '  Merging submeshes...'
    CALL merge_all_submeshes( submesh, orientation)
    IF (debug_mesh_creation .AND. par%master) WRITE(0,*) '  Creating final mesh...'
    CALL create_final_mesh_from_merged_submesh( submesh, region%mesh_new)

    IF (par%master) THEN
      WRITE(0,'(A)')                '   Finished creating final mesh.'
      WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh_new%nV
      WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh_new%nTri
      WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh_new%resolution_min/1000._dp, ' - ', region%mesh_new%resolution_max/1000._dp, ' km'
    END IF
    
    n2 = par%mem%n
    CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE create_new_mesh
  
  SUBROUTINE widen_high_res_zones( ice, mesh, time)
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
    
      mask_widened = ice%mask_margin_a
      
      DO vi = 1, mesh%nV
        IF (ice%mask_margin_a( vi)==1) THEN
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            D = NORM2( mesh%V( vi,:) - mesh%V( vc,:))
            IF (D < C%res_max_margin * 1000._dp * widening_factor) THEN
              mask_widened( vc) = 1
            END IF
          END DO
        END IF
      END DO
      
      ice%mask_margin_a = mask_widened
    
  ! Grounding line
  ! ==============
    
      mask_widened = ice%mask_gl_a
      
      DO vi = 1, mesh%nV
        IF (ice%mask_gl_a( vi)==1) THEN
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            D = NORM2( mesh%V( vi,:) - mesh%V( vc,:))
            IF (D < C%res_max_gl * 1000._dp * widening_factor) THEN
              mask_widened( vc) = 1
            END IF
          END DO
        END IF
      END DO
      
      ice%mask_gl_a = mask_widened
    
  ! Calving front
  ! =============
    
      mask_widened = ice%mask_cf_a
      
      DO vi = 1, mesh%nV
        IF (ice%mask_cf_a( vi)==1) THEN
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            D = NORM2( mesh%V( vi,:) - mesh%V(vc,:))
            IF (D < C%res_max_cf * 1000._dp * widening_factor) THEN
              mask_widened( vc) = 1
            END IF
          END DO
        END IF
      END DO
      
      ice%mask_cf_a = mask_widened
      
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE widen_high_res_zones
  
  ! == Extended Ruppert's algorithm
  SUBROUTINE refine_mesh( mesh, mesh_old, ice, refgeo_PD)
    ! Refine a mesh. Single-core, but multiple submeshes can be done in parallel on different cores.
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh_old
    TYPE(type_ice_model),       INTENT(IN)        :: ice
    TYPE(type_reference_geometry), INTENT(IN)     :: refgeo_PD
    
    ! Local variables
    INTEGER                                       :: ti
    LOGICAL                                       :: IsGood, FinishedRefining, DoExtendMemory
    REAL(dp), DIMENSION(2)                        :: p
    
    FinishedRefining = .FALSE.
    DoExtendMemory   = .FALSE.
    
    CALL extend_submesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
    
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
        CALL is_good_triangle(mesh, ti, mesh_old, ice, refgeo_PD, IsGood)
              
        IF (IsGood) THEN
          ! Remove this triangle from the stack
          mesh%RefMap(ti) = 0
          mesh%RefStack(mesh%RefStackN) = 0
          mesh%RefStackN = mesh%RefStackN - 1
        ELSE
          ! Spit this triangle, add the affected triangles to the stack
          p = mesh%Tricc(ti,:)
          CALL split_triangle( mesh, ti, p)
        END IF
        
        ! If we're reaching the memory limit, stop refining and extend the memory.
        IF (mesh%nV > mesh%nV_mem - 10) THEN
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
        CALL extend_submesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF
    
    END DO ! DO WHILE (.NOT. FinishedRefining)
  
  END SUBROUTINE refine_mesh
  
  ! == Determine if mesh updating is needed
  SUBROUTINE determine_mesh_fitness( mesh, ice, fitness)
    ! Determine how "fit" the current mesh is.
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    TYPE(type_ice_model),       INTENT(INOUT)     :: ice
    
    ! Output variables
    REAL(dp),                   INTENT(OUT)       :: fitness
    
    ! Local variables
    INTEGER                                       :: ti, i, status(MPI_STATUS_SIZE)
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
    
    DO ti = mesh%ti1, mesh%ti2
     
      ! Triangle vertex indices
      v1 = mesh%Tri( ti,1)
      v2 = mesh%Tri( ti,2)
      v3 = mesh%Tri( ti,3)
    
      ! Triangle vertex coordinates
      p = mesh%V( v1,:)
      q = mesh%V( v2,:)
      r = mesh%V( v3,:)

      ! Triangle legs
      pq = p-q
      qr = q-r
      rp = r-p
      
      ! Longest triangle leg
      dmax = MAXVAL([SQRT(pq(1)**2+pq(2)**2), SQRT(qr(1)**2+qr(2)**2), SQRT(rp(1)**2+rp(2)**2)])
      
      IF (ice%mask_coast_a( v1)==1 .OR. ice%mask_coast_a( v2)==1 .OR. ice%mask_coast_a( v3)==1) THEN
        ncoast = ncoast + 1
        lcoast = lcoast + dmax
        IF (dmax > C%res_max_coast*2.0_dp*1000._dp*res_tol) THEN
          nucoast = nucoast + 1
          lucoast = lucoast + dmax
        END IF
      END IF
      
      IF (ice%mask_margin_a( v1)==1 .OR. ice%mask_margin_a( v2)==1 .OR. ice%mask_margin_a( v3)==1) THEN
        nmargin = nmargin + 1
        lmargin = lmargin + dmax
        IF (dmax > C%res_max_margin*2.0_dp*1000._dp*res_tol) THEN
          numargin = numargin + 1
          lumargin = lumargin + dmax
        END IF
      END IF
      IF (ice%mask_gl_a( v1)==1 .OR. ice%mask_gl_a( v2)==1 .OR. ice%mask_gl_a( v3)==1) THEN
        ngl = ngl + 1
        lgl = lgl + dmax
        IF (dmax > C%res_max_gl*2.0_dp*1000._dp*res_tol) THEN
          nugl = nugl + 1
          lugl = lugl + dmax
        END IF
      END IF
      IF (ice%mask_cf_a( v1)==1 .OR. ice%mask_cf_a( v2)==1 .OR. ice%mask_cf_a( v3)==1) THEN
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
    
  END SUBROUTINE determine_mesh_fitness

END MODULE mesh_update_module
