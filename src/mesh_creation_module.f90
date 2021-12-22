MODULE mesh_creation_module

  ! Routines for creating the first mesh from a collection of forcing data on a Cartesian grid.

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
  USE parallel_module,                 ONLY: allocate_shared_dist_int_0D, allocate_shared_dist_dp_0D, &
                                             allocate_shared_dist_int_1D, allocate_shared_dist_dp_1D, &
                                             allocate_shared_dist_int_2D, allocate_shared_dist_dp_2D, &
                                             allocate_shared_dist_int_3D, allocate_shared_dist_dp_3D, &
                                             allocate_shared_dist_bool_1D
  USE data_types_module,               ONLY: type_model_region, type_mesh, type_reference_geometry
  USE mesh_help_functions_module,      ONLY: find_connection_widths, find_triangle_areas, find_Voronoi_cell_areas, get_lat_lon_coordinates, &
                                             determine_mesh_resolution, write_mesh_to_screen, merge_vertices, switch_vertices, redo_Tri_edge_indices, check_mesh, &
                                             cart_bilinear_dp, cart_bilinear_int, max_cart_over_triangle_int, max_cart_over_triangle_dp, cross2, &
                                             min_cart_over_triangle_int, sum_cart_over_triangle_dp, is_in_triangle, segment_intersection, is_walltowall, &
                                             find_POI_xy_coordinates, find_Voronoi_cell_geometric_centres, partition_domain_regular, find_POI_vertices_and_weights, &
                                             update_triangle_circumcenter, write_mesh_to_text_file, is_encroached_upon, is_boundary_segment, &
                                             calc_triangle_geometric_centres
  USE mesh_memory_module,              ONLY: allocate_mesh_primary, allocate_submesh_primary, extend_mesh_primary, extend_submesh_primary, &
                                             crop_mesh_primary, crop_submesh_primary, allocate_mesh_secondary, &
                                             deallocate_submesh_primary, move_data_from_submesh_to_mesh, share_submesh_access
  USE mesh_Delaunay_module,            ONLY: split_triangle, split_segment, flip_triangle_pairs, move_vertex
  USE mesh_operators_module,           ONLY: calc_matrix_operators_mesh
  USE mesh_ArakawaC_module,            ONLY: make_Ac_mesh
  USE mesh_five_colour_module,         ONLY: calculate_five_colouring_AaAc
  
  IMPLICIT NONE
  
  LOGICAL :: debug_mesh_creation = .FALSE.
  
  CONTAINS
  
  ! === Check whether or not a triangle meets all the fitness criteria.
  ! If you want to change the rules for mesh creation, this is where to do it.
  SUBROUTINE is_good_triangle( mesh, ti, refgeo_init, is_good)
    ! Check if triangle ti of the mesh is Good
    ! A triangle is not Good if:
    !   - its smallest internal angle is too small
    !   - its 2nd order surface deviation (=max(curvature)*typical_length) is too large
    !   - its area exceeds the limits based on ice velocity, grounding line or calving front
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                    INTENT(IN)        :: ti
    TYPE(type_reference_geometry),INTENT(IN)      :: refgeo_init
    LOGICAL,                    INTENT(OUT)       :: is_good

    INTEGER                                       :: vp,vq,vr
    REAL(dp), DIMENSION(2)                        :: p, q, r, POI
    REAL(dp), DIMENSION(2)                        :: pq, qr, rp
    LOGICAL                                       :: isso
    REAL(dp)                                      :: dmax
    INTEGER                                       :: n
    REAL(dp)                                      :: trisumel, mean_curvature, dz
    INTEGER                                       :: trinel
    REAL(dp), PARAMETER                           :: Hb_lo = 500._dp
    REAL(dp), PARAMETER                           :: Hb_hi = 1500._dp  
    REAL(dp)                                      :: Hb_max, w_Hb, lr_lo, lr_hi, r_crit  
    INTEGER                                       :: min_mask_int, max_mask_int
    REAL(dp)                                      :: mean_mask_dp
    LOGICAL                                       :: contains_ice, contains_nonice, contains_margin, contains_gl, contains_cf, contains_coast

    is_good = .TRUE.
    
    ! First check if the basic triangle geometry meets Ruppert's criteria
    ! ===================================================================
    
    CALL is_good_triangle_geo_only( mesh, ti, isso)
    IF (.NOT. is_good) RETURN

    ! Find length of longest triangle leg (for resolution checks)
    ! ===========================================================
    
    vp = mesh%Tri(ti,1)
    vq = mesh%Tri(ti,2)
    vr = mesh%Tri(ti,3)
    
    p  = mesh%V(vp,:)
    q  = mesh%V(vq,:)
    r  = mesh%V(vr,:) 

    ! Triangle legs
    pq = p-q
    qr = q-r
    rp = r-p

    ! Longest leg
    dmax = MAXVAL([NORM2(pq), NORM2(qr), NORM2(rp)])
    
    ! Coarsest allowed resolution
    ! ===========================
    
    IF (dmax > mesh%res_max * 1.5_dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF
    
    ! Finest allowed resolution
    ! =========================
    
    IF (dmax < mesh%res_min * 1.5_dp * 1000._dp) THEN
      is_good = .TRUE.
      RETURN
    END IF
    
    ! Resolution at points of interest
    ! ================================
    
    DO n = 1, mesh%nPOI
      POI = mesh%POI_XY_coordinates(n,:)
      IF (is_in_triangle( p, q, r, POI) .AND. dmax > mesh%POI_resolutions(n) * 1.5_dp * 1000._dp) THEN
        is_good = .FALSE.
        RETURN
      END IF
    END DO

    ! Determine what's inside the triangle
    ! ====================================
    
    contains_ice    = .FALSE.
    contains_nonice = .FALSE.
    contains_margin = .FALSE.
    contains_gl     = .FALSE.
    contains_cf     = .FALSE.
    contains_coast  = .FALSE.
    
    CALL min_cart_over_triangle_int( p, q, r, refgeo_init%mask_ice, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, min_mask_int, trinel)
    CALL max_cart_over_triangle_int( p, q, r, refgeo_init%mask_ice, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, max_mask_int, trinel)    
    IF (trinel>0) THEN
      IF (max_mask_int==1) contains_ice    = .TRUE.
      IF (min_mask_int==0) contains_nonice = .TRUE.
    ELSE
      CALL cart_bilinear_int( refgeo_init%mask_ice, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, mean_mask_dp)
      IF (mean_mask_dp>0.1_dp) contains_ice    = .TRUE.
      IF (mean_mask_dp<0.9_dp) contains_nonice = .TRUE.
    END IF
    IF (contains_ice .AND. contains_nonice) contains_margin = .TRUE.
    
    CALL max_cart_over_triangle_int(p,q,r, refgeo_init%mask_gl, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, max_mask_int, trinel)    
    IF (trinel>0) THEN
      IF (max_mask_int==1) contains_gl = .TRUE.
    ELSE
      CALL cart_bilinear_int( refgeo_init%mask_gl, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, mean_mask_dp)
      IF (mean_mask_dp>0.1_dp) contains_gl = .TRUE.
    END IF
    
    CALL max_cart_over_triangle_int(p,q,r, refgeo_init%mask_cf, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, max_mask_int, trinel)    
    IF (trinel>0) THEN
      IF (max_mask_int==1) contains_cf = .TRUE.
    ELSE
      CALL cart_bilinear_int( refgeo_init%mask_cf, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, mean_mask_dp)
      IF (mean_mask_dp>0.1_dp) contains_cf = .TRUE.
    END IF
    
    CALL max_cart_over_triangle_int(p,q,r, refgeo_init%mask_coast, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, max_mask_int, trinel)    
    IF (trinel>0) THEN
      IF (max_mask_int==1) contains_coast = .TRUE.
    ELSE
      CALL cart_bilinear_int( refgeo_init%mask_coast, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, mean_mask_dp)
      IF (mean_mask_dp>0.1_dp) contains_coast = .TRUE.
    END IF

    ! Second-order surface deviation (curvature times size)
    ! =====================================================
    
    CALL sum_cart_over_triangle_dp(p,q,r, refgeo_init%surf_curv, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, trisumel, trinel)
    IF (trinel>4) THEN
      mean_curvature = trisumel / trinel
    ELSE
      CALL cart_bilinear_dp( refgeo_init%surf_curv, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, mean_curvature)
    END IF
    dz = 0.5_dp * mean_curvature * dmax**2

    IF (contains_ice .AND. dz > mesh%dz_max_ice) THEN
      is_good = .FALSE.
     RETURN
    END IF

    ! Special area resolution - ice margin, grounding line, calving front
    ! ===================================================================
    
    ! Coastline
    IF (contains_coast .AND. dmax > mesh%res_max_coast * 2._dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF
    
    ! Ice margin
    IF (contains_margin .AND. dmax > mesh%res_max_margin * 2._dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF
    
    ! Grounding line
    IF (contains_gl .AND. dmax > mesh%res_max_gl * 2._dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF

    ! Calving front
    IF (contains_cf .AND. dmax > mesh%res_max_cf * 2._dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN
    END IF
    
    ! Ice-free bed topography (higher res for mountains so inception is captured better)
    ! ==================================================================================
    
    CALL max_cart_over_triangle_dp(p,q,r, refgeo_init%Hb_grid, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, Hb_max,trinel)
    IF (trinel==0) THEN
      CALL cart_bilinear_dp( refgeo_init%Hb_grid, refgeo_init%grid%x, refgeo_init%grid%y, refgeo_init%grid%nx, refgeo_init%grid%ny, (p+q+r)/3._dp, Hb_max)
    END IF
    
    lr_lo = LOG( mesh%res_max)
    lr_hi = LOG( mesh%res_max_mountain)
    
    w_Hb = MIN(1._dp, MAX(0._dp, (Hb_max - Hb_lo) / (Hb_hi - Hb_lo)))    
    r_crit = EXP( (w_Hb * lr_hi) + ((1._dp - w_Hb) * lr_lo))
    
    IF (contains_nonice .AND. dmax > r_crit * 2._dp * 1000._dp) THEN
      is_good = .FALSE.
      RETURN      
    END IF

  END SUBROUTINE is_good_triangle
  SUBROUTINE is_good_triangle_geo_only( mesh, ti, is_good)
    ! Check if triangle ti of the mesh is Good
    ! A triangle is not Good if:
    !   - its smallest internal angle is too small
    !   - its 2nd order surface deviation (=max(curvature)*typical_length) is too large
    !   - its area exceeds the limits based on ice velocity, grounding line or calving front
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                    INTENT(IN)        :: ti

    LOGICAL,                    INTENT(OUT)       :: is_good

    INTEGER                                       :: vp,vq,vr
    REAL(dp), DIMENSION(2)                        :: p,q,r
    REAL(dp), DIMENSION(2)                        :: pq,qr,rp
    REAL(dp)                                      :: ap,aq,ar,alpha
    LOGICAL                                       :: isso
    
    is_good = .TRUE.

    ! Triangle geometry (the basis of the original version of Rupperts Algorithm)
    ! ===========================================================================
    
    vp = mesh%Tri(ti,1)
    vq = mesh%Tri(ti,2)
    vr = mesh%Tri(ti,3)
    
    p  = mesh%V(vp,:)
    q  = mesh%V(vq,:)
    r  = mesh%V(vr,:)  

    ! Triangle legs
    pq = p-q
    qr = q-r
    rp = r-p

    ! Internal angles
    ap = ACOS(-(rp(1)*pq(1) + rp(2)*pq(2))/(NORM2(rp)*NORM2(pq)))
    aq = ACOS(-(pq(1)*qr(1) + pq(2)*qr(2))/(NORM2(pq)*NORM2(qr)))
    ar = ACOS(-(rp(1)*qr(1) + rp(2)*qr(2))/(NORM2(rp)*NORM2(qr)))

    ! Smallest internal angle
    alpha = MINVAL([ap, aq, ar])

    IF (alpha < mesh%alpha_min) THEN
      is_good = .FALSE.
      RETURN
    END IF
    
    ! If its an edge triangle, check if the third vertex encroaches on the edge segment
    IF (is_boundary_segment( mesh, vp, vq)) THEN
      CALL is_encroached_upon( mesh, vp, vq, isso)
      IF (isso) THEN
        is_good = .FALSE.
        RETURN
      END IF
    ELSEIF (is_boundary_segment( mesh, vq, vr)) THEN
      CALL is_encroached_upon( mesh, vq, vr, isso)
      IF (isso) THEN
        is_good = .FALSE.
        RETURN
      END IF
    ELSEIF (is_boundary_segment( mesh, vr, vp)) THEN
      CALL is_encroached_upon( mesh, vr, vp, isso)
      IF (isso) THEN
        is_good = .FALSE.
        RETURN
      END IF
    END IF
    
    ! Forbid "wall to wall" triangles (i.e. triangles with vertices lying on opposite domain boundaries)
    IF (is_walltowall( mesh, ti)) THEN
      is_good = .FALSE.
      RETURN
    END IF
    
  END SUBROUTINE is_good_triangle_geo_only
    
  ! == Mesh creation routines ==
  SUBROUTINE create_mesh_from_cart_data( region)
    ! Create the first mesh, using the data from the initial file to force the resolution.
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    ! Local variables
    CHARACTER(LEN=64), PARAMETER                  :: routine_name = 'create_mesh_from_cart_data'
    INTEGER                                       :: n1, n2
    INTEGER                                       :: orientation
    TYPE(type_mesh)                               :: submesh
    REAL(dp)                                      :: xmin, xmax, ymin, ymax
    REAL(dp)                                      :: res_min_inc
    CHARACTER(LEN=2)                              :: str_processid
    
    n1 = par%mem%n
    
    IF (par%master) WRITE(0,*) '  Creating the first mesh...'
    
    ! Orientation of domain partitioning: east-west for GRL, north-south everywhere else
    IF (region%name == 'GRL') THEN
      orientation = 1
    ELSE
      orientation = 0
    END IF
        
    ! Determine the domain of this process' submesh.
    IF (orientation == 0) THEN
      CALL partition_domain_regular( MINVAL(region%refgeo_init%grid%x), MAXVAL(region%refgeo_init%grid%x), par%i, par%n, xmin, xmax)
      ymin = MINVAL(region%refgeo_init%grid%y)
      ymax = MAXVAL(region%refgeo_init%grid%y)
    ELSE
      CALL partition_domain_regular( MINVAL(region%refgeo_init%grid%y), MAXVAL(region%refgeo_init%grid%y), par%i, par%n, ymin, ymax)
      xmin = MINVAL(region%refgeo_init%grid%x)
      xmax = MAXVAL(region%refgeo_init%grid%x)
    END IF
    
    ! Allocate memory and initialise a dummy mesh
    CALL allocate_submesh_primary( submesh, region%name, 10, 20, C%nconmax)    
    CALL initialise_dummy_mesh(    submesh, xmin, xmax, ymin, ymax)
    CALL perturb_dummy_mesh(       submesh, 0)
    
    ! Exception for the EISMINT schematic tests, which start with a flat bedrock and no ice, leading to a very coarse mesh.
    ! After the first time step, a thin ice sheet forms in the area with positive SMB, so that the margin lies in an area
    ! with very coarse resolution. This triggers a mesh update, which results in a mesh with a very wide high-resolution band.
    ! Prevent this by making the very first mesh slightly finer than strictly needed.
    IF (C%do_benchmark_experiment) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6' ) THEN
        submesh%res_max = MAX(C%res_min * 2._dp, 16._dp)
      ELSEIF (C%choice_benchmark_experiment == 'Bueler'.OR. &
              C%choice_benchmark_experiment == 'Halfar') THEN
        ! No need for an exception here, as this one starts with a (small) ice sheet
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
              C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        ! No need for an exception here either, as the initial state has a coastline.
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in create_mesh_from_cart_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF
        
    ! Refine the process submesh with incrementally increasing resolution, aligning with neighbouring
    ! submeshes after every step, until we reach the final desired resolution
    
    res_min_inc = C%res_max * 2._dp
    
    DO WHILE (res_min_inc > C%res_min)
      
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
      CALL refine_mesh( submesh, region%refgeo_init)
      
      ! Align with neighbouring submeshes
      CALL align_all_submeshes( submesh, orientation)
      
      ! Split any new triangles (added during alignment) that are too sharp
      CALL refine_submesh_geo_only( submesh)
      
      ! Smooth the submesh using Lloyd' algorithm
      CALL Lloyds_algorithm_single_iteration_submesh( submesh)
      
      ! After the last refinement step, apply Lloyds algorithm two more times, because we can.
      IF (res_min_inc <= C%res_min) THEN
      CALL Lloyds_algorithm_single_iteration_submesh( submesh)
      CALL Lloyds_algorithm_single_iteration_submesh( submesh)
      END IF
      
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
    CALL create_final_mesh_from_merged_submesh( submesh, region%mesh)

    IF (par%master) THEN
      WRITE(0,'(A)')                '   Finished creating final mesh.'
      WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh%nV
      WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh%nTri
      WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh%resolution_min/1000._dp, ' - ', region%mesh%resolution_max/1000._dp, ' km'
    END IF
    
    n2 = par%mem%n
    CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE create_mesh_from_cart_data
  
  ! == Extended and basic Ruppert's algorithm
  SUBROUTINE refine_mesh( mesh, refgeo_init)
    ! Refine a mesh. Single-core, but multiple submeshes can be done in parallel on different cores.
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    TYPE(type_reference_geometry),INTENT(IN)      :: refgeo_init
    
    ! Local variables
    INTEGER                                       :: ti
    REAL(dp), DIMENSION(2)                        :: p
    LOGICAL                                       :: IsGood, FinishedRefining, DoExtendMemory
    
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
    
    DO WHILE (.NOT. FinishedRefining)
    
      ! Refine the mesh until it's done, or until it's almost out of memory.
      ! ====================================================================
      
      DoExtendMemory = .FALSE.
          
      DO WHILE (mesh%RefStackN > 0)
    
        ! Check the last triangle list in the RefineStack. If it's
        ! Bad, refine it and add the affected triangles to the RefineStack.
      
        ti = mesh%RefStack( mesh%RefStackN)       
        CALL is_good_triangle( mesh, ti, refgeo_init, IsGood)
              
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
  SUBROUTINE refine_mesh_geo_only( mesh)
    ! Refine a mesh, using only triangle geometry as a condition (so really the original version of Ruppert's algorithm)
    ! Meant to be run on the final, merged mesh - called by all processes, but the work is only done by the Master.
    ! Must be called by all to be able to call ExtendMeshMemory if necessary.
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: ti
    REAL(dp), DIMENSION(2)                        :: p
    LOGICAL                                       :: IsGood, FinishedRefining, DoExtendMemory
    
    ! List all triangles for checking
    IF (par%master) THEN
      mesh%RefMap    = 0
      mesh%RefStack  = 0
      mesh%RefStackN = 0
      DO ti = 1, mesh%nTri
        mesh%RefMap(ti)               = 1
        mesh%RefStackN                = mesh%RefStackN + 1
        mesh%RefStack(mesh%RefStackN) = ti
      END DO
    END IF ! IF (par%master) THEN
    CALL sync
    
    FinishedRefining = .FALSE.
    
    DO WHILE (.NOT. FinishedRefining)
    
      ! Refine the mesh until it's done, or until it's almost out of memory.
      ! ====================================================================
      
      DoExtendMemory = .FALSE.
          
      IF (par%master) THEN
      DO WHILE (mesh%RefStackN > 0)
    
        ! Check the last triangle list in the RefineStack. If it's
        ! Bad, refine it and add the affected triangles to the RefineStack.
      
        ti = mesh%RefStack( mesh%RefStackN)       
        CALL is_good_triangle_geo_only( mesh, ti, IsGood)
              
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
      END IF ! IF (par%master)
      
      FinishedRefining = .FALSE.
      IF (mesh%RefStackN == 0) FinishedRefining = .TRUE.
      
      ! Check if all processes finished refining. If so, exit.
      CALL MPI_BCAST( FinishedRefining, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)      
      IF (FinishedRefining) EXIT
      
      ! Check if any process needs to extend their memory.
      CALL MPI_BCAST( DoExtendMemory, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      IF (DoExtendMemory) THEN
        ! By extending the memory to mesh%nV + 1000, we ensure that processes that 
        ! have already finished refining do not keep adding useless extra memory.
        CALL extend_mesh_primary( mesh, mesh%nV + 1000, mesh%nTri + 2000)
      END IF
    
    END DO !DO WHILE (.NOT. FinishedRefining)
  
  END SUBROUTINE refine_mesh_geo_only
  SUBROUTINE refine_submesh_geo_only( mesh)
    ! Refine a mesh. Single-core, but multiple submeshes can be done in parallel on different cores.
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: ti
    REAL(dp), DIMENSION(2)                        :: p
    LOGICAL                                       :: IsGood, FinishedRefining, DoExtendMemory
    
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
    
    DO WHILE (.NOT. FinishedRefining)
    
      ! Refine the mesh until it's done, or until it's almost out of memory.
      ! ====================================================================
      
      DoExtendMemory = .FALSE.
          
      DO WHILE (mesh%RefStackN > 0)
    
        ! Check the last triangle list in the RefineStack. If it's
        ! Bad, refine it and add the affected triangles to the RefineStack.
      
        ti = mesh%RefStack( mesh%RefStackN)         
        CALL is_good_triangle_geo_only( mesh, ti, IsGood)
              
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
  
  END SUBROUTINE refine_submesh_geo_only
  
  ! == Lloyd's algorithm for "smoothing" a (sub)mesh
  SUBROUTINE Lloyds_algorithm_single_iteration_submesh( mesh)
    ! Lloyd's algorithm: move all vertices to the geometric centers of their Voronoi cells, and update the triangulation.
    ! This "smooths" the mesh, reducing resolution gradients and widening internal angles, thus making it more
    ! suitable for numerical methods (particularly the SSA).
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: vi, ci, cip1
    REAL(dp)                                      :: VorTriA, sumVorTriA
    REAL(dp), DIMENSION(2)                        :: pa, pb, pc, VorGC
    
    ! Move all non-boundary vertices to their Voronoi cell geometric centre
    DO vi = 1, mesh%nV
    
      ! Leave boundary vertices where they are
      IF (mesh%edge_index( vi) > 0) CYCLE
        
      ! Find the geometric centre of this vertex' Voronoi cell
      VorGC      = 0._dp
      sumVorTriA = 0._dp
      
      DO ci = 1, mesh%nC( vi)
      
        cip1 = ci + 1
        IF (cip1 > mesh%nC( vi)) cip1 = 1
        
        pa = mesh%V( vi,:)
        pb = mesh%V( mesh%C( vi,ci  ),:)
        pc = mesh%V( mesh%C( vi,cip1),:)
        
        VorTriA = cross2( pb - pa, pc - pa)
        
        VorGC = VorGC + VorTriA * (pa + pb + pc) / 3._dp
        sumVorTriA   = sumVorTriA   + VorTriA
        
      END DO ! DO ci = 1, mesh%nC( vi)
      
      VorGC = VorGC / sumVorTriA
        
      ! Move the vertex
      CALL move_vertex( mesh, vi, VorGC)
    
    END DO
    CALL sync
    
  END SUBROUTINE Lloyds_algorithm_single_iteration_submesh
  
  ! == Align and merge submeshes created by parallel processes
  SUBROUTINE align_all_submeshes( submesh, orientation)
    ! Align (= ensure all seam vertices are shared) all adjacent submeshes.
  
    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: submesh
    INTEGER,                    INTENT(IN)        :: orientation
    
    INTEGER, DIMENSION(:,:  ), ALLOCATABLE        :: alignlist
    INTEGER                                       :: nalign, i_left, i_right, i, nVl_east, nVr_west
    
    ! No need for this if we're running on a single core
    IF (par%n == 1) RETURN
    
    ! Since each submesh can only be aligned with one neighbour at a time, and each one has
    ! at most two neighbours, we need two passes.
    
    ALLOCATE( alignlist( par%n, 2))
    
    ! == Pass one: even
    ! =================
  
    alignlist = 0
    nalign    = 0
    
    DO i = 0, par%n-2, 2
        
      nalign = nalign + 1
      alignlist( nalign,:) = [i, i+1]
      
    END DO
        
    ! Determine if we're participating as Left, Right or Passive
    ! (Passive still needs to run through the code, to participate in ExtendSubmesh calls)    
    i_left  = -1
    i_right = -1
    DO i = 1, nalign
      IF ( alignlist(i,1) == par%i .OR. alignlist(i,2) == par%i) THEN
        i_left  = alignlist(i,1)
        i_right = alignlist(i,2)
      END IF
    END DO
    
    IF (nalign > 0) THEN
      CALL align_submeshes( submesh, orientation, i_left, i_right, nVl_east, nVr_west)
    END IF
    
  ! == Pass two: odd
  ! ================
  
    alignlist = 0
    nalign    = 0
    
    DO i = 1, par%n-2, 2
      
      nalign = nalign + 1
      alignlist(nalign,:) = [i, i+1]
      
    END DO
        
    ! Determine if we're participating as Left, Right or Passive
    ! (Passive still needs to run through the code, to participate in ExtendSubmesh calls)    
    i_left  = -1
    i_right = -1
    DO i = 1, nalign
      IF ( alignlist(i,1) == par%i .OR. alignlist(i,2) == par%i) THEN
        i_left  = alignlist(i,1)
        i_right = alignlist(i,2)
      END IF
    END DO
    
    IF (nalign > 0) THEN
      CALL align_submeshes( submesh, orientation, i_left, i_right, nVl_east, nVr_west)
    END IF
    
  ! Clean up after yourself
  ! =======================
    
    DEALLOCATE( alignlist)
    
  END SUBROUTINE align_all_submeshes
  SUBROUTINE align_submeshes( submesh, orientation, p_left, p_right, nVl_east, nVr_west)
    ! Align: make sure two adjacent submeshes share all their boundary vertices
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: submesh
    INTEGER,                    INTENT(IN)        :: orientation  ! 0 = eastwest, 1 = northsouth
    INTEGER,                    INTENT(IN)        :: p_left, p_right
    INTEGER,                    INTENT(OUT)       :: nVl_east, nVr_west
    
    INTEGER                                       :: status(MPI_STATUS_SIZE)
    
    INTEGER                                       :: nVl_tot, nVr_tot, nV_extra
    INTEGER,  DIMENSION(:  ), ALLOCATABLE         :: Vil_east, Vir_west, Vi_dummy
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Vl_east, Vr_west, V_dummy
    INTEGER                                       :: vi_prev, vi
    LOGICAL                                       :: FoundNorth
    
    INTEGER,  DIMENSION(:,:), ALLOCATABLE         :: T
    INTEGER                                       :: nT
    INTEGER                                       :: vil, vir
    LOGICAL,  DIMENSION(:  ), ALLOCATABLE         :: FoundAsMatchl, FoundAsMatchr
    REAL(dp), DIMENSION(2  )                      :: p_new, pc, pu, pl
    INTEGER                                       :: vilc, vilu, vill, virc, viru, virl
    
    INTEGER                                       :: vi2
    
!    IF (par%i == p_left .OR. par%i == p_right) THEN
!      WRITE(0,'(A,I2,A,I2,A,I2)') '  align_submeshes - process ', par%i, ': aligning mesh ', p_left, ' with mesh ', p_right
!    ELSE
!      WRITE(0,'(A,I2,A)')         '  align_submeshes - process ', par%i, ': passively running align_submeshes'
!    END IF
!    CALL sync
        
! Create the list of boundary vertices
! ====================================

    ! Allocate dummy memory in non-participating processes (to prevent compiler warnings only)
    IF (.NOT. (par%i == p_left .OR. par%i == p_right)) THEN
      ALLOCATE( FoundAsMatchl(1))
      ALLOCATE( FoundAsMatchr(1))
      ALLOCATE( T(          1,1))
      ALLOCATE( Vil_east(     1))
      ALLOCATE( Vl_east(    1,1))
      ALLOCATE( Vir_west(     1))
      ALLOCATE( Vr_west(    1,1))
      nT = 0
    END IF
    
    IF (par%i == p_left) THEN
      
      ALLOCATE( Vi_dummy( submesh%nV   ))
      ALLOCATE( V_dummy(  submesh%nV, 2))
      
      IF (orientation == 0) THEN
      
        Vi_dummy(1)  = 2
        V_dummy(1,:) = submesh%V(2,:)
        nVl_east     = 1
        nVl_tot      = submesh%nV
        
        FoundNorth = .FALSE.
        vi_prev    = 2
        
        DO WHILE (.NOT. FoundNorth)
          vi = submesh%C(vi_prev,1)
          nVl_east = nVl_east + 1
          Vi_dummy( nVl_east  ) = vi
          V_dummy(  nVl_east,:) = submesh%V(vi,:)
          vi_prev = vi
          IF ( vi == 3) FoundNorth = .TRUE.
        END DO ! DO WHILE (.NOT. FoundNorth)
      
      ELSEIF (orientation == 1) THEN
      
        Vi_dummy(1)  = 3
        V_dummy(1,:) = submesh%V(3,:)
        nVl_east     = 1
        nVl_tot      = submesh%nV
        
        FoundNorth = .FALSE.
        vi_prev    = 3
        
        DO WHILE (.NOT. FoundNorth)
          vi = submesh%C(vi_prev,1)
          nVl_east = nVl_east + 1
          Vi_dummy( nVl_east  ) = vi
          V_dummy(  nVl_east,:) = submesh%V(vi,:)
          vi_prev = vi
          IF ( vi == 4) FoundNorth = .TRUE.
        END DO ! DO WHILE (.NOT. FoundNorth)
        
      END IF
      
      ALLOCATE( Vil_east( nVl_east   ))
      ALLOCATE( Vl_east(  nVl_east, 2))
      Vil_east = Vi_dummy( 1:nVl_east  )
      Vl_east  = V_dummy(  1:nVl_east,:)
      DEALLOCATE( Vi_dummy)
      DEALLOCATE( V_dummy)
      
    ELSEIF (par%i == p_right) THEN
      
      ALLOCATE( Vi_dummy( submesh%nV   ))
      ALLOCATE( V_dummy(  submesh%nV, 2))
      
      IF (orientation == 0) THEN
      
        Vi_dummy(1)  = 1
        V_dummy(1,:) = submesh%V(1,:)
        nVr_west     = 1
        nVr_tot      = submesh%nV
        
        FoundNorth = .FALSE.
        vi_prev    = 1
        
        DO WHILE (.NOT. FoundNorth)
          vi = submesh%C(vi_prev, submesh%nC(vi_prev))
          nVr_west = nVr_west + 1
          Vi_dummy( nVr_west  ) = vi
          V_dummy(  nVr_west,:) = submesh%V(vi,:)
          vi_prev = vi
          IF ( vi == 4) FoundNorth = .TRUE.
        END DO ! DO WHILE (.NOT. FoundNorth)
      
      ELSEIF (orientation == 1) THEN
      
        Vi_dummy(1)  = 2
        V_dummy(1,:) = submesh%V(2,:)
        nVr_west     = 1
        nVr_tot      = submesh%nV
        
        FoundNorth = .FALSE.
        vi_prev    = 2
        
        DO WHILE (.NOT. FoundNorth)
          vi = submesh%C(vi_prev, submesh%nC(vi_prev))
          nVr_west = nVr_west + 1
          Vi_dummy( nVr_west  ) = vi
          V_dummy(  nVr_west,:) = submesh%V(vi,:)
          vi_prev = vi
          IF ( vi == 1) FoundNorth = .TRUE.
        END DO ! DO WHILE (.NOT. FoundNorth)
      
      END IF
      
      ALLOCATE( Vir_west( nVr_west   ))
      ALLOCATE( Vr_west(  nVr_west, 2))
      Vir_west = Vi_dummy( 1:nVr_west  )
      Vr_west  = V_dummy(  1:nVr_west,:)
      DEALLOCATE( Vi_dummy)
      DEALLOCATE( V_dummy)
      
    END IF ! IF (par%i == p_left) THEN
    
! Exchange lists of boundary vertices
! ===================================
    
    IF (par%i == p_left) THEN
    
      ! Receive list of boundary vertices from p_right
      CALL MPI_RECV( nVr_west,  1,          MPI_INTEGER,          p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr) 
      CALL MPI_RECV( nVr_tot,   1,          MPI_INTEGER,          p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)      
      ALLOCATE( Vir_west( nVr_west   ))
      ALLOCATE( Vr_west(  nVr_west, 2))
      CALL MPI_RECV(  Vir_west, nVr_west,   MPI_INTEGER,          p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      CALL MPI_RECV(  Vr_west,  nVr_west*2, MPI_DOUBLE_PRECISION, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
    
      ! Send list of boundary vertices to p_left
      CALL MPI_SEND( nVl_east, 1,          MPI_INTEGER,          p_right, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_SEND( nVl_tot,  1,          MPI_INTEGER,          p_right, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_SEND( Vil_east, nVl_east,   MPI_INTEGER,          p_right, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_SEND( Vl_east,  nVl_east*2, MPI_DOUBLE_PRECISION, p_right, 0, MPI_COMM_WORLD, ierr)
      
    ELSEIF (par%i == p_right) THEN
    
      ! Send list of boundary vertices to p_left
      CALL MPI_SEND( nVr_west, 1,          MPI_INTEGER,          p_left, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_SEND( nVr_tot,  1,          MPI_INTEGER,          p_left, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_SEND( Vir_west, nVr_west,   MPI_INTEGER,          p_left, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_SEND( Vr_west,  nVr_west*2, MPI_DOUBLE_PRECISION, p_left, 0, MPI_COMM_WORLD, ierr)
      
      ! Receive list of boundary vertices from p_left
      CALL MPI_RECV( nVl_east,  1,          MPI_INTEGER,          p_left, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr) 
      CALL MPI_RECV( nVl_tot,   1,          MPI_INTEGER,          p_left, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)      
      ALLOCATE( Vil_east( nVl_east   ))
      ALLOCATE( Vl_east(  nVl_east, 2))
      CALL MPI_RECV(  Vil_east, nVl_east,   MPI_INTEGER,          p_left, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      CALL MPI_RECV(  Vl_east,  nVl_east*2, MPI_DOUBLE_PRECISION, p_left, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      
    END IF ! IF (par%i == p_left) THEN
    
! Create the boundary translation table T
! =======================================
    
    IF (par%i == p_left .OR. par%i == p_right) THEN
      
      ! Create T
      ALLOCATE( T( nVl_east + nVr_west, 2))
      ALLOCATE( FoundAsMatchl( nVl_east))
      ALLOCATE( FoundAsMatchr( nVr_west))
      
      IF (orientation == 0) THEN
      
        T(1,:) = [2, 1]
        T(2,:) = [3, 4]
        nT = 2
        
        FoundAsMatchr = .FALSE.
        FoundAsMatchl = .FALSE.
        
        DO vi = 2, nVl_east - 1
          vil = Vil_east( vi)
          DO vi2 = 2, nVr_west - 1
            IF ( FoundAsMatchr( vi2)) CYCLE
            vir = Vir_west( vi2)
            IF ( ABS( Vl_east( vi,2) - Vr_west( vi2,2)) < submesh%tol_dist) THEN
              FoundAsMatchl(vi)  = .TRUE.
              FoundAsMatchr(vi2) = .TRUE.
              nT = nT+1
              T(nT,:) = [vil, vir]
              EXIT
            END IF
          END DO
        END DO
        DO vi2 = 2, nVr_west - 1
          vir = Vir_west(vi2)
          DO vi = 2, nVl_east - 1
            IF (FoundAsMatchl( vi)) CYCLE
            vil = Vil_east(vi)
            IF ( ABS( Vl_east( vi,2) - Vr_west( vi2,2)) < submesh%tol_dist) THEN
              FoundAsMatchl(vi)  = .TRUE.
              FoundAsMatchr(vi2) = .TRUE.
              nT = nT+1
              T(nT,:) = [vil, vir]
              EXIT
            END IF
          END DO
        END DO
      
      ELSEIF (orientation == 1) THEN
      
        T(1,:) = [3, 2]
        T(2,:) = [4, 1]
        nT = 2
        
        FoundAsMatchr = .FALSE.
        FoundAsMatchl = .FALSE.
                
        DO vi = 2, nVl_east - 1
          vil = Vil_east( vi)
          DO vi2 = 2, nVr_west - 1
            IF ( FoundAsMatchr( vi2)) CYCLE
            vir = Vir_west( vi2)
            IF ( ABS( Vl_east( vi,1) - Vr_west( vi2,1)) < submesh%tol_dist) THEN
              FoundAsMatchl(vi)  = .TRUE.
              FoundAsMatchr(vi2) = .TRUE.
              nT = nT+1
              T(nT,:) = [vil, vir]
              EXIT
            END IF
          END DO
        END DO
        DO vi2 = 2, nVr_west - 1
          vir = Vir_west(vi2)
          DO vi = 2, nVl_east - 1
            IF (FoundAsMatchl( vi)) CYCLE
            vil = Vil_east(vi)
            IF ( ABS( Vl_east( vi,1) - Vr_west( vi2,1)) < submesh%tol_dist) THEN
              FoundAsMatchl(vi)  = .TRUE.
              FoundAsMatchr(vi2) = .TRUE.
              nT = nT+1
              T(nT,:) = [vil, vir]
              EXIT
            END IF
          END DO
        END DO
      
      END IF ! IF (orientation == 0) THEN
      
    END IF ! IF (par%i == p_left .OR. par%i == p_right) THEN
    
! Extend submesh memory to accomodate the extra vertices
! ======================================================
    
    IF (par%i == p_left) THEN    
      nV_extra = nVr_west - nT
    ELSEIF (par%i == p_right) THEN
      nV_extra = nVl_east - nT
    ELSE
      nV_extra = 0
    END IF
    
    CALL extend_submesh_primary( submesh, submesh%nV + nV_extra, submesh%nTri + 2 * nV_extra) 
    
! Add the extra vertices to the submesh
! =====================================
    
    IF (par%i == p_left) THEN
      
      DO vi2 = 2, nVr_west-1
        IF (FoundAsMatchr(vi2)) CYCLE

        vir = Vir_west( vi2)
        vil = submesh%nV+1

        ! Find the vertices spanning the segment that needs to be split
        p_new = Vr_west( vi2,:)
        
        DO vi = 1, nVl_east

          vilc = Vil_east( vi)
          vilu = submesh%C( vilc, 1)
          vill = submesh%C( vilc, submesh%nC( vilc))

          pc   = submesh%V(vilc,:)
          pu   = submesh%V(vilu,:)
          pl   = submesh%V(vill,:)
          
          IF (orientation == 0) THEN

            IF ( pu(2) > p_new(2) .AND. pc(2) < p_new(2)) THEN
              nT = nT+1
              T(nT,:) = [vil, vir]
              CALL split_segment( submesh, vilc, vilu, p_new)
              EXIT
            ELSEIF ( pc(2) > p_new(2) .AND. pl(2) < p_new(2)) THEN
              nT = nT+1
              T(nT,:) = [vil, vir]
              CALL split_segment( submesh, vill, vilc, p_new)
              EXIT
            END IF
          
          ELSEIF (orientation == 1) THEN

            IF ( pu(1) < p_new(1) .AND. pc(1) > p_new(1)) THEN
              nT = nT+1
              T(nT,:) = [vil, vir]
              CALL split_segment( submesh, vilc, vilu, p_new)
              EXIT
            ELSEIF ( pc(1) < p_new(1) .AND. pl(1) > p_new(1)) THEN
              nT = nT+1
              T(nT,:) = [vil, vir]
              CALL split_segment( submesh, vill, vilc, p_new)
              EXIT
            END IF
          
          END IF ! IF (orientation == 0) THEN
                    
        END DO ! DO vi = 1:nVl_east
      END DO ! DO vi2 = 1, nVr_west
      
    ELSEIF (par%i == p_right) THEN
      
      DO vi = 2, nVl_east-1
        IF (FoundAsMatchl( vi)) CYCLE

        vil = Vil_east( vi)
        vir = submesh%nV+1

        ! Find the vertices spanning the segment that needs to be split
        p_new = Vl_east( vi,:)
        
        DO vi2 = 1, nVr_west

          virc = Vir_west( vi2)
          viru = submesh%C( virc, submesh%nC( virc))
          virl = submesh%C( virc, 1)

          pc   = submesh%V(virc,:)
          pu   = submesh%V(viru,:)
          pl   = submesh%V(virl,:)
          
          IF (orientation == 0) THEN

            IF ( pu(2) > p_new(2) .AND. pc(2) < p_new(2)) THEN
              nT = nT+1
              T(nT,:) = [vil, vir]
              CALL split_segment( submesh, virc, viru, p_new)
              EXIT
            ELSEIF ( pc(2) > p_new(2) .AND. pl(2) < p_new(2)) THEN
              nT = nT+1
              T(nT,:) = [vil, vir]
              CALL split_segment( submesh, virl, virc, p_new)
              EXIT
            END IF
          
          ELSEIF (orientation == 1) THEN

            IF ( pu(1) < p_new(1) .AND. pc(1) > p_new(1)) THEN
              nT = nT+1
              T(nT,:) = [vil, vir]
              CALL split_segment( submesh, virc, viru, p_new)
              EXIT
            ELSEIF ( pc(1) < p_new(1) .AND. pl(1) > p_new(1)) THEN
              nT = nT+1
              T(nT,:) = [vil, vir]
              CALL split_segment( submesh, virl, virc, p_new)
              EXIT
            END IF
          
          END IF ! IF (orientation == 0) THEN
          
        END DO ! DO vi2 = 1, nVr_west        
      END DO ! DO vi = 1, nVl_east
      
    END IF ! IF (par%i == p_left) THEN
    
    ! Crop submesh memory
    CALL crop_submesh_primary( submesh) 
    
    ! Clean up after yourself!
    DEALLOCATE( Vil_east     )
    DEALLOCATE( Vir_west     )   
    DEALLOCATE( Vl_east      )
    DEALLOCATE( Vr_west      )
    DEALLOCATE( T            )
    DEALLOCATE( FoundAsMatchl)
    DEALLOCATE( FoundAsMatchr)
    
!    WRITE(0,'(A,I2,A)') '  align_submeshes - process ', par%i, ': finished'
    
  END SUBROUTINE align_submeshes
  SUBROUTINE merge_all_submeshes( submesh, orientation)
    ! Iteratively merge the submeshes created by the different processes.
  
    TYPE(type_mesh),            INTENT(INOUT)     :: submesh
    INTEGER,                    INTENT(IN)        :: orientation
    
    INTEGER, DIMENSION(:,:), ALLOCATABLE          :: mergelist
    INTEGER                                       :: nmerge, merge_it, n_merge_it, i
    
    ! No need for this if we're running on a single core
    IF (par%n == 1) RETURN
    
    ALLOCATE( mergelist( par%n, 2))
    
    ! Determine number of required merging iterations
    n_merge_it = 1
    DO WHILE (2**n_merge_it < par%n)
      n_merge_it = n_merge_it + 1
    END DO
    
    ! Iteratively merge meshes
    DO merge_it = 1, n_merge_it
      
      IF (debug_mesh_creation .AND. par%master) WRITE(0,*) '  Merging submeshes: iteration ', merge_it
    
      mergelist = 0
      nmerge    = 0
      
      DO i = 0, par%n-2, 2**merge_it
        IF (i + (2**(merge_it-1)) < par%n) THEN
          nmerge = nmerge + 1
          mergelist(nmerge,:) = [i, i + (2**(merge_it-1))]
        END IF
      END DO
      
      CALL merge_submeshes( submesh, mergelist, nmerge, orientation)
    
    END DO
    
  ! Clean up after yourself
  ! =======================
    
    DEALLOCATE( mergelist)
  
  END SUBROUTINE merge_all_submeshes
  SUBROUTINE merge_submeshes( submesh, mergelist, nmerge, orientation)
    ! Merge the submeshes created by proc1 and proc2
    
    IMPLICIT NONE
  
    ! Input variables
    TYPE(type_mesh),            INTENT(INOUT)     :: submesh
    INTEGER, DIMENSION(:,:),    INTENT(IN)        :: mergelist
    INTEGER,                    INTENT(IN)        :: nmerge
    INTEGER,                    INTENT(IN)        :: orientation  ! 0 = eastwest, 1 = northsouth
    
    ! Local variables
    TYPE(type_mesh)                               :: submesh_right
    
    INTEGER                                       :: status(MPI_STATUS_SIZE)
    
    INTEGER                                       :: p_left, p_right, i
    
    INTEGER,  DIMENSION(:,:), ALLOCATABLE         :: T
    INTEGER                                       :: nT, ti, nVl_east, nVr_west
    LOGICAL                                       :: FoundNorth
    INTEGER                                       :: vi, vil, vir, vil_prev, vir_prev
    
    INTEGER                                       :: nVl, nVr, nVtot, nTril, nTrir, nTritot, RefStackNl, RefStackNr, RefStackNtot
    INTEGER                                       :: ci, iti, n, vi1, vi2, nf, ti1, ti2
    LOGICAL                                       :: did_flip
    
    INTEGER                                       :: cip1
    REAL(dp)                                      :: VorTriA, sumVorTriA
    REAL(dp), DIMENSION(2)                        :: pa, pb, pc, VorGC
        
    ! Determine if we're participating as Left, Right or Passive
    ! (Passive still needs to run through the code, to participate in ExtendSubmesh calls)
    
    p_left  = -1
    p_right = -1
    DO i = 1, nmerge
      IF ( mergelist(i,1) == par%i .OR. mergelist(i,2) == par%i) THEN
        p_left  = mergelist(i,1)
        p_right = mergelist(i,2)
      END IF
    END DO
    
    IF (par%i == p_left .OR. par%i == p_right) THEN
      IF (debug_mesh_creation) WRITE(0,'(A,I2,A,I2,A,I2)') ' merge_submeshes - process ', par%i, ': merging mesh ', p_left, ' with mesh ', p_right
    ELSE
      IF (debug_mesh_creation) WRITE(0,'(A,I2,A)')         ' merge_submeshes - process ', par%i, ': passively running merge_submeshes'
    END IF
    CALL sync
   
  ! Align the two submeshes (make sure all their shared boundary vertices match)
  ! ============================================================================
     
    CALL align_submeshes( submesh, orientation, p_left, p_right, nVl_east, nVr_west)
    
  ! Extend memory to accomodate data from submesh_right
  ! ===================================================

    ! First communicate size of aligned meshes, so p_left can extend their memory to accomodate data from p_right
    ! (must be done before sharing memory access, since ExtendSubmesh realloates the memory, resulting in new addresses)

    IF (par%i == p_left) THEN
      CALL MPI_RECV( nVr,   1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      CALL MPI_RECV( nTrir, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      nVl     = submesh%nV
      nTril   = submesh%nTri
      nVtot   = nVl + nVr
      nTritot = nTril + nTrir
    ELSEIF (par%i == p_right) THEN
      CALL MPI_SEND( submesh%nV,   1, MPI_INTEGER, p_left, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_SEND( submesh%nTri, 1, MPI_INTEGER, p_left, 0, MPI_COMM_WORLD, ierr)
    END IF

    IF (par%i == p_left) THEN
      CALL extend_submesh_primary( submesh, nVtot, nTritot)
    ELSE
      CALL extend_submesh_primary( submesh, submesh%nV, submesh%nTri)
    END IF
    
  ! Give p_left access to submesh memory from p_right
  ! =================================================  
    
    IF (par%i == p_left .OR. par%i == p_right) THEN
      IF (debug_mesh_creation .AND. par%i == p_left)  WRITE(0,'(A,I2,A,I2)') ' merge_submeshes - process ', par%i, ': gaining   access to submesh memory from process ', p_right
      IF (debug_mesh_creation .AND. par%i == p_right) WRITE(0,'(A,I2,A,I2)') ' merge_submeshes - process ', par%i, ': providing access to submesh memory to   process ', p_left
      CALL share_submesh_access( p_left, p_right, submesh, submesh_right)
    END IF
    
  ! Merge the data
  ! ==============

    IF (par%i == p_left) THEN
    
  ! Recalculate T with the extra vertices, again sorthed south-north
  ! (must be done to know which vertices to merge. Can be done much
  ! more easily now that we have access to the data from submesh_right)
  ! ===================================================================
  
      IF (debug_mesh_creation) WRITE(0,'(A,I2,A,I2)') ' merge_submeshes - process ', par%i, ': recalculating T'
      
      ALLOCATE( T( nVl_east + nVr_west, 2))
      
      IF (orientation == 0) THEN 
        
        T = 0
        T(1, :) = [2, 1]
        nT = 1
        
        vil_prev = 2
        vir_prev = 1
        FoundNorth = .FALSE.
        
        DO WHIlE (.NOT. FoundNorth)
          vil = submesh%C(       vil_prev, 1)
          vir = submesh_right%C( vir_prev, submesh_right%nC( vir_prev))
          nT = nT+1
          T(nT,:) = [vil, vir]
          vil_prev = vil
          vir_prev = vir
          IF (vil == 3) FoundNorth = .TRUE.
        END DO
      
      ELSEIF (orientation == 1) THEN
        
        T = 0
        T(1, :) = [3, 2]
        nT = 1
        
        vil_prev = 3
        vir_prev = 2
        FoundNorth = .FALSE.
        
        DO WHIlE (.NOT. FoundNorth)
          vil = submesh%C(       vil_prev, 1)
          vir = submesh_right%C( vir_prev, submesh_right%nC( vir_prev))
          nT = nT+1
          T(nT,:) = [vil, vir]
          vil_prev = vil
          vir_prev = vir
          IF (vil == 4) FoundNorth = .TRUE.
        END DO
      
      END IF ! IF (orientation == 0) THEN
      
    ! Add the data from submesh_right to this process' submesh
    ! ======================================================== 
  
      IF (debug_mesh_creation) WRITE(0,'(A,I2,A,I2)') ' merge_submeshes - process ', par%i, ': adding data from p_right to p_left'     
  
      nVl          = submesh%nV
      nVr          = submesh_right%nV
      nVtot        = nVl + nVr
      
      nTril        = submesh%nTri
      nTrir        = submesh_right%nTri
      nTritot      = nTril + nTrir
      
      RefStackNl   = submesh%RefStackN
      RefStackNr   = submesh_right%RefStackN
      RefStackNtot = RefStackNl + RefStackNr
      
      submesh%V(              nVl   +1:nVtot  ,:) = submesh_right%V(              1:nVr,  :)
      submesh%nC(             nVl   +1:nVtot    ) = submesh_right%nC(             1:nVr    )
      submesh%C(              nVl   +1:nVtot  ,:) = submesh_right%C(              1:nVr,  :)
      submesh%niTri(          nVl   +1:nVtot    ) = submesh_right%niTri(          1:nVr    )
      submesh%iTri(           nVl   +1:nVtot  ,:) = submesh_right%iTri(           1:nVr,  :)
      submesh%edge_index(     nVl   +1:nVtot    ) = submesh_right%edge_index(     1:nVr    )

      submesh%Tri(            nTril +1:nTritot,:) = submesh_right%Tri(            1:nTrir,:)
      submesh%Tricc(          nTril +1:nTritot,:) = submesh_right%Tricc(          1:nTrir,:)
      submesh%TriC(           nTril +1:nTritot,:) = submesh_right%TriC(           1:nTrir,:)
      submesh%Tri_edge_index( nTril +1:nTritot  ) = submesh_right%Tri_edge_index( 1:nTrir  )

      submesh%RefMap(         nTril +1:nTritot  ) = submesh_right%RefMap(         1:nTrir  )
      
      submesh%RefStack(RefStackNl+1:RefStackNl+RefStackNr) = submesh_right%RefStack(1:RefStackNr)
      submesh%RefStackN = RefStackNl + RefStackNr

      submesh%nV   = nVtot
      submesh%nTri = nTritot
            
    ! Increase the vertex and triangle count of data from submesh_right to start where submesh_left stops
    ! ===================================================================================================
  
      DO vi = nVl+1, nVtot
        DO ci = 1, submesh%nC(vi)
          submesh%C(vi,ci) = submesh%C(vi,ci) + nVl
        END DO
        DO iti = 1, submesh%niTri(vi)
          submesh%iTri(vi,iti) = submesh%iTri(vi,iti) + nTril
        END DO
      END DO
      DO ti = nTril+1, nTritot
        submesh%Tri(ti,:) = submesh%Tri(ti,:) + nVl
        DO n = 1, 3
          IF (submesh%TriC(ti,n)>0) submesh%TriC(ti,n) = submesh%TriC(ti,n) + nTril
        END DO
      END DO

      ! Do the same for the boundary translation table
      DO ti = 1, nT
        T(ti,2) = T(ti,2) + nVl
      END DO

      ! And for the Refinement Stack
      DO n = RefStackNl+1, RefStackNl+RefStackNr
        submesh%RefStack(n) = submesh%RefStack(n) + nTril
      END DO

    ! Merge the shared vertices
    ! =========================
  
      IF (debug_mesh_creation) WRITE(0,'(A,I2,A,I2)') ' merge_submeshes - process ', par%i, ': merging shared vertices'  
  
      DO ti = 1, nT
        vil = T(ti,1)
        vir = T(ti,2)
        CALL merge_vertices( submesh, nVl, nVr, nTril, nTrir, T, nT, vil, vir, orientation)
      END DO
      
    ! Update domain boundaries
    ! ========================
      
      IF (orientation == 0) THEN
        submesh%xmin = submesh%xmin
        submesh%xmax = submesh_right%xmax        
      ELSEIF (orientation == 1) THEN
        submesh%ymin = submesh%ymin
        submesh%ymax = submesh_right%ymax   
      END IF

    ! Update Tri_edge_index
    ! =====================
  
      IF (debug_mesh_creation) WRITE(0,'(A,I2,A,I2)') ' merge_submeshes - process ', par%i, ': redoing triangle edge indices'  
    
      CALL redo_Tri_edge_indices( submesh)

    ! Make sure vertices 1, 2, 3, 4 are again at the corners
    ! ======================================================
  
      IF (debug_mesh_creation) WRITE(0,'(A,I2,A,I2)') ' merge_submeshes - process ', par%i, ': resetting corner vertices'  
    
      IF (orientation == 0) THEN
      
        vi1 = 2
        vi2 = nVl+1
        CALL switch_vertices( submesh, vi1, vi2)
        vi1 = 3
        vi2 = nVl+2
        CALL switch_vertices( submesh, vi1, vi2)
      
      ELSEIF (orientation == 1) THEN
      
        vi1 = 3
        vi2 = nVl+1
        CALL switch_vertices( submesh, vi1, vi2)
        vi1 = 4
        vi2 = nVl+2
        CALL switch_vertices( submesh, vi1, vi2)
        
      END IF
      
    ! Check if any seam triangles require flipping
    ! ============================================
  
      IF (debug_mesh_creation) WRITE(0,'(A,I2,A,I2)') ' merge_submeshes - process ', par%i, ': updating Delaunay triangulation'  
    
      DO ti = 1, nT
      
        vi = T( ti,1)
        
        nf = 0
  
        DO iti = 1, submesh%niTri( vi)
          ti1 = submesh%iTri( vi,iti)
          DO n = 1, 3
            ti2 = submesh%TriC( ti1,n)
            IF (ti2 > 0) THEN              
              nf = nf + 1          
              submesh%Triflip(nf,:) = [ti1,ti2]
            END IF
          END DO
        END DO
  
        ! Flip triangle pairs
        DO WHILE (nf > 0)
          CALL flip_triangle_pairs( submesh, nf, did_flip)
        END DO
        
      END DO ! DO ti = 1, nT
      
    ! Lloyd's algorithm: move seam vertices to their Voronoi cell geometric centres
    ! =============================================================================
  
      IF (debug_mesh_creation) WRITE(0,'(A,I2,A,I2)') ' merge_submeshes - process ', par%i, ': apply Lloyds algorithm to seam'
      
      DO ti = 1, nT
        
        vi = T( ti,1)
        
        ! Skip the two boundary vertices
        IF (submesh%edge_index( vi) > 0) CYCLE
        
        ! Find the geometric centre of this vertex' Voronoi cell
        VorGC      = 0._dp
        sumVorTriA = 0._dp
        
        DO ci = 1, submesh%nC( vi)
        
          cip1 = ci + 1
          IF (cip1 > submesh%nC( vi)) cip1 = 1
          
          pa = submesh%V( vi,:)
          pb = submesh%V( submesh%C( vi,ci  ),:)
          pc = submesh%V( submesh%C( vi,cip1),:)
          
          VorTriA = cross2( pb - pa, pc - pa)
          
          VorGC = VorGC + VorTriA * (pa + pb + pc) / 3._dp
          sumVorTriA   = sumVorTriA   + VorTriA
          
        END DO ! DO ci = 1, submesh%nC( vi)
        
        VorGC = VorGC / sumVorTriA
        
        ! Move the vertex
        CALL move_vertex( submesh, vi, VorGC)
        
      END DO ! DO ti = 1, nT
      
    ! Clean up after yourself
    ! =======================
    
      DEALLOCATE(T)
      
    END IF ! IF (par%i == p_left) THEN
    
    CALL sync
    
  END SUBROUTINE merge_submeshes
  
  ! == Once merging is finished, finalise the mesh
  SUBROUTINE create_final_mesh_from_merged_submesh( submesh, mesh)
    
    IMPLICIT NONE
  
    ! In/output variables
    TYPE(type_mesh),            INTENT(INOUT)     :: submesh
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables
    INTEGER                                       :: nV, nTri
    
    ! Communicate final merged mesh size to all processes
    ! ===================================================
    
    nV = submesh%nV
    nTri = submesh%nTri
    CALL MPI_BCAST( nV,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( nTri, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  
    ! Copy data from the final merged submesh, deallocate all the submeshes,
    ! do one final refining pass for the new triangles that may be too sharp.
    ! =======================================================================

    CALL allocate_mesh_primary( mesh, submesh%region_name, nV + 1000, nTri + 2000, submesh%nC_mem)
    IF (par%master) CALL move_data_from_submesh_to_mesh( mesh, submesh)  
    CALL deallocate_submesh_primary( submesh)
    CALL refine_mesh_geo_only( mesh)
    CALL crop_mesh_primary(    mesh)

    ! Finish up - mesh metadata and extra info
    ! ========================================
                
    ! Determine vertex and triangle domains
    CALL partition_list( mesh%nV,   par%i, par%n, mesh%vi1, mesh%vi2)
    CALL partition_list( mesh%nTri, par%i, par%n, mesh%ti1, mesh%ti2)
    
    ! Calculate extra mesh data
    CALL allocate_mesh_secondary(             mesh)
    CALL calc_triangle_geometric_centres(     mesh)
    CALL find_Voronoi_cell_areas(             mesh)
    CALL get_lat_lon_coordinates(             mesh)
    CALL find_triangle_areas(                 mesh)
    CALL find_connection_widths(              mesh)
    CALL make_Ac_mesh(                        mesh)
    CALL calc_matrix_operators_mesh(          mesh)
    CALL determine_mesh_resolution(           mesh)
    IF (par%master) CALL find_POI_xy_coordinates( mesh)
    CALL sync
    CALL find_POI_vertices_and_weights(       mesh)
    CALL find_Voronoi_cell_geometric_centres( mesh)
    CALL create_transect(                     mesh)
    CALL calculate_five_colouring_AaAc(       mesh)
    
    CALL check_mesh( mesh)
    
  END SUBROUTINE create_final_mesh_from_merged_submesh
  
  ! == Create the list of vertex indices and weights used in making transects
  SUBROUTINE create_transect( mesh)
    ! Create a transect along the x-axis, through the centre of the model domain.
    ! Useful for benchmark experiments, not so much for realistic simulations.
    
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: p_start, p_end
    INTEGER                                       :: vi_cur, vj_cur, ti_cur, vi_next, vj_next, ti_next, vi_end, vj_end, ti_end, ti
    INTEGER                                       :: nV_transect
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE       :: vi_transect
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: w_transect
    INTEGER                                       :: n, n1
    REAL(dp), DIMENSION(2)                        :: pi_next, pj_next, llis
    LOGICAL                                       :: do_cross
    INTEGER                                       :: iti, n2, n3
    
    ! Allocate temporary memory for the list of transect vertex pairs
    ! (since we don't know in advance how many there will be)
    nV_transect = CEILING( 2._dp * (mesh%xmax - mesh%xmin) / mesh%resolution_min)
    ALLOCATE( vi_transect( nV_transect,2))
    ALLOCATE( w_transect(  nV_transect,2))
    nV_transect = 0
    vi_transect = 0
    
    ! Done only by the Master
    IF (par%master) THEN
      
      ! The start and end points of the transect
      p_start = [mesh%xmin, (mesh%ymin + mesh%ymax) / 2._dp]
      p_end   = [mesh%xmax, (mesh%ymin + mesh%ymax) / 2._dp]
      
      ! Find the pair of vertices on whose connection p_start lies, and the
      ! triangle whose side is made up by that connection
      vi_cur = 0
      vj_cur = 1
      DO WHILE (mesh%V( vj_cur,2) < p_start(2))
        vi_cur = vj_cur
        vj_cur = mesh%C( vi_cur, mesh%nC( vi_cur))
      END DO
      ti_cur = mesh%iTri( vi_cur, mesh%niTri( vi_cur))
      
      ! Exception for Greenland: sometimes the process domain border can intersect with the transect,
      ! so that it overlaps with lines everywhere. In this case, move it slightly.
      IF     (p_start(2) == mesh%V( vi_cur,2)) THEN
        p_start(2) = p_start(2) + 100._dp
        p_end(  2) = p_end(  2) + 100._dp
      ELSEIF (p_start(2) == mesh%V( vj_cur,2)) THEN
        p_start(2) = p_start(2) - 100._dp
        p_end(  2) = p_end(  2) - 100._dp
      END IF
      
      ! List this as the first pair of transect vertices
      nV_transect = 1
      vi_transect( 1,:) = [vi_cur,vj_cur]
      w_transect(1,:) = [NORM2( mesh%V( vj_cur,:) - p_start), NORM2( mesh%V( vi_cur,:) - p_start)] / NORM2( mesh%V( vj_cur,:) - mesh%V( vi_cur,:))
      
      ! Find the pair of vertices on whose connection p_end lies, and the
      ! triangle whose side is made up by that connection
      vi_end = 0
      vj_end = 2
      DO WHILE (mesh%V( vj_end,2) < p_end(2))
        vi_end = vj_end
        vj_end = mesh%C( vi_end, 1)
      END DO
      ti_end = mesh%iTri( vi_end, 1)
      
      ! Trace the transect through the mesh
      DO WHILE (ti_cur /= ti_end)
        
        ! Find out where the transect exits the current triangle
        ti_next = 0
        DO n = 1, 3
          
          n1 = n + 1
          IF (n1==4) n1 = 1
          
          vi_next = mesh%Tri( ti_cur,n)
          vj_next = mesh%Tri( ti_cur,n1)
          
          IF ((vi_next == vi_cur .AND. vj_next == vj_cur) .OR. (vi_next == vj_cur .AND. vj_next == vi_cur)) CYCLE
          
          pi_next = mesh%V( vi_next,:)
          pj_next = mesh%V( vj_next,:)
          
          CALL segment_intersection( p_start, p_end, pi_next, pj_next, llis, do_cross, mesh%tol_dist)
          
          IF (do_cross) THEN
            ! Find the next triangle
            DO iti = 1, mesh%niTri( vi_next)
              ti = mesh%iTri( vi_next, iti)
              DO n2 = 1, 3
                n3 = n2 + 1
                IF (n3 == 4) n3 = 1
                IF (((mesh%Tri( ti,n2) == vj_next .AND. mesh%Tri( ti,n3) == vi_next) .OR. &
                     (mesh%Tri( ti,n2) == vi_next .AND. mesh%Tri( ti,n3) == vj_next)) .AND. &
                    ti /= ti_cur) THEN
                  ti_next = ti
                  EXIT
                END IF
              END DO
            END DO ! DO iti = 1, mesh%niTri( vi_next)
          
            EXIT
          END IF ! IF (do_cross) THEN
          
        END DO ! DO n = 1, 3
        
        ! Check if we managed to find the crossing
        IF (ti_next == 0) THEN
          WRITE(0,*) '  create_transect - ERROR: couldnt find next triangle along transect!'
          STOP
        END IF
        
        ! Add this new vertex pair to the list
        nV_transect = nV_transect + 1
        vi_transect( nV_transect,:) = [vi_next, vj_next]
        w_transect(  nV_transect,:) = [NORM2( mesh%V( vj_next,:) - llis), NORM2( mesh%V( vi_next,:) - llis)] / NORM2( mesh%V( vj_next,:) - mesh%V( vi_next,:))
        
        ! Cycle the vertex pairs
        vi_cur = vi_next
        vj_cur = vj_next
        ti_cur = ti_next
        
      END DO ! DO WHILE (ti_cur /= ti_end)
      
      ! Add the final pair of vertices to the transect
      nV_transect = nV_transect + 1
      vi_transect( nV_transect,:) = [vi_end, vj_end]
      w_transect(  nV_transect,:) = [NORM2( mesh%V( vj_end,:) - p_end), NORM2( mesh%V( vi_end,:) - p_end)] / NORM2( mesh%V( vj_end,:) - mesh%V( vi_end,:))
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate shared memory, copy list of transect vertex pairs
    CALl allocate_shared_int_0D( mesh%nV_transect, mesh%wnV_transect)
    IF (par%master) mesh%nV_transect = nV_transect
    CALL sync
    CALL allocate_shared_int_2D( mesh%nV_transect, 2, mesh%vi_transect, mesh%wvi_transect)
    CALL allocate_shared_dp_2D(  mesh%nV_transect, 2, mesh%w_transect,  mesh%ww_transect )
    IF (par%master) mesh%vi_transect = vi_transect( 1:mesh%nV_transect,:)
    IF (par%master) mesh%w_transect  = w_transect(  1:mesh%nV_transect,:)
    CALL sync
    
    ! Clean up after yourself
    DEALLOCATE( vi_transect)
    DEALLOCATE( w_transect)
    
  END SUBROUTINE create_transect
  
  ! == Initialise a five-vertex dummy mesh
  SUBROUTINE initialise_dummy_mesh( mesh, xmin, xmax, ymin, ymax)
    ! Initialises a 5-vertex, 4-triangle "dummy"  mesh:
    !
    !   v4 - - - - - - - - v3   V          nC     C             niTri   iTri          edge_index
    !   | \              / |    -1 -1      3      2  5  4        2      1  4            6
    !   |  \            /  |     1 -1      3      3  5  1        2      2  1            4
    !   |   \    t3    /   |     1  1      3      4  5  2        2      3  2            2
    !   |    \        /    |    -1  1      3      1  5  3        2      4  3            8
    !   |     \      /     |     0  0      4      1  2  3  4     4      1  2  3  4      0
    !   |      \    /      |
    !   |       \  /       |    Tri           TriC
    !   |  t4    v5    t2  |    1  2  5      2  4  0
    !   |       /  \       |    2  3  5      3  1  0
    !   |      /    \      |    3  4  5      4  2  0
    !   |     /      \     |    4  1  5      1  3  0
    !   |    /        \    |
    !   |   /    t1    \   |
    !   |  /            \  |
    !   | /              \ |
    !   v1 - - - - - - - - v2
    
    IMPLICIT NONE

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    REAL(dp),                   INTENT(IN)        :: xmin    ! X and Y ranges of the newly initialised (dummy) mesh
    REAL(dp),                   INTENT(IN)        :: xmax
    REAL(dp),                   INTENT(IN)        :: ymin
    REAL(dp),                   INTENT(IN)        :: ymax

    ! Local variables
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp

    ! Meta properties
    mesh%xmin                 = xmin    ! Boundaries of the square domain.
    mesh%xmax                 = xmax
    mesh%ymin                 = ymin
    mesh%ymax                 = ymax
    mesh%alpha_min            = C%alpha_min
    mesh%dz_max_ice           = C%dz_max_ice
    mesh%res_max              = C%res_max
    mesh%res_max_margin       = C%res_max_margin
    mesh%res_max_gl           = C%res_max_gl
    mesh%res_max_cf           = C%res_max_cf
    mesh%res_max_mountain     = C%res_max_mountain
    mesh%res_max_coast        = C%res_max_coast
    mesh%res_min              = C%res_min
    
    ! Horizontal distance of tolerance. Used for several small routines - points that lie within
    ! this distance of each other (vertices, triangle circumcenters, etc.) are treated as identical.
    mesh%tol_dist   = ((mesh%xmax - mesh%xmin) + (mesh%ymax-mesh%ymin)) * tol / 2._dp
    
    ! Points of interest
    CALL find_POI_xy_coordinates( mesh)

    ! The four corners, plus one central vertex.
    mesh%nV           = 5

    mesh%V            = 0._dp
    mesh%V(1,:)       = [xmin, ymin]
    mesh%V(2,:)       = [xmax, ymin]
    mesh%V(3,:)       = [xmax, ymax]
    mesh%V(4,:)       = [xmin, ymax]
    mesh%V(5,:)       = [(xmin+xmax)/2, (ymin+ymax)/2]

    mesh%edge_index      = 0
    mesh%edge_index(1:5) = [6, 4, 2, 8, 0]

    mesh%nC          = 0
    mesh%nC(1:5)     = [3, 3, 3, 3, 4]

    mesh%C           = 0
    mesh%C(1,1:4)    = [2, 5, 4, 0]
    mesh%C(2,1:4)    = [3, 5, 1, 0]
    mesh%C(3,1:4)    = [4, 5, 2, 0]
    mesh%C(4,1:4)    = [1, 5, 3, 0]
    mesh%C(5,1:4)    = [1, 2, 3, 4]

    mesh%niTri       = 0
    mesh%niTri(1:5)  = [2, 2, 2, 2, 4]

    mesh%iTri        = 0
    mesh%iTri(1,1:4) = [1, 4, 0, 0]
    mesh%iTri(2,1:4) = [2, 1, 0, 0]
    mesh%iTri(3,1:4) = [3, 2, 0, 0]
    mesh%iTri(4,1:4) = [4, 3, 0, 0]
    mesh%iTri(5,1:4) = [1, 2, 3, 4]

    mesh%nTri         = 4

    mesh%Tri          = 0
    mesh%Tri(1,:)     = [1, 2, 5]
    mesh%Tri(2,:)     = [2, 3, 5]
    mesh%Tri(3,:)     = [3, 4, 5]
    mesh%Tri(4,:)     = [4, 1, 5]

    mesh%Tri_edge_index      = 0
    mesh%Tri_edge_index(1:4) = [5, 3, 1, 7]

    mesh%TriC        = 0
    mesh%TriC(1,:)   = [2, 4, 0]
    mesh%TriC(2,:)   = [3, 1, 0]
    mesh%TriC(3,:)   = [4, 2, 0]
    mesh%TriC(4,:)   = [1, 3, 0]

    mesh%Triflip     = 0

    mesh%TriCC = 0._dp
    CALL update_triangle_circumcenter( mesh, 1)
    CALL update_triangle_circumcenter( mesh, 2)
    CALL update_triangle_circumcenter( mesh, 3)
    CALL update_triangle_circumcenter( mesh, 4)
  
    ! Map and stack used for refining
    mesh%RefMap    = 0
    mesh%RefStack  = 0
    mesh%RefStackN = 0

  END SUBROUTINE initialise_dummy_mesh
  SUBROUTINE perturb_dummy_mesh( mesh, perturb_dir_global)
    ! "Perturb" the five-vertex dummy mesh; slightly offset the centre vertex,
    ! and split the four edges just next to their midpoints. This ensures that
    ! any new triangles created during mesh refinement are never cocircular.
    
    USE parameters_module, ONLY: pi
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: perturb_dir_global
    
    ! Local variables:
    INTEGER                                       :: perturb_dir_local, i, vi1, vi2
    REAL(dp), DIMENSION(2)                        :: p
    REAL(dp)                                      :: dx, dy
    
    perturb_dir_local = perturb_dir_global
    DO i = 1, par%i
      perturb_dir_local = 1 - perturb_dir_local
    END DO
    
    IF (perturb_dir_local == 0) THEN
      dx = (mesh%xmax - mesh%xmin) *  pi     / 1000._dp ! Offset in x-direction by ~0.3 % of domain width
      dy = (mesh%ymax - mesh%ymin) * (pi**2) / 1000._dp ! Offset in y-direction by ~1   % of domain width
    ELSEIF (perturb_dir_local == 1) THEN
      dx = (mesh%xmin - mesh%xmax) *  pi     / 1000._dp ! Offset in x-direction by ~0.3 % of domain width
      dy = (mesh%ymin - mesh%ymax) * (pi**2) / 1000._dp ! Offset in y-direction by ~1   % of domain width
    ELSE
      WRITE(0,*) ' ERROR: mesh perturbation direction can only be 0 or 1!'
      STOP
    END IF
    
    ! Save perturbation direction
    mesh%perturb_dir = perturb_dir_local
    
    ! Offset the center vertex
    mesh%V( 5,1) = mesh%V( 5,1) + dx / 2._dp
    mesh%V( 5,2) = mesh%V( 5,2) + dy / 2._dp
    
    ! Update triangle circumcenters
    CALL update_triangle_circumcenter( mesh, 1)
    CALL update_triangle_circumcenter( mesh, 2)
    CALL update_triangle_circumcenter( mesh, 3)
    CALL update_triangle_circumcenter( mesh, 4)
    
    ! Split the southern edge
    vi1 = 1
    vi2 = 2
    p = [(mesh%xmax + mesh%xmin) / 2._dp + dx, mesh%ymin]
    CALL split_segment( mesh, vi1, vi2, p)
    
    ! Split the eastern edge
    vi1 = 2
    vi2 = 3
    p = [mesh%xmax, (mesh%ymax + mesh%ymin) / 2._dp + dy]
    CALL split_segment( mesh, vi1, vi2, p)
    
    ! Split the northern edge
    vi1 = 3
    vi2 = 4
    p = [(mesh%xmax + mesh%xmin) / 2._dp - dx, mesh%ymax]
    CALL split_segment( mesh, vi1, vi2, p)
    
    ! Split the western edge
    vi1 = 4
    vi2 = 1
    p = [mesh%xmin, (mesh%ymax + mesh%ymin) / 2._dp - dy]
    CALL split_segment( mesh, vi1, vi2, p)
    
  END SUBROUTINE perturb_dummy_mesh

END MODULE mesh_creation_module
