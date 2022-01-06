MODULE general_ice_model_data_module

  ! Only "secondary geometry" right now: masks, surface elevation, TAF

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
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  
  ! Import specific functionality
  USE data_types_module,               ONLY: type_model_region, type_mesh, type_ice_model
  USE utilities_module,                ONLY: is_floating, surface_elevation, thickness_above_floatation
  USE mesh_help_functions_module,      ONLY: find_triangle_area
  USE mesh_operators_module,           ONLY: map_a_to_ac_2D, map_ac_to_bb_2D, map_a_to_b_2D, map_a_to_c_2D, map_c_to_a_2D, map_bb_to_ac_2D

  IMPLICIT NONE

CONTAINS
  
  ! Routines for calculating general ice model data - Hs, masks, ice physical properties
  SUBROUTINE update_general_ice_model_data( mesh, ice)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables
    CHARACTER(LEN=64), PARAMETER                       :: routine_name = 'update_general_ice_model_data'
    INTEGER                                            :: n1, n2
    INTEGER                                            :: vi
    
    n1 = par%mem%n
    
    ! Calculate surface elevation and thickness above floatation
    DO vi = mesh%vi1, mesh%vi2
      ice%Hs_a(  vi) = surface_elevation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))
      ice%TAF_a( vi) = thickness_above_floatation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))
    END DO
    CALL sync
    
    ! Determine masks
    CALL determine_masks( mesh, ice)
    
    n2 = par%mem%n
    !CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE update_general_ice_model_data
  SUBROUTINE determine_masks( mesh, ice)
    ! Determine the different masks, on both the Aa and the Ac mesh
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
  
    INTEGER                                       :: vi, ci, vc
    
    ! Start out with land everywhere, fill in the rest based on input.
    ice%mask_land_a(   mesh%vi1:mesh%vi2) = 1
    ice%mask_ocean_a(  mesh%vi1:mesh%vi2) = 0
    ice%mask_lake_a(   mesh%vi1:mesh%vi2) = 0
    ice%mask_ice_a(    mesh%vi1:mesh%vi2) = 0
    ice%mask_sheet_a(  mesh%vi1:mesh%vi2) = 0
    ice%mask_shelf_a(  mesh%vi1:mesh%vi2) = 0
    ice%mask_coast_a(  mesh%vi1:mesh%vi2) = 0
    ice%mask_margin_a( mesh%vi1:mesh%vi2) = 0
    ice%mask_gl_a(     mesh%vi1:mesh%vi2) = 0
    ice%mask_cf_a(     mesh%vi1:mesh%vi2) = 0
    ice%mask_a(        mesh%vi1:mesh%vi2) = C%type_land
    CALL sync
  
    DO vi = mesh%vi1, mesh%vi2
      
      ! Determine ocean (both open and shelf-covered)
      IF (is_floating( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))) THEN
        ice%mask_ocean_a( vi) = 1
        ice%mask_land_a(  vi) = 0
        ice%mask_a(       vi) = C%type_ocean
      END IF
    
      ! Determine ice
      IF (ice%Hi_a( vi) > 0._dp) THEN
        ice%mask_ice_a( vi)  = 1
      END IF
      
      ! Determine sheet
      IF (ice%mask_ice_a( vi) == 1 .AND. ice%mask_land_a( vi) == 1) THEN
        ice%mask_sheet_a( vi) = 1
        ice%mask_a(       vi) = C%type_sheet
      END IF
    
      ! Determine shelf
      IF (ice%mask_ice_a( vi) == 1 .AND. ice%mask_ocean_a( vi) == 1) THEN
        ice%mask_shelf_a( vi) = 1
        ice%mask_a(       vi) = C%type_shelf
      END IF
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
  
    ! Determine coast, grounding line and calving front
    DO vi = mesh%vi1, mesh%vi2
    
      IF (ice%mask_land_a( vi) == 1) THEN  
        ! Land bordering ocean equals coastline
        
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ocean_a( vc) == 1) THEN
            ice%mask_a( vi) = C%type_coast
            ice%mask_coast_a( vi) =  1
          END IF
        END DO
      END IF
    
      IF (ice%mask_ice_a( vi) == 1) THEN  
        ! Ice bordering non-ice equals margin
        
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ice_a( vc) == 0) THEN
            ice%mask_a( vi) = C%type_margin
            ice%mask_margin_a( vi) =  1
          END IF
        END DO
      END IF
    
      IF (ice%mask_sheet_a( vi) == 1) THEN  
        ! Sheet bordering shelf equals groundingline
        
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_shelf_a( vc) == 1) THEN
            ice%mask_a( vi) = C%type_groundingline
            ice%mask_gl_a( vi) = 1
          END IF
        END DO
      END IF
  
      IF (ice%mask_ice_a( vi) == 1) THEN  
        ! Ice (sheet or shelf) bordering open ocean equals calvingfront
        
        DO ci = 1, mesh%nC(vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ocean_a( vc) == 1 .AND. ice%mask_ice_a( vc) == 0) THEN
            ice%mask_a( vi) = C%type_calvingfront
            ice%mask_cf_a( vi) = 1
          END IF
        END DO
        
      END IF
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
  
  END SUBROUTINE determine_masks
  
! == Routines for calculating sub-grid grounded fractions
  SUBROUTINE determine_grounded_fractions_ac( mesh, ice)
    ! Determine the grounded fraction of next-to-grounding-line pixels on the ac-grid
    ! (used for determining basal friction in the DIVA)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER                ::  TAF_b
    INTEGER                                            :: wTAF_b
    INTEGER                                            :: avi, vi, aci, ci, vj, iti, iti2, ti1, ti2, vl, vr, ti, n1, n2, n3, til, tir
    REAL(dp)                                           :: TAF_max, TAF_min
    REAL(dp), DIMENSION(2)                             :: va, ccb1, ccb2
    REAL(dp)                                           :: TAFa, TAFb, TAFc, A_vor, A_tri_tot, A_tri_grnd, A_grnd
    REAL(dp), DIMENSION(:    ), POINTER                ::  f_grnd_bb
    INTEGER                                            :: wf_grnd_bb
    
    ! Map thickness-above-floatation to the b-grid
    CALL allocate_shared_dp_1D( mesh%nTri, TAF_b, wTAF_b)
    CALL map_a_to_b_2D(  mesh, ice%TAF_a, TAF_b)
    
  ! ===== On the vertices =====
  ! ===========================
  
    DO vi = mesh%vi1, mesh%vi2
      
      avi = vi
      
      ! Skip border vertices
      IF (mesh%edge_index( vi) > 0) THEN
        ice%f_grnd_ac( avi) = 0._dp
        CYCLE
      END IF
      
      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = -1E6_dp
      TAF_min =  1E6_dp
      
      TAF_max = MAX( TAF_max, ice%TAF_a( vi))
      TAF_min = MIN( TAF_min, ice%TAF_a( vi))
      
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        TAF_max = MAX( TAF_max, ice%TAF_a( vj))
        TAF_min = MIN( TAF_min, ice%TAF_a( vj))
      END DO
      
      ! If the entire local neighbourhood is grounded, the answer is trivial
      IF (TAF_min >= 0._dp) THEN
        ice%f_grnd_ac( avi) = 1._dp
        CYCLE
      END IF
      
      ! If the entire local neighbourhood is floating, the answer is trivial
      IF (TAF_max <= 0._dp) THEN
        ice%f_grnd_ac( avi) = 0._dp
        CYCLE
      END IF
      
      ! The local neighbourhood contains both grounded and floating vertices.
      A_vor  = 0._dp
      A_grnd = 0._dp
      
      va   = mesh%V( vi,:)
      TAFa = ice%TAF_a( vi)
      
      DO iti = 1, mesh%niTri( vi)
        
        iti2 = iti + 1
        IF (iti == mesh%niTri( vi)) iti2 = 1
        
        ti1 = mesh%iTri( vi,iti )
        ti2 = mesh%iTri( vi,iti2)
        
        ccb1 = mesh%Tricc( ti1,:)
        ccb2 = mesh%Tricc( ti2,:)
        
        TAFb = TAF_b( ti1)
        TAFc = TAF_b( ti2)
        
        ! Determine total area of, and grounded area within, this subtriangle
        CALL determine_grounded_fraction_triangle( va, ccb1, ccb2, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
        
        A_vor  = A_vor  + A_tri_tot
        A_grnd = A_grnd + A_tri_grnd
        
      END DO ! DO iati = 1, mesh%niTriAaAc( avi)
      
      ! Calculate the grounded fraction of this Voronoi cell
      ice%f_grnd_ac( avi) = A_grnd / A_vor
      
    END DO
    CALL sync
    
  ! ===== On the edges =====
  ! ========================
  
    DO aci = mesh%ci1, mesh%ci2
      
      avi = aci + mesh%nV
      
      ! Skip border edges
      IF (mesh%edge_index_ac( aci) > 0) THEN
        ice%f_grnd_ac( avi) = 0._dp
        CYCLE
      END IF
      
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      vl = mesh%Aci( aci,3)
      vr = mesh%Aci( aci,4)
      
      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = MAXVAL( [ice%TAF_a( vi), ice%TAF_a( vj), ice%TAF_a( vl), ice%TAF_a( vr)])
      TAF_min = MINVAL( [ice%TAF_a( vi), ice%TAF_a( vj), ice%TAF_a( vl), ice%TAF_a( vr)])
      
      ! If the entire local neighbourhood is grounded, the answer is trivial
      IF (TAF_min >= 0._dp) THEN
        ice%f_grnd_ac( avi) = 1._dp
        CYCLE
      END IF
      
      ! If the entire local neighbourhood is floating, the answer is trivial
      IF (TAF_max <= 0._dp) THEN
        ice%f_grnd_ac( avi) = 0._dp
        CYCLE
      END IF
      
      ! The local neighbourhood contains both grounded and floating vertices.
      
      ! Find the two triangles adjacent to this edge
      til = 0
      tir = 0
      DO iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        DO n1 = 1, 3
          n2 = n1 + 1
          IF (n2 == 4) n2 = 1
          n3 = n2 + 1
          IF (n3 == 4) n3 = 1
          IF (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj .AND. mesh%Tri( ti,n3) == vl) til = ti
          IF (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vr .AND. mesh%Tri( ti,n3) == vj) tir = ti
        END DO
      END DO
      IF (til == 0 .OR. tir == 0) THEN
        WRITE(0,*) 'determine_grounded_fractions_ac - whaa!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
      ! Find grounded fractions of the two sub-triangles
      A_vor  = 0._dp
      A_grnd = 0._dp
      
      ! Lower subtriangle: [vi,tir,til]
      va   = mesh%V( vi,:)
      ccb1 = mesh%Tricc( tir,:)
      ccb2 = mesh%Tricc( til,:)
      
      TAFa = ice%TAF_a( vi)
      TAFb = TAF_b( tir)
      TAFc = TAF_b( til)
      
      CALL determine_grounded_fraction_triangle( va, ccb1, ccb2, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
      
      A_vor  = A_vor  + A_tri_tot
      A_grnd = A_grnd + A_tri_grnd
      
      ! Upper subtriangle: [vj,til,tir]
      va   = mesh%V( vj,:)
      ccb1 = mesh%Tricc( til,:)
      ccb2 = mesh%Tricc( tir,:)
      
      TAFa = ice%TAF_a( vj)
      TAFb = TAF_b( til)
      TAFc = TAF_b( tir)
      
      CALL determine_grounded_fraction_triangle( va, ccb1, ccb2, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
      
      A_vor  = A_vor  + A_tri_tot
      A_grnd = A_grnd + A_tri_grnd
      
      ! Calculate the grounded fraction of this Voronoi cell
      ice%f_grnd_ac( avi) = A_grnd / A_vor
      
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wTAF_b)
    
    ! DENK DROM
    CALL allocate_shared_dp_1D( mesh%nTriAaAc, f_grnd_bb, wf_grnd_bb)
    CALL map_ac_to_bb_2D( mesh, ice%f_grnd_ac, f_grnd_bb)
    CALL map_bb_to_ac_2D( mesh, f_grnd_bb, ice%f_grnd_ac)
    CALL deallocate_shared( wf_grnd_bb)
    
  END SUBROUTINE determine_grounded_fractions_ac
  SUBROUTINE determine_grounded_fractions_ac_old( mesh, ice)
    ! Determine the grounded fraction of next-to-grounding-line pixels on the ac-grid
    ! (used for determining basal friction in the DIVA)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    REAL(dp), DIMENSION(:    ), POINTER                ::  TAF_ac,  TAF_bb
    INTEGER                                            :: wTAF_ac, wTAF_bb
    INTEGER                                            :: avi, ci, avj, iati, iati2, ati1, ati2
    REAL(dp)                                           :: TAF_max, TAF_min
    REAL(dp), DIMENSION(2)                             :: vac, ccbb1, ccbb2
    REAL(dp)                                           :: TAFa, TAFb, TAFc, A_vor, A_tri_tot, A_grnd, A_tri_grnd
    
    ! Map thickness-above-floatation to the ac- and bb-grids
    CALL allocate_shared_dp_1D( mesh%nVAaAc  , TAF_ac, wTAF_ac)
    CALL allocate_shared_dp_1D( mesh%nTriAaAc, TAF_bb, wTAF_bb)
    CALL map_a_to_ac_2D(  mesh, ice%TAF_a, TAF_ac)
    CALL map_ac_to_bb_2D( mesh, TAF_ac, TAF_bb)
    
    DO avi = mesh%avi1, mesh%avi2
      
      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = -1E6_dp
      TAF_min =  1E6_dp
      
      TAF_max = MAX( TAF_max, TAF_ac( avi))
      TAF_min = MIN( TAF_min, TAF_ac( avi))
      
      DO ci = 1, mesh%nCAaAc( avi)
        avj = mesh%CAaAc( avi,ci)
        TAF_max = MAX( TAF_max, TAF_ac( avj))
        TAF_min = MIN( TAF_min, TAF_ac( avj))
      END DO
      
      ! If the entire local neighbourhood is grounded, the answer is trivial
      IF (TAF_min >= 0._dp) THEN
        ice%f_grnd_ac( avi) = 1._dp
        CYCLE
      END IF
      
      ! If the entire local neighbourhood is floating, the answer is trivial
      IF (TAF_max <= 0._dp) THEN
        ice%f_grnd_ac( avi) = 0._dp
        CYCLE
      END IF
      
      ! The local neighbourhood contains both grounded and floating vertices.
      A_vor  = 0._dp
      A_grnd = 0._dp
      
      vac  = mesh%VAaAc( avi,:)
      TAFa = TAF_ac( avi)
      
      DO iati = 1, mesh%niTriAaAc( avi)
        
        iati2 = iati + 1
        IF (iati == mesh%niTriAaAc( avi)) iati2 = 1
        
        ati1 = mesh%iTriAaAc( avi,iati )
        ati2 = mesh%iTriAaAc( avi,iati2)
        
        ccbb1 = mesh%TriccAaAc( ati1,:)
        ccbb2 = mesh%TriccAaAc( ati2,:)
        
        TAFb = TAF_bb( ati1)
        TAFc = TAF_bb( ati2)
        
        ! Determine total area of, and grounded area within, this subtriangle
        CALL determine_grounded_fraction_triangle( vac, ccbb1, ccbb2, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
        
        A_vor  = A_vor  + A_tri_tot
        A_grnd = A_grnd + A_tri_grnd
        
      END DO ! DO iati = 1, mesh%niTriAaAc( avi)
      
      ! Calculate the grounded fraction of this Voronoi cell
      ice%f_grnd_ac( avi) = A_grnd / A_vor
      
    END DO ! DO avi = mesh%avi1, mesh%avi2
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wTAF_ac)
    CALL deallocate_shared( wTAF_bb)
    
    
    ! DENK DROM
    IF (par%master) THEN
      debug%dp_2D_ac_01 = ice%f_grnd_ac
      CALL write_to_debug_file
    END IF
    CALL sync
    !CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    
  END SUBROUTINE determine_grounded_fractions_ac_old
  SUBROUTINE determine_grounded_fraction_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
    ! Determine the grounded area of the triangle [va,vb,vc], where the thickness-above-floatation is given at all three corners
    
    IMPLICIT NONE
    
    ! In- and output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
    REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
    REAL(dp),                            INTENT(OUT)   :: A_tri_tot, A_tri_grnd
    
    ! Local variables:
    REAL(dp)                                           :: A_flt
    
    ! Determine total area of this subtriangle
    CALL find_triangle_area( va, vb, vc, A_tri_tot)
        
    IF     (TAFa >= 0._dp .AND. TAFb >= 0._dp .AND. TAFc >= 0._dp) THEN
      ! If all three corners are grounded, the answer is trivial
      A_tri_grnd = A_tri_tot
    ELSEIF (TAFa <= 0._dp .AND. TAFb <= 0._dp .AND. TAFc <= 0._dp) THEN
      ! If all three corners are floating, the answer is trivial
      A_tri_grnd = 0._dp
    ELSE
      ! At least one corner is grounded and at least one corner is floating
      
      IF     (TAFa >= 0._dp .AND. TAFb <= 0._dp .AND. TAFc <= 0._dp) THEN
        ! a is grounded, b and c are floating
        CALL determine_grounded_fraction_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
      ELSEIF (TAFa <= 0._dp .AND. TAFb >= 0._dp .AND. TAFc <= 0._dp) THEN
        ! b is grounded, a and c are floating
        CALL determine_grounded_fraction_triangle_1grnd_2flt( vb, vc, va, TAFb, TAFc, TAFa, A_tri_grnd)
      ELSEIF (TAFa <= 0._dp .AND. TAFb <= 0._dp .AND. TAFc >= 0._dp) THEN
        ! c is grounded, a and b are floating
        CALL determine_grounded_fraction_triangle_1grnd_2flt( vc, va, vb, TAFc, TAFa, TAFb, A_tri_grnd)
      ELSEIF (TAFa <= 0._dp .AND. TAFb >= 0._dp .AND. TAFc >= 0._dp) THEN
        ! a is floating, b and c are grounded
        CALL determine_grounded_fraction_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      ELSEIF (TAFa >= 0._dp .AND. TAFb <= 0._dp .AND. TAFc >= 0._dp) THEN
        ! b is floating, c and a are grounded
        CALL determine_grounded_fraction_triangle_1flt_2grnd( vb, vc, va, TAFb, TAFc, TAFa, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      ELSEIF (TAFa >= 0._dp .AND. TAFb >= 0._dp .AND. TAFc <= 0._dp) THEN
        ! c is floating, a and b are grounded
        CALL determine_grounded_fraction_triangle_1flt_2grnd( vc, va, vb, TAFc, TAFa, TAFb, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      ELSE
        A_tri_grnd = 0._dp
        WRITE(0,*) 'determine_grounded_fraction_triangle - ERROR: TAF = [', TAFa, ',', TAFb, ',', TAFc, ']'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END IF
    
  END SUBROUTINE determine_grounded_fraction_triangle
  SUBROUTINE determine_grounded_fraction_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
    ! Determine the grounded area of the triangle [va,vb,vc], where vertex a is grounded
    ! and b and c are floating
    
    IMPLICIT NONE
    
    ! In- and output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
    REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
    REAL(dp),                            INTENT(OUT)   :: A_tri_grnd
    
    ! Local variables:
    REAL(dp)                                           :: lambda_ab, lambda_ac
    REAL(dp), DIMENSION(2)                             :: pab, pac
    
    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)
    
    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)
    
    CALL find_triangle_area( va, pab, pac, A_tri_grnd)
    
  END SUBROUTINE determine_grounded_fraction_triangle_1grnd_2flt
  SUBROUTINE determine_grounded_fraction_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_tri_flt)
    ! Determine the grounded area of the triangle [va,vb,vc], where vertex a is floating
    ! and b and c are grounded
    
    IMPLICIT NONE
    
    ! In- and output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
    REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
    REAL(dp),                            INTENT(OUT)   :: A_tri_flt
    
    ! Local variables:
    REAL(dp)                                           :: lambda_ab, lambda_ac
    REAL(dp), DIMENSION(2)                             :: pab, pac
    
    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)
    
    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)
    
    CALL find_triangle_area( va, pab, pac, A_tri_flt)
    
  END SUBROUTINE determine_grounded_fraction_triangle_1flt_2grnd
  
! == The no-ice mask, to prevent ice growth in certain areas
  SUBROUTINE initialise_mask_noice( region, mesh)
    ! Mask a certain area where no ice is allowed to grow. This is used to "remove"
    ! Greenland from NAM and EAS, and Ellesmere Island from GRL.
    ! 
    ! Also used to define calving fronts in certain idealised-geometry experiments
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Initialise
    region%mask_noice( mesh%vi1:mesh%vi2) = 0
    CALL sync
    
    IF     (region%name == 'NAM') THEN
      ! Define a no-ice mask for North America
    
      IF     (C%choice_mask_noice_NAM == 'none') THEN
        ! No no-ice mask is defined for North America
      ELSEIF (C%choice_mask_noice_NAM == 'NAM_remove_GRL') THEN
        ! Prevent ice growth in the Greenlandic part of the North America domain
        CALL initialise_mask_noice_NAM_remove_GRL( mesh, region%mask_noice)
      ELSE
        IF (par%master) WRITE(0,*) 'initialise_mask_noice - ERROR: unknown choice_mask_noice_NAM "', TRIM(C%choice_mask_noice_NAM), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSEIF (region%name == 'EAS') THEN
      ! Define a no-ice mask for Eurasia
    
      IF     (C%choice_mask_noice_EAS == 'none') THEN
        ! No no-ice mask is defined for Eurasia
      ELSEIF (C%choice_mask_noice_EAS == 'EAS_remove_GRL') THEN
        ! Prevent ice growth in the Greenlandic part of the Eurasia domain
        CALL initialise_mask_noice_EAS_remove_GRL( mesh, region%mask_noice)
      ELSE
        IF (par%master) WRITE(0,*) 'initialise_mask_noice - ERROR: unknown choice_mask_noice_EAS "', TRIM(C%choice_mask_noice_EAS), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSEIF (region%name == 'GRL') THEN
      ! Define a no-ice mask for Greenland
    
      IF     (C%choice_mask_noice_GRL == 'none') THEN
        ! No no-ice mask is defined for Greenland
      ELSEIF (C%choice_mask_noice_GRL == 'GRL_remove_Ellesmere') THEN
        ! Prevent ice growth in the Ellesmere Island part of the Greenland domain
        CALL initialise_mask_noice_GRL_remove_Ellesmere( mesh, region%mask_noice)
      ELSE
        IF (par%master) WRITE(0,*) 'initialise_mask_noice - ERROR: unknown choice_mask_noice_GRL "', TRIM(C%choice_mask_noice_GRL), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSEIF (region%name == 'ANT') THEN
      ! Define a no-ice mask for Antarctica, or for an idealised-geometry experiment
    
      IF     (C%choice_mask_noice_ANT == 'none') THEN
        ! No no-ice mask is defined for Antarctica
      ELSEIF (C%choice_mask_noice_ANT == 'MISMIP_mod') THEN
        ! Confine ice to the circular shelf around the cone-shaped island of the MISMIP_mod idealised geometry
        CALL initialise_mask_noice_MISMIP_mod( mesh, region%mask_noice)
      ELSEIF (C%choice_mask_noice_ANT == 'MISMIP+') THEN
        ! Enforce the static calving front at x = 640 km in the MISMIP+ idealised geometry
        CALL initialise_mask_noice_MISMIPplus( mesh, region%mask_noice)
      ELSE
        IF (par%master) WRITE(0,*) 'initialise_mask_noice - ERROR: unknown choice_mask_noice_ANT "', TRIM(C%choice_mask_noice_ANT), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END IF
  
  END SUBROUTINE initialise_mask_noice
  SUBROUTINE initialise_mask_noice_NAM_remove_GRL( mesh, mask_noice)
    ! Prevent ice growth in the Greenlandic part of the North America domain
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: mask_noice
  
    ! Local variables:
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(2)                             :: pa, pb
    REAL(dp)                                           :: yl_ab
      
    pa = [ 490000._dp, 1530000._dp]
    pb = [2030000._dp,  570000._dp]
    
    DO vi = mesh%vi1, mesh%vi2
      yl_ab = pa(2) + (mesh%V( vi,1) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
      IF (mesh%V( vi,2) > yl_ab .AND. mesh%V( vi,1) > pa(1) .AND. mesh%V( vi,2) > pb(2)) THEN
        mask_noice( vi) = 1
      ELSE
        mask_noice( vi) = 0
      END IF
    END DO
    CALL sync
  
  END SUBROUTINE initialise_mask_noice_NAM_remove_GRL
  SUBROUTINE initialise_mask_noice_EAS_remove_GRL( mesh, mask_noice)
    ! Prevent ice growth in the Greenlandic part of the Eurasia domain
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: mask_noice
  
    ! Local variables:
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc, pd
    REAL(dp)                                           :: yl_ab, yl_bc, yl_cd
    
    pa = [-2900000._dp, 1300000._dp]
    pb = [-1895000._dp,  900000._dp]
    pc = [ -835000._dp, 1135000._dp]
    pd = [ -400000._dp, 1855000._dp]
    
    DO vi = mesh%vi1, mesh%vi2
      yl_ab = pa(2) + (mesh%V( vi,1) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
      yl_bc = pb(2) + (mesh%V( vi,1) - pb(1))*(pc(2)-pb(2))/(pc(1)-pb(1))
      yl_cd = pc(2) + (mesh%V( vi,1) - pc(1))*(pd(2)-pc(2))/(pd(1)-pc(1))
      IF ((mesh%V( vi,1) <  pa(1) .AND. mesh%V( vi,2) > pa(2)) .OR. &
          (mesh%V( vi,1) >= pa(1) .AND. mesh%V( vi,1) < pb(1) .AND. mesh%V( vi,2) > yl_ab) .OR. &
          (mesh%V( vi,1) >= pb(1) .AND. mesh%V( vi,1) < pc(1) .AND. mesh%V( vi,2) > yl_bc) .OR. &
          (mesh%V( vi,1) >= pc(1) .AND. mesh%V( vi,1) < pd(1) .AND. mesh%V( vi,2) > yl_cd)) THEN
        mask_noice( vi) = 1
      ELSE
        mask_noice( vi) = 0
      END IF
    END DO
    CALL sync
  
  END SUBROUTINE initialise_mask_noice_EAS_remove_GRL
  SUBROUTINE initialise_mask_noice_GRL_remove_Ellesmere( mesh, mask_noice)
    ! Prevent ice growth in the Ellesmere Island part of the Greenland domain
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: mask_noice
  
    ! Local variables:
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(2)                             :: pa, pb
    REAL(dp)                                           :: yl_ab
      
    pa = [-750000._dp,  900000._dp]
    pb = [-250000._dp, 1250000._dp]
    
    DO vi = mesh%vi1, mesh%vi2
      yl_ab = pa(2) + (mesh%V( vi,1) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
      IF (mesh%V( vi,2) > pa(2) .AND. mesh%V( vi,2) > yl_ab .AND. mesh%V( vi,1) < pb(1)) THEN
        mask_noice( vi) = 1
      ELSE
        mask_noice( vi) = 0
      END IF
    END DO
    CALL sync
  
  END SUBROUTINE initialise_mask_noice_GRL_remove_Ellesmere
  SUBROUTINE initialise_mask_noice_MISMIP_mod( mesh, mask_noice)
    ! Confine ice to the circular shelf around the cone-shaped island of the MISMIP_mod idealised-geometry experiment
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: mask_noice
  
    ! Local variables:
    INTEGER                                            :: vi
        
    ! Create a nice circular ice shelf
    DO vi = mesh%vi1, mesh%vi2
      IF (SQRT(mesh%V( vi,1)**2 + mesh%V( vi,2)**2) > mesh%xmax * 0.95_dp) THEN
        mask_noice( vi) = 1
      ELSE
        mask_noice( vi) = 0
      END IF
    END DO
    CALL sync
  
  END SUBROUTINE initialise_mask_noice_MISMIP_mod
  SUBROUTINE initialise_mask_noice_MISMIPplus( mesh, mask_noice)
    ! Enforce the static calving front at x = 640 km in the MISMIP+ idealised geometry
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: mask_noice
  
    ! Local variables:
    INTEGER                                            :: vi
        
    DO vi = mesh%vi1, mesh%vi2
      ! NOTE: because UFEMISM wants to centre the domain at x=0, the front now lies at x = 240 km
      IF (mesh%V( vi,1) > 240000._dp) THEN
        mask_noice( vi) = 1
      ELSE
        mask_noice( vi) = 0
      END IF
    END DO
    CALL sync
  
  END SUBROUTINE initialise_mask_noice_MISMIPplus

END MODULE general_ice_model_data_module
