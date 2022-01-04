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
  USE data_types_module,               ONLY: type_mesh, type_ice_model
  USE utilities_module,                ONLY: is_floating, surface_elevation, thickness_above_floatation
  USE mesh_help_functions_module,      ONLY: find_triangle_area
  USE mesh_operators_module,           ONLY: map_a_to_ac_2D, map_ac_to_bb_2D

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
        ! Ice (sheet or shelf) bordering ocean equals calvingfront
        
        DO ci = 1, mesh%nC(vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ocean_a( vc) == 1) THEN
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
    REAL(dp), DIMENSION(:    ), POINTER                ::  TAF_ac,  TAF_bb
    INTEGER                                            :: wTAF_ac, wTAF_bb
    INTEGER                                            :: avi, ci, avj, iati, iati2, ati1, ati2
    REAL(dp)                                           :: TAF_max, TAF_min
    REAL(dp), DIMENSION(2)                             :: vac, ccbb1, ccbb2
    REAL(dp)                                           :: TAFa, TAFb, TAFc, A_vor_ac, A_subtri_bb, A_grnd_ac, A_grnd_bb
    
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
      A_vor_ac  = 0._dp
      A_grnd_ac = 0._dp
      
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
        CALL determine_grounded_fraction_subtri_bb( vac, ccbb1, ccbb2, TAFa, TAFb, TAFc, A_subtri_bb, A_grnd_bb)
        
        A_vor_ac  = A_vor_ac  + A_subtri_bb
        A_grnd_ac = A_grnd_ac + A_grnd_bb
        
      END DO ! DO iati = 1, mesh%niTriAaAc( avi)
      
      ! Calculate the grounded fraction of this Voronoi cell
      ice%f_grnd_ac( avi) = A_grnd_ac / A_vor_ac
      
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
    
  END SUBROUTINE determine_grounded_fractions_ac
  SUBROUTINE determine_grounded_fraction_subtri_bb( va, vb, vc, TAFa, TAFb, TAFc, A_subtri, A_grnd)
    ! Determine the grounded area of the triangle [va,vb,vc], where the thickness-above-floatation is given at all three corners
    
    IMPLICIT NONE
    
    ! In- and output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
    REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
    REAL(dp),                            INTENT(OUT)   :: A_subtri, A_grnd
    
    ! Local variables:
    REAL(dp)                                           :: A_flt
    
    ! Determine total area of this subtriangle
    CALL find_triangle_area( va, vb, vc, A_subtri)
        
    IF     (TAFa >= 0._dp .AND. TAFb >= 0._dp .AND. TAFc >= 0._dp) THEN
      ! If all three corners are grounded, the answer is trivial
      A_grnd = A_subtri
    ELSEIF (TAFa <= 0._dp .AND. TAFb <= 0._dp .AND. TAFc <= 0._dp) THEN
      ! If all three corners are floating, the answer is trivial
      A_grnd = 0._dp
    ELSE
      ! At least one corner is grounded and at least one corner is floating
      
      IF     (TAFa >= 0._dp .AND. TAFb <= 0._dp .AND. TAFc <= 0._dp) THEN
        ! a is grounded, b and c are floating
        CALL determine_grounded_fraction_subtri_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_grnd)
      ELSEIF (TAFa <= 0._dp .AND. TAFb >= 0._dp .AND. TAFc <= 0._dp) THEN
        ! b is grounded, a and c are floating
        CALL determine_grounded_fraction_subtri_1grnd_2flt( vb, vc, va, TAFb, TAFc, TAFa, A_grnd)
      ELSEIF (TAFa <= 0._dp .AND. TAFb <= 0._dp .AND. TAFc >= 0._dp) THEN
        ! c is grounded, a and b are floating
        CALL determine_grounded_fraction_subtri_1grnd_2flt( vc, va, vb, TAFc, TAFa, TAFb, A_grnd)
      ELSEIF (TAFa <= 0._dp .AND. TAFb >= 0._dp .AND. TAFc >= 0._dp) THEN
        ! a is floating, b and c are grounded
        CALL determine_grounded_fraction_subtri_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_flt)
        A_grnd = A_subtri - A_flt
      ELSEIF (TAFa >= 0._dp .AND. TAFb <= 0._dp .AND. TAFc >= 0._dp) THEN
        ! b is floating, c and a are grounded
        CALL determine_grounded_fraction_subtri_1flt_2grnd( vb, vc, va, TAFb, TAFc, TAFa, A_flt)
        A_grnd = A_subtri - A_flt
      ELSEIF (TAFa >= 0._dp .AND. TAFb >= 0._dp .AND. TAFc <= 0._dp) THEN
        ! c is floating, a and b are grounded
        CALL determine_grounded_fraction_subtri_1flt_2grnd( vc, va, vb, TAFc, TAFa, TAFb, A_flt)
        A_grnd = A_subtri - A_flt
      ELSE
        A_grnd = 0._dp
        WRITE(0,*) 'determine_grounded_fraction_subtri_bb - ERROR: TAF = [', TAFa, ',', TAFb, ',', TAFc, ']'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END IF
    
  END SUBROUTINE determine_grounded_fraction_subtri_bb
  SUBROUTINE determine_grounded_fraction_subtri_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_grnd)
    ! Determine the grounded area of the triangle [va,vb,vc], where vertex a is grounded
    ! and b and c are floating
    
    IMPLICIT NONE
    
    ! In- and output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
    REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
    REAL(dp),                            INTENT(OUT)   :: A_grnd
    
    ! Local variables:
    REAL(dp)                                           :: lambda_ab, lambda_ac
    REAL(dp), DIMENSION(2)                             :: pab, pac
    
    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)
    
    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)
    
    CALL find_triangle_area( va, pab, pac, A_grnd)
    
  END SUBROUTINE determine_grounded_fraction_subtri_1grnd_2flt
  SUBROUTINE determine_grounded_fraction_subtri_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_flt)
    ! Determine the grounded area of the triangle [va,vb,vc], where vertex a is floating
    ! and b and c are grounded
    
    IMPLICIT NONE
    
    ! In- and output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
    REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
    REAL(dp),                            INTENT(OUT)   :: A_flt
    
    ! Local variables:
    REAL(dp)                                           :: lambda_ab, lambda_ac
    REAL(dp), DIMENSION(2)                             :: pab, pac
    
    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)
    
    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)
    
    CALL find_triangle_area( va, pab, pac, A_flt)
    
  END SUBROUTINE determine_grounded_fraction_subtri_1flt_2grnd

END MODULE general_ice_model_data_module
