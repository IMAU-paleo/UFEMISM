MODULE general_ice_model_data_module

  ! Only "secondary geometry" right now: masks, surface elevation, TAF

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list

  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_model_region
  USE utilities_module,                ONLY: is_floating, surface_elevation, thickness_above_floatation
  use mpi_module,                      only: allgather_array
  ! USE mesh_help_functions_module,      ONLY: find_triangle_area
  ! USE mesh_operators_module,           ONLY: map_a_to_b_2D

  IMPLICIT NONE

CONTAINS
  
  ! Routines for calculating general ice model data - Hs, masks, ice physical properties
  SUBROUTINE update_general_ice_model_data( mesh, ice)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_general_ice_model_data'
    INTEGER                                            :: vi
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Calculate surface elevation and thickness above floatation
    DO vi = mesh%vi1, mesh%vi2
      ice%Hs_a(  vi) = surface_elevation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))
      ice%TAF_a( vi) = thickness_above_floatation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))
    END DO
    
    ! Determine masks
    CALL determine_masks( mesh, ice)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE update_general_ice_model_data
  SUBROUTINE determine_masks( mesh, ice)
    ! Determine the different masks, on both the Aa and the Ac mesh
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_masks'
    INTEGER                                            :: vi, ci, vc
    real(dp), dimension(:), allocatable                :: Hi_a, Hb_a, SL_a
    
    ! Add routine to path
    CALL init_routine( routine_name)

    
    ! Get necessary information
    allocate(Hi_a(1:mesh%nV))
    allocate(Hb_a(1:mesh%nV))
    allocate(SL_a(1:mesh%nV))
    Hi_a(mesh%vi1:mesh%vi2) = ice%Hi_a
    Hb_a(mesh%vi1:mesh%vi2) = ice%Hb_a
    SL_a(mesh%vi1:mesh%vi2) = ice%SL_a
    call allgather_array(Hi_a)
    call allgather_array(Hb_a)
    call allgather_array(SL_a)
    

    ! Start out with land everywhere, fill in the rest based on input.
    ice%mask_land_a   = 1
    ice%mask_lake_a   = 0
    ice%mask_ocean_a  = 0
    ice%mask_ice_a    = 0
    ice%mask_shelf_a  = 0
    ice%mask_sheet_a  = 0
    ice%mask_coast_a  = 0
    ice%mask_margin_a = 0
    ice%mask_gl_a     = 0
    ice%mask_cf_a     = 0
    ice%mask_a        = C%type_land
  
    DO vi = mesh%vi1, mesh%vi2
      
      ! Determine ocean (both open and shelf-covered)
      IF (is_floating( Hi_a( vi), Hb_a( vi), SL_a( vi))) THEN
        ice%mask_ocean_a( vi) = 1
        ice%mask_land_a(  vi) = 0
        ice%mask_a(       vi) = C%type_ocean
      END IF
    
      ! Determine ice
      IF (Hi_a( vi) > 0._dp) THEN
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
      
    END DO ! DO vi = 1, mesh%nV
  
    ! Determine coast, grounding line and calving front
    DO vi = mesh%vi1, mesh%vi2
    
      IF (ice%mask_land_a( vi) == 1) THEN  
        ! Land bordering ocean equals coastline
        
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          ! if neighbour is sea
          IF (is_floating( Hi_a( vc), Hb_a( vc), SL_a( vc))) THEN
            ice%mask_a( vi) = C%type_coast
            ice%mask_coast_a( vi) =  1
          END IF
        END DO
      END IF
    
      IF (ice%mask_ice_a( vi) == 1) THEN  
        ! Ice bordering non-ice equals margin
        
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          ! if neighbour has no ice
          IF (Hi_a( vc) <= 0._dp) THEN
            ice%mask_a( vi) = C%type_margin
            ice%mask_margin_a( vi) =  1
          END IF
        END DO
      END IF
    
      IF (ice%mask_sheet_a( vi) == 1) THEN  
        ! Sheet bordering shelf equals groundingline
        
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          ! if neighbour has ocean and ice
          IF (is_floating( Hi_a( vc), Hb_a( vc), SL_a( vc)) .and. Hi_a( vc) > 0._dp) THEN
            ice%mask_a( vi) = C%type_groundingline
            ice%mask_gl_a( vi) = 1
          END IF
        END DO
      END IF
  
      IF (ice%mask_ice_a( vi) == 1) THEN  
        ! Ice (sheet or shelf) bordering open ocean equals calvingfront
        
        DO ci = 1, mesh%nC(vi)
          vc = mesh%C( vi,ci)
          ! if neighbour has ocean and no ice
          IF (is_floating( Hi_a( vc), Hb_a( vc), SL_a( vc)) .and. Hi_a( vc) <= 0._dp) THEN
            ice%mask_a( vi) = C%type_calvingfront
            ice%mask_cf_a( vi) = 1
          END IF
        END DO
        
      END IF
    END DO ! DO vi = 1, mesh%nV

    deallocate(Hi_a)
    deallocate(Hb_a)
    deallocate(SL_a)

    call allgather_array(ice%mask_land_a)
    call allgather_array(ice%mask_lake_a  )
    call allgather_array(ice%mask_ocean_a )
    call allgather_array(ice%mask_ice_a   )
    call allgather_array(ice%mask_shelf_a )
    call allgather_array(ice%mask_sheet_a )
    call allgather_array(ice%mask_coast_a )
    call allgather_array(ice%mask_margin_a)
    call allgather_array(ice%mask_gl_a    )
    call allgather_array(ice%mask_cf_a    )
    call allgather_array(ice%mask_a       )

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE determine_masks
  
! ! == Routines for calculating sub-grid grounded fractions
!   SUBROUTINE determine_grounded_fractions( mesh, ice)
!     ! Determine the grounded fractions of all grid cells
    
!     IMPLICIT NONE
    
!     ! In- and output variables
!     TYPE(type_mesh),                     INTENT(IN)    :: mesh
!     TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
!     ! Local variables:
!     CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_grounded_fractions'
    
!     ! Add routine to path
!     CALL init_routine( routine_name)
    
!     CALL determine_grounded_fractions_a( mesh, ice)
!     CALL determine_grounded_fractions_b( mesh, ice)
    
!     ! Finalise routine path
!     CALL finalise_routine( routine_name)
    
!   END SUBROUTINE determine_grounded_fractions
!   SUBROUTINE determine_grounded_fractions_a( mesh, ice)
!     ! Determine the grounded fractions of all grid cells on the a-grid
    
!     IMPLICIT NONE
    
!     ! In- and output variables
!     TYPE(type_mesh),                     INTENT(IN)    :: mesh
!     TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
!     ! Local variables:
!     CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_grounded_fractions_a'
!     REAL(dp), DIMENSION(:    ), POINTER                ::  TAF_b
!     INTEGER                                            :: wTAF_b
!     INTEGER                                            :: vi, ci, vj, iti, iti2, ti1, ti2
!     REAL(dp)                                           :: TAF_max, TAF_min
!     REAL(dp), DIMENSION(2)                             :: va, ccb1, ccb2
!     REAL(dp)                                           :: TAFa, TAFb, TAFc, A_vor, A_tri_tot, A_tri_grnd, A_grnd
    
!     ! Add routine to path
!     CALL init_routine( routine_name)
    
!     ! Map thickness-above-floatation to the b-grid
!     CALL allocate_shared_dp_1D( mesh%nTri, TAF_b, wTAF_b)
!     CALL map_a_to_b_2D(  mesh, ice%TAF_a, TAF_b)
  
!     DO vi = mesh%vi1, mesh%vi2
      
!       ! Skip border vertices
!       IF (mesh%edge_index( vi) > 0) THEN
!         ice%f_grnd_a( vi) = 0._dp
!         CYCLE
!       END IF
      
!       ! Determine maximum and minimum TAF of the local neighbourhood
!       TAF_max = -1E6_dp
!       TAF_min =  1E6_dp
      
!       TAF_max = MAX( TAF_max, ice%TAF_a( vi))
!       TAF_min = MIN( TAF_min, ice%TAF_a( vi))
      
!       DO ci = 1, mesh%nC( vi)
!         vj = mesh%C( vi,ci)
!         TAF_max = MAX( TAF_max, ice%TAF_a( vj))
!         TAF_min = MIN( TAF_min, ice%TAF_a( vj))
!       END DO
      
!       ! If the entire local neighbourhood is grounded, the answer is trivial
!       IF (TAF_min >= 0._dp) THEN
!         ice%f_grnd_a( vi) = 1._dp
!         CYCLE
!       END IF
      
!       ! If the entire local neighbourhood is floating, the answer is trivial
!       IF (TAF_max <= 0._dp) THEN
!         ice%f_grnd_a( vi) = 0._dp
!         CYCLE
!       END IF
      
!       ! The local neighbourhood contains both grounded and floating vertices.
!       A_vor  = 0._dp
!       A_grnd = 0._dp
      
!       va   = mesh%V( vi,:)
!       TAFa = ice%TAF_a( vi)
      
!       DO iti = 1, mesh%niTri( vi)
        
!         iti2 = iti + 1
!         IF (iti == mesh%niTri( vi)) iti2 = 1
        
!         ti1 = mesh%iTri( vi,iti )
!         ti2 = mesh%iTri( vi,iti2)
        
!         ccb1 = mesh%Tricc( ti1,:)
!         ccb2 = mesh%Tricc( ti2,:)
        
!         TAFb = TAF_b( ti1)
!         TAFc = TAF_b( ti2)
        
!         ! Determine total area of, and grounded area within, this subtriangle
!         CALL determine_grounded_area_triangle( va, ccb1, ccb2, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
        
!         A_vor  = A_vor  + A_tri_tot
!         A_grnd = A_grnd + A_tri_grnd
        
!       END DO ! DO iati = 1, mesh%niTriAaAc( avi)
      
!       ! Calculate the grounded fraction of this Voronoi cell
!       ice%f_grnd_a( vi) = A_grnd / A_vor
      
!     END DO
!     CALL sync
    
!     ! Clean up after yourself
!     CALL deallocate_shared( wTAF_b)
    
!     ! Finalise routine path
!     CALL finalise_routine( routine_name)
    
!   END SUBROUTINE determine_grounded_fractions_a
!   SUBROUTINE determine_grounded_fractions_b( mesh, ice)
!     ! Determine the grounded fractions of all grid cells on the b-grid
    
!     IMPLICIT NONE
    
!     ! In- and output variables
!     TYPE(type_mesh),                     INTENT(IN)    :: mesh
!     TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
!     ! Local variables:
!     CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_grounded_fractions_b'
!     INTEGER                                            :: ti, via, vib, vic
!     REAL(dp)                                           :: TAF_max, TAF_min
!     REAL(dp), DIMENSION(2)                             :: va, vb, vc
!     REAL(dp)                                           :: TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd
    
!     ! Add routine to path
!     CALL init_routine( routine_name)
  
!     DO ti = mesh%ti1, mesh%ti2
      
!       via = mesh%Tri( ti,1)
!       vib = mesh%Tri( ti,2)
!       vic = mesh%Tri( ti,3)
      
!       ! Determine maximum and minimum TAF of the local neighbourhood
!       TAF_max = MAXVAL([ ice%TAF_a( via), ice%TAF_a( vib), ice%TAF_a( vic)])
!       TAF_min = MINVAL([ ice%TAF_a( via), ice%TAF_a( vib), ice%TAF_a( vic)])
      
!       ! If the entire local neighbourhood is grounded, the answer is trivial
!       IF (TAF_min >= 0._dp) THEN
!         ice%f_grnd_b( ti) = 1._dp
!         CYCLE
!       END IF
      
!       ! If the entire local neighbourhood is floating, the answer is trivial
!       IF (TAF_max <= 0._dp) THEN
!         ice%f_grnd_b( ti) = 0._dp
!         CYCLE
!       END IF
      
!       ! The local neighbourhood contains both grounded and floating vertices.
      
!       va   = mesh%V( via,:)
!       vb   = mesh%V( vib,:)
!       vc   = mesh%V( vic,:)
        
!       TAFa = ice%TAF_a( via)
!       TAFb = ice%TAF_a( vib)
!       TAFc = ice%TAF_a( vic)
        
!       ! Determine total area of, and grounded area within, this subtriangle
!       CALL determine_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
      
!       ! Calculate the grounded fraction of this Voronoi cell
!       ice%f_grnd_b( ti) = A_tri_grnd / A_tri_tot
      
!     END DO
!     CALL sync
    
!     ! Finalise routine path
!     CALL finalise_routine( routine_name)
    
!   END SUBROUTINE determine_grounded_fractions_b
!   SUBROUTINE determine_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
!     ! Determine the grounded area of the triangle [va,vb,vc], where the thickness-above-floatation is given at all three corners
    
!     IMPLICIT NONE
    
!     ! In- and output variables
!     REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
!     REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
!     REAL(dp),                            INTENT(OUT)   :: A_tri_tot, A_tri_grnd
    
!     ! Local variables:
!     REAL(dp)                                           :: A_flt
    
!     ! Determine total area of this subtriangle
!     CALL find_triangle_area( va, vb, vc, A_tri_tot)
        
!     IF     (TAFa >= 0._dp .AND. TAFb >= 0._dp .AND. TAFc >= 0._dp) THEN
!       ! If all three corners are grounded, the answer is trivial
!       A_tri_grnd = A_tri_tot
!     ELSEIF (TAFa <= 0._dp .AND. TAFb <= 0._dp .AND. TAFc <= 0._dp) THEN
!       ! If all three corners are floating, the answer is trivial
!       A_tri_grnd = 0._dp
!     ELSE
!       ! At least one corner is grounded and at least one corner is floating
      
!       IF     (TAFa >= 0._dp .AND. TAFb <= 0._dp .AND. TAFc <= 0._dp) THEN
!         ! a is grounded, b and c are floating
!         CALL determine_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
!       ELSEIF (TAFa <= 0._dp .AND. TAFb >= 0._dp .AND. TAFc <= 0._dp) THEN
!         ! b is grounded, a and c are floating
!         CALL determine_grounded_area_triangle_1grnd_2flt( vb, vc, va, TAFb, TAFc, TAFa, A_tri_grnd)
!       ELSEIF (TAFa <= 0._dp .AND. TAFb <= 0._dp .AND. TAFc >= 0._dp) THEN
!         ! c is grounded, a and b are floating
!         CALL determine_grounded_area_triangle_1grnd_2flt( vc, va, vb, TAFc, TAFa, TAFb, A_tri_grnd)
!       ELSEIF (TAFa <= 0._dp .AND. TAFb >= 0._dp .AND. TAFc >= 0._dp) THEN
!         ! a is floating, b and c are grounded
!         CALL determine_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_flt)
!         A_tri_grnd = A_tri_tot - A_flt
!       ELSEIF (TAFa >= 0._dp .AND. TAFb <= 0._dp .AND. TAFc >= 0._dp) THEN
!         ! b is floating, c and a are grounded
!         CALL determine_grounded_area_triangle_1flt_2grnd( vb, vc, va, TAFb, TAFc, TAFa, A_flt)
!         A_tri_grnd = A_tri_tot - A_flt
!       ELSEIF (TAFa >= 0._dp .AND. TAFb >= 0._dp .AND. TAFc <= 0._dp) THEN
!         ! c is floating, a and b are grounded
!         CALL determine_grounded_area_triangle_1flt_2grnd( vc, va, vb, TAFc, TAFa, TAFb, A_flt)
!         A_tri_grnd = A_tri_tot - A_flt
!       ELSE
!         A_tri_grnd = 0._dp
!         CALL crash('TAF = [{dp_01},{dp_02},{dp_03}]', dp_01 = TAFa, dp_02 = TAFb, dp_03 = TAFc)
!       END IF
      
!     END IF
    
!   END SUBROUTINE determine_grounded_area_triangle
!   SUBROUTINE determine_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
!     ! Determine the grounded area of the triangle [va,vb,vc], where vertex a is grounded
!     ! and b and c are floating
    
!     IMPLICIT NONE
    
!     ! In- and output variables
!     REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
!     REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
!     REAL(dp),                            INTENT(OUT)   :: A_tri_grnd
    
!     ! Local variables:
!     REAL(dp)                                           :: lambda_ab, lambda_ac
!     REAL(dp), DIMENSION(2)                             :: pab, pac
    
!     lambda_ab = TAFa / (TAFa - TAFb)
!     pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)
    
!     lambda_ac = TAFa / (TAFa - TAFc)
!     pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)
    
!     CALL find_triangle_area( va, pab, pac, A_tri_grnd)
    
!   END SUBROUTINE determine_grounded_area_triangle_1grnd_2flt
!   SUBROUTINE determine_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_tri_flt)
!     ! Determine the grounded area of the triangle [va,vb,vc], where vertex a is floating
!     ! and b and c are grounded
    
!     IMPLICIT NONE
    
!     ! In- and output variables
!     REAL(dp), DIMENSION(2),              INTENT(IN)    :: va, vb, vc
!     REAL(dp),                            INTENT(IN)    :: TAFa, TAFb, TAFc
!     REAL(dp),                            INTENT(OUT)   :: A_tri_flt
    
!     ! Local variables:
!     REAL(dp)                                           :: lambda_ab, lambda_ac
!     REAL(dp), DIMENSION(2)                             :: pab, pac
    
!     lambda_ab = TAFa / (TAFa - TAFb)
!     pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)
    
!     lambda_ac = TAFa / (TAFa - TAFc)
!     pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)
    
!     CALL find_triangle_area( va, pab, pac, A_tri_flt)
    
!   END SUBROUTINE determine_grounded_area_triangle_1flt_2grnd
  
! ! == The no-ice mask, to prevent ice growth in certain areas
  SUBROUTINE initialise_mask_noice( region, mesh)
    ! Mask a certain area where no ice is allowed to grow. This is used to "remove"
    ! Greenland from NAM and EAS, and Ellesmere Island from GRL.
    !
    ! Also used to define calving fronts in certain idealised-geometry experiments
    
    IMPLICIT NONE
  
    ! In- and output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
  
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_mask_noice'
  
    ! Add routine to path
    CALL init_routine( routine_name)
  
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
        CALL crash('unknown choice_mask_noice_NAM "' // TRIM( C%choice_mask_noice_NAM) // '"!')
      END IF
    
    ELSEIF (region%name == 'EAS') THEN
      ! Define a no-ice mask for Eurasia
  
      IF     (C%choice_mask_noice_EAS == 'none') THEN
        ! No no-ice mask is defined for Eurasia
      ELSEIF (C%choice_mask_noice_EAS == 'EAS_remove_GRL') THEN
        ! Prevent ice growth in the Greenlandic part of the Eurasia domain
        CALL initialise_mask_noice_EAS_remove_GRL( mesh, region%mask_noice)
      ELSE
        CALL crash('unknown choice_mask_noice_EAS "' // TRIM( C%choice_mask_noice_EAS) // '"!')
      END IF
    
    ELSEIF (region%name == 'GRL') THEN
      ! Define a no-ice mask for Greenland
  
      IF     (C%choice_mask_noice_GRL == 'none') THEN
        ! No no-ice mask is defined for Greenland
      ELSEIF (C%choice_mask_noice_GRL == 'GRL_remove_Ellesmere') THEN
        ! Prevent ice growth in the Ellesmere Island part of the Greenland domain
        CALL initialise_mask_noice_GRL_remove_Ellesmere( mesh, region%mask_noice)
      ELSE
        CALL crash('unknown choice_mask_noice_GRL "' // TRIM( C%choice_mask_noice_GRL) // '"!')
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
        CALL crash('unknown choice_mask_noice_ANT "' // TRIM( C%choice_mask_noice_ANT) // '"!')
      END IF
    
    END IF
  
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_mask_noice
  SUBROUTINE initialise_mask_noice_NAM_remove_GRL( mesh, mask_noice)
    ! Prevent ice growth in the Greenlandic part of the North America domain
    
    IMPLICIT NONE
  
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_mask_noice_NAM_remove_GRL'
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(2)                             :: pa, pb
    REAL(dp)                                           :: yl_ab
  
    ! Add routine to path
    CALL init_routine( routine_name)
    
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
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_mask_noice_NAM_remove_GRL
  SUBROUTINE initialise_mask_noice_EAS_remove_GRL( mesh, mask_noice)
    ! Prevent ice growth in the Greenlandic part of the Eurasia domain
    
    IMPLICIT NONE
  
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_mask_noice_EAS_remove_GRL'
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc, pd
    REAL(dp)                                           :: yl_ab, yl_bc, yl_cd
  
    ! Add routine to path
    CALL init_routine( routine_name)
  
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
  
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_mask_noice_EAS_remove_GRL
  SUBROUTINE initialise_mask_noice_GRL_remove_Ellesmere( mesh, mask_noice)
    ! Prevent ice growth in the Ellesmere Island part of the Greenland domain
    
    IMPLICIT NONE
  
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_mask_noice_GRL_remove_Ellesmere'
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(2)                             :: pa, pb
    REAL(dp)                                           :: yl_ab
  
    ! Add routine to path
    CALL init_routine( routine_name)
    
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
  
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_mask_noice_GRL_remove_Ellesmere
  SUBROUTINE initialise_mask_noice_MISMIP_mod( mesh, mask_noice)
    ! Confine ice to the circular shelf around the cone-shaped island of the MISMIP_mod idealised-geometry experiment
    
    IMPLICIT NONE
  
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_mask_noice_MISMIP_mod'
    INTEGER                                            :: vi
  
    ! Add routine to path
    CALL init_routine( routine_name)
      
    ! Create a nice circular ice shelf
    DO vi = mesh%vi1, mesh%vi2
      IF (SQRT(mesh%V( vi,1)**2 + mesh%V( vi,2)**2) > mesh%xmax * 0.95_dp) THEN
        mask_noice( vi) = 1
      ELSE
        mask_noice( vi) = 0
      END IF
    END DO

  
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_mask_noice_MISMIP_mod
  SUBROUTINE initialise_mask_noice_MISMIPplus( mesh, mask_noice)
    ! Enforce the static calving front at x = 640 km in the MISMIP+ idealised geometry
    
    IMPLICIT NONE
  
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_mask_noice_MISMIPplus'
    INTEGER                                            :: vi
  
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO vi = mesh%vi1, mesh%vi2
      ! NOTE: because UFEMISM wants to centre the domain at x=0, the front now lies at x = 240 km
      IF (mesh%V( vi,1) > 240000._dp) THEN
        mask_noice( vi) = 1
      ELSE
        mask_noice( vi) = 0
      END IF
    END DO
  
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_mask_noice_MISMIPplus

  ! == Routines for defining ice drainage basins from an external polygon file
  SUBROUTINE initialise_basins( mesh, basin_ID, nbasins, region_name)
    ! Define the ice basins mask from an external text file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),          INTENT(INOUT) :: basin_ID
    INTEGER,                             INTENT(INOUT) :: nbasins
    CHARACTER(LEN=3)                                   :: region_name

    ! Local variables:
    INTEGER                                            :: i,j,vi,vj,bi
    CHARACTER(LEN=256)                                 :: choice_basin_scheme
    CHARACTER(LEN=256)                                 :: filename_basins
    LOGICAL                                            :: recognised_file
    INTEGER                                            :: n_skip
    INTEGER                                            :: n_header_lines, n_vertices
    REAL(dp), DIMENSION(:    ), POINTER                :: Vlat, Vlon, Vx, Vy
    INTEGER,  DIMENSION(:    ), POINTER                :: Vid
    INTEGER                                            :: wVlat, wVlon, wVx, wVy, wVid
    REAL(dp)                                           :: VID_dp
    INTEGER                                            :: ios, dummy
    INTEGER                                            :: vi1, vi2
    REAL(dp), DIMENSION(:,:  ), POINTER                :: poly_bi
    INTEGER                                            :: wpoly_bi
    REAL(dp)                                           :: xmin,xmax,ymin,ymax
    REAL(dp), DIMENSION(2)                             :: p
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  basin_ID_loc,  basin_ID_ext_loc
    INTEGER                                            :: wbasin_ID_loc, wbasin_ID_ext_loc
    INTEGER                                            :: n_front
    INTEGER,  DIMENSION(:,:  ), POINTER                :: ij_front
    INTEGER                                            :: wij_front
    INTEGER                                            :: ii,jj,k,i_nearest,j_nearest
    REAL(dp)                                           :: dist,dist_min

    ! Determine what to do for this region
    choice_basin_scheme = 'none'
    filename_basins     = ''

    IF (region_name == 'NAM') THEN
      choice_basin_scheme = C%choice_basin_scheme_NAM
      filename_basins     = C%filename_basins_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_basin_scheme = C%choice_basin_scheme_EAS
      filename_basins     = C%filename_basins_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_basin_scheme = C%choice_basin_scheme_GRL
      filename_basins     = C%filename_basins_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_basin_scheme = C%choice_basin_scheme_ANT
      filename_basins     = C%filename_basins_ANT
    END IF

    IF     (choice_basin_scheme == 'none') THEN
      ! No basins are defined (i.e. the whole region is one big, big basin)

      basin_ID( mesh%vi1:mesh%vi2) = 1
      IF (par%master) nbasins = 1
      CALL sync

    ELSEIF (choice_basin_scheme == 'file') THEN
      ! Define basins from an external text file describing the polygons

      IF (par%master) WRITE(0,*) '  Reading basins for ', TRIM(region_name), ' from file "', TRIM(filename_basins), '"'

      ! WIP wall
      IF (par%master) WRITE(0,*) 'This subroutine (initialise_basins) does not work for multiple basins yet.'
      IF (par%master) WRITE(0,*) 'Use 1 basin (option "choice_basin_scheme=none"), or feel free to implement it yoself before running the model :)'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

    !   IF (region_name == 'ANT') THEN
    !     ! Antarctica: ant_full_drainagesystem_polygons.txt
    !     ! Can be downloaded from: https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems
    !     ! A text file with 7 header lines, followed by three columns of data:
    !     !    Lat, Lon, basin ID

    !   ! ===== Check if this is really the file we're reading =====
    !   ! ==========================================================

    !     recognised_file = .FALSE.
    !     n_header_lines  = 0
    !     n_vertices      = 0
    !     n_skip          = 0

    !     DO i = 1, 256-36
    !       IF (filename_basins(i:i+35) == 'ant_full_drainagesystem_polygons.txt') THEN
    !         recognised_file = .TRUE.
    !         n_header_lines  = 7
    !         n_skip          = 5 ! Since the polygon file is at a ridiculously high resolution, downscale it a bit for efficiency
    !         n_vertices      = CEILING( REAL(901322,dp) / REAL(n_skip,dp))
    !       END IF
    !     END DO

    !     IF ((.NOT. recognised_file) .AND. par%master) THEN
    !       WRITE(0,*) ''
    !       WRITE(0,*) ' ===== '
    !       WRITE(0,*) 'initialise_basins - WARNING: for Antarctica, we expect the file "ant_full_drainagesystem_polygons.txt"'
    !       WRITE(0,*) '                             This can be downloaded from https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems'
    !       WRITE(0,*) '                             If another file is used, make sure that it has 7 header lines, or alternatively just change the code in initialise_basins'
    !       WRITE(0,*) ' ===== '
    !       WRITE(0,*) ''
    !     END IF

    !     ! Allocate shared memory
    !     CALL allocate_shared_dp_1D(  n_vertices, Vlat, wVlat)
    !     CALL allocate_shared_dp_1D(  n_vertices, Vlon, wVlon)
    !     CALL allocate_shared_dp_1D(  n_vertices, Vx,   wVx  )
    !     CALL allocate_shared_dp_1D(  n_vertices, Vy,   wVy  )
    !     CALL allocate_shared_int_1D( n_vertices, VID,  wVID )

    !   ! ===== Read the file =====
    !   ! =========================

    !     IF (par%master) THEN

    !       ! Open the file
    !       OPEN( UNIT = 1337, FILE=filename_basins, ACTION='READ')

    !       ! Skip the header lines
    !       DO i = 1, n_header_lines
    !         READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
    !       END DO

    !       ! Read the actual data
    !       DO vi = 1, n_vertices
    !         READ( UNIT = 1337, FMT=*, IOSTAT=ios) Vlat( vi), Vlon( vi), VID( vi)
    !         IF (ios /= 0) THEN
    !           WRITE(0,*) ' initialise_basins - ERROR: length of text file "', TRIM( filename_basins), '" does not match n_vertices = ', n_vertices
    !           CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    !         END IF

    !         DO vj = 1, n_skip-1
    !           READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
    !         END DO
    !       END DO

    !       ! Close the file
    !       CLOSE( UNIT = 1337)

    !     END IF ! IF (par%master) THEN
    !     CALL sync

    !     ! Project [lat,lon] to [x,y]
    !     CALL partition_list( n_vertices, par%i, par%n, vi1, vi2)
    !     DO vi = vi1, vi2
    !       CALL oblique_sg_projection( Vlon( vi), Vlat( vi), grid%lambda_M, grid%phi_M, grid%alpha_stereo, Vx( vi), Vy( vi))
    !     END DO

    !   ELSEIF (region_name == 'GRL') THEN
    !     ! Greenland: grndrainagesystems_ekholm.txt
    !     ! Can be downloaded from: https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems
    !     ! A text file with 7 header lines, followed by three columns of data:
    !     !    basin ID, Lat, Lon   (NOTE: here basin ID is a floating-point rather than an integer!)

    !   ! ===== Check if this is really the file we're reading =====
    !   ! ==========================================================

    !     recognised_file = .FALSE.
    !     n_header_lines  = 0
    !     n_vertices      = 0
    !     n_skip          = 0

    !     DO i = 1, 256-29
    !       IF (filename_basins(i:i+28) == 'grndrainagesystems_ekholm.txt') THEN
    !         recognised_file = .TRUE.
    !         n_header_lines  = 7
    !         n_skip          = 5 ! Since the polygon file is at a ridiculously high resolution, downscale it a bit for efficiency
    !         n_vertices      = CEILING( REAL(272965,dp) / REAL(n_skip,dp))
    !       END IF
    !     END DO

    !     IF ((.NOT. recognised_file) .AND. par%master) THEN
    !       WRITE(0,*) ''
    !       WRITE(0,*) ' ===== '
    !       WRITE(0,*) 'initialise_basins - WARNING: for Greenland, we expect the file "grndrainagesystems_ekholm.txt"'
    !       WRITE(0,*) '                             This can be downloaded from https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems'
    !       WRITE(0,*) '                             If another file is used, make sure that it has 7 header lines, or alternatively just change the code in initialise_basins'
    !       WRITE(0,*) ' ===== '
    !       WRITE(0,*) ''
    !     END IF

    !     ! Allocate shared memory
    !     CALL allocate_shared_dp_1D(  n_vertices, Vlat,   wVlat  )
    !     CALL allocate_shared_dp_1D(  n_vertices, Vlon,   wVlon  )
    !     CALL allocate_shared_dp_1D(  n_vertices, Vx,     wVx    )
    !     CALL allocate_shared_dp_1D(  n_vertices, Vy,     wVy    )
    !     CALL allocate_shared_int_1D( n_vertices, VID,    wVID   )

    !   ! ===== Read the file =====
    !   ! =========================

    !     IF (par%master) THEN

    !       ! Open the file
    !       OPEN( UNIT = 1337, FILE=filename_basins, ACTION='READ')

    !       ! Skip the header lines
    !       DO i = 1, n_header_lines
    !         READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
    !       END DO

    !       ! Read the actual data
    !       DO vi = 1, n_vertices
    !         READ( UNIT = 1337, FMT=*, IOSTAT=ios) VID_dp, Vlat( vi), Vlon( vi)
    !         IF (ios /= 0) THEN
    !           WRITE(0,*) ' initialise_basins - ERROR: length of text file "', TRIM( filename_basins), '" does not match n_vertices = ', n_vertices
    !           CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    !         END IF

    !         DO vj = 1, n_skip-1
    !           READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
    !         END DO

    !         ! Convert basin ID from floating point to integer
    !         IF     ( viD_dp == 1.1_dp) THEN
    !           VID( vi) = 1
    !         ELSEIF ( viD_dp == 1.2_dp) THEN
    !           VID( vi) = 2
    !         ELSEIF ( viD_dp == 1.3_dp) THEN
    !           VID( vi) = 3
    !         ELSEIF ( viD_dp == 1.4_dp) THEN
    !           VID( vi) = 4
    !         ELSEIF ( viD_dp == 2.1_dp) THEN
    !           VID( vi) = 5
    !         ELSEIF ( viD_dp == 2.2_dp) THEN
    !           VID( vi) = 6
    !         ELSEIF ( viD_dp == 3.1_dp) THEN
    !           VID( vi) = 7
    !         ELSEIF ( viD_dp == 3.2_dp) THEN
    !           VID( vi) = 8
    !         ELSEIF ( viD_dp == 3.3_dp) THEN
    !           VID( vi) = 9
    !         ELSEIF ( viD_dp == 4.1_dp) THEN
    !           VID( vi) = 10
    !         ELSEIF ( viD_dp == 4.2_dp) THEN
    !           VID( vi) = 11
    !         ELSEIF ( viD_dp == 4.3_dp) THEN
    !           VID( vi) = 12
    !         ELSEIF ( viD_dp == 5.0_dp) THEN
    !           VID( vi) = 13
    !         ELSEIF ( viD_dp == 6.1_dp) THEN
    !           VID( vi) = 14
    !         ELSEIF ( viD_dp == 6.2_dp) THEN
    !           VID( vi) = 15
    !         ELSEIF ( viD_dp == 7.1_dp) THEN
    !           VID( vi) = 16
    !         ELSEIF ( viD_dp == 7.2_dp) THEN
    !           VID( vi) = 17
    !         ELSEIF ( viD_dp == 8.1_dp) THEN
    !           VID( vi) = 18
    !         ELSEIF ( viD_dp == 8.2_dp) THEN
    !           VID( vi) = 19
    !         ELSE
    !           WRITE(0,*) 'initialise_basins - ERROR: unrecognised floating-point basin ID in file "', TRIM( filename_basins), '"!'
    !           CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    !         END IF

    !       END DO

    !       ! Close the file
    !       CLOSE( UNIT = 1337)

    !     END IF ! IF (par%master) THEN
    !     CALL sync

    !     ! Project [lat,lon] to [x,y]
    !     CALL partition_list( n_vertices, par%i, par%n, vi1, vi2)
    !     DO vi = vi1, vi2
    !       CALL oblique_sg_projection( Vlon( vi), Vlat( vi), grid%lambda_M, grid%phi_M, grid%alpha_stereo, Vx( vi), Vy( vi))
    !     END DO

    !   END IF ! IF (region_name == 'ANT') THEN

    ! ! ===== Fill in the basins =====
    ! ! ==============================

    !   ! Allocate shared memory
    !   CALL allocate_shared_int_2D( grid%ny, grid%nx, basin_ID_loc, wbasin_ID_loc)

    !   ! Determine number of basins
    !   IF (par%master) nbasins = MAXVAL( VID)
    !   CALL sync

    !   DO bi = 1, nbasins

    !     ! Find range of vertices [vi1,vi2] for this basin
    !     vi1 = 1
    !     DO WHILE (VID( vi1) /= bi)
    !       vi1 = vi1 + 1
    !     END DO
    !     vi2 = vi1
    !     DO WHILE (VID( vi2) == bi .AND. vi2 < n_vertices)
    !       vi2 = vi2 + 1
    !     END DO
    !     vi2 = vi2 - 1

    !     ! Copy these vertices to a single array
    !     CALL allocate_shared_dp_2D( vi2+1-vi1, 2, poly_bi, wpoly_bi)
    !     IF (par%master) THEN
    !       poly_bi( :,1) = Vx( vi1:vi2)
    !       poly_bi( :,2) = Vy( vi1:vi2)
    !     END IF
    !     CALL sync

    !     ! Determine maximum polygon extent, for quick checks
    !     xmin = MINVAL( poly_bi(:,1))
    !     xmax = MAXVAL( poly_bi(:,1))
    !     ymin = MINVAL( poly_bi(:,2))
    !     ymax = MAXVAL( poly_bi(:,2))

    !     ! Check which grid cells lie inside the polygon spanned by these vertices
    !     DO i = grid%i1, grid%i2
    !     DO j = 1, grid%ny
    !       p = [grid%x( i), grid%y( j)]

    !       ! Quick test
    !       IF (p(1) < xmin .OR. p(1) > xmax .OR. &
    !           p(2) < ymin .OR. p(2) > ymax) THEN
    !         ! p cannot lie in the polygon, don't bother checking
    !       ElSE
    !         IF (is_in_polygon( poly_bi, p)) basin_ID_loc( j,i) = bi
    !       END IF

    !     END DO
    !     END DO
    !     CALL sync

    !     ! Clean up this basin's polygon
    !     CALL deallocate_shared( wpoly_bi)

    !   END DO ! DO bi = 1, nbasins

    ! ! ===== Extend basins into the ocean =====
    ! ! ========================================

    !   ! Allocate shared memory
    !   CALL allocate_shared_int_2D( grid%ny, grid%nx, basin_ID_ext_loc, wbasin_ID_ext_loc)

    !   ! Copy data
    !   basin_ID_ext_loc( :,grid%i1:grid%i2) = basin_ID_loc( :,grid%i1:grid%i2)
    !   CALL sync

    !   ! Compile list of ice-front pixels and their basin ID's
    !   IF (par%master) THEN
    !     n_front = 0
    !     DO i = 2, grid%nx-1
    !     DO j = 2, grid%ny-1
    !       IF (basin_ID_loc( j,i) > 0) THEN
    !         IF (MINVAL( basin_ID_loc( j-1:j+1,i-1:i+1)) == 0) THEN
    !           n_front = n_front + 1
    !         END IF
    !       END IF
    !     END DO
    !     END DO
    !   END IF
    !   CALL MPI_BCAST( n_front, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    !   CALL allocate_shared_int_2D( n_front, 2, ij_front, wij_front)

    !   IF (par%master) THEN
    !     k = 0
    !     DO i = 2, grid%nx-1
    !     DO j = 2, grid%ny-1
    !       IF (basin_ID_loc( j,i) > 0) THEN
    !         IF (MINVAL( basin_ID_loc( j-1:j+1,i-1:i+1)) == 0) THEN
    !           k = k + 1
    !           ij_front( k,:) = [i,j]
    !         END IF
    !       END IF
    !     END DO
    !     END DO
    !   END IF
    !   CALL sync

    !   ! For all non-assigned grid cells, find the nearest front cell and copy that basin ID
    !   DO i = grid%i1, grid%i2
    !   DO j = 1, grid%ny

    !     IF (basin_ID_ext_loc( j,i) == 0) THEN

    !       ! Go over all front cells, find the nearest
    !       dist_min  = REAL(MAX(grid%nx,grid%ny),dp) * grid%dx
    !       i_nearest = 0
    !       j_nearest = 0
    !       DO k = 1, n_front
    !         ii = ij_front( k,1)
    !         jj = ij_front( k,2)
    !         dist = SQRT( REAL(ii-i,dp)**2 + REAL(jj-j,dp)**2) * grid%dx
    !         IF (dist < dist_min) THEN
    !           dist_min = dist
    !           i_nearest = ii
    !           j_nearest = jj
    !         END IF
    !       END DO ! DO k = 1, n_front

    !       ! Assign basin ID of nearest front cell
    !       basin_ID_ext_loc( j,i) = basin_ID_loc( j_nearest, i_nearest)

    !     END IF

    !   END DO
    !   END DO
    !   CALL sync

    ! ! ===== Copy final result to the ice structure =====
    ! ! ==================================================

    !   basin_ID( :,grid%i1:grid%i2) = basin_ID_ext_loc( :,grid%i1:grid%i2)
    !   CALL sync

    !   ! Clean up after yourself
    !   CALL deallocate_shared( wVlat            )
    !   CALL deallocate_shared( wVlon            )
    !   CALL deallocate_shared( wVx              )
    !   CALL deallocate_shared( wVy              )
    !   CALL deallocate_shared( wVID             )
    !   CALL deallocate_shared( wbasin_ID_loc    )
    !   CALL deallocate_shared( wbasin_ID_ext_loc)
    !   CALL deallocate_shared( wij_front        )

    ! ! ==== If so specified, merge certain ice basins =====
    ! ! ====================================================

    !   IF     (region_name == 'ANT' .AND. C%do_merge_basins_ANT) THEN
    !     CALL merge_basins_ANT( grid, basin_ID, nbasins)
    !   ELSEIF (region_name == 'GRL' .AND. C%do_merge_basins_GRL) THEN
    !     CALL merge_basins_GRL( grid, basin_ID, nbasins)
    !   END IF

    ELSE
      IF (par%master) WRITE(0,*) 'initialise_basins - ERROR: unknown choice_basin_scheme "', TRIM(choice_basin_scheme), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

  END SUBROUTINE initialise_basins

END MODULE general_ice_model_data_module
