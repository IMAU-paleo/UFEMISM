MODULE general_ice_model_data_module
 !
 ! Only "secondary geometry" right now: masks, surface elevation, TAF

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
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
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             is_floating, thickness_above_floatation, surface_elevation, &
                                             oblique_sg_projection, is_in_polygon
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE data_types_module,               ONLY: type_model_region, type_mesh, type_ice_model, type_grid
  USE mesh_help_functions_module,      ONLY: find_triangle_area
  USE mesh_operators_module,           ONLY: map_a_to_b_2D

  IMPLICIT NONE

CONTAINS

! ===== General ice model data =====
! ==================================

  SUBROUTINE update_general_ice_model_data( mesh, ice)
    ! Update masks, surface elevation, and thickness above floatation

    USE parameters_module, ONLY: ice_density, seawater_density

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
      IF (ice%mask_land_a( vi) == 1) THEN
        ice%dHs_dt_a( vi) = ice%dHb_dt_a( vi) + ice%dHi_dt_a( vi)
      ELSE
        ice%dHs_dt_a( vi) = ice%dHi_dt_a( vi) * (1._dp - ice_density / seawater_density)
      END IF
    END DO
    CALL sync

    ! Determine masks
    CALL determine_masks( mesh, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_general_ice_model_data

! ===== Masks =====
! =================

  SUBROUTINE determine_masks( mesh, ice)
    ! Determine the different masks

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_masks'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL determine_masks_landocean(   mesh, ice)
    CALL determine_masks_ice(         mesh, ice)
    CALL determine_masks_transitions( mesh, ice)
    CALL determine_masks_total(       mesh, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_masks

  SUBROUTINE determine_masks_landocean( mesh, ice)
    ! Determine the different masks
    !
    ! Determine the land/ocean masks (where "land" can also include ice-covered land,
    ! and "ocean" can also include shelf-covered ocean, but not marine grounded ice)
    !
    ! If the C%do_ocean_floodfill option is .TRUE., then lakes (i.e. floating ice on
    ! "no-ocean" points) are also identified as part of the flood fill

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_masks_landocean'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise: start with land everywhere.
    ice%mask_land_a(  mesh%vi1:mesh%vi2) = 1
    ice%mask_ocean_a( mesh%vi1:mesh%vi2) = 0
    ice%mask_lake_a(  mesh%vi1:mesh%vi2) = 0
    CALL sync

    ! Determine land/ocean masks
    IF (C%do_ocean_floodfill) THEN

      ! Use a flood-fill algorith to identify connected ocean
      CALL ocean_floodfill( mesh, ice)

      ! Now that we know the ocean, mask the lakes
      DO vi = mesh%vi1, mesh%vi2
        IF (is_floating( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi)) .AND. ice%mask_land_a( vi) == 1) THEN
          ! This element floats and it's on land, thus a lake
          ice%mask_lake_a( vi) = 1
        END IF
      END DO
      CALL sync

    ELSE

      DO vi = mesh%vi1, mesh%vi2
        IF (is_floating( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))) THEN
          ice%mask_land_a(  vi) = 0
          ice%mask_ocean_a( vi) = 1
        END IF
      END DO
      CALL sync

    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_masks_landocean

  SUBROUTINE determine_masks_ice( mesh, ice)
    ! Determine the different masks
    !
    ! Determine the ice/sheet/shelf masks

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_masks_ice'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise: no ice everywhere
    ice%mask_ice_a(   mesh%vi1:mesh%vi2) = 0
    ice%mask_sheet_a( mesh%vi1:mesh%vi2) = 0
    ice%mask_shelf_a( mesh%vi1:mesh%vi2) = 0
    CALL sync

    ! Determine ice/sheet/shelf masks
    DO vi = mesh%vi1, mesh%vi2

      ! Ice
      IF (ice%Hi_a( vi) > TINY(ice%Hi_a( 1))) THEN
        ice%mask_ice_a(  vi) = 1
      END IF

      ! Sheet
      IF (ice%mask_ice_a( vi) == 1 .AND. ice%mask_land_a( vi) == 1) THEN
        ice%mask_sheet_a( vi) = 1
      END IF

      ! Shelf
      IF (ice%mask_ice_a( vi) == 1 .AND. ice%mask_ocean_a( vi) == 1) THEN
        ice%mask_shelf_a( vi) = 1
      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_masks_ice

  SUBROUTINE determine_masks_transitions( mesh, ice)
    ! Determine the different masks
    !
    ! Determine the "transitional" masks (coast, ice margin, grounding line, calving front)

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_masks_transitions'
    INTEGER                                            :: vi, ci, vc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    ice%mask_coast_a(  mesh%vi1:mesh%vi2) = 0
    ice%mask_margin_a( mesh%vi1:mesh%vi2) = 0
    ice%mask_gl_a(     mesh%vi1:mesh%vi2) = 0
    ice%mask_cf_a(     mesh%vi1:mesh%vi2) = 0
    CALL sync

    ! Determine coast, grounding line and calving front
    DO vi = mesh%vi1, mesh%vi2

      IF (ice%mask_land_a( vi) == 1 .AND. ice%mask_ice_a( vi) == 0) THEN
        ! Ice-free land bordering ocean equals coastline
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ocean_a( vc) == 1) THEN
            ice%mask_coast_a( vi) =  1
          END IF
        END DO
      END IF

      IF (ice%mask_ice_a( vi) == 1) THEN
        ! Ice bordering non-ice equals margin
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ice_a( vc) == 0) THEN
            ice%mask_margin_a( vi) =  1
          END IF
        END DO
      END IF

      IF (ice%mask_sheet_a( vi) == 1) THEN
        ! Sheet bordering shelf equals grounding line
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_shelf_a( vc) == 1) THEN
            ice%mask_gl_a( vi) =  1
          END IF
        END DO
      END IF

      IF (ice%mask_shelf_a( vi) == 1) THEN
        ! Shelf bordering sheet equals floating side of grounding line
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_sheet_a( vc) == 1) THEN
            ice%mask_glf_a( vi) =  1
          END IF
        END DO
      END IF

      IF (ice%mask_ice_a( vi) == 1) THEN
        ! Ice (sheet or shelf) bordering open ocean equals calving front
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ocean_a( vc) == 1 .AND. ice%mask_ice_a( vc) == 0) THEN
            ice%mask_cf_a( vi) =  1
          END IF
        END DO
      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_masks_transitions

  SUBROUTINE determine_masks_total( mesh, ice)
    ! Determine the different masks
    !
    ! Determine the "total" mask (used only for writing to output!)

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_masks_total'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Land/ocean
      IF   (ice%mask_land_a( vi) == 1) THEN
        ice%mask_a( vi) = C%type_land
      ELSE
        ice%mask_a( vi) = C%type_ocean
      END IF

      ! Coast
      IF (ice%mask_coast_a( vi) == 1) THEN
        ice%mask_a( vi) = C%type_coast
      END IF

      ! Sheet/shelf
      IF     (ice%mask_sheet_a( vi) == 1) THEN
        ice%mask_a( vi) = C%type_sheet
      ELSEIF (ice%mask_shelf_a( vi) == 1) THEN
        ice%mask_a( vi) = C%type_shelf
      END IF

      ! Lakes

      IF (ice%mask_lake_a( vi) == 1) THEN
        ice%mask_a( vi) = C%type_lake
      END IF

      ! Ice margin / grounding line / calving front
      IF (ice%mask_margin_a( vi) == 1) THEN
        ice%mask_a( vi) = C%type_margin
      END IF
      IF (ice%mask_gl_a( vi) == 1) THEN
        ice%mask_a( vi) = C%type_groundingline
      END IF
      IF (ice%mask_cf_a( vi) == 1) THEN
        ice%mask_a( vi) = C%type_calvingfront
      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_masks_total

! ===== Ocean flood-fill =====
! ============================

  SUBROUTINE ocean_floodfill( mesh, ice)
    ! Use a simple floodfill algorithm to determine the ocean mask,
    ! to prevent the formation of (pro-/sub-glacial) lakes

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'ocean_floodfill'
    INTEGER                                            :: vi, ci, vc
    INTEGER, DIMENSION(:    ), POINTER                 :: map
    INTEGER, DIMENSION(:    ), POINTER                 :: stack
    INTEGER                                            :: wmap, wstack
    INTEGER                                            :: stackN

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory for the map, because even though the flood-fill is done only
    ! by the Master, the other processes need access to the resulting filled map.
    CALL allocate_shared_int_1D( mesh%nV,   map,   wmap)
    CALL allocate_shared_int_1D( mesh%nV, stack, wstack)

    ! No easy way to parallelise flood-fill, just let the Master do it
    IF (par%master) THEN

      ice%mask_land_a  = 1
      ice%mask_ocean_a = 0
      ice%mask_a       = C%type_land

      map    = 0
      stack  = 0
      stackN = 0

      ! Let the ocean flow in from the domain edges
      DO vi = 1, mesh%nV
        IF ( .NOT. mesh%edge_index( vi) == 0) THEN
          stackN = stackN + 1
          stack( stackN) = vi
          map( vi) = 1
        END IF
      END DO

      ! Flood fill
      DO WHILE (stackN > 0)

        ! Inspect the last element of the stack
        vi = stack( stackN)

        ! Remove it from the stack
        stack( stackN) = 0
        stackN = stackN-1

        IF (is_floating( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))) THEN
          ! This element is ocean

          ! Mark it as such on the map
          map( vi) = 2
          ice%mask_ocean_a( vi) = 1
          ice%mask_land_a(  vi) = 0

          ! Add its neighbours to the stack
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            IF (map( vc)==0) THEN
              ! This neighbour isn't yet mapped or stacked
              map( vc) = 1
              stackN = stackN + 1
              stack( stackN) = vc
            END IF
          END DO

        ELSE
          ! This element is land
        END IF

      END DO ! WHILE (stackN > 0)

    END IF ! (par%master)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wmap)
    CALL deallocate_shared( wstack)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE ocean_floodfill

! ===== Sub-grid grounded fractions =====
! =======================================

  SUBROUTINE determine_grounded_fractions( mesh, ice)
    ! Determine the grounded fractions of all grid cells

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_grounded_fractions'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL determine_grounded_fractions_a( mesh, ice)
    CALL determine_grounded_fractions_b( mesh, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_grounded_fractions

  SUBROUTINE determine_grounded_fractions_a( mesh, ice)
    ! Determine the grounded fractions of all grid cells on the a-grid

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_grounded_fractions_a'
    REAL(dp), DIMENSION(:    ), POINTER                ::  TAF_b
    INTEGER                                            :: wTAF_b
    INTEGER                                            :: vi, ci, vj, iti, iti2, ti1, ti2
    REAL(dp)                                           :: TAF_max, TAF_min
    REAL(dp), DIMENSION(2)                             :: va, ccb1, ccb2
    REAL(dp)                                           :: TAFa, TAFb, TAFc, A_vor, A_tri_tot, A_tri_grnd, A_grnd

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Map thickness-above-floatation to the b-grid
    CALL allocate_shared_dp_1D( mesh%nTri, TAF_b, wTAF_b)
    CALL map_a_to_b_2D(  mesh, ice%TAF_a, TAF_b)

    DO vi = mesh%vi1, mesh%vi2

      ! Skip border vertices
      IF (mesh%edge_index( vi) > 0) THEN
        ice%f_grnd_a( vi) = 0._dp
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
        ice%f_grnd_a( vi) = 1._dp
        CYCLE
      END IF

      ! If the entire local neighbourhood is floating, the answer is trivial
      IF (TAF_max <= 0._dp) THEN
        ice%f_grnd_a( vi) = 0._dp
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
        CALL determine_grounded_area_triangle( va, ccb1, ccb2, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

        A_vor  = A_vor  + A_tri_tot
        A_grnd = A_grnd + A_tri_grnd

      END DO ! DO iati = 1, mesh%niTriAaAc( avi)

      ! Calculate the grounded fraction of this Voronoi cell
      ice%f_grnd_a( vi) = A_grnd / A_vor

    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wTAF_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_grounded_fractions_a

  SUBROUTINE determine_grounded_fractions_b( mesh, ice)
    ! Determine the grounded fractions of all grid cells on the b-grid

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_grounded_fractions_b'
    INTEGER                                            :: ti, via, vib, vic
    REAL(dp)                                           :: TAF_max, TAF_min
    REAL(dp), DIMENSION(2)                             :: va, vb, vc
    REAL(dp)                                           :: TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd

    ! Add routine to path
    CALL init_routine( routine_name)

    DO ti = mesh%ti1, mesh%ti2

      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = MAXVAL([ ice%TAF_a( via), ice%TAF_a( vib), ice%TAF_a( vic)])
      TAF_min = MINVAL([ ice%TAF_a( via), ice%TAF_a( vib), ice%TAF_a( vic)])

      ! If the entire local neighbourhood is grounded, the answer is trivial
      IF (TAF_min >= 0._dp) THEN
        ice%f_grnd_b( ti) = 1._dp
        CYCLE
      END IF

      ! If the entire local neighbourhood is floating, the answer is trivial
      IF (TAF_max <= 0._dp) THEN
        ice%f_grnd_b( ti) = 0._dp
        CYCLE
      END IF

      ! The local neighbourhood contains both grounded and floating vertices.

      va   = mesh%V( via,:)
      vb   = mesh%V( vib,:)
      vc   = mesh%V( vic,:)

      TAFa = ice%TAF_a( via)
      TAFb = ice%TAF_a( vib)
      TAFc = ice%TAF_a( vic)

      ! Determine total area of, and grounded area within, this subtriangle
      CALL determine_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

      ! Calculate the grounded fraction of this Voronoi cell
      ice%f_grnd_b( ti) = A_tri_grnd / A_tri_tot

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_grounded_fractions_b

  SUBROUTINE determine_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
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
        CALL determine_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
      ELSEIF (TAFa <= 0._dp .AND. TAFb >= 0._dp .AND. TAFc <= 0._dp) THEN
        ! b is grounded, a and c are floating
        CALL determine_grounded_area_triangle_1grnd_2flt( vb, vc, va, TAFb, TAFc, TAFa, A_tri_grnd)
      ELSEIF (TAFa <= 0._dp .AND. TAFb <= 0._dp .AND. TAFc >= 0._dp) THEN
        ! c is grounded, a and b are floating
        CALL determine_grounded_area_triangle_1grnd_2flt( vc, va, vb, TAFc, TAFa, TAFb, A_tri_grnd)
      ELSEIF (TAFa <= 0._dp .AND. TAFb >= 0._dp .AND. TAFc >= 0._dp) THEN
        ! a is floating, b and c are grounded
        CALL determine_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      ELSEIF (TAFa >= 0._dp .AND. TAFb <= 0._dp .AND. TAFc >= 0._dp) THEN
        ! b is floating, c and a are grounded
        CALL determine_grounded_area_triangle_1flt_2grnd( vb, vc, va, TAFb, TAFc, TAFa, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      ELSEIF (TAFa >= 0._dp .AND. TAFb >= 0._dp .AND. TAFc <= 0._dp) THEN
        ! c is floating, a and b are grounded
        CALL determine_grounded_area_triangle_1flt_2grnd( vc, va, vb, TAFc, TAFa, TAFb, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      ELSE
        A_tri_grnd = 0._dp
        CALL crash('TAF = [{dp_01},{dp_02},{dp_03}]', dp_01 = TAFa, dp_02 = TAFb, dp_03 = TAFc)
      END IF

    END IF

  END SUBROUTINE determine_grounded_area_triangle

  SUBROUTINE determine_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
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

  END SUBROUTINE determine_grounded_area_triangle_1grnd_2flt

  SUBROUTINE determine_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_tri_flt)
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

  END SUBROUTINE determine_grounded_area_triangle_1flt_2grnd

! ===== Sub-grid ice-filled fraction and effective ice thickness at the calving front =====
! =========================================================================================

  SUBROUTINE determine_floating_margin_fraction( mesh, ice)
    ! Determine the ice-filled fraction and effective ice thickness of floating margin pixels

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_floating_margin_fraction'
    INTEGER                                            :: vi, ci, vc
    LOGICAL                                            :: has_noncf_neighbours
    REAL(dp)                                           :: Hi_neighbour_max

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise
      IF (ice%mask_ice_a( vi) == 1) THEN
        ice%float_margin_frac_a( vi) = 1._dp
        ice%Hi_eff_cf_a(         vi) = ice%Hi_a( vi)
      ELSE
        ice%float_margin_frac_a( vi) = 0._dp
        ice%Hi_eff_cf_a(         vi) = 0._dp
      END IF

      IF (ice%mask_cf_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 1) THEN

        ! First check if any non-calving-front neighbours actually exist
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ice_a( vc) == 1 .AND. ice%mask_cf_a( vc) == 0) THEN
            has_noncf_neighbours = .TRUE.
          END IF
        END DO

        ! If not, then the floating fraction is defined as 1
        IF (.NOT. has_noncf_neighbours) THEN
          ice%float_margin_frac_a( vi) = 1._dp
          ice%Hi_eff_cf_a(         vi) = ice%Hi_a( vi)
          CYCLE
        END IF

        ! If so, find the ice thickness the thickest non-calving-front neighbour
        Hi_neighbour_max = 0._dp
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ice_a( vc) == 1 .AND. ice%mask_cf_a( vc) == 0) THEN
            Hi_neighbour_max = MAX( Hi_neighbour_max, ice%Hi_a( vc))
          END IF
        END DO

        ! If the thickest non-calving-front neighbour has thinner ice, define the fraction as 1
        IF (Hi_neighbour_max < ice%Hi_a( vi)) THEN
          ice%float_margin_frac_a( vi) = 1._dp
          ice%Hi_eff_cf_a(         vi) = ice%Hi_a( vi)
          CYCLE
        END IF

        ! Calculate ice-filled fraction
        ice%float_margin_frac_a( vi) = ice%Hi_a( vi) / Hi_neighbour_max
        ice%Hi_eff_cf_a(         vi) = Hi_neighbour_max

      END IF

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%float_margin_frac_a, 'ice%float_margin_frac_a')
    CALL check_for_NaN_dp_1D( ice%Hi_eff_cf_a        , 'ice%Hi_eff_cf_a'        )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_floating_margin_fraction

! ===== The no-ice mask, to prevent ice growth in certain areas =====
! ===================================================================

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
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_mask_noice_EAS_remove_GRL

  SUBROUTINE initialise_mask_noice_GRL_remove_Ellesmere( mesh, mask_noice)
    ! Prevent ice growth in the Ellesmere Island part of the Greenland domain

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),               INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),    INTENT(OUT)   :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'initialise_mask_noice_GRL_remove_Ellesmere'
    INTEGER                                      :: vi
    REAL(dp), DIMENSION(2)                       :: pa_latlon, pb_latlon
    REAL(dp)                                     :: xa,ya,xb,yb
    REAL(dp), DIMENSION(2)                       :: pa, pb
    REAL(dp)                                     :: yl_ab

    ! Add routine to path
    CALL init_routine( routine_name)

    ! The two endpoints in lat,lon
    pa_latlon = [76.74_dp, -74.79_dp]
    pb_latlon = [82.19_dp, -60.00_dp]

    ! The two endpoints in x,y
    CALL oblique_sg_projection( pa_latlon(2), pa_latlon(1), mesh%lambda_M, mesh%phi_M, mesh%alpha_stereo, xa, ya)
    CALL oblique_sg_projection( pb_latlon(2), pb_latlon(1), mesh%lambda_M, mesh%phi_M, mesh%alpha_stereo, xb, yb)

    pa = [xa,ya]
    pb = [xb,yb]

    DO vi = mesh%vi1, mesh%vi2
      yl_ab = pa(2) + (mesh%V( vi,1) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
      IF (mesh%V( vi,2) > pa(2) .AND. mesh%V( vi,2) > yl_ab .AND. mesh%V( vi,1) < pb(1)) THEN
        mask_noice( vi) = 1
      ELSE
        mask_noice( vi) = 0
      END IF
    END DO
    CALL sync

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
    CALL sync

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
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_mask_noice_MISMIPplus

! ===== Ice drainage basins from an external polygon file (mesh version) =====
! ============================================================================

  SUBROUTINE initialise_basins( mesh, basin_ID, nbasins, region_name)
    ! Define the ice basins mask from an external text file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),          INTENT(INOUT) :: basin_ID
    INTEGER,                             INTENT(INOUT) :: nbasins
    CHARACTER(LEN=3)                                   :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basins'
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
    INTEGER,  DIMENSION(:    ), POINTER                ::  basin_ID_loc,  basin_ID_ext_loc
    INTEGER                                            :: wbasin_ID_loc, wbasin_ID_ext_loc
    INTEGER                                            :: n_front
    INTEGER,  DIMENSION(:    ), POINTER                :: ij_front
    INTEGER                                            :: wij_front
    INTEGER                                            :: ii,jj,k,i_nearest,j_nearest
    REAL(dp)                                           :: dist,dist_min

    ! Add routine to path
    CALL init_routine( routine_name)

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

      IF (par%master) WRITE(0,*) '  Reading basins for ', TRIM(region_name), ' from file "', TRIM(filename_basins), '"...'

      IF (region_name == 'ANT') THEN
        ! Antarctica: ant_full_drainagesystem_polygons.txt
        ! Can be downloaded from: https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems
        ! A text file with 7 header lines, followed by three columns of data:
        !    Lat, Lon, basin ID

      ! ===== Check if this is really the file we're reading =====
      ! ==========================================================

        recognised_file = .FALSE.
        n_header_lines  = 0
        n_vertices      = 0
        n_skip          = 0

        DO i = 1, 256-36
          IF (filename_basins(i:i+35) == 'ant_full_drainagesystem_polygons.txt') THEN
            recognised_file = .TRUE.
            n_header_lines  = 7
            n_skip          = 5 ! Since the polygon file is at a ridiculously high resolution, downscale it a bit for efficiency
            n_vertices      = CEILING( REAL(901322,dp) / REAL(n_skip,dp))
          END IF
        END DO

        IF ((.NOT. recognised_file) .AND. par%master) THEN
          CALL warning('for Antarctica, we expect the file "ant_full_drainagesystem_polygons.txt". ' // &
                       'This can be downloaded from https: // earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems. ' // &
                       'If another file is used, make sure that it has 7 header lines, or alternatively just change the code in initialise_basins.')
        END IF

        ! Allocate shared memory
        CALL allocate_shared_dp_1D(  n_vertices, Vlat, wVlat)
        CALL allocate_shared_dp_1D(  n_vertices, Vlon, wVlon)
        CALL allocate_shared_dp_1D(  n_vertices, Vx,   wVx  )
        CALL allocate_shared_dp_1D(  n_vertices, Vy,   wVy  )
        CALL allocate_shared_int_1D( n_vertices, VID,  wVID )

      ! ===== Read the file =====
      ! =========================

        IF (par%master) THEN

          ! Open the file
          OPEN( UNIT = 1337, FILE=filename_basins, ACTION='READ')

          ! Skip the header lines
          DO i = 1, n_header_lines
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
          END DO

          ! Read the actual data
          DO vi = 1, n_vertices
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) Vlat( vi), Vlon( vi), VID( vi)
            IF (ios /= 0) THEN
              CALL crash('length of text file "' // TRIM( filename_basins) // '" does not match n_vertices = {int_01}!', int_01 = n_vertices)
            END IF

            DO vj = 1, n_skip-1
              READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
            END DO
          END DO

          ! Close the file
          CLOSE( UNIT = 1337)

        END IF ! IF (par%master) THEN
        CALL sync

        ! Project [lat,lon] to [x,y]
        CALL partition_list( n_vertices, par%i, par%n, vi1, vi2)
        DO vi = vi1, vi2
          CALL oblique_sg_projection( Vlon( vi), Vlat( vi), mesh%lambda_M, mesh%phi_M, mesh%alpha_stereo, Vx( vi), Vy( vi))
        END DO

      ELSEIF (region_name == 'GRL') THEN
        ! Greenland: grndrainagesystems_ekholm.txt
        ! Can be downloaded from: https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems
        ! A text file with 7 header lines, followed by three columns of data:
        !    basin ID, Lat, Lon   (NOTE: here basin ID is a floating-point rather than an integer!)

      ! ===== Check if this is really the file we're reading =====
      ! ==========================================================

        recognised_file = .FALSE.
        n_header_lines  = 0
        n_vertices      = 0
        n_skip          = 0

        DO i = 1, 256-29
          IF (filename_basins(i:i+28) == 'grndrainagesystems_ekholm.txt') THEN
            recognised_file = .TRUE.
            n_header_lines  = 7
            n_skip          = 5 ! Since the polygon file is at a ridiculously high resolution, downscale it a bit for efficiency
            n_vertices      = CEILING( REAL(272965,dp) / REAL(n_skip,dp))
          END IF
        END DO

        IF ((.NOT. recognised_file) .AND. par%master) THEN
          CALL warning('for Greenland, we expect the file "grndrainagesystems_ekholm.txt". ' // &
                       'This can be downloaded from https: // earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems. ' // &
                       'If another file is used, make sure that it has 7 header lines, or alternatively just change the code in initialise_basins.')
        END IF

        ! Allocate shared memory
        CALL allocate_shared_dp_1D(  n_vertices, Vlat,   wVlat  )
        CALL allocate_shared_dp_1D(  n_vertices, Vlon,   wVlon  )
        CALL allocate_shared_dp_1D(  n_vertices, Vx,     wVx    )
        CALL allocate_shared_dp_1D(  n_vertices, Vy,     wVy    )
        CALL allocate_shared_int_1D( n_vertices, VID,    wVID   )

      ! ===== Read the file =====
      ! =========================

        IF (par%master) THEN

          ! Open the file
          OPEN( UNIT = 1337, FILE=filename_basins, ACTION='READ')

          ! Skip the header lines
          DO i = 1, n_header_lines
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
          END DO

          ! Read the actual data
          DO vi = 1, n_vertices
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) VID_dp, Vlat( vi), Vlon( vi)
            IF (ios /= 0) THEN
              CALL crash('length of text file "' // TRIM( filename_basins) // '" does not match n_vertices = {int_01}!', int_01 = n_vertices)
            END IF

            DO vj = 1, n_skip-1
              READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
            END DO

            ! Convert basin ID from floating point to integer
            IF     ( viD_dp == 1.1_dp) THEN
              VID( vi) = 1
            ELSEIF ( viD_dp == 1.2_dp) THEN
              VID( vi) = 2
            ELSEIF ( viD_dp == 1.3_dp) THEN
              VID( vi) = 3
            ELSEIF ( viD_dp == 1.4_dp) THEN
              VID( vi) = 4
            ELSEIF ( viD_dp == 2.1_dp) THEN
              VID( vi) = 5
            ELSEIF ( viD_dp == 2.2_dp) THEN
              VID( vi) = 6
            ELSEIF ( viD_dp == 3.1_dp) THEN
              VID( vi) = 7
            ELSEIF ( viD_dp == 3.2_dp) THEN
              VID( vi) = 8
            ELSEIF ( viD_dp == 3.3_dp) THEN
              VID( vi) = 9
            ELSEIF ( viD_dp == 4.1_dp) THEN
              VID( vi) = 10
            ELSEIF ( viD_dp == 4.2_dp) THEN
              VID( vi) = 11
            ELSEIF ( viD_dp == 4.3_dp) THEN
              VID( vi) = 12
            ELSEIF ( viD_dp == 5.0_dp) THEN
              VID( vi) = 13
            ELSEIF ( viD_dp == 6.1_dp) THEN
              VID( vi) = 14
            ELSEIF ( viD_dp == 6.2_dp) THEN
              VID( vi) = 15
            ELSEIF ( viD_dp == 7.1_dp) THEN
              VID( vi) = 16
            ELSEIF ( viD_dp == 7.2_dp) THEN
              VID( vi) = 17
            ELSEIF ( viD_dp == 8.1_dp) THEN
              VID( vi) = 18
            ELSEIF ( viD_dp == 8.2_dp) THEN
              VID( vi) = 19
            ELSE
              WRITE(0,*) 'initialise_basins - ERROR: unrecognised floating-point basin ID in file "', TRIM( filename_basins), '"!'
              CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
            END IF

          END DO

          ! Close the file
          CLOSE( UNIT = 1337)

        END IF ! IF (par%master) THEN
        CALL sync

        ! Project [lat,lon] to [x,y]
        CALL partition_list( n_vertices, par%i, par%n, vi1, vi2)
        DO vi = vi1, vi2
          CALL oblique_sg_projection( Vlon( vi), Vlat( vi), mesh%lambda_M, mesh%phi_M, mesh%alpha_stereo, Vx( vi), Vy( vi))
        END DO

      END IF ! IF (region_name == 'ANT') THEN

    ! ===== Fill in the basins =====
    ! ==============================

      ! Allocate shared memory
      CALL allocate_shared_int_1D( mesh%nV, basin_ID_loc, wbasin_ID_loc)

      ! Determine number of basins
      IF (par%master) nbasins = MAXVAL( VID)
      CALL sync

      DO bi = 1, nbasins

        ! Find range of vertices [vi1,vi2] for this basin
        vi1 = 1
        DO WHILE (VID( vi1) /= bi)
          vi1 = vi1 + 1
        END DO
        vi2 = vi1
        DO WHILE (VID( vi2) == bi .AND. vi2 < n_vertices)
          vi2 = vi2 + 1
        END DO
        vi2 = vi2 - 1

        ! Copy these vertices to a single array
        CALL allocate_shared_dp_2D( vi2+1-vi1, 2, poly_bi, wpoly_bi)
        IF (par%master) THEN
          poly_bi( :,1) = Vx( vi1:vi2)
          poly_bi( :,2) = Vy( vi1:vi2)
        END IF
        CALL sync

        ! Determine maximum polygon extent, for quick checks
        xmin = MINVAL( poly_bi(:,1))
        xmax = MAXVAL( poly_bi(:,1))
        ymin = MINVAL( poly_bi(:,2))
        ymax = MAXVAL( poly_bi(:,2))

        ! Check which grid cells lie inside the polygon spanned by these vertices
        DO vi = mesh%vi1, mesh%vi2
          p = [mesh%V(vi,1), mesh%V(vi,2)]

          ! Quick test
          IF (p(1) < xmin .OR. p(1) > xmax .OR. &
              p(2) < ymin .OR. p(2) > ymax) THEN
            ! p cannot lie in the polygon, don't bother checking
          ElSE
            IF (is_in_polygon( poly_bi, p)) basin_ID_loc( vi) = bi
          END IF

        END DO
        CALL sync

        ! Clean up this basin's polygon
        CALL deallocate_shared( wpoly_bi)

      END DO ! DO bi = 1, nbasins

    ! ===== Extend basins into the ocean =====
    ! ========================================

      ! Allocate shared memory
      CALL allocate_shared_int_1D( mesh%nV, basin_ID_ext_loc, wbasin_ID_ext_loc)

      ! Copy data
      basin_ID_ext_loc( mesh%vi1:mesh%vi2) = basin_ID_loc( mesh%vi1:mesh%vi2)
      CALL sync

      ! Compile list of ice-front pixels and their basin ID's
      IF (par%master) THEN
        n_front = 0
        DO vi = 1, mesh%nV
          IF (basin_ID_loc( vi) > 0) THEN
            ! Check if any connected vertex doesn't have a basin ID yet
            ! (trimming triling 0s in mesh%C)
            IF (MINVAL( basin_ID_loc(mesh%C(vi,1:mesh%nC(vi)))) == 0) THEN
              n_front = n_front + 1
            END IF
          END IF
        END DO
      END IF
      CALL MPI_BCAST( n_front, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      CALL allocate_shared_int_1D( n_front, ij_front, wij_front)

      IF (par%master) THEN
        k = 0
        DO vi = 1, mesh%nV
          IF (basin_ID_loc( vi) > 0) THEN
            IF (MINVAL( basin_ID_loc( mesh%C(vi,1:mesh%nC(vi)))) == 0) THEN
              k = k + 1
              ij_front( k) = vi
            END IF
          END IF
        END DO
      END IF
      CALL sync

      ! For all non-assigned grid cells, find the nearest front cell and copy that basin ID
      DO vi = mesh%vi1, mesh%vi2

        IF (basin_ID_ext_loc( vi) == 0) THEN

          ! Go over all front cells, find the nearest
          dist_min  = REAL(mesh%nV,dp) * C%res_max * 1000.0_dp
          i_nearest = 0
          DO k = 1, n_front
            ii = ij_front( k)
            dist = NORM2( mesh%V( vi,:) - mesh%V(ii,:))
            IF (dist < dist_min) THEN
              dist_min = dist
              i_nearest = ii
            END IF
          END DO ! DO k = 1, n_front

          ! Assign basin ID of nearest front cell
          basin_ID_ext_loc( vi) = basin_ID_loc( i_nearest)

        END IF

      END DO
      CALL sync

    ! ===== Copy final result to the ice structure =====
    ! ==================================================

      basin_ID( mesh%vi1:mesh%vi2) = basin_ID_ext_loc( mesh%vi1:mesh%vi2)
      CALL sync

      ! Clean up after yourself
      CALL deallocate_shared( wVlat            )
      CALL deallocate_shared( wVlon            )
      CALL deallocate_shared( wVx              )
      CALL deallocate_shared( wVy              )
      CALL deallocate_shared( wVID             )
      CALL deallocate_shared( wbasin_ID_loc    )
      CALL deallocate_shared( wbasin_ID_ext_loc)
      CALL deallocate_shared( wij_front        )

    ! ==== If so specified, merge certain ice basins =====
    ! ====================================================

      IF     (region_name == 'ANT' .AND. C%do_merge_basins_ANT) THEN
        CALL merge_basins_ANT( mesh, basin_ID, nbasins)
      ELSEIF (region_name == 'GRL' .AND. C%do_merge_basins_GRL) THEN
        CALL merge_basins_GRL( mesh, basin_ID, nbasins)
      END IF

    ELSE
      CALL crash('unknown choice_basin_scheme "' // TRIM(choice_basin_scheme) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basins

  SUBROUTINE merge_basins_ANT(  mesh, basin_ID, nbasins)
    ! Merge some of the Antarctic basins to match the more common definitions

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),          INTENT(INOUT) :: basin_ID
    INTEGER,                             INTENT(INOUT) :: nbasins

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'merge_basins_ANT'
    INTEGER                                            :: vi,k
    INTEGER,  DIMENSION( 27)                           :: T
    INTEGER,  DIMENSION(:    ), POINTER                ::  basin_ID_old
    INTEGER                                            :: wbasin_ID_old

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_1D( mesh%nV, basin_ID_old, wbasin_ID_old)

    ! Copy data
    basin_ID_old( mesh%vi1:mesh%vi2) = basin_ID( mesh%vi1:mesh%vi2)
    CALL sync

    ! Define the merging table
    T(  1) = 1
    T(  2) = 1
    T(  3) = 1
    T(  4) = 2
    T(  5) = 3
    T(  6) = 4
    T(  7) = 5
    T(  8) = 6
    T(  9) = 6
    T( 10) = 6
    T( 11) = 6
    T( 12) = 7
    T( 13) = 8
    T( 14) = 9
    T( 15) = 10
    T( 16) = 11
    T( 17) = 11
    T( 18) = 11
    T( 19) = 11
    T( 20) = 12
    T( 21) = 13
    T( 22) = 13
    T( 23) = 13
    T( 24) = 14
    T( 25) = 15
    T( 26) = 16
    T( 27) = 17

    ! Perform the merge
    DO vi = mesh%vi1, mesh%vi2
      DO k = 1, 27
        IF (basin_ID_old( vi) == k) basin_ID( vi) = T( k)
      END DO
    END DO
    IF (par%master) nbasins = 17
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wbasin_ID_old)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE merge_basins_ANT

  SUBROUTINE merge_basins_GRL(  mesh, basin_ID, nbasins)
    ! Merge some of the Greenland basins to match the more common definitions

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(:    ),          INTENT(INOUT) :: basin_ID
    INTEGER,                             INTENT(INOUT) :: nbasins

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'merge_basins_GRL'
    INTEGER                                            :: vi,k
    INTEGER,  DIMENSION( 19)                           :: T
    INTEGER,  DIMENSION(:    ), POINTER                ::  basin_ID_old
    INTEGER                                            :: wbasin_ID_old

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_1D( mesh%nV, basin_ID_old, wbasin_ID_old)

    ! Copy data
    basin_ID_old( mesh%vi1:mesh%vi2) = basin_ID( mesh%vi1:mesh%vi2)
    CALL sync

    ! Define the merging table
    T(  1) = 1
    T(  2) = 1
    T(  3) = 1
    T(  4) = 1
    T(  5) = 2
    T(  6) = 2
    T(  7) = 3
    T(  8) = 3
    T(  9) = 3
    T( 10) = 4
    T( 11) = 4
    T( 12) = 4
    T( 13) = 5
    T( 14) = 6
    T( 15) = 6
    T( 16) = 7
    T( 17) = 7
    T( 18) = 8
    T( 19) = 8

    ! Perform the merge
    DO vi = mesh%vi1, mesh%vi2
      DO k = 1, 19
        IF (basin_ID_old( vi) == k) basin_ID( vi) = T( k)
      END DO
    END DO
    IF (par%master) nbasins = 8
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wbasin_ID_old)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE merge_basins_GRL

! ===== Ice drainage basins from an external polygon file (grid version) =====
! ============================================================================

  SUBROUTINE initialise_basins_grid( grid, basin_ID, nbasins, region_name)
    ! Define the ice basins mask from an external text file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(INOUT) :: basin_ID
    INTEGER,                             INTENT(INOUT) :: nbasins
    CHARACTER(LEN=3)                                   :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basins_grid'
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

    ! Add routine to path
    CALL init_routine( routine_name)

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

      basin_ID( grid%i1:grid%i2,:) = 1
      IF (par%master) nbasins = 1
      CALL sync

    ELSEIF (choice_basin_scheme == 'file') THEN
      ! Define basins from an external text file describing the polygons

      IF (par%master) WRITE(0,*) '    Reading basins for ', TRIM(region_name), ' from file "', TRIM(filename_basins), '"...'

      IF (region_name == 'ANT') THEN
        ! Antarctica: ant_full_drainagesystem_polygons.txt
        ! Can be downloaded from: https: // earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems
        ! A text file with 7 header lines, followed by three columns of data:
        !    Lat, Lon, basin ID

      ! ===== Check if this is really the file we're reading =====
      ! ==========================================================

        recognised_file = .FALSE.
        n_header_lines  = 0
        n_vertices      = 0
        n_skip          = 0

        DO i = 1, 256-36
          IF (filename_basins(i:i+35) == 'ant_full_drainagesystem_polygons.txt') THEN
            recognised_file = .TRUE.
            n_header_lines  = 7
            n_skip          = 5 ! Since the polygon file is at a ridiculously high resolution, downscale it a bit for efficiency
            n_vertices      = CEILING( REAL(901322,dp) / REAL(n_skip,dp))
          END IF
        END DO

        IF ((.NOT. recognised_file) .AND. par%master) THEN
          CALL warning('for Antarctica, we expect the file "ant_full_drainagesystem_polygons.txt". ' // &
                       'This can be downloaded from https: // earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems. ' // &
                       'If another file is used, make sure that it has 7 header lines, or alternatively just change the code in initialise_basins.')
        END IF

        ! Allocate shared memory
        CALL allocate_shared_dp_1D(  n_vertices, Vlat, wVlat)
        CALL allocate_shared_dp_1D(  n_vertices, Vlon, wVlon)
        CALL allocate_shared_dp_1D(  n_vertices, Vx,   wVx  )
        CALL allocate_shared_dp_1D(  n_vertices, Vy,   wVy  )
        CALL allocate_shared_int_1D( n_vertices, VID,  wVID )

      ! ===== Read the file =====
      ! =========================

        IF (par%master) THEN

          ! Open the file
          OPEN( UNIT = 1337, FILE=filename_basins, ACTION='READ')

          ! Skip the header lines
          DO i = 1, n_header_lines
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
          END DO

          ! Read the actual data
          DO vi = 1, n_vertices
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) Vlat( vi), Vlon( vi), VID( vi)
            IF (ios /= 0) THEN
              CALL crash('length of text file "' // TRIM( filename_basins) // '" does not match n_vertices = {int_01}!', int_01 = n_vertices)
            END IF

            DO vj = 1, n_skip-1
              READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
            END DO
          END DO

          ! Close the file
          CLOSE( UNIT = 1337)

        END IF ! IF (par%master) THEN
        CALL sync

        ! Project [lat,lon] to [x,y]
        CALL partition_list( n_vertices, par%i, par%n, vi1, vi2)
        DO vi = vi1, vi2
          CALL oblique_sg_projection( Vlon( vi), Vlat( vi), grid%lambda_M, grid%phi_M, grid%alpha_stereo, Vx( vi), Vy( vi))
        END DO

      ELSEIF (region_name == 'GRL') THEN
        ! Greenland: grndrainagesystems_ekholm.txt
        ! Can be downloaded from: https: // earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems
        ! A text file with 7 header lines, followed by three columns of data:
        !    basin ID, Lat, Lon   (NOTE: here basin ID is a floating-point rather than an integer!)

      ! ===== Check if this is really the file we're reading =====
      ! ==========================================================

        recognised_file = .FALSE.
        n_header_lines  = 0
        n_vertices      = 0
        n_skip          = 0

        DO i = 1, 256-29
          IF (filename_basins(i:i+28) == 'grndrainagesystems_ekholm.txt') THEN
            recognised_file = .TRUE.
            n_header_lines  = 7
            n_skip          = 5 ! Since the polygon file is at a ridiculously high resolution, downscale it a bit for efficiency
            n_vertices      = CEILING( REAL(272965,dp) / REAL(n_skip,dp))
          END IF
        END DO

        IF ((.NOT. recognised_file) .AND. par%master) THEN
          CALL warning('for Greenland, we expect the file "grndrainagesystems_ekholm.txt". ' // &
                       'This can be downloaded from https: // earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems. ' // &
                       'If another file is used, make sure that it has 7 header lines, or alternatively just change the code in initialise_basins.')
        END IF

        ! Allocate shared memory
        CALL allocate_shared_dp_1D(  n_vertices, Vlat,   wVlat  )
        CALL allocate_shared_dp_1D(  n_vertices, Vlon,   wVlon  )
        CALL allocate_shared_dp_1D(  n_vertices, Vx,     wVx    )
        CALL allocate_shared_dp_1D(  n_vertices, Vy,     wVy    )
        CALL allocate_shared_int_1D( n_vertices, VID,    wVID   )

      ! ===== Read the file =====
      ! =========================

        IF (par%master) THEN

          ! Open the file
          OPEN( UNIT = 1337, FILE=filename_basins, ACTION='READ')

          ! Skip the header lines
          DO i = 1, n_header_lines
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
          END DO

          ! Read the actual data
          DO vi = 1, n_vertices
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) VID_dp, Vlat( vi), Vlon( vi)
            IF (ios /= 0) THEN
              CALL crash('length of text file "' // TRIM( filename_basins) // '" does not match n_vertices = {int_01}!', int_01 = n_vertices)
            END IF

            DO vj = 1, n_skip-1
              READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
            END DO

            ! Convert basin ID from floating point to integer
            IF     ( viD_dp == 1.1_dp) THEN
              VID( vi) = 1
            ELSEIF ( viD_dp == 1.2_dp) THEN
              VID( vi) = 2
            ELSEIF ( viD_dp == 1.3_dp) THEN
              VID( vi) = 3
            ELSEIF ( viD_dp == 1.4_dp) THEN
              VID( vi) = 4
            ELSEIF ( viD_dp == 2.1_dp) THEN
              VID( vi) = 5
            ELSEIF ( viD_dp == 2.2_dp) THEN
              VID( vi) = 6
            ELSEIF ( viD_dp == 3.1_dp) THEN
              VID( vi) = 7
            ELSEIF ( viD_dp == 3.2_dp) THEN
              VID( vi) = 8
            ELSEIF ( viD_dp == 3.3_dp) THEN
              VID( vi) = 9
            ELSEIF ( viD_dp == 4.1_dp) THEN
              VID( vi) = 10
            ELSEIF ( viD_dp == 4.2_dp) THEN
              VID( vi) = 11
            ELSEIF ( viD_dp == 4.3_dp) THEN
              VID( vi) = 12
            ELSEIF ( viD_dp == 5.0_dp) THEN
              VID( vi) = 13
            ELSEIF ( viD_dp == 6.1_dp) THEN
              VID( vi) = 14
            ELSEIF ( viD_dp == 6.2_dp) THEN
              VID( vi) = 15
            ELSEIF ( viD_dp == 7.1_dp) THEN
              VID( vi) = 16
            ELSEIF ( viD_dp == 7.2_dp) THEN
              VID( vi) = 17
            ELSEIF ( viD_dp == 8.1_dp) THEN
              VID( vi) = 18
            ELSEIF ( viD_dp == 8.2_dp) THEN
              VID( vi) = 19
            ELSE
              CALL crash('unrecognised floating-point basin ID in file "' // TRIM( filename_basins) // '"!')
            END IF

          END DO

          ! Close the file
          CLOSE( UNIT = 1337)

        END IF ! IF (par%master) THEN
        CALL sync

        ! Project [lat,lon] to [x,y]
        CALL partition_list( n_vertices, par%i, par%n, vi1, vi2)
        DO vi = vi1, vi2
          CALL oblique_sg_projection( Vlon( vi), Vlat( vi), grid%lambda_M, grid%phi_M, grid%alpha_stereo, Vx( vi), Vy( vi))
        END DO

      END IF ! IF (region_name == 'ANT') THEN

    ! ===== Fill in the basins =====
    ! ==============================

      IF (par%master) WRITE(0,*) '    Filling in the basins...'

      ! Allocate shared memory
      CALL allocate_shared_int_2D( grid%nx, grid%ny, basin_ID_loc, wbasin_ID_loc)

      ! Determine number of basins
      IF (par%master) nbasins = MAXVAL( VID)
      CALL sync

      DO bi = 1, nbasins

        ! Find range of vertices [vi1,vi2] for this basin
        vi1 = 1
        DO WHILE (VID( vi1) /= bi)
          vi1 = vi1 + 1
        END DO
        vi2 = vi1
        DO WHILE (VID( vi2) == bi .AND. vi2 < n_vertices)
          vi2 = vi2 + 1
        END DO
        vi2 = vi2 - 1

        ! Copy these vertices to a single array
        CALL allocate_shared_dp_2D( vi2+1-vi1, 2, poly_bi, wpoly_bi)
        IF (par%master) THEN
          poly_bi( :,1) = Vx( vi1:vi2)
          poly_bi( :,2) = Vy( vi1:vi2)
        END IF
        CALL sync

        ! Determine maximum polygon extent, for quick checks
        xmin = MINVAL( poly_bi(:,1))
        xmax = MAXVAL( poly_bi(:,1))
        ymin = MINVAL( poly_bi(:,2))
        ymax = MAXVAL( poly_bi(:,2))

        ! Check which grid cells lie inside the polygon spanned by these vertices
        DO j = 1, grid%ny
        DO i = grid%i1, grid%i2
          p = [grid%x( i), grid%y( j)]

          ! Quick test
          IF (p(1) < xmin .OR. p(1) > xmax .OR. &
              p(2) < ymin .OR. p(2) > ymax) THEN
            ! p cannot lie in the polygon, don't bother checking
          ElSE
            IF (is_in_polygon( poly_bi, p)) basin_ID_loc( i,j) = bi
          END IF

        END DO
        END DO
        CALL sync

        ! Clean up this basin's polygon
        CALL deallocate_shared( wpoly_bi)

      END DO ! DO bi = 1, nbasins

    ! ===== Extend basins into the ocean =====
    ! ========================================

      IF (par%master) WRITE(0,*) '    Extending basins into the ocean...'

      ! Allocate shared memory
      CALL allocate_shared_int_2D( grid%nx, grid%ny, basin_ID_ext_loc, wbasin_ID_ext_loc)

      ! Copy data
      basin_ID_ext_loc( grid%i1:grid%i2,:) = basin_ID_loc( grid%i1:grid%i2,:)
      CALL sync

      ! Compile list of ice-front pixels and their basin ID's
      IF (par%master) THEN
        n_front = 0
        DO j = 2, grid%ny-1
        DO i = 2, grid%nx-1
          IF (basin_ID_loc( i,j) > 0) THEN
            IF (MINVAL( basin_ID_loc( i-1:i+1,j-1:j+1)) == 0) THEN
              n_front = n_front + 1
            END IF
          END IF
        END DO
        END DO
      END IF
      CALL MPI_BCAST( n_front, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      CALL allocate_shared_int_2D( n_front, 2, ij_front, wij_front)

      IF (par%master) THEN
        k = 0
        DO j = 2, grid%ny-1
        DO i = 2, grid%nx-1
          IF (basin_ID_loc( i,j) > 0) THEN
            IF (MINVAL( basin_ID_loc( i-1:i+1,j-1:j+1)) == 0) THEN
              k = k + 1
              ij_front( k,:) = [i,j]
            END IF
          END IF
        END DO
        END DO
      END IF
      CALL sync

      ! For all non-assigned grid cells, find the nearest front cell and copy that basin ID
      DO j = 1, grid%ny
      DO i = grid%i1, grid%i2

        IF (basin_ID_ext_loc( i,j) == 0) THEN

          ! Go over all front cells, find the nearest
          dist_min  = REAL(MAX(grid%nx,grid%ny),dp) * grid%dx
          i_nearest = 0
          j_nearest = 0
          DO k = 1, n_front
            ii = ij_front( k,1)
            jj = ij_front( k,2)
            dist = SQRT( REAL(ii-i,dp)**2 + REAL(jj-j,dp)**2) * grid%dx
            IF (dist < dist_min) THEN
              dist_min = dist
              i_nearest = ii
              j_nearest = jj
            END IF
          END DO ! DO k = 1, n_front

          ! Assign basin ID of nearest front cell
          basin_ID_ext_loc( i,j) = basin_ID_loc( i_nearest, j_nearest)

        END IF

      END DO
      END DO
      CALL sync

    ! ===== Copy final result to the ice structure =====
    ! ==================================================

      basin_ID( grid%i1:grid%i2,:) = basin_ID_ext_loc( grid%i1:grid%i2,:)
      CALL sync

      ! Clean up after yourself
      CALL deallocate_shared( wVlat            )
      CALL deallocate_shared( wVlon            )
      CALL deallocate_shared( wVx              )
      CALL deallocate_shared( wVy              )
      CALL deallocate_shared( wVID             )
      CALL deallocate_shared( wbasin_ID_loc    )
      CALL deallocate_shared( wbasin_ID_ext_loc)
      CALL deallocate_shared( wij_front        )

    ! ==== If so specified, merge certain ice basins =====
    ! ====================================================

      IF (par%master .AND. (C%do_merge_basins_ANT .OR. &
                            C%do_merge_basins_GRL)) THEN
        WRITE(0,*) '    Merging specific basins...'
      END IF

      IF     (region_name == 'ANT' .AND. C%do_merge_basins_ANT) THEN
        CALL merge_basins_ANT_grid( grid, basin_ID, nbasins)
      ELSEIF (region_name == 'GRL' .AND. C%do_merge_basins_GRL) THEN
        CALL merge_basins_GRL_grid( grid, basin_ID, nbasins)
      END IF

    ELSE
      CALL crash('unknown choice_basin_scheme "' // TRIM(choice_basin_scheme) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basins_grid

  SUBROUTINE merge_basins_ANT_grid(  grid, basin_ID, nbasins)
    ! Merge some of the Antarctic basins to match the more common definitions

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(INOUT) :: basin_ID
    INTEGER,                             INTENT(INOUT) :: nbasins

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'merge_basins_ANT_grid'
    INTEGER                                            :: i,j,k
    INTEGER,  DIMENSION( 27)                           :: T
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  basin_ID_old
    INTEGER                                            :: wbasin_ID_old

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%nx, grid%ny, basin_ID_old, wbasin_ID_old)

    ! Copy data
    basin_ID_old( grid%i1:grid%i2,:) = basin_ID( grid%i1:grid%i2,:)
    CALL sync

    ! Define the merging table
    T(  1) = 1
    T(  2) = 1
    T(  3) = 1
    T(  4) = 2
    T(  5) = 3
    T(  6) = 4
    T(  7) = 5
    T(  8) = 6
    T(  9) = 6
    T( 10) = 6
    T( 11) = 6
    T( 12) = 7
    T( 13) = 8
    T( 14) = 9
    T( 15) = 10
    T( 16) = 11
    T( 17) = 11
    T( 18) = 11
    T( 19) = 11
    T( 20) = 12
    T( 21) = 13
    T( 22) = 13
    T( 23) = 13
    T( 24) = 14
    T( 25) = 15
    T( 26) = 16
    T( 27) = 17

    ! Perform the merge
    DO j = 1, grid%ny
    DO i = grid%i1, grid%i2
      DO k = 1, 27
        IF (basin_ID_old( i,j) == k) basin_ID( i,j) = T( k)
      END DO
    END DO
    END DO
    IF (par%master) nbasins = 17
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wbasin_ID_old)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE merge_basins_ANT_grid

  SUBROUTINE merge_basins_GRL_grid(  grid, basin_ID, nbasins)
    ! Merge some of the Greenland basins to match the more common definitions

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(INOUT) :: basin_ID
    INTEGER,                             INTENT(INOUT) :: nbasins

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'merge_basins_GRL_grid'
    INTEGER                                            :: i,j,k
    INTEGER,  DIMENSION( 19)                           :: T
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  basin_ID_old
    INTEGER                                            :: wbasin_ID_old

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%nx, grid%ny, basin_ID_old, wbasin_ID_old)

    ! Copy data
    basin_ID_old( grid%i1:grid%i2,:) = basin_ID( grid%i1:grid%i2,:)
    CALL sync

    ! Define the merging table
    T(  1) = 1
    T(  2) = 1
    T(  3) = 1
    T(  4) = 1
    T(  5) = 2
    T(  6) = 2
    T(  7) = 3
    T(  8) = 3
    T(  9) = 3
    T( 10) = 4
    T( 11) = 4
    T( 12) = 4
    T( 13) = 5
    T( 14) = 6
    T( 15) = 6
    T( 16) = 7
    T( 17) = 7
    T( 18) = 8
    T( 19) = 8

    ! Perform the merge
    DO j = 1, grid%ny
    DO i = grid%i1, grid%i2
      DO k = 1, 19
        IF (basin_ID_old( i,j) == k) basin_ID( i,j) = T( k)
      END DO
    END DO
    END DO
    IF (par%master) nbasins = 8
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wbasin_ID_old)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE merge_basins_GRL_grid

END MODULE general_ice_model_data_module