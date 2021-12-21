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

END MODULE general_ice_model_data_module
