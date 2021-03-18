MODULE general_ice_model_data_module

  USE mpi
  USE configuration_module,          ONLY: dp, C  
  USE parallel_module,               ONLY: par, sync, &
                                           allocate_shared_int_0D, allocate_shared_dp_0D, &
                                           allocate_shared_int_1D, allocate_shared_dp_1D, &
                                           allocate_shared_int_2D, allocate_shared_dp_2D, &
                                           allocate_shared_int_3D, allocate_shared_dp_3D, &
                                           deallocate_shared
  USE data_types_module,             ONLY: type_mesh, type_ice_model
  USE netcdf_module,                 ONLY: debug, write_to_debug_file
  USE mesh_ArakawaC_module,          ONLY: map_Aa_to_Ac, map_Aa_to_Ac_3D, get_mesh_derivatives_Ac
  USE mesh_derivatives_module,       ONLY: get_mesh_derivatives, get_mesh_curvatures
  USE zeta_module,                   ONLY: vertical_average, vertical_integrate
  USE parameters_module,             ONLY: seawater_density, ice_density 

  IMPLICIT NONE

CONTAINS
  
  ! Routines for calculating general ice model data - Hs, masks, ice physical properties
  SUBROUTINE update_general_ice_model_data( mesh, ice, time)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables
    INTEGER                                            :: vi, aci
 
    ! Map basic ice data to Ac mesh
    CALL map_Aa_to_Ac(    mesh, ice%Hi, ice%Hi_Ac)
    CALL map_Aa_to_Ac(    mesh, ice%Hb, ice%Hb_Ac)
    CALL map_Aa_to_Ac(    mesh, ice%SL, ice%SL_Ac)
    CALL map_Aa_to_Ac_3D( mesh, ice%Ti, ice%Ti_Ac)
    
    ! Update Hs
    DO vi = mesh%v1, mesh%v2
      ice%Hs(    vi ) = ice%Hi(    vi ) + MAX(ice%SL(    vi ) - ice_density / seawater_density * ice%Hi(    vi ), ice%Hb(    vi ))
    END DO
    CALL sync
    DO aci = mesh%ac1, mesh%ac2
      ice%Hs_Ac( aci) = ice%Hi_Ac( aci) + MAX(ice%SL_Ac( aci) - ice_density / seawater_density * ice%Hi_Ac( aci), ice%Hb_Ac( aci))
    END DO
    CALL sync
    
    ! Update dHs_dt
    ice%dHs_dt( mesh%v1:mesh%v2) = ice%dHb_dt( mesh%v1:mesh%v2) + ice%dHi_dt( mesh%v1:mesh%v2)
    CALL sync
    
    ! Determine masks
    CALL determine_masks( mesh, ice)
    
    ! Calculate surface slopes on Aa and Ac mesh
    CALL get_mesh_derivatives(    mesh, ice%Hi, ice%dHi_dx,    ice%dHi_dy)
    CALL get_mesh_derivatives(    mesh, ice%Hs, ice%dHs_dx,    ice%dHs_dy)
    CALL get_mesh_derivatives_Ac( mesh, ice%Hi, ice%dHi_dx_Ac, ice%dHi_dy_Ac, ice%dHi_dp_Ac, ice%dHi_do_Ac)
    CALL get_mesh_derivatives_Ac( mesh, ice%Hb, ice%dHb_dx_Ac, ice%dHb_dy_Ac, ice%dHb_dp_Ac, ice%dHb_do_Ac)
    CALL get_mesh_derivatives_Ac( mesh, ice%Hs, ice%dHs_dx_Ac, ice%dHs_dy_Ac, ice%dHs_dp_Ac, ice%dHs_do_Ac)
    CALL get_mesh_derivatives_Ac( mesh, ice%SL, ice%dSL_dx_Ac, ice%dSL_dy_Ac, ice%dSL_dp_Ac, ice%dSL_do_Ac)
    
    ! Use a different surface slope for shelves, to make sure the SSA gets the
    ! correct slopes around the discontinuity that is the grounding line.
    DO vi = mesh%v1, mesh%v2
      IF (ice%mask_sheet( vi)==1) THEN
        ice%dHs_dx_shelf( vi) = ice%dHs_dx( vi)
        ice%dHs_dy_shelf( vi) = ice%dHs_dy( vi)
      ELSE
        ice%dHs_dx_shelf( vi) = (1._dp - ice_density / seawater_density) * ice%dHi_dx( vi)
        ice%dHs_dy_shelf( vi) = (1._dp - ice_density / seawater_density) * ice%dHi_dy( vi)
      END IF
    END DO
    DO aci = mesh%ac1, mesh%ac2
      IF (ice%mask_sheet_Ac( aci)==1) THEN
        ice%dHs_dx_shelf_Ac( aci) = ice%dHs_dx_Ac( aci)
        ice%dHs_dy_shelf_Ac( aci) = ice%dHs_dy_Ac( aci)
      ELSE
        ice%dHs_dx_shelf_Ac( aci) = (1._dp - ice_density / seawater_density) * ice%dHi_dx_Ac( aci)
        ice%dHs_dy_shelf_Ac( aci) = (1._dp - ice_density / seawater_density) * ice%dHi_dy_Ac( aci)
      END IF
    END DO
    CALL sync
        
    ! Calculate physical properties on both Aa and Ac mesh
    CALL ice_physical_properties( mesh, ice, time)
    
  END SUBROUTINE update_general_ice_model_data
  SUBROUTINE determine_masks( mesh, ice)
    ! Determine the different masks, on both the Aa and the Ac mesh
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
  
    INTEGER                                       :: vi, ci, vc, aci, vj

    INTEGER, PARAMETER                            :: type_land           = 0
    INTEGER, PARAMETER                            :: type_ocean          = 1
    INTEGER, PARAMETER                            :: type_lake           = 2
    INTEGER, PARAMETER                            :: type_sheet          = 3
    INTEGER, PARAMETER                            :: type_shelf          = 4
    INTEGER, PARAMETER                            :: type_coast          = 5
    INTEGER, PARAMETER                            :: type_margin         = 6
    INTEGER, PARAMETER                            :: type_groundingline  = 7
    INTEGER, PARAMETER                            :: type_calvingfront   = 8
    
    ! Start out with land everywhere, fill in the rest based on input.
    ice%mask_land(      mesh%v1:mesh%v2  ) = 1
    ice%mask_land_Ac(   mesh%ac1:mesh%ac2) = 1
    ice%mask_ocean(     mesh%v1:mesh%v2  ) = 0
    ice%mask_ocean_Ac(  mesh%ac1:mesh%ac2) = 0
    ice%mask_lake(      mesh%v1:mesh%v2  ) = 0
    ice%mask_lake_Ac(   mesh%ac1:mesh%ac2) = 0
    ice%mask_ice(       mesh%v1:mesh%v2  ) = 0
    ice%mask_ice_Ac(    mesh%ac1:mesh%ac2) = 0
    ice%mask_sheet(     mesh%v1:mesh%v2  ) = 0
    ice%mask_sheet_Ac(  mesh%ac1:mesh%ac2) = 0
    ice%mask_shelf(     mesh%v1:mesh%v2  ) = 0
    ice%mask_shelf_Ac(  mesh%ac1:mesh%ac2) = 0
    ice%mask_coast(     mesh%v1:mesh%v2  ) = 0
    ice%mask_coast_Ac(  mesh%ac1:mesh%ac2) = 0
    ice%mask_margin(    mesh%v1:mesh%v2  ) = 0
    ice%mask_margin_Ac( mesh%ac1:mesh%ac2) = 0
    ice%mask_gl(        mesh%v1:mesh%v2  ) = 0
    ice%mask_gl_Ac(     mesh%ac1:mesh%ac2) = 0
    ice%mask_cf(        mesh%v1:mesh%v2  ) = 0
    ice%mask_cf_Ac(     mesh%ac1:mesh%ac2) = 0
    ice%mask(           mesh%v1:mesh%v2  ) = type_land
    ice%mask_Ac(        mesh%ac1:mesh%ac2) = type_land
    CALL sync
    
    ! Start out with land everywhere, fill in the rest based on input.
    CALL sync
  
    ! First on the Aa mesh
    ! ====================
  
    DO vi = mesh%v1, mesh%v2
    
      ! Determine ice
      IF (ice%Hi( vi) > 0._dp) THEN
        ice%mask_ice( vi)  = 1
        ice%mask_land( vi) = 0
      END IF
      
      ! Determine sheet
      IF (ice%mask_ice( vi) == 1 .AND. (ice%Hi( vi) > (ice%SL( vi) - ice%Hb( vi)) * seawater_density/ice_density)) THEN
        ice%mask_sheet( vi) = 1
        ice%mask( vi)       = type_sheet
      END IF
    
      ! Determine shelf
      IF (ice%mask_ice( vi) == 1 .AND. ice%mask_sheet( vi) == 0) THEN
        ice%mask_shelf( vi) = 1
        ice%mask( vi)       = type_shelf
      END IF
      
      ! Determine open ocean
      IF ((ice%Hb( vi) < ice%SL( vi)) .AND. (ice%Hi( vi) == 0._dp)) THEN
        ice%mask_ocean( vi) = 1
        ice%mask_land( vi)  = 0
        ice%mask( vi)       = type_ocean
      END IF
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
  
    ! Determine coast, grounding line and calving front
    DO vi = mesh%v1, mesh%v2
    
      IF (ice%mask_land( vi) == 1) THEN  
        ! Land bordering ocean equals coastline
        
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ocean( vc) == 1) THEN
            ice%mask( vi) = type_coast
            ice%mask_coast( vi) =  1
          END IF
        END DO
      END IF
    
      IF (ice%mask_ice( vi) == 1) THEN  
        ! Ice bordering non-ice equals margin
        
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ice( vc) == 0) THEN
            ice%mask( vi) = type_margin
            ice%mask_margin( vi) =  1
          END IF
        END DO
      END IF
    
      IF (ice%mask_sheet( vi) == 1) THEN  
        ! Sheet bordering shelf equals groundingline
        
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_shelf( vc) == 1) THEN
            ice%mask( vi) = type_groundingline
            ice%mask_gl( vi) = 1
          END IF
        END DO
      END IF
  
      IF (ice%mask_ice( vi) == 1) THEN  
        ! Ice (sheet or shelf) bordering ocean equals calvingfront
        
        DO ci = 1, mesh%nC(vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_ocean( vc) == 1) THEN
            ice%mask( vi) = type_calvingfront
            ice%mask_cf( vi) = 1
          END IF
        END DO
        
      END IF
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! Then on the Ac mesh
    ! ===================
  
    ! Determine ocean (open ocean + shelf)
    DO aci = mesh%ac1, mesh%ac2
      IF ((ice%Hb_Ac( aci) < ice%SL_Ac( aci)) .AND. (ice%Hi_Ac( aci) < (ice%SL_Ac( aci) - ice%Hb_Ac( aci)) * seawater_density/ice_density)) THEN
        ice%mask_ocean_Ac( aci) = 1
        ice%mask_Ac( aci)       = type_ocean
      END IF
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync
  
    ! Determine shelf
    DO aci = mesh%ac1, mesh%ac2
      IF (ice%mask_ocean_Ac( aci) == 1 .AND. ice%Hi_Ac( aci) > 0._dp) THEN
        ice%mask_shelf_Ac( aci) = 1
        ice%mask_ice_Ac( aci)   = 1
        ice%mask_Ac( aci)       = type_shelf
      END IF
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync
  
    ! Determine sheet
    DO aci = mesh%ac1, mesh%ac2
      IF (ice%mask_ocean_Ac( aci) == 0 .AND. ice%Hi_Ac( aci) > 0._dp) THEN
        ice%mask_sheet_Ac( aci) = 1
        ice%mask_ice_Ac( aci)   = 1
        ice%mask_Ac( aci)       = type_sheet
      END IF
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync
  
    ! Determine coast, margin, grounding line and calving front
    DO aci = mesh%ac1, mesh%ac2
    
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      
      IF ((ice%mask_land( vi) == 1 .AND. ice%mask_ocean( vj) == 1) .OR. &
          (ice%mask_land( vj) == 1 .AND. ice%mask_ocean( vi) == 1)) THEN
        ! Land bordering ocean equals coast
        ice%mask_Ac( aci) = type_coast
        ice%mask_coast_Ac( aci) = 1
      END IF
      
      IF ((ice%mask_ice( vi) == 1 .AND. ice%mask_ice( vj) == 0) .OR. &
          (ice%mask_ice( vj) == 1 .AND. ice%mask_ice( vi) == 0)) THEN
        ! Ice bordering non-ice equals margin
        ice%mask_Ac( aci) = type_margin
        ice%mask_margin_Ac( aci) = 1
      END IF
      
      IF ((ice%mask_sheet( vi) == 1 .AND. ice%mask_shelf( vj) == 1) .OR. &
          (ice%mask_sheet( vj) == 1 .AND. ice%mask_shelf( vi) == 1)) THEN
        ! Sheet bordering shelf equals groundingline
        ice%mask_Ac( aci) = type_groundingline
        ice%mask_gl_Ac( aci) = 1
      END IF
      
      IF ((ice%mask_ice( vi) == 1 .AND. ice%mask_shelf( vj) == 0 .AND. ice%mask_ocean( vj) == 1) .OR. &
          (ice%mask_ice( vj) == 1 .AND. ice%mask_shelf( vi) == 0 .AND. ice%mask_ocean( vi) == 1)) THEN
        ! Ice bordering open ocean equals calvingfront
        ice%mask_Ac( aci) = type_calvingfront
        ice%mask_cf_Ac( aci) = 1
      END IF
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync
  
  END SUBROUTINE determine_masks
  SUBROUTINE ice_physical_properties( mesh, ice, time)
    ! Calculate the pressure melting point, flow parameter, specific heat and thermal conductivity of the ice.
      
    USE parameters_module, ONLY: T0, CC, SMT, sec_per_year
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: time
  
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: vi, ci, k
    REAL(dp)                                           :: Ti_mean     ! Mean ice temperature at the shelf [K]
    REAL(dp), DIMENSION(C%nZ)                          :: prof
    
    REAL(dp)                                           :: A_flow_MISMIP
    
    REAL(dp), PARAMETER                                :: A_low_temp  = 1.14E-05_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: A_high_temp = 5.47E+10_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_low_temp  = 6.0E+04_dp    ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_high_temp = 13.9E+04_dp   ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: R_gas       = 8.314_dp      ! Gas constant [J mol^-1 K^-1]
    
    ! If we're doing one of the EISMINT experiments, use fixed values
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test') THEN
          
        ice%A_flow(         mesh%v1:mesh%v2,   :) = 1.0E-16_dp
        ice%A_flow_Ac(      mesh%ac1:mesh%ac2, :) = 1.0E-16_dp
        ice%A_flow_mean(    mesh%v1:mesh%v2     ) = 1.0E-16_dp
        ice%A_flow_mean_Ac( mesh%ac1:mesh%ac2   ) = 1.0E-16_dp
        ice%Ki(             mesh%v1:mesh%v2,   :) = 2.1_dp * sec_per_year
        ice%Cpi(            mesh%v1:mesh%v2,   :) = 2009._dp
        
        DO vi = mesh%v1, mesh%v2
          DO k = 1, C%nZ
            ice%Ti_pmp( vi,k) = T0 - (C%zeta(k) * ice%Hi( vi) * 8.7E-04_dp)
          END DO
        END DO    
        CALL sync
        
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
      
        A_flow_MISMIP = 1.0E-16_dp
        IF     (time < 25000._dp) THEN
          A_flow_MISMIP = 1.0E-16_dp
        ELSEIF (time < 50000._dp) THEN
          A_flow_MISMIP = 1.0E-17_dp
        ELSEIF (time < 75000._dp) THEN
          A_flow_MISMIP = 1.0E-16_dp
        END IF
        
        ice%A_flow(         mesh%v1:mesh%v2,   :) = A_flow_MISMIP
        ice%A_flow_Ac(      mesh%ac1:mesh%ac2, :) = A_flow_MISMIP
        ice%A_flow_mean(    mesh%v1:mesh%v2     ) = A_flow_MISMIP
        ice%A_flow_mean_Ac( mesh%ac1:mesh%ac2   ) = A_flow_MISMIP
        CALL sync
        
        RETURN
        
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in ice_physical_properties!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
  
    ! First on the Aa mesh
    ! ====================
    
    ice%A_flow(         mesh%v1:mesh%v2,   :) = 0._dp
    ice%A_flow_Ac(      mesh%ac1:mesh%ac2, :) = 0._dp
    ice%A_flow_mean(    mesh%v1:mesh%v2     ) = 0._dp
    ice%A_flow_mean_Ac( mesh%ac1:mesh%ac2   ) = 0._dp
    ice%Ti_pmp(         mesh%v1:mesh%v2,   :) = 0._dp
    ice%Cpi(            mesh%v1:mesh%v2,   :) = 0._dp
    ice%Ki(             mesh%v1:mesh%v2,   :) = 0._dp
    CALL sync
    
    DO vi = mesh%v1, mesh%v2
      ! Calculate the pressure melting point temperature (= the maximum temperature) for each depth (see equation (11.2)):
      ice%Ti_pmp( vi,:) = T0 - CC * ice%Hi( vi) * C%zeta
      
      DO k = 1, C%nZ 

        ! Calculation of the flow parameter at the sheet and groundline as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        IF (ice%Ti( vi,k) < 263.15_dp) THEN
          ice%A_flow( vi,k) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti( vi,k)))  
        ELSE
          ice%A_flow( vi,k) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti( vi,k)))  
        END IF
           
        ! Calculation of the parameterization of the specific heat capacity of ice, based on Pounder (1965):
        ice%Cpi( vi,k) = 2115.3_dp + 7.79293_dp * (ice%Ti( vi,k) - T0)  ! See equation (11.9)
           
        ! Calculation of the parameterization of the thermal conductivity of ice, based on Ritz (1987):
        ice%Ki( vi,k)  = 3.101E+08_dp * EXP(-0.0057_dp * ice%Ti( vi,k)) ! See equation (11.5), Huybrechts (4.40)
      END DO 

      IF (ice%mask_sheet( vi) == 1) THEN
        ! Calculation of the vertical average flow parameter at the sheet and groundline
        prof = ice%A_flow( vi,:)
        ice%A_flow_mean( vi) = vertical_average(prof)
      ELSE
        ! Calculation of the flow parameter at the shelf as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        Ti_mean = (ice%Ti( vi,1) + SMT) / 2._dp
        IF (Ti_mean < 263.15_dp) THEN
          ice%A_flow_mean( vi) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * Ti_mean))  
        ELSE
          ice%A_flow_mean( vi) = A_high_temp * EXP(-Q_high_temp / (R_gas * Ti_mean))  
        END IF
      END IF  
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
  
    ! Then on the Ac mesh
    ! ===================
    
    DO ci = mesh%ac1, mesh%ac2
      
      DO k = 1, C%nZ 

        ! Calculation of the flow parameter at the sheet and groundline as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        IF (ice%Ti_Ac( ci,k) < 263.15_dp) THEN
          ice%A_flow_Ac( ci,k) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti_Ac( ci,k)))  
        ELSE
          ice%A_flow_Ac( ci,k) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti_Ac( ci,k)))  
        END IF
        
      END DO 

      IF (ice%mask_sheet_Ac( ci) == 1) THEN
        ! Calculation of the vertical average flow parameter at the sheet and groundline
        prof = ice%A_flow_Ac( ci,:)
        ice%A_flow_mean_Ac( ci) = vertical_average(prof)
      ELSE
        ! Calculation of the flow parameter at the shelf as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        Ti_mean = (ice%Ti_Ac( ci,1) + SMT) / 2._dp
        IF (Ti_mean < 263.15_dp) THEN
          ice%A_flow_mean_Ac( ci) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * Ti_mean))  
        ELSE
          ice%A_flow_mean_Ac( ci) = A_high_temp * EXP(-Q_high_temp / (R_gas * Ti_mean))  
        END IF
      END IF  
      
    END DO ! DO ci = mesh%ac1, mesh%ac2
    CALL sync
    
  END SUBROUTINE ice_physical_properties
  
  FUNCTION is_floating( Hi, Hb, SL) RESULT( isso)
    ! The flotation criterion
      
    IMPLICIT NONE
    
    REAL(dp),                            INTENT(IN)    :: Hi, Hb, SL
    LOGICAL                                            :: isso
    
    isso = .FALSE.
    IF (Hi < (SL - Hb) * seawater_density/ice_density) isso = .TRUE.
    
  END FUNCTION is_floating

END MODULE general_ice_model_data_module
