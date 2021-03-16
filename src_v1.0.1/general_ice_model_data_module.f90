MODULE general_ice_model_data_module

  USE mpi
  USE configuration_module,        ONLY: dp, C  
  USE parallel_module,             ONLY: par, sync, &
                                         allocate_shared_memory_int_0D, allocate_shared_memory_dp_0D, &
                                         allocate_shared_memory_int_1D, allocate_shared_memory_dp_1D, &
                                         allocate_shared_memory_int_2D, allocate_shared_memory_dp_2D, &
                                         allocate_shared_memory_int_3D, allocate_shared_memory_dp_3D, &
                                         deallocate_shared_memory
  USE data_types_module,           ONLY: type_mesh, type_ice_model
  USE mesh_ArakawaC_module,        ONLY: MapAaToAc, MapAaToAc_3D, GetAcMeshDerivatives
  USE mesh_derivatives_module,     ONLY: GetMeshDerivatives, GetMeshCurvatures
  USE zeta_module,                 ONLY: vertical_average, vertical_integrate
  USE parameters_module,           ONLY: seawater_density, ice_density 

  IMPLICIT NONE

CONTAINS
  
  ! Routines for calculating general ice model data - Hs, masks, ice physical properties
  SUBROUTINE update_general_ice_model_data( ice, mesh)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    
    ! Local variables
    INTEGER                                            :: vi, ci
    
    ! Update Hs
    DO vi = mesh%v1, mesh%v2
      ice%Hs(vi) = ice%Hi(vi) + MAX(ice%sealevel(vi) - ice_density / seawater_density * ice%Hi(vi), ice%Hb(vi))
    END DO
    CALL sync
    
    ! Update dHs_dt
    ice%dHs_dt(mesh%v1:mesh%v2) = ice%dHb_dt(mesh%v1:mesh%v2) + ice%dHi_dt(mesh%v1:mesh%v2)
    CALL sync
 
    ! Map basic ice data to Ac mesh
    CALL MapAaToAc(    mesh, ice%Hi,       ice%Hi_ac          )
    CALL MapAaToAc(    mesh, ice%Hb,       ice%Hb_ac          )
    CALL MapAaToAc(    mesh, ice%sealevel, ice%sealevel_ac    )
    CALL MapAaToAc_3D( mesh, ice%Ti,       ice%Ti_ac          )
    
    DO ci = mesh%c1, mesh%c2
      ice%Hs_ac(ci) = ice%Hi_ac(ci) + MAX(ice%sealevel_ac(ci) - ice_density / seawater_density * ice%Hi_ac(ci), ice%Hb_ac(ci))
    END DO
    CALL sync
    
    ! Determine masks
    CALL determine_masks(ice, mesh)
    
    ! Calculate surface slopes on Aa and Ac mesh
    CALL GetMeshDerivatives(   mesh, ice%Hi, ice%dHi_dx,    ice%dHi_dy)
    CALL GetMeshDerivatives(   mesh, ice%Hs, ice%dHs_dx,    ice%dHs_dy)
    CALL GetAcMeshDerivatives( mesh, ice%Hs, ice%dHs_dx_ac, ice%dHs_dy_ac, ice%dHs_dp_ac, ice%dHs_do_ac)
    
    ! Use a different surface slope for shelves, to make sure the SSA gets the
    ! correct slopes around the discontinuity that is the grounding line.
    DO vi = mesh%v1, mesh%v2
      IF (ice%mask_sheet(vi)==1) THEN
        ice%dHs_dx_shelf(vi) = ice%dHs_dx(vi)
        ice%dHs_dy_shelf(vi) = ice%dHs_dy(vi)
      ELSE
        ice%dHs_dx_shelf(vi) = (1._dp - ice_density / seawater_density) * ice%dHi_dx(vi)
        ice%dHs_dy_shelf(vi) = (1._dp - ice_density / seawater_density) * ice%dHi_dy(vi)
      END IF
    END DO
    CALL sync
        
    ! Calculate physical properties on both Aa and Ac mesh
    CALL ice_physical_properties( ice, mesh)
    
  END SUBROUTINE update_general_ice_model_data
  SUBROUTINE determine_masks( ice, mesh)
    ! Determine the different masks, on both the Aa and the Ac mesh
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh  
  
    INTEGER                                       :: vi, cc, vn, ci, vj

    INTEGER, PARAMETER                            :: type_land           = 0
    INTEGER, PARAMETER                            :: type_ocean          = 1
    INTEGER, PARAMETER                            :: type_lake           = 2
    INTEGER, PARAMETER                            :: type_sheet          = 3
    INTEGER, PARAMETER                            :: type_shelf          = 4
    INTEGER, PARAMETER                            :: type_coast          = 5
    INTEGER, PARAMETER                            :: type_margin         = 6
    INTEGER, PARAMETER                            :: type_groundingline  = 7
    INTEGER, PARAMETER                            :: type_calvingfront   = 8
  
    ! First on the Aa mesh
    ! ====================
    
    ! Start out with land everywhere, fill in the rest based on input.
    ice%mask(              mesh%v1:mesh%v2) = type_land
    ice%mask_land(         mesh%v1:mesh%v2) = 1
    ice%mask_ocean(        mesh%v1:mesh%v2) = 0
    ice%mask_lake(         mesh%v1:mesh%v2) = 0
    ice%mask_ice(          mesh%v1:mesh%v2) = 0
    ice%mask_sheet(        mesh%v1:mesh%v2) = 0
    ice%mask_shelf(        mesh%v1:mesh%v2) = 0
    ice%mask_coast(        mesh%v1:mesh%v2) = 0
    ice%mask_margin(       mesh%v1:mesh%v2) = 0
    ice%mask_groundingline(mesh%v1:mesh%v2) = 0
    ice%mask_calvingfront( mesh%v1:mesh%v2) = 0
    CALL sync
  
    DO vi = mesh%v1, mesh%v2
    
      ! Determine ice
      IF (ice%Hi(vi) > 0._dp) THEN
        ice%mask_ice(vi)  = 1
        ice%mask_land(vi) = 0
      END IF
      
      ! Determine sheet
      IF (ice%mask_ice(vi) == 1 .AND. (ice%Hi(vi) > (ice%sealevel(vi) - ice%Hb(vi)) * seawater_density/ice_density)) THEN
        ice%mask_sheet(vi) = 1
        ice%mask(vi)       = type_sheet
      END IF
    
      ! Determine shelf
      IF (ice%mask_ice(vi) == 1 .AND. ice%mask_sheet(vi) == 0) THEN
        ice%mask_shelf(vi) = 1
        ice%mask(vi)       = type_shelf
      END IF
      
      ! Determine open ocean
      IF ((ice%Hb(vi) < ice%sealevel(vi)) .AND. (ice%Hi(vi) == 0._dp)) THEN
        ice%mask_ocean(vi) = 1
        ice%mask_land(vi)  = 0
        ice%mask(vi)       = type_ocean
      END IF
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
  
    ! Determine coast, grounding line and calving front
    DO vi = mesh%v1, mesh%v2
    
      IF (ice%mask_land(vi) == 1) THEN  
        ! Land bordering ocean equals coastline
        
        DO cc = 1, mesh%nC(vi)
          vn = mesh%C(vi,cc)
          IF (ice%mask_ocean(vn) == 1) THEN
            ice%mask(vi) = type_coast
            ice%mask_coast(vi) =  1
          END IF
        END DO
      END IF
    
      IF (ice%mask_ice(vi) == 1) THEN  
        ! Ice bordering non-ice equals margin
        
        DO cc = 1, mesh%nC(vi)
          vn = mesh%C(vi,cc)
          IF (ice%mask_ice(vn) == 0) THEN
            ice%mask(vi) = type_margin
            ice%mask_margin(vi) =  1
          END IF
        END DO
      END IF
    
      IF (ice%mask_sheet(vi) == 1) THEN  
        ! Sheet bordering shelf equals groundingline
        
        DO cc = 1, mesh%nC(vi)
          vn = mesh%C(vi,cc)
          IF (ice%mask_shelf(vn) == 1) THEN
            ice%mask(vi) = type_groundingline
            ice%mask_groundingline(vi) = 1
          END IF
        END DO
      END IF
  
      IF (ice%mask_ice(vi) == 1) THEN  
        ! Ice (sheet or shelf) bordering ocean equals calvingfront
        
        DO cc = 1, mesh%nC(vi)
          vn = mesh%C(vi,cc)
          IF (ice%mask_ocean(vn) == 1) THEN
            ice%mask(vi) = type_calvingfront
            ice%mask_calvingfront(vi) = 1
          END IF
        END DO
        
      END IF
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! Then on the Ac mesh
    ! ===================
    
    ! Start out with land everywhere, fill in the rest based on input.
    ice%mask_ac(               mesh%c1:mesh%c2) = type_land
    ice%mask_ice_ac(           mesh%c1:mesh%c2) = 0
    ice%mask_sheet_ac(         mesh%c1:mesh%c2) = 0
    ice%mask_shelf_ac(         mesh%c1:mesh%c2) = 0
    ice%mask_ocean_ac(         mesh%c1:mesh%c2) = 0
    ice%mask_groundingline_ac( mesh%c1:mesh%c2) = 0
    ice%mask_calvingfront_ac(  mesh%c1:mesh%c2) = 0
    CALL sync
  
    ! Determine ocean (open ocean + shelf)
    DO ci = mesh%c1, mesh%c2
      IF ((ice%Hb_ac(ci) < ice%sealevel_ac(ci)) .AND. (ice%Hi_ac(ci) < (ice%sealevel_ac(ci) - ice%Hb_ac(ci)) * seawater_density/ice_density)) THEN
        ice%mask_ocean_ac(ci) = 1
        ice%mask_ac(ci)       = type_ocean
      END IF
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
  
    ! Determine shelf
    DO ci = mesh%c1, mesh%c2
      IF (ice%mask_ocean_ac(ci) == 1 .AND. ice%Hi_ac(ci) > 0._dp) THEN
        ice%mask_shelf_ac(ci) = 1
        ice%mask_ice_ac(ci)   = 1
        ice%mask_ac(ci)       = type_shelf
      END IF
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
  
    ! Determine sheet
    DO ci = mesh%c1, mesh%c2
      IF (ice%mask_ocean_ac(ci) == 0 .AND. ice%Hi_ac(ci) > 0._dp) THEN
        ice%mask_sheet_ac(ci) = 1
        ice%mask_ice_ac(ci)   = 1
        ice%mask_ac(ci)       = type_sheet
      END IF
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
  
    ! Determine grounding line and calving front
    DO ci = mesh%c1, mesh%c2
      vi = mesh%Aci(ci,1)
      vj = mesh%Aci(ci,2)
      
      IF ((ice%mask_sheet(vi)==1 .AND. ice%mask_shelf(vj)==1) .OR. &
          (ice%mask_sheet(vj)==1 .AND. ice%mask_shelf(vj)==1)) THEN
        ! Sheet bordering shelf equals groundingline
        ice%mask_ac(ci) = type_groundingline
        ice%mask_groundingline_ac(ci) = 1
      END IF
      
      IF ((ice%mask_ice(vi)==1 .AND. ice%mask_shelf(vj)==0 .AND. ice%mask_ocean(vj)==1) .OR. &
          (ice%mask_ice(vj)==1 .AND. ice%mask_shelf(vi)==0 .AND. ice%mask_ocean(vi)==1)) THEN
        ! Ice bordering open ocean equals calvingfront
        ice%mask_ac(ci) = type_calvingfront
        ice%mask_calvingfront_ac(ci) = 1
      END IF
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
  
  END SUBROUTINE determine_masks
  SUBROUTINE ice_physical_properties( ice, mesh)
    ! Calculate the pressure melting point, flow parameter, specific heat and thermal conductivity of the ice.
      
    USE parameters_module, ONLY: T0, CC, SMT
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh     
  
    ! Local variables:
    INTEGER                                          :: vi, ci, k
    REAL(dp)                                         :: Ti_mean     ! Mean ice temperature at the shelf [K]
    REAL(dp), DIMENSION(C%nZ)                        :: prof
    REAL(dp), PARAMETER         :: A_low_temp  = 1.14E-05_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER         :: A_high_temp = 5.47E+10_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER         :: Q_low_temp  = 6.0E+04_dp    ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER         :: Q_high_temp = 13.9E+04_dp   ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER         :: R_gas       = 8.314_dp      ! Gas constant [J mol^-1 K^-1]
    
    ! If we're doing one of the EISMINT experiments, use fixed values
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler') THEN
          
        ice%A_flow(         mesh%v1:mesh%v2, :) = 1.0E-16_dp
        ice%A_flow_mean(    mesh%v1:mesh%v2   ) = 1.0E-16_dp
        ice%A_flow_ac(      mesh%c1:mesh%c2, :) = 1.0E-16_dp
        ice%A_flow_mean_ac( mesh%c1:mesh%c2   ) = 1.0E-16_dp
        ice%Ki(             mesh%v1:mesh%v2, :) = 2.1_dp * C%sec_per_year
        ice%Ki_ac(          mesh%c1:mesh%c2, :) = 2.1_dp * C%sec_per_year
        ice%Cpi(            mesh%v1:mesh%v2, :) = 2009._dp
        ice%Cpi_ac(         mesh%c1:mesh%c2, :) = 2009._dp
        
        DO vi = mesh%v1, mesh%v2
          DO k = 1, C%nZ
            ice%Ti_pmp(    vi,k) = T0 - (C%zeta(k) * ice%Hi(   vi) * 8.7E-04_dp)
          END DO
        END DO
        DO ci = mesh%c1, mesh%c2
          DO k = 1, C%nZ
            ice%Ti_pmp_ac( ci,k) = T0 - (C%zeta(k) * ice%Hi_ac(ci) * 8.7E-04_dp)
          END DO
        END DO      
        CALL sync
        
        RETURN
        
      ELSE
        WRITE(0,*) '  ERROR: "', C%choice_benchmark_experiment, '" is not a valid choice_benchmark_experiment!'
        STOP
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
  
    ! First on the Aa mesh
    ! ====================
    
    ! All variables are initialized at zero, Ti_pmp is caclculated everywhere
    ice%A_flow(      mesh%v1:mesh%v2, :) = 0._dp
    ice%A_flow_mean( mesh%v1:mesh%v2   ) = 0._dp
    ice%Cpi(         mesh%v1:mesh%v2, :) = 0._dp
    ice%Ki(          mesh%v1:mesh%v2, :) = 0._dp
    CALL sync
    
    DO vi = mesh%v1, mesh%v2
      ! Calculate the pressure melting point temperature (= the maximum temperature) for each depth (see equation (11.2)):
      ice%Ti_pmp(vi,:) = T0 - CC * ice%Hi(vi) * C%zeta(:)
      
      DO k = 1, C%NZ 

        ! Calculation of the flow parameter at the sheet and groundline as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        IF (ice%Ti(vi,k) < 263.15_dp) THEN
          ice%A_flow(vi,k) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti(vi,k)))  
        ELSE
          ice%A_flow(vi,k) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti(vi,k)))  
        END IF
           
        ! Calculation of the parameterization of the specific heat capacity of ice, based on Pounder (1965):
        ice%Cpi(vi,k) = 2115.3_dp + 7.79293_dp * (ice%Ti(vi,k) - T0)  ! See equation (11.9)
           
        ! Calculation of the parameterization of the thermal conductivity of ice, based on Ritz (1987):
        ice%Ki(vi,k)  = 3.101E+08_dp * EXP(-0.0057_dp * ice%Ti(vi,k)) ! See equation (11.5), Huybrechts (4.40)
      END DO 

      IF (ice%mask_sheet(vi) == 1) THEN
        ! Calculation of the vertical average flow parameter at the sheet and groundline
        prof = ice%A_flow(vi,:)
        ice%A_flow_mean(vi) = vertical_average(prof)
      ELSE
        ! Calculation of the flow parameter at the shelf as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        Ti_mean = (ice%Ti(vi,1) + SMT) / 2._dp
        IF (Ti_mean < 263.15_dp) THEN
          ice%A_flow_mean(vi) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * Ti_mean))  
        ELSE
          ice%A_flow_mean(vi) = A_high_temp * EXP(-Q_high_temp / (R_gas * Ti_mean))  
        END IF
      END IF  
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
  
    ! Then on the Ac mesh
    ! ===================
    
    ! All variables are initialized at zero, Ti_pmp is caclculated everywhere
    ice%A_flow_ac(      mesh%c1:mesh%c2, :) = 0._dp
    ice%A_flow_mean_ac( mesh%c1:mesh%c2   ) = 0._dp
    ice%Cpi_ac(         mesh%c1:mesh%c2, :) = 0._dp
    ice%Ki_ac(          mesh%c1:mesh%c2, :) = 0._dp
    CALL sync
    
    DO ci = mesh%c1, mesh%c2
      ! Calculate the pressure melting point temperature (= the maximum temperature) for each depth (see equation (11.2)):
      ice%Ti_pmp_ac(ci,:) = T0 - CC * ice%Hi_ac(ci) * C%zeta(:)
      
      DO k = 1, C%NZ 

        ! Calculation of the flow parameter at the sheet and groundline as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        IF (ice%Ti_ac(ci,k) < 263.15_dp) THEN
          ice%A_flow_ac(ci,k) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti_ac(ci,k)))  
        ELSE
          ice%A_flow_ac(ci,k) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti_ac(ci,k)))  
        END IF
         
        ! Calculation of the parameterization of the specific heat capacity of ice, based on Pounder (1965):
        ice%Cpi_ac(ci,k) = 2115.3_dp + 7.79293_dp * (ice%Ti_ac(ci,k) - T0)  ! See equation (11.9)
           
        ! Calculation of the parameterization of the thermal conductivity of ice, based on Ritz (1987):
        ice%Ki_ac(ci,k)  = 3.101E+08_dp * EXP(-0.0057_dp * ice%Ti_ac(ci,k)) ! See equation (11.5), Huybrechts (4.40)
      END DO 

      IF (ice%mask_sheet_ac(ci) == 1) THEN
        ! Calculation of the vertical average flow parameter at the sheet and groundline
        prof = ice%A_flow_ac(ci,:)
        ice%A_flow_mean_Ac(ci) = vertical_average(prof)
      ELSE
        ! Calculation of the flow parameter at the shelf as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        Ti_mean = (ice%Ti_ac(ci,1) + SMT) / 2._dp
        IF (Ti_mean < 263.15_dp) THEN
          ice%A_flow_mean_ac(ci) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * Ti_mean))  
        ELSE
          ice%A_flow_mean_ac(ci) = A_high_temp * EXP(-Q_high_temp / (R_gas * Ti_mean))  
        END IF
      END IF  
      
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
    
  END SUBROUTINE ice_physical_properties 
  SUBROUTINE basal_yield_stress( ice, mesh)
    ! Calculate the basal yield stress, which is used to determine the basal stress, used 
    ! for sliding. Using the parameterisations given by Martin et al. (TCD: 2010) 
    ! As used in the PISM-PIK model
    
    USE parameters_module, ONLY : ice_density, grav
    
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    
    ! local variables
    REAL(dp), DIMENSION(mesh%nV)                        :: lambda_p      ! scaling of the pore water pressure
    REAL(dp), DIMENSION(mesh%nV)                        :: P_wat         ! pore water pressure [Pa]
    REAL(dp), DIMENSION(mesh%nV)                        :: phi_fric      ! the friction angle (degrees)
    
    REAL(dp), PARAMETER                                 :: pf1   = -1000._dp
    REAL(dp), PARAMETER                                 :: pf2   = 0._dp
    REAL(dp), PARAMETER                                 :: p_min = 5._dp
    REAL(dp), PARAMETER                                 :: p_max = 20._dp 
    
    INTEGER                                             :: vi
    
    DO vi = mesh%v1, mesh%v2
      ! the pore water pressure is scaled with a bedrock height dependend parameterisation
      ! Equation (13) in Martin et al. (2011)
      lambda_p(vi) = MAX(0._dp, MIN(1._dp, (1._dp - (ice%Hb(vi) - ice%sealevel(vi)) / 1000._dp)))
      
      ! The pore water pressure, equation (12) in Martin et al. (2011)
      P_wat(vi) = 0.96_dp * lambda_p(vi) * ice_density * grav
      
      ! The friction angle, used for the yield stress, equation (11) in Martin et al. (2011)
      phi_fric(vi) = MAX(p_min, MIN(p_max, (p_min + (p_max - p_min) * (1._dp + (ice%Hb(vi) - pf2) / (pf2 - pf1))) ))
      
      IF (ice%mask_ocean(vi) == 1 .OR. ice%mask_shelf(vi) == 1) THEN
        ! set the yield stress to zero for ocean or shelf grid points
        ice%tau_yield(vi) = 0._dp
      ELSE
        ! The yield stress, equation (10) in Martin et al. (2011)
        ice%tau_yield(vi) = TAN(C%deg2rad * phi_fric(vi)) * (ice_density * grav - P_wat(vi))
      END IF
    END DO
    CALL sync
   
  END SUBROUTINE basal_yield_stress

END MODULE general_ice_model_data_module
