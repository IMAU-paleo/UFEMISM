MODULE BMB_module

  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync, &
                                         allocate_shared_memory_int_0D, allocate_shared_memory_dp_0D, &
                                         allocate_shared_memory_int_1D, allocate_shared_memory_dp_1D, &
                                         allocate_shared_memory_int_2D, allocate_shared_memory_dp_2D, &
                                         allocate_shared_memory_int_3D, allocate_shared_memory_dp_3D, &
                                         deallocate_shared_memory
  USE data_types_module,           ONLY: type_mesh, type_ice_model, type_climate_model, type_BMB_model, &
                                         type_remapping
  USE mesh_help_functions_module,  ONLY: FindVoronoiCellVertices
  USE forcing_module,              ONLY: forcing
  USE parameters_module,           ONLY: T0, L_fusion, seawater_density, ice_density
  USE mesh_mapping_module,           ONLY: Remap_field_dp, Remap_field_dp_3D, Remap_field_dp_monthly, &
                                           Reallocate_field_dp, Reallocate_field_dp_3D, Reallocate_field_int 

  IMPLICIT NONE
    
CONTAINS

  ! Run the SMB model on the region mesh
  SUBROUTINE run_BMB_model( mesh, ice, climate, region_time, BMB)
    ! Calculate basal melt
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_climate_model),            INTENT(IN)    :: climate
    REAL(dp),                            INTENT(IN)    :: region_time
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables
    INTEGER                                            :: vi
    
    REAL(dp)                                           :: weight                      ! weighting index for changes in ocean melt
    REAL(dp)                                           :: T_ocean                     ! temperature beneath the  shelves [Celcius]
    REAL(dp)                                           :: M_deep                      ! melt for deep-ocean areas [m/year]
    REAL(dp)                                           :: M_expo                      ! melt for exposed shelf [m/year]
    REAL(dp)                                           :: z_deep                      ! deep-ocean weighting
    REAL(dp)                                           :: z_expo                      ! exposed-shelf weighting
    
    REAL(dp)                                           :: T_freeze                    ! Freezing temperature at the base of the shelf (Celcius)
    REAL(dp)                                           :: S_flux                      ! melt flux from ocean to shelf
    
    REAL(dp), PARAMETER                                :: cp0        = 3974._dp       ! specific heat capacity of the ocean mixed layer (J kg-1 K-1) 
    REAL(dp), PARAMETER                                :: F_melt     = 5.0E-03_dp     ! Melt parameter (m s-1) (PISM-PIK original is 5E-03)
    REAL(dp), PARAMETER                                :: gamma_T    = 1.0E-04_dp     ! Thermal exchange velocity (m s-1)
    
    REAL(dp), PARAMETER                                :: T_ocean_CD = -5._dp         ! cold period temperature of the ocean beneath the shelves [Celcius]
    REAL(dp), PARAMETER                                :: T_ocean_PD = -1.7_dp        ! present day temperature of the ocean beneath the shelves [Celcius]
    REAL(dp), PARAMETER                                :: T_ocean_WM =  2._dp         ! warm period temperature of the ocean beneath the shelves [Celcius]
            
    REAL(dp), PARAMETER                                :: M_deep_CD =  2._dp          ! cold period value for deep-ocean areas [m/year]
    REAL(dp), PARAMETER                                :: M_deep_PD =  5._dp          ! present day value for deep-ocean areas [m/year]
    REAL(dp), PARAMETER                                :: M_deep_WM = 10._dp          ! warm period value for deep-ocean areas [m/year]    

    REAL(dp), PARAMETER                                :: M_expo_CD =  0._dp          ! cold period value for exposed areas [m/year]
    REAL(dp), PARAMETER                                :: M_expo_PD =  3._dp          ! present day value for exposed areas [m/year]
    REAL(dp), PARAMETER                                :: M_expo_WM =  6._dp          ! warm period value for exposed areas [m/year]
    
    ! Check if we need to apply any special benchmark experiment SMB
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler') THEN
          
        BMB%BMB( mesh%v1:mesh%v2) = 0._dp
        CALL sync
        RETURN
        
      ELSE
        WRITE(0,*) '  ERROR: "', C%choice_benchmark_experiment, '" is not a valid choice_benchmark_experiment!'
        STOP
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
      
    ! Initialise everything at zero
    BMB%BMB(              mesh%v1:mesh%v2) = 0._dp
    BMB%BMB_sheet(        mesh%v1:mesh%v2) = 0._dp
    BMB%BMB_shelf(        mesh%v1:mesh%v2) = 0._dp
    CALL sync  
    
    ! Find how much ocean water intrudes underneath the shelf
    CALL FindOceanIntrusion(mesh, ice, BMB)
    
!    weight = MAX(0._dp, MIN(2._dp, 1._dp + gT_offset/12._dp + MAX(0._dp, (insolation_data%Jun65N - 462.29_dp)/40._dp)))
    weight = 1._dp
    
    IF  ( weight >= 0._dp .AND. weight < 1._dp) THEN
      T_ocean = weight * T_ocean_PD + (1._dp - weight) * T_ocean_CD
      M_deep  = weight * M_deep_PD  + (1._dp - weight) * M_deep_CD
      M_expo  = weight * M_expo_PD  + (1._dp - weight) * M_expo_CD
    ELSE IF ( weight >= 1._dp .AND. weight <= 2._dp) THEN
      T_ocean = (2._dp - weight) * T_ocean_PD + (weight - 1._dp) * T_ocean_WM
      M_deep  = (2._dp - weight) * M_deep_PD  + (weight - 1._dp) * M_deep_WM
      M_expo  = (2._dp - weight) * M_expo_PD  + (weight - 1._dp) * M_expo_WM
    ELSE 
      T_ocean = T_ocean_PD
      M_deep  = M_deep_PD
      M_expo  = M_expo_PD
    END IF
    
    ! Calculate the two weighting factors for the deep-ocean and exposed shelves
    DO vi = mesh%v1, mesh%v2
      IF (.NOT. (ice%mask_shelf(vi)==1 .OR. ice%mask_ocean(vi)==1)) CYCLE
      
      ! Freezing temperature at the bottom of the ice shelves
      T_freeze = 0.0939_dp - 0.057_dp * 35._dp - 7.64E-04_dp * ice%Hi(vi) * ice_density / seawater_density

      ! Basal melting, mass flux, from shelf to ocean (Martin, TCD, 2010) - only melt values, when T_ocean > T_freeze.
      S_flux   = seawater_density * cp0 * C%sec_per_year * gamma_T * F_melt * & 
                      (T_ocean - T_freeze) / (L_fusion * ice_density)
      z_deep = MAX(0._dp, MIN(1._dp, (ice%sealevel(vi) - ice%Hb(vi) - 1200._dp)/200._dp))
      z_expo = MAX(0._dp, MIN(1._dp, BMB%ocean_intrusion(vi)))
        
      BMB%BMB_shelf(vi) = (z_deep - 1._dp) * ( (1._dp - z_expo) * S_flux + z_expo * M_expo) - z_deep * M_deep
      
    END DO
    CALL sync    
    
    ! Add sheet and shelf melt together
    BMB%BMB(mesh%v1:mesh%v2) = BMB%BMB_sheet(mesh%v1:mesh%v2) + BMB%BMB_shelf(mesh%v1:mesh%v2)
    CALL sync
          
  END SUBROUTINE run_BMB_model
  
  ! Determine the "ocean intrusion" coefficient
  SUBROUTINE FindOceanIntrusion( mesh, ice, BMB)
    ! Determine an "ocean intrusion" coefficient; 1 for open water, decreasing to 0
    ! as we go further underneath the shelf. Calculated using a simple diffusion scheme.
    ! Much faster than the subtended angle parameterisation.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    INTEGER,  PARAMETER                                :: n_it   = 25
    REAL(dp), PARAMETER                                :: c_diff = 0.15_dp
    INTEGER                                            :: it, vi, ci, vj
    REAL(dp)                                           :: dist, F
    REAL(dp), DIMENSION(mesh%nV)                       :: ocean_intrusion_old
        
    BMB%ocean_intrusion(mesh%v1:mesh%v2) = 0._dp
    DO vi = mesh%v1, mesh%v2
      IF (ice%mask_ocean(vi)==1) BMB%ocean_intrusion(vi) = 1._dp
    END DO
    CALL sync

    DO it = 1, n_it

      ocean_intrusion_old = BMB%ocean_intrusion
      CALL sync

      DO vi = mesh%v1, mesh%v2
        IF (ice%mask_shelf(vi)==0) CYCLE

        F = 0._dp

        DO ci = 1, mesh%nC(vi)
          vj = mesh%C(vi,ci)
          
          ! Only diffuse inside shelf and at shelf-ocean interface
          IF (.NOT. (ice%mask_shelf(vj)==1 .OR. ice%mask_ocean(vj)==1)) CYCLE

          dist = SQRT( (mesh%V(vi,1)-mesh%V(vj,1))**2 + (mesh%V(vi,2)-mesh%V(vj,2))**2);
          F = F - ((ocean_intrusion_old(vi) - ocean_intrusion_old(vj)) * mesh%Cw(vi,ci) * c_diff / dist)
        END DO

        BMB%ocean_intrusion(vi) = MAX(0._dp, MIN(1._dp, BMB%ocean_intrusion(vi) + F))

      END DO ! DO vi = mesh%v1, mesh%v2
      CALL sync
 
    END DO ! DO it = 1, nit
    
  END SUBROUTINE FindOceanIntrusion
  
  ! Administration: allocation, initialisation, and remapping
  SUBROUTINE initialise_BMB_model( mesh, BMB)
    ! Allocate memory for the data fields of the SMB model.
    
    IMPLICIT NONE
    
    !In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Allocate memory
    CALL allocate_shared_memory_dp_1D( mesh%nV,     BMB%BMB,              BMB%wBMB             )
    CALL allocate_shared_memory_dp_1D( mesh%nV,     BMB%BMB_sheet,        BMB%wBMB_sheet       )
    CALL allocate_shared_memory_dp_1D( mesh%nV,     BMB%BMB_shelf,        BMB%wBMB_shelf       )
    CALL allocate_shared_memory_dp_1D( mesh%nV,     BMB%sub_angle,        BMB%wsub_angle       )
    CALL allocate_shared_memory_dp_1D( mesh%nV,     BMB%dist_open,        BMB%wdist_open       )
    CALL allocate_shared_memory_dp_1D( mesh%nV,     BMB%ocean_intrusion,  BMB%wocean_intrusion )
      
    BMB%BMB(              mesh%v1:mesh%v2) = 0._dp
    BMB%BMB_sheet(        mesh%v1:mesh%v2) = 0._dp
    BMB%BMB_shelf(        mesh%v1:mesh%v2) = 0._dp
    BMB%sub_angle(        mesh%v1:mesh%v2) = 0._dp
    BMB%dist_open(        mesh%v1:mesh%v2) = 0._dp
    BMB%ocean_intrusion(  mesh%v1:mesh%v2) = 0._dp
      
  END SUBROUTINE initialise_BMB_model
  SUBROUTINE remap_BMB_model( mesh_old, mesh_new, map, BMB)
    ! Remap or reallocate all the data fields
  
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
        
    ! All BMB fields are remapped, since they are not updated in every timestep.
    CALL Remap_field_dp(         mesh_old, mesh_new, map, BMB%BMB,                 BMB%wBMB,                 'trilin'        )
    CALL Remap_field_dp(         mesh_old, mesh_new, map, BMB%BMB_sheet,           BMB%wBMB_sheet,           'trilin'        )
    CALL Remap_field_dp(         mesh_old, mesh_new, map, BMB%BMB_shelf,           BMB%wBMB_shelf,           'trilin'        )
    CALL Remap_field_dp(         mesh_old, mesh_new, map, BMB%sub_angle,           BMB%wsub_angle,           'trilin'        )
    CALL Remap_field_dp(         mesh_old, mesh_new, map, BMB%dist_open,           BMB%wdist_open,           'trilin'        )
    CALL Remap_field_dp(         mesh_old, mesh_new, map, BMB%ocean_intrusion,     BMB%wocean_intrusion,     'trilin'        ) 
    
  END SUBROUTINE remap_BMB_model
  
END MODULE BMB_module
