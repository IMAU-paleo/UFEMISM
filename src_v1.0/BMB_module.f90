MODULE BMB_module

  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync, &
                                         allocate_shared_memory_int_0D, allocate_shared_memory_dp_0D, &
                                         allocate_shared_memory_int_1D, allocate_shared_memory_dp_1D, &
                                         allocate_shared_memory_int_2D, allocate_shared_memory_dp_2D, &
                                         allocate_shared_memory_int_3D, allocate_shared_memory_dp_3D, &
                                         deallocate_shared_memory
  USE data_types_module,           ONLY: type_mesh, type_ice_model, type_climate_model, type_init_data_fields, type_BMB_model
  USE mesh_help_functions_module,  ONLY: FindVoronoiCellVertices
  USE forcing_module,              ONLY: forcing
  USE parameters_module,           ONLY: T0, L_fusion, seawater_density, ice_density 

  IMPLICIT NONE
  
  REAL(dp), PARAMETER :: albedo_water        = 0.1_dp
  REAL(dp), PARAMETER :: albedo_soil         = 0.2_dp
  REAL(dp), PARAMETER :: albedo_ice          = 0.5_dp
  REAL(dp), PARAMETER :: albedo_snow         = 0.85_dp
  REAL(dp), PARAMETER :: initial_snow_depth  = 0.1_dp
    
CONTAINS

  ! Run the SMB model on the region mesh
  SUBROUTINE run_BMB_model(mesh, ice, climate, region_time, BMB)
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
    
    ! If we're doing one of the EISMINT experiments, set basal mass balance to zero
    IF (C%do_eismint_experiment .OR. C%do_mismip3d_experiment) THEN
      BMB%BMB(mesh%v1:mesh%v2) = 0._dp
      CALL sync
      RETURN
    END IF
      
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
  
  SUBROUTINE run_BMB_model_old(mesh, ice, climate, region_time, BMB)
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
    
    ! If we're doing one of the EISMINT experiments, set basal mass balance to zero
    IF (C%do_eismint_experiment) THEN
      BMB%BMB(mesh%v1:mesh%v2) = 0._dp
      CALL sync
      RETURN
    END IF
      
    ! Initialise everything at zero
    BMB%BMB(              mesh%v1:mesh%v2) = 0._dp
    BMB%BMB_sheet(        mesh%v1:mesh%v2) = 0._dp
    BMB%BMB_shelf(        mesh%v1:mesh%v2) = 0._dp
    
    ! Find the subtended angles and distances to the open ocean for the BMB parameterisation
    CALL FindSubtendedAnglesAndDistances(mesh, ice, BMB)
    
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
      z_expo = MAX(0._dp, MIN(1._dp, (BMB%sub_angle(vi) - 80._dp)/30._dp)) * EXP(-BMB%dist_open(vi)/100000._dp)
        
      BMB%BMB_shelf(vi) = (z_deep - 1._dp) * ( (1._dp - z_expo) * S_flux + z_expo * M_expo) - z_deep * M_deep
      
    END DO
    CALL sync    
    
    ! Add sheet and shelf melt together
    BMB%BMB(mesh%v1:mesh%v2) = BMB%BMB_sheet(mesh%v1:mesh%v2) + BMB%BMB_shelf(mesh%v1:mesh%v2)
    
    ! DENK DROM
    !BMB%BMB(mesh%v1:mesh%v2) = BMB%sub_angle(mesh%v1:mesh%v2)
    !BMB%BMB(mesh%v1:mesh%v2) = BMB%dist_open(mesh%v1:mesh%v2)
          
  END SUBROUTINE run_BMB_model_old
  
  ! Initialise the BMB model (allocating shared memory)
  SUBROUTINE initialise_BMB_model(mesh, BMB)
    ! Allocate memory for the data fields of the SMB model.
    
    IMPLICIT NONE
    
    !In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables
    INTEGER                                            :: vi
    REAL(dp)                                           :: T_annual
    
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
  
  ! Finding the "subtended angle" and distance to open ocean for the sub-shelf melt parameterisation
  SUBROUTINE FindSubtendedAnglesAndDistances(mesh, ice, BMB)
    ! The idea: perform a flood-fill style search outward from each vertex, bounded by non-ocean
    ! (either land or marine ice sheet) and open ocean (defined as water deeper than 1500 m)
    ! List all the shoreline and open ocean border vertices we find. Then, make a list of all
    ! the surrounding angles (with a 1 degree resolution) where we find open ocean, and another
    ! one for where we find land. Also keep track of how close the land/ocean is for each angle.
    ! Integrate the angles where the ocean is closer than the land to find the subtended angle.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables
    INTEGER,  DIMENSION(mesh%nV)                       :: ocean_map,       land_map
    INTEGER,  DIMENSION(mesh%nV)                       :: ocean_stack,     land_stack
    INTEGER                                            :: ocean_stackN,    land_stackN
    INTEGER,  DIMENSION(mesh%nV)                       :: map, stack
    INTEGER                                            :: stackN
    INTEGER,  DIMENSION(360)                           :: ocean_angles,    land_angles
    REAL(dp), DIMENSION(360)                           :: ocean_distances, land_distances
    INTEGER                                            :: vi, ci, vc, ii, il, iu, vii, vio, v
    REAL(dp), DIMENSION(:,:), ALLOCATABLE              :: Vor
    INTEGER                                            :: nVor
    REAL(dp)                                           :: D, theta, thetamin, thetamax
    LOGICAL                                            :: is_shelf_or_shallow
    
    REAL(dp), PARAMETER                                :: R_search = 250._dp * 1000._dp
    REAL(dp), PARAMETER                                :: open_ocean_depth = 1200._dp
        
    ALLOCATE(Vor(mesh%nconmax+2,2))

    ! Initialise
    DO vi = mesh%v1, mesh%v2    
      IF (ice%mask_ocean(vi)==1) THEN
        BMB%sub_angle(vi) = 360._dp
        BMB%dist_open(vi) = 0._dp
      ELSEIF (ice%mask_shelf(vi)==1) THEN
        BMB%sub_angle(vi) = 360._dp
        BMB%dist_open(vi) = R_search
      ELSE
        BMB%sub_angle(vi) = 0._dp
        BMB%dist_open(vi) = R_search
      END IF
    END DO
    
    DO vi = mesh%v1, mesh%v2
        
      ! Only treat shelf vertices
      IF (ice%mask_shelf(vi)==0) CYCLE
              
      ! Clean maps and stacks
      ocean_map       = 0
      ocean_stackN    = 0      
      land_map        = 0
      land_stackN     = 0      
      map             = 0
      stackN          = 0      
      ocean_angles    = 0
      ocean_distances = (mesh%xmax-mesh%xmin) + (mesh%ymax-mesh%ymin)
      land_angles     = 0
      land_distances  = (mesh%xmax-mesh%xmin) + (mesh%ymax-mesh%ymin)
            
      ! == Search outward floodfill-style, find the borders of land and open ocean
      ! surrounding this shallow ocean vertex

      map(vi)  = 1
      stack(1) = vi
      stackN   = 1

      DO WHILE (stackN > 0)

        ! Get vertex index from the stack
        vii = stack(stackN)

        ! Remove this vertex from the stack
        stack(stackN) = 0
        stackN = stackN-1

        ! Only search within prescribed radius
        D = SQRT( (mesh%V(vi,1)-mesh%V(vii,1))**2 + (mesh%V(vi,2)-mesh%V(vii,2))**2)
        IF (D > R_search) CYCLE

        ! Add all shallow-ocean neighbours to the stack
        ! Add all ocean border and land border vertices to their respective maps and stacks
        DO ci = 1, mesh%nC(vii)
          vc = mesh%C(vii,ci)
          
          IF (ice%mask_shelf(vc)==1 .AND. map(vc) == 0) THEN
            map(vc)       = 1
            stackN        = stackN + 1
            stack(stackN) = vc
          ELSE
            IF ((ice%mask_land(vc)==1 .OR. ice%mask_sheet(vc)==1) .AND. land_map(vc) == 0) THEN
              land_map(vc)              = 1
              land_stackN               = land_stackN+1
              land_stack(land_stackN)   = vc;
            ELSEIF (ice%mask_ocean(vc)==1 .AND. ocean_map(vc) == 0) THEN
              ocean_map(vc)             = 1
              ocean_stackN              = ocean_stackN+1
              ocean_stack(ocean_stackN) = vc;
            END IF
          END IF
          
        END DO ! DO ci = 1, mesh%nC(vii)

      END DO !  DO WHILE (stackN > 0)

      ! == For every ocean border vertex, add its angles to the list
      DO vii = 1, ocean_stackN
      
        vio = ocean_stack(vii)
        CALL FindVoronoiCellVertices(mesh, vio, Vor, nVor)

        thetamin =  400._dp
        thetamax = -400._dp

        DO v = 1, nVor
          theta = ATAN2( Vor(v,2) - mesh%V(vi,2), Vor(v,1) - mesh%V(vi,1)) * C%rad2deg
          thetamin = MIN(thetamin, theta)
          thetamax = MAX(thetamax, theta)
          IF (thetamax-thetamin > 180._dp) thetamax = thetamax - 360._dp
          il = MAX(1,   FLOOR(   thetamin+180._dp))
          iu = MIN(360, CEILING( thetamax+180._dp))
          D = SQRT( (mesh%V(vi,1)-Vor(v,1))**2 + (mesh%V(vi,2)-Vor(v,2))**2)
          DO ii = il, iu
            ocean_angles(   ii) = 1
            ocean_distances(ii) = MIN(ocean_distances(ii),D)
          END DO
        END DO
        
      END DO ! DO vii = 1, ocean_stackN

      ! == For every land border vertex, add its angles to the list
      DO vii = 1, land_stackN
      
        vio = land_stack(vii)
        CALL FindVoronoiCellVertices(mesh, vio, Vor, nVor)

        thetamin =  400._dp
        thetamax = -400._dp

        DO v = 1, nVor
          theta = ATAN2( Vor(v,2) - mesh%V(vi,2), Vor(v,1) - mesh%V(vi,1)) * C%rad2deg
          thetamin = MIN(thetamin, theta)
          thetamax = MAX(thetamax, theta)
          IF (thetamax-thetamin > 180._dp) thetamax = thetamax - 360._dp
          il = MAX(1,   FLOOR(   thetamin+180._dp))
          iu = MIN(360, CEILING( thetamax+180._dp))
          D = SQRT( (mesh%V(vi,1)-Vor(v,1))**2 + (mesh%V(vi,2)-Vor(v,2))**2)
          DO ii = il, iu
            land_angles(   ii) = 1
            land_distances(ii) = MIN(land_distances(ii),D)
          END DO
        END DO
        
      END DO ! DO vii = 1, land_stackN

      ! == Find number of angle bins where the open ocean is closer than the coast      
      DO ii = 1, 360
        IF (land_angles(ii)==1 .AND. land_distances(ii) < ocean_distances(ii)) THEN
          BMB%sub_angle(vi) = BMB%sub_angle(vi) - 1._dp
        END IF
        BMB%dist_open(vi) = MIN( BMB%dist_open(vi), ocean_distances(ii))
      END DO
    
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync  
    
    DEALLOCATE(Vor)
    
  END SUBROUTINE FindSubtendedAnglesAndDistances  
  ! Determine the "ocean intrusion" coefficient
  SUBROUTINE FindOceanIntrusion(mesh, ice, BMB)
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
  
END MODULE BMB_module
