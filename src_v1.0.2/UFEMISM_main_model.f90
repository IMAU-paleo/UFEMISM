MODULE UFEMISM_main_model

  USE mpi
  USE parallel_module,               ONLY: par, sync, &
                                           allocate_shared_int_0D, allocate_shared_dp_0D, &
                                           allocate_shared_int_1D, allocate_shared_dp_1D, &
                                           allocate_shared_int_2D, allocate_shared_dp_2D, &
                                           allocate_shared_int_3D, allocate_shared_dp_3D, &
                                           deallocate_shared
  USE data_types_module,             ONLY: type_model_region, type_mesh, type_ice_model, type_PD_data_fields, type_init_data_fields, &
                                           type_climate_model, type_climate_matrix, type_SMB_model, type_BMB_model, type_remapping
  USE configuration_module,          ONLY: dp, C  
  USE parameters_module,             ONLY: seawater_density, ice_density, T0
  USE reference_fields_module,       ONLY: initialise_PD_data_fields, initialise_init_data_fields_cart, deallocate_init_data, &
                                           initialise_init_data_fields_mesh, map_PD_data_to_mesh, map_init_data_to_mesh
  USE mesh_memory_module,            ONLY: deallocate_mesh_all, allocate_mesh_secondary
  USE mesh_help_functions_module,    ONLY: get_lat_lon_coordinates, find_Voronoi_cell_areas, find_triangle_areas, find_connection_widths, &
                                           determine_mesh_resolution
  USE mesh_creation_module,          ONLY: create_mesh_from_cart_data
  USE mesh_mapping_module,           ONLY: create_remapping_arrays, deallocate_remapping_arrays, reallocate_field_dp, reallocate_field_int
  USE mesh_ArakawaC_module,          ONLY: make_Ac_mesh
  USE mesh_derivatives_module,       ONLY: get_neighbour_functions
  USE mesh_update_module,            ONLY: determine_mesh_fitness, create_new_mesh
  USE netcdf_module,                 ONLY: create_output_file, write_to_output_file
  USE general_ice_model_data_module, ONLY: ice_physical_properties, update_general_ice_model_data
  USE ice_dynamics_module,           ONLY: initialise_ice_model, remap_ice_model, solve_SIA, solve_SSA, calculate_ice_thickness_change
  USE thermodynamics_module,         ONLY: initialise_ice_temperature, update_ice_temperature
  USE climate_module,                ONLY: initialise_climate_model, remap_climate_model, run_climate_model
  USE SMB_module,                    ONLY: initialise_SMB_model,     remap_SMB_model,     run_SMB_model
  USE BMB_module,                    ONLY: initialise_BMB_model,     remap_BMB_model,     run_BMB_model

  IMPLICIT NONE

CONTAINS

  SUBROUTINE run_model( region, matrix, t_end)
    ! Run the model until t_end (usually a 100 years further)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    TYPE(type_climate_matrix),  INTENT(IN)        :: matrix
    REAL(dp),                   INTENT(IN)        :: t_end
    
    ! Local variables
    REAL(dp)                                      :: meshfitness
    REAL(dp)                                      :: t1, t2
    
    IF (par%master) WRITE (0,'(A,A,A,A,A,F9.3,A,F9.3,A)') '  Running model region ', region%name, ' (', TRIM(region%long_name), & 
                                                          ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' ky'
               
    ! Computation time tracking
    region%tcomp_total          = 0._dp
    region%tcomp_mesh           = 0._dp
    region%tcomp_SIA            = 0._dp
    region%tcomp_SSA            = 0._dp
    region%tcomp_thermo         = 0._dp
    region%tcomp_climate        = 0._dp
    region%tcomp_output         = 0._dp
       
    t1 = MPI_WTIME()
    t2 = 0._dp
    
    ! Write to text output at t=0
    IF (par%master .AND. region%time == C%start_time_of_run) THEN
      CALL write_text_output(region)
    END IF
                            
    ! The main model time loop
    ! ========================
    
    DO WHILE (region%time < t_end)
          
      ! Update ice thickness based on ice velocity and mass balance
      IF (par%master) t2 = MPI_WTIME()
      CALL calculate_ice_thickness_change( region%mesh, region%ice, region%SMB, region%BMB, region%dt)
      IF (par%master) region%tcomp_SIA = region%tcomp_SIA + MPI_WTIME() - t2
      
      ! Check if the mesh needs to be updated
      IF (par%master) t2 = MPI_WTIME()
      meshfitness = 1._dp
      IF (region%time > region%t1_mesh) THEN
        CALL determine_mesh_fitness(region%mesh, region%ice, meshfitness)        
        ! Update the time when mesh fitness was last checked
        region%t0_mesh = region%time
        region%t1_mesh = region%time + C%dt_mesh_min
      END IF
      IF (par%master) region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2 
    
      ! If required, update the mesh
      IF (meshfitness < C%mesh_fitness_threshold) THEN
    !  IF (.FALSE.) THEN
    !  IF (.TRUE.) THEN    
    
        IF (par%master) t2 = MPI_WTIME()
    
        ! Create a new mesh
        CALL create_new_mesh( region)      
        
        ! Remap model data from the old mesh to the new one
        CALL remap_model( region%mesh, region%mesh_new, region%PD, region%ice, matrix, region%climate, region%SMB, region%BMB)
        
        ! Deallocate the old mesh
        CALL deallocate_mesh_all( region%mesh)
        region%mesh = region%mesh_new
        
        IF (par%master) region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2
        
        ! A new mesh requires a new output file
        region%output_file_exists = .FALSE.
        
        ! Run all model components again after updating the mesh
        region%do_solve_SIA      = .TRUE.
        region%do_solve_SSA      = .TRUE.
        region%do_climate        = .TRUE.
        region%do_SMB            = .TRUE.
        region%do_BMB            = .TRUE.
        region%do_thermodynamics = .TRUE.
      END IF
      
      ! Update general ice model data - Hs, masks, ice physical properties, on both Aa and Ac mesh
      IF (par%master) t2 = MPI_WTIME()
      CALL update_general_ice_model_data( region%mesh, region%ice, region%time)
      IF (par%master) region%tcomp_SIA = region%tcomp_SIA + MPI_WTIME() - t2
      
    ! Ice dynamics
    ! ============
      
      ! Solve the SIA
      IF (region%do_solve_SIA) THEN
        IF (par%master) t2 = MPI_WTIME()
        CALL solve_SIA( region%mesh, region%ice)
        IF (par%master) region%tcomp_SIA = region%tcomp_SIA + MPI_WTIME() - t2
        region%t0_SIA = region%time
      END IF
      
      ! Solve the SSA
      IF (region%do_solve_SSA) THEN
        IF (par%master) t2 = MPI_WTIME()
        CALL solve_SSA( region%mesh, region%ice)
        IF (par%master) region%tcomp_SSA = region%tcomp_SSA + MPI_WTIME() - t2
        region%t0_SSA = region%time
      END IF
      
    ! Climate , SMB and BMB
    ! =====================
    
      IF (par%master) t2 = MPI_WTIME()
            
      ! Run the climate model
      IF (region%do_climate) THEN
        CALL run_climate_model( region%mesh, region%ice, region%climate, region%time)
        region%t0_climate = region%time
      END IF
    
      ! Run the SMB model
      IF (region%do_SMB) THEN
        CALL run_SMB_model( region%mesh, region%ice, region%climate, region%time, region%SMB)
        region%t0_SMB = region%time
      END IF
    
      ! Run the BMB model
      IF (region%do_BMB) THEN
        CALL run_BMB_model( region%mesh, region%ice, region%BMB)
        region%t0_BMB = region%time
      END IF
      
      IF (par%master) region%tcomp_climate = region%tcomp_climate + MPI_WTIME() - t2
      
    ! Thermodynamics
    ! ==============
      
      ! Solve thermodynamics
      IF (region%do_thermodynamics) THEN
        IF (par%master) t2 = MPI_WTIME()
        CALL update_ice_temperature(region%mesh, region%ice, region%climate, region%SMB)
        IF (par%master) region%tcomp_thermo = region%tcomp_thermo + MPI_WTIME() - t2
        region%t0_thermo = region%time
      END IF
      
    ! Time step and output
    ! ====================
                            
      ! Write output
      IF (region%do_write_output) THEN
        IF (par%master) t2 = MPI_WTIME()
        ! If the mesh has been updated, create a new NetCDF file
        IF (.NOT. region%output_file_exists) THEN
          IF (par%master) CALL create_output_file(region)
          CALL sync
          region%output_file_exists = .TRUE.
        END IF  
        IF (par%master) CALL write_to_output_file(region)
        region%t0_output = region%time
        IF (par%master) region%tcomp_output = region%tcomp_output + MPI_WTIME() - t2
      END IF
      
      ! Calculate time steps, determine actions and advance region time
      IF (par%master) t2 = MPI_WTIME()
      CALL determine_timesteps_and_actions(region, t_end)
      IF (par%master) region%tcomp_SIA = region%tcomp_SIA + MPI_WTIME() - t2
      
      ! DENK DROM
     ! region%time = t_end
    
    END DO ! DO WHILE (region%time < t_end)
    
    ! Write to NetCDF output one last time at the end of the simulation
    IF (region%time == C%end_time_of_run) THEN
      ! If the mesh has been updated, create a new NetCDF file
      IF (.NOT. region%output_file_exists) THEN
        IF (par%master) CALL create_output_file(region)
        CALL sync
        region%output_file_exists = .TRUE.
      END IF      
      IF (par%master) CALL write_to_output_file(region)
    END IF 
    
    ! Determine total ice sheet area, volume, volume-above-flotation and GMSL contribution,
    ! used for writing to text output and in the inverse routine
    CALL calculate_icesheet_volume_and_area(region)
    
    ! Write to text output
    IF (par%master) THEN
      region%tcomp_total = MPI_WTIME() - t1
      CALL write_text_output(region)
    END IF
    
    IF (par%master) WRITE(0,*) ''
    
  END SUBROUTINE run_model
  
  ! Map data from an old to a new mesh
  SUBROUTINE remap_model( mesh_old, mesh_new, PD, ice, matrix, climate, SMB, BMB)
    ! Map all the model data from the old mesh to the new mesh.
    
    USE climate_module,              ONLY: map_subclimate_to_mesh
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_old
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_new
    TYPE(type_PD_data_fields),           INTENT(INOUT) :: PD
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_climate_matrix),           INTENT(IN)    :: matrix
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables
    TYPE(type_remapping)                               :: map
    
    IF (par%master) WRITE(0,*) '  Mapping model data to the new mesh...'
    
! Remap essential data
! ====================
    
    ! Calculate the mapping arrays
    CALL create_remapping_arrays( mesh_old, mesh_new, map)
    
    ! Remap all the submodels
    CALL remap_ice_model(     mesh_old, mesh_new, map, ice)
    CALL remap_climate_model( mesh_old, mesh_new, map, climate, matrix)
    CALL remap_SMB_model(     mesh_old, mesh_new, map, SMB)
    CALL remap_BMB_model(     mesh_old, mesh_new, map, BMB)
            
    ! Reallocate and remap reference PD data    
    CALL reallocate_field_dp(  mesh_new%nV, PD%Hi,   PD%wHi  )
    CALL reallocate_field_dp(  mesh_new%nV, PD%Hb,   PD%wHb  )
    CALL reallocate_field_dp(  mesh_new%nV, PD%Hs,   PD%wHs  )
    CALL reallocate_field_int( mesh_new%nV, PD%mask, PD%wmask)
    CALL map_PD_data_to_mesh(mesh_new, PD)
    
    ! Hb is more accurate when mapped from the PD reference field
    ice%Hb(mesh_new%v1:mesh_new%v2) = PD%Hb(mesh_new%v1:mesh_new%v2)
    
    ! Redefine Hs
    ice%Hs(mesh_new%v1:mesh_new%v2) = ice%Hb(mesh_new%v1:mesh_new%v2) + ice%Hi(mesh_new%v1:mesh_new%v2)
            
    ! Deallocate shared memory for the mapping arrays    
    CALL deallocate_remapping_arrays(map)
      
  END SUBROUTINE remap_model
  
  ! Calculate critical time steps from the SIA and SSA, and determine which actions should be taken
  SUBROUTINE determine_timesteps_and_actions( region, t_end)
    ! This subroutine calculates the SIA time step with the CFK-criterium
    ! Methods are similar to those in ANICE, which are based on PISM (see also Bueler  et al., 2007)
    
    ! It also determines whether the SIA velocities, SSA velocities or thermodynamics should be updated,
    ! or whether results should be written to the output file
    
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end
    
    ! Local variables:
    INTEGER                                     :: ci, vi, vj, k, ierr
    REAL(dp)                                    :: dist
    REAL(dp)                                    :: dt_D_2D,     dt_D_2D_min
    REAL(dp)                                    :: dt_V_3D_SIA, dt_V_3D_SIA_min
    REAL(dp)                                    :: dt_V_2D_SSA, dt_V_2D_SSA_min
    REAL(dp)                                    :: t_next_action
    
    REAL(dp), PARAMETER                         :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.
    
    ! Calculate the critical time steps resulting from the SIA and SSA velocities
    ! ===========================================================================

    dt_D_2D_min     = 1000._dp
    dt_V_3D_SIA_min = 1000._dp
    dt_V_2D_SSA_min = 1000._dp
    
    dt_D_2D         = 0._dp
    dt_V_3D_SIA     = 0._dp
    dt_V_2D_SSA     = 0._dp

    ! As in PISM we check the three 'domains' of diffusivity: the SIA, the grounded 3-D SIA,
    ! and the shelf velocities (SSA). Equations are based on Bueler et al. (JoG, 2007)
    
    DO ci = region%mesh%ac1, region%mesh%ac2      
      vi = region%mesh%Aci(ci,1)
      vj = region%mesh%Aci(ci,2)
      dist = SQRT((region%mesh%V(vj,1)-region%mesh%V(vi,1))**2 + (region%mesh%V(vj,2)-region%mesh%V(vi,2))**2)
      dt_D_2D = dist**2 / (-6._dp * C%pi * (region%ice%D_SIA_ac(ci) - 1E-09))
      dt_D_2D_min = MIN(dt_D_2D, dt_D_2D_min) 
      
      dt_V_2D_SSA = dist / (ABS(region%ice%U_SSA(vi)) + ABS(region%ice%V_SSA(vi)))
      dt_V_2D_SSA_min = MIN(dt_V_2D_SSA, dt_V_2D_SSA_min)   
      dt_V_2D_SSA = dist / (ABS(region%ice%U_SSA(vj)) + ABS(region%ice%V_SSA(vj)))
      dt_V_2D_SSA_min = MIN(dt_V_2D_SSA, dt_V_2D_SSA_min)  
    END DO
    
    DO vi = region%mesh%v1, region%mesh%v2
      dt_V_2D_SSA = SQRT(region%mesh%A(vi)/C%pi) / (ABS(region%ice%U_SSA(vi)) + ABS(region%ice%V_SSA(vi)))
      dt_V_2D_SSA_min = MIN(dt_V_2D_SSA, dt_V_2D_SSA_min)
      
      DO k = 1, C%nZ
        dt_V_3D_SIA = SQRT(region%mesh%A(vi)/C%pi) / (ABS(region%ice%U_3D(vi,k)) + ABS(region%ice%V_3D(vi,k)))
        dt_V_3D_SIA_min = MIN(dt_V_3D_SIA, dt_V_3D_SIA_min)
      END DO
    END DO
     
    ! Collect smallest values across all processes    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_D_2D_min,     1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_V_2D_SSA_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_V_3D_SIA_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
     
    ! Correction factor. Using the exact result of the CFK criterion result in some small not-quite-unstable oscillations
    ! in the solution; better to take a slightly larger time step to prevent this.
    dt_D_2D_min     = dt_D_2D_min     * dt_correction_factor
    dt_V_2D_SSA_min = dt_V_2D_SSA_min * dt_correction_factor
    dt_V_3D_SIA_min = dt_V_3D_SIA_min * dt_correction_factor
    
    ! Calculate ice dynamics critical time step - only for writing to screen, the actual time step can be smaller
    ! if we need to do thermodynamics or write to output...
    region%dt = MIN( MIN( MIN( dt_D_2D_min, dt_V_2D_SSA_min), dt_V_3D_SIA_min), C%dt_max)      
    IF (ABS(1._dp - region%dt / region%dt_prev) > 0.1_dp) THEN
      region%dt_prev = region%dt
      IF (par%master) WRITE(0,'(A,F7.4,A)') '   Critical time step: ', region%dt, ' yr'
    END IF

    ! Determine which action should be taken next
    ! ===========================================
    
    region%dt_SIA = MIN(C%dt_max, MIN(dt_D_2D_min, dt_V_3D_SIA_min))
    region%dt_SSA = MIN(C%dt_max, dt_V_2D_SSA_min)
    
    region%t1_SIA      = region%t0_SIA      + region%dt_SIA
    region%t1_SSA      = region%t0_SSA      + region%dt_SSA
    region%t1_thermo   = region%t0_thermo   + C%dt_thermo
    region%t1_climate  = region%t0_climate  + C%dt_climate
    region%t1_SMB      = region%t0_SMB      + C%dt_SMB
    region%t1_BMB      = region%t0_BMB      + C%dt_BMB
    region%t1_output   = region%t0_output   + C%dt_output
    
    t_next_action = MINVAL( [region%t1_SIA, region%t1_SSA, region%t1_thermo, region%t1_climate, region%t1_SMB, region%t1_BMB, region%t1_output])
    
    region%dt = t_next_action - region%time
    
    region%do_solve_SIA       = .FALSE.
    region%do_solve_SSA       = .FALSE.
    region%do_thermodynamics  = .FALSE.
    region%do_climate         = .FALSE.
    region%do_SMB             = .FALSE.
    region%do_BMB             = .FALSE.
    region%do_write_output    = .FALSE.
    
    IF (t_next_action == region%t1_SIA     ) region%do_solve_SIA      = .TRUE.
    IF (t_next_action == region%t1_SSA     ) region%do_solve_SSA      = .TRUE.
    IF (t_next_action == region%t1_thermo  ) region%do_thermodynamics = .TRUE.
    IF (t_next_action == region%t1_climate ) region%do_climate        = .TRUE.
    IF (t_next_action == region%t1_SMB     ) region%do_SMB            = .TRUE.
    IF (t_next_action == region%t1_BMB     ) region%do_BMB            = .TRUE.
    IF (t_next_action == region%t1_output  ) region%do_write_output   = .TRUE.
    
    ! If the next action occurs after the coupling interval, crop the time step and set all "do" logicals to TRUE
    ! (except for writing output, that one really should happen only at the prescribed interval)
    IF (t_next_action >= t_end) THEN
      region%dt = t_end - region%time
      region%do_solve_SIA       = .TRUE.
      region%do_solve_SSA       = .TRUE.
      region%do_thermodynamics  = .TRUE.   
      region%do_climate         = .TRUE.
      region%do_SMB             = .TRUE.
      region%do_BMB             = .TRUE.  
    END IF
            
    ! Advance region time
    region%time = region%time + region%dt
    
    CALL sync
        
  END SUBROUTINE determine_timesteps_and_actions
  
  ! Calculate this region's ice sheet's volume and area
  SUBROUTINE calculate_icesheet_volume_and_area( region)
    
    USE parameters_module,           ONLY: ocean_area, seawater_density, ice_density
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    INTEGER                                       :: vi, p
    REAL(dp)                                      :: ice_area
    REAL(dp)                                      :: ice_volume
    REAL(dp)                                      :: ice_volume_above_flotation
    
    ice_area                   = 0._dp
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp
    
    IF (par%master) THEN
      region%ice_area                   = 0._dp
      region%ice_volume                 = 0._dp
      region%ice_volume_above_flotation = 0._dp
    END IF
    CALL sync
    
    ! Calculate ice area and volume for processor domain
    DO vi = region%mesh%v1, region%mesh%v2
      IF (region%ice%mask_ice( vi) == 1) THEN
        ice_volume = ice_volume + (region%ice%Hi( vi) * region%mesh%A( vi) * seawater_density / (ice_density * ocean_area))
        ice_area   = ice_area   + (region%mesh%A( vi)) * 1.0E-06_dp
        
        IF (region%ice%Hi( vi) > (region%ice%SL( vi) - region%ice%Hb( vi)) * ice_density / seawater_density) THEN
          ice_volume_above_flotation = ice_volume_above_flotation + &
            ((region%ice%Hi( vi) - (MAX(0._dp, region%ice%SL( vi) - region%ice%Hb( vi)) * ice_density / seawater_density)) &
              * region%mesh%A( vi) * seawater_density / (ice_density * ocean_area))
        END IF
      END IF     
    END DO ! DO vi = region%mesh%v1, region%mesh%v2
    
    ! Add integrated values from the different processes
    DO p = 0, par%n-1
      IF (p == par%i) THEN
        region%ice_area                   = region%ice_area                   + ice_area
        region%ice_volume                 = region%ice_volume                 + ice_volume
        region%ice_volume_above_flotation = region%ice_volume_above_flotation + ice_volume_above_flotation
      END IF
      CALL sync
    END DO
    
    ! Calculate GMSL contribution
    IF (par%master) THEN
      region%GMSL_contribution = region%ice_volume_above_flotation - region%ice_volume_above_flotation_PD
    END IF
    CALL sync
    
  END SUBROUTINE calculate_icesheet_volume_and_area
  SUBROUTINE calculate_PD_sealevel_contribution( region)
    
    USE parameters_module,           ONLY: ocean_area, seawater_density, ice_density
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    INTEGER                                       :: i, j, p
    REAL(dp)                                      :: dx, dy
    REAL(dp)                                      :: ice_volume_above_flotation
    
    dx = ABS(region%PD%x(2) - region%PD%x(1))
    dy = ABS(region%PD%y(2) - region%PD%y(1))
    
    ice_volume_above_flotation = 0._dp
    
    IF (par%master) THEN
      region%ice_volume_above_flotation_PD = 0._dp
    END IF
    CALL sync
    
    ! Calculate ice volume for processor domain
    DO i = region%PD%i1, region%PD%i2
    DO j = 1, region%PD%ny
        
      IF (region%PD%Hi_cart(i,j) > (0._dp - region%PD%Hb_cart(i,j)) * ice_density / seawater_density) THEN
        ice_volume_above_flotation = ice_volume_above_flotation + &
          ((region%PD%Hi_cart(i,j) - (MAX(0._dp, 0._dp - region%PD%Hb_cart(i,j)) * ice_density / seawater_density)) &
            * dx * dy * seawater_density / (ice_density * ocean_area))
      END IF
      
    END DO
    END DO ! DO i = region%PD%i1, region%PD%i2
    
    ! Add integrated values from the different processes
    DO p = 0, par%n-1
      IF (p == par%i) THEN
        region%ice_volume_above_flotation_PD = region%ice_volume_above_flotation_PD + ice_volume_above_flotation
      END IF
      CALL sync
    END DO
    
  END SUBROUTINE calculate_PD_sealevel_contribution
  
  ! Create and write to this region's text output files - both the region-wide one,
  ! and the ones for the different Points-Of-Interest
  SUBROUTINE create_text_output_files( region)
    ! Creates the following text output files:
    !   time_log_REG.txt             - a log of how much computation time the different model parts take
    !   general_output_REG.txt       - some general info - ice sheet volume, average surface temperature, total mass balance, etc.
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(IN)        :: region
    
    CHARACTER(LEN=256)                            :: filename
    CHARACTER(LEN=3)                              :: ns
    CHARACTER(LEN=1)                              :: ns3
    CHARACTER(LEN=2)                              :: ns2
    INTEGER                                       :: n, k
        
    ! Time log
    ! ========
    
    filename = TRIM(C%output_dir) // 'aa_time_log_' // region%name // '.txt'
    OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
    
    WRITE(UNIT = 1337, FMT = '(A)') 'Time log for region ' // TRIM(region%long_name)
    WRITE(UNIT = 1337, FMT = '(A)') 'Computation time (in seconds) required by each model component'
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A)') '     Time  vertices     total       mesh        SIA        SSA     thermo    climate     output'
    
    CLOSE(UNIT = 1337)
        
    ! General output
    ! ==============
    
    filename = TRIM(C%output_dir) // 'aa_general_output_' // region%name // '.txt'
    OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
    
    WRITE(UNIT = 1337, FMT = '(A)') 'General output for region ' // TRIM(region%long_name)
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A)') ' Columns in order:'
    WRITE(UNIT = 1337, FMT = '(A)') '   1)  Model time                  (years) '
    WRITE(UNIT = 1337, FMT = '(A)') '   2)  Ice volume                  (meter sea level equivalent)'
    WRITE(UNIT = 1337, FMT = '(A)') '   3)  Ice volume above flotation  (meter sea level equivalent)'
    WRITE(UNIT = 1337, FMT = '(A)') '   4)  Ice area                    (km^2)'
    WRITE(UNIT = 1337, FMT = '(A)') '   5)  Mean surface temperature    (Kelvin)'
    WRITE(UNIT = 1337, FMT = '(A)') '   6)  Total snowfall     over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '   7)  Total rainfall     over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '   8)  Total melt         over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '   9)  Total refreezing   over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '  10)  Total runoff       over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '  11)  Total SMB          over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '  12)  Total BMB          over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '  13)  Total mass balance over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '  14)  Grounding line x position (for MISMIP benchmark experiments)'
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A)') '     Time     Ice  Ice-af     Ice-area     T2m       Snow       Rain       Melt   Refreeze     Runoff        SMB        BMB         MB       x_GL'
    
    CLOSE(UNIT = 1337)
    
    ! Point-of-interest output
    ! ========================
    
    DO n = 1, region%mesh%nPOI
    
      IF (n<10) THEN
        WRITE(ns3,'(I1)') n
        ns(1:2) = '00'
        ns(3:3) = ns3
      ELSEIF (n<100) THEN
        WRITE(ns2,'(I2)') n
        ns(1:1) = '0'
        ns(2:3) = ns3
      ELSE
        WRITE(ns,'(I3)') n
      END IF
    
      filename = TRIM(C%output_dir) // 'ab_POI_' // TRIM(ns) // '_data.txt'
      OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
    
      WRITE(UNIT = 1337, FMT = '(A)') 'Relevant data for Point of Interest ' // TRIM(ns)
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A,F10.2)') '  Lat = ', region%mesh%POI_coordinates(n,1)
      WRITE(UNIT = 1337, FMT = '(A,F10.2)') '  Lon = ', region%mesh%POI_coordinates(n,2)
      WRITE(UNIT = 1337, FMT = '(A,F10.2)') '  x   = ', region%mesh%POI_XY_coordinates(n,1)
      WRITE(UNIT = 1337, FMT = '(A,F10.2)') '  y   = ', region%mesh%POI_XY_coordinates(n,2)
      WRITE(UNIT = 1337, FMT = '(A,F10.2)') '  res = ', region%mesh%POI_resolutions(n)
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)') '  zeta = '
      DO k = 1, C%nZ
        WRITE(UNIT = 1337, FMT = '(F22.16)') C%zeta(k)
      END DO
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)') '     Time         Hi         Hb         Hs         Ti'
      
      CLOSE(UNIT = 1337)
    
    END DO
    
  END SUBROUTINE create_text_output_files
  SUBROUTINE write_text_output( region)
    ! Write data to the following text output files:
    !   time_log_REG.txt             - a log of how much computation time the different model parts take
    !   general_output_REG.txt       - some general info - ice sheet volume, average surface temperature, total mass balance, etc.
    
    USE parameters_module,           ONLY: ocean_area, seawater_density, ice_density
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(IN)        :: region
    
    CHARACTER(LEN=256)                            :: filename
    CHARACTER(LEN=3)                              :: ns
    CHARACTER(LEN=1)                              :: ns3
    CHARACTER(LEN=2)                              :: ns2
    INTEGER                                       :: n
    INTEGER                                       :: vi, m, k, aci, vj
    REAL(dp)                                      :: T2m_mean
    REAL(dp)                                      :: total_snowfall
    REAL(dp)                                      :: total_rainfall
    REAL(dp)                                      :: total_melt
    REAL(dp)                                      :: total_refreezing
    REAL(dp)                                      :: total_runoff
    REAL(dp)                                      :: total_SMB
    REAL(dp)                                      :: total_BMB
    REAL(dp)                                      :: total_MB
    REAL(dp)                                      :: TAFi, TAFj, lambda, x_GL
    INTEGER                                       :: n_GL
    
    INTEGER                                       :: vi1, vi2, vi3
    REAL(dp)                                      :: w1, w2, w3
    REAL(dp)                                      :: Hi_POI, Hb_POI, Hs_POI
    REAL(dp), DIMENSION(C%nZ)                     :: Ti_POI
        
  ! Time log
  ! ========
    
    filename = TRIM(C%output_dir) // 'aa_time_log_' // region%name // '.txt'
    OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
    
    WRITE(UNIT = 1337, FMT = '(F10.1,I9,F11.3,F11.3,F11.3,F11.3,F11.3,F11.3,F11.3,F11.3)') region%time, region%mesh%nV, &
      region%tcomp_total, region%tcomp_mesh, region%tcomp_SIA, region%tcomp_SSA, region%tcomp_thermo, region%tcomp_climate, region%tcomp_output
    
    CLOSE(UNIT = 1337)
    
  ! General output
  ! ==============
    
    T2m_mean                   = 0._dp
    total_snowfall             = 0._dp
    total_rainfall             = 0._dp
    total_melt                 = 0._dp
    total_refreezing           = 0._dp
    total_runoff               = 0._dp
    total_SMB                  = 0._dp
    total_BMB                  = 0._dp
    total_MB                   = 0._dp
    
    DO vi = 1, region%mesh%nV
      IF (region%ice%Hi(vi) > 0._dp) THEN
        
        total_BMB = total_BMB + (region%BMB%BMB(vi) * region%mesh%A(vi) / 1E9_dp)
          
        DO m = 1, 12
          total_snowfall   = total_snowfall   + (region%SMB%Snowfall(  vi,m) * region%mesh%A(vi) / 1E9_dp)
          total_rainfall   = total_rainfall   + (region%SMB%Rainfall(  vi,m) * region%mesh%A(vi) / 1E9_dp)
          total_melt       = total_melt       + (region%SMB%Melt(      vi,m) * region%mesh%A(vi) / 1E9_dp)
          total_refreezing = total_refreezing + (region%SMB%Refreezing(vi,m) * region%mesh%A(vi) / 1E9_dp)
          total_runoff     = total_runoff     + (region%SMB%Runoff(    vi,m) * region%mesh%A(vi) / 1E9_dp)
          total_SMB        = total_SMB        + (region%SMB%SMB(       vi,m) * region%mesh%A(vi) / 1E9_dp)
        END DO
        
      END IF

      T2m_mean = T2m_mean + SUM(region%climate%applied%T2m(vi,:)) * region%mesh%A(vi) / (12._dp * (region%mesh%xmax - region%mesh%xmin) * (region%mesh%ymax - region%mesh%ymin))
      
    END DO
    
    total_MB = total_SMB + total_BMB
    
    ! Average x-position of grounding line
    x_GL = 0._dp
    n_GL = 0
    DO aci = 1, region%mesh%nAc
      IF (region%ice%mask_gl_Ac( aci) == 1) THEN
        n_GL = n_GL + 1
        
        vi = region%mesh%Aci( aci,1)
        vj = region%mesh%Aci( aci,2)
        
        ! Find interpolated GL position
        TAFi = region%ice%Hi( vi) - ((region%ice%SL( vi) - region%ice%Hb( vi)) * (seawater_density / ice_density))
        TAFj = region%ice%Hi( vj) - ((region%ice%SL( vj) - region%ice%Hb( vj)) * (seawater_density / ice_density))
        lambda = TAFi / (TAFi - TAFj)
        
        x_GL = x_GL + (NORM2(region%mesh%V( vi,:)) * (1._dp - lambda)) + (NORM2(region%mesh%V( vj,:)) * lambda)
      END IF
    END DO
    x_GL = x_GL / n_GL
    
    filename = TRIM(C%output_dir) // 'aa_general_output_' // region%name // '.txt'
    OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
    
    WRITE(UNIT = 1337, FMT = '(F10.1,2F8.2,F13.2,F8.2,9F11.2)') region%time, &
      region%ice_volume, region%ice_volume_above_flotation, region%ice_area, T2m_mean, &
      total_snowfall, total_rainfall, total_melt, total_refreezing, total_runoff, total_SMB, total_BMB, total_MB, x_GL
    
    CLOSE(UNIT = 1337)
    
  ! Point-of-interest output
  ! ========================
    
    DO n = 1, region%mesh%nPOI
    
      vi1 = region%mesh%POI_vi(n,1)
      vi2 = region%mesh%POI_vi(n,2)
      vi3 = region%mesh%POI_vi(n,3)
    
      w1  = region%mesh%POI_w(n,1)
      w2  = region%mesh%POI_w(n,2)
      w3  = region%mesh%POI_w(n,3)
      
      Hi_POI = (region%ice%Hi(vi1  ) * w1) + (region%ice%Hi(vi2  ) * w2) + (region%ice%Hi(vi3  ) * w3)
      Hb_POI = (region%ice%Hb(vi1  ) * w1) + (region%ice%Hb(vi2  ) * w2) + (region%ice%Hb(vi3  ) * w3)
      Hs_POI = (region%ice%Hs(vi1  ) * w1) + (region%ice%Hs(vi2  ) * w2) + (region%ice%Hs(vi3  ) * w3)
      Ti_POI = (region%ice%Ti(vi1,:) * w1) + (region%ice%Ti(vi2,:) * w2) + (region%ice%Ti(vi3,:) * w3)
    
      IF (n<10) THEN
        WRITE(ns3,'(I1)') n
        ns(1:2) = '00'
        ns(3:3) = ns3
      ELSEIF (n<100) THEN
        WRITE(ns2,'(I2)') n
        ns(1:1) = '0'
        ns(2:3) = ns3
      ELSE
        WRITE(ns,'(I3)') n
      END IF
    
      filename = TRIM(C%output_dir) // 'ab_POI_' // TRIM(ns) // '_data.txt'
      OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
    
      WRITE(UNIT = 1337, FMT = '(F10.1,3F11.2)', ADVANCE='NO') region%time, Hi_POI, Hb_POI, Hs_POI
      DO k = 1, C%nZ
        WRITE(UNIT = 1337, FMT = '(F11.2)', ADVANCE='NO') Ti_POI(k)
      END DO
      WRITE(UNIT = 1337, FMT = '(A)') ''
      
      CLOSE(UNIT = 1337)
    
    END DO
    
  END SUBROUTINE write_text_output
  
  ! Initialise the entire model region - read initial and PD data, create the mesh,
  ! initialise the ice dynamics, climate and SMB sub models
  SUBROUTINE initialise_model( region, name, matrix)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    CHARACTER(LEN=3),           INTENT(IN)        :: name
    TYPE(type_climate_matrix),  INTENT(IN)        :: matrix
    
    ! Local variables
    CHARACTER(LEN=4)                                   :: filetype_init
              
    ! Basic initialisation
    ! ====================
    
    ! Region name
    region%name      = name    
    IF (region%name == 'NAM') THEN
      region%long_name = 'North America'
      filetype_init    = C%filetype_init_NAM
    ELSE IF (region%name == 'EAS') THEN
      region%long_name = 'Eurasia'
      filetype_init    = C%filetype_init_EAS
    ELSE IF (region%name == 'GRL') THEN
      region%long_name = 'Greenland'
      filetype_init    = C%filetype_init_GRL
    ELSE IF (region%name == 'ANT') THEN
      region%long_name = 'Antarctica'
      filetype_init    = C%filetype_init_ANT
    END IF
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Initialising model region ', region%name, ' (', TRIM(region%long_name), ')...'
    
    ! Initialise last update times
    ! ============================
    
    region%time               = C%start_time_of_run
    
    region%t0_SIA             = C%start_time_of_run
    region%t0_SSA             = C%start_time_of_run
    region%t0_thermo          = C%start_time_of_run
    region%t0_climate         = C%start_time_of_run
    region%t0_SMB             = C%start_time_of_run
    region%t0_BMB             = C%start_time_of_run
    region%t0_output          = C%start_time_of_run
    region%t0_mesh            = C%start_time_of_run
    
    region%t1_SIA             = C%start_time_of_run
    region%t1_SSA             = C%start_time_of_run
    region%t1_thermo          = C%start_time_of_run
    region%t1_climate         = C%start_time_of_run
    region%t1_SMB             = C%start_time_of_run
    region%t1_BMB             = C%start_time_of_run
    region%t1_output          = C%start_time_of_run
    region%t1_mesh            = C%start_time_of_run + C%dt_mesh_min ! So that the mesh won't be updated immediately after starting the run.
    
    region%dt                 = 0._dp
    region%dt_prev            = 1000._dp
    
    region%dt_SIA             = 0._dp
    region%dt_SSA             = 0._dp
    
    region%do_solve_SIA       = .TRUE.
    region%do_solve_SSA       = .TRUE.
    region%do_thermodynamics  = .TRUE.
    region%do_write_output    = .TRUE.
    
    ! Allocate shared memory for the ice sheet metadata (area, volume, GMSL contribution)    
    CALL allocate_shared_dp_0D( region%ice_area,                      region%wice_area                      )
    CALL allocate_shared_dp_0D( region%ice_volume,                    region%wice_volume                    )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation,    region%wice_volume_above_flotation    )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation_PD, region%wice_volume_above_flotation_PD )
    CALL allocate_shared_dp_0D( region%GMSL_contribution,             region%wGMSL_contribution             )
    
    ! ===== PD and init reference data fields =====
    ! =============================================
    
    IF (par%master) WRITE(0,*) '  Reading reference data...'
    ! Memory is allocated in all processes using allocate_shared, data is only read from the NetCDF files by the master.
    CALL initialise_PD_data_fields(  region%PD,   region%name)
    CALL calculate_PD_sealevel_contribution(region)
    
    IF (filetype_init == 'cart') THEN
      CALL initialise_init_data_fields_cart(region%init, region%name)
    ELSE
      CALL initialise_init_data_fields_mesh(region%init, region%name, region%mesh)
    END IF
    
    ! ===== The mesh =====
    ! ====================
    
    IF (filetype_init == 'cart') THEN
      IF (par%master) WRITE(0,*) '  Creating the first mesh...'
      CALL create_mesh_from_cart_data(region)
    ELSE
      ! basic mesh data has already been read from the NetCDF file, calculate additional data.
      region%mesh%region_name = region%name
      CALL allocate_mesh_secondary(  region%mesh)       
      CALL get_lat_lon_coordinates(    region%mesh)
      CALL find_Voronoi_cell_areas(    region%mesh)
      CALL find_triangle_areas(       region%mesh)
      CALL find_connection_widths(    region%mesh)
      CALL make_Ac_mesh(              region%mesh)
      CALL get_neighbour_functions(   region%mesh)
      CALL determine_mesh_resolution( region%mesh)
    END IF
    CALL sync
    
    ! ===== Map PD and init reference data fields to the mesh =====
    ! =============================================================
    
    IF (par%master) WRITE(0,*) '  Mapping reference data to the mesh...'
    CALL map_PD_data_to_mesh(   region%mesh, region%PD)
    
    IF (filetype_init == 'cart') THEN
      CALL map_init_data_to_mesh( region%mesh, region%init)
    END IF
       
    ! ===== Output file =====
    ! =======================
  
    IF (par%master) CALL create_output_file(region)
    region%output_file_exists = .TRUE.
    CALL sync
        
    ! ===== The climate model =====
    ! =============================    
    
    IF (par%master) WRITE (0,*) '  Initialising climate model...'
    CALL initialise_climate_model(region%climate, matrix, region%mesh)
    CALL sync    
    
    ! ===== The SMB model =====
    ! =========================    
    
    IF (par%master) WRITE (0,*) '  Initialising SMB model...'
    CALL initialise_SMB_model(region%mesh, region%init, filetype_init, region%SMB)
    CALL sync     
    
    ! ===== The BMB model =====
    ! =========================    
    
    IF (par%master) WRITE (0,*) '  Initialising BMB model...'
    CALL initialise_BMB_model(region%mesh, region%BMB)    
    CALL sync   
  
    ! ===== The ice dynamics model
    ! ============================
    
    IF (par%master) WRITE(0,*) '  Initialising ice dynamics model...'
    CALL initialise_ice_model( region%mesh, region%ice, region%init, filetype_init)
    
    ! Run the climate, BMB and SMB models once, to get the correct surface temperature field for the ice temperature initialisation,
    ! and so that all the relevant fields already have sensible values in the first time frame of the output file.
    CALL run_climate_model( region%mesh, region%ice, region%climate, C%start_time_of_run)
    CALL run_SMB_model(     region%mesh, region%ice, region%climate, C%start_time_of_run, region%SMB)
    CALL run_BMB_model(     region%mesh, region%ice, region%BMB)
    
    ! Initialise the ice temperature profile
    CALL initialise_ice_temperature( region%mesh, region%ice, region%init, filetype_init, region%climate)
    
    ! Calculate physical properties again, now with the initialised temperature profile, determine the masks and slopes
    CALL update_general_ice_model_data( region%mesh, region%ice, C%start_time_of_run)
    
    ! Calculate ice sheet metadata (volume, area, GMSL contribution) for writing to the first line of the output file
    CALL calculate_icesheet_volume_and_area( region)
    
    IF (par%master) WRITE (0,*) ' Finished initialising model region ', region%name, '.'
    
    ! ===== ASCII output (computation time tracking, general output)
    IF (par%master) CALL create_text_output_files( region)
    
  END SUBROUTINE initialise_model

END MODULE UFEMISM_main_model
