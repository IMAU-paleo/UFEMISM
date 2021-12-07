MODULE UFEMISM_main_model

  USE mpi
  USE configuration_module,          ONLY: dp, C
  USE parallel_module,               ONLY: par, sync, ierr, cerr, write_to_memory_log, &
                                           allocate_shared_int_0D, allocate_shared_dp_0D, &
                                           allocate_shared_int_1D, allocate_shared_dp_1D, &
                                           allocate_shared_int_2D, allocate_shared_dp_2D, &
                                           allocate_shared_int_3D, allocate_shared_dp_3D, &
                                           deallocate_shared
  USE data_types_module,             ONLY: type_model_region, type_mesh, type_grid, type_climate_matrix, type_remapping
  USE parameters_module,             ONLY: seawater_density, ice_density
  USE reference_fields_module,       ONLY: initialise_PD_data_fields, initialise_init_data_fields, &
                                           map_PD_data_to_mesh, map_init_data_to_mesh
  USE mesh_memory_module,            ONLY: deallocate_mesh_all
  USE mesh_help_functions_module,    ONLY: inverse_oblique_sg_projection, partition_list
  USE mesh_creation_module,          ONLY: create_mesh_from_cart_data
  USE mesh_mapping_module,           ONLY: create_remapping_arrays, deallocate_remapping_arrays, reallocate_field_dp, &
                                           create_remapping_arrays_mesh_grid, deallocate_remapping_arrays_mesh_grid
  USE mesh_update_module,            ONLY: determine_mesh_fitness, create_new_mesh
  USE netcdf_module,                 ONLY: debug, write_to_debug_file, create_output_files, write_to_output_files, &
                                           initialise_debug_fields, associate_debug_fields, reallocate_debug_fields, create_debug_file
  USE restart_module,                ONLY: read_mesh_from_restart_file, read_init_data_from_restart_file
  USE general_ice_model_data_module, ONLY: ice_physical_properties, update_general_ice_model_data
  USE ice_dynamics_module,           ONLY: initialise_ice_model, remap_ice_model, solve_SIA, solve_SSA, calculate_ice_thickness_change
  USE thermodynamics_module,         ONLY: initialise_ice_temperature, update_ice_temperature
  USE climate_module,                ONLY: initialise_climate_model,  remap_climate_model,  run_climate_model, map_subclimate_to_mesh
  USE SMB_module,                    ONLY: initialise_SMB_model,      remap_SMB_model,      run_SMB_model
  USE BMB_module,                    ONLY: initialise_BMB_model,      remap_BMB_model,      run_BMB_model
  USE isotopes_module,               ONLY: initialise_isotopes_model, remap_isotopes_model, run_isotopes_model, calculate_reference_isotopes
  USE bedrock_ELRA_module,           ONLY: initialise_ELRA_model,     remap_ELRA_model,     run_ELRA_model

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
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE (0,'(A,A,A,A,A,F9.3,A,F9.3,A)') '  Running model region ', region%name, ' (', TRIM(region%long_name), & 
                                                          ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'
    
    ! Set the intermediary pointers in "debug" to this region's debug data fields
    CALL associate_debug_fields(  region)
               
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
    
      ! Update bedrock elevation
      IF (C%choice_GIA_model == 'ELRA') THEN
        CALL run_ELRA_model( region)
      ELSE
        WRITE(0,*) '  ERROR - choice_GIA_model "', C%choice_GIA_model, '" not implemented in run_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
          
      ! Update ice thickness based on ice velocity and mass balance
      IF (par%master) t2 = MPI_WTIME()
      CALL calculate_ice_thickness_change( region%mesh, region%ice, region%SMB, region%BMB, region%dt, region%mask_noice)
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
        CALL run_model_update_mesh( region, matrix)
        IF (par%master) region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2
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
        CALL run_climate_model( region%mesh, region%ice, region%SMB, region%climate, region%time, region%name, region%grid_smooth)
        region%t0_climate = region%time
      END IF
    
      ! Run the SMB model
      IF (region%do_SMB) THEN
        CALL run_SMB_model( region%mesh, region%ice, region%climate%applied, region%time, region%SMB, region%mask_noice)
        region%t0_SMB = region%time
      END IF
    
      ! Run the BMB model
      IF (region%do_BMB) THEN
        CALL run_BMB_model( region%mesh, region%ice, region%climate%applied, region%BMB, region%name)
        region%t0_BMB = region%time
      END IF
      
      IF (par%master) region%tcomp_climate = region%tcomp_climate + MPI_WTIME() - t2
      
    ! Thermodynamics
    ! ==============
      
      ! Solve thermodynamics
      IF (region%do_thermodynamics) THEN
        IF (par%master) t2 = MPI_WTIME()
        CALL update_ice_temperature( region%mesh, region%ice, region%climate, region%SMB)
        IF (par%master) region%tcomp_thermo = region%tcomp_thermo + MPI_WTIME() - t2
        region%t0_thermo = region%time
      END IF
      
      ! Surface, basal, and vertically averaged velocities (only diagnostically)
      region%ice%U_surf( region%mesh%v1:region%mesh%v2) = region%ice%U_3D( region%mesh%v1:region%mesh%v2,1)
      region%ice%V_surf( region%mesh%v1:region%mesh%v2) = region%ice%V_3D( region%mesh%v1:region%mesh%v2,1)
      region%ice%U_base( region%mesh%v1:region%mesh%v2) = region%ice%U_3D( region%mesh%v1:region%mesh%v2,C%nZ)
      region%ice%V_base( region%mesh%v1:region%mesh%v2) = region%ice%V_3D( region%mesh%v1:region%mesh%v2,C%nZ)
      region%ice%U_vav(  region%mesh%v1:region%mesh%v2) = region%ice%U_SIA( region%mesh%v1:region%mesh%v2) + &
                                                          region%ice%U_SSA( region%mesh%v1:region%mesh%v2)
      region%ice%V_vav(  region%mesh%v1:region%mesh%v2) = region%ice%V_SIA( region%mesh%v1:region%mesh%v2) + &
                                                          region%ice%V_SSA( region%mesh%v1:region%mesh%v2)
      CALL sync
      
    ! Isotopes
    ! ========
    
      CALL run_isotopes_model( region)
      
    ! Time step and output
    ! ====================
                            
      ! Write output
      IF (region%do_write_output) THEN
        IF (par%master) t2 = MPI_WTIME()
        ! If the mesh has been updated, create a new NetCDF file
        IF (.NOT. region%output_file_exists) THEN
          CALL create_output_files( region)
          CALL sync
          region%output_file_exists = .TRUE.
        END IF  
        CALL write_to_output_files( region)
        region%t0_output = region%time
        IF (par%master) region%tcomp_output = region%tcomp_output + MPI_WTIME() - t2
      END IF
      
      ! Calculate time steps, determine actions and advance region time
      IF (par%master) t2 = MPI_WTIME()
      CALL determine_timesteps_and_actions( region, t_end)
      IF (par%master) region%tcomp_SIA = region%tcomp_SIA + MPI_WTIME() - t2
      
      ! DENK DROM
      !region%time = t_end
    
    END DO ! DO WHILE (region%time < t_end)
    
    ! Write to NetCDF output one last time at the end of the simulation
    IF (region%time == C%end_time_of_run) THEN
      ! If the mesh has been updated, create a new NetCDF file
      IF (.NOT. region%output_file_exists) THEN
        CALL create_output_files( region)
        CALL sync
        region%output_file_exists = .TRUE.
      END IF      
      CALL write_to_output_files( region)
    END IF 
    
    ! Determine total ice sheet area, volume, volume-above-flotation and GMSL contribution,
    ! used for writing to text output and in the inverse routine
    CALL calculate_icesheet_volume_and_area( region)
    
    ! Write to text output
    IF (par%master) THEN
      region%tcomp_total = MPI_WTIME() - t1
      CALL write_text_output( region)
    END IF
    
  END SUBROUTINE run_model
  
  ! Update the mesh and everything else that needs it
  SUBROUTINE run_model_update_mesh( region, matrix)
    ! Perform a mesh update: create a new mesh based on the current modelled ice-sheet
    ! geometry, map all the model data from the old to the new mesh. Deallocate the
    ! old mesh, update the output files and the square grid maps.
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_climate_matrix),           INTENT(IN)    :: matrix
    
    ! Local variables
    TYPE(type_remapping)                               :: map
    
    ! Create a new mesh
    CALL create_new_mesh( region)
    
    IF (par%master) WRITE(0,*) '  Mapping model data to the new mesh...'
    
    ! Calculate the mapping arrays
    CALL create_remapping_arrays( region%mesh, region%mesh_new, map)
    
    ! Update the square grid mapping arrays (needed for remapping/reinitialising the climate model)
    CALL deallocate_remapping_arrays_mesh_grid(              region%grid_output)
    CALL deallocate_remapping_arrays_mesh_grid(              region%grid_GIA   )
    CALL deallocate_remapping_arrays_mesh_grid(              region%grid_smooth)
    CALL create_remapping_arrays_mesh_grid( region%mesh_new, region%grid_GIA   )
    CALL create_remapping_arrays_mesh_grid( region%mesh_new, region%grid_output)
    CALL create_remapping_arrays_mesh_grid( region%mesh_new, region%grid_smooth)
    
    ! Recalculate the "no ice" mask (also needed for the climate model)
    CALL deallocate_shared(      region%wmask_noice)
    CALL allocate_shared_int_1D( region%mesh_new%nV, region%mask_noice, region%wmask_noice)
    CALL initialise_mask_noice(  region, region%mesh_new)
            
    ! Reallocate and remap reference PD data 
    CALL deallocate_shared(  region%PD%wHi  )
    CALL deallocate_shared(  region%PD%wHb  )
    CALL deallocate_shared(  region%PD%wHs  )
    CALL map_PD_data_to_mesh(  region%mesh_new,    region%PD)
    
    ! Remap all the submodels
    CALL remap_ELRA_model(     region%mesh, region%mesh_new, map, region%ice, region%PD, region%grid_GIA)
    CALL remap_ice_model(      region%mesh, region%mesh_new, map, region%ice, region%PD)
    CALL remap_climate_model(  region%mesh, region%mesh_new, map, region%climate, matrix, region%PD, region%grid_smooth, region%mask_noice, region%name)
    CALL remap_SMB_model(      region%mesh, region%mesh_new, map, region%SMB)
    CALL remap_BMB_model(      region%mesh, region%mesh_new, map, region%BMB)
    CALL remap_isotopes_model( region%mesh, region%mesh_new, map, region)
            
    ! Deallocate shared memory for the mapping arrays
    CALL deallocate_remapping_arrays( map)
    
    ! Deallocate the old mesh, bind the region%mesh pointers to the new mesh.
    CALL deallocate_mesh_all( region%mesh)
    region%mesh = region%mesh_new
    
    ! When the next output is written, new output files must be created.
    region%output_file_exists = .FALSE.
    
    ! Reallocate the debug fields for the new mesh, create a new debug file
    CALL reallocate_debug_fields( region)
    CALL create_debug_file(       region)
    
    ! Recalculate the reference precipitation isotope content (must be done after region%mesh has been cycled)
    CALL calculate_reference_isotopes( region)
    
    ! Run all model components again after updating the mesh
    region%do_solve_SIA      = .TRUE.
    region%do_solve_SSA      = .TRUE.
    region%do_climate        = .TRUE.
    region%do_SMB            = .TRUE.
    region%do_BMB            = .TRUE.
    region%do_ELRA           = .TRUE.
    region%do_thermodynamics = .TRUE.
      
  END SUBROUTINE run_model_update_mesh
  
  ! Initialise the entire model region - read initial and PD data, create the mesh,
  ! initialise the ice dynamics, climate and SMB sub models
  SUBROUTINE initialise_model( region, name, matrix)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    CHARACTER(LEN=3),           INTENT(IN)        :: name
    TYPE(type_climate_matrix),  INTENT(IN)        :: matrix
    
    ! Local variables
    CHARACTER(LEN=64), PARAMETER                  :: routine_name = 'initialise_model'
    INTEGER                                       :: n1, n2
    
    n1 = par%mem%n
              
    ! Basic initialisation
    ! ====================
    
    ! Region name
    region%name      = name    
    IF (region%name == 'NAM') THEN
      region%long_name = 'North America'
    ELSE IF (region%name == 'EAS') THEN
      region%long_name = 'Eurasia'
    ELSE IF (region%name == 'GRL') THEN
      region%long_name = 'Greenland'
    ELSE IF (region%name == 'ANT') THEN
      region%long_name = 'Antarctica'
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
    region%t0_ELRA            = C%start_time_of_run
    
    region%t1_SIA             = C%start_time_of_run
    region%t1_SSA             = C%start_time_of_run
    region%t1_thermo          = C%start_time_of_run + C%dt_thermo
    region%t1_climate         = C%start_time_of_run
    region%t1_SMB             = C%start_time_of_run
    region%t1_BMB             = C%start_time_of_run
    region%t1_output          = C%start_time_of_run
    region%t1_mesh            = C%start_time_of_run + C%dt_mesh_min ! So that the mesh won't be updated immediately after starting the run.
    region%t1_ELRA            = C%start_time_of_run
    
    region%dt                 = 0._dp
    region%dt_prev            = 1000._dp
    
    region%dt_SIA             = 0._dp
    region%dt_SSA             = 0._dp
    
    region%do_solve_SIA       = .TRUE.
    region%do_solve_SSA       = .TRUE.
    region%do_thermodynamics  = .FALSE.
    region%do_climate         = .TRUE.
    region%do_SMB             = .TRUE.
    region%do_BMB             = .TRUE.
    region%do_ELRA            = .TRUE.
    region%do_write_output    = .TRUE.
    
    ! Allocate shared memory for the ice sheet metadata (area, volume, GMSL contribution)    
    CALL allocate_shared_dp_0D( region%ice_area                     , region%wice_area                     )
    CALL allocate_shared_dp_0D( region%ice_volume                   , region%wice_volume                   )
    CALL allocate_shared_dp_0D( region%ice_volume_PD                , region%wice_volume_PD                )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation   , region%wice_volume_above_flotation   )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation_PD, region%wice_volume_above_flotation_PD)
    CALL allocate_shared_dp_0D( region%GMSL_contribution            , region%wGMSL_contribution            )
    CALL allocate_shared_dp_0D( region%mean_isotope_content         , region%wmean_isotope_content         )
    CALL allocate_shared_dp_0D( region%mean_isotope_content_PD      , region%wmean_isotope_content_PD      )
    CALL allocate_shared_dp_0D( region%d18O_contribution            , region%wd18O_contribution            )
    CALL allocate_shared_dp_0D( region%d18O_contribution_PD         , region%wd18O_contribution_PD         )
    
    ! ===== PD and init reference data fields =====
    ! =============================================
    
    CALL initialise_PD_data_fields(   region%PD,   region%name)
    CALL calculate_PD_sealevel_contribution( region)
    
    ! If we're restarting from a previous run, no need to read an initial file.
    IF (.NOT. C%is_restart) CALL initialise_init_data_fields( region%init, region%name)
    
    ! ===== The mesh =====
    ! ====================
    
    IF (C%is_restart) THEN
      CALL read_mesh_from_restart_file( region)
    ELSE
      CALL create_mesh_from_cart_data( region)
    END IF
    
    ! The different square grids
    ! ==========================
    
    CALL initialise_model_square_grid( region, region%grid_output, C%dx_grid_output)
    CALL initialise_model_square_grid( region, region%grid_GIA,    C%dx_grid_GIA   )
    CALL initialise_model_square_grid( region, region%grid_smooth, C%dx_grid_smooth)
    
    ! ===== Initialise dummy fields for debugging
    ! ===========================================
    
    CALL initialise_debug_fields( region)
    
    ! ===== Map PD and init reference data fields to the mesh =====
    ! =============================================================
    
    CALL map_PD_data_to_mesh(   region%mesh, region%PD  )
    
    IF (.NOT. C%is_restart) THEN
      CALL map_init_data_to_mesh( region%mesh, region%init)
    ELSE
      ! Read initial data from the restart file
      CALL read_init_data_from_restart_file( region)
    END IF
       
    ! ===== Output file =====
    ! =======================
  
    IF (par%master) CALL create_output_files( region)
    CALL associate_debug_fields(  region)
    region%output_file_exists = .TRUE.
    CALL sync
    
    ! ===== The "no ice" mask
    ! =======================

    CALL allocate_shared_int_1D( region%mesh%nV, region%mask_noice, region%wmask_noice)
    CALL initialise_mask_noice( region, region%mesh)
        
    ! ===== The climate model =====
    ! =============================    
    
    CALL initialise_climate_model( region%climate, matrix, region%mesh, region%PD, region%name, region%mask_noice, region%grid_smooth)
    
    ! ===== The SMB model =====
    ! =========================    
    
    CALL initialise_SMB_model( region%mesh, region%init, region%SMB, region%name)
    
    ! ===== The BMB model =====
    ! =========================    
    
    CALL initialise_BMB_model( region%mesh, region%BMB, region%name) 
  
    ! ===== The ice dynamics model
    ! ============================
    
    CALL initialise_ice_model( region%mesh, region%ice, region%init)
    
    ! Run the climate, BMB and SMB models once, to get the correct surface temperature field for the ice temperature initialisation,
    ! and so that all the relevant fields already have sensible values in the first time frame of the output file.
    CALL run_climate_model( region%mesh, region%ice, region%SMB, region%climate, C%start_time_of_run, region%name, region%grid_smooth)
    CALL run_SMB_model(     region%mesh, region%ice, region%climate%applied, C%start_time_of_run, region%SMB, region%mask_noice)
    CALL run_BMB_model(     region%mesh, region%ice, region%climate%applied, region%BMB, region%name)
    
    ! Initialise the ice temperature profile
    CALL initialise_ice_temperature( region%mesh, region%ice, region%init, region%climate)
    
    ! Calculate physical properties again, now with the initialised temperature profile, determine the masks and slopes
    CALL update_general_ice_model_data( region%mesh, region%ice, C%start_time_of_run)
    
    ! Calculate ice sheet metadata (volume, area, GMSL contribution) for writing to the first line of the output file
    CALL calculate_icesheet_volume_and_area( region)
    
    ! ===== The isotopes model =====
    ! ==============================
   
    CALL initialise_isotopes_model( region)
    
    ! ===== The ELRA GIA model =====
    ! ==============================
    
    IF (C%choice_GIA_model == 'ELRA') THEN
      CALL initialise_ELRA_model( region%mesh, region%grid_GIA, region%ice, region%PD)
    ELSE
      WRITE(0,*) '  ERROR - choice_GIA_model "', C%choice_GIA_model, '" not implemented in initialise_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    IF (par%master) WRITE (0,*) ' Finished initialising model region ', region%name, '.'
    
    ! ===== ASCII output (computation time tracking, general output)
    IF (par%master) CALL create_text_output_files( region)
    
    n2 = par%mem%n
    !CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE initialise_model
  SUBROUTINE initialise_model_square_grid( region, grid, dx)
    ! Initialise a regular square grid enveloping this model region
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    TYPE(type_grid),            INTENT(INOUT)     :: grid
    REAL(dp),                   INTENT(IN)        :: dx
    
    ! Local variables:
    CHARACTER(LEN=64), PARAMETER                  :: routine_name = 'initialise_model_square_grid'
    INTEGER                                       :: n1, n2
    REAL(dp)                                      :: xmid, ymid
    INTEGER                                       :: nsx, nsy, i, j
    
    n1 = par%mem%n
    
    nsx = 0
    nsy = 0
    
    CALL allocate_shared_int_0D( grid%nx,   grid%wnx  )    
    CALL allocate_shared_int_0D( grid%ny,   grid%wny  )
    CALL allocate_shared_dp_0D(  grid%dx,   grid%wdx  )
    CALL allocate_shared_dp_0D(  grid%xmin, grid%wxmin)
    CALL allocate_shared_dp_0D(  grid%xmax, grid%wxmax)
    CALL allocate_shared_dp_0D(  grid%ymin, grid%wymin)
    CALL allocate_shared_dp_0D(  grid%ymax, grid%wymax)
    
    IF (par%master) THEN
    
      ! Resolution
      grid%dx = dx
    
      ! Determine the center of the model domain
      xmid = (region%mesh%xmin + region%mesh%xmax) / 2._dp
      ymid = (region%mesh%ymin + region%mesh%ymax) / 2._dp
      
      nsx = CEILING((region%mesh%xmax - xmid - grid%dx/2._dp) / grid%dx)
      nsy = CEILING((region%mesh%ymax - ymid - grid%dx/2._dp) / grid%dx)
      
      grid%nx = 2*nsx + 1
      grid%ny = 2*nsy + 1
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    CALL allocate_shared_dp_1D( grid%nx, grid%x, grid%wx)
    CALL allocate_shared_dp_1D( grid%ny, grid%y, grid%wy)
    
    IF (par%master) THEN
    
      DO i = 1, grid%nx
        grid%x( i) = -nsx*grid%dx + (i-1)*grid%dx
      END DO
      DO j = 1, grid%ny
        grid%y( j) = -nsy*grid%dx + (j-1)*grid%dx
      END DO
      
      grid%xmin = grid%x(1      )
      grid%xmax = grid%x(grid%nx)
      grid%ymin = grid%y(1      )
      grid%ymax = grid%y(grid%ny)
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Assign range to each processor
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)
    
    ! Calculate lat-lon coordinates
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, grid%lat, grid%wlat)
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, grid%lon, grid%wlon)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL inverse_oblique_sg_projection( grid%x( i), grid%y( j), region%mesh%lambda_M, region%mesh%phi_M, region%mesh%alpha_stereo, grid%lon( i,j), grid%lat( i,j))
    END DO
    END DO
    CALL sync
    
    ! Calculate mapping arrays between the mesh and the grid
    CALL create_remapping_arrays_mesh_grid( region%mesh, grid)
    
    n2 = par%mem%n
    CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE initialise_model_square_grid
  SUBROUTINE initialise_mask_noice( region, mesh)
    ! Mask a certain area where no ice is allowed to grow. This is used to "remove"
    ! Greenland from NAM and EAS, and Ellesmere Island from GRL.
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    TYPE(type_mesh),                 INTENT(IN)        :: mesh
    
    ! Local variables:
    CHARACTER(LEN=64), PARAMETER                       :: routine_name = 'initialise_mask_noice'
    INTEGER                                            :: n1, n2
    INTEGER                                            :: vi
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc, pd
    REAL(dp)                                           :: yl_ab, yl_bc, yl_cd
    
    n1 = par%mem%n
    
    region%mask_noice( mesh%v1:mesh%v2) = 0
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_mask_remove_land!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    IF (region%name == 'NAM') THEN
      ! North America: remove Greenland
      
      pa = [ 490000._dp, 1530000._dp]
      pb = [2030000._dp,  570000._dp]
      
      DO vi = mesh%v1, mesh%v2
        yl_ab = pa(2) + (mesh%V( vi,1) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
        IF (mesh%V( vi,2) > yl_ab .AND. mesh%V( vi,1) > pa(1) .AND. mesh%V( vi,2) > pb(2)) THEN
          region%mask_noice( vi) = 1
        END IF
      END DO
      CALL sync
    
    ELSEIF (region%name == 'EAS') THEN
      ! Eurasia: remove Greenland
    
      pa = [-2900000._dp, 1300000._dp]
      pb = [-1895000._dp,  900000._dp]
      pc = [ -835000._dp, 1135000._dp]
      pd = [ -400000._dp, 1855000._dp]
      
      DO vi = mesh%v1, mesh%v2
        yl_ab = pa(2) + (mesh%V( vi,1) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
        yl_bc = pb(2) + (mesh%V( vi,1) - pb(1))*(pc(2)-pb(2))/(pc(1)-pb(1))
        yl_cd = pc(2) + (mesh%V( vi,1) - pc(1))*(pd(2)-pc(2))/(pd(1)-pc(1))
        IF ((mesh%V( vi,1) <  pa(1) .AND. mesh%V( vi,2) > pa(2)) .OR. &
            (mesh%V( vi,1) >= pa(1) .AND. mesh%V( vi,1) < pb(1) .AND. mesh%V( vi,2) > yl_ab) .OR. &
            (mesh%V( vi,1) >= pb(1) .AND. mesh%V( vi,1) < pc(1) .AND. mesh%V( vi,2) > yl_bc) .OR. &
            (mesh%V( vi,1) >= pc(1) .AND. mesh%V( vi,1) < pd(1) .AND. mesh%V( vi,2) > yl_cd)) THEN
          region%mask_noice( vi) = 1
        END IF
      END DO
      CALL sync
      
    ELSEIF (region%name == 'GRL') THEN
      ! Greenland: remove Ellesmere island
      
      pa = [-750000._dp,  900000._dp]
      pb = [-250000._dp, 1250000._dp]
      
      DO vi = mesh%v1, mesh%v2
        yl_ab = pa(2) + (mesh%V( vi,1) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
        IF (mesh%V( vi,2) > pa(2) .AND. mesh%V( vi,2) > yl_ab .AND. mesh%V( vi,1) < pb(1)) THEN
          region%mask_noice( vi) = 1
        END IF
      END DO
      CALL sync
      
    ELSEIF (region%name == 'ANT') THEN
      ! Antarctica: no changes needed
      
    END IF ! IF (region%name == 'NAM') THEN
    
    n2 = par%mem%n
    CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE initialise_mask_noice
  
  ! Calculate critical time steps from the SIA and SSA, and determine which actions should be taken
  SUBROUTINE determine_timesteps_and_actions( region, t_end)
    ! This subroutine calculates the SIA time step with the CFK-criterium
    ! Methods are similar to those in ANICE, which are based on PISM (see also Bueler  et al., 2007)
    
    ! It also determines whether the SIA velocities, SSA velocities or thermodynamics should be updated,
    ! or whether results should be written to the output file
    
    USE parameters_module, ONLY: pi
    
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end
    
    ! Local variables:
    INTEGER                                     :: ci, vi, vj, k
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
      dt_D_2D = dist**2 / (-6._dp * pi * (region%ice%D_SIA_ac(ci) - 1E-09))
      dt_D_2D_min = MIN(dt_D_2D, dt_D_2D_min) 
      
      dt_V_2D_SSA = dist / (ABS(region%ice%U_SSA(vi)) + ABS(region%ice%V_SSA(vi)))
      dt_V_2D_SSA_min = MIN(dt_V_2D_SSA, dt_V_2D_SSA_min)   
      dt_V_2D_SSA = dist / (ABS(region%ice%U_SSA(vj)) + ABS(region%ice%V_SSA(vj)))
      dt_V_2D_SSA_min = MIN(dt_V_2D_SSA, dt_V_2D_SSA_min)  
    END DO
    
    DO vi = region%mesh%v1, region%mesh%v2
      dt_V_2D_SSA = SQRT(region%mesh%A(vi)/pi) / (ABS(region%ice%U_SSA(vi)) + ABS(region%ice%V_SSA(vi)))
      dt_V_2D_SSA_min = MIN(dt_V_2D_SSA, dt_V_2D_SSA_min)
      
      DO k = 1, C%nZ
        dt_V_3D_SIA = SQRT(region%mesh%A(vi)/pi) / (ABS(region%ice%U_3D(vi,k)) + ABS(region%ice%V_3D(vi,k)))
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
    region%t1_ELRA     = region%t0_ELRA     + C%dt_bedrock_ELRA
    region%t1_output   = region%t0_output   + C%dt_output
    
    t_next_action = MINVAL( [region%t1_SIA, region%t1_SSA, region%t1_thermo, region%t1_climate, region%t1_SMB, region%t1_BMB, region%t1_ELRA, region%t1_output])
    
    region%dt = t_next_action - region%time
    
    region%do_solve_SIA       = .FALSE.
    region%do_solve_SSA       = .FALSE.
    region%do_thermodynamics  = .FALSE.
    region%do_climate         = .FALSE.
    region%do_SMB             = .FALSE.
    region%do_BMB             = .FALSE.
    region%do_ELRA            = .FALSE.
    region%do_write_output    = .FALSE.
    
    IF (t_next_action == region%t1_SIA     ) region%do_solve_SIA      = .TRUE.
    IF (t_next_action == region%t1_SSA     ) region%do_solve_SSA      = .TRUE.
    IF (t_next_action == region%t1_thermo  ) region%do_thermodynamics = .TRUE.
    IF (t_next_action == region%t1_climate ) region%do_climate        = .TRUE.
    IF (t_next_action == region%t1_SMB     ) region%do_SMB            = .TRUE.
    IF (t_next_action == region%t1_BMB     ) region%do_BMB            = .TRUE.
    IF (t_next_action == region%t1_ELRA    ) region%do_ELRA           = .TRUE.
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
    
    INTEGER                                       :: vi
    REAL(dp)                                      :: ice_area, ice_volume, thickness_above_flotation, ice_volume_above_flotation
    
    ice_area                   = 0._dp
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp
    
    ! Calculate ice area and volume for processor domain
    DO vi = region%mesh%v1, region%mesh%v2
    
      IF (region%ice%mask_ice( vi) == 1) THEN
        ice_volume = ice_volume + (region%ice%Hi( vi) * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area))
        ice_area   = ice_area   + region%mesh%A( vi) * 1.0E-06_dp ! [km^3]
    
        ! Thickness above flotation
        thickness_above_flotation = MAX(0._dp, region%ice%Hi( vi) - MAX(0._dp, (region%ice%SL( vi) - region%ice%Hb( vi)) * (seawater_density / ice_density)))
        
        ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area)
      END IF     
      
    END DO
    CALL sync
    
    CALL MPI_REDUCE( ice_area,                   region%ice_area,                   1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume,                 region%ice_volume,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    ! Calculate GMSL contribution
    IF (par%master) region%GMSL_contribution = -1._dp * (region%ice_volume_above_flotation - region%ice_volume_above_flotation_PD)
    CALL sync
    
  END SUBROUTINE calculate_icesheet_volume_and_area
  SUBROUTINE calculate_PD_sealevel_contribution( region)
    
    USE parameters_module,           ONLY: ocean_area, seawater_density, ice_density
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    INTEGER                                       :: i, j
    REAL(dp)                                      :: ice_volume, thickness_above_flotation, ice_volume_above_flotation
    
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp
    
    DO i = region%PD%grid%i1, region%PD%grid%i2
    DO j = 1, region%PD%grid%ny
    
      ! Thickness above flotation
      IF (region%PD%Hi_grid( i,j) > 0._dp) THEN
        thickness_above_flotation = MAX(0._dp, region%PD%Hi_grid( i,j) - MAX(0._dp, (0._dp - region%PD%Hb_grid( i,j)) * (seawater_density / ice_density)))
      ELSE
        thickness_above_flotation = 0._dp
      END IF
      
      ! Ice volume (above flotation) in m.s.l.e
      ice_volume                 = ice_volume                 + region%PD%Hi_grid( i,j)   * region%PD%grid%dx * region%PD%grid%dx * ice_density / (seawater_density * ocean_area)
      ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation * region%PD%grid%dx * region%PD%grid%dx * ice_density / (seawater_density * ocean_area)
      
    END DO
    END DO
    
    CALL MPI_REDUCE( ice_volume                , region%ice_volume_PD,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation_PD, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
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
    
      filename = TRIM(C%output_dir) // 'ab_POI_' // TRIM(region%name) // '_' // TRIM(ns) // '_data.txt'
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
    
      filename = TRIM(C%output_dir) // 'ab_POI_' // TRIM(region%name) // '_' // TRIM(ns) // '_data.txt'
      OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
    
      WRITE(UNIT = 1337, FMT = '(F10.1,3F11.2)', ADVANCE='NO') region%time, Hi_POI, Hb_POI, Hs_POI
      DO k = 1, C%nZ
        WRITE(UNIT = 1337, FMT = '(F11.2)', ADVANCE='NO') Ti_POI(k)
      END DO
      WRITE(UNIT = 1337, FMT = '(A)') ''
      
      CLOSE(UNIT = 1337)
    
    END DO
    
  END SUBROUTINE write_text_output

END MODULE UFEMISM_main_model
