MODULE UFEMISM_main_model

  ! The main regional ice-sheet model

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
  USE data_types_module,               ONLY: type_model_region, type_mesh, type_grid, type_climate_matrix, type_remapping
  USE reference_fields_module,         ONLY: initialise_PD_data_fields, initialise_init_data_fields, &
                                             map_PD_data_to_mesh, map_init_data_to_mesh
  USE mesh_memory_module,              ONLY: deallocate_mesh_all
  USE mesh_help_functions_module,      ONLY: inverse_oblique_sg_projection
  USE mesh_creation_module,            ONLY: create_mesh_from_cart_data
  USE mesh_mapping_module,             ONLY: create_remapping_arrays, deallocate_remapping_arrays, &
                                             create_remapping_arrays_mesh_grid, deallocate_remapping_arrays_mesh_grid
  USE mesh_update_module,              ONLY: determine_mesh_fitness, create_new_mesh
  USE netcdf_module,                   ONLY: create_output_files, write_to_output_files, initialise_debug_fields, &
                                             associate_debug_fields, reallocate_debug_fields, create_debug_file, write_CSR_matrix_to_NetCDF
  USE restart_module,                  ONLY: read_mesh_from_restart_file, read_init_data_from_restart_file
  USE general_ice_model_data_module,   ONLY: update_general_ice_model_data
  USE ice_dynamics_module,             ONLY: initialise_ice_model,      remap_ice_model,      run_ice_model, update_ice_thickness
  USE thermodynamics_module,           ONLY: initialise_ice_temperature,                      run_thermo_model, calc_ice_rheology
  USE climate_module,                  ONLY: initialise_climate_model,  remap_climate_model,  run_climate_model, map_subclimate_to_mesh
  USE SMB_module,                      ONLY: initialise_SMB_model,      remap_SMB_model,      run_SMB_model
  USE BMB_module,                      ONLY: initialise_BMB_model,      remap_BMB_model,      run_BMB_model
  USE isotopes_module,                 ONLY: initialise_isotopes_model, remap_isotopes_model, run_isotopes_model, calculate_reference_isotopes
  USE bedrock_ELRA_module,             ONLY: initialise_ELRA_model,     remap_ELRA_model,     run_ELRA_model

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
    INTEGER                                       :: it
    REAL(dp)                                      :: meshfitness
    REAL(dp)                                      :: t1, t2
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE (0,'(A,A,A,A,A,F9.3,A,F9.3,A)') '  Running model region ', region%name, ' (', TRIM(region%long_name), & 
                                                          ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'
    
    ! Set the intermediary pointers in "debug" to this region's debug data fields
    CALL associate_debug_fields(  region)
               
    ! Computation time tracking
    region%tcomp_total          = 0._dp
    region%tcomp_ice            = 0._dp
    region%tcomp_thermo         = 0._dp
    region%tcomp_climate        = 0._dp
    region%tcomp_GIA            = 0._dp
    region%tcomp_mesh           = 0._dp
       
    t1 = MPI_WTIME()
    t2 = 0._dp
                            
  ! ====================================
  ! ===== The main model time loop =====
  ! ====================================
    
    it = 0
    DO WHILE (region%time < t_end)
      it = it + 1
      
    ! GIA
    ! ===
    
      t1 = MPI_WTIME()
      IF     (C%choice_GIA_model == 'none') THEN
        ! Nothing to be done
      ELSEIF (C%choice_GIA_model == 'ELRA') THEN
        CALL run_ELRA_model( region)
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR - choice_GIA_model "', C%choice_GIA_model, '" not implemented in run_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_GIA = region%tcomp_GIA + t2 - t1
      
    ! Ice dynamics
    ! ============
    
      ! Calculate ice velocities and the resulting change in ice geometry
      ! NOTE: geometry is not updated yet; this happens at the end of the time loop
      t1 = MPI_WTIME()
      CALL run_ice_model( region, t_end)
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_ice = region%tcomp_ice + t2 - t1
      
    ! Mesh update
    ! ===========
      
      ! Check if the mesh needs to be updated
      IF (par%master) t2 = MPI_WTIME()
      meshfitness = 1._dp
      IF (region%time > region%t_last_mesh) THEN
        CALL determine_mesh_fitness(region%mesh, region%ice, meshfitness)    
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
      
    ! Climate , SMB and BMB
    ! =====================
    
      t1 = MPI_WTIME()
            
      ! Run the climate model
      IF (region%do_climate) THEN
        CALL run_climate_model( region%mesh, region%ice, region%SMB, region%climate, region%time, region%name, region%grid_smooth)
      END IF     
    
      ! Run the SMB model
      IF (region%do_SMB) THEN
        CALL run_SMB_model( region%mesh, region%ice, region%climate%applied, region%time, region%SMB, region%mask_noice)
      END IF
    
      ! Run the BMB model
      IF (region%do_BMB) THEN
        CALL run_BMB_model( region%mesh, region%ice, region%climate%applied, region%BMB, region%name)
      END IF
      
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_climate = region%tcomp_climate + t2 - t1
      
    ! Thermodynamics
    ! ==============
    
      t1 = MPI_WTIME()
      CALL run_thermo_model( region%mesh, region%ice, region%climate%applied, region%SMB, region%time, do_solve_heat_equation = region%do_thermo)
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_thermo = region%tcomp_thermo + t2 - t1
      
    ! Isotopes
    ! ========
    
      CALL run_isotopes_model( region)
      
    ! Time step and output
    ! ====================
                            
      ! Write output
      IF (region%do_output) THEN
        ! If the mesh has been updated, create a new NetCDF file
        IF (.NOT. region%output_file_exists) THEN
          CALL create_output_files( region)
          CALL sync
          region%output_file_exists = .TRUE.
        END IF  
        CALL write_to_output_files( region)
      END IF

      ! Update ice geometry and advance region time
      CALL update_ice_thickness( region%mesh, region%ice)
      IF (par%master) region%time = region%time + region%dt
      CALL sync
      
      ! DENK DROM
      !region%time = t_end
    
    END DO ! DO WHILE (region%time < t_end)
                            
  ! ===========================================
  ! ===== End of the main model time loop =====
  ! ===========================================
    
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
    region%do_SIA     = .TRUE.
    region%do_SSA     = .TRUE.
    region%do_DIVA    = .TRUE.
    region%do_thermo  = .TRUE.
    region%do_climate = .TRUE.
    region%do_SMB     = .TRUE.
    region%do_BMB     = .TRUE.
    region%do_output  = .TRUE.
    region%do_ELRA    = .TRUE.
      
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
    
    ! ===== Allocate memory for timers and scalars =====
    ! ==================================================
    
    CALL allocate_region_timers_and_scalars( region)
    
    ! ===== PD and init reference data fields =====
    ! =============================================
    
    CALL initialise_PD_data_fields( region%PD, region%name)
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
    
    
    ! DENK DROM
    CALL test_matrix_operators_mesh(      region%mesh)
    CALL solve_modified_Laplace_equation_c( region%mesh)
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    
    
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
    
    ! Run the climate and SMB models once, to get the correct surface temperature+SMB fields for the ice temperature initialisation
    CALL run_climate_model( region%mesh, region%ice, region%SMB, region%climate, C%start_time_of_run, region%name, region%grid_smooth)
    CALL run_SMB_model(     region%mesh, region%ice, region%climate%applied, C%start_time_of_run, region%SMB, region%mask_noice)
    
    ! Initialise the ice temperature profile
    CALL initialise_ice_temperature( region%mesh, region%ice, region%climate%applied, region%SMB, region%name)
    
    ! Initialise the rheology
    CALL calc_ice_rheology( region%mesh, region%ice, C%start_time_of_run)
    
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
    
!    ! ===== ASCII output (computation time tracking, general output)
!    IF (par%master) CALL create_text_output_files( region)
    
    n2 = par%mem%n
    !CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE initialise_model
  SUBROUTINE allocate_region_timers_and_scalars( region)
    ! Allocate shared memory for this region's timers (used for the asynchronous coupling between the
    ! ice dynamics and the secondary model components), and for the scalars (integrated ice volume and
    ! area, SMB components, computation times, etc.)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    ! Timers and time steps
    ! =====================
    
    CALL allocate_shared_dp_0D(   region%time,             region%wtime            )
    CALL allocate_shared_dp_0D(   region%dt,               region%wdt              )
    CALL allocate_shared_dp_0D(   region%dt_prev,          region%wdt_prev         )
    CALL allocate_shared_dp_0D(   region%dt_crit_SIA,      region%wdt_crit_SIA     )
    CALL allocate_shared_dp_0D(   region%dt_crit_SSA,      region%wdt_crit_SSA     )
    CALL allocate_shared_dp_0D(   region%dt_crit_ice,      region%wdt_crit_ice     )
    CALL allocate_shared_dp_0D(   region%dt_crit_ice_prev, region%wdt_crit_ice_prev)
    
    CALL allocate_shared_dp_0D(   region%t_last_mesh,      region%wt_last_mesh     )
    CALL allocate_shared_dp_0D(   region%t_next_mesh,      region%wt_next_mesh     )
    CALL allocate_shared_bool_0D( region%do_mesh,          region%wdo_mesh         ) 
    
    CALL allocate_shared_dp_0D(   region%t_last_SIA,       region%wt_last_SIA      )
    CALL allocate_shared_dp_0D(   region%t_next_SIA,       region%wt_next_SIA      )
    CALL allocate_shared_bool_0D( region%do_SIA,           region%wdo_SIA          )
    
    CALL allocate_shared_dp_0D(   region%t_last_SSA,       region%wt_last_SSA      )
    CALL allocate_shared_dp_0D(   region%t_next_SSA,       region%wt_next_SSA      )
    CALL allocate_shared_bool_0D( region%do_SSA,           region%wdo_SSA          )
    
    CALL allocate_shared_dp_0D(   region%t_last_DIVA,      region%wt_last_DIVA     )
    CALL allocate_shared_dp_0D(   region%t_next_DIVA,      region%wt_next_DIVA     )
    CALL allocate_shared_bool_0D( region%do_DIVA,          region%wdo_DIVA         )
    
    CALL allocate_shared_dp_0D(   region%t_last_thermo,    region%wt_last_thermo   )
    CALL allocate_shared_dp_0D(   region%t_next_thermo,    region%wt_next_thermo   )
    CALL allocate_shared_bool_0D( region%do_thermo,        region%wdo_thermo       )
    
    CALL allocate_shared_dp_0D(   region%t_last_climate,   region%wt_last_climate  )
    CALL allocate_shared_dp_0D(   region%t_next_climate,   region%wt_next_climate  )
    CALL allocate_shared_bool_0D( region%do_climate,       region%wdo_climate      )
    
    CALL allocate_shared_dp_0D(   region%t_last_ocean,     region%wt_last_ocean    )
    CALL allocate_shared_dp_0D(   region%t_next_ocean,     region%wt_next_ocean    )
    CALL allocate_shared_bool_0D( region%do_ocean,         region%wdo_ocean        )
    
    CALL allocate_shared_dp_0D(   region%t_last_SMB,       region%wt_last_SMB      )
    CALL allocate_shared_dp_0D(   region%t_next_SMB,       region%wt_next_SMB      )
    CALL allocate_shared_bool_0D( region%do_SMB,           region%wdo_SMB          )
    
    CALL allocate_shared_dp_0D(   region%t_last_BMB,       region%wt_last_BMB      )
    CALL allocate_shared_dp_0D(   region%t_next_BMB,       region%wt_next_BMB      )
    CALL allocate_shared_bool_0D( region%do_BMB,           region%wdo_BMB          )
    
    CALL allocate_shared_dp_0D(   region%t_last_ELRA,      region%wt_last_ELRA     )
    CALL allocate_shared_dp_0D(   region%t_next_ELRA,      region%wt_next_ELRA     )
    CALL allocate_shared_bool_0D( region%do_ELRA,          region%wdo_ELRA         )
    
    CALL allocate_shared_dp_0D(   region%t_last_output,    region%wt_last_output   )
    CALL allocate_shared_dp_0D(   region%t_next_output,    region%wt_next_output   )
    CALL allocate_shared_bool_0D( region%do_output,        region%wdo_output       )
    
    IF (par%master) THEN
      region%time           = C%start_time_of_run
      region%dt             = C%dt_min
      region%dt_prev        = C%dt_min
      
      region%t_last_mesh    = C%start_time_of_run
      region%t_next_mesh    = C%start_time_of_run + C%dt_mesh_min
      region%do_mesh        = .FALSE.
      
      region%t_last_SIA     = C%start_time_of_run
      region%t_next_SIA     = C%start_time_of_run
      region%do_SIA         = .TRUE.
      
      region%t_last_SSA     = C%start_time_of_run
      region%t_next_SSA     = C%start_time_of_run
      region%do_SSA         = .TRUE.
      
      region%t_last_DIVA    = C%start_time_of_run
      region%t_next_DIVA    = C%start_time_of_run
      region%do_DIVA        = .TRUE.
      
      region%t_last_thermo  = C%start_time_of_run
      region%t_next_thermo  = C%start_time_of_run + C%dt_thermo
      region%do_thermo      = .FALSE.
      
      region%t_last_climate = C%start_time_of_run
      region%t_next_climate = C%start_time_of_run
      region%do_climate     = .TRUE.
      
      region%t_last_ocean   = C%start_time_of_run
      region%t_next_ocean   = C%start_time_of_run
      region%do_ocean       = .TRUE.
      
      region%t_last_SMB     = C%start_time_of_run
      region%t_next_SMB     = C%start_time_of_run
      region%do_SMB         = .TRUE.
      
      region%t_last_BMB     = C%start_time_of_run
      region%t_next_BMB     = C%start_time_of_run
      region%do_BMB         = .TRUE.
      
      region%t_last_ELRA    = C%start_time_of_run
      region%t_next_ELRA    = C%start_time_of_run
      region%do_ELRA        = .TRUE.
      
      region%t_last_output  = C%start_time_of_run
      region%t_next_output  = C%start_time_of_run
      region%do_output      = .TRUE.
    END IF
    
    ! ===== Scalars =====
    ! ===================
    
    ! Ice-sheet volume and area
    CALL allocate_shared_dp_0D( region%ice_area                     , region%wice_area                     )
    CALL allocate_shared_dp_0D( region%ice_volume                   , region%wice_volume                   )
    CALL allocate_shared_dp_0D( region%ice_volume_PD                , region%wice_volume_PD                )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation   , region%wice_volume_above_flotation   )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation_PD, region%wice_volume_above_flotation_PD)
    
    ! Regionally integrated SMB components
    CALL allocate_shared_dp_0D( region%int_T2m                      , region%wint_T2m                      )
    CALL allocate_shared_dp_0D( region%int_snowfall                 , region%wint_snowfall                 )
    CALL allocate_shared_dp_0D( region%int_rainfall                 , region%wint_rainfall                 )
    CALL allocate_shared_dp_0D( region%int_melt                     , region%wint_melt                     )
    CALL allocate_shared_dp_0D( region%int_refreezing               , region%wint_refreezing               )
    CALL allocate_shared_dp_0D( region%int_runoff                   , region%wint_runoff                   )
    CALL allocate_shared_dp_0D( region%int_SMB                      , region%wint_SMB                      )
    CALL allocate_shared_dp_0D( region%int_BMB                      , region%wint_BMB                      )
    CALL allocate_shared_dp_0D( region%int_MB                       , region%wint_MB                       )
    
    ! Englacial isotope content
    CALL allocate_shared_dp_0D( region%GMSL_contribution            , region%wGMSL_contribution            )
    CALL allocate_shared_dp_0D( region%mean_isotope_content         , region%wmean_isotope_content         )
    CALL allocate_shared_dp_0D( region%mean_isotope_content_PD      , region%wmean_isotope_content_PD      )
    CALL allocate_shared_dp_0D( region%d18O_contribution            , region%wd18O_contribution            )
    CALL allocate_shared_dp_0D( region%d18O_contribution_PD         , region%wd18O_contribution_PD         )
    
    ! Computation times
    CALL allocate_shared_dp_0D( region%tcomp_total                  , region%wtcomp_total                  )
    CALL allocate_shared_dp_0D( region%tcomp_ice                    , region%wtcomp_ice                    )
    CALL allocate_shared_dp_0D( region%tcomp_thermo                 , region%wtcomp_thermo                 )
    CALL allocate_shared_dp_0D( region%tcomp_climate                , region%wtcomp_climate                )
    CALL allocate_shared_dp_0D( region%tcomp_GIA                    , region%wtcomp_GIA                    )
    CALL allocate_shared_dp_0D( region%tcomp_mesh                   , region%wtcomp_mesh                   )
    
  END SUBROUTINE allocate_region_timers_and_scalars
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
    
      IF (region%ice%mask_ice_a( vi) == 1) THEN
        ice_volume = ice_volume + (region%ice%Hi_a( vi) * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area))
        ice_area   = ice_area   + region%mesh%A( vi) * 1.0E-06_dp ! [km^3]
    
        ! Thickness above flotation
        thickness_above_flotation = MAX(0._dp, region%ice%Hi_a( vi) - MAX(0._dp, (region%ice%SL_a( vi) - region%ice%Hb_a( vi)) * (seawater_density / ice_density)))
        
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
  
!  ! Create and write to this region's text output files - both the region-wide one,
!  ! and the ones for the different Points-Of-Interest
!  SUBROUTINE create_text_output_files( region)
!    ! Creates the following text output files:
!    !   time_log_REG.txt             - a log of how much computation time the different model parts take
!    !   general_output_REG.txt       - some general info - ice sheet volume, average surface temperature, total mass balance, etc.
!  
!    IMPLICIT NONE  
!    
!    TYPE(type_model_region),    INTENT(IN)        :: region
!    
!    CHARACTER(LEN=256)                            :: filename
!    CHARACTER(LEN=3)                              :: ns
!    CHARACTER(LEN=1)                              :: ns3
!    CHARACTER(LEN=2)                              :: ns2
!    INTEGER                                       :: n, k
!        
!    ! Time log
!    ! ========
!    
!    filename = TRIM(C%output_dir) // 'aa_time_log_' // region%name // '.txt'
!    OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
!    
!    WRITE(UNIT = 1337, FMT = '(A)') 'Time log for region ' // TRIM(region%long_name)
!    WRITE(UNIT = 1337, FMT = '(A)') 'Computation time (in seconds) required by each model component'
!    WRITE(UNIT = 1337, FMT = '(A)') ''
!    WRITE(UNIT = 1337, FMT = '(A)') '     Time  vertices     total       mesh        SIA        SSA     thermo    climate     output'
!    
!    CLOSE(UNIT = 1337)
!        
!    ! General output
!    ! ==============
!    
!    filename = TRIM(C%output_dir) // 'aa_general_output_' // region%name // '.txt'
!    OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
!    
!    WRITE(UNIT = 1337, FMT = '(A)') 'General output for region ' // TRIM(region%long_name)
!    WRITE(UNIT = 1337, FMT = '(A)') ''
!    WRITE(UNIT = 1337, FMT = '(A)') ' Columns in order:'
!    WRITE(UNIT = 1337, FMT = '(A)') '   1)  Model time                  (years) '
!    WRITE(UNIT = 1337, FMT = '(A)') '   2)  Ice volume                  (meter sea level equivalent)'
!    WRITE(UNIT = 1337, FMT = '(A)') '   3)  Ice volume above flotation  (meter sea level equivalent)'
!    WRITE(UNIT = 1337, FMT = '(A)') '   4)  Ice area                    (km^2)'
!    WRITE(UNIT = 1337, FMT = '(A)') '   5)  Mean surface temperature    (Kelvin)'
!    WRITE(UNIT = 1337, FMT = '(A)') '   6)  Total snowfall     over ice (Gton/y)'
!    WRITE(UNIT = 1337, FMT = '(A)') '   7)  Total rainfall     over ice (Gton/y)'
!    WRITE(UNIT = 1337, FMT = '(A)') '   8)  Total melt         over ice (Gton/y)'
!    WRITE(UNIT = 1337, FMT = '(A)') '   9)  Total refreezing   over ice (Gton/y)'
!    WRITE(UNIT = 1337, FMT = '(A)') '  10)  Total runoff       over ice (Gton/y)'
!    WRITE(UNIT = 1337, FMT = '(A)') '  11)  Total SMB          over ice (Gton/y)'
!    WRITE(UNIT = 1337, FMT = '(A)') '  12)  Total BMB          over ice (Gton/y)'
!    WRITE(UNIT = 1337, FMT = '(A)') '  13)  Total mass balance over ice (Gton/y)'
!    WRITE(UNIT = 1337, FMT = '(A)') '  14)  Grounding line x position (for MISMIP benchmark experiments)'
!    WRITE(UNIT = 1337, FMT = '(A)') ''
!    WRITE(UNIT = 1337, FMT = '(A)') '     Time     Ice  Ice-af     Ice-area     T2m       Snow       Rain       Melt   Refreeze     Runoff        SMB        BMB         MB       x_GL'
!    
!    CLOSE(UNIT = 1337)
!    
!    ! Point-of-interest output
!    ! ========================
!    
!    DO n = 1, region%mesh%nPOI
!    
!      IF (n<10) THEN
!        WRITE(ns3,'(I1)') n
!        ns(1:2) = '00'
!        ns(3:3) = ns3
!      ELSEIF (n<100) THEN
!        WRITE(ns2,'(I2)') n
!        ns(1:1) = '0'
!        ns(2:3) = ns3
!      ELSE
!        WRITE(ns,'(I3)') n
!      END IF
!    
!      filename = TRIM(C%output_dir) // 'ab_POI_' // TRIM(region%name) // '_' // TRIM(ns) // '_data.txt'
!      OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
!    
!      WRITE(UNIT = 1337, FMT = '(A)') 'Relevant data for Point of Interest ' // TRIM(ns)
!      WRITE(UNIT = 1337, FMT = '(A)') ''
!      WRITE(UNIT = 1337, FMT = '(A,F10.2)') '  Lat = ', region%mesh%POI_coordinates(n,1)
!      WRITE(UNIT = 1337, FMT = '(A,F10.2)') '  Lon = ', region%mesh%POI_coordinates(n,2)
!      WRITE(UNIT = 1337, FMT = '(A,F10.2)') '  x   = ', region%mesh%POI_XY_coordinates(n,1)
!      WRITE(UNIT = 1337, FMT = '(A,F10.2)') '  y   = ', region%mesh%POI_XY_coordinates(n,2)
!      WRITE(UNIT = 1337, FMT = '(A,F10.2)') '  res = ', region%mesh%POI_resolutions(n)
!      WRITE(UNIT = 1337, FMT = '(A)') ''
!      WRITE(UNIT = 1337, FMT = '(A)') '  zeta = '
!      DO k = 1, C%nZ
!        WRITE(UNIT = 1337, FMT = '(F22.16)') C%zeta(k)
!      END DO
!      WRITE(UNIT = 1337, FMT = '(A)') ''
!      WRITE(UNIT = 1337, FMT = '(A)') '     Time         Hi         Hb         Hs         Ti'
!      
!      CLOSE(UNIT = 1337)
!    
!    END DO
!    
!  END SUBROUTINE create_text_output_files
!  SUBROUTINE write_text_output( region)
!    ! Write data to the following text output files:
!    !   time_log_REG.txt             - a log of how much computation time the different model parts take
!    !   general_output_REG.txt       - some general info - ice sheet volume, average surface temperature, total mass balance, etc.
!    
!    USE parameters_module,           ONLY: ocean_area, seawater_density, ice_density
!  
!    IMPLICIT NONE  
!    
!    TYPE(type_model_region),    INTENT(IN)        :: region
!    
!    CHARACTER(LEN=256)                            :: filename
!    CHARACTER(LEN=3)                              :: ns
!    CHARACTER(LEN=1)                              :: ns3
!    CHARACTER(LEN=2)                              :: ns2
!    INTEGER                                       :: n
!    INTEGER                                       :: vi, m, k, aci, vj
!    REAL(dp)                                      :: T2m_mean
!    REAL(dp)                                      :: total_snowfall
!    REAL(dp)                                      :: total_rainfall
!    REAL(dp)                                      :: total_melt
!    REAL(dp)                                      :: total_refreezing
!    REAL(dp)                                      :: total_runoff
!    REAL(dp)                                      :: total_SMB
!    REAL(dp)                                      :: total_BMB
!    REAL(dp)                                      :: total_MB
!    REAL(dp)                                      :: TAFi, TAFj, lambda, x_GL
!    INTEGER                                       :: n_GL
!    
!    INTEGER                                       :: vi1, vi2, vi3
!    REAL(dp)                                      :: w1, w2, w3
!    REAL(dp)                                      :: Hi_POI, Hb_POI, Hs_POI
!    REAL(dp), DIMENSION(C%nZ)                     :: Ti_POI
!        
!  ! Time log
!  ! ========
!    
!    filename = TRIM(C%output_dir) // 'aa_time_log_' // region%name // '.txt'
!    OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
!    
!    WRITE(UNIT = 1337, FMT = '(F10.1,I9,F11.3,F11.3,F11.3,F11.3,F11.3,F11.3,F11.3,F11.3)') region%time, region%mesh%nV, &
!      region%tcomp_total, region%tcomp_mesh, region%tcomp_SIA, region%tcomp_SSA, region%tcomp_thermo, region%tcomp_climate, region%tcomp_output
!    
!    CLOSE(UNIT = 1337)
!    
!  ! General output
!  ! ==============
!    
!    T2m_mean                   = 0._dp
!    total_snowfall             = 0._dp
!    total_rainfall             = 0._dp
!    total_melt                 = 0._dp
!    total_refreezing           = 0._dp
!    total_runoff               = 0._dp
!    total_SMB                  = 0._dp
!    total_BMB                  = 0._dp
!    total_MB                   = 0._dp
!    
!    DO vi = 1, region%mesh%nV
!      IF (region%ice%Hi(vi) > 0._dp) THEN
!        
!        total_BMB = total_BMB + (region%BMB%BMB(vi) * region%mesh%A(vi) / 1E9_dp)
!          
!        DO m = 1, 12
!          total_snowfall   = total_snowfall   + (region%SMB%Snowfall(  vi,m) * region%mesh%A(vi) / 1E9_dp)
!          total_rainfall   = total_rainfall   + (region%SMB%Rainfall(  vi,m) * region%mesh%A(vi) / 1E9_dp)
!          total_melt       = total_melt       + (region%SMB%Melt(      vi,m) * region%mesh%A(vi) / 1E9_dp)
!          total_refreezing = total_refreezing + (region%SMB%Refreezing(vi,m) * region%mesh%A(vi) / 1E9_dp)
!          total_runoff     = total_runoff     + (region%SMB%Runoff(    vi,m) * region%mesh%A(vi) / 1E9_dp)
!          total_SMB        = total_SMB        + (region%SMB%SMB(       vi,m) * region%mesh%A(vi) / 1E9_dp)
!        END DO
!        
!      END IF
!
!      T2m_mean = T2m_mean + SUM(region%climate%applied%T2m(vi,:)) * region%mesh%A(vi) / (12._dp * (region%mesh%xmax - region%mesh%xmin) * (region%mesh%ymax - region%mesh%ymin))
!      
!    END DO
!    
!    total_MB = total_SMB + total_BMB
!    
!    ! Average x-position of grounding line
!    x_GL = 0._dp
!    n_GL = 0
!    DO aci = 1, region%mesh%nAc
!      IF (region%ice%mask_gl_Ac( aci) == 1) THEN
!        n_GL = n_GL + 1
!        
!        vi = region%mesh%Aci( aci,1)
!        vj = region%mesh%Aci( aci,2)
!        
!        ! Find interpolated GL position
!        TAFi = region%ice%Hi( vi) - ((region%ice%SL( vi) - region%ice%Hb( vi)) * (seawater_density / ice_density))
!        TAFj = region%ice%Hi( vj) - ((region%ice%SL( vj) - region%ice%Hb( vj)) * (seawater_density / ice_density))
!        lambda = TAFi / (TAFi - TAFj)
!        
!        x_GL = x_GL + (NORM2(region%mesh%V( vi,:)) * (1._dp - lambda)) + (NORM2(region%mesh%V( vj,:)) * lambda)
!      END IF
!    END DO
!    x_GL = x_GL / n_GL
!    
!    filename = TRIM(C%output_dir) // 'aa_general_output_' // region%name // '.txt'
!    OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
!    
!    WRITE(UNIT = 1337, FMT = '(F10.1,2F8.2,F13.2,F8.2,9F11.2)') region%time, &
!      region%ice_volume, region%ice_volume_above_flotation, region%ice_area, T2m_mean, &
!      total_snowfall, total_rainfall, total_melt, total_refreezing, total_runoff, total_SMB, total_BMB, total_MB, x_GL
!    
!    CLOSE(UNIT = 1337)
!    
!  ! Point-of-interest output
!  ! ========================
!    
!    DO n = 1, region%mesh%nPOI
!    
!      vi1 = region%mesh%POI_vi(n,1)
!      vi2 = region%mesh%POI_vi(n,2)
!      vi3 = region%mesh%POI_vi(n,3)
!    
!      w1  = region%mesh%POI_w(n,1)
!      w2  = region%mesh%POI_w(n,2)
!      w3  = region%mesh%POI_w(n,3)
!      
!      Hi_POI = (region%ice%Hi(vi1  ) * w1) + (region%ice%Hi(vi2  ) * w2) + (region%ice%Hi(vi3  ) * w3)
!      Hb_POI = (region%ice%Hb(vi1  ) * w1) + (region%ice%Hb(vi2  ) * w2) + (region%ice%Hb(vi3  ) * w3)
!      Hs_POI = (region%ice%Hs(vi1  ) * w1) + (region%ice%Hs(vi2  ) * w2) + (region%ice%Hs(vi3  ) * w3)
!      Ti_POI = (region%ice%Ti(vi1,:) * w1) + (region%ice%Ti(vi2,:) * w2) + (region%ice%Ti(vi3,:) * w3)
!    
!      IF (n<10) THEN
!        WRITE(ns3,'(I1)') n
!        ns(1:2) = '00'
!        ns(3:3) = ns3
!      ELSEIF (n<100) THEN
!        WRITE(ns2,'(I2)') n
!        ns(1:1) = '0'
!        ns(2:3) = ns3
!      ELSE
!        WRITE(ns,'(I3)') n
!      END IF
!    
!      filename = TRIM(C%output_dir) // 'ab_POI_' // TRIM(region%name) // '_' // TRIM(ns) // '_data.txt'
!      OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
!    
!      WRITE(UNIT = 1337, FMT = '(F10.1,3F11.2)', ADVANCE='NO') region%time, Hi_POI, Hb_POI, Hs_POI
!      DO k = 1, C%nZ
!        WRITE(UNIT = 1337, FMT = '(F11.2)', ADVANCE='NO') Ti_POI(k)
!      END DO
!      WRITE(UNIT = 1337, FMT = '(A)') ''
!      
!      CLOSE(UNIT = 1337)
!    
!    END DO
!    
!  END SUBROUTINE write_text_output
  
! == Test the matrix operators
  SUBROUTINE test_matrix_operators_mesh( mesh)
    ! Test all the matrix operators
    
    USE mesh_operators_module
    USE utilities_module
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, ti, aci
    REAL(dp)                                           :: x, y, d, ddx, ddy
    
    ! Exact solutions
    REAL(dp), DIMENSION(:    ), POINTER                :: d_a_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: d_c_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_a_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_c_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_a_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_b_ex
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_c_ex
    INTEGER :: wd_a_ex, wd_b_ex, wd_c_ex, wddx_a_ex, wddx_b_ex, wddx_c_ex, wddy_a_ex, wddy_b_ex, wddy_c_ex
    
    ! Discretised approximations
    REAL(dp), DIMENSION(:    ), POINTER                :: d_aa
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ab
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ba
    REAL(dp), DIMENSION(:    ), POINTER                :: d_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: d_bc
    REAL(dp), DIMENSION(:    ), POINTER                :: d_ca
    REAL(dp), DIMENSION(:    ), POINTER                :: d_cb
    REAL(dp), DIMENSION(:    ), POINTER                :: d_cc
    INTEGER :: wd_aa, wd_ab, wd_ac, wd_ba, wd_bb, wd_bc, wd_ca, wd_cb, wd_cc
    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_aa
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ab
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ba
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_bc
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_ca
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_cb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddx_cc
    INTEGER :: wddx_aa, wddx_ab, wddx_ac, wddx_ba, wddx_bb, wddx_bc, wddx_ca, wddx_cb, wddx_cc
    
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_aa
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ab
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ac
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ba
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_bb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_bc
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_ca
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_cb
    REAL(dp), DIMENSION(:    ), POINTER                :: ddy_cc
    INTEGER :: wddy_aa, wddy_ab, wddy_ac, wddy_ba, wddy_bb, wddy_bc, wddy_ca, wddy_cb, wddy_cc
        
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' =================================='
    IF (par%master) WRITE(0,*) ' == Testing the matrix operators =='
    IF (par%master) WRITE(0,*) ' =================================='
    IF (par%master) WRITE(0,*) ''
    CALL sync
    
  ! Write all the matrix operators to files, for double-checking stuff in Matlab
  ! ============================================================================
    
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_a_b, 'M_map_a_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_a_c, 'M_map_a_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_b_a, 'M_map_b_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_b_c, 'M_map_b_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_c_a, 'M_map_c_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_map_c_b, 'M_map_c_b.nc')
    
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_a_a, 'M_ddx_a_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_a_b, 'M_ddx_a_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_a_c, 'M_ddx_a_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_b_a, 'M_ddx_b_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_b_b, 'M_ddx_b_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_b_c, 'M_ddx_b_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_c_a, 'M_ddx_c_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_c_b, 'M_ddx_c_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddx_c_c, 'M_ddx_c_c.nc')
    
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_a_a, 'M_ddy_a_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_a_b, 'M_ddy_a_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_a_c, 'M_ddy_a_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_b_a, 'M_ddy_b_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_b_b, 'M_ddy_b_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_b_c, 'M_ddy_b_c.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_c_a, 'M_ddy_c_a.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_c_b, 'M_ddy_c_b.nc')
    CALL write_CSR_matrix_to_NetCDF( mesh%M_ddy_c_c, 'M_ddy_c_c.nc')
    
  ! Allocate shared memory
  ! ======================
    
    ! Exact solutions
    CALL allocate_shared_dp_1D( mesh%nV,   d_a_ex,   wd_a_ex  )
    CALL allocate_shared_dp_1D( mesh%nTri, d_b_ex,   wd_b_ex  )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_c_ex,   wd_c_ex  )
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_a_ex, wddx_a_ex)
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_b_ex, wddx_b_ex)
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_c_ex, wddx_c_ex)
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_a_ex, wddy_a_ex)
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_b_ex, wddy_b_ex)
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_c_ex, wddy_c_ex)
    
    ! Discretised approximations
    CALL allocate_shared_dp_1D( mesh%nV,   d_aa,     wd_aa    )
    CALL allocate_shared_dp_1D( mesh%nTri, d_ab,     wd_ab    )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_ac,     wd_ac    )
    CALL allocate_shared_dp_1D( mesh%nV,   d_ba,     wd_ba    )
    CALL allocate_shared_dp_1D( mesh%nTri, d_bb,     wd_bb    )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_bc,     wd_bc    )
    CALL allocate_shared_dp_1D( mesh%nV,   d_ca,     wd_ca    )
    CALL allocate_shared_dp_1D( mesh%nTri, d_cb,     wd_cb    )
    CALL allocate_shared_dp_1D( mesh%nAc,  d_cc,     wd_cc    )
    
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_aa,   wddx_aa  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_ab,   wddx_ab  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_ac,   wddx_ac  )
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_ba,   wddx_ba  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_bb,   wddx_bb  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_bc,   wddx_bc  )
    CALL allocate_shared_dp_1D( mesh%nV,   ddx_ca,   wddx_ca  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddx_cb,   wddx_cb  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddx_cc,   wddx_cc  )
    
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_aa,   wddy_aa  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_ab,   wddy_ab  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_ac,   wddy_ac  )
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_ba,   wddy_ba  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_bb,   wddy_bb  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_bc,   wddy_bc  )
    CALL allocate_shared_dp_1D( mesh%nV,   ddy_ca,   wddy_ca  )
    CALL allocate_shared_dp_1D( mesh%nTri, ddy_cb,   wddy_cb  )
    CALL allocate_shared_dp_1D( mesh%nAc,  ddy_cc,   wddy_cc  )
    
  ! Calculate exact solutions
  ! =========================
    
    DO vi = mesh%v1, mesh%v2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      CALL test_matrix_operators_mesh_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy)
      d_a_ex(   vi) = d
      ddx_a_ex( vi) = ddx
      ddy_a_ex( vi) = ddy
    END DO
    DO ti = mesh%t1, mesh%t2
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      CALL test_matrix_operators_mesh_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy)
      d_b_ex(   ti) = d
      ddx_b_ex( ti) = ddx
      ddy_b_ex( ti) = ddy
    END DO
    DO aci = mesh%ac1, mesh%ac2
      x = mesh%VAc( aci,1)
      y = mesh%VAc( aci,2)
      CALL test_matrix_operators_mesh_test_function( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, d, ddx, ddy)
      d_c_ex(   aci) = d
      ddx_c_ex( aci) = ddx
      ddy_c_ex( aci) = ddy
    END DO
    CALL sync
    
  ! Calculate discretised approximations
  ! ====================================
  
    ! Mapping
    d_aa( mesh%v1:mesh%v2) = d_a_ex( mesh%v1:mesh%v2)
    CALL map_a_to_b_2D( mesh, d_a_ex, d_ab)
    CALL map_a_to_c_2D( mesh, d_a_ex, d_ac)
    
    CALL map_b_to_a_2D( mesh, d_b_ex, d_ba)
    d_bb( mesh%t1:mesh%t2) = d_b_ex( mesh%t1:mesh%t2)
    CALL map_b_to_c_2D( mesh, d_b_ex, d_bc)
    
    CALL map_c_to_a_2D( mesh, d_c_ex, d_ca)
    CALL map_c_to_b_2D( mesh, d_c_ex, d_cb)
    d_cc( mesh%ac1:mesh%ac2) = d_c_ex( mesh%ac1:mesh%ac2)
    
    ! d/dx
    CALL ddx_a_to_a_2D( mesh, d_a_ex, ddx_aa)
    CALL ddx_a_to_b_2D( mesh, d_a_ex, ddx_ab)
    CALL ddx_a_to_c_2D( mesh, d_a_ex, ddx_ac)
    
    CALL ddx_b_to_a_2D( mesh, d_b_ex, ddx_ba)
    CALL ddx_b_to_b_2D( mesh, d_b_ex, ddx_bb)
    CALL ddx_b_to_c_2D( mesh, d_b_ex, ddx_bc)
    
    CALL ddx_c_to_a_2D( mesh, d_c_ex, ddx_ca)
    CALL ddx_c_to_b_2D( mesh, d_c_ex, ddx_cb)
    CALL ddx_c_to_c_2D( mesh, d_c_ex, ddx_cc)
    
    ! d/dy
    CALL ddy_a_to_a_2D( mesh, d_a_ex, ddy_aa)
    CALL ddy_a_to_b_2D( mesh, d_a_ex, ddy_ab)
    CALL ddy_a_to_c_2D( mesh, d_a_ex, ddy_ac)
    
    CALL ddy_b_to_a_2D( mesh, d_b_ex, ddy_ba)
    CALL ddy_b_to_b_2D( mesh, d_b_ex, ddy_bb)
    CALL ddy_b_to_c_2D( mesh, d_b_ex, ddy_bc)
    
    CALL ddy_c_to_a_2D( mesh, d_c_ex, ddy_ca)
    CALL ddy_c_to_b_2D( mesh, d_c_ex, ddy_cb)
    CALL ddy_c_to_c_2D( mesh, d_c_ex, ddy_cc)
    
  ! Write to debug file
  ! ===================
  
    IF (par%master) THEN
    
      ! a-grid (vertex)
      debug%dp_2D_a_01 = d_aa   - d_a_ex
      debug%dp_2D_a_02 = d_ba   - d_a_ex
      debug%dp_2D_a_03 = d_ca   - d_a_ex
      debug%dp_2D_a_04 = ddx_aa - ddx_a_ex
      debug%dp_2D_a_05 = ddx_ba - ddx_a_ex
      debug%dp_2D_a_06 = ddx_ca - ddx_a_ex
      debug%dp_2D_a_07 = ddy_aa - ddy_a_ex
      debug%dp_2D_a_08 = ddy_ba - ddy_a_ex
      debug%dp_2D_a_09 = ddy_ca - ddy_a_ex
    
      ! b-grid (triangle)
      debug%dp_2D_b_01 = d_ab   - d_b_ex
      debug%dp_2D_b_02 = d_bb   - d_b_ex
      debug%dp_2D_b_03 = d_cb   - d_b_ex
      debug%dp_2D_b_04 = ddx_ab - ddx_b_ex
      debug%dp_2D_b_05 = ddx_bb - ddx_b_ex
      debug%dp_2D_b_06 = ddx_cb - ddx_b_ex
      debug%dp_2D_b_07 = ddy_ab - ddy_b_ex
      debug%dp_2D_b_08 = ddy_bb - ddy_b_ex
      debug%dp_2D_b_09 = ddy_cb - ddy_b_ex
    
      ! c-grid (edge)
      debug%dp_2D_c_01 = d_ac   - d_c_ex
      debug%dp_2D_c_02 = d_bc   - d_c_ex
      debug%dp_2D_c_03 = d_cc   - d_c_ex
      debug%dp_2D_c_04 = ddx_ac - ddx_c_ex
      debug%dp_2D_c_05 = ddx_bc - ddx_c_ex
      debug%dp_2D_c_06 = ddx_cc - ddx_c_ex
      debug%dp_2D_c_07 = ddy_ac - ddy_c_ex
      debug%dp_2D_c_08 = ddy_bc - ddy_c_ex
      debug%dp_2D_c_09 = ddy_cc - ddy_c_ex
      
      CALL write_to_debug_file
      
    END IF
    CALL sync
    
  ! Clean up after yourself
  ! =======================
    
    CALL deallocate_shared( wd_a_ex  )
    CALL deallocate_shared( wd_b_ex  )
    CALL deallocate_shared( wd_c_ex  )
    CALL deallocate_shared( wddx_a_ex)
    CALL deallocate_shared( wddx_b_ex)
    CALL deallocate_shared( wddx_c_ex)
    CALL deallocate_shared( wddy_a_ex)
    CALL deallocate_shared( wddy_b_ex)
    CALL deallocate_shared( wddy_c_ex)
    
    CALL deallocate_shared( wd_aa    )
    CALL deallocate_shared( wd_ab    )
    CALL deallocate_shared( wd_ac    )
    CALL deallocate_shared( wd_ba    )
    CALL deallocate_shared( wd_bb    )
    CALL deallocate_shared( wd_bc    )
    CALL deallocate_shared( wd_ca    )
    CALL deallocate_shared( wd_cb    )
    CALL deallocate_shared( wd_cc    )
    
    CALL deallocate_shared( wddx_aa  )
    CALL deallocate_shared( wddx_ab  )
    CALL deallocate_shared( wddx_ac  )
    CALL deallocate_shared( wddx_ba  )
    CALL deallocate_shared( wddx_bb  )
    CALL deallocate_shared( wddx_bc  )
    CALL deallocate_shared( wddx_ca  )
    CALL deallocate_shared( wddx_cb  )
    CALL deallocate_shared( wddx_cc  )
    
    CALL deallocate_shared( wddy_aa  )
    CALL deallocate_shared( wddy_ab  )
    CALL deallocate_shared( wddy_ac  )
    CALL deallocate_shared( wddy_ba  )
    CALL deallocate_shared( wddy_bb  )
    CALL deallocate_shared( wddy_bc  )
    CALL deallocate_shared( wddy_ca  )
    CALL deallocate_shared( wddy_cb  )
    CALL deallocate_shared( wddy_cc  )
        
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' ==========================================='
    IF (par%master) WRITE(0,*) ' == Finished testing the matrix operators =='
    IF (par%master) WRITE(0,*) ' ==========================================='
    IF (par%master) WRITE(0,*) ''
    CALL sync
    
  END SUBROUTINE test_matrix_operators_mesh
  SUBROUTINE test_matrix_operators_mesh_test_function( x, y, xmin, xmax, ymin, ymax, d, ddx, ddy)
    ! Simple function for testing all the matrix operators
      
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: x
    REAL(dp),                            INTENT(IN)    :: y
    REAL(dp),                            INTENT(IN)    :: xmin
    REAL(dp),                            INTENT(IN)    :: xmax
    REAL(dp),                            INTENT(IN)    :: ymin
    REAL(dp),                            INTENT(IN)    :: ymax
    REAL(dp),                            INTENT(OUT)   :: d
    REAL(dp),                            INTENT(OUT)   :: ddx
    REAL(dp),                            INTENT(OUT)   :: ddy
    
    ! Local variables:
    REAL(dp)                                           :: xp, yp
    
    xp = 2._dp * pi * (x - xmin) / (xmax - xmin)
    yp = 2._dp * pi * (y - ymin) / (ymax - ymin)
    
    d   = SIN( xp) * SIN( yp)
    ddx = COS( xp) * SIN( yp) * 2._dp * pi / (xmax - xmin)
    ddy = SIN( xp) * COS( yp) * 2._dp * pi / (xmax - xmin)
        
  END SUBROUTINE test_matrix_operators_mesh_test_function
  SUBROUTINE solve_modified_Laplace_equation( mesh)
    ! Test the discretisation by solving (a modified version of) the Laplace equation:
    !
    ! d/dx ( N * df/dx) + d/dx ( N * df/dy) = 0
    !
    ! (Modified to include a spatially variable stiffness, making it much more similar to
    ! the SSA/DIVA, and necessitating the use of a staggered grid/mesh to solve it)
    
    USE mesh_operators_module
    USE utilities_module
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, ti
    REAL(dp), DIMENSION(:    ), POINTER                ::  f_a,  N_b,  b_a
    INTEGER                                            :: wf_a, wN_b, wb_a
    TYPE(type_sparse_matrix_CSR)                       :: M_a, BC_a, mN_b, mNddx_b
    REAL(dp)                                           :: x, y, xp, yp
    REAL(dp), PARAMETER                                :: amp = 5._dp
    INTEGER                                            :: nnz_BC
        
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' ==========================================='
    IF (par%master) WRITE(0,*) ' == Solving the modified Laplace equation =='
    IF (par%master) WRITE(0,*) ' ==========================================='
    IF (par%master) WRITE(0,*) ''
    CALL sync
    
  ! Solve the equation on the a-grid (vertex)
  ! =========================================
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV,   f_a, wf_a)
    CALL allocate_shared_dp_1D( mesh%nV,   b_a, wb_a)
    CALL allocate_shared_dp_1D( mesh%nTri, N_b, wN_b)
    
    ! Set up the spatially variable stiffness N on the b-grid (triangle)
    DO ti = mesh%t1, mesh%t2
      x = mesh%TriGC( ti,1)
      y = mesh%TriGC( ti,2)
      xp = (x - mesh%xmin) / (mesh%xmax - mesh%xmin)
      yp = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
      N_b( ti) = 1._dp + EXP(amp * SIN( 2._dp * pi * xp) * SIN(2._dp * pi * yp))
    END DO
    CALL sync
    
    ! Convert stiffness from vector to diagonal matrix
    CALL convert_vector_to_diag_CSR( N_b, mN_b)
    
    ! Create the operator matrix M
    CALL multiply_matrix_matrix_CSR( mN_b          , mesh%M_ddx_a_b, mNddx_b)
    CALL multiply_matrix_matrix_CSR( mesh%M_ddx_b_a, mNddx_b       , M_a    )
    
    ! Create the boundary conditions matrix BC
    
    ! Find maximum number of non-zero entries
    nnz_BC = 0
    DO vi = mesh%v1, mesh%v2
      IF (mesh%edge_index( vi) > 0) THEN
        nnz_BC = nnz_BC + 1 + mesh%nC( vi)
      END IF
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nnz_BC, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)   
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( BC_a, mesh%nV, mesh%nV, nnz_BC)
    
    ! Fill in values
    BC_a%nnz = 0
    BC_a%ptr = 1
    DO vi = mesh%v1, mesh%v2
    
      IF (mesh%edge_index( vi) == 0) THEN
        ! Free vertex: no boundary conditions apply
      
      ELSEIF (mesh%edge_index( vi) == 6 .OR. mesh%edge_index( vi) == 7 .OR. mesh%edge_index( vi) == 8) THEN
        ! West: bump
        
        ! Matrix
        BC_a%nnz  = BC_a%nnz+1
        BC_a%index( BC_a%nnz) = vi
        BC_a%val(   BC_a%nnz) = 1._dp
        
        ! Right-hand side vector
        y = mesh%V( vi,2)
        y = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
        b_a( vi) = 0.5_dp * (1._dp - COS( 2._dp * pi * y))
        
      ELSE
        ! Otherwise: zero
        
        ! Matrix
        BC_a%nnz  = BC_a%nnz+1
        BC_a%index( BC_a%nnz) = vi
        BC_a%val(   BC_a%nnz) = 1._dp
        
        ! Right-hand side vector
        b_a( vi) = 0._dp
        
      END IF
      
      ! Finalise matrix row
      BC_a%ptr( vi+1) = BC_a%nnz + 1
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! Finalise matrix
    CALL finalise_matrix_CSR_dist( BC_a, mesh%v1, mesh%v2)
    
    ! Add boundary conditions to the operator matrix
    CALL overwrite_rows_CSR( M_a, BC_a)
    
    ! Solve the equation
    CALL solve_matrix_equation_CSR( M_a, b_a, f_a, C%DIVA_choice_matrix_solver, &
      SOR_nit = 5000, SOR_tol = 0.0001_dp, SOR_omega = 1.3_dp, &
      PETSc_rtol = C%DIVA_PETSc_rtol, PETSc_abstol = C%DIVA_PETSc_abstol)
      
    ! Write result to debug file
    IF (par%master) THEN
      debug%dp_2D_a_01 = f_a
      CALL write_to_debug_file
    END IF
    CALL sync
        
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' ===================================================='
    IF (par%master) WRITE(0,*) ' == Finished solving the modified Laplace equation =='
    IF (par%master) WRITE(0,*) ' ===================================================='
    IF (par%master) WRITE(0,*) ''
    CALL sync
    
  END SUBROUTINE solve_modified_Laplace_equation
  SUBROUTINE solve_modified_Laplace_equation_c( mesh)
    ! Test the discretisation by solving (a modified version of) the Laplace equation:
    !
    ! d/dx ( N * df/dx) + d/dx ( N * df/dy) = 0
    !
    ! (Modified to include a spatially variable stiffness, making it much more similar to
    ! the SSA/DIVA, and necessitating the use of a staggered grid/mesh to solve it)
    
    USE mesh_operators_module
    USE utilities_module
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    
    ! Local variables:
    INTEGER                                            :: vi, aci
    REAL(dp), DIMENSION(:    ), POINTER                ::  f_c,  N_a,  b_c
    INTEGER                                            :: wf_c, wN_a, wb_c
    TYPE(type_sparse_matrix_CSR)                       :: M_c, BC_c, mN_a, mNddx_a
    REAL(dp)                                           :: x, y, xp, yp
    REAL(dp), PARAMETER                                :: amp = 5._dp
    INTEGER                                            :: nnz_BC
        
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' ======================================================='
    IF (par%master) WRITE(0,*) ' == Solving the modified Laplace equation (STAGGERED) =='
    IF (par%master) WRITE(0,*) ' ======================================================='
    IF (par%master) WRITE(0,*) ''
    CALL sync
    
  ! Solve the equation on the a-grid (vertex)
  ! =========================================
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nAc, f_c, wf_c)
    CALL allocate_shared_dp_1D( mesh%nAc, b_c, wb_c)
    CALL allocate_shared_dp_1D( mesh%nV,  N_a, wN_a)
    
    ! Set up the spatially variable stiffness N on the a-grid (vertex)
    DO vi = mesh%v1, mesh%v2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      xp = (x - mesh%xmin) / (mesh%xmax - mesh%xmin)
      yp = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
      N_a( vi) = 1._dp + EXP(amp * SIN( 2._dp * pi * xp) * SIN(2._dp * pi * yp))
    END DO
    CALL sync
    
    ! Convert stiffness from vector to diagonal matrix
    CALL convert_vector_to_diag_CSR( N_a, mN_a)
    
    ! Create the operator matrix M
    CALL multiply_matrix_matrix_CSR( mN_a          , mesh%M_ddx_c_a, mNddx_a)
    CALL multiply_matrix_matrix_CSR( mesh%M_ddx_a_c, mNddx_a       , M_c    )
    
    ! Create the boundary conditions matrix BC
    
    ! Find maximum number of non-zero entries
    nnz_BC = 0
    DO aci = mesh%ac1, mesh%ac2
      IF (mesh%edge_index_Ac( aci) > 0) THEN
        nnz_BC = nnz_BC + 1 + mesh%nC_mem
      END IF
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, nnz_BC, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)   
    
    ! Allocate distributed shared memory
    CALL allocate_matrix_CSR_dist( BC_c, mesh%nAc, mesh%nAc, nnz_BC)
    
    ! Fill in values
    BC_c%nnz = 0
    BC_c%ptr = 1
    DO aci = mesh%ac1, mesh%ac2
    
      IF (mesh%edge_index_Ac( aci) == 0) THEN
        ! Free vertex: no boundary conditions apply
      
      ELSEIF (mesh%edge_index_Ac( aci) == 6 .OR. mesh%edge_index_Ac( aci) == 7 .OR. mesh%edge_index_Ac( aci) == 8) THEN
        ! West: bump
        
        ! Matrix
        BC_c%nnz  = BC_c%nnz+1
        BC_c%index( BC_c%nnz) = aci
        BC_c%val(   BC_c%nnz) = 1._dp
        
        ! Right-hand side vector
        y = mesh%VAc( aci,2)
        y = (y - mesh%ymin) / (mesh%ymax - mesh%ymin)
        b_c( aci) = 0.5_dp * (1._dp - COS( 2._dp * pi * y))
        
      ELSE
        ! Otherwise: zero
        
        ! Matrix
        BC_c%nnz  = BC_c%nnz+1
        BC_c%index( BC_c%nnz) = aci
        BC_c%val(   BC_c%nnz) = 1._dp
        
        ! Right-hand side vector
        b_c( aci) = 0._dp
        
      END IF
      
      ! Finalise matrix row
      BC_c%ptr( aci+1) = BC_c%nnz + 1
      
    END DO ! DO aci = mesh%ac1, mesh%ac2
    CALL sync
    
    ! Finalise matrix
    CALL finalise_matrix_CSR_dist( BC_c, mesh%ac1, mesh%ac2)
    
    ! Add boundary conditions to the operator matrix
    CALL overwrite_rows_CSR( M_c, BC_c)
    
    ! Solve the equation
    CALL solve_matrix_equation_CSR( M_c, b_c, f_c, C%DIVA_choice_matrix_solver, &
      SOR_nit = 1, SOR_tol = 0.0001_dp, SOR_omega = 1.3_dp, &
      PETSc_rtol = C%DIVA_PETSc_rtol, PETSc_abstol = C%DIVA_PETSc_abstol)
      
    ! Write result to debug file
    IF (par%master) THEN
      debug%dp_2D_c_01 = f_c
      CALL write_to_debug_file
    END IF
    CALL sync
        
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' ================================================================'
    IF (par%master) WRITE(0,*) ' == Finished solving the modified Laplace equation (STAGGERED) =='
    IF (par%master) WRITE(0,*) ' ================================================================'
    IF (par%master) WRITE(0,*) ''
    CALL sync
    
  END SUBROUTINE solve_modified_Laplace_equation_c

END MODULE UFEMISM_main_model
