MODULE UFEMISM_main_model

  ! The main regional ice-sheet model

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parameters_module
  USE petsc_module,                    ONLY: perr
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
  USE data_types_module,               ONLY: type_model_region, type_mesh, type_grid, type_climate_matrix, type_remapping_mesh_mesh
  USE reference_fields_module,         ONLY: initialise_reference_geometries, map_reference_geometries_to_mesh
  USE mesh_memory_module,              ONLY: deallocate_mesh_all
  USE mesh_help_functions_module,      ONLY: inverse_oblique_sg_projection
  USE mesh_creation_module,            ONLY: create_mesh_from_cart_data
  USE mesh_mapping_module,             ONLY: calc_remapping_operators_mesh_mesh, deallocate_remapping_operators_mesh_mesh, &
                                             calc_remapping_operator_mesh2grid, deallocate_remapping_operators_mesh2grid, &
                                             calc_remapping_operator_grid2mesh, deallocate_remapping_operators_grid2mesh
  USE mesh_update_module,              ONLY: determine_mesh_fitness, create_new_mesh
  USE netcdf_module,                   ONLY: create_output_files, write_to_output_files, initialise_debug_fields, &
                                             associate_debug_fields, reallocate_debug_fields, create_debug_file, write_PETSc_matrix_to_NetCDF
  USE restart_module,                  ONLY: read_mesh_from_restart_file, read_init_data_from_restart_file
  USE general_ice_model_data_module,   ONLY: update_general_ice_model_data, initialise_mask_noice
  USE ice_dynamics_module,             ONLY: initialise_ice_model,      remap_ice_model,      run_ice_model, update_ice_thickness
  USE thermodynamics_module,           ONLY: initialise_ice_temperature,                      run_thermo_model, calc_ice_rheology
  USE climate_module,                  ONLY: initialise_climate_model,  remap_climate_model,  run_climate_model, map_subclimate_to_mesh
  USE SMB_module,                      ONLY: initialise_SMB_model,      remap_SMB_model,      run_SMB_model
  USE BMB_module,                      ONLY: initialise_BMB_model,      remap_BMB_model,      run_BMB_model
  USE isotopes_module,                 ONLY: initialise_isotopes_model, remap_isotopes_model, run_isotopes_model, calculate_reference_isotopes
  USE bedrock_ELRA_module,             ONLY: initialise_ELRA_model,     remap_ELRA_model,     run_ELRA_model
  USE tests_and_checks_module,         ONLY: run_all_matrix_tests

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
      region%tcomp_GIA = region%tcomp_GIA + t2 - t1
      
    ! Mesh update
    ! ===========
      
      ! Check if the mesh needs to be updated
      t2 = MPI_WTIME()
      meshfitness = 1._dp
      IF (region%time > region%t_last_mesh + C%dt_mesh_min) THEN
        CALL determine_mesh_fitness(region%mesh, region%ice, meshfitness)    
      END IF
      region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2 
    
      ! If required, update the mesh
      IF (meshfitness < C%mesh_fitness_threshold) THEN
    !  IF (.FALSE.) THEN
    !  IF (.TRUE.) THEN 
        region%t_last_mesh = region%time
        t2 = MPI_WTIME()
        CALL run_model_update_mesh( region, matrix)
        region%tcomp_mesh = region%tcomp_mesh + MPI_WTIME() - t2
      END IF
      
    ! Ice dynamics
    ! ============
    
      ! Calculate ice velocities and the resulting change in ice geometry
      ! NOTE: geometry is not updated yet; this happens at the end of the time loop
      t1 = MPI_WTIME()
      CALL run_ice_model( region, t_end)
      t2 = MPI_WTIME()
      region%tcomp_ice = region%tcomp_ice + t2 - t1
      
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
      region%tcomp_climate = region%tcomp_climate + t2 - t1
      
    ! Thermodynamics
    ! ==============
    
      t1 = MPI_WTIME()
      CALL run_thermo_model( region%mesh, region%ice, region%climate%applied, region%SMB, region%time, do_solve_heat_equation = region%do_thermo)
      t2 = MPI_WTIME()
      region%tcomp_thermo = region%tcomp_thermo + t2 - t1
      
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
      region%time = region%time + region%dt
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
    TYPE(type_remapping_mesh_mesh)                     :: map
    
    ! Create a new mesh
    CALL create_new_mesh( region)
    
    ! Update the mapping operators between the new mesh and the fixed square grids
    CALL deallocate_remapping_operators_mesh2grid(            region%grid_output)
    CALL deallocate_remapping_operators_mesh2grid(            region%grid_GIA   )
    CALL deallocate_remapping_operators_mesh2grid(            region%grid_smooth)
    
    CALL deallocate_remapping_operators_grid2mesh(            region%grid_output)
    CALL deallocate_remapping_operators_grid2mesh(            region%grid_GIA   )
    CALL deallocate_remapping_operators_grid2mesh(            region%grid_smooth)
    
    CALL calc_remapping_operator_mesh2grid( region%mesh_new, region%grid_output)
    CALL calc_remapping_operator_mesh2grid( region%mesh_new, region%grid_GIA   )
    CALL calc_remapping_operator_mesh2grid( region%mesh_new, region%grid_smooth)
    
    CALL calc_remapping_operator_grid2mesh( region%grid_output, region%mesh_new)
    CALL calc_remapping_operator_grid2mesh( region%grid_GIA   , region%mesh_new)
    CALL calc_remapping_operator_grid2mesh( region%grid_smooth, region%mesh_new)
    
    ! Calculate the mapping arrays
    CALL calc_remapping_operators_mesh_mesh( region%mesh, region%mesh_new, map)
    
    ! Recalculate the "no ice" mask (also needed for the climate model)
    deallocate(region%mask_noice)
    allocate(region%mask_noice(region%mesh_new%vi1:region%mesh_new%vi2))

    CALL initialise_mask_noice(  region, region%mesh_new)
            
    ! Reallocate and remap reference geometries
    CALL deallocate_shared(  region%refgeo_init%wHi )
    CALL deallocate_shared(  region%refgeo_init%wHb )
    CALL deallocate_shared(  region%refgeo_init%wHs )
    CALL deallocate_shared(  region%refgeo_PD%wHi   )
    CALL deallocate_shared(  region%refgeo_PD%wHb   )
    CALL deallocate_shared(  region%refgeo_PD%wHs   )
    CALL deallocate_shared(  region%refgeo_GIAeq%wHi)
    CALL deallocate_shared(  region%refgeo_GIAeq%wHb)
    CALL deallocate_shared(  region%refgeo_GIAeq%wHs)
    CALL map_reference_geometries_to_mesh( region, region%mesh_new)
    
    ! Remap all the submodels
    CALL remap_ELRA_model(     region%mesh, region%mesh_new, map, region%ice, region%refgeo_PD, region%grid_GIA)
    CALL remap_ice_model(      region%mesh, region%mesh_new, map, region%ice, region%refgeo_PD, region%time)
    CALL remap_climate_model(  region%mesh, region%mesh_new, map, region%climate, matrix, region%refgeo_PD, region%grid_smooth, region%mask_noice, region%name)
    CALL remap_SMB_model(      region%mesh, region%mesh_new, map, region%SMB)
    CALL remap_BMB_model(      region%mesh, region%mesh_new, map, region%BMB)
    CALL remap_isotopes_model( region%mesh, region%mesh_new, map, region)
            
    ! Deallocate shared memory for the mapping arrays
    CALL deallocate_remapping_operators_mesh_mesh( map)
    
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
    region%t_next_SIA     = region%time
    region%t_next_SSA     = region%time
    region%t_next_DIVA    = region%time
    region%t_next_thermo  = region%time
    region%t_next_climate = region%time
    region%t_next_SMB     = region%time
    region%t_next_BMB     = region%time
    region%t_next_ELRA    = region%time
    
    region%do_SIA         = .TRUE.
    region%do_SSA         = .TRUE.
    region%do_DIVA        = .TRUE.
    region%do_thermo      = .TRUE.
    region%do_climate     = .TRUE.
    region%do_SMB         = .TRUE.
    region%do_BMB         = .TRUE.
    region%do_ELRA        = .TRUE.
      
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
    
    CALL initialise_reference_geometries( region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq, region%name)
    
    ! ===== The mesh =====
    ! ====================
    
    CALL create_mesh_from_cart_data( region)
    
    ! ===== Map reference geometries to the mesh =====
    ! ================================================
    
    CALL map_reference_geometries_to_mesh( region, region%mesh)
    
    ! The different square grids
    ! ==========================
    
    IF (par%master) WRITE(0,*) '  Initialising square grids for output, GIA, and data smoothing...'
    CALL initialise_model_square_grid( region, region%grid_output, C%dx_grid_output)
    CALL initialise_model_square_grid( region, region%grid_GIA,    C%dx_grid_GIA   )
    CALL initialise_model_square_grid( region, region%grid_smooth, C%dx_grid_smooth)
    
    ! ===== Initialise dummy fields for debugging
    ! ===========================================
    
    IF (par%master) WRITE(0,*) '  Initialising debug fields...'
    CALL initialise_debug_fields( region)
       
    ! ===== Output files =====
    ! =======================
  
    CALL create_output_files(    region)
    CALL associate_debug_fields( region)
    region%output_file_exists = .TRUE.
    CALL sync
    
    ! ===== The "no ice" mask
    ! =======================

    allocate(region%mask_noice(region%mesh%vi1:region%mesh%vi2))
    CALL initialise_mask_noice( region, region%mesh)
  
    ! ===== The ice dynamics model
    ! ============================
    
    CALL initialise_ice_model( region%mesh, region%ice, region%refgeo_init)
        
    ! ===== The climate model =====
    ! =============================    
    
    CALL initialise_climate_model( region%climate, matrix, region%mesh, region%refgeo_PD, region%name, region%mask_noice, region%grid_smooth)
    
    ! ===== The SMB model =====
    ! =========================    
    
    CALL initialise_SMB_model( region%mesh, region%ice, region%SMB, region%name)
    
    ! ===== The BMB model =====
    ! =========================    
    
    CALL initialise_BMB_model( region%mesh, region%BMB, region%name) 
    
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
    
    IF     (C%choice_GIA_model == 'none') THEN
      ! Nothing to be done
    ELSEIF (C%choice_GIA_model == 'ELRA') THEN
      CALL initialise_ELRA_model( region%mesh, region%grid_GIA, region%ice, region%refgeo_PD)
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
    INTEGER                                       :: nsx, nsy, i, j, n
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp
    
    n1 = par%mem%n
    
    nsx = 0
    nsy = 0
    
    CALL allocate_shared_int_0D(                   grid%nx          , grid%wnx          )
    CALL allocate_shared_int_0D(                   grid%ny          , grid%wny          )
    CALL allocate_shared_dp_0D(                    grid%dx          , grid%wdx          )
    CALL allocate_shared_dp_0D(                    grid%xmin        , grid%wxmin        )
    CALL allocate_shared_dp_0D(                    grid%xmax        , grid%wxmax        )
    CALL allocate_shared_dp_0D(                    grid%ymin        , grid%wymin        )
    CALL allocate_shared_dp_0D(                    grid%ymax        , grid%wymax        )
    
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
      
    ! Tolerance; points lying within this distance of each other are treated as identical
    CALL allocate_shared_dp_0D( grid%tol_dist, grid%wtol_dist)
    IF (par%master) grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp
    CALL sync
    
    ! Set up grid-to-vector translation tables
    CALL allocate_shared_int_0D(                   grid%n           , grid%wn           )
    IF (par%master) grid%n  = grid%nx * grid%ny
    CALL sync
    CALL allocate_shared_int_2D( grid%nx, grid%ny, grid%ij2n        , grid%wij2n        )
    CALL allocate_shared_int_2D( grid%n , 2      , grid%n2ij        , grid%wn2ij        )
    IF (par%master) THEN
      n = 0
      DO i = 1, grid%nx
        IF (MOD(i,2) == 1) THEN
          DO j = 1, grid%ny
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = grid%ny, 1, -1
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO
    END IF
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
    CALL calc_remapping_operator_mesh2grid( region%mesh, grid)
    CALL calc_remapping_operator_grid2mesh( grid, region%mesh)
    
    n2 = par%mem%n
    CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE initialise_model_square_grid
  
  ! Calculate this region's ice sheet's volume and area
  SUBROUTINE calculate_icesheet_volume_and_area( region)
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    INTEGER                                       :: vi
    REAL(dp)                                      :: ice_area, ice_volume, thickness_above_flotation, ice_volume_above_flotation
    REAL(dp), dimension(3)                        :: buffer
    
    ice_area                   = 0._dp
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp
    
    ! Calculate ice area and volume for processor domain
    DO vi = region%mesh%vi1, region%mesh%vi2
    
      IF (region%ice%mask_ice_a( vi) == 1) THEN
        ice_volume = ice_volume + (region%ice%Hi_a( vi) * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area))
        ice_area   = ice_area   + region%mesh%A( vi) * 1.0E-06_dp ! [km^3]
    
        ! Thickness above flotation
        thickness_above_flotation = MAX(0._dp, region%ice%Hi_a( vi) - MAX(0._dp, (region%ice%SL_a( vi) - region%ice%Hb_a( vi)) * (seawater_density / ice_density)))
        
        ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area)
      END IF     
      
    END DO
    
    buffer(1) = ice_area
    buffer(2) = ice_volume
    buffer(3) = ice_volume_above_flotation

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, buffer, 3, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
    
    region% ice_area = buffer(1)
    region% ice_volume = buffer(2)
    region% ice_volume_above_flotation = buffer(3) 

    ! Calculate GMSL contribution
    region%GMSL_contribution = -1._dp * (region%ice_volume_above_flotation - region%ice_volume_above_flotation_PD)
    
  END SUBROUTINE calculate_icesheet_volume_and_area
  SUBROUTINE calculate_PD_sealevel_contribution( region)
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    INTEGER                                       :: i, j
    REAL(dp)                                      :: ice_volume, thickness_above_flotation, ice_volume_above_flotation
    
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp
    
    DO i = region%refgeo_PD%grid%i1, region%refgeo_PD%grid%i2
    DO j = 1, region%refgeo_PD%grid%ny
    
      ! Thickness above flotation
      IF (region%refgeo_PD%Hi_grid( i,j) > 0._dp) THEN
        thickness_above_flotation = MAX(0._dp, region%refgeo_PD%Hi_grid( i,j) - MAX(0._dp, (0._dp - region%refgeo_PD%Hb_grid( i,j)) * (seawater_density / ice_density)))
      ELSE
        thickness_above_flotation = 0._dp
      END IF
      
      ! Ice volume (above flotation) in m.s.l.e
      ice_volume                 = ice_volume                 + region%refgeo_PD%Hi_grid( i,j)   * region%refgeo_PD%grid%dx * region%refgeo_PD%grid%dx * ice_density / (seawater_density * ocean_area)
      ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation        * region%refgeo_PD%grid%dx * region%refgeo_PD%grid%dx * ice_density / (seawater_density * ocean_area)
      
    END DO
    END DO
    
    CALL MPI_ALLREDUCE( ice_volume                , region%ice_volume_PD,                 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation_PD, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
    
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
!      T2m_mean = T2m_mean + SUM(region%climate%applied%ti2m(vi,:)) * region%mesh%A(vi) / (12._dp * (region%mesh%xmax - region%mesh%xmin) * (region%mesh%ymax - region%mesh%ymin))
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

END MODULE UFEMISM_main_model
