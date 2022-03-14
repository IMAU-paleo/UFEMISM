PROGRAM UFEMISM_program
! The Utrecht FinitE voluMe Ice Sheet Model (UFEMISM), by Tijn Berends, 2019.
! Institute for Marine and Atmospheric Research Utrecht (IMAU)
!
! e-mail: c.j.berends@uu.nl
!
! After some initialisation work (starting the program on multiple cores using MPI_INIT,
! reading the config file, creating an output folder, etc.), this program runs the four
! copies of the ice-sheet model (North America, Eurasia, Greenland and Antarctica).
! These are all run individually for 100 years (the "coupling interval"), after which
! control is passed back to this program. At this point, the sea-level model SELEN will be
! called (not implemented yet), some global output data is calculated and written to the
! output file, and the coupling loop is run again.
!
! The four ice-sheet models are four instances of the "model_region" data type (declared in
! the data_types module), which is accepted as an argument by the "run_model" subroutine.
!
! Some general notes:
! - Model data is arranged into several large structures, all of which are declared in the 
!   data_types module, and USEd by the different model subroutines. This prevents dependency
!   problems during compiling, and makes the modules very clean and easy to read.
! - The general rule is that any kind of data that is a property only of the adaptive mesh,
!   not of the ice-sheet model in particular (vertex coordinates, neighbour functions, etc.),
!   is stored in the "mesh" data type. Likewise, any subroutines that perform operations on
!   the mesh which are not exclusive to the ice-sheet model (such as calculation of derivatives
!   or remapping) are contained in the different "mesh_XXXX" modules.
! - Because of the adaptive mesh, memory for the different model data fields needs to be
!   reallocated when the mesh is updated. Since this concerns data fields that are a property
!   of the ice-sheet model components, rather than the mesh itself, this is done in the
!   "remap_COMPONENT" routines contained in the different model component modules.
  
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE petscksp
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER
  USE petsc_module,                ONLY: perr
  USE configuration_module,        ONLY: dp, C, routine_path, crash, warning, initialise_model_configuration, write_total_model_time_to_screen, &
                                         reset_resource_tracker
  USE parallel_module,             ONLY: par, sync, ierr, cerr, initialise_parallelisation
  USE data_types_module,           ONLY: type_model_region, type_climate_matrix, type_netcdf_resource_tracker
  USE forcing_module,              ONLY: forcing, initialise_insolation_data, update_insolation_data, initialise_CO2_record, update_CO2_at_model_time, &
                                         initialise_d18O_record, update_d18O_at_model_time, initialise_d18O_data, update_global_mean_temperature_change_history, &
                                         calculate_modelled_d18O, initialise_inverse_routine_data, inverse_routine_global_temperature_offset, inverse_routine_CO2, &
                                         initialise_geothermal_heat_flux
  USE climate_module,              ONLY: initialise_climate_matrix
  USE zeta_module,                 ONLY: initialise_zeta_discretisation
  USE global_text_output_module,   ONLY: create_text_output_files, write_text_output
  USE UFEMISM_main_model,          ONLY: initialise_model, run_model
  USE netcdf_module,               ONLY: create_resource_tracking_file, write_to_resource_tracking_file

  IMPLICIT NONE
  
  CHARACTER(LEN=256), PARAMETER          :: version_number = '1.2'
  
  ! The four model regions
  TYPE(type_model_region)                :: NAM, EAS, GRL, ANT
  
  ! The global climate matrix
  TYPE(type_climate_matrix)              :: matrix
  
  ! Coupling timer
  REAL(dp)                               :: t_coupling, t_end_models
  REAL(dp)                               :: GMSL_NAM, GMSL_EAS, GMSL_GRL, GMSL_ANT, GMSL_glob
  
  ! Computation time tracking
  TYPE(type_netcdf_resource_tracker)     :: resources
  REAL(dp)                               :: tstart, tstop
  
  ! ======================================================================================
  
  routine_path = 'UFEMISM_program'
  
  ! Initialise MPI and PETSc
  CALL initialise_parallelisation
  CALL PetscInitialize( PETSC_NULL_CHARACTER, perr)
  
  IF (par%master) WRITE(0,*) ''
  IF (par%master) WRITE(0,*) '=================================================='
  IF (par%master) WRITE(0,'(A,A,A,I3,A)') ' ===== Running UFEMISM v', TRIM(version_number), ' on ', par%n, ' cores ====='
  IF (par%master) WRITE(0,*) '=================================================='
  IF (par%master) WRITE(0,*) ''
  
  tstart = MPI_WTIME()
    
  ! Set up the model configuration from the provided config file(s) and create an output directory
  ! ==============================================================================================
  
  CALL initialise_model_configuration( version_number)
    
  ! ===== Initialise parameters for the vertical scaled coordinate transformation =====
  ! (the same for all model regions, so stored in the "C" structure)
  ! ===================================================================================
  
  CALL initialise_zeta_discretisation
  
  ! ===== Initialise forcing data =====
  ! ===================================
  
  CALL initialise_d18O_data
  CALL initialise_insolation_data
  CALL initialise_CO2_record
  CALL initialise_d18O_record
  CALL initialise_inverse_routine_data
  CALL initialise_geothermal_heat_flux
  
  ! ===== Create the resource tracking output file =====
  ! ====================================================
  
  CALL create_resource_tracking_file( resources)
  
  ! ===== Initialise the climate matrix =====
  ! =========================================
  
  CALL initialise_climate_matrix(matrix)  

  ! ===== Initialise the model regions ======
  ! =========================================
  
  IF (C%do_NAM) CALL initialise_model( NAM, 'NAM', matrix)
  IF (C%do_EAS) CALL initialise_model( EAS, 'EAS', matrix)
  IF (C%do_GRL) CALL initialise_model( GRL, 'GRL', matrix)
  IF (C%do_ANT) CALL initialise_model( ANT, 'ANT', matrix)
    
  ! Determine GMSL contributions of all simulated ice sheets
  GMSL_NAM = 0._dp
  GMSL_EAS = 0._dp
  GMSL_GRL = 0._dp
  GMSL_ANT = 0._dp
  IF (C%do_NAM) GMSL_NAM = NAM%GMSL_contribution
  IF (C%do_EAS) GMSL_EAS = EAS%GMSL_contribution
  IF (C%do_GRL) GMSL_GRL = GRL%GMSL_contribution
  IF (C%do_ANT) GMSL_ANT = ANT%GMSL_contribution
  GMSL_glob = GMSL_NAM + GMSL_EAS + GMSL_GRL + GMSL_ANT 
  
  ! Determine d18O contributions of all simulated ice sheets
  !CALL update_global_mean_temperature_change_history( NAM, EAS, GRL, ANT)
  CALL calculate_modelled_d18O( NAM, EAS, GRL, ANT)
    
  ! Write global data at t=0 to output file
  CALL write_text_output( &
    C%start_time_of_run,               &  ! time
    GMSL_glob,                         &  ! global mean sea level
    forcing%CO2_obs,                   &  ! observed CO2  from prescribed record (if any)
    forcing%CO2_mod,                   &  ! modelled CO2     (if any)
    forcing%d18O_obs,                  &  ! observed d18O from prescribed record (if any)
    forcing%d18O_mod,                  &  ! modelled d18O    (always)
    forcing%d18O_from_ice_volume_mod,  &  ! contribution to modelled d18O from ice volume
    forcing%d18O_from_temperature_mod, &  !     ""            ""          ""   deep-sea temperature change
    GMSL_NAM,                          &  ! contribution to GMSL from North America
    GMSL_EAS,                          &  ! contribution to GMSL from Eurasia
    GMSL_GRL,                          &  ! contribution to GMSL from Greenland
    GMSL_ANT,                          &  ! contribution to GMSL from Antarctica
    forcing%d18O_NAM,                  &  ! mean isotope content of North America
    forcing%d18O_EAS,                  &  ! mean isotope content of Eurasia
    forcing%d18O_GRL,                  &  ! mean isotope content of Greenland
    forcing%d18O_ANT,                  &  ! mean isotope content of Antarctica
    forcing%dT_glob,                   &  ! global mean surface temperature anomaly
    forcing%dT_deepwater               )  ! deep-water temperature anomaly
  
  ! ===== The big time loop =====
  ! =============================
  
  t_coupling = C%start_time_of_run
  
  DO WHILE (t_coupling < C%end_time_of_run)
  
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,'(A,F9.3,A)') ' Coupling model: t = ', t_coupling/1000._dp, ' kyr'
  
    ! Update global insolation forcing, CO2, and d18O at the current model time
    CALL update_insolation_data(    t_coupling)
    CALL update_CO2_at_model_time(  t_coupling)
    CALL update_d18O_at_model_time( t_coupling)
    
    ! Update regional sea level (needs to be moved to separate subroutine at some point!)
    IF (C%choice_sealevel_model == 'fixed') THEN
      IF (C%do_NAM) NAM%ice%SL_a( NAM%mesh%vi1:NAM%mesh%vi2) = C%fixed_sealevel
      IF (C%do_EAS) EAS%ice%SL_a( EAS%mesh%vi1:EAS%mesh%vi2) = C%fixed_sealevel
      IF (C%do_GRL) GRL%ice%SL_a( GRL%mesh%vi1:GRL%mesh%vi2) = C%fixed_sealevel
      IF (C%do_ANT) ANT%ice%SL_a( ANT%mesh%vi1:ANT%mesh%vi2) = C%fixed_sealevel
    ELSEIF (C%choice_sealevel_model == 'eustatic') THEN
      IF (C%do_NAM) NAM%ice%SL_a( NAM%mesh%vi1:NAM%mesh%vi2) = GMSL_glob
      IF (C%do_EAS) EAS%ice%SL_a( EAS%mesh%vi1:EAS%mesh%vi2) = GMSL_glob
      IF (C%do_GRL) GRL%ice%SL_a( GRL%mesh%vi1:GRL%mesh%vi2) = GMSL_glob
      IF (C%do_ANT) ANT%ice%SL_a( ANT%mesh%vi1:ANT%mesh%vi2) = GMSL_glob
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_sealevel_model "', TRIM(C%choice_sealevel_model), '" not implemented in IMAU_ICE_program!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Run all four model regions for 100 years
    t_end_models = MIN(C%end_time_of_run, t_coupling + C%dt_coupling)
    
    IF (C%do_NAM) CALL run_model( NAM, matrix, t_end_models)
    IF (C%do_EAS) CALL run_model( EAS, matrix, t_end_models)
    IF (C%do_GRL) CALL run_model( GRL, matrix, t_end_models)
    IF (C%do_ANT) CALL run_model( ANT, matrix, t_end_models)
    
    ! Advance coupling time
    t_coupling = t_end_models
    
    ! Determine GMSL contributions
    GMSL_NAM = 0._dp
    GMSL_EAS = 0._dp
    GMSL_GRL = 0._dp
    GMSL_ANT = 0._dp
    IF (C%do_NAM) GMSL_NAM = NAM%GMSL_contribution
    IF (C%do_EAS) GMSL_EAS = EAS%GMSL_contribution
    IF (C%do_GRL) GMSL_GRL = GRL%GMSL_contribution
    IF (C%do_ANT) GMSL_ANT = ANT%GMSL_contribution
    GMSL_glob = GMSL_NAM + GMSL_EAS + GMSL_GRL + GMSL_ANT
  
    ! Calculate contributions to benthic d18O from the different ice sheets
    CALL update_global_mean_temperature_change_history( NAM, EAS, GRL, ANT)
    CALL calculate_modelled_d18O( NAM, EAS, GRL, ANT)
    
    ! If applicable, call the inverse routine to update the climate forcing parameter
    IF     (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CALL inverse_routine_global_temperature_offset
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL inverse_routine_CO2
    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      ! No inverse routine is used in these forcing methods
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in IMAU_ICE_program!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Write global data to output file
    CALL write_text_output( &
      t_coupling,  &  ! time
      GMSL_glob,   &  ! global mean sea level
      forcing%CO2_obs,                   &  ! observed CO2  from prescribed record (if any)
      forcing%CO2_mod,                   &  ! modelled CO2   (if any)
      forcing%d18O_obs,                  &  ! observed d18O from prescribed record (if any)
      forcing%d18O_mod,                  &  ! modelled d18O  (always)
      forcing%d18O_from_ice_volume_mod,  &  ! contribution to modelled d18O from ice volume
      forcing%d18O_from_temperature_mod, &  !     ""            ""          ""   deep-sea temperature change
      GMSL_NAM,    &  ! contribution to GMSL from North America
      GMSL_EAS,    &  ! contribution to GMSL from Eurasia
      GMSL_GRL,    &  ! contribution to GMSL from Greenland
      GMSL_ANT,    &  ! contribution to GMSL from Antarctica
      forcing%d18O_NAM,                  &  ! mean isotope content of North America
      forcing%d18O_EAS,                  &  ! mean isotope content of Eurasia
      forcing%d18O_GRL,                  &  ! mean isotope content of Greenland
      forcing%d18O_ANT,                  &  ! mean isotope content of Antarctica
      forcing%dT_glob,                   &  ! global mean surface temperature anomaly
      forcing%dT_deepwater               )  ! deep-water temperature anomaly
  
    ! Write resource use to the resource tracking file
    CALL write_to_resource_tracking_file( resources, t_coupling)
    CALL reset_resource_tracker
      
  END DO ! DO WHILE (t_coupling < C%end_time_of_run)
  
! ====================================
! ===== End of the big time loop =====
! ==================================== 
  
  ! Write total elapsed time to screen
  tstop = MPI_WTIME()
  IF (par%master) CALL write_total_model_time_to_screen( tstart, tstop)
  CALL sync
  
  ! Finalise MPI and PETSc
  CALL PetscFinalize( perr)
  CALL MPI_FINALIZE( ierr)
    
END PROGRAM UFEMISM_program
