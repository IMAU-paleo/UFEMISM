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

! ===== USE modules =====
! =======================

  USE mpi
  USE petscksp
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER
  USE petsc_module,                ONLY: perr
  USE configuration_module,        ONLY: dp, C, routine_path, crash, warning, initialise_model_configuration, write_total_model_time_to_screen, &
                                         reset_resource_tracker
  USE parallel_module,             ONLY: par, sync, ierr, cerr, initialise_parallelisation
  USE data_types_module,           ONLY: type_model_region, type_netcdf_resource_tracker, &
                                         type_climate_matrix_global, type_ocean_matrix_global, &
                                         type_SELEN_global
  USE forcing_module,              ONLY: forcing, initialise_global_forcing, update_global_forcing, &
                                         update_global_mean_temperature_change_history, calculate_modelled_d18O
  USE climate_module,              ONLY: initialise_climate_model_global
  USE ocean_module,                ONLY: initialise_ocean_model_global, initialise_ocean_vertical_grid
  ! USE SELEN_main_module,           ONLY: initialise_SELEN, run_SELEN
  USE zeta_module,                 ONLY: initialise_zeta_discretisation
  USE text_output_module,          ONLY: create_global_text_output, write_global_text_output
  USE UFEMISM_main_model,          ONLY: initialise_model, run_model
  USE netcdf_module,               ONLY: create_resource_tracking_file, write_to_resource_tracking_file
  USE general_sea_level_module,    ONLY: determine_GMSL_contributions, update_regional_sea_level

! ===== Main variables =====
! ==========================

  IMPLICIT NONE

  CHARACTER(LEN=256), PARAMETER          :: version_number = '1.2'

  ! The four model regions
  TYPE(type_model_region)                :: NAM, EAS, GRL, ANT

  ! The global climate and ocean matrices
  TYPE(type_climate_matrix_global)       :: climate_matrix_global
  TYPE(type_ocean_matrix_global)         :: ocean_matrix_global

  ! SELEN
  TYPE(type_SELEN_global)                :: SELEN
  REAL(dp)                               :: ocean_area
  REAL(dp)                               :: ocean_depth

  ! Coupling
  REAL(dp)                               :: t_coupling, t_end_models
  REAL(dp)                               :: GMSL_NAM, GMSL_EAS, GMSL_GRL, GMSL_ANT, GMSL_glob

  ! Computation time tracking
  TYPE(type_netcdf_resource_tracker)     :: resources
  REAL(dp)                               :: tstart, tstop, t1, tcomp_loop

! ===== START =====
! =================

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
  t1     = MPI_WTIME()

  ! == Model set-up
  ! ===============

  CALL initialise_model_configuration( version_number)

  ! == Vertical scaled coordinate transformation
  ! ============================================

  CALL initialise_zeta_discretisation

  ! == Initialise global forcing data
  ! =================================

  CALL initialise_global_forcing

  ! == Create the resource tracking output file
  ! ===========================================

  CALL create_resource_tracking_file( resources)

  ! == Initialise the climate matrix
  ! ================================

  CALL initialise_climate_model_global( climate_matrix_global)

  ! == Initialise the ocean matrix
  ! ==============================

  CALL initialise_ocean_vertical_grid
  CALL initialise_ocean_model_global( ocean_matrix_global)

  ! == Initialise the model regions
  ! ===============================

  IF (C%do_NAM) CALL initialise_model( NAM, 'NAM', climate_matrix_global, ocean_matrix_global)
  IF (C%do_EAS) CALL initialise_model( EAS, 'EAS', climate_matrix_global, ocean_matrix_global)
  IF (C%do_GRL) CALL initialise_model( GRL, 'GRL', climate_matrix_global, ocean_matrix_global)
  IF (C%do_ANT) CALL initialise_model( ANT, 'ANT', climate_matrix_global, ocean_matrix_global)

  ! == Initial GMSL contribution
  ! ============================

  CALL determine_GMSL_contributions( GMSL_NAM, GMSL_EAS, GMSL_GRL, GMSL_ANT, GMSL_glob, NAM, EAS, GRL, ANT)

  ! == Initial d18O contribution
  ! ============================

  IF ( (.NOT. C%choice_forcing_method == 'none') .AND. &
              C%do_calculate_benthic_d18O ) THEN
      CALL update_global_mean_temperature_change_history( NAM, EAS, GRL, ANT)
      CALL calculate_modelled_d18O( NAM, EAS, GRL, ANT)
  END IF

  ! == Initialise SELEN
  ! ===================

  ! IF (C%choice_GIA_model == 'SELEN' .OR. C%choice_sealevel_model == 'SELEN') THEN
  !   CALL initialise_SELEN( SELEN, NAM, EAS, GRL, ANT, version_number)
  ! END IF

  ! ! Timers and run-at-start switch
  ! IF (C%SELEN_run_at_t_start) THEN
  !   SELEN%t0_SLE = C%start_time_of_run - C%dt_SELEN
  !   SELEN%t1_SLE = C%start_time_of_run
  ! ELSE
  !   SELEN%t0_SLE = C%start_time_of_run
  !   SELEN%t1_SLE = C%start_time_of_run + C%dt_SELEN
  ! END IF

  ! == Initial global output
  ! ========================

  ! Create file and header
  CALL create_global_text_output

  ! Write global data at t=0 to output file
  CALL write_global_text_output( C%start_time_of_run, GMSL_glob, GMSL_NAM, GMSL_EAS, GMSL_GRL, GMSL_ANT, forcing)

! ===== The big time loop =====
! =============================

  t_coupling = C%start_time_of_run

  DO WHILE (t_coupling < C%end_time_of_run)

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,'(A,F9.3,A)') ' Coupling model: t = ', t_coupling/1000._dp, ' kyr'

    ! == Global forcing update (pre regional runs)
    ! ============================================

    ! Update global insolation forcing, CO2, and d18O at the current model time
    CALL update_global_forcing( NAM, EAS, GRL, ANT, t_coupling, switch = 'pre')

    ! == SELEN
    ! ========

    ! ! Solve the SLE
    ! IF (t_coupling >= SELEN%t1_SLE .AND. (C%choice_GIA_model == 'SELEN' .OR. C%choice_sealevel_model == 'SELEN')) THEN
    !   CALL run_SELEN( SELEN, NAM, EAS, GRL, ANT, t_coupling, ocean_area, ocean_depth)
    !   SELEN%t0_SLE = t_coupling
    !   SELEN%t1_SLE = t_coupling + C%dt_SELEN
    ! END IF

    ! == Regional sea level update
    ! ============================

    CALL update_regional_sea_level( NAM, EAS, GRL, ANT, GMSL_glob, t_coupling)

    ! == Regional model runs
    ! ======================

    ! Run all four model regions for C%dt_coupling years
    t_end_models = MIN(C%end_time_of_run, t_coupling + C%dt_coupling)

    IF (C%do_NAM) CALL run_model( NAM, climate_matrix_global, t_end_models)
    IF (C%do_EAS) CALL run_model( EAS, climate_matrix_global, t_end_models)
    IF (C%do_GRL) CALL run_model( GRL, climate_matrix_global, t_end_models)
    IF (C%do_ANT) CALL run_model( ANT, climate_matrix_global, t_end_models)

    ! Advance coupling time
    t_coupling = t_end_models

    ! == Global sea level update
    ! ==========================

    ! Determine GMSL contributions from all simulated ice sheets
    CALL determine_GMSL_contributions( GMSL_NAM, GMSL_EAS, GMSL_GRL, GMSL_ANT, GMSL_glob, NAM, EAS, GRL, ANT)

    ! == Global forcing update (post regional runs)
    ! =============================================

    ! Calculate contributions to benthic d18O from the different ice sheets and,
    ! if applicable, call the inverse routine to update the climate forcing parameter
    CALL update_global_forcing( NAM, EAS, GRL, ANT, t_coupling, switch = 'post')

    ! == Global output
    ! ================

    ! Write global data to output file
    CALL write_global_text_output( t_coupling, GMSL_glob, GMSL_NAM, GMSL_EAS, GMSL_GRL, GMSL_ANT, forcing)

    ! == Resource tracking output
    ! ===========================

    tcomp_loop = MPI_WTIME() - t1
    CALL write_to_resource_tracking_file( resources, t_coupling, tcomp_loop)
    t1 = MPI_WTIME()
    CALL reset_resource_tracker

  END DO

! ===== END =====
! ===============

  ! == Total elapsed time
  ! =====================

  tstop = MPI_WTIME()
  IF (par%master) CALL write_total_model_time_to_screen( tstart, tstop)
  CALL sync

  !== Finalise MPI and PETSc
  !=========================

  CALL PetscFinalize( perr)
  CALL MPI_FINALIZE( ierr)

END PROGRAM UFEMISM_program