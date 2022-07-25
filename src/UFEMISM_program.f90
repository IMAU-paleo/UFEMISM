PROGRAM UFEMISM_program
! The Utrecht FinitE voluMe Ice Sheet Model (UFEMISM),
  ! by Tijn Berends and Jorjo Bernales, 2019-2022.
  ! Institute for Marine and Atmospheric Research Utrecht (IMAU)
  !
  ! e-mail: c.j.berends@uu.nl / j.a.bernalesconcha@uu.nl
  !
  ! Model optimisation and distributed-memory version thanks
  ! to Victor Azizi, at the Netherlands eScience Center.
  !
  ! After some initialisation work (e.g. starting the program on
  ! multiple cores using Open MPI, reading the config file, creating
  ! an output folder, etc.), this program runs up to five copies of
  ! the ice-sheet model (North America, Eurasia, Greenland, Antarctica,
  ! and Patagonia). These are all run individually for a prescribed
  ! number of years (the "coupling interval") through the subroutines
  ! in the UFEMISM_main_model module, after which control is passed
  ! back to this program. At this point, the sea-level model SELEN is
  ! optionally called, some global output data is calculated and
  ! written to output files, and the coupling loop is run again.
  !
  ! The five ice-sheet models are five instances of the "model_region"
  ! data type (declared in the data_types module), which is accepted
  ! as an argument by the "run_model" subroutine.
  !
  ! Some general notes:
  ! - Model data are arranged into several large structures, all of
  !   which are declared in the data_types module, and USEd by the
  !   different model subroutines. This prevents dependency problems
  !   during compilation, and as a result the modules look very clean
  !   and easy to read.
  ! - The general rule is that any kind of data that is a property
  !   only of, e.g., the adaptive mesh (vertex coordinates, cells,
  !   etc.) and not of the ice-sheet model in particular, will be
  !   stored in the "mesh" data type. Likewise, any subroutines that
  !   perform operations on the mesh which are not exclusive to the
  !   ice-sheet model (such as calculation of derivatives or remapping)
  !   are contained in the different "mesh_XXXX" modules.
  ! - Because of the adaptive mesh, memory for the different model data
  !   fields needs to be reallocated when the mesh is updated. Since
  !   this concerns data fields that are a property of the ice-sheet
  !   model components, rather than the mesh itself, this is done in
  !   the "remap_COMPONENT" routines contained in the different model
  !   component modules.

#include <petsc/finclude/petscksp.h>
  USE mpi
  USE petscksp
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER
  USE petsc_module,                ONLY: perr
  USE configuration_module,        ONLY: dp, routine_path, write_total_model_time_to_screen, initialise_model_configuration, &
                                         C, crash, warning, reset_resource_tracker
  USE parallel_module,             ONLY: initialise_parallelisation, par, sync, ierr
  USE data_types_module,           ONLY: type_netcdf_resource_tracker, type_model_region
  USE forcing_module,              ONLY: initialise_global_forcing

  USE zeta_module,                 ONLY: initialise_zeta_discretisation
  USE UFEMISM_main_model,          ONLY: initialise_model, run_model
  USE netcdf_module,               ONLY: create_resource_tracking_file, write_to_resource_tracking_file

  IMPLICIT NONE

  CHARACTER(LEN=256), PARAMETER          :: version_number = '0.1'

  ! The four model regions
  TYPE(type_model_region)                :: NAM, EAS, GRL, ANT

  ! Coupling timer
  REAL(dp)                               :: t_coupling, t_end_models

  ! Computation time tracking
  TYPE(type_netcdf_resource_tracker)     :: resources
  REAL(dp)                               :: tstart, tstop, t1, tcomp_loop

  ! ======================================================================================

  routine_path = 'UFEMISM_program'

  ! Initialise MPI and PETSc
  CALL initialise_parallelisation
  CALL PetscInitialize( PETSC_NULL_CHARACTER, perr)

  IF (par%master) WRITE(0,*) ''
  IF (par%master) WRITE(0,*) '=================================================='
  IF (par%master) WRITE(0,'(A,A,A,I3,A)') ' ===== Running MINIMISM v', TRIM(version_number), ' on ', par%n, ' cores ====='
  IF (par%master) WRITE(0,*) '=================================================='
  IF (par%master) WRITE(0,*) ''

  tstart = MPI_WTIME()
  t1     = MPI_WTIME()

  ! Set up the model configuration from the provided config file(s) and create an output directory
  ! ==============================================================================================

  CALL initialise_model_configuration( version_number)

  ! ===== Initialise parameters for the vertical scaled coordinate transformation =====
  ! (the same for all model regions, so stored in the "C" structure)
  ! ===================================================================================

  CALL initialise_zeta_discretisation

  ! ===== Initialise global forcing data (d18O, CO2, insolation, geothermal heat flux) =====
  ! ========================================================================================

  CALL initialise_global_forcing

  ! ===== Create the resource tracking output file =====
  ! ====================================================

  CALL create_resource_tracking_file( resources)

  ! ===== Initialise the model regions ======
  ! =========================================

  IF (C%do_NAM) CALL initialise_model( NAM, 'NAM')
  IF (C%do_EAS) CALL initialise_model( EAS, 'EAS')
  IF (C%do_GRL) CALL initialise_model( GRL, 'GRL')
  IF (C%do_ANT) CALL initialise_model( ANT, 'ANT')

! =============================
! ===== The big time loop =====
! =============================

  t_coupling = C%start_time_of_run

  DO WHILE (t_coupling < C%end_time_of_run)

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,'(A,F9.3,A)') ' Coupling model: t = ', t_coupling/1000._dp, ' kyr'

    ! Run all four model regions for 100 years
    t_end_models = MIN(C%end_time_of_run, t_coupling + C%dt_coupling)

    IF (C%do_NAM) CALL run_model( NAM, t_end_models)
    IF (C%do_EAS) CALL run_model( EAS, t_end_models)
    IF (C%do_GRL) CALL run_model( GRL, t_end_models)
    IF (C%do_ANT) CALL run_model( ANT, t_end_models)

    ! Advance coupling time
    t_coupling = t_end_models

    ! Write resource use to the resource tracking file
    tcomp_loop = MPI_WTIME() - t1
    CALL write_to_resource_tracking_file( resources, t_coupling, tcomp_loop)
    t1 = MPI_WTIME()
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
