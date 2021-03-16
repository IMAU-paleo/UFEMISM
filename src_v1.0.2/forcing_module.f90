MODULE forcing_module
 
  USE mpi
  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync, &
                                         allocate_shared_int_0D, allocate_shared_dp_0D, &
                                         allocate_shared_int_1D, allocate_shared_dp_1D, &
                                         allocate_shared_int_2D, allocate_shared_dp_2D, &
                                         allocate_shared_int_3D, allocate_shared_dp_3D, &
                                         deallocate_shared
  USE data_types_module,           ONLY: type_forcing_data, type_mesh, type_SMB_model
  USE netcdf_module,               ONLY: inquire_insolation_data_file, read_insolation_data_file_time_lat, read_insolation_data_file

  IMPLICIT NONE
  
  ! The data structure containing model forcing data - CO2 record, d18O record, (global) insolation record
  ! Updated at every coupling time step. For Q_TOA, only the two timeframes in the file enveloping the coupling time
  ! are read from the NetCDF file, so that during the next model loops, actual Q_TOA can be calculated by interpolating between them.
  
  TYPE(type_forcing_data), SAVE :: forcing
    
CONTAINS

  SUBROUTINE update_insolation_forcing_data( t_coupling)
    ! Read the text file containing the insolation forcing data. Only read the time frames enveloping the current
    ! coupling timestep to save on memory usage. Only done by master.
    
    ! NOTE: assumes time in forcing file is in kyr
    
    IMPLICIT NONE

    REAL(dp),                            INTENT(IN)    :: t_coupling
    
    ! Local variables
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: ti0, ti1
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        RETURN
      ELSE 
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_forcing_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    IF (par%master) WRITE(0,*) ' Updating insolation forcing...'
    
    ! Initialise at zero
    IF (par%master) THEN
      forcing%ins_Q_TOA0 = 0._dp
      forcing%ins_Q_TOA1 = 0._dp
    END IF
    CALL sync
    
    ! Check if data for model time is available
    IF (t_coupling < forcing%ins_time(1)) THEN
      IF (par%master) WRITE(0,*) '  ERROR: insolation data only available between ', MINVAL(forcing%ins_time), ' y and ', MAXVAL(forcing%ins_time), ' y'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Find time indices to be read
    IF (par%master) THEN
      IF (t_coupling <= forcing%ins_time( forcing%ins_nyears)) THEN
        IF (par%master) WRITE(0,*) ''
        ti1 = 1
        DO WHILE (forcing%ins_time(ti1) < t_coupling)
          ti1 = ti1 + 1
        END DO
        ti0 = ti1 - 1
        
        forcing%ins_t0 = forcing%ins_time(ti0)
        forcing%ins_t1 = forcing%ins_time(ti1)
      ELSE
        IF (par%master) WRITE(0,*) '  WARNING: using constant PD insolation for future projections!'
        IF (par%master) WRITE(0,*) ''
        ti0 = forcing%ins_nyears
        ti1 = forcing%ins_nyears
        
        forcing%ins_t0 = forcing%ins_time(ti0) - 1._dp
        forcing%ins_t1 = forcing%ins_time(ti1)
      END IF
    END IF ! IF (par%master) THEN
        
    ! Read new insolation fields from the NetCDF file
    IF (par%master) CALL read_insolation_data_file( forcing, ti0, ti1)
    CALL sync
    
  END SUBROUTINE update_insolation_forcing_data

  SUBROUTINE initialise_forcing_data
    ! Allocate shared memory for the forcing data fields
    
    IMPLICIT NONE
    
    INTEGER :: cerr, ierr
        
    ! The times at which we have insolation fields from Laskar, between which we'll interpolate
    ! to find the insolation at model time (ins_t0 < model_time < ins_t1)
    
    CALL allocate_shared_dp_0D( forcing%ins_t0, forcing%wins_t0)
    CALL allocate_shared_dp_0D( forcing%ins_t1, forcing%wins_t1)
    
    IF (par%master) THEN
      forcing%ins_t0 = C%start_time_of_run
      forcing%ins_t1 = C%end_time_of_run
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        RETURN
      ELSE 
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_forcing_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Inquire into the insolation forcing netcdf file    
    CALL allocate_shared_int_0D( forcing%ins_nyears, forcing%wins_nyears)
    CALL allocate_shared_int_0D( forcing%ins_nlat,   forcing%wins_nlat  )
    
    forcing%netcdf%filename = C%filename_insolation
    
    IF (par%master) CALL inquire_insolation_data_file( forcing)
    CALL sync
    
    ! Insolation    
    CALL allocate_shared_dp_1D( forcing%ins_nyears,   forcing%ins_time,    forcing%wins_time   )
    CALL allocate_shared_dp_1D( forcing%ins_nlat,     forcing%ins_lat,     forcing%wins_lat    )
    CALL allocate_shared_dp_2D( forcing%ins_nlat, 12, forcing%ins_Q_TOA0,  forcing%wins_Q_TOA0 )
    CALL allocate_shared_dp_2D( forcing%ins_nlat, 12, forcing%ins_Q_TOA1,  forcing%wins_Q_TOA1 )
    
    forcing%ins_t0 = C%start_time_of_run
    forcing%ins_t1 = C%start_time_of_run + 50._dp
    
    IF (par%master) CALL read_insolation_data_file_time_lat( forcing)
    CALL sync
    
  END SUBROUTINE initialise_forcing_data

END MODULE forcing_module