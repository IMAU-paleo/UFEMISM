MODULE forcing_module

  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync, &
                                         allocate_shared_memory_int_0D, allocate_shared_memory_dp_0D, &
                                         allocate_shared_memory_int_1D, allocate_shared_memory_dp_1D, &
                                         allocate_shared_memory_int_2D, allocate_shared_memory_dp_2D, &
                                         allocate_shared_memory_int_3D, allocate_shared_memory_dp_3D, &
                                         deallocate_shared_memory
  USE data_types_module,           ONLY: type_forcing_data, type_mesh, type_SMB_model
  USE netcdf_module,               ONLY: inquire_insolation_data_file, read_insolation_data_file_time_lat, read_insolation_data_file

  IMPLICIT NONE
  
  ! The data structure containing model forcing data - CO2 record, d18O record, (global) insolation record
  ! Updated at every coupling time step, reading the two timeframes in the file enveloping the coupling timestep,
  ! so that during the next model loops, actual Q_TOA can be calculated by interpolating between them.
  
  TYPE(type_forcing_data), SAVE :: forcing
    
CONTAINS

  SUBROUTINE update_insolation_forcing_data(t_coupling)
    ! Read the text file containing the insolation forcing data. Only read the time frames enveloping the current
    ! coupling timestep to save on memory usage. Only done by master.
    
    ! NOTE: assumes time in forcing file is in ky
    
    IMPLICIT NONE

    REAL(dp),                            INTENT(IN)    :: t_coupling
    
    ! Not needed for benchmark experiments
    IF (C%do_eismint_experiment) RETURN 
        
    CALL read_insolation_data_file(forcing, t_coupling)
    
  END SUBROUTINE update_insolation_forcing_data

  SUBROUTINE initialise_forcing_data
    ! Allocate shared memory for the forcing data fields
    
    IMPLICIT NONE   
    
    INTEGER                                            :: i
        
    CALL allocate_shared_memory_dp_0D(                      forcing%ins_t0,      forcing%wins_t0     )
    CALL allocate_shared_memory_dp_0D(                      forcing%ins_t1,      forcing%wins_t1     )
    IF (par%master) THEN
      forcing%ins_t0 = C%start_time_of_run
      forcing%ins_t1 = C%end_time_of_run
    END IF
    CALL sync
    
    ! Not needed for benchmark experiments
    IF (C%do_eismint_experiment) RETURN 
    
    ! Inquire into the insolation forcing netcdf file
    forcing%ins_netcdf%filename = C%filename_insolation
    DO i = 0, par%n-1
      IF (i==par%i) THEN
        CALL inquire_insolation_data_file(forcing)
      END IF
      CALL sync
    END DO
    
    ! Insolation    
    CALL allocate_shared_memory_dp_1D(forcing%ins_nyears,   forcing%ins_time,    forcing%wins_time   )
    CALL allocate_shared_memory_dp_1D(forcing%ins_nlat,     forcing%ins_lat,     forcing%wins_lat    )
    CALL allocate_shared_memory_dp_2D(forcing%ins_nlat, 12, forcing%ins_Q_TOA0,  forcing%wins_Q_TOA0 )
    CALL allocate_shared_memory_dp_2D(forcing%ins_nlat, 12, forcing%ins_Q_TOA1,  forcing%wins_Q_TOA1 )
    
    forcing%ins_t0 = C%start_time_of_run
    forcing%ins_t1 = C%start_time_of_run + 50._dp
    
    If (par%master) THEN
      CALL read_insolation_data_file_time_lat(forcing) 
    END IF
    CALL sync
    
  END SUBROUTINE initialise_forcing_data

END MODULE forcing_module
