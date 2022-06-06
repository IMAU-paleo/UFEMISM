MODULE global_text_output_module

  USE mpi
  USE configuration_module,        ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,             ONLY: par, sync, ierr, cerr
  USE data_types_module,           ONLY: type_forcing_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE create_text_output_files

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_text_output_files'
    CHARACTER(LEN=256)                                 :: filename

    IF (.NOT. par%master) RETURN

    ! Add routine to path
    CALL init_routine( routine_name)

    ! The general output file
    ! =======================

    filename = TRIM(C%output_dir) // 'aa_general_output.txt'
    OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')

    WRITE(UNIT = 1337, FMT = '(A)') '% UFEMISM global output data'
    WRITE(UNIT = 1337, FMT = '(A)') '%'
    WRITE(UNIT = 1337, FMT = '(A)') '% Time     : in yr, so LGM occurs at -21000'
    WRITE(UNIT = 1337, FMT = '(A)') '% sealevel : global mean sea level in m w.r.t. PD, so a sea-level drop shows up as a negative number'
    WRITE(UNIT = 1337, FMT = '(A)') '% CO2_obs  : observed CO2          from a prescribed record (set to zero if no record is prescribed)'
    WRITE(UNIT = 1337, FMT = '(A)') '% CO2_mod  : modelled CO2          from the inverse routine (set to zero if the inverse routine is not used)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_obs : observed benthic d18O from a prescribed record (set to zero if no record is prescribed)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_mod : modelled benthic d18O from the inverse routine (set to zero if the inverse routine is not used)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_ice : contribution to benthic d18O from global ice volume'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_Tdw : contribution to benthic d18O from deep-water temperature change'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_NAM   : global mean sea level contribution from the North American ice sheet (= ice volume above flotation / ocean area)'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_EAS   : global mean sea level contribution from the Eurasian       ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_GRL   : global mean sea level contribution from the Greenland      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_ANT   : global mean sea level contribution from the Antarctic      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_NAM : contribution to benthic d18O       from the North American ice sheet (= mean isotope content * sea-level equivalent ice volume / mean ocean depth)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_EAS : contribution to benthic d18O       from the Eurasian       ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_GRL : contribution to benthic d18O       from the Greenland      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_ANT : contribution to benthic d18O       from the Antarctic      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% dT_glob  : global mean annual surface temperature change (scaled to sea-level)'
    WRITE(UNIT = 1337, FMT = '(A)') '% dT_dw    : deep-water temperature anomaly'
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A,A)') '     Time     sealevel    CO2_obs    CO2_mod   d18O_obs   d18O_mod   d18O_ice   d18O_Tdw     ', &
                                      'SL_NAM     SL_EAS     SL_GRL     SL_ANT   d18O_NAM   d18O_EAS   d18O_GRL   d18O_ANT    dT_glob      dT_dw'

    CLOSE(UNIT = 1337)

    ! The memory use log
    ! ==================

    IF (C%do_write_memory_tracker) THEN

      filename = TRIM(C%output_dir) // 'aa_memory_use_log.txt'
      OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')

      WRITE(UNIT = 1337, FMT = '(A)') '% UFEMISM memory use log'
      WRITE(UNIT = 1337, FMT = '(A)') ''

      CLOSE(UNIT = 1337)

    END IF ! IF (C%do_write_memory_tracker) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_text_output_files
  SUBROUTINE write_text_output( time, SL_glob, SL_NAM, SL_EAS, SL_GRL, SL_ANT, forcing)
    ! Write data to global output file

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                        INTENT(IN)        :: time, SL_glob
    REAL(dp),                        INTENT(IN)        :: SL_NAM, SL_EAS, SL_GRL, SL_ANT
    TYPE(type_forcing_data),         INTENT(IN)        :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_text_output'
    CHARACTER(LEN=256)                                 :: filename

    IF (.NOT. par%master) RETURN

    ! Add routine to path
    CALL init_routine( routine_name)

    ! The general output file
    ! =======================

    filename = TRIM(C%output_dir) // 'aa_general_output.txt'
    OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')

    WRITE(UNIT = 1337, FMT = '(19F11.2)') &
      time,                               &   ! 1  - time
      SL_glob,                            &   ! 2  - global mean sea level
      forcing%CO2_obs,                    &   ! 3  - observed CO2  from prescribed record (if any)
      forcing%CO2_mod,                    &   ! 4  - modelled CO2     (if any)
      forcing%d18O_obs,                   &   ! 5  - observed d18O from prescribed record (if any)
      forcing%d18O_mod,                   &   ! 6  - modelled d18O    (always)
      forcing%d18O_from_ice_volume_mod,   &   ! 7  - contribution to modelled d18O from ice volume
      forcing%d18O_from_temperature_mod,  &   ! 8  -     ""            ""          ""   deep-sea temperature change
      SL_NAM,                             &   ! 9  - contribution to GMSL from North America
      SL_EAS,                             &   ! 10 - contribution to GMSL from Eurasia
      SL_GRL,                             &   ! 11 - contribution to GMSL from Greenland
      SL_ANT,                             &   ! 12 - contribution to GMSL from Antarctica
      forcing%d18O_NAM,                   &   ! 13 - mean isotope content of North America
      forcing%d18O_EAS,                   &   ! 14 - mean isotope content of Eurasia
      forcing%d18O_GRL,                   &   ! 15 - mean isotope content of Greenland
      forcing%d18O_ANT,                   &   ! 16 - mean isotope content of Antarctica
      forcing%dT_glob,                    &   ! 17 - global mean surface temperature anomaly
      forcing%dT_deepwater                    ! 18 - deep-water temperature anomaly

    CLOSE(UNIT = 1337)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_text_output

END MODULE global_text_output_module
