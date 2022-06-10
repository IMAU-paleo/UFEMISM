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
    WRITE(UNIT = 1337, FMT = '(A)') '% Time         : in yr, so LGM occurs at -21000'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_obs       : global mean sea level from an observational or synthetic record, in m w.r.t. PD, so a sea-level drop shows up as a negative number'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_mod       : global mean sea level contribution from all ice sheets'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_NAM       : global mean sea level contribution from the North American ice sheet (= ice volume above flotation / ocean area)'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_EAS       : global mean sea level contribution from the Eurasian       ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_GRL       : global mean sea level contribution from the Greenland      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_ANT       : global mean sea level contribution from the Antarctic      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% CO2_obs      : observed CO2          from a prescribed record (set to zero if no record is prescribed)'
    WRITE(UNIT = 1337, FMT = '(A)') '% CO2_mod      : modelled CO2          from the inverse routine (set to zero if the inverse routine is not used)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_obs     : observed benthic d18O from a prescribed record (set to zero if no record is prescribed)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_mod     : modelled benthic d18O from the inverse routine (set to zero if the inverse routine is not used)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_ice     : contribution to benthic d18O from global ice volume'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_Tdw     : contribution to benthic d18O from deep-water temperature change'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_NAM     : contribution to benthic d18O       from the North American ice sheet (= mean isotope content * sea-level equivalent ice volume / mean ocean depth)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_EAS     : contribution to benthic d18O       from the Eurasian       ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_GRL     : contribution to benthic d18O       from the Greenland      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_ANT     : contribution to benthic d18O       from the Antarctic      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% dT_glob      : global mean annual surface temperature change (scaled to sea-level)'
    WRITE(UNIT = 1337, FMT = '(A)') '% dT_dw        : deep-water temperature anomaly'
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A,A)') '       Time     SL_obs     SL_mod     SL_NAM     SL_EAS     SL_GRL     SL_ANT    CO2_obs    CO2_mod   d18O_obs', &
                                      '   d18O_mod   d18O_ice   d18O_Tdw   d18O_NAM   d18O_EAS   d18O_GRL   d18O_ANT    dT_glob      dT_dw'

    CLOSE(UNIT = 1337)

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
    REAL(dp)                                           :: sealevel_obs

    IF (.NOT. par%master) RETURN

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Deal with uninitialised variables
    IF (C%choice_sealevel_model == 'prescribed') THEN
       sealevel_obs = forcing%sealevel_obs
    ELSEIF (C%choice_sealevel_model == 'fixed') THEN
       sealevel_obs = C%fixed_sealevel
     ELSE
       sealevel_obs = 0._dp
    END IF

    ! The general output file
    ! =======================

    filename = TRIM(C%output_dir) // 'aa_general_output.txt'
    OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')

    WRITE(UNIT = 1337, FMT = '(19F11.2)') &
      time,                               &   ! 1  - time
      sealevel_obs,                       &   ! 2  - observed sea level from prescribed record (if any)
      SL_glob,                            &   ! 3  - contribution to GMSL from all ice sheets
      SL_NAM,                             &   ! 4  - contribution to GMSL from North America
      SL_EAS,                             &   ! 5  - contribution to GMSL from Eurasia
      SL_GRL,                             &   ! 6  - contribution to GMSL from Greenland
      SL_ANT,                             &   ! 7  - contribution to GMSL from Antarctica
      forcing%CO2_obs,                    &   ! 8  - observed CO2  from prescribed record (if any)
      forcing%CO2_mod,                    &   ! 9  - modelled CO2     (if any)
      forcing%d18O_obs,                   &   ! 10 - observed d18O from prescribed record (if any)
      forcing%d18O_mod,                   &   ! 11 - modelled d18O    (always)
      forcing%d18O_from_ice_volume_mod,   &   ! 12 - contribution to modelled d18O from ice volume
      forcing%d18O_from_temperature_mod,  &   ! 13 -     ""            ""          ""   deep-sea temperature change
      forcing%d18O_NAM,                   &   ! 14 - mean isotope content of North America
      forcing%d18O_EAS,                   &   ! 15 - mean isotope content of Eurasia
      forcing%d18O_GRL,                   &   ! 16 - mean isotope content of Greenland
      forcing%d18O_ANT,                   &   ! 17 - mean isotope content of Antarctica
      forcing%dT_glob,                    &   ! 18 - global mean surface temperature anomaly
      forcing%dT_deepwater                    ! 19 - deep-water temperature anomaly

    CLOSE(UNIT = 1337)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_text_output

END MODULE global_text_output_module
