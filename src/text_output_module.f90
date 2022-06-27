MODULE text_output_module

  USE mpi
  USE configuration_module,        ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,             ONLY: par, sync, ierr, cerr
  USE data_types_module,           ONLY: type_forcing_data, type_model_region

  IMPLICIT NONE

CONTAINS

! ===== Global output =====
! =========================

  ! Create the global text output file
  SUBROUTINE create_global_text_output

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

  END SUBROUTINE create_global_text_output

  ! Write scalar data to the global text output file
  SUBROUTINE write_global_text_output( time, SL_glob, SL_NAM, SL_EAS, SL_GRL, SL_ANT, forcing)
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

  END SUBROUTINE write_global_text_output

! ===== Regional output =====
! ===========================

  ! Create regional text output files (regional and POIs)
  SUBROUTINE create_regional_text_output( region)
    ! Creates the following text output files:
    !   time_log_REG.txt             - a log of how much computation time the different model parts take
    !   general_output_REG.txt       - some general info - ice sheet volume, average surface temperature, total mass balance, etc.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),    INTENT(IN)        :: region

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_regional_text_output'
    CHARACTER(LEN=256)                            :: filename
    CHARACTER(LEN=3)                              :: ns
    CHARACTER(LEN=1)                              :: ns3
    CHARACTER(LEN=2)                              :: ns2
    INTEGER                                       :: n, k

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! General output
    ! ==============

    filename = TRIM(C%output_dir) // 'ab_general_output_' // region%name // '.txt'
    OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')

    WRITE(UNIT = 1337, FMT = '(A)') 'General output for region ' // TRIM(region%long_name)
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A)') ' Columns in order:'
    WRITE(UNIT = 1337, FMT = '(A)') '   1)  Model time                  (years) '
    WRITE(UNIT = 1337, FMT = '(A)') '   2)  Ice volume                  (meter sea level equivalent)'
    WRITE(UNIT = 1337, FMT = '(A)') '   3)  Ice volume above flotation  (meter sea level equivalent)'
    WRITE(UNIT = 1337, FMT = '(A)') '   4)  Ice area                    (km^2)'
    WRITE(UNIT = 1337, FMT = '(A)') '   5)  Mean surface temperature    (degrees Celsius)'
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
    WRITE(UNIT = 1337, FMT = '(A)') '      Time     Volume  Volume_AF         Area     T2m       Snow       Rain       Melt   Refreeze     Runoff        SMB        BMB         MB       x_GL'

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

      filename = TRIM(C%output_dir) // 'ac_POI_' // TRIM(region%name) // '_' // TRIM(ns) // '_data.txt'
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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_regional_text_output

  ! Write scalar data to regional text output files (regional and POIs)
  SUBROUTINE write_regional_text_output( region)
    ! Write data to the following text output files:
    !   time_log_REG.txt             - a log of how much computation time the different model parts take
    !   general_output_REG.txt       - some general info - ice sheet volume, average surface temperature, total mass balance, etc.

    USE parameters_module, ONLY: ocean_area, seawater_density, ice_density

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),    INTENT(IN)        :: region

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_regional_text_output'
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

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! General output
    ! ==============

    T2m_mean                   = -273.5_dp
    total_snowfall             = 0._dp
    total_rainfall             = 0._dp
    total_melt                 = 0._dp
    total_refreezing           = 0._dp
    total_runoff               = 0._dp
    total_SMB                  = 0._dp
    total_BMB                  = 0._dp
    total_MB                   = 0._dp

    DO vi = 1, region%mesh%nV
      IF (region%ice%mask_ice_a( vi)==0) THEN

        total_BMB = total_BMB + (region%BMB%BMB(vi) * region%mesh%A(vi) / 1E9_dp) * ice_density / 1000._dp ! m3ie -> m3we -> Gt

        DO m = 1, 12
          total_snowfall   = total_snowfall   + (region%SMB%Snowfall(  vi,m) * region%mesh%A(vi) / 1E9_dp) ! Already in water equivalent
          total_rainfall   = total_rainfall   + (region%SMB%Rainfall(  vi,m) * region%mesh%A(vi) / 1E9_dp) ! Already in water equivalent
          total_melt       = total_melt       + (region%SMB%Melt(      vi,m) * region%mesh%A(vi) / 1E9_dp) ! Already in water equivalent
          total_refreezing = total_refreezing + (region%SMB%Refreezing(vi,m) * region%mesh%A(vi) / 1E9_dp) ! Already in water equivalent
          total_runoff     = total_runoff     + (region%SMB%Runoff(    vi,m) * region%mesh%A(vi) / 1E9_dp) ! Already in water equivalent
          total_SMB        = total_SMB        + (region%SMB%SMB(       vi,m) * region%mesh%A(vi) / 1E9_dp) * ice_density / 1000._dp
        END DO

      END IF

      T2m_mean = T2m_mean + SUM(region%climate_matrix%applied%T2m(vi,:)) * region%mesh%A(vi) &
                            / (12._dp * (region%mesh%xmax - region%mesh%xmin) * (region%mesh%ymax - region%mesh%ymin))

    END DO

    total_MB = total_SMB + total_BMB

    ! Average x-position of grounding line
    x_GL = 0._dp
    ! n_GL = 0
    ! DO aci = 1, region%mesh%nAc
    !   IF (region%ice%mask_gl_Ac( aci) == 1) THEN
    !     n_GL = n_GL + 1

    !     vi = region%mesh%Aci( aci,1)
    !     vj = region%mesh%Aci( aci,2)

    !     ! Find interpolated GL position
    !     TAFi = region%ice%Hi( vi) - ((region%ice%SL( vi) - region%ice%Hb( vi)) * (seawater_density / ice_density))
    !     TAFj = region%ice%Hi( vj) - ((region%ice%SL( vj) - region%ice%Hb( vj)) * (seawater_density / ice_density))
    !     lambda = TAFi / (TAFi - TAFj)

    !     x_GL = x_GL + (NORM2(region%mesh%V( vi,:)) * (1._dp - lambda)) + (NORM2(region%mesh%V( vj,:)) * lambda)
    !   END IF
    ! END DO
    ! x_GL = x_GL / n_GL

    filename = TRIM(C%output_dir) // 'ab_general_output_' // region%name // '.txt'
    OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')

    WRITE(UNIT = 1337, FMT = '(F10.1,2F11.2,F13.2,F8.2,9F11.2)') region%time, &
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

      Hi_POI = (region%ice%Hi_a(vi1  ) * w1) + (region%ice%Hi_a(vi2  ) * w2) + (region%ice%Hi_a(vi3  ) * w3)
      Hb_POI = (region%ice%Hb_a(vi1  ) * w1) + (region%ice%Hb_a(vi2  ) * w2) + (region%ice%Hb_a(vi3  ) * w3)
      Hs_POI = (region%ice%Hs_a(vi1  ) * w1) + (region%ice%Hs_a(vi2  ) * w2) + (region%ice%Hs_a(vi3  ) * w3)
      Ti_POI = (region%ice%Ti_a(vi1,:) * w1) + (region%ice%Ti_a(vi2,:) * w2) + (region%ice%Ti_a(vi3,:) * w3)

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

      filename = TRIM(C%output_dir) // 'ac_POI_' // TRIM(region%name) // '_' // TRIM(ns) // '_data.txt'
      OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')

      WRITE(UNIT = 1337, FMT = '(F10.1,3F11.2)', ADVANCE='NO') region%time, Hi_POI, Hb_POI, Hs_POI
      DO k = 1, C%nZ
        WRITE(UNIT = 1337, FMT = '(F11.2)', ADVANCE='NO') Ti_POI(k)
      END DO
      WRITE(UNIT = 1337, FMT = '(A)') ''

      CLOSE(UNIT = 1337)

    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_regional_text_output

END MODULE text_output_module
