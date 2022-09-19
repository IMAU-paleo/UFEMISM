module forcing_module
  ! Contains all the routines for reading and calculating the model forcing (CO2, d18O, insolation, sea level),
  ! as well as the "forcing" structure which stores all the results from these routines (so that they
  ! can be accessed from all four ice-sheet models and the coupling routine).

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module, only : dp, C, routine_path, init_routine, &
                                   finalise_routine, crash, warning
  use parallel_module,      only : par, sync, ierr, cerr
  use data_types_module,    only : type_forcing_data, type_mesh
  use netcdf_module,        only : inquire_insolation_file, read_insolation_file_time_lat, &
                                   read_insolation_file_timeframes, inquire_geothermal_heat_flux_file, &
                                   read_geothermal_heat_flux_file

  implicit none

  ! The data structure containing model forcing data - CO2 record, d18O record, (global) insolation record
  ! Updated at every coupling time step. For Q_TOA, only the two timeframes in the file enveloping the coupling time
  ! are read from the NetCDF file, so that during the next model loops, actual Q_TOA can be calculated by interpolating between them.

  type(type_forcing_data), save :: forcing

contains

! ===== Main routines =====
! =========================

  subroutine initialise_global_forcing
    ! Initialise global forcing data (d18O, CO2,
    ! insolation, geothermal heat flux, sea level)

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_global_forcing'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_forcing_method)

      case ('none')
        ! No need to do anything

      case ('CO2_direct')
        ! Global climate calculated based on a prescribed CO2 record
        call initialise_CO2_record

      case default
        ! Unknown option
        call crash('unknown choice_forcing_method"' // &
                    TRIM(C%choice_forcing_method) // '"!')

    end select

    ! Insolation
    call initialise_insolation_data

    ! Geothermal heat flux
    call initialise_geothermal_heat_flux_global

    ! Sea level
    if (C%choice_sealevel_model == 'prescribed') then
      ! Global sea level calculated based on a prescribed sea-level record
      call initialise_sealevel_record
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_global_forcing

! ===== Prescribed CO2 record =====
! =================================

  subroutine update_CO2_at_model_time( time)
    ! Interpolate the data in forcing%CO2 to find the
    ! value at the queried time. If time lies outside
    ! the range of forcing%CO2_time, return the first/
    ! last value. NOTE: assumes time is listed in kyr
    ! (so LGM would be -21.0).

    implicit none

    ! In/output variables:
    real(dp),                      intent(in) :: time

    ! Local variables:
    character(len=256), parameter             :: routine_name = 'update_CO2_at_model_time'
    INTEGER                                   :: il, iu
    REAL(dp)                                  :: wl, wu

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if model time is within CO2 record

    if (time / 1000._dp < minval( forcing%CO2_time)) then
      ! Model time before start of CO2 record; using constant extrapolation
      forcing%CO2_obs = forcing%CO2_record( 1)

    elseif (time / 1000._dp > maxval( forcing%CO2_time)) then
      ! Model time beyond end of CO2 record; using constant extrapolation
      forcing%CO2_obs = forcing%CO2_record( C%CO2_record_length)

    else
      ! Model time within record

      ! Find the two record times enclosing model time
      iu = 1
      do while (forcing%CO2_time(iu) * 1000._dp < time)
        iu = iu+1
      end do
      il = iu-1

      ! Compute interpolation weights
      wl = (forcing%CO2_time(iu) * 1000._dp - time) / &
           ((forcing%CO2_time(iu)-forcing%CO2_time(il)) * 1000._dp)
      wu = 1._dp - wl

      ! Interpolate CO2 record to model time
      forcing%CO2_obs = forcing%CO2_record(il) * wl + &
                        forcing%CO2_record(iu) * wu
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_CO2_at_model_time

  subroutine initialise_CO2_record
    ! Read the CO2 record specified in C%filename_CO2_record.
    ! Assumes this is an ASCII text file with at least two columns
    ! (time in kyr and CO2 in ppmv) and the number of rows being
    ! equal to C%CO2_record_length.
    ! NOTE: assumes time is listed in kyr (so LGM would be -21.0).

    implicit none

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_CO2_record'
    integer                       :: i, ios, fp

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory to take the data
    allocate(forcing%CO2_time  ( C%CO2_record_length ))
    allocate(forcing%CO2_record( C%CO2_record_length ))

    if (par%master) then
      write(*,"(A)") ''
      write(*,"(3A)") ' Initialising CO2 record from ', &
                        TRIM(C%filename_CO2_record), '...'
    end if
    call sync

    ! === Read data ===
    ! =================

    ! Read CO2 record (time and values) from specified text file
    open(  newunit=fp, file=C%filename_CO2_record, action='READ')

    do i = 1, C%CO2_record_length
      read( unit=fp, fmt=*, iostat=ios) forcing%CO2_time(i), &
                                        forcing%CO2_record(i)

      if (ios /= 0) then
        call crash('length of text file "' // &
                    TRIM(C%filename_CO2_record) // &
                    '" does not match C%CO2_record_length!')
      end if

    end do

    close( unit=fp)

    ! === Update CO2 ===
    ! ==================

    ! Set the value for the current (starting) model time
    call update_CO2_at_model_time( C%start_time_of_run)

    ! === Finalisation ===
    ! ====================

    ! Safety
    if (par%master) then
      if (C%start_time_of_run/1000._dp < forcing%CO2_time(1)) then
        call warning('Model time starts before start of CO2 record;' // &
                     ' constant extrapolation will be used in that case!')
      end if

      if (C%end_time_of_run/1000._dp > forcing%CO2_time(C%CO2_record_length)) then
        call warning('Model time will reach beyond end of CO2 record;' // &
                     ' constant extrapolation will be used in that case!')
      end if
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_CO2_record

! ===== Insolation =====
! ======================

  subroutine get_insolation_at_time( mesh, time, Q_TOA)
    ! Get monthly insolation at time t on the regional grid

    implicit none

    ! In/output variables:
    type(type_mesh),                           intent(in)  :: mesh
    real(dp),                                  intent(in)  :: time
    real(dp), dimension(mesh%vi1:mesh%vi2,12), intent(out) :: Q_TOA

    ! Local variables:
    character(len=256), parameter                          :: routine_name = 'get_insolation_at_time'
    real(dp)                                               :: time_applied
    integer                                                :: vi, m, ilat_l, ilat_u
    real(dp)                                               :: wt0, wt1, wlat_l, wlat_u
    real(dp), dimension(:,:), allocatable                  :: Q_TOA_int

    ! Add routine to path
    call init_routine( routine_name)

    time_applied = 0._dp

    select case (C%choice_insolation_forcing)

      case ('none')
        ! Do nothing and return
        call finalise_routine( routine_name)
        return

      case ('static')
        ! Use the same insolation always
        time_applied = C%static_insolation_time

      case ('realistic')
        ! Following a data record
        time_applied = time

      case default
        ! Unknown option
        call crash('unknown choice_insolation_forcing "' // &
                    TRIM( C%choice_insolation_forcing) // '"!')

    end select

    ! Allocate timeframe-interpolated lat-month-only insolation
    allocate( Q_TOA_int (1:forcing%ins_nlat, 12) )

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    if (time_applied < forcing%ins_t0 .or. time_applied > forcing%ins_t1) then
      call update_insolation_timeframes_from_file( time_applied)
    end if

    ! Calculate timeframe interpolation
    ! weights + safety net
    if (time_applied < forcing%ins_t0) then
      ! If applied time is still before the first
      ! timeframe (out of data record)
      wt0 = 1.0_dp
      wt1 = 0.0_dp
    elseif (time_applied > forcing%ins_t1) then
      ! If applied time is still after the second
      ! timeframe (out of data record)
      wt0 = 0.0_dp
      wt1 = 1.0_dp
    else
      ! If applied time is within the timeframes
      wt0 = (forcing%ins_t1 - time_applied) / &
            (forcing%ins_t1 - forcing%ins_t0)
      wt1 = 1._dp - wt0
    end if

    ! Interpolate the two timeframes
    Q_TOA_int = wt0 * forcing%ins_Q_TOA0 + &
                wt1 * forcing%ins_Q_TOA1

    ! Map the timeframe-interpolated lat-month-only
    ! insolation to the model mesh
    do vi = mesh%vi1, mesh%vi2

      ilat_l = FLOOR(mesh%lat( vi) + 91)
      ilat_u = ilat_l + 1

      wlat_l = forcing%ins_lat(ilat_u) - mesh%lat( vi)
      wlat_u = 1._dp - wlat_l

      do m = 1, 12
        Q_TOA( vi, m) = wlat_l * Q_TOA_int( ilat_l,m) + &
                        wlat_u * Q_TOA_int( ilat_u,m)
      end do

    end do

    ! Clean up after yourself
    deallocate( Q_TOA_int)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine get_insolation_at_time

  subroutine get_insolation_at_time_month_and_lat( time, month, lat, Q_TOA)
    ! Get monthly insolation at time t, month m and latitude l on the regional grid

    implicit none

    ! In/output variables:
    real(dp),                      intent(in)  :: time
    integer,                       intent(in)  :: month
    real(dp),                      intent(in)  :: lat
    real(dp),                      intent(out) :: Q_TOA

    ! Local variables:
    character(len=256), parameter              :: routine_name = 'get_insolation_at_time_month_and_lat'
    real(dp)                                   :: time_applied
    integer                                    :: ilat_l,ilat_u
    real(dp)                                   :: wt0, wt1, wlat_l, wlat_u

    ! Add routine to path
    call init_routine( routine_name)

    time_applied = 0._dp

    select case (C%choice_insolation_forcing)

      case ('none')
        ! Do nothing and return
        call finalise_routine( routine_name)
        return

      case ('static')
        ! Use the same insolation always
        time_applied = C%static_insolation_time

      case ('realistic')
        ! Following a data record
        time_applied = time

      case default
        ! Unknown option
        call crash('unknown choice_insolation_forcing "' // &
                    TRIM( C%choice_insolation_forcing) // '"!')

    end select

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    if (time_applied < forcing%ins_t0 .and. time_applied > forcing%ins_t1) then
      call update_insolation_timeframes_from_file( time_applied)
    end if

    ! Calculate timeframe interpolation weights
    wt0 = (forcing%ins_t1 - time) / (forcing%ins_t1 - forcing%ins_t0)
    wt1 = 1._dp - wt0

    ! Get value at month m and latitude l
    ilat_l = floor(lat + 91)
    ilat_u = ilat_l + 1

    wlat_l = (forcing%ins_lat( ilat_u) - lat) / &
             (forcing%ins_lat( ilat_u) - forcing%ins_lat( ilat_l))
    wlat_u = 1._dp - wlat_l

    Q_TOA = wt0 * wlat_l * forcing%ins_Q_TOA0( ilat_l,month) + &
            wt0 * wlat_u * forcing%ins_Q_TOA0( ilat_u,month) + &
            wt1 * wlat_l * forcing%ins_Q_TOA1( ilat_l,month) + &
            wt1 * wlat_u * forcing%ins_Q_TOA1( ilat_u,month)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine get_insolation_at_time_month_and_lat

  subroutine update_insolation_timeframes_from_file( time)
    ! Read the NetCDF file containing the insolation forcing data. Only read the time frames enveloping the current
    ! coupling timestep to save on memory usage. Only done by master.

    ! NOTE: assumes time in forcing file is in kyr

    implicit none

    real(dp),                      intent(in) :: time

    ! Local variables:
    character(len=256), parameter             :: routine_name = 'update_insolation_timeframes_from_file'
    integer                                   :: ti0, ti1

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_insolation_forcing)

      case ('none')
        ! Do nothing and return
        call finalise_routine( routine_name)
        return

      case ('static','realistic')
        ! Update insolation

        ! Find time indices to be read
        if (time <= forcing%ins_time( forcing%ins_nyears)) then
          ti1 = 1
          do while (forcing%ins_time(ti1) < time)
            ti1 = ti1 + 1
          end do
          ti0 = ti1 - 1

          forcing%ins_t0 = forcing%ins_time(ti0)
          forcing%ins_t1 = forcing%ins_time(ti1)
        else
          ! Constant PD insolation for future projections
          ti0 = forcing%ins_nyears
          ti1 = forcing%ins_nyears

          forcing%ins_t0 = forcing%ins_time(ti0) - 1._dp
          forcing%ins_t1 = forcing%ins_time(ti1)
        end if

        ! Read new insolation fields from the NetCDF file
        call read_insolation_file_timeframes( forcing, ti0, ti1)

      case default
        ! Unknown option
        call crash('unknown choice_insolation_forcing "' // &
                    TRIM( C%choice_insolation_forcing) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_insolation_timeframes_from_file

  subroutine initialise_insolation_data
    ! Allocate shared memory for the forcing data fields

    implicit none

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_insolation_data'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Read data ===
    ! =================

    select case (C%choice_insolation_forcing)

      case ('none')
        ! Do nothing; no insolation included

      case ('static','realistic')
        ! Initialise insolation

        if (par%master) then
          write(*,"(3A)") ' Initialising insolation data from ', &
                            trim(C%filename_insolation), '...'
        end if
        call sync

        ! Name of file containing record
        forcing%netcdf_ins%filename = C%filename_insolation

        ! Initialise to impossible times to force a first read
        forcing%ins_t0 = C%start_time_of_run - 100._dp
        forcing%ins_t1 = C%start_time_of_run - 90._dp

        ! Inquire into the insolation forcing netcdf file
        call inquire_insolation_file( forcing)

        ! Allocation based on inquired data
        allocate( forcing%ins_time   (1:forcing%ins_nyears)   )
        allocate( forcing%ins_lat    (1:forcing%ins_nlat)     )
        allocate( forcing%ins_Q_TOA0 (1:forcing%ins_nlat, 12) )
        allocate( forcing%ins_Q_TOA1 (1:forcing%ins_nlat, 12) )

        ! Read time and latitude data
        call read_insolation_file_time_lat( forcing)

        ! Time-out-of-record warning/crash
        if (par%master) then
          if (C%start_time_of_run < forcing%ins_time(1)) then
            call crash('Model time starts before start of insolation' // &
                       ' record; the model will crash lol')
          end if
          if (C%end_time_of_run > forcing%ins_time(forcing%ins_nyears)) then
            call warning('Model time will reach beyond end of insolation record;' // &
                         ' constant extrapolation will be used in that case!')
          end if
        end if
        call sync

      case default
        ! Unknown option
        call crash('unknown choice_insolation_forcing "' // &
                    trim( C%choice_insolation_forcing) // '"!')

    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  END SUBROUTINE initialise_insolation_data

! ===== Geothermal heat flux =====
! ================================

  subroutine initialise_geothermal_heat_flux_global
    ! Initialise global geothermal heat flux data

    implicit none

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_geothermal_heat_flux_global'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_geothermal_heat_flux)

      case ('constant')
        ! Use constant value, no need to read a file

      case ('spatial')
        ! Use a spatially variable geothermal heat fux
        ! read from the specified NetCDF file

        if (par%master) then
          write(*,"(3A)") ' Initialising geothermal heat flux data from ', &
                            TRIM(C%filename_geothermal_heat_flux), '...'
        end if
        call sync

        ! Name of file containing data
        forcing%netcdf_ghf%filename = C%filename_geothermal_heat_flux

        ! Read size of data field
        call inquire_geothermal_heat_flux_file( forcing)

        ! Allocation based on inquired data
        allocate( forcing%grid_ghf%lon (1:forcing%grid_ghf%nlon)                          )
        allocate( forcing%grid_ghf%lat (1:forcing%grid_ghf%nlat)                          )
        allocate( forcing%ghf_ghf      (1:forcing%grid_ghf%nlon, 1:forcing%grid_ghf%nlat) )

        ! Read GHF and lat/lon data
        call read_geothermal_heat_flux_file( forcing)

      case default
        ! Unknown option
        call crash('unknown choice_geothermal_heat_flux "' // &
                    TRIM( C%choice_geothermal_heat_flux) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_geothermal_heat_flux_global

! ===== Sea level records =====
! =============================

  subroutine initialise_sealevel_record
    ! Read the sea level record specified in C%filename_sealevel_record.
    ! Assumes this is an ASCII text file with at least two columns (time
    ! in kyr and sea level in m) and the number of rows being equal to
    ! C%sealevel_record_length
    ! NOTE: assumes time is listed in YEARS (so LGM would be -21000.0)

    implicit none

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_sealevel_record'
    integer                       :: i,ios, fp

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory to take the data
    allocate(forcing%sealevel_time  ( C%sealevel_record_length ))
    allocate(forcing%sealevel_record( C%sealevel_record_length ))

    if (par%master) then
      write(*,"(3A)") ' Initialising sea level record from ', &
                        TRIM(C%filename_sealevel_record), '...'
    end if
    call sync

    ! === Read data ===
    ! =================

    ! Read CO2 record (time and values) from specified text file
    open(  newunit=fp, file=C%filename_sealevel_record, action='READ')

    do i = 1, C%sealevel_record_length
      read( unit=fp, fmt=*, iostat=ios) forcing%sealevel_time(i), &
                                        forcing%sealevel_record(i)

      if (ios /= 0) then
        call crash('length of text file "' // &
                    TRIM(C%filename_sealevel_record) // &
                    '" does not match C%sealevel_record_length!')
      end if

    end do

    close( unit=fp)

    ! === Update sea level ===
    ! ========================

    ! Set the value for the current (starting) model time
    call update_sealevel_at_model_time( C%start_time_of_run)

    ! === Finalisation ===
    ! ====================

    if (par%master) then
      if (C%start_time_of_run < forcing%sealevel_time(1)) then
        call warning('Model time starts before start of sea level record;' // &
                     ' constant extrapolation will be used in that case!')
      end if

      if (C%end_time_of_run > forcing%sealevel_time(C%sealevel_record_length)) then
        call warning('Model time will reach beyond end of sea level record;' // &
                     ' constant extrapolation will be used in that case!')
      end if
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_sealevel_record

  subroutine update_sealevel_at_model_time( time)
    ! Interpolate the data in forcing%sealevel to find the
    ! value at the queried time. If time lies outside the
    ! range of forcing%sealevel_time, return the first/last value.
    ! NOTE: assumes time is listed in YEARS (so LGM would be -21000.0)

    implicit none

    ! In/output variables:
    real(dp),                      intent(in) :: time

    ! Local variables:
    character(len=256), parameter             :: routine_name = 'update_sealevel_record_at_model_time'
    integer                                   :: il, iu
    real(dp)                                  :: wl, wu

    ! Add routine to path
    call init_routine( routine_name)

    if (time <= minval( forcing%sealevel_time)) then
      ! Model time before start of sea level record; using constant extrapolation
      forcing%sealevel_obs = forcing%sealevel_record( 1)

    elseif (time >= maxval( forcing%sealevel_time)) then
      ! Model time beyond end of sea level record; using constant extrapolation
      forcing%sealevel_obs = forcing%sealevel_record( C%sealevel_record_length)

    else

      iu = 1
      do while (forcing%sealevel_time(iu) < time)
        iu = iu+1
      end do
      il = iu-1

      wl = (forcing%sealevel_time(iu) - time) / ((forcing%sealevel_time(iu)-forcing%sealevel_time(il)))
      wu = 1._dp - wl

      forcing%sealevel_obs = forcing%sealevel_record(il) * wl + forcing%sealevel_record(iu) * wu

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_sealevel_at_model_time

end module forcing_module
