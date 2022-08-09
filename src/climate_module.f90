module climate_module
  ! Contains all the routines related to the global and regional climate forcing

! ===== Preamble =====
! ====================

  use configuration_module, only : dp, C, init_routine, finalise_routine, crash
  use parallel_module,      only : par, partition_list
  use data_types_module,    only : type_climate_matrix_global, type_climate_snapshot_global
  use netcdf_module,        only : inquire_PD_obs_global_climate_file, read_PD_obs_global_climate_file

  implicit none

contains

! ===== Main routines =====
! =========================

  subroutine initialise_climate_model_global( climate_matrix)
    ! Initialise the global climate model

    implicit none

    ! In/output variables:
    type(type_climate_matrix_global), intent(inout) :: climate_matrix

    ! Local variables:
    character(len=256), parameter                   :: routine_name = 'initialise_climate_model_global'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write(*,'(3A)') ' Initialising global climate model "', trim(C%choice_climate_model), '"...'
    end if

    ! Pick selected method
    select case (C%choice_climate_model)

      case('none')
        ! No need to do anything

      case('PD_obs')
        ! Keep the climate fixed to present-day observed conditions
        call initialise_climate_model_global_PD_obs( climate_matrix%PD_obs)

      case('matrix_warm_cold')
        ! Allocate all global snapshots used in the warm/cold climate matrix
        call crash(TRIM(C%choice_climate_model) // ' not implemented yet...')

      case default
        ! Unknown option
        call crash('unknown choice_climate_model"' // TRIM(C%choice_climate_model) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_climate_model_global

! ===== Observational PD climate =====
! ====================================

  subroutine initialise_climate_model_global_PD_obs( PD_obs)
    ! Initialise the observational present-day global climate model

    implicit none

    ! In/output variables:
    type(type_climate_snapshot_global), intent(inout) :: PD_obs

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'initialise_climate_model_global_PD_obs'

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise the present-day observed global climate
    PD_obs%name = 'PD_obs'
    PD_obs%netcdf%filename = C%filename_PD_obs_climate

    ! Inquire data from the NetCDF file + get grid size (nlon, nlat)
    call inquire_PD_obs_global_climate_file( PD_obs)

    ! Allocate memory
    allocate( PD_obs%lon     (1:PD_obs%nlon)                    )
    allocate( PD_obs%lat     (1:PD_obs%nlat)                    )
    allocate( PD_obs%Hs      (1:PD_obs%nlon, 1:PD_obs%nlat)     )
    allocate( PD_obs%T2m     (1:PD_obs%nlon, 1:PD_obs%nlat, 12) )
    allocate( PD_obs%Precip  (1:PD_obs%nlon, 1:PD_obs%nlat, 12) )
    allocate( PD_obs%Wind_WE (1:PD_obs%nlon, 1:PD_obs%nlat, 12) )
    allocate( PD_obs%Wind_SN (1:PD_obs%nlon, 1:PD_obs%nlat, 12) )

    if (par%master) then
      write(*,'(3A)') '  Reading PD observed climate data from file ', &
                         TRIM(PD_obs%netcdf%filename), '...'
    end if

    ! Read data from the NetCDF file
    call read_PD_obs_global_climate_file( PD_obs)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_climate_model_global_PD_obs

end module