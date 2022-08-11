module ocean_module

! ===== Preamble =====
! ====================

  use configuration_module, only : dp, C, init_routine, finalise_routine, crash
  use parallel_module,      only : par, partition_list
  use data_types_module,    only : type_ocean_matrix_global, type_ocean_snapshot_global
  use netcdf_module,        only : inquire_PD_obs_global_ocean_file, read_PD_obs_global_ocean_file
  use utilities_module,     only : remap_cons_2nd_order_1D

contains

! ===== Main routines =====
! =========================

  subroutine initialise_ocean_model_global( ocean_matrix)
    ! Initialise the global ocean model

    implicit none

    ! In/output variables:
    type(type_ocean_matrix_global), intent(inout) :: ocean_matrix

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'initialise_ocean_model_global'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write (*,'(3A)') ' Initialising global ocean model "', TRIM(C%choice_ocean_model), '"...'
    end if

    ! Initialise ocean vertical grid
    CALL initialise_ocean_vertical_grid

    ! Pick selected method
    select case(C%choice_ocean_model)

      case('none')
        ! No need to do anything

      case('PD_obs')
        ! Keep the ocean fixed to present-day observed conditions
        call initialise_ocean_model_global_PD_obs( ocean_matrix%PD_obs)

      case('matrix_warm_cold')
        ! Allocate all global snapshots used in the warm/cold ocean matrix
        call crash(TRIM(C%choice_ocean_model) // ' not implemented yet...')

      case default
        ! Unknown option
        CALL crash('unknown choice_ocean_model "' // TRIM( C%choice_ocean_model) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_global

! ===== Observational PD ocean =====
! ==================================

  subroutine initialise_ocean_model_global_PD_obs( PD_obs)
    ! Initialise the present-day observed global ocean snapshot

    implicit none

    ! In/output variables:
    type(type_ocean_snapshot_global), intent(inout) :: PD_obs

    ! Local variables:
    character(len=256), parameter                   :: routine_name = 'initialise_ocean_model_global_PD_obs'

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise the present-day observed global climate
    PD_obs%name = 'PD_obs'
    PD_obs%netcdf%filename = C%filename_PD_obs_ocean

    ! Inquire data from the NetCDF file + get grid size (nlon, nlat, nz)
    call inquire_PD_obs_global_ocean_file( PD_obs)

    ! Allocate memory
    allocate( PD_obs%lon         (1:PD_obs%nlon)                                       )
    allocate( PD_obs%lat         (1:PD_obs%nlat)                                       )
    allocate( PD_obs%z_ocean_raw (1:PD_obs%nz_ocean_raw)                               )
    allocate( PD_obs%T_ocean_raw (1:PD_obs%nlon, 1:PD_obs%nlat, 1:PD_obs%nz_ocean_raw) )
    allocate( PD_obs%S_ocean_raw (1:PD_obs%nlon, 1:PD_obs%nlat, 1:PD_obs%nz_ocean_raw) )

    if (par%master) then
      write(*,'(3A)') '  Reading PD observed ocean data from file ', &
                         TRIM(PD_obs%netcdf%filename), '...'
    end if

    ! Read data from the NetCDF file
    call read_PD_obs_global_ocean_file( PD_obs)

    ! Map the data to the desired vertical grid
    call map_global_ocean_data_to_UFEMISM_vertical_grid( PD_obs)

    ! Deallocate raw ocean data
    deallocate( PD_obs%z_ocean_raw )
    deallocate( PD_obs%T_ocean_raw )
    deallocate( PD_obs%S_ocean_raw )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_global_PD_obs

! ===== Regridding ======
! =======================

  subroutine initialise_ocean_vertical_grid
    ! Set up the vertical grid used for ocean data

    implicit none

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_ocean_vertical_grid'

    ! Add routine to path
    call init_routine( routine_name)

    ! Pick selected method
    select case(C%choice_ocean_vertical_grid)

      case('regular')
        ! Regular vertical grid
        call initialise_ocean_vertical_grid_regular

      case default
        ! Unknown option
        call crash('unknown choice_ocean_vertical_grid "' // TRIM( C%choice_ocean_vertical_grid) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_vertical_grid

  subroutine initialise_ocean_vertical_grid_regular
    ! Set up the vertical grid used for ocean data - regular grid

    implicit none

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_ocean_vertical_grid_regular'
    INTEGER                       :: k

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the number of vertical layers to be used
    C%nz_ocean = 1 + floor( C%ocean_vertical_grid_max_depth / C%ocean_regular_grid_dz)

    ! Allocate memory
    allocate( C%z_ocean( C%nz_ocean))

    ! Fill in the values
    do k = 1, C%nz_ocean
      C%z_ocean( k) = real(k-1,dp) * C%ocean_regular_grid_dz
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_vertical_grid_regular

  subroutine map_global_ocean_data_to_UFEMISM_vertical_grid( ocean_glob)
    ! Map global 3-D ocean temperature and salinity from whatever vertical grid
    ! it has been provided on to the vertical grid used by UFEMISM

    implicit none

    ! Input variables:
    type(type_ocean_snapshot_global), intent(inout) :: ocean_glob

    ! Local variables:
    character(len=256), parameter      :: routine_name = 'map_global_ocean_data_to_UFEMISM_vertical_grid'
    integer                            :: i,j,k
    integer, dimension(:), allocatable :: z_mask_old, z_mask_new
    real(dp)                           :: z_floor
    real(dp)                           :: NaN

    ! Add routine to path
    call init_routine( routine_name)

    ! A trick
    NaN = -1._dp
    NaN = sqrt( NaN)

    ! Allocate shared memory for the global 3-D ocean temperature and salinity
    allocate(ocean_glob%T_ocean (1:ocean_glob%nlon, 1:ocean_glob%nlat, 1:C%nz_ocean))
    allocate(ocean_glob%S_ocean (1:ocean_glob%nlon, 1:ocean_glob%nlat, 1:C%nz_ocean))

    ! Use "masked" 2-nd order conservative 1-D remapping
    allocate( z_mask_old( 1:ocean_glob%nz_ocean_raw))
    allocate( z_mask_new( 1:C%nz_ocean             ))

    do j = 1, ocean_glob%nlat
    do i = 1, ocean_glob%nlon

      ! Determine local depth of the ocean floor, fill in both data masks
      if (ocean_glob%T_ocean_raw( i,j,ocean_glob%nz_ocean_raw) == &
          ocean_glob%T_ocean_raw( i,j,ocean_glob%nz_ocean_raw)) then

        ! Ocean floor lies below the vertical limit of the provided data
        z_mask_old = 1
        z_floor = ocean_glob%z_ocean_raw( ocean_glob%nz_ocean_raw) + &
                  (ocean_glob%z_ocean_raw( 2) - ocean_glob%z_ocean_raw( 1))

      elseif (ocean_glob%T_ocean_raw( i,j,1) /= ocean_glob%T_ocean_raw( i,j,1)) then

        ! This grid cell isn't ocean at all
        z_mask_old = 0
        z_floor    = 0._dp

      else

        z_mask_old = 1
        k = ocean_glob%nz_ocean_raw

        do while (ocean_glob%T_ocean_raw( i,j,k) /= ocean_glob%T_ocean_raw( i,j,k))
          z_mask_old( k) = 0
          z_floor = ocean_glob%z_ocean_raw( k)
          k = k - 1
        end do

      end if

      z_mask_new = 0

      do k = 1, C%nz_ocean
        if (C%z_ocean( k) < z_floor) then
          z_mask_new = 1
        end if
      end do

      ! Regrid vertical column
      call remap_cons_2nd_order_1D( ocean_glob%z_ocean_raw, z_mask_old, &
                                    ocean_glob%T_ocean_raw( i,j,:), &
                                    C%z_ocean, z_mask_new, &
                                    ocean_glob%T_ocean( i,j,:))

      call remap_cons_2nd_order_1D( ocean_glob%z_ocean_raw, z_mask_old, &
                                    ocean_glob%S_ocean_raw( i,j,:), &
                                    C%z_ocean, z_mask_new, &
                                    ocean_glob%S_ocean( i,j,:))

      ! Fill masked values with NaN
      do k = 1, C%nz_ocean
        if (z_mask_new( k) == 0) then
          ocean_glob%T_ocean( i,j,k) = NaN
          ocean_glob%S_ocean( i,j,k) = NaN
        end if
      end do

    end do
    end do

    ! Clean up after yourself
    deallocate( z_mask_old)
    deallocate( z_mask_new)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_global_ocean_data_to_UFEMISM_vertical_grid

end module