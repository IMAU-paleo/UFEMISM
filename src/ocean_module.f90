module ocean_module

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module, only : dp, C, init_routine, finalise_routine, crash
  use parallel_module,      only : par, partition_list, cerr, ierr
  use data_types_module,    only : type_ocean_matrix_global, type_ocean_snapshot_global, &
                                   type_ocean_matrix_regional, type_ocean_snapshot_regional, &
                                   type_model_region, type_mesh, type_highres_ocean_data
  use netcdf_module,        only : inquire_PD_obs_global_ocean_file, read_PD_obs_global_ocean_file, &
                                   inquire_hires_geometry_file, read_hires_geometry_file
  use utilities_module,     only : remap_cons_2nd_order_1D, inverse_oblique_sg_projection, &
                                   map_glob_to_grid_3D, extrapolate_Gaussian_floodfill
  use mesh_mapping_module,  only : calc_remapping_operator_mesh2grid, map_mesh2grid_2D, &
                                   calc_remapping_operator_grid2mesh, map_grid2mesh_3D, &
                                   deallocate_remapping_operators_mesh2grid
  use utilities_module,     only : surface_elevation

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
      write (*,'(3A)') ' Initialising global ocean model "', trim(C%choice_ocean_model), '"...'
    end if

    ! Initialise ocean vertical grid
    call initialise_ocean_vertical_grid

    ! Pick selected method
    select case(C%choice_ocean_model)

      case('none')
        ! No need to do anything

      case('PD_obs')
        ! Keep the ocean fixed to present-day observed conditions
        call initialise_ocean_model_global_PD_obs( ocean_matrix%PD_obs)

      case('matrix_warm_cold')
        ! Allocate all global snapshots used in the warm/cold ocean matrix
        call crash(trim(C%choice_ocean_model) // ' not implemented yet...')

      case default
        ! Unknown option
        call crash('unknown choice_ocean_model "' // trim( C%choice_ocean_model) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_global

  subroutine initialise_ocean_model_regional( region, ocean_matrix_global)
    ! Initialise the regional ocean model

    implicit none

    ! In/output variables:
    type(type_model_region),        intent(inout) :: region
    type(type_ocean_matrix_global), intent(in)    :: ocean_matrix_global

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'initialise_ocean_model_regional'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write (*,'(3A)') '  Initialising regional ocean model "', trim(C%choice_ocean_model), '"...'
    end if

    select case (C%choice_ocean_model)

    case ('none')
      ! No need to do anything

    case ('PD_obs')
      ! Keep the ocean fixed to present-day observed conditions
      call initialise_ocean_model_PD_obs_regional( region, ocean_matrix_global)

    case ('matrix_warm_cold')
      ! Allocate all global snapshots used in the warm/cold ocean matrix
      call crash(trim(C%choice_ocean_model) // ' not implemented yet...')

    case default
      ! Unknown option
      call crash('unknown choice_ocean_model"' // trim(C%choice_ocean_model) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected=37)

  end subroutine initialise_ocean_model_regional

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
                         trim(PD_obs%netcdf%filename), '...'
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

  subroutine initialise_ocean_model_PD_obs_regional( region, ocean_matrix_global)
    ! Initialise the regional ocean model
    !
    ! Just use the present-day observed ocean data

    implicit none

    ! In/output variables:
    type(type_model_region),        intent(inout) :: region
    type(type_ocean_matrix_global), intent(in)    :: ocean_matrix_global

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'initialise_ocean_model_PD_obs_regional'
    integer                                       :: vi,k

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate all the snapshots
    call allocate_ocean_snapshot_regional( region%mesh, region%ocean_matrix%PD_obs,  name = 'PD_obs' )
    call allocate_ocean_snapshot_regional( region%mesh, region%ocean_matrix%applied, name = 'applied')

    ! Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid,
    ! extrapolate mapped ocean data to cover the entire 3D domain, and finally
    ! map to the actual ice model resolution
    call get_extrapolated_ocean_data( region, ocean_matrix_global%PD_obs, region%ocean_matrix%PD_obs, C%filename_PD_obs_ocean)

    do k = 1, C%nz_ocean
    do vi = region%mesh%vi1, region%mesh%vi2

      ! PD_obs doesn't have a bias-corrected version
      region%ocean_matrix%PD_obs%T_ocean_corr_ext(  vi,k) = region%ocean_matrix%PD_obs%T_ocean_ext( vi,k)
      region%ocean_matrix%PD_obs%S_ocean_corr_ext(  vi,k) = region%ocean_matrix%PD_obs%S_ocean_ext( vi,k)

      ! Initialise applied ocean forcing with present-day observations
      region%ocean_matrix%applied%T_ocean(          vi,k) = region%ocean_matrix%PD_obs%T_ocean(          vi,k)
      region%ocean_matrix%applied%T_ocean_ext(      vi,k) = region%ocean_matrix%PD_obs%T_ocean_ext(      vi,k)
      region%ocean_matrix%applied%T_ocean_corr_ext( vi,k) = region%ocean_matrix%PD_obs%T_ocean_corr_ext( vi,k)
      region%ocean_matrix%applied%S_ocean(          vi,k) = region%ocean_matrix%PD_obs%S_ocean(          vi,k)
      region%ocean_matrix%applied%S_ocean_ext(      vi,k) = region%ocean_matrix%PD_obs%S_ocean_ext(      vi,k)
      region%ocean_matrix%applied%S_ocean_corr_ext( vi,k) = region%ocean_matrix%PD_obs%S_ocean_corr_ext( vi,k)

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected=14)

  end subroutine initialise_ocean_model_PD_obs_regional

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
        call crash('unknown choice_ocean_vertical_grid "' // trim( C%choice_ocean_vertical_grid) // '"!')

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

! ===== Extrapolation =====
! =========================

  subroutine get_extrapolated_ocean_data( region, ocean_glob, ocean_reg, filename_ocean_glob)
    ! Check if extrapolated ocean files for the current ice model setting exist. If so,
    ! read those. If not, perform the extrapolation and save the results to a new netCDF file.

    ! When creating a set of extrapolated files, a header file is created that describes
    ! the ice model settings for which those files were created. We check all existing header
    ! files, if any of them match the current settings, we read the extrapolated files listed
    ! there. If none of them match, we create a set of extrapolated files (and the accompanying
    ! header file) from scratch.

    implicit none

    ! In/output variables:
    type(type_model_region),            intent(inout) :: region
    type(type_ocean_snapshot_global),   intent(in)    :: ocean_glob
    type(type_ocean_snapshot_regional), intent(inout) :: ocean_reg
    character(len=256),                 intent(in)    :: filename_ocean_glob

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'get_extrapolated_ocean_data'
    logical                                           :: foundmatch
    type(type_highres_ocean_data)                     :: hires
    integer                                           :: i,j,n
    real(dp), parameter                               :: tol = 1E-9_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! ===== Initial check =====
    ! =========================

    ! First, check if any existing ocean header matches the current ice model set-up.
    call check_for_matching_ocean_header( region, filename_ocean_glob, foundmatch, ocean_reg%hires_ocean_foldername)

    ! If a valid preprocessed file exists, read data from there. If not, perform
    ! the preprocessing and save the result to a file to save on future work

    ! if (foundmatch) then

    !   if (par%master) WRITE(0,*) '   Found valid extrapolated ocean data in folder "', trim( ocean_reg%hires_ocean_foldername), '"'
    !   call get_hires_ocean_data_from_file( region%mesh, hires, ocean_reg%hires_ocean_foldername)

    ! else

      ! No header fitting the current ice model set-up was found. Create a new one describing
      ! the current set-up, and generate extrapolated ocean data files from scratch.
      if (par%master) then
        write(*,"(3A)") '   Creating new extrapolated ocean data in folder "', trim( ocean_reg%hires_ocean_foldername), '"...'
      end if

      call map_and_extrapolate_hires_ocean_data( region, ocean_glob, hires)
      ! call write_hires_extrapolated_ocean_data_to_file( hires, filename_ocean_glob, ocean_reg%hires_ocean_foldername)

    ! end if

    ! ===== Map extrapolated data from the high-resolution grid to the actual ice-model mesh =====
    ! ============================================================================================

    if (par%master) then
      write(*,"(A)") '    Mapping high-resolution extrapolated ocean data unto the ice-model mesh...'
    end if

    ! Tolerance; points lying within this distance of each other are treated as identical
    hires%grid%tol_dist = ((hires%grid%xmax - hires%grid%xmin) + (hires%grid%ymax - hires%grid%ymin)) * tol / 2._dp

    ! Set up grid-to-vector translation tables
    hires%grid%n  = hires%grid%nx * hires%grid%ny

    allocate( hires%grid%ij2n (hires%grid%nx, hires%grid%ny))
    allocate( hires%grid%n2ij (hires%grid%n , 2            ))

    if (par%master) then
      n = 0
      do i = 1, hires%grid%nx
        if (mod(i,2) == 1) then
          do j = 1, hires%grid%ny
            n = n+1
            hires%grid%ij2n( i,j) = n
            hires%grid%n2ij( n,:) = [i,j]
          end do
        else
          do j = hires%grid%ny, 1, -1
            n = n+1
            hires%grid%ij2n( i,j) = n
            hires%grid%n2ij( n,:) = [i,j]
          end do
        end if
      end do
    end if

    ! Map high-resolution ocean data to UFEMISM's mesh
    call calc_remapping_operator_grid2mesh( hires%grid, region%mesh)
    call map_grid2mesh_3D( hires%grid, region%mesh, hires%T_ocean, ocean_reg%T_ocean_ext)
    call map_grid2mesh_3D( hires%grid, region%mesh, hires%S_ocean, ocean_reg%S_ocean_ext)
    call deallocate_remapping_operators_mesh2grid( hires%grid)

    ! Clean up after yourself
    deallocate( hires%grid%x              )
    deallocate( hires%grid%y              )
    deallocate( hires%grid%lat            )
    deallocate( hires%grid%lon            )
    deallocate( hires%grid%ij2n           )
    deallocate( hires%grid%n2ij           )
    deallocate( hires%T_ocean             )
    deallocate( hires%S_ocean             )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine get_extrapolated_ocean_data

  subroutine check_for_matching_ocean_header( region, filename_ocean_glob, foundmatch, hires_foldername)
    ! Inspect all the folder inside the "extrapolated_ocean_files" folder, read their header files,
    ! and see if any of them match the current model settings. If so, return the name of the folder
    ! where the matching header was found. If not, return the name of the folder where a new header
    ! should be created.

    implicit none

    ! In/output variables:
    type(type_model_region), intent(in)  :: region
    character(LEN=256),      intent(in)  :: filename_ocean_glob
    logical,                 intent(out) :: foundmatch
    character(LEN=256),      intent(out) :: hires_foldername

    ! Local variables:
    character(len=256), parameter        :: routine_name = 'check_for_matching_ocean_header'
    integer                              :: cerr, ierr
    character(len=256)                   :: header_filename
    logical                              :: header_exists
    integer                              :: folder_i

    character(len=256)                   :: original_ocean_filename_read
    character(len=256)                   :: choice_ocean_vertical_grid_read
    integer                              :: nz_ocean_read
    real(dp)                             :: ocean_vertical_grid_max_depth_read
    real(dp)                             :: ocean_extrap_res_read
    real(dp)                             :: ocean_extrap_Gauss_sigma_read
    real(dp)                             :: lambda_M_read
    real(dp)                             :: phi_M_read
    real(dp)                             :: alpha_stereo_read

    ! Add routine to path
    call init_routine( routine_name)

    folder_i = 1
    do while (folder_i < 1000)
      ! Generate a foldername to inspect

      if     (folder_i < 10)   then
        WRITE( hires_foldername,'(A,A,A,A,I1)') trim(C%ocean_extrap_dir), '/', region%name, '_00', folder_i
      elseif (folder_i < 100)  then
        WRITE( hires_foldername,'(A,A,A,A,I2)') trim(C%ocean_extrap_dir), '/', region%name, '_0',  folder_i
      elseif (folder_i < 1000) then
        WRITE( hires_foldername,'(A,A,A,A,I3)') trim(C%ocean_extrap_dir), '/', region%name, '_',   folder_i
      else
        call crash('tried a thousand folders!')
      end if

      ! Check if a header in this folder exists. If not, then we've inspected all existing headers
      ! without finding the good one, so we must generate the extrapolated ocean files from scratch.
      header_filename = trim( hires_foldername) // '/' // 'header.txt'

      inquire( file = header_filename, exist = header_exists)

      if (.not. header_exists) then
        ! No more headers exist to be inspected. Return foundmatch = .false.

        foundmatch = .false.
        exit

      else
        ! If the header exists, read it and see if it fits the current ice-model set-up.

        call read_ocean_header( header_filename, original_ocean_filename_read, &
                                choice_ocean_vertical_grid_read, nz_ocean_read, &
                                ocean_vertical_grid_max_depth_read, ocean_extrap_res_read, &
                                ocean_extrap_Gauss_sigma_read, lambda_M_read, &
                                phi_M_read, alpha_stereo_read)

        if ( trim(original_ocean_filename_read)    == trim(filename_ocean_glob)             .and. &
             trim(choice_ocean_vertical_grid_read) == trim(choice_ocean_vertical_grid_read) .and. &
             nz_ocean_read                         == C%nz_ocean                            .and. &
             ocean_vertical_grid_max_depth_read    == C%ocean_vertical_grid_max_depth       .and. &
             ocean_extrap_res_read                 == C%ocean_extrap_res                    .and. &
             ocean_extrap_Gauss_sigma_read         == C%ocean_extrap_Gauss_sigma            .and. &
             lambda_M_read                         == region%mesh%lambda_M                  .and. &
             phi_M_read                            == region%mesh%phi_M                     .and. &
             alpha_stereo_read                     == region%mesh%alpha_stereo) then
          ! This header matches the current model set-up!

          foundmatch = .true.
          exit

        else
          ! This header doesn't match the current model set-up. Try the next one.

          folder_i = folder_i + 1

        end if
      end if

    end do ! while (header_i < 1000)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_for_matching_ocean_header

  subroutine read_ocean_header( header_filename, original_ocean_filename, &
                                choice_ocean_vertical_grid, nz_ocean, &
                                ocean_vertical_grid_max_depth, &
                                ocean_extrap_res, ocean_extrap_Gauss_sigma, &
                                lambda_M, phi_M, alpha_stereo)
    ! Read a header file listing the model settings that were used to create a high-resolution extrapolated ocean data file

    implicit none

    ! In/output variables:
    character(len=256), intent(in)  :: header_filename
    character(len=256), intent(out) :: original_ocean_filename
    character(len=256), intent(out) :: choice_ocean_vertical_grid
    integer,            intent(out) :: nz_ocean
    real(dp),           intent(out) :: ocean_vertical_grid_max_depth
    real(dp),           intent(out) :: ocean_extrap_res
    real(dp),           intent(out) :: ocean_extrap_Gauss_sigma
    real(dp),           intent(out) :: lambda_M
    real(dp),           intent(out) :: phi_M
    real(dp),           intent(out) :: alpha_stereo

    ! Local variables:
    character(len=256), parameter   :: routine_name = 'read_ocean_header'
    integer                         :: ios, cerr, ierr

    ! The NAMELIST that's used to read the external header file.
    namelist /HEADER/original_ocean_filename,       &
                     choice_ocean_vertical_grid,    &
                     nz_ocean,                      &
                     ocean_vertical_grid_max_depth, &
                     ocean_extrap_res,              &
                     ocean_extrap_Gauss_sigma,      &
                     lambda_M,                      &
                     phi_M,                         &
                     alpha_stereo

    ! Add routine to path
    call init_routine( routine_name)

    open( unit = 29, file = trim(header_filename), status='old', action='read', iostat=ios)

    if (ios /= 0) then
      call crash('could not open "' // trim(header_filename) // '"!')
    end if

    ! In the following statement the entire configuration file is read, using the namelist (NML=HEADER)
    read(  unit = 29, nml = HEADER, iostat = ios)
    close( unit = 29)

    if (ios /= 0) then
      call crash('could not read "' // trim(header_filename) // '"!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_ocean_header

  subroutine map_and_extrapolate_hires_ocean_data( region, ocean_glob, hires)
    ! Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid,
    ! extrapolate mapped ocean data to cover the entire 3D domain, and finally
    ! map to the actual ice model resolution.

    implicit none

    ! In/output variables:
    type(type_model_region),          intent(inout) :: region
    type(type_ocean_snapshot_global), intent(in)    :: ocean_glob
    type(type_highres_ocean_data),    intent(inout) :: hires

    ! Local variables:
    character(len=256), parameter                   :: routine_name = 'map_and_extrapolate_hires_ocean_data'
    integer                                         :: vi,i,j
    real(dp), dimension(:  ), allocatable           :: basin_ID_dp_lores
    real(dp), dimension(:,:), allocatable           :: basin_ID_dp_hires,  basin_ID_dp_hires_ext
    integer                                         :: ii,jj,n
    logical                                         :: foundit
    real(dp), PARAMETER                             :: tol = 1E-9_dp
    character(len=256)                              :: choice_basin_scheme
    character(len=256)                              :: filename_basins

    ! Add routine to path
    call init_routine( routine_name)

    ! ===== High-resolution grid =====
    ! ================================

    ! Determine which file to use for this region
    if     (region%name == 'NAM') then
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_NAM
      choice_basin_scheme       = C%choice_basin_scheme_NAM
      filename_basins           = C%filename_basins_NAM
    elseif (region%name == 'EAS') then
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_EAS
      choice_basin_scheme       = C%choice_basin_scheme_EAS
      filename_basins           = C%filename_basins_EAS
    elseif (region%name == 'GRL') then
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_GRL
      choice_basin_scheme       = C%choice_basin_scheme_GRL
      filename_basins           = C%filename_basins_GRL
    elseif (region%name == 'ANT') then
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_ANT
      choice_basin_scheme       = C%choice_basin_scheme_ANT
      filename_basins           = C%filename_basins_ANT
    end if

    ! Check if the NetCDF file has all the required dimensions and variables
    call inquire_hires_geometry_file( hires)

    ! Allocate shared memory for x,y and the actual data
    allocate( hires%grid%x (1:hires%grid%nx))
    allocate( hires%grid%y (1:hires%grid%ny))
    allocate( hires%Hi     (1:hires%grid%nx, 1:hires%grid%nx))
    allocate( hires%Hb     (1:hires%grid%nx, 1:hires%grid%nx))

    ! Read the data from the NetCDF file
    if (par%master) then
      write(*,"(3A)") '    Reading high-resolution geometry for ocean extrapolation from file "', &
                           trim( hires%netcdf_geo%filename), '"...'

    end if
    call read_hires_geometry_file( hires)

    ! Projection parameters are of course identical to those used for this ice model region
    hires%grid%lambda_M     = region%mesh%lambda_M
    hires%grid%phi_M        = region%mesh%phi_M
    hires%grid%alpha_stereo = region%mesh%alpha_stereo

    ! But the resolution is different
    hires%grid%dx           = hires%grid%x( 2) - hires%grid%x( 1)
    hires%grid%xmin         = hires%grid%x( 1            )
    hires%grid%xmax         = hires%grid%x( hires%grid%nx)
    hires%grid%ymin         = hires%grid%y( 1            )
    hires%grid%ymax         = hires%grid%y( hires%grid%ny)

    ! Check if this is the resolution we want
    if (hires%grid%dx /= C%ocean_extrap_res) then
      call crash('high-resolution geometry file "' // &
                  TRIM( hires%netcdf_geo%filename) // &
                  '" has a different resolution from C%ocean_extrap_res = {dp_01}', &
                  dp_01 = C%ocean_extrap_res)
    end if

    ! Tolerance; points lying within this distance of each other are treated as identical
    hires%grid%tol_dist = ((hires%grid%xmax - hires%grid%xmin) + (hires%grid%ymax - hires%grid%ymin)) * tol / 2._dp

    ! Set up grid-to-vector translation tables
    hires%grid%n = hires%grid%nx * hires%grid%ny

    allocate( hires%grid%ij2n (1:hires%grid%nx, 1:hires%grid%ny))
    allocate( hires%grid%n2ij (hires%grid%n, 2))

    n = 0
    do i = 1, hires%grid%nx
      if (mod(i,2) == 1) then
        do j = 1, hires%grid%ny
          n = n+1
          hires%grid%ij2n( i,j) = n
          hires%grid%n2ij( n,:) = [i,j]
        end do
      else
        do j = hires%grid%ny, 1, -1
          n = n+1
          hires%grid%ij2n( i,j) = n
          hires%grid%n2ij( n,:) = [i,j]
        end do
      end if
    end do

    ! Assign range to each processor
    call partition_list( hires%grid%nx, par%i, par%n, hires%grid%i1, hires%grid%i2)
    call partition_list( hires%grid%ny, par%i, par%n, hires%grid%j1, hires%grid%j2)

    ! Lat,lon coordinates
    allocate( hires%grid%lat (1:hires%grid%nx, 1:hires%grid%ny))
    allocate( hires%grid%lon (1:hires%grid%nx, 1:hires%grid%ny))

    do i = 1, hires%grid%nx
    do j = 1, hires%grid%ny
      call inverse_oblique_sg_projection( hires%grid%x( i), hires%grid%y( j), &
                                          hires%grid%lambda_M, hires%grid%phi_M, &
                                          hires%grid%alpha_stereo, &
                                          hires%grid%lon( i,j), hires%grid%lat( i,j))
    end do
    end do

    ! ===== Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid =====
    ! ================================================================================================

    if (par%master) then
      write(*,"(A,F4.1,A)") '     Mapping ocean data from the global lat/lon-grid to the ', &
                                  hires%grid%dx / 1000._dp, ' km regional x/y-grid...'
    end if

    ! Allocate shared memory for high-resolution ocean data
    allocate( hires%T_ocean (hires%grid%nx, hires%grid%ny, C%nz_ocean))
    allocate( hires%S_ocean (hires%grid%nx, hires%grid%ny, C%nz_ocean))

    ! Map the data from the global lon/lat-grid to the high-resolution regional x/y-grid
    call map_glob_to_grid_3D( ocean_glob%nlat, ocean_glob%nlon, &
                              ocean_glob%lat, ocean_glob%lon, &
                              hires%grid, ocean_glob%T_ocean, hires%T_ocean)
    call map_glob_to_grid_3D( ocean_glob%nlat, ocean_glob%nlon, &
                              ocean_glob%lat, ocean_glob%lon, &
                              hires%grid, ocean_glob%S_ocean, hires%S_ocean)

    ! ===== Perform the extrapolation on the high-resolution grid =====
    ! =================================================================

    IF (par%master) then
      write(*,'(A,F4.1,A)') '     Defining ice basins on the ', &
                                  hires%grid%dx / 1000._dp, &
                                  ' km regional x/y-grid...'
    end if

    ! Allocate shared memory for ice basins on the high-resolution grid
    allocate( hires%basin_ID (1:hires%grid%nx, 1:hires%grid%ny))

    hires%nbasins = region%ice%nbasins

    ! Instead of doing the "proper" basin definition on high resolution (which is insanely slow),
    ! just downscale the basin ID field from the ice model (using some tricks to get accurate values near the boundaries)

    ! Allocate shared memory
    allocate( basin_ID_dp_lores     (1:region%mesh%nV))
    allocate( basin_ID_dp_hires     (1:hires%grid%nx, 1:hires%grid%ny))
    allocate( basin_ID_dp_hires_ext (1:hires%grid%nx, 1:hires%grid%ny))

    select case (choice_basin_scheme)

      case ('none')
        ! No basins are defined (i.e. the whole region is one big, big basin)

        do j = 1, hires%grid%ny
        do i = hires%grid%i1, hires%grid%i2
          hires%basin_ID( i,j) = 1
        end do
        end do

      case ('file')
        ! Convert basin ID field to double precision (for remapping)
        do vi = region%mesh%vi1, region%mesh%vi2
          basin_ID_dp_lores( vi) = REAL( region%ice%basin_ID( vi), dp)
        end do

        ! Map double-precision basin ID from ice-model mesh to high-resolution grid
        call calc_remapping_operator_mesh2grid( region%mesh, hires%grid)
        call map_mesh2grid_2D( region%mesh, hires%grid, basin_ID_dp_lores, basin_ID_dp_hires)
        call deallocate_remapping_operators_mesh2grid( hires%grid)

        ! Remove all near-boundary cells
        do j = 1, hires%grid%ny
        do i = hires%grid%i1, hires%grid%i2
          if (modulo( basin_ID_dp_hires( i,j), 1._dp) > 0.01_dp) then
            basin_ID_dp_hires( i,j) = -1._dp
          end if
        end do
        end do

        ! For those, use extrapolation instead
        basin_ID_dp_hires_ext( hires%grid%i1:hires%grid%i2,:) = basin_ID_dp_hires( hires%grid%i1:hires%grid%i2,:)

        do j = 1, hires%grid%ny
        do i = hires%grid%i1, hires%grid%i2

          if (basin_ID_dp_hires_ext( i,j) == -1._dp) then

            n = 0
            foundit = .false.
            do while (.not. foundit)

              n = n+1

              ! Take the value of the nearest non-boundary cell
              do jj = MAX(1,j-n), MIN(hires%grid%ny,j+n)
              do ii = MAX(1,i-n), MIN(hires%grid%nx,i+n)
                if (basin_ID_dp_hires( ii,jj) > -1._dp) then
                  basin_ID_dp_hires_ext( i,j) = basin_ID_dp_hires( ii,jj)
                  foundit = .true.
                  exit
                end if
              end do
              if (foundit) then
                exit
              end if
              end do

              ! Safety
              if (n > MAX(hires%grid%nx, hires%grid%ny)) then
                call crash('map_and_extrapolate_ocean_data - ERROR: basin ID downscaling got stuck!')
              end if

            end do ! while (.not. foundit)

          end if ! (basin_ID_dp_hires_ext( i,j) == -1._dp)

        end do
        end do

        ! Convert hi-resolution basin ID field back to integer precision
        do j = 1, hires%grid%ny
        do i = hires%grid%i1, hires%grid%i2
          hires%basin_ID( i,j) = nint( basin_ID_dp_hires_ext( i,j))
        end do
        end do

      case default
        ! Unknown case
        call crash('unknown choice_basin_scheme "' // TRIM(choice_basin_scheme) // '"!')

    end select

    ! Clean up after yourself
    deallocate( basin_ID_dp_lores)
    deallocate( basin_ID_dp_hires)
    deallocate( basin_ID_dp_hires_ext)

    ! Perform the extrapolation on the high-resolution grid
    if (par%master) then
      write(*,'(A,F4.1,A)') '     Performing ocean data extrapolation on the ', &
                                  hires%grid%dx / 1000._dp, ' km regional x/y-grid...'
    end if

    call extend_regional_ocean_data_to_cover_domain( hires)

    ! Clean up fields that were needed only for the extrapolation
    deallocate( hires%Hi      )
    deallocate( hires%Hb      )
    deallocate( hires%basin_ID)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected=20)

  end subroutine map_and_extrapolate_hires_ocean_data

  subroutine extend_regional_ocean_data_to_cover_domain( hires)
    ! Extend global ocean data over the whole grid, based on the procedure outlined in
    ! Jourdain, N. C., Asay-Davis, X., Hattermann, T., Straneo, F., Seroussi, H., Little, C. M., & Nowicki, S. (2020).
    ! A protocol for calculating basal melt rates in the ISMIP6 Antarctic ice sheet projections. The Cryosphere, 14(9), 3111-3134.

    implicit none

    ! In/output variables:
    type(type_highres_ocean_data), intent(inout) :: hires

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'extend_regional_ocean_data_to_cover_domain'
    integer                                      :: i,j,k
    real(dp)                                     :: NaN, Hs, z_bedrock, z_icebase, z
    integer,  dimension(:,:,:), allocatable      :: mask_wetdry, mask_hasdata
    integer                                      :: k1,k2,bi
    integer,  dimension(:,:  ), allocatable      :: mask, mask_filled
    real(dp), dimension(:,:  ), allocatable      :: d_T, d_S
    real(dp), dimension(:,:,:), allocatable      :: T_ocean_ext,  S_ocean_ext
    logical,  parameter                          :: verbose = .FALSE.

    ! Add routine to path
    call init_routine( routine_name)

    ! Useful trick
    NaN = -1._dp
    NaN = sqrt( NaN)

    ! Allocate shared memory
    allocate( T_ocean_ext (hires%grid%nx, hires%grid%ny, C%nz_ocean))
    allocate( S_ocean_ext (hires%grid%nx, hires%grid%ny, C%nz_ocean))

    ! Initialise the extrapolated product with the provided ocean data
    do k = 1, C%nz_ocean
    do j = 1, hires%grid%ny
    do i = hires%grid%i1, hires%grid%i2
      T_ocean_ext( i,j,k) = hires%T_ocean( i,j,k)
      S_ocean_ext( i,j,k) = hires%S_ocean( i,j,k)
    end do
    end do
    end do

    ! Define the two masks needed for the four extrapolation steps:
    !
    !  - mask_wetdry:
    !      1 = actually    wet (i.e. open ocean, sub-shelf cavity, above sea floor and beneath ice base)
    !      2 = potentially wet (i.e. grounded marine ice, above sea floor)
    !      3 =             dry (i.e. beneath bedrock                     )
    !
    !  - mask_hasdata:
    !      0 = has no data
    !      1 = has data provided
    !      2 = has data extrapolated

    allocate( mask_wetdry  (hires%grid%nx, hires%grid%ny, C%nz_ocean))
    allocate( mask_hasdata (hires%grid%nx, hires%grid%ny, C%nz_ocean))

    do j = 1, hires%grid%ny
    do i = hires%grid%i1, hires%grid%i2

      Hs = surface_elevation( hires%Hi( i,j), hires%Hb( i,j), 0._dp)
      z_bedrock = hires%Hb( i,j)
      z_icebase = Hs - hires%Hi( i,j)

      do k = 1, C%nz_ocean

        z = -C%z_ocean( k)

        ! mask_wetdry
        if (z < z_bedrock) then
          ! This 3D hires%grid box is beneath the bedrock surface, so dry
          mask_wetdry( i,j,k) = 3
        else
          ! This 3D hires%grid box is above the bedrock surface so at least potentially wet
          if (z < z_icebase) then
            ! This 3D hires%grid box is above the bedrock surface and below the ice base, so it is actually wet
            mask_wetdry( i,j,k) = 1
          else
            ! This 3D hires%grid box is above the bedrock surface and above the ice base, so it is potentially wet (i.e. inside grounded marine ice)
            mask_wetdry( i,j,k) = 2
          end if
        end if

        ! mask_hasdata
        if (hires%T_ocean( i,j,k) /= hires%T_ocean( i,j,k)) then
          ! This 3D hires%grid box has no data (yet)
          mask_hasdata( i,j,k) = 0
        else
          ! Data is already provided for this 3D hires%grid box
          mask_hasdata( i,j,k) = 1
        end if

      end do

    end do
    end do

    ! ================================================================
    ! ===== Step 1: horizontal extrapolation into shelf cavities =====
    ! ================================================================

    if (par%master .and. verbose) then
      write(*,"(A)") '    extend_regional_ocean_data_to_cover_domain - step 1'
    end if

    ! Here, we start with the ocean data as provided (i.e. only for open ocean), and
    ! perform a horizontal extrapolation (so for each vertical layer separately) into
    ! the shelf cavities. Only "actually wet" 3D grid boxes are allowed to be filled,
    ! so the fill is limited by both bedrock sills and grounded ice.

    ! Allocate memory for the mask and data field of a single extrapolation step
    allocate( mask(        hires%grid%nx, hires%grid%ny))
    allocate( mask_filled( hires%grid%nx, hires%grid%ny))
    allocate( d_T(         hires%grid%nx, hires%grid%ny))
    allocate( d_S(         hires%grid%nx, hires%grid%ny))

    ! Parallelised by partitioning the vertical domain
    call partition_list( C%nz_ocean, par%i, par%n, k1, k2)

    do k = k1, k2

      ! Extrapolate per basin
      do bi = 1, hires%nbasins

        if (verbose) then
          write(*,'(A,I2,A,I3,A,I3,A,I3,A,I3)') '        process ', par%i, &
                                                ': vertical layer ', k, '/', C%nz_ocean, &
                                                ', basin ', bi, '/', hires%nbasins
        end if

        ! Define the mask and initial data fields for this particular flood-fill
        ! (i.e. this vertical layer and this basin)
        mask        = 0
        d_T         = NaN
        d_S         = NaN
        do j = 1, hires%grid%ny
        do i = 1, hires%grid%nx
          if (hires%basin_ID( i,j) == bi) then
            if (mask_hasdata( i,j,k) == 1) then
              ! This is where the source data comes from
              mask( i,j) = 2
              d_T(  i,j) = T_ocean_ext( i,j,k)
              d_S(  i,j) = S_ocean_ext( i,j,k)
            elseif (mask_hasdata( i,j,k) == 0 .AND. mask_wetdry( i,j,k) == 1) then
              ! This is where we're supposed to fill it in
              mask( i,j) = 1
            end if
          end if
        end do
        end do

        ! Perform the flood-fill-based Gaussian extrapolation
        call extrapolate_Gaussian_floodfill( hires%grid, mask, d_T, C%ocean_extrap_Gauss_sigma, mask_filled)
        call extrapolate_Gaussian_floodfill( hires%grid, mask, d_S, C%ocean_extrap_Gauss_sigma, mask_filled)

        ! Copy extrapolated data to the data structure

        do j = 1, hires%grid%ny
        do i = 1, hires%grid%nx
          if (mask_filled( i,j) == 1) then
            T_ocean_ext(  i,j,k) = d_T( i,j)
            S_ocean_ext(  i,j,k) = d_S( i,j)
            mask_hasdata( i,j,k) = 2
          end if
        end do
        end do

      end do ! DO bi = 1, ice%nbasins

    end do ! DO k = k1, k2

    ! Clean up after yourself
    deallocate( mask       )
    deallocate( mask_filled)
    deallocate( d_T        )
    deallocate( d_S        )

    ! ===========================================================================
    ! ===== Step 2: vertical extrapolation into sill-blocked shelf cavities =====
    ! ===========================================================================

    if (par%master .and. verbose) then
      write(*,"(A)") '    extend_regional_ocean_data_to_cover_domain - step 2'
    end if

    ! Here, we start with the ocean data that has been horizontally extrapolated into
    ! the shelf cavities, allowing for bedrock topography to block the fill, so that
    ! for example the lower parts of the Filchner-Ronne and Ross cavities have not yet
    ! been filled. We now extrapolate the data vertically from the filled parts to
    ! fill those parts of the cavities. Barring any really weird geometry, the entire
    ! cavities will now be filled.

    do j = 1, hires%grid%ny
    do i = hires%grid%i1, hires%grid%i2

      ! Move down through the vertical column
      do k = 2, C%nz_ocean
        ! If this grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        if (mask_wetdry( i,j,k) == 1 .and. mask_hasdata( i,j,k) == 0) then
          ! This 3D grid box is wet but has no data
          if (mask_hasdata( i,j,k-1) == 1 .or. mask_hasdata( i,j,k-1) == 2) then
            ! The one above it has data; copy data
            mask_hasdata( i,j,k) = 2
            T_ocean_ext( i,j,k) = T_ocean_ext( i,j,k-1)
            S_ocean_ext( i,j,k) = S_ocean_ext( i,j,k-1)
          end if
        end if
      end do

    end do
    end do

    ! ===============================================================
    ! ===== Step 3: vertical extrapolation into ice and bedrock =====
    ! ===============================================================

    if (par%master .and. verbose) then
      write(*,"(A)") '    extend_regional_ocean_data_to_cover_domain - step 3'
    end if

    ! Extrapolate data vertically into 3D grid boxes that are occupied by ice
    ! or bedrock (since they might turn into ocean at some point during a simulation)

    do j = 1, hires%grid%ny
    do i = hires%grid%i1, hires%grid%i2

      ! Move down through the vertical column
      do k = 2, C%nz_ocean
        ! If this grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        if (mask_hasdata( i,j,k) == 0) then
          ! This 3D grid box is wet but has no data
          if (mask_hasdata( i,j,k-1) == 1 .or. mask_hasdata( i,j,k-1) == 2) then
            ! The one above it has data; copy data
            mask_hasdata( i,j,k) = 2
            T_ocean_ext( i,j,k) = T_ocean_ext( i,j,k-1)
            S_ocean_ext( i,j,k) = S_ocean_ext( i,j,k-1)
          end if
        end if
      end do

      ! Move up through the vertical column
      do k = C%nz_ocean-1, 1, -1
        ! If this grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        if (mask_hasdata( i,j,k) == 0) then
          ! This 3D grid box is wet but has no data
          if (mask_hasdata( i,j,k+1) == 1 .or. mask_hasdata( i,j,k+1) == 2) then
            ! The one above it has data; copy data
            mask_hasdata( i,j,k) = 2
            T_ocean_ext(  i,j,k) = T_ocean_ext( i,j,k+1)
            S_ocean_ext(  i,j,k) = S_ocean_ext( i,j,k+1)
          end if
        end if
      end do

    end do
    end do

    ! =================================================================
    ! ===== Step 4: horizontal extrapolation into ice and bedrock =====
    ! =================================================================

    if (par%master .and. verbose) then
      write(*,"(A)") '    extend_regional_ocean_data_to_cover_domain - step 4'
    end if

    ! In the last step, extrapolate data horizontally into 3D
    ! grid boxes that are occupied by ice or bedrock

    ! Allocate memory for the mask and data field of a single extrapolation step
    allocate( mask(        hires%grid%nx, hires%grid%ny))
    allocate( mask_filled( hires%grid%nx, hires%grid%ny))
    allocate( d_T(         hires%grid%nx, hires%grid%ny))
    allocate( d_S(         hires%grid%nx, hires%grid%ny))

    ! Parallelised by partitioning the vertical domain
    call partition_list( C%nz_ocean, par%i, par%n, k1, k2)

    do k = k1, k2

      ! Extrapolate per basin
      do bi = 1, hires%nbasins

        if (verbose) then
          write(*,'(A,I2,A,I3,A,I3,A,I3,A,I3)') '        process ', par%i, &
                                                ': vertical layer ', k, &
                                                '/', C%nz_ocean, &
                                                ', basin ', bi, &
                                                '/', hires%nbasins
        end if

        ! Define the mask and initial data fields for this particular flood-fill
        ! (i.e. this vertical layer and this basin)
        mask        = 0
        d_T         = NaN
        d_S         = NaN
        do j = 1, hires%grid%ny
        do i = 1, hires%grid%nx
          if (hires%basin_ID( i,j) == bi) then
            if (mask_hasdata( i,j,k) == 1 .or. mask_hasdata( i,j,k) == 2) then
              ! This is where the source data comes from
              mask( i,j) = 2
              d_T(  i,j) = T_ocean_ext( i,j,k)
              d_S(  i,j) = S_ocean_ext( i,j,k)
            elseif (mask_hasdata( i,j,k) == 0) then
              ! This is where we're supposed to fill it in
              mask( i,j) = 1
            end if
          end if
        end do
        end do

        ! Perform the flood-fill-based Gaussian extrapolation

        call extrapolate_Gaussian_floodfill( hires%grid, mask, d_T, C%ocean_extrap_Gauss_sigma, mask_filled)
        call extrapolate_Gaussian_floodfill( hires%grid, mask, d_S, C%ocean_extrap_Gauss_sigma, mask_filled)

        ! Copy extrapolated data to the data structure
        do j = 1, hires%grid%ny
        do i = 1, hires%grid%nx
          if (mask_filled( i,j) == 1) then
            T_ocean_ext(  i,j,k) = d_T( i,j)
            S_ocean_ext(  i,j,k) = d_S( i,j)
            mask_hasdata( i,j,k) = 2
          end if
        end do
        end do

      end do ! bi = 1, ice%nbasins

      ! One more pass without considering basins (sometimes, a handful of isolated "basin enclaves" can
      ! occur at high resolution, which will not be filled when using the basin-constrained flood-fill)

      ! Define the mask and initial data fields for this particular flood-fill
      ! (i.e. this vertical layer and this basin)
      mask        = 0
      d_T         = NaN
      d_S         = NaN

      do j = 1, hires%grid%ny
      do i = 1, hires%grid%nx
        if (mask_hasdata( i,j,k) == 1 .or. mask_hasdata( i,j,k) == 2) then
          ! This is where the source data comes from
          mask( i,j) = 2
          d_T(  i,j) = T_ocean_ext( i,j,k)
          d_S(  i,j) = S_ocean_ext( i,j,k)
        elseif (mask_hasdata( i,j,k) == 0) then
          ! This is where we're supposed to fill it in
          mask( i,j) = 1
        end if
      end do
      end do

      ! Perform the flood-fill-based Gaussian extrapolation
      call extrapolate_Gaussian_floodfill( hires%grid, mask, d_T, C%ocean_extrap_Gauss_sigma, mask_filled)
      call extrapolate_Gaussian_floodfill( hires%grid, mask, d_S, C%ocean_extrap_Gauss_sigma, mask_filled)

      ! Copy extrapolated data to the data structure
      do j = 1, hires%grid%ny
      do i = 1, hires%grid%nx
        if (mask_filled( i,j) == 1) then
          T_ocean_ext(  i,j,k) = d_T( i,j)
          S_ocean_ext(  i,j,k) = d_S( i,j)
          mask_hasdata( i,j,k) = 2
        end if
      end do
      end do

      ! Check if all pixels have now been filled
      do j = 1, hires%grid%ny
      do i = 1, hires%grid%nx
        if (T_ocean_ext( i,j,k) /= T_ocean_ext( i,j,k)) then
          write(0,*) '    extend_regional_ocean_data_to_cover_domain - ERROR: unfilled pixels remains at [i,j] = [', i,',',j,']'
          call MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        end if
      end do
      end do

    end do ! k = k1, k2

    ! Copy data to the hires structure
    do k = 1, C%nz_ocean
    do j = 1, hires%grid%ny
    do i = hires%grid%i1, hires%grid%i2
      hires%T_ocean( i,j,k) = T_ocean_ext( i,j,k)
      hires%S_ocean( i,j,k) = S_ocean_ext( i,j,k)
    end do
    end do
    end do

    ! Clean up after yourself
    deallocate( mask        )
    deallocate( mask_filled )
    deallocate( d_T         )
    deallocate( d_S         )
    deallocate( mask_wetdry )
    deallocate( mask_hasdata)
    deallocate( T_ocean_ext )
    deallocate( S_ocean_ext )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extend_regional_ocean_data_to_cover_domain

! ===== Administration =====
! ==========================

  subroutine allocate_ocean_snapshot_regional( mesh, ocean, name)
    ! Allocate shared memory for a regional ocean snapshot

    implicit none

    ! In/output variables:
    type(type_mesh),                    intent(in)    :: mesh
    type(type_ocean_snapshot_regional), intent(inout) :: ocean
    character(len=3),                   intent(in)    :: name

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'allocate_ocean_snapshot_regional'

    ! Add routine to path
    call init_routine( routine_name)

    ocean%name = name

    ! Allocate shared memory
    allocate( ocean%T_ocean          (mesh%vi1:mesh%vi2, C%nz_ocean) )
    allocate( ocean%S_ocean          (mesh%vi1:mesh%vi2, C%nz_ocean) )
    allocate( ocean%T_ocean_ext      (mesh%vi1:mesh%vi2, C%nz_ocean) )
    allocate( ocean%S_ocean_ext      (mesh%vi1:mesh%vi2, C%nz_ocean) )
    allocate( ocean%T_ocean_corr_ext (mesh%vi1:mesh%vi2, C%nz_ocean) )
    allocate( ocean%S_ocean_corr_ext (mesh%vi1:mesh%vi2, C%nz_ocean) )

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected=6)

  end subroutine allocate_ocean_snapshot_regional

end module