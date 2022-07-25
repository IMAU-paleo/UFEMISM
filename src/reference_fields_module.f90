MODULE reference_fields_module
  ! Contains the routines for setting up the three "reference geometries":
  ! - refgeo_PD:     present-day, used to calculate sea-level contribution, isotope change, and more
  ! - refgeo_init:   initial, used to initialise the simulation
  ! - refgeo_GIA_eq: GIA equilibrium, used for the GIA model

#include <petsc/finclude/petscksp.h>

! ===== USEs =====
! ================

  use mpi
  use petscksp
  use configuration_module,            only: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module
  use petsc_module,                    only: perr
  use parallel_module,                 only: par, sync, ierr, cerr, partition_list
  use data_types_module,               only: type_reference_geometry, type_grid, &
                                             type_model_region, type_mesh
  use netcdf_module,                   only: inquire_reference_geometry_file, read_reference_geometry_file
  use mesh_mapping_module,             only: calc_remapping_operator_grid2mesh, map_grid2mesh_2D, map_grid2mesh_2D_partial, deallocate_remapping_operators_grid2mesh
  use utilities_module,                only: check_for_NaN_dp_2D, remove_Lake_Vostok, &
                                             is_floating, surface_elevation
  ! use netcdf_module,                   only: debug, write_to_debug_file

! ===== Preamble =====
! ====================

  implicit none

contains

  ! Initialise all three reference geometries
  subroutine initialise_reference_geometries( refgeo_init, refgeo_PD, refgeo_GIAeq, region_name)
    ! Initialise all three reference geometries

    implicit none

    ! In/output variables:
    type(type_reference_geometry),  intent(inout) :: refgeo_init
    type(type_reference_geometry),  intent(inout) :: refgeo_PD
    type(type_reference_geometry),  intent(inout) :: refgeo_GIAeq
    character(len=3),               intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'initialise_reference_geometries'
    character(len=256)                            :: choice_refgeo_init, choice_refgeo_PD, choice_refgeo_GIAeq
    character(len=256)                            :: filename_refgeo_init, filename_refgeo_PD, filename_refgeo_GIAeq
    real(dp)                                      :: time_to_restart_from

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine parameters for this region
    if     (region_name == 'NAM') then
      choice_refgeo_init    = C%choice_refgeo_init_NAM
      choice_refgeo_PD      = C%choice_refgeo_PD_NAM
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_NAM
      filename_refgeo_init  = C%filename_refgeo_init_NAM
      filename_refgeo_PD    = C%filename_refgeo_PD_NAM
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_NAM
      time_to_restart_from  = C%time_to_restart_from_NAM
    elseif (region_name == 'EAS') then
      choice_refgeo_init    = C%choice_refgeo_init_EAS
      choice_refgeo_PD      = C%choice_refgeo_PD_EAS
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_EAS
      filename_refgeo_init  = C%filename_refgeo_init_EAS
      filename_refgeo_PD    = C%filename_refgeo_PD_EAS
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_EAS
      time_to_restart_from  = C%time_to_restart_from_EAS
    elseif (region_name == 'GRL') then
      choice_refgeo_init    = C%choice_refgeo_init_GRL
      choice_refgeo_PD      = C%choice_refgeo_PD_GRL
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_GRL
      filename_refgeo_init  = C%filename_refgeo_init_GRL
      filename_refgeo_PD    = C%filename_refgeo_PD_GRL
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_GRL
      time_to_restart_from  = C%time_to_restart_from_GRL
    elseif (region_name == 'ANT') then
      choice_refgeo_init    = C%choice_refgeo_init_ANT
      choice_refgeo_PD      = C%choice_refgeo_PD_ANT
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_ANT
      filename_refgeo_init  = C%filename_refgeo_init_ANT
      filename_refgeo_PD    = C%filename_refgeo_PD_ANT
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_ANT
      time_to_restart_from  = C%time_to_restart_from_ANT
    end if

    ! Initial ice-sheet geometry
    ! ==========================

    ! if     (choice_refgeo_init == 'idealised') then
    !   if (par%master) WRITE(0,*) '  Initialising initial         reference geometry from idealised case "', TRIM(C%choice_refgeo_init_idealised), '"...'
    !   call initialise_reference_geometry_idealised_grid( refgeo_init, C%choice_refgeo_init_idealised, region_name, C%dx_refgeo_init_idealised)
    ! elseif (choice_refgeo_init == 'realistic') then
      if (par%master) then
        write(*,"(3A)") '  Initialising initial         reference geometry from file ', TRIM( filename_refgeo_init), '...'
      end if
      call initialise_reference_geometry_from_file( refgeo_init, filename_refgeo_init, region_name)
    ! else
    !   call crash('unknown choice_refgeo_init "' // TRIM( choice_refgeo_init) // '"!')
    ! end if

    ! Present-day ice-sheet geometry
    ! ==============================

    ! if     (choice_refgeo_PD == 'idealised') then
    !   if (par%master) WRITE(0,*) '  Initialising present-day     reference geometry from idealised case "', TRIM(C%choice_refgeo_PD_idealised), '"...'
    !   call initialise_reference_geometry_idealised_grid( refgeo_PD, C%choice_refgeo_PD_idealised, region_name, C%dx_refgeo_PD_idealised)
    ! elseif (choice_refgeo_PD == 'realistic') then
      if (par%master) then
        write(*,"(3A)") '  Initialising present-day     reference geometry from file ', TRIM( filename_refgeo_PD), '...'
      end if
      call initialise_reference_geometry_from_file( refgeo_PD, filename_refgeo_PD, region_name)
    ! else
    !   call crash('unknown choice_refgeo_PD "' // TRIM( choice_refgeo_PD) // '"!')
    ! end if

    ! ! GIA equilibrium ice-sheet geometry
    ! ! ==================================

    ! if     (choice_refgeo_GIAeq == 'idealised') then
    !   if (par%master) WRITE(0,*) '  Initialising GIA equilibrium reference geometry from idealised case "', TRIM(C%choice_refgeo_GIAeq_idealised), '"...'
    !   call initialise_reference_geometry_idealised_grid( refgeo_GIAeq, C%choice_refgeo_GIAeq_idealised, region_name, C%dx_refgeo_GIAeq_idealised)
    ! elseif (choice_refgeo_GIAeq == 'realistic') then
      if (par%master) then
        write(*,"(3A)") '  Initialising GIA equilibrium reference geometry from file ', TRIM( filename_refgeo_GIAeq), '...'
      end if
      call initialise_reference_geometry_from_file( refgeo_GIAeq, filename_refgeo_GIAeq, region_name)
    ! else
    !   call crash('unknown choice_refgeo_GIAeq "' // TRIM( choice_refgeo_GIAeq) // '"!')
    ! end if

    ! Fill in secondary data for the reference geometry (used to force mesh creation)
    call calc_reference_geometry_secondary_data( refgeo_init%grid , refgeo_init )
    call calc_reference_geometry_secondary_data( refgeo_PD%grid   , refgeo_PD   )
    call calc_reference_geometry_secondary_data( refgeo_GIAeq%grid, refgeo_GIAeq)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = 78)

  end subroutine initialise_reference_geometries

  ! Initialise a reference geometry with data from a (timeless) NetCDF file (e.g. BedMachine)
  subroutine initialise_reference_geometry_from_file( refgeo, filename_refgeo, region_name)
    ! Initialise a reference geometry with data from a NetCDF file

    implicit none

    ! In/output variables:
    type(type_reference_geometry),  intent(inout) :: refgeo
    character(len=256),             intent(in)    :: filename_refgeo
    character(len=3),               intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_from_file'
    integer                                       :: i,j,n
    real(dp), parameter                           :: tol = 1E-9_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire if all the required fields are present in the specified NetCDF file,
    ! and determine the dimensions of the memory to be allocated.
    refgeo%netcdf%filename = filename_refgeo
    call inquire_reference_geometry_file( refgeo)

    ! Assign range to each processor
    call partition_list( refgeo%grid%nx, par%i, par%n, refgeo%grid%i1, refgeo%grid%i2)
    call partition_list( refgeo%grid%ny, par%i, par%n, refgeo%grid%j1, refgeo%grid%j2)

    ! Allocate memory for raw data
    allocate( refgeo%grid%x (  refgeo%grid%nx ))
    allocate( refgeo%grid%y (  refgeo%grid%ny ))

    allocate( refgeo%Hi_grid( refgeo%grid%nx, refgeo%grid%ny ))
    allocate( refgeo%Hb_grid( refgeo%grid%nx, refgeo%grid%ny ))
    allocate( refgeo%Hs_grid( refgeo%grid%nx, refgeo%grid%ny ))

    ! Read data from input file
    call read_reference_geometry_file( refgeo)

    ! Fill in secondary grid parameters
    refgeo%grid%dx = refgeo%grid%x( 2) - refgeo%grid%x( 1)
    refgeo%grid%xmin = refgeo%grid%x( 1             )
    refgeo%grid%xmax = refgeo%grid%x( refgeo%grid%nx)
    refgeo%grid%ymin = refgeo%grid%y( 1             )
    refgeo%grid%ymax = refgeo%grid%y( refgeo%grid%ny)

    ! Tolerance; points lying within this distance of each other are treated as identical
    refgeo%grid%tol_dist = ((refgeo%grid%xmax - refgeo%grid%xmin) + (refgeo%grid%ymax - refgeo%grid%ymin)) * tol / 2._dp

    ! Set up grid-to-vector translation tables
    refgeo%grid%n  = refgeo%grid%nx * refgeo%grid%ny
    allocate( refgeo%grid%ij2n (refgeo%grid%nx, refgeo%grid%ny))
    allocate( refgeo%grid%n2ij (refgeo%grid%n , 2 ))
    n = 0
    do i = 1, refgeo%grid%nx
      if (MOD(i,2) == 1) then
        do j = 1, refgeo%grid%ny
          n = n+1
          refgeo%grid%ij2n( i,j) = n
          refgeo%grid%n2ij( n,:) = [i,j]
        end do
      else
        do j = refgeo%grid%ny, 1, -1
          n = n+1
          refgeo%grid%ij2n( i,j) = n
          refgeo%grid%n2ij( n,:) = [i,j]
        end do
      end if
    end do

    ! Safety
    call check_for_NaN_dp_2D( refgeo%Hi_grid, 'refgeo%Hi_grid')
    call check_for_NaN_dp_2D( refgeo%Hb_grid, 'refgeo%Hb_grid')
    call check_for_NaN_dp_2D( refgeo%Hs_grid, 'refgeo%Hs_grid')

    ! Remove Lake Vostok from Antarctica (because it's annoying)
    if (region_name == 'ANT'.AND. C%remove_Lake_Vostok) then
      call remove_Lake_Vostok( refgeo%grid%x, refgeo%grid%y, refgeo%Hi_grid, refgeo%Hb_grid, refgeo%Hs_grid)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = 16)

  end subroutine initialise_reference_geometry_from_file

!   ! Initialise a reference geometry according to an idealised world on the initial grid
!   subroutine initialise_reference_geometry_idealised_grid( refgeo, choice_refgeo_idealised, region_name, dx)
!     ! Initialise a reference geometry according to an idealised world

!     implicit none

!     ! In/output variables:
!     type(type_reference_geometry),  intent(inout) :: refgeo
!     character(len=256),             intent(in)    :: choice_refgeo_idealised
!     character(len=3),               intent(in)    :: region_name
!     real(dp),                       intent(in)    :: dx

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid'

!     ! Add routine to path
!     call init_routine( routine_name)

!     ! Set up a grid
!     call setup_idealised_geometry_grid( refgeo%grid, region_name, dx)

!     ! Allocate shared memory
!     call allocate_shared_dp_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%Hi_grid, refgeo%wHi_grid)
!     call allocate_shared_dp_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%Hb_grid, refgeo%wHb_grid)
!     call allocate_shared_dp_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%Hs_grid, refgeo%wHs_grid)

!     if     (choice_refgeo_idealised == 'flatearth') then
!       ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments
!       call initialise_reference_geometry_idealised_grid_flatearth(     refgeo%grid, refgeo)
!     elseif (choice_refgeo_idealised == 'Halfar') then
!       ! The Halfar dome solution at t = 0
!       call initialise_reference_geometry_idealised_grid_Halfar(        refgeo%grid, refgeo)
!     elseif (choice_refgeo_idealised == 'Bueler') then
!       ! The Bueler dome solution at t = 0
!       call initialise_reference_geometry_idealised_grid_Bueler(        refgeo%grid, refgeo)
!     elseif (choice_refgeo_idealised == 'SSA_icestream') then
!       ! The SSA_icestream infinite slab on a flat slope
!       call initialise_reference_geometry_idealised_grid_SSA_icestream( refgeo%grid, refgeo)
!     elseif (choice_refgeo_idealised == 'MISMIP_mod') then
!       ! The MISMIP_mod cone-shaped island
!       call initialise_reference_geometry_idealised_grid_MISMIP_mod(    refgeo%grid, refgeo)
!     elseif (choice_refgeo_idealised == 'ISMIP_HOM_A') then
!       ! The ISMIP-HOM A bumpy slope
!       call initialise_reference_geometry_idealised_grid_ISMIP_HOM_A(   refgeo%grid, refgeo)
!     elseif (choice_refgeo_idealised == 'ISMIP_HOM_B') then
!       ! The ISMIP-HOM B bumpy slope
!       call initialise_reference_geometry_idealised_grid_ISMIP_HOM_B(   refgeo%grid, refgeo)
!     elseif (choice_refgeo_idealised == 'ISMIP_HOM_C' .OR. &
!             choice_refgeo_idealised == 'ISMIP_HOM_D') then
!       ! The ISMIP-HOM C/D bumpy slope
!       call initialise_reference_geometry_idealised_grid_ISMIP_HOM_CD(  refgeo%grid, refgeo)
!     elseif (choice_refgeo_idealised == 'ISMIP_HOM_E') then
!       ! The ISMIP-HOM E Glacier d'Arolla geometry
!       call initialise_reference_geometry_idealised_grid_ISMIP_HOM_E(   refgeo%grid, refgeo)
!     elseif (choice_refgeo_idealised == 'ISMIP_HOM_F') then
!       ! The ISMIP-HOM A bumpy slope
!       call initialise_reference_geometry_idealised_grid_ISMIP_HOM_F(   refgeo%grid, refgeo)
!     elseif (choice_refgeo_idealised == 'MISMIP+') then
!       ! The MISMIP+ fjord geometry
!       call initialise_reference_geometry_idealised_grid_MISMIPplus(    refgeo%grid, refgeo)
!     else
!       call crash('unknown choice_refgeo_idealised "' // TRIM( choice_refgeo_idealised) // '"!')
!     end if

!     ! Finalise routine path
!     call finalise_routine( routine_name, n_extra_windows_expected = 16)

!   end subroutine initialise_reference_geometry_idealised_grid
!   subroutine initialise_reference_geometry_idealised_grid_flatearth(     grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_flatearth'
!     integer                                       :: i,j

!     ! Add routine to path
!     call init_routine( routine_name)

!     do i = grid%i1, grid%i2
!     do j = 1, grid%ny
!       refgeo%Hi_grid( i,j) = 0._dp
!       refgeo%Hb_grid( i,j) = 0._dp
!       refgeo%Hs_grid( i,j) = 0._dp
!     end do
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_flatearth
!   subroutine initialise_reference_geometry_idealised_grid_Halfar(        grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The Halfar dome solution at t = 0

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_Halfar'
!     integer                                       :: i,j

!     ! Add routine to path
!     call init_routine( routine_name)

!     do i = grid%i1, grid%i2
!     do j = 1, grid%ny
!       refgeo%Hi_grid( i,j) = Halfar_solution( grid%x( i), grid%y( j), C%start_time_of_run)
!       refgeo%Hb_grid( i,j) = 0._dp
!       refgeo%Hs_grid( i,j) = refgeo%Hi_grid( i,j)
!     end do
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_Halfar
!   subroutine initialise_reference_geometry_idealised_grid_Bueler(        grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The Bueler dome solution at t = 0

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_Bueler'
!     integer                                       :: i,j

!     ! Add routine to path
!     call init_routine( routine_name)

!     do i = grid%i1, grid%i2
!     do j = 1, grid%ny
!       refgeo%Hi_grid( i,j) = Bueler_solution( grid%x( i), grid%y( j), C%start_time_of_run)
!       refgeo%Hb_grid( i,j) = 0._dp
!       refgeo%Hs_grid( i,j) = refgeo%Hi_grid( i,j)
!     end do
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_Bueler
!   subroutine initialise_reference_geometry_idealised_grid_SSA_icestream( grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The SSA_icestream infinite slab on a flat slope

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_SSA_icestream'
!     integer                                       :: i,j

!     ! Add routine to path
!     call init_routine( routine_name)

!     do i = grid%i1, grid%i2
!     do j = 1, grid%ny
!       refgeo%Hi_grid( i,j) = 2000._dp
!       refgeo%Hb_grid( i,j) = -0.001_dp * grid%x( i)
!       refgeo%Hs_grid( i,j) = refgeo%Hb_grid( i,j) + refgeo%Hi_grid( i,j)
!     end do
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_SSA_icestream
!   subroutine initialise_reference_geometry_idealised_grid_MISMIP_mod(    grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The MISMIP_mod cone-shaped island

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_MISMIP_mod'
!     integer                                       :: i,j

!     ! Add routine to path
!     call init_routine( routine_name)

!     do i = grid%i1, grid%i2
!     do j = 1, grid%ny

!       ! Create a nice circular ice shelf
!       if (SQRT(grid%x( i)**2 + grid%y( j)**2) < grid%xmax * 0.95_dp) then
!         refgeo%Hi_grid( i,j) = 100._dp
!       else
!         refgeo%Hi_grid( i,j) = 0._dp
!       end if

!       refgeo%Hb_grid( i,j) = 720._dp - 778.5_dp * SQRT( grid%x(i)**2 + grid%y(j)**2)/ 750000._dp
!       refgeo%Hs_grid( i,j) = surface_elevation( refgeo%Hi_grid( i,j), refgeo%Hb_grid( i,j), 0._dp)
!     end do
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_MISMIP_mod
!   subroutine initialise_reference_geometry_idealised_grid_ISMIP_HOM_A(   grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The ISMIP-HOM A bumpy slope

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_ISMIP_HOM_A'
!     integer                                       :: i,j

!     ! Add routine to path
!     call init_routine( routine_name)

!     do i = grid%i1, grid%i2
!     do j = 1, grid%ny
!       refgeo%Hs_grid( i,j) = 2000._dp - grid%x( i) * TAN( 0.5_dp * pi / 180._dp)
!       refgeo%Hb_grid( i,j) = refgeo%Hs_grid( i,j) - 1000._dp + 500._dp * SIN( grid%x( i) * 2._dp * pi / C%ISMIP_HOM_L) * SIN( grid%y( j) * 2._dp * pi / C%ISMIP_HOM_L)
!       refgeo%Hi_grid( i,j) = refgeo%Hs_grid( i,j) - refgeo%Hb_grid( i,j)
!     end do
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_ISMIP_HOM_A
!   subroutine initialise_reference_geometry_idealised_grid_ISMIP_HOM_B(   grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The ISMIP-HOM B bumpy slope

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_ISMIP_HOM_B'
!     integer                                       :: i,j

!     ! Add routine to path
!     call init_routine( routine_name)

!     do i = grid%i1, grid%i2
!     do j = 1, grid%ny
!       refgeo%Hs_grid( i,j) = 2000._dp - grid%x( i) * TAN( 0.5_dp * pi / 180._dp)
!       refgeo%Hb_grid( i,j) = refgeo%Hs_grid( i,j) - 1000._dp + 500._dp * SIN( grid%x(i) * 2._dp * pi / C%ISMIP_HOM_L)
!       refgeo%Hi_grid( i,j) = refgeo%Hs_grid( i,j) - refgeo%Hb_grid( i,j)
!     end do
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_ISMIP_HOM_B
!   subroutine initialise_reference_geometry_idealised_grid_ISMIP_HOM_CD(  grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The ISMIP-HOM C/D bumpy slope

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_ISMIP_HOM_CD'
!     integer                                       :: i,j

!     ! Add routine to path
!     call init_routine( routine_name)

!     do i = grid%i1, grid%i2
!     do j = 1, grid%ny
!       refgeo%Hs_grid( i,j) = 2000._dp - grid%x( i) * TAN( 0.1_dp * pi / 180._dp)
!       refgeo%Hb_grid( i,j) = refgeo%Hs_grid( i,j) - 1000._dp
!       refgeo%Hi_grid( i,j) = refgeo%Hs_grid( i,j) - refgeo%Hb_grid( i,j)
!     end do
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_ISMIP_HOM_CD
!   subroutine initialise_reference_geometry_idealised_grid_ISMIP_HOM_E(   grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The ISMIP-HOM E Glacier d'Arolla geometry

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_ISMIP_HOM_E'
!     integer                                       :: i,j
!     real(dp)                                      :: x,Hs,Hb
!     integer                                       :: ios,slides

!     ! Add routine to path
!     call init_routine( routine_name)

!     ! To prevent compiler warnings from unused variables
!     i = grid%nx

!     ! Read data from external file
!     if (par%master) then

!       OPEN( UNIT = 1337, FILE=C%ISMIP_HOM_E_Arolla_filename, ACTION = 'READ')
!       do i = 1, 51
!         READ( UNIT = 1337, FMT=*, IOSTAT=ios) x, Hb, Hs, slides
!         do j = 1, refgeo%grid%ny
!           refgeo%Hb_grid( i,j) = Hb
!           refgeo%Hi_grid( i,j) = Hs - Hb
!           refgeo%Hs_grid( i,j) = Hs
!         end do
!         if (ios /= 0) then
!           call crash(' length of text file "' // TRIM( C%ISMIP_HOM_E_Arolla_filename) // '" should be 51 lines!')
!         end if
!       end do
!       CLOSE( UNIT  = 1337)

!     end if ! if (par%master) then
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_ISMIP_HOM_E
!   subroutine initialise_reference_geometry_idealised_grid_ISMIP_HOM_F(   grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The ISMIP-HOM A bumpy slope

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_ISMIP_HOM_F'
!     integer                                       :: i,j

!     real(dp), parameter                           :: H0    = 1000._dp
!     real(dp), parameter                           :: a0    = 100._dp
!     real(dp), parameter                           :: sigma = 10000._dp

!     ! Add routine to path
!     call init_routine( routine_name)

!     do i = grid%i1, grid%i2
!     do j = 1, grid%ny
!       refgeo%Hs_grid( i,j) = 5000._dp - grid%x( i) * TAN( 3._dp * pi / 180._dp)
!       refgeo%Hb_grid( i,j) = refgeo%Hs_grid( i,j) - H0 + a0 * EXP( -((grid%x( i) - 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                                        + a0 * EXP( -((grid%x( i) - 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                                        + a0 * EXP( -((grid%x( i) - 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                                        + a0 * EXP( -((grid%x( i) - 0._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                                        + a0 * EXP( -((grid%x( i) - 0._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                                        + a0 * EXP( -((grid%x( i) - 0._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                                        + a0 * EXP( -((grid%x( i) + 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                                        + a0 * EXP( -((grid%x( i) + 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                                        + a0 * EXP( -((grid%x( i) + 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2)
!       refgeo%Hi_grid( i,j) = refgeo%Hs_grid( i,j) - refgeo%Hb_grid( i,j)
!     end do
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_ISMIP_HOM_F
!   subroutine initialise_reference_geometry_idealised_grid_MISMIPplus(    grid, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The MISMIpplus fjord geometry

!     implicit none

!     ! In/output variables:
!     type(type_grid),                intent(in)    :: grid
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_grid_MISMIPplus'
!     integer                                       :: i,j

!     real(dp)                                      :: x,y,xtilde,Bx,By
!     real(dp), parameter                           :: B0     = -150._dp
!     real(dp), parameter                           :: B2     = -728.8_dp
!     real(dp), parameter                           :: B4     = 343.91_dp
!     real(dp), parameter                           :: B6     = -50.57_dp
!     real(dp), parameter                           :: xbar   = 300000._dp
!     real(dp), parameter                           :: fc     = 4000._dp
!     real(dp), parameter                           :: dc     = 500._dp
!     real(dp), parameter                           :: wc     = 24000._dp
!     real(dp), parameter                           :: zbdeep = -720._dp

!     ! Add routine to path
!     call init_routine( routine_name)

!     do i = grid%i1, grid%i2
!     do j = 1, grid%ny
!       x = grid%x( i) + 400000._dp
!       y = -40000._dp +  80000._dp * real(j-1,dp) / real(grid%ny-1,dp)
!       xtilde = x / xbar
!       Bx = B0 + (B2 * xtilde**2._dp) + (B4 * xtilde**4._dp) + (B6 * xtilde**6._dp)
!       By = (dc / (1 + EXP(-2._dp*(y - wc)/fc))) + &
!            (dc / (1 + EXP( 2._dp*(y + wc)/fc)))
!       refgeo%Hi_grid( i,j) = 100._dp
!       refgeo%Hb_grid( i,j) = MAX( Bx + By, zbdeep)
!       refgeo%Hs_grid( i,j) = surface_elevation( refgeo%Hi_grid( i,j), refgeo%Hb_grid( i,j), 0._dp)
!     end do
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_grid_MISMIPplus

!   ! Set up a square grid for the idealised reference geometry (since it is now not provided externally)
!   subroutine setup_idealised_geometry_grid( grid, region_name, dx)
!     ! Set up a square grid for the idealised reference geometry (since it is now not provided externally)

!     implicit none

!     ! In/output variables:
!     type(type_grid),                 intent(inout)     :: grid
!     character(len=3),                intent(in)        :: region_name
!     real(dp),                        intent(in)        :: dx

!     ! Local variables:
!     character(len=256), parameter                     :: routine_name = 'setup_idealised_geometry_grid'
!     integer                                            :: nsx, nsy, i, j, n
!     real(dp)                                           :: xmin, xmax, ymin, ymax, xmid, ymid
!     real(dp), parameter                                :: tol = 1E-9_dp

!     ! Add routine to path
!     call init_routine( routine_name)

!     ! Assign dummy values to suppress compiler warnings
!     xmin = 0._dp
!     xmax = 0._dp
!     ymin = 0._dp
!     ymax = 0._dp
!     xmid = 0._dp
!     ymid = 0._dp
!     nsx  = 0
!     nsy  = 0

!     ! Allocate shared memory
!     call allocate_shared_int_0D(                   grid%nx          , grid%wnx          )
!     call allocate_shared_int_0D(                   grid%ny          , grid%wny          )
!     call allocate_shared_dp_0D(                    grid%dx          , grid%wdx          )
!     call allocate_shared_dp_0D(                    grid%xmin        , grid%wxmin        )
!     call allocate_shared_dp_0D(                    grid%xmax        , grid%wxmax        )
!     call allocate_shared_dp_0D(                    grid%ymin        , grid%wymin        )
!     call allocate_shared_dp_0D(                    grid%ymax        , grid%wymax        )

!     ! Let the Master do the work
!     if (par%master) then

!       grid%dx = dx

!       ! Resolution, domain size, and projection parameters for this region are determined from the config
!       if     (region_name == 'NAM') then
!         xmin = C%xmin_NAM
!         xmax = C%xmax_NAM
!         ymin = C%ymin_NAM
!         ymax = C%ymax_NAM
!       elseif (region_name == 'EAS') then
!         xmin = C%xmin_EAS
!         xmax = C%xmax_EAS
!         ymin = C%ymin_EAS
!         ymax = C%ymax_EAS
!       elseif (region_name == 'GRL') then
!         xmin = C%xmin_GRL
!         xmax = C%xmax_GRL
!         ymin = C%ymin_GRL
!         ymax = C%ymax_GRL
!       elseif (region_name == 'ANT') then
!         xmin = C%xmin_ANT
!         xmax = C%xmax_ANT
!         ymin = C%ymin_ANT
!         ymax = C%ymax_ANT
!       end if

!       ! Determine the number of grid cells we can fit in this domain
!       xmid = (xmax + xmin) / 2._dp
!       ymid = (ymax + ymin) / 2._dp
!       nsx  = FLOOR( (xmax - xmid) / grid%dx)
!       nsy  = FLOOR( (ymax - ymid) / grid%dx)

!       ! Small exceptions for very weird benchmark experiments
!       if (C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'SSA_icestream') nsx = 3
!       if (C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'ISMIP_HOM_E')   nsx = 25

!       grid%nx = 1 + 2*nsx
!       grid%ny = 1 + 2*nsy

!     end if ! if (par%master) then
!     call sync

!     ! Assign range to each processor
!     call partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
!     call partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

!     ! Allocate shared memory for x and y
!     call allocate_shared_dp_1D( grid%nx, grid%x, grid%wx)
!     call allocate_shared_dp_1D( grid%ny, grid%y, grid%wy)

!     ! Fill in x and y
!     if (par%master) then
!       do i = 1, grid%nx
!         grid%x( i) = -nsx*grid%dx + (i-1)*grid%dx
!       end do
!       do j = 1, grid%ny
!         grid%y( j) = -nsy*grid%dx + (j-1)*grid%dx
!       end do

!       grid%xmin = MINVAL(grid%x)
!       grid%xmax = MAXVAL(grid%x)
!       grid%ymin = MINVAL(grid%y)
!       grid%ymax = MAXVAL(grid%y)
!     end if ! if (par%master) then
!     call sync

!     ! Tolerance; points lying within this distance of each other are treated as identical
!     call allocate_shared_dp_0D( grid%tol_dist, grid%wtol_dist)
!     if (par%master) grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

!     ! Set up grid-to-vector translation tables
!     call allocate_shared_int_0D(                   grid%n           , grid%wn           )
!     if (par%master) grid%n  = grid%nx * grid%ny
!     call sync
!     call allocate_shared_int_2D( grid%nx, grid%ny, grid%ij2n        , grid%wij2n        )
!     call allocate_shared_int_2D( grid%n , 2      , grid%n2ij        , grid%wn2ij        )
!     if (par%master) then
!       n = 0
!       do i = 1, grid%nx
!         if (MOD(i,2) == 1) then
!           do j = 1, grid%ny
!             n = n+1
!             grid%ij2n( i,j) = n
!             grid%n2ij( n,:) = [i,j]
!           end do
!         else
!           do j = grid%ny, 1, -1
!             n = n+1
!             grid%ij2n( i,j) = n
!             grid%n2ij( n,:) = [i,j]
!           end do
!         end if
!       end do
!     end if
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name, n_extra_windows_expected = 13)

!   end subroutine setup_idealised_geometry_grid

  ! Fill in secondary data for the reference geometry (used to force mesh creation)
  subroutine calc_reference_geometry_secondary_data( grid, refgeo)
    ! Fill in secondary data for the reference geometry (used to force mesh creation)

    implicit none

    ! In/output variables:
    type(type_grid),                intent(in)    :: grid
    type(type_reference_geometry),  intent(inout) :: refgeo

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'calc_reference_geometry_secondary_data'
    integer                                       :: i,j,ii,jj
    real(dp)                                      :: d2Hs_dx2, d2Hs_dxdy, d2Hs_dy2

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( refgeo%surf_curv  ( grid%nx, grid%ny), source=0.0_dp)
    allocate( refgeo%mask_land  ( grid%nx, grid%ny), source=0)
    allocate( refgeo%mask_ocean ( grid%nx, grid%ny), source=0)
    allocate( refgeo%mask_ice   ( grid%nx, grid%ny), source=0)
    allocate( refgeo%mask_sheet ( grid%nx, grid%ny), source=0)
    allocate( refgeo%mask_shelf ( grid%nx, grid%ny), source=0)
    allocate( refgeo%mask_margin( grid%nx, grid%ny), source=0)
    allocate( refgeo%mask_gl    ( grid%nx, grid%ny), source=0)
    allocate( refgeo%mask_cf    ( grid%nx, grid%ny), source=0)
    allocate( refgeo%mask_coast ( grid%nx, grid%ny), source=0)

    ! Calculate surface curvature
    do i = 1, grid%nx
    do j = 1, grid%ny

      if (i == 1 .OR. i == grid%nx .OR. j == 1 .OR. j == grid%ny) then
        d2Hs_dx2  = 0._dp
        d2Hs_dxdy = 0._dp
        d2Hs_dy2  = 0._dp
      else
        d2Hs_dx2  = (refgeo%Hs_grid( i+1,j) + refgeo%Hs_grid( i-1,j) - 2._dp * refgeo%Hs_grid( i,j)) / (grid%dx**2)
        d2Hs_dxdy = (refgeo%Hs_grid( i+1,j+1) + refgeo%Hs_grid( i-1,j-1) - refgeo%Hs_grid( i+1,j-1) - refgeo%Hs_grid( i-1,j+1)) / (4._dp * grid%dx * grid%dx)
        d2Hs_dy2  = (refgeo%Hs_grid( i,j+1) + refgeo%Hs_grid( i,j-1) - 2._dp * refgeo%Hs_grid( i,j)) / (grid%dx**2)
      end if

      refgeo%surf_curv( i,j) = MAX( -1E-6_dp, MIN( 1E-6_dp, SQRT( d2Hs_dx2**2 + d2Hs_dxdy**2 + d2Hs_dy2**2)))

    end do
    end do

    ! Fill in masks

    ! Land/ocean
    do i = 1, grid%nx
    do j = 1, grid%ny

      if (is_floating( refgeo%Hi_grid( i,j), refgeo%Hb_grid( i,j), 0._dp)) then
        refgeo%mask_land(  i,j) = 1
      else
        refgeo%mask_ocean( i,j) = 1
      end if

    end do
    end do

    ! Ice/sheet/shelf
    do i = 1, grid%nx
    do j = 1, grid%ny

      if (refgeo%Hi_grid( i,j) > 0._dp) then

        refgeo%mask_ice(  i,j) = 1

        if (refgeo%mask_land( i,j) == 1) then
          refgeo%mask_sheet( i,j) = 1
        else
          refgeo%mask_shelf( i,j) = 1
        end if

      end if

    end do
    end do

    ! Transitions (margin, grounding line, calving front, coastline)
    do i = 1, grid%nx
    do j = 1, grid%ny

      ! Ice next to non-ice equals ice margin
      if (refgeo%mask_ice( i,j) == 1) then
        do ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        do jj = MAX( 1, j-1), MIN( grid%ny, j+1)
          if (refgeo%mask_ice( ii,jj) == 0) then
            refgeo%mask_margin( i,j) = 1
          end if
        end do
        end do
      end if

      ! Sheet next to shelf equals grounding line
      if (refgeo%mask_sheet( i,j) == 1) then
        do ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        do jj = MAX( 1, j-1), MIN( grid%ny, j+1)
          if (refgeo%mask_shelf( ii,jj) == 1) then
            refgeo%mask_gl( i,j) = 1
          end if
        end do
        end do
      end if

      ! Ice next to open ocean equals calving front
      if (refgeo%mask_ice( i,j) == 1) then
        do ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        do jj = MAX( 1, j-1), MIN( grid%ny, j+1)
          if (refgeo%mask_ocean( ii,jj) == 1 .AND. refgeo%mask_ice( ii,jj) == 0) then
            refgeo%mask_cf( i,j) = 1
          end if
        end do
        end do
      end if

      ! Dry land next to open ocean equals coastline
      if (refgeo%mask_land( i,j) == 1) then
        do ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        do jj = MAX( 1, j-1), MIN( grid%ny, j+1)
          if (refgeo%mask_ocean( ii,jj) == 1 .AND. refgeo%mask_ice( ii,jj) == 0) then
            refgeo%mask_coast( i,j) = 1
          end if
        end do
        end do
      end if

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = 10)

  end subroutine calc_reference_geometry_secondary_data

  ! ! Initialise a reference geometry according to an idealised world on the model mesh
  ! subroutine initialise_reference_geometry_idealised_mesh( mesh, refgeo, choice_refgeo_idealised)
  !   ! Initialise a reference geometry according to an idealised world

  !   implicit none

  !   ! In/output variables:
  !   type(type_mesh),                intent(in)    :: mesh
  !   type(type_reference_geometry),  intent(inout) :: refgeo
  !   character(len=256),             intent(in)    :: choice_refgeo_idealised

  !   ! Local variables:
  !   character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh'

  !   ! Add routine to path
  !   call init_routine( routine_name)

  !   ! Allocate shared memory
  !   call allocate_shared_dp_1D( mesh%nV, refgeo%Hi, refgeo%wHi)
  !   call allocate_shared_dp_1D( mesh%nV, refgeo%Hb, refgeo%wHb)
  !   call allocate_shared_dp_1D( mesh%nV, refgeo%Hs, refgeo%wHs)

  !   if     (choice_refgeo_idealised == 'flatearth') then
  !     ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments
  !     call initialise_reference_geometry_idealised_mesh_flatearth(     mesh, refgeo)
  !   elseif (choice_refgeo_idealised == 'Halfar') then
  !     ! The Halfar dome solution at t = 0
  !     call initialise_reference_geometry_idealised_mesh_Halfar(        mesh, refgeo)
  !   elseif (choice_refgeo_idealised == 'Bueler') then
  !     ! The Bueler dome solution at t = 0
  !     call initialise_reference_geometry_idealised_mesh_Bueler(        mesh, refgeo)
  !   elseif (choice_refgeo_idealised == 'SSA_icestream') then
  !     ! The SSA_icestream infinite slab on a flat slope
  !     call initialise_reference_geometry_idealised_mesh_SSA_icestream( mesh, refgeo)
  !   elseif (choice_refgeo_idealised == 'MISMIP_mod') then
  !     ! The MISMIP_mod cone-shaped island
  !     call initialise_reference_geometry_idealised_mesh_MISMIP_mod(    mesh, refgeo)
  !   elseif (choice_refgeo_idealised == 'ISMIP_HOM_A') then
  !     ! The ISMIP-HOM A bumpy slope
  !     call initialise_reference_geometry_idealised_mesh_ISMIP_HOM_A(   mesh, refgeo)
  !   elseif (choice_refgeo_idealised == 'ISMIP_HOM_B') then
  !     ! The ISMIP-HOM B bumpy slope
  !     call initialise_reference_geometry_idealised_mesh_ISMIP_HOM_B(   mesh, refgeo)
  !   elseif (choice_refgeo_idealised == 'ISMIP_HOM_C' .OR. &
  !           choice_refgeo_idealised == 'ISMIP_HOM_D') then
  !     ! The ISMIP-HOM C/D bumpy slope
  !     call initialise_reference_geometry_idealised_mesh_ISMIP_HOM_CD(  mesh, refgeo)
  !   elseif (choice_refgeo_idealised == 'ISMIP_HOM_E') then
  !     ! The ISMIP-HOM E Glacier d'Arolla geometry
  !     call initialise_reference_geometry_idealised_mesh_ISMIP_HOM_E(   mesh, refgeo)
  !   elseif (choice_refgeo_idealised == 'ISMIP_HOM_F') then
  !     ! The ISMIP-HOM A bumpy slope
  !     call initialise_reference_geometry_idealised_mesh_ISMIP_HOM_F(   mesh, refgeo)
  !   elseif (choice_refgeo_idealised == 'MISMIP+') then
  !     ! The MISMIP+ fjord geometry
  !     call initialise_reference_geometry_idealised_mesh_MISMIPplus(    mesh, refgeo)
  !   else
  !     call crash('unknown choice_refgeo_idealised "' // TRIM( choice_refgeo_idealised) // '"!')
  !   end if

  !   ! Finalise routine path
  !   call finalise_routine( routine_name, n_extra_windows_expected = 3)

  ! end subroutine initialise_reference_geometry_idealised_mesh
!   subroutine initialise_reference_geometry_idealised_mesh_flatearth(     mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_flatearth'
!     integer                                       :: vi

!     ! Add routine to path
!     call init_routine( routine_name)

!     do vi = mesh%vi1, mesh%vi2
!       refgeo%Hi( vi) = 0._dp
!       refgeo%Hb( vi) = 0._dp
!       refgeo%Hs( vi) = 0._dp
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_flatearth
!   subroutine initialise_reference_geometry_idealised_mesh_Halfar(        mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The Halfar dome solution at t = 0

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_Halfar'
!     integer                                       :: vi

!     ! Add routine to path
!     call init_routine( routine_name)

!     do vi = mesh%vi1, mesh%vi2
!       refgeo%Hi( vi) = Halfar_solution( mesh%V( vi,1), mesh%V( vi,2), C%start_time_of_run)
!       refgeo%Hb( vi) = 0._dp
!       refgeo%Hs( vi) = refgeo%Hi( vi)
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_Halfar
!   subroutine initialise_reference_geometry_idealised_mesh_Bueler(        mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The Bueler dome solution at t = 0

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_Bueler'
!     integer                                       :: vi

!     ! Add routine to path
!     call init_routine( routine_name)

!     do vi = mesh%vi1, mesh%vi2
!       refgeo%Hi( vi) = Bueler_solution( mesh%V( vi,1), mesh%V( vi,2), C%start_time_of_run)
!       refgeo%Hb( vi) = 0._dp
!       refgeo%Hs( vi) = refgeo%Hi( vi)
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_Bueler
!   subroutine initialise_reference_geometry_idealised_mesh_SSA_icestream( mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The SSA_icestream infinite slab on a flat slope

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_SSA_icestream'
!     integer                                       :: vi

!     ! Add routine to path
!     call init_routine( routine_name)

!     do vi = mesh%vi1, mesh%vi2
!       refgeo%Hi( vi) = 2000._dp
!       refgeo%Hb( vi) = -0.001_dp * mesh%V( vi,1)
!       refgeo%Hs( vi) = refgeo%Hb( vi) + refgeo%Hi( vi)
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_SSA_icestream
!   subroutine initialise_reference_geometry_idealised_mesh_MISMIP_mod(    mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The MISMIP_mod cone-shaped island

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_MISMIP_mod'
!     integer                                       :: vi

!     ! Add routine to path
!     call init_routine( routine_name)

!     do vi = mesh%vi1, mesh%vi2

!       ! Create a nice circular ice shelf
!       if (SQRT(mesh%V( vi,1)**2 + mesh%V( vi,2)**2) < mesh%xmax * 0.95_dp) then
!         refgeo%Hi( vi) = 100._dp
!       else
!         refgeo%Hi( vi) = 0._dp
!       end if

!       refgeo%Hb( vi) = 720._dp - 778.5_dp * SQRT( mesh%V( vi,1)**2 + mesh%V( vi,2)**2)/ 750000._dp
!       refgeo%Hs( vi) = surface_elevation( refgeo%Hi( vi), refgeo%Hb( vi), 0._dp)
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_MISMIP_mod
!   subroutine initialise_reference_geometry_idealised_mesh_ISMIP_HOM_A(   mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The ISMIP-HOM A bumpy slope

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_ISMIP_HOM_A'
!     integer                                       :: vi

!     ! Add routine to path
!     call init_routine( routine_name)

!     do vi = mesh%vi1, mesh%vi2
!       refgeo%Hs( vi) = 2000._dp - mesh%V( vi,1) * TAN( 0.5_dp * pi / 180._dp)
!       refgeo%Hb( vi) = refgeo%Hs( vi) - 1000._dp + 500._dp * SIN( mesh%V( vi,1) * 2._dp * pi / C%ISMIP_HOM_L) * SIN( mesh%V( vi,2) * 2._dp * pi / C%ISMIP_HOM_L)
!       refgeo%Hi( vi) = refgeo%Hs( vi) - refgeo%Hb( vi)
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_ISMIP_HOM_A
!   subroutine initialise_reference_geometry_idealised_mesh_ISMIP_HOM_B(   mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The ISMIP-HOM B bumpy slope

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_ISMIP_HOM_B'
!     integer                                       :: vi

!     ! Add routine to path
!     call init_routine( routine_name)

!     do vi = mesh%vi1, mesh%vi2
!       refgeo%Hs( vi) = 2000._dp - mesh%V( vi,1) * TAN( 0.5_dp * pi / 180._dp)
!       refgeo%Hb( vi) = refgeo%Hs( vi) - 1000._dp + 500._dp * SIN( mesh%V( vi,1) * 2._dp * pi / C%ISMIP_HOM_L)
!       refgeo%Hi( vi) = refgeo%Hs( vi) - refgeo%Hb( vi)
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_ISMIP_HOM_B
!   subroutine initialise_reference_geometry_idealised_mesh_ISMIP_HOM_CD(  mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The ISMIP-HOM C/D bumpy slope

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_ISMIP_HOM_CD'
!     integer                                       :: vi

!     ! Add routine to path
!     call init_routine( routine_name)

!     do vi = mesh%vi1, mesh%vi2
!       refgeo%Hs( vi) = 2000._dp - mesh%V( vi,1) * TAN( 0.1_dp * pi / 180._dp)
!       refgeo%Hb( vi) = refgeo%Hs( vi) - 1000._dp
!       refgeo%Hi( vi) = refgeo%Hs( vi) - refgeo%Hb( vi)
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_ISMIP_HOM_CD
!   subroutine initialise_reference_geometry_idealised_mesh_ISMIP_HOM_E(   mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The ISMIP-HOM E Glacier d'Arolla geometry

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_ISMIP_HOM_E'

!     ! DENK DROM
!     real(dp) :: dp_dummy

!     ! Add routine to path
!     call init_routine( routine_name)

!     ! DENK DROM
!     dp_dummy = mesh%V( 1,1)
!     dp_dummy = refgeo%Hi( 1)
!     call crash('FIXME!')

! !    ! Local variables
! !    integer                                       :: vi
! !    real(dp)                                      :: x,Hs,Hb
! !    integer                                       :: ios,slides
! !
! !    ! To prevent compiler warnings from unused variables
! !    i = grid%nx
! !
! !    ! Read data from external file
! !    if (par%master) then
! !
! !      OPEN( UNIT = 1337, FILE=C%ISMIP_HOM_E_Arolla_filename, ACTION = 'READ')
! !      do i = 1, 51
! !        READ( UNIT = 1337, FMT=*, IOSTAT=ios) x, Hb, Hs, slides
! !        do j = 1, refgeo%grid%ny
! !          refgeo%Hb( vi) = Hb
! !          refgeo%Hi( vi) = Hs - Hb
! !          refgeo%Hs( vi) = Hs
! !        end do
! !        if (ios /= 0) then
! !          WRITE(0,*) ' initialise_reference_geometry_idealised_mesh_ISMIP_HOM_E - ERROR: length of text file "', TRIM(C%ISMIP_HOM_E_Arolla_filename), '" should be 51 lines!'
! !          call MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
! !        end if
! !      end do
! !      CLOSE( UNIT  = 1337)
! !
! !    end if ! if (par%master) then
! !    call sync
! !
! !    ! Finalise routine path
! !    call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_ISMIP_HOM_E
!   subroutine initialise_reference_geometry_idealised_mesh_ISMIP_HOM_F(   mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The ISMIP-HOM A bumpy slope

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_ISMIP_HOM_F'
!     integer                                       :: vi

!     real(dp), parameter                           :: H0    = 1000._dp
!     real(dp), parameter                           :: a0    = 100._dp
!     real(dp), parameter                           :: sigma = 10000._dp

!     ! Add routine to path
!     call init_routine( routine_name)

!     do vi = mesh%vi1, mesh%vi2
!       refgeo%Hs( vi) = 5000._dp - mesh%V( vi,1) * TAN( 3._dp * pi / 180._dp)
!       refgeo%Hb( vi) = refgeo%Hs( vi) - H0 + a0 * EXP( -((mesh%V( vi,1) - 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                            + a0 * EXP( -((mesh%V( vi,1) - 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                            + a0 * EXP( -((mesh%V( vi,1) - 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                            + a0 * EXP( -((mesh%V( vi,1) - 0._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                            + a0 * EXP( -((mesh%V( vi,1) - 0._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                            + a0 * EXP( -((mesh%V( vi,1) - 0._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                            + a0 * EXP( -((mesh%V( vi,1) + 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                            + a0 * EXP( -((mesh%V( vi,1) + 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
!                                            + a0 * EXP( -((mesh%V( vi,1) + 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2)
!       refgeo%Hi( vi) = refgeo%Hs( vi) - refgeo%Hb( vi)
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_ISMIP_HOM_F
!   subroutine initialise_reference_geometry_idealised_mesh_MISMIPplus(    mesh, refgeo)
!     ! Initialise reference geometry according to an idealised world
!     !
!     ! The MISMIpplus fjord geometry

!     implicit none

!     ! In/output variables:
!     type(type_mesh),                intent(in)    :: mesh
!     type(type_reference_geometry),  intent(inout) :: refgeo

!     ! Local variables:
!     character(len=256), parameter                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_MISMIPplus'
!     integer                                       :: vi

!     real(dp)                                      :: x,y,xtilde,Bx,By
!     real(dp), parameter                           :: B0     = -150._dp
!     real(dp), parameter                           :: B2     = -728.8_dp
!     real(dp), parameter                           :: B4     = 343.91_dp
!     real(dp), parameter                           :: B6     = -50.57_dp
!     real(dp), parameter                           :: xbar   = 300000._dp
!     real(dp), parameter                           :: fc     = 4000._dp
!     real(dp), parameter                           :: dc     = 500._dp
!     real(dp), parameter                           :: wc     = 24000._dp
!     real(dp), parameter                           :: zbdeep = -720._dp

!     ! Add routine to path
!     call init_routine( routine_name)

!     do vi = mesh%vi1, mesh%vi2
!       x = mesh%V( vi,1) + 400000._dp
!       y = mesh%V( vi,2)
!       xtilde = x / xbar
!       Bx = B0 + (B2 * xtilde**2._dp) + (B4 * xtilde**4._dp) + (B6 * xtilde**6._dp)
!       By = (dc / (1 + EXP(-2._dp*(y - wc)/fc))) + &
!            (dc / (1 + EXP( 2._dp*(y + wc)/fc)))
!       refgeo%Hi( vi) = 100._dp
!       refgeo%Hb( vi) = MAX( Bx + By, zbdeep)
!       refgeo%Hs( vi) = surface_elevation( refgeo%Hi( vi), refgeo%Hb( vi), 0._dp)
!     end do
!     call sync

!     ! Finalise routine path
!     call finalise_routine( routine_name)

!   end subroutine initialise_reference_geometry_idealised_mesh_MISMIPplus

  ! Map init and PD references data from their supplied grids to the model mesh
  subroutine map_reference_geometries_to_mesh( region, mesh)
    ! Map the initial, present-day, and GIAeq reference geometries from their original
    ! square grids to the model mesh.
    !
    ! Since calculating remapping operators can take some time, and since the three
    ! square grids are often identical, some time can be saved by calculating the
    ! operators only once.

    implicit none

    ! Input and output variables
    type(type_model_region),        intent(inout) :: region
    type(type_mesh),                intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'map_reference_geometries_to_mesh'
    character(len=256)                            :: choice_refgeo_init, choice_refgeo_PD, choice_refgeo_GIAeq
    logical                                       :: did_remap_init, did_remap_PD, did_remap_GIAeq, do_reuse_init_map, do_reuse_PD_map

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine parameters for this region
    if     (region%name == 'NAM') then
      choice_refgeo_init    = C%choice_refgeo_init_NAM
      choice_refgeo_PD      = C%choice_refgeo_PD_NAM
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_NAM
    elseif (region%name == 'EAS') then
      choice_refgeo_init    = C%choice_refgeo_init_EAS
      choice_refgeo_PD      = C%choice_refgeo_PD_EAS
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_EAS
    elseif (region%name == 'GRL') then
      choice_refgeo_init    = C%choice_refgeo_init_GRL
      choice_refgeo_PD      = C%choice_refgeo_PD_GRL
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_GRL
    elseif (region%name == 'ANT') then
      choice_refgeo_init    = C%choice_refgeo_init_ANT
      choice_refgeo_PD      = C%choice_refgeo_PD_ANT
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_ANT
    end if

    ! Initial ice-sheet geometry
    ! ==========================

    if (par%master) then
      write(*,"(A)") '  Mapping initial reference geometry to the mesh...'
    end if

    did_remap_init = .FALSE.

    ! if     (choice_refgeo_init == 'idealised') then
    !   ! For idealised geometries, calculate the exact solution directly on the mesh instead of remapping it from the grid

    !   call initialise_reference_geometry_idealised_mesh( mesh, region%refgeo_init, C%choice_refgeo_init_idealised)

    ! elseif (choice_refgeo_init == 'realistic') then
    !   ! For realistic geometries, remap the data from the grid

      call calc_remapping_operator_grid2mesh( region%refgeo_init%grid, mesh)
      call map_reference_geometry_to_mesh( mesh, region%refgeo_init)
      did_remap_init = .TRUE.

    ! else
    !   call crash('unknown choice_refgeo_init "' // TRIM( choice_refgeo_init) // '"!')
    ! end if

    ! Present-day ice-sheet geometry
    ! ==============================

    if (par%master) then
      write(*,"(A)") '  Mapping present-day reference geometry to the mesh...'
    end if

    did_remap_PD = .FALSE.

    ! if     (choice_refgeo_PD == 'idealised') then
    !   ! For idealised geometries, calculate the exact solution directly on the mesh instead of remapping it from the grid

    !   call initialise_reference_geometry_idealised_mesh( mesh, region%refgeo_PD, C%choice_refgeo_PD_idealised)

    ! elseif (choice_refgeo_PD == 'realistic') then
      ! For realistic geometries, remap the data from the grid

      ! Check if we can re-use the remapping arrays from the initial geometry
      do_reuse_init_map = .FALSE.
      if (did_remap_init) then
        if (region%refgeo_PD%grid%xmin == region%refgeo_init%grid%xmin .AND. &
            region%refgeo_PD%grid%xmax == region%refgeo_init%grid%xmax .AND. &
            region%refgeo_PD%grid%ymin == region%refgeo_init%grid%ymin .AND. &
            region%refgeo_PD%grid%ymax == region%refgeo_init%grid%ymax .AND. &
            region%refgeo_PD%grid%dx   == region%refgeo_init%grid%dx   .AND. &
            region%refgeo_PD%grid%nx   == region%refgeo_init%grid%nx   .AND. &
            region%refgeo_PD%grid%ny   == region%refgeo_init%grid%ny) then
          do_reuse_init_map = .TRUE.
        end if
      end if
      if (do_reuse_init_map) then
        call MatDuplicate( region%refgeo_init%grid%M_map_grid2mesh, MAT_COPY_VALUES, region%refgeo_PD%grid%M_map_grid2mesh, perr)
      else
        call calc_remapping_operator_grid2mesh( region%refgeo_PD%grid, mesh)
      end if

      call map_reference_geometry_to_mesh( mesh, region%refgeo_PD)
      did_remap_PD = .TRUE.

    ! else
    !   call crash('unknown choice_refgeo_PD "' // TRIM( choice_refgeo_PD) // '"!')
    ! end if

    ! GIA equilibrium ice-sheet geometry
    ! ==================================

    if (par%master) then
      write(*,"(A)") '  Mapping GIA equilibrium reference geometry to the mesh...'
    end if

    did_remap_GIAeq = .FALSE.

    ! if     (choice_refgeo_GIAeq == 'idealised') then
    !   ! For idealised geometries, calculate the exact solution directly on the mesh instead of remapping it from the grid

    !   call initialise_reference_geometry_idealised_mesh( mesh, region%refgeo_GIAeq, C%choice_refgeo_GIAeq_idealised)

    ! elseif (choice_refgeo_GIAeq == 'realistic') then
      ! For realistic geometries, remap the data from the grid

      ! Check if we can re-use the remapping arrays from the initial/PD geometry
      do_reuse_init_map = .FALSE.
      if (did_remap_init) then
        if (region%refgeo_GIAeq%grid%xmin == region%refgeo_init%grid%xmin .AND. &
            region%refgeo_GIAeq%grid%xmax == region%refgeo_init%grid%xmax .AND. &
            region%refgeo_GIAeq%grid%ymin == region%refgeo_init%grid%ymin .AND. &
            region%refgeo_GIAeq%grid%ymax == region%refgeo_init%grid%ymax .AND. &
            region%refgeo_GIAeq%grid%dx   == region%refgeo_init%grid%dx   .AND. &
            region%refgeo_GIAeq%grid%nx   == region%refgeo_init%grid%nx   .AND. &
            region%refgeo_GIAeq%grid%ny   == region%refgeo_init%grid%ny) then
          do_reuse_init_map = .TRUE.
        end if
      end if
      do_reuse_PD_map = .FALSE.
      if (did_remap_init) then
        if (region%refgeo_GIAeq%grid%xmin == region%refgeo_PD%grid%xmin .AND. &
            region%refgeo_GIAeq%grid%xmax == region%refgeo_PD%grid%xmax .AND. &
            region%refgeo_GIAeq%grid%ymin == region%refgeo_PD%grid%ymin .AND. &
            region%refgeo_GIAeq%grid%ymax == region%refgeo_PD%grid%ymax .AND. &
            region%refgeo_GIAeq%grid%dx   == region%refgeo_PD%grid%dx   .AND. &
            region%refgeo_GIAeq%grid%nx   == region%refgeo_PD%grid%nx   .AND. &
            region%refgeo_GIAeq%grid%ny   == region%refgeo_PD%grid%ny) then
          do_reuse_PD_map = .TRUE.
        end if
      end if
      if     (do_reuse_init_map) then
        call MatDuplicate( region%refgeo_init%grid%M_map_grid2mesh, MAT_COPY_VALUES, region%refgeo_GIAeq%grid%M_map_grid2mesh, perr)
      elseif (do_reuse_PD_map) then
        call MatDuplicate( region%refgeo_PD%grid%M_map_grid2mesh  , MAT_COPY_VALUES, region%refgeo_GIAeq%grid%M_map_grid2mesh, perr)
      else
        call calc_remapping_operator_grid2mesh( region%refgeo_GIAeq%grid, mesh)
      end if

      call map_reference_geometry_to_mesh( mesh, region%refgeo_GIAeq)
      did_remap_GIAeq = .TRUE.

    ! else
    !   call crash('unknown choice_refgeo_GIAeq "' // TRIM( choice_refgeo_GIAeq) // '"!')
    ! end if

    ! Clean up after yourself
    if (did_remap_init ) call deallocate_remapping_operators_grid2mesh( region%refgeo_init%grid )
    if (did_remap_PD   ) call deallocate_remapping_operators_grid2mesh( region%refgeo_PD%grid   )
    if (did_remap_GIAeq) call deallocate_remapping_operators_grid2mesh( region%refgeo_GIAeq%grid)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = 9)

  end subroutine map_reference_geometries_to_mesh
  subroutine map_reference_geometry_to_mesh( mesh, refgeo)
    ! Map data for a single reference geometry from its original square grid to the model mesh.

    implicit none

    ! Input and output variables
    type(type_mesh),                intent(inout) :: mesh
    type(type_reference_geometry),  intent(inout) :: refgeo

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'map_reference_geometry_to_mesh'
    integer                                       :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Map PD data to the mesh
    call map_grid2mesh_2D_partial( refgeo%grid, mesh, refgeo%Hi_grid, refgeo%Hi)
    call map_grid2mesh_2D_partial( refgeo%grid, mesh, refgeo%Hb_grid, refgeo%Hb)

    do vi = mesh%vi1, mesh%vi2
      refgeo%Hs( vi) = surface_elevation( refgeo%Hi( vi), refgeo%Hb( vi), 0._dp)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_reference_geometry_to_mesh

!   ! Analytical solutions used to initialise some benchmark experiments
!   FUNCTION Halfar_solution( x, y, t) RESULT(H)
!     ! Describes an ice-sheet at time t (in years) conforming to the Halfar similarity function
!     ! with dome thickness H0 and margin radius R0 at t0. Used to initialise the model
!     ! for the Halfar solution test run

!     implicit none

!     ! Input variables
!     real(dp), intent(in) :: x  ! x coordinate [m]
!     real(dp), intent(in) :: y  ! y coordinate [m]
!     real(dp), intent(in) :: t  ! Time from t0 [years]

!     ! Result
!     real(dp)             :: H  ! Ice thickness at [x,y] at t=0 [m]

!     ! Local variables
!     real(dp) :: A_flow, rho, g, Gamma, t0, r, f1, f2, f3, tp

!     real(dp), parameter :: H0 = 5000._dp   ! Ice dome thickness at t=0 [m]
!     real(dp), parameter :: R0 = 300000._dp ! Ice margin radius  at t=0 [m]

!     A_flow  = 1E-16_dp
!     rho     = 910._dp
!     g       = 9.81_dp

!     Gamma = (2._dp / 5._dp) * (A_flow / sec_per_year) * (rho * g)**3._dp
!     t0 = 1._dp / (18._dp * Gamma) * (7._dp/4._dp)**3._dp * (R0**4._dp)/(H0**7._dp)

!     tp = (t * sec_per_year) + t0

!     r = SQRT(x**2._dp + y**2._dp)

!     f1 = (t0/tp)**(1._dp/9._dp)
!     f2 = (t0/tp)**(1._dp/18._dp)
!     f3 = (r/R0)

!     H = H0 * f1 * MAX(0._dp, (1._dp - (f2*f3)**(4._dp/3._dp)))**(3._dp/7._dp)

!   END FUNCTION Halfar_solution
!   FUNCTION Bueler_solution( x, y, t) RESULT(H)
!     ! Describes an ice-sheet at time t (in years) conforming to the Bueler solution
!     ! with dome thickness H0 and margin radius R0 at t0, with a surface mass balance
!     ! determined by lambda. Used to intialise the model for the Bueler solution test run

!     implicit none

!     ! Input variables
!     real(dp), intent(in) :: x       ! x coordinate [m]
!     real(dp), intent(in) :: y       ! y coordinate [m]
!     real(dp), intent(in) :: t       ! Time from t0 [years]

!     ! Result
!     real(dp)             :: H  ! Ice thickness at [x,y] at t=0 [m]

!     ! Local variables
!     real(dp) :: A_flow, rho, g, n, alpha, beta, Gamma, f1, f2, t0, tp, f3, f4

!     real(dp), parameter :: H0     = 3000._dp    ! Ice dome thickness at t=0 [m]
!     real(dp), parameter :: R0     = 500000._dp  ! Ice margin radius  at t=0 [m]
!     real(dp), parameter :: lambda = 5.0_dp      ! Mass balance parameter

!     A_flow  = 1E-16_dp
!     rho     = 910._dp
!     g       = 9.81_dp
!     n       = 3._dp

!     alpha = (2._dp - (n+1._dp)*lambda) / ((5._dp*n)+3._dp)
!     beta  = (1._dp + ((2._dp*n)+1._dp)*lambda) / ((5._dp*n)+3._dp)
!     Gamma = 2._dp/5._dp * (A_flow/sec_per_year) * (rho * g)**n

!     f1 = ((2._dp*n)+1)/(n+1._dp)
!     f2 = (R0**(n+1._dp))/(H0**((2._dp*n)+1._dp))
!     t0 = (beta / Gamma) * (f1**n) * f2

!     !tp = (t * sec_per_year) + t0; % Actual equation needs t in seconds from zero , but we want to supply t in years from t0
!     tp = t * sec_per_year

!     f1 = (tp / t0)**(-alpha)
!     f2 = (tp / t0)**(-beta)
!     f3 = SQRT( (x**2._dp) + (y**2._dp) )/R0
!     f4 = MAX(0._dp, 1._dp - (f2*f3)**((n+1._dp)/n))
!     H = H0 * f1 * f4**(n/((2._dp*n)+1._dp))

!     !M = (lambda / tp) * H * sec_per_year

!   END FUNCTION Bueler_solution

END MODULE reference_fields_module
