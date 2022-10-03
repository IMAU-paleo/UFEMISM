module reference_fields_module
  ! Contains the routines for setting up the three "reference geometries":
  ! - refgeo_PD:     present-day, used to calculate sea-level contribution, isotope change, and more
  ! - refgeo_init:   initial, used to initialise the simulation
  ! - refgeo_GIA_eq: GIA equilibrium, used for the GIA model

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

  use mpi
  use petscksp
  use configuration_module, only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use petsc_module,         only : perr
  use parallel_module,      only : par, sync, ierr, cerr, partition_list
  use data_types_module,    only : type_reference_geometry, type_grid, &
                                   type_model_region, type_mesh
  use netcdf_module,        only : inquire_reference_geometry_file, read_reference_geometry_file
  use mesh_mapping_module,  only : calc_remapping_operator_grid2mesh, map_grid2mesh_2D, &
                                   map_grid2mesh_2D_partial, deallocate_remapping_operators_grid2mesh
  use utilities_module,     only : check_for_NaN_dp_2D, remove_Lake_Vostok, &
                                   is_floating, surface_elevation, oblique_sg_projection

  implicit none

contains

! ===== Initialisation =====
! ==========================

  ! Initialise all three reference geometries
  subroutine initialise_reference_geometries( refgeo_init, refgeo_PD, refgeo_GIAeq, region_name)
    ! Initialise all three reference geometries

    implicit none

    ! In/output variables:
    type(type_reference_geometry), intent(inout) :: refgeo_init
    type(type_reference_geometry), intent(inout) :: refgeo_PD
    type(type_reference_geometry), intent(inout) :: refgeo_GIAeq
    character(len=3),              intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'initialise_reference_geometries'
    character(len=256)                           :: choice_refgeo_init, choice_refgeo_PD, choice_refgeo_GIAeq
    character(len=256)                           :: filename_refgeo_init, filename_refgeo_PD, filename_refgeo_GIAeq
    real(dp)                                     :: time_to_restart_from

    ! Initialisation
    ! ==============

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

    ! Present-day ice-sheet geometry
    ! ==============================

    select case (choice_refgeo_PD)

      case ('idealised')
        ! Idealised geometry
        if (par%master) then
          write(*,'(3A)') '  Initialising present-day     reference geometry from idealised case "', &
                             trim(C%choice_refgeo_PD_idealised), '"...'
        end if
        call sync

        ! WIP
        call crash('choice_refgeo_PD - idealised not implemented yet...')

      case ('realistic')
        ! Realistic geometry
        if (par%master) then
          write(*,"(3A)") '  Initialising present-day     reference geometry from file ', &
                             trim( filename_refgeo_PD), '...'
        end if
        call sync

        ! Initialise PD geometry from a file
        call initialise_reference_geometry_from_file( refgeo_PD, filename_refgeo_PD, region_name)

      case default
        ! Unkown option
        call crash('unknown choice_refgeo_PD "' // trim( choice_refgeo_PD) // '"!')

    end select

    ! Initial ice-sheet geometry
    ! ==========================

    select case (choice_refgeo_init)

      case ('idealised')
        ! Idealised geometry

        if (par%master) then
          write(*,'(3A)') '  Initialising initial         reference geometry from idealised case "', &
                             trim(C%choice_refgeo_init_idealised), '"...'
        end if
        call sync

        ! WIP
        call crash('choice_refgeo_init - idealised not implemented yet...')

      case ('realistic')
        ! Realistic geometry

        if (par%master) then
          write(*,"(3A)") '  Initialising initial         reference geometry from file ', &
                             trim( filename_refgeo_init), '...'
        end if
        call sync

        ! Initialise initial geometry from a file
        call initialise_reference_geometry_from_file( refgeo_init, filename_refgeo_init, region_name)

      case ('restart')
        ! Geometry from a previous simulation

        ! WIP
        call crash('choice_refgeo_init - restart not implemented yet...')

      case default
        ! Unkown option
        call crash('unknown choice_refgeo_init "' // trim( choice_refgeo_init) // '"!')

    end select

    ! GIA equilibrium ice-sheet geometry
    ! ==================================

    select case (choice_refgeo_GIAeq)

      case ('idealised')
        ! Idealised geometry

        if (par%master) then
          write(*,'(3A)') '  Initialising GIA equilibrium reference geometry from idealised case "', &
                             trim(C%choice_refgeo_GIAeq_idealised), '"...'
        end if
        call sync

        ! WIP
        call crash('choice_refgeo_GIAeq - idealised not implemented yet...')

      case ('realistic')
        ! Realistic geometry
        if (par%master) then
          write(*,"(3A)") '  Initialising GIA equilibrium reference geometry from file ', &
                             trim( filename_refgeo_GIAeq), '...'
        end if
        call sync

        ! Initialise GIA-equilibrium geometry from a file
        call initialise_reference_geometry_from_file( refgeo_GIAeq, filename_refgeo_GIAeq, region_name)

      case default
        ! Unkown option
        call crash('unknown choice_refgeo_init "' // trim( choice_refgeo_GIAeq) // '"!')

    end select

    ! Secondary data
    ! ==============

    ! Fill in secondary data for the reference geometry (used to force mesh creation)
    call calc_reference_geometry_secondary_data( refgeo_init%grid , refgeo_init )
    call calc_reference_geometry_secondary_data( refgeo_PD%grid   , refgeo_PD   )
    call calc_reference_geometry_secondary_data( refgeo_GIAeq%grid, refgeo_GIAeq)

    ! Finalisation
    ! ============

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometries

! ===== Realistic reference geometries (grid) =====
! =================================================

  subroutine initialise_reference_geometry_from_file( refgeo, filename_refgeo, region_name)
    ! Initialise a reference geometry with data from a NetCDF file

    implicit none

    ! In/output variables:
    type(type_reference_geometry), intent(inout) :: refgeo
    character(len=256),            intent(in)    :: filename_refgeo
    character(len=3),              intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'initialise_reference_geometry_from_file'
    integer                                      :: i,j,n
    real(dp), parameter                          :: tol = 1E-9_dp

    ! Add routine to path
    call init_routine( routine_name)

    if     (region_name == 'NAM') then
      refgeo%grid%lambda_M     = C%lambda_M_NAM
      refgeo%grid%phi_M        = C%phi_M_NAM
      refgeo%grid%alpha_stereo = C%alpha_stereo_NAM
    elseif (region_name == 'EAS') then
      refgeo%grid%lambda_M     = C%lambda_M_EAS
      refgeo%grid%phi_M        = C%phi_M_EAS
      refgeo%grid%alpha_stereo = C%alpha_stereo_EAS
    elseif (region_name == 'GRL') then
      refgeo%grid%lambda_M     = C%lambda_M_GRL
      refgeo%grid%phi_M        = C%phi_M_GRL
      refgeo%grid%alpha_stereo = C%alpha_stereo_GRL
    elseif (region_name == 'ANT') then
      refgeo%grid%lambda_M     = C%lambda_M_ANT
      refgeo%grid%phi_M        = C%phi_M_ANT
      refgeo%grid%alpha_stereo = C%alpha_stereo_ANT
    end if

    ! Inquire if all the required fields are present in the specified NetCDF file,
    ! and determine the dimensions of the memory to be allocated.
    refgeo%netcdf%filename = filename_refgeo
    call inquire_reference_geometry_file( refgeo)

    ! Assign range to each processor
    call partition_list( refgeo%grid%nx, par%i, par%n, refgeo%grid%i1, refgeo%grid%i2)
    call partition_list( refgeo%grid%ny, par%i, par%n, refgeo%grid%j1, refgeo%grid%j2)

    ! Allocate memory for raw data
    allocate( refgeo%grid%x ( refgeo%grid%nx ))
    allocate( refgeo%grid%y ( refgeo%grid%ny ))

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
    if (region_name == 'ANT'.and. C%remove_Lake_Vostok) then
      call remove_Lake_Vostok( refgeo%grid%x, refgeo%grid%y, refgeo%Hi_grid, refgeo%Hb_grid, refgeo%Hs_grid)
    end if

    ! Remove ice based on the no-ice masks (grid versions)
    call apply_mask_noice_grid( refgeo, region_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_reference_geometry_from_file

! ===== Secondary data =====
! ==========================

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

      if (i == 1 .or. i == grid%nx .or. j == 1 .or. j == grid%ny) then
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
    call finalise_routine( routine_name)

  end subroutine calc_reference_geometry_secondary_data

! ===== Mesh mapping =====
! ========================

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
    type(type_model_region),       intent(inout) :: region
    type(type_mesh),               intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'map_reference_geometries_to_mesh'
    character(len=256)                           :: choice_refgeo_init, choice_refgeo_PD, choice_refgeo_GIAeq
    logical                                      :: did_remap_init, did_remap_PD, did_remap_GIAeq
    logical                                      :: do_reuse_init_map, do_reuse_PD_map

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
    call sync

    did_remap_init = .false.

    select case (choice_refgeo_init)

      case ('idealised')
        ! Idealised geometry

        ! Calculate the exact solution directly on the mesh
        call crash('choice_refgeo_init - idealised not implemented yet...')

      case ('realistic')
        ! Realistic geometry

        ! Remap data from the grid
        call calc_remapping_operator_grid2mesh( region%refgeo_init%grid, mesh)
        call map_reference_geometry_to_mesh( mesh, region%refgeo_init)
        did_remap_init = .true.

      case ('restart')
        ! Geometry from a previous simulation

        ! WIP
        call crash('choice_refgeo_init - restart not implemented yet...')

        if (region%time == C%start_time_of_run) then
          ! Initialisation of the model

          if (par%master) then
            write(*,'(A)') '   Copying initial reference geometry from the restart mesh...'
          end if
          call sync

          ! ! Get initial geometry from restart data
          ! region%refgeo_init%Hi( mesh%vi1:mesh%vi2) = region%restart%Hi( mesh%vi1:mesh%vi2)
          ! region%refgeo_init%Hb( mesh%vi1:mesh%vi2) = region%restart%Hb( mesh%vi1:mesh%vi2)
          ! region%refgeo_init%Hs( mesh%vi1:mesh%vi2) = region%restart%Hs( mesh%vi1:mesh%vi2)

        else
          ! Mesh update

          if (par%master) then
            write(*,'(A)') '   Mapping initial reference geometry to the mesh...'
          end if
          call sync

          ! Remap initial geometry data from the grid
          call calc_remapping_operator_grid2mesh( region%refgeo_init%grid, mesh)
          call map_reference_geometry_to_mesh( mesh, region%refgeo_init)
          did_remap_init = .true.

      end IF

      case default
        ! Unknown option
        call crash('unknown choice_refgeo_init "' // trim( choice_refgeo_init) // '"!')

    end select

    ! Present-day ice-sheet geometry
    ! ==============================

    if (par%master) then
      write(*,"(A)") '  Mapping present-day reference geometry to the mesh...'
    end if
    call sync

    did_remap_PD = .false.

    select case (choice_refgeo_PD)

      case ('idealised')
        ! Idealised geometry

        ! Calculate the exact solution directly on the mesh
        call crash('choice_refgeo_PD - idealised not implemented yet...')

      case ('realistic')
        ! Realistic geometry

        ! Check if we can re-use the remapping arrays from the initial geometry
        do_reuse_init_map = .false.
        if (did_remap_init) then
          if (region%refgeo_PD%grid%xmin == region%refgeo_init%grid%xmin .and. &
              region%refgeo_PD%grid%xmax == region%refgeo_init%grid%xmax .and. &
              region%refgeo_PD%grid%ymin == region%refgeo_init%grid%ymin .and. &
              region%refgeo_PD%grid%ymax == region%refgeo_init%grid%ymax .and. &
              region%refgeo_PD%grid%dx   == region%refgeo_init%grid%dx   .and. &
              region%refgeo_PD%grid%nx   == region%refgeo_init%grid%nx   .and. &
              region%refgeo_PD%grid%ny   == region%refgeo_init%grid%ny) then
            do_reuse_init_map = .true.
          end if
        end if

        if (do_reuse_init_map) then
          call MatDuplicate( region%refgeo_init%grid%M_map_grid2mesh, MAT_COPY_VALUES, &
                             region%refgeo_PD%grid%M_map_grid2mesh, perr)
        else
          call calc_remapping_operator_grid2mesh( region%refgeo_PD%grid, mesh)
        end if

        ! Remap data from the grid
        call map_reference_geometry_to_mesh( mesh, region%refgeo_PD)
        did_remap_PD = .true.

      case default
        ! Unknown option
        call crash('unknown choice_refgeo_PD "' // trim( choice_refgeo_PD) // '"!')

    end select

    ! GIA equilibrium ice-sheet geometry
    ! ==================================

    if (par%master) then
      write(*,"(A)") '  Mapping GIA equilibrium reference geometry to the mesh...'
    end if
    call sync

    did_remap_GIAeq = .false.

    select case (choice_refgeo_GIAeq)

      case ('idealised')
        ! Idealised geometry

        ! Calculate the exact solution directly on the mesh
        call crash('choice_refgeo_GIAeq - idealised not implemented yet...')

      case ('realistic')
        ! Realistic geometry

        ! Check if we can re-use the remapping arrays from the initial geometry
        do_reuse_init_map = .false.
        if (did_remap_init) then
          if (region%refgeo_GIAeq%grid%xmin == region%refgeo_init%grid%xmin .and. &
              region%refgeo_GIAeq%grid%xmax == region%refgeo_init%grid%xmax .and. &
              region%refgeo_GIAeq%grid%ymin == region%refgeo_init%grid%ymin .and. &
              region%refgeo_GIAeq%grid%ymax == region%refgeo_init%grid%ymax .and. &
              region%refgeo_GIAeq%grid%dx   == region%refgeo_init%grid%dx   .and. &
              region%refgeo_GIAeq%grid%nx   == region%refgeo_init%grid%nx   .and. &
              region%refgeo_GIAeq%grid%ny   == region%refgeo_init%grid%ny) then
            do_reuse_init_map = .true.
          end if
        end if

        ! Check if we can re-use the remapping arrays from the PD geometry
        do_reuse_PD_map = .false.
        if (did_remap_PD) then
          if (region%refgeo_GIAeq%grid%xmin == region%refgeo_PD%grid%xmin .and. &
              region%refgeo_GIAeq%grid%xmax == region%refgeo_PD%grid%xmax .and. &
              region%refgeo_GIAeq%grid%ymin == region%refgeo_PD%grid%ymin .and. &
              region%refgeo_GIAeq%grid%ymax == region%refgeo_PD%grid%ymax .and. &
              region%refgeo_GIAeq%grid%dx   == region%refgeo_PD%grid%dx   .and. &
              region%refgeo_GIAeq%grid%nx   == region%refgeo_PD%grid%nx   .and. &
              region%refgeo_GIAeq%grid%ny   == region%refgeo_PD%grid%ny) then
            do_reuse_PD_map = .true.
          end if
        end if

        ! Get remapping arrays
        if (do_reuse_init_map) then
          call MatDuplicate( region%refgeo_init%grid%M_map_grid2mesh, MAT_COPY_VALUES, &
                             region%refgeo_GIAeq%grid%M_map_grid2mesh, perr)
        elseif (do_reuse_PD_map) then
          call MatDuplicate( region%refgeo_PD%grid%M_map_grid2mesh,   MAT_COPY_VALUES, &
                             region%refgeo_GIAeq%grid%M_map_grid2mesh, perr)
        else
          call calc_remapping_operator_grid2mesh( region%refgeo_GIAeq%grid, mesh)
        end if

        ! Remap data from the grid
        call map_reference_geometry_to_mesh( mesh, region%refgeo_GIAeq)
        did_remap_GIAeq = .true.

      case default
        ! Unknown option
        call crash('unknown choice_refgeo_GIAeq "' // trim( choice_refgeo_GIAeq) // '"!')

    end select

    ! Finalisation
    ! ============

    ! Clean up after yourself
    if (did_remap_init ) call deallocate_remapping_operators_grid2mesh( region%refgeo_init%grid )
    if (did_remap_PD   ) call deallocate_remapping_operators_grid2mesh( region%refgeo_PD%grid   )
    if (did_remap_GIAeq) call deallocate_remapping_operators_grid2mesh( region%refgeo_GIAeq%grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

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

! ===== Clean-up of realistic geometries (grid) =====
! ===================================================

  subroutine apply_mask_noice_grid( refgeo, region_name)
    ! Remove ice from a certain area. This is used to remove
    ! Greenland from NAM and EAS, and Ellesmere Island from GRL.

    implicit none

    ! In- and output variables
    type(type_reference_geometry), intent(inout) :: refgeo
    character(len=3),              intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'apply_mask_noice_grid'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Apply no-ice mask ===
    ! =========================

    if     (region_name == 'NAM') THEN
      ! Clean up the NAM domain

      select case (C%choice_mask_noice_NAM)

        case ('none')
          ! No no-ice mask is defined for North America

        case ('NAM_remove_GRL')
          ! WIP
          call crash('No-ice mask for "' // region_name // '" not implemented yet!')
          ! Prevent ice growth in the Greenlandic part of the North America domain
          call apply_mask_noice_NAM_remove_GRL_grid( refgeo)

        case default
          ! Unknown case
          call crash('unknown choice_mask_noice_NAM "' // trim( C%choice_mask_noice_NAM) // '"!')

      end select

    elseif (region_name == 'EAS') then
      ! Clean up the EAS domain

      select case (C%choice_mask_noice_EAS)

        case ('none')
          ! No no-ice mask is defined for Eurasia

        case ('EAS_remove_GRL')
          ! WIP
          call crash('No-ice mask for "' // region_name // '" not implemented yet!')
          ! Prevent ice growth in the Greenlandic part of the Eurasian domain
          call apply_mask_noice_EAS_remove_GRL_grid( refgeo)

        case default
          ! Unknown case
          call crash('unknown choice_mask_noice_EAS "' // trim( C%choice_mask_noice_EAS) // '"!')

      end select

    elseif (region_name == 'GRL') then
      ! Clean up the GRL domain

      select case (C%choice_mask_noice_GRL)

        case ('none')
          ! No no-ice mask is defined for Greenland

        case ('GRL_remove_Ellesmere')
          ! Prevent ice growth in the Ellesmere Island part of the Greenland domain
          call apply_mask_noice_GRL_remove_Ellesmere_grid( refgeo)

        case default
          ! Unknown case
          call crash('unknown choice_mask_noice_GRL "' // trim( C%choice_mask_noice_GRL) // '"!')

      end select

    elseif (region_name == 'ANT') then
      ! Clean up the ANT domain

      select case (C%choice_mask_noice_ANT)

        case ('none')
          ! No no-ice mask is defined for Antarctica

        case default
          ! Unknown case
          call crash('unknown choice_mask_noice_ANT "' // trim( C%choice_mask_noice_ANT) // '"!')

      end select

    end if

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_mask_noice_grid

  subroutine apply_mask_noice_NAM_remove_GRL_grid( refgeo)
    ! Remove ice from the Greenlandic part of the North America domain

    implicit none

    ! In- and output variables
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(LEN=256), parameter                :: routine_name = 'apply_mask_noice_NAM_remove_GRL_grid'
    integer                                      :: i,j
    real(dp), dimension(2)                       :: pa, pb
    real(dp)                                     :: yl_ab

    ! Add routine to path
    call init_routine( routine_name)

    pa = [ 490000._dp, 1530000._dp]
    pb = [2030000._dp,  570000._dp]

    do i = 1, refgeo%grid%nx

      yl_ab = pa(2) + (refgeo%grid%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))

      do j = 1, refgeo%grid%ny

        if (refgeo%grid%y(j) > yl_ab .and. refgeo%grid%x(i) > pa(1) .and. refgeo%grid%y(j) > pb(2)) then

          refgeo%Hi_grid( i,j) = 0._dp

        end if

      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_mask_noice_NAM_remove_GRL_grid

  subroutine apply_mask_noice_EAS_remove_GRL_grid( refgeo)
    ! Remove ice from the Greenlandic part of the Eurasia domain

    implicit none

    ! In- and output variables
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'apply_mask_noice_EAS_remove_GRL_grid'
    integer                                      :: i,j
    real(dp), dimension(2)                       :: pa, pb, pc, pd
    real(dp)                                     :: yl_ab, yl_bc, yl_cd

    ! Add routine to path
    call init_routine( routine_name)

    pa = [-2900000._dp, 1300000._dp]
    pb = [-1895000._dp,  900000._dp]
    pc = [ -835000._dp, 1135000._dp]
    pd = [ -400000._dp, 1855000._dp]

    do i = 1, refgeo%grid%nx

      yl_ab = pa(2) + (refgeo%grid%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
      yl_bc = pb(2) + (refgeo%grid%x(i) - pb(1))*(pc(2)-pb(2))/(pc(1)-pb(1))
      yl_cd = pc(2) + (refgeo%grid%x(i) - pc(1))*(pd(2)-pc(2))/(pd(1)-pc(1))

      do j = 1, refgeo%grid%ny

        if ((refgeo%grid%x(i) <  pa(1) .and. refgeo%grid%y(j) > pa(2)) .or. &
            (refgeo%grid%x(i) >= pa(1) .and. refgeo%grid%x(i) < pb(1) .and. refgeo%grid%y(j) > yl_ab) .or. &
            (refgeo%grid%x(i) >= pb(1) .and. refgeo%grid%x(i) < pc(1) .and. refgeo%grid%y(j) > yl_bc) .or. &
            (refgeo%grid%x(i) >= pc(1) .and. refgeo%grid%x(i) < pd(1) .and. refgeo%grid%y(j) > yl_cd)) then

          refgeo%Hi_grid( i,j) = 0._dp

        end if

      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_mask_noice_EAS_remove_GRL_grid

  subroutine apply_mask_noice_GRL_remove_Ellesmere_grid( refgeo)
    ! Remove ice from the Ellesmere Island part of the Greenland domain

    implicit none

    ! In- and output variables
    type(type_reference_geometry), intent(inout) :: refgeo

    ! Local variables:
    character(LEN=256), parameter                :: routine_name = 'apply_mask_noice_GRL_remove_Ellesmere_grid'
    integer                                      :: i,j
    real(dp), dimension(2)                       :: pa_latlon, pb_latlon
    real(dp)                                     :: xa,ya,xb,yb
    real(dp), dimension(2)                       :: pa, pb
    real(dp)                                     :: yl_ab

    ! Add routine to path
    call init_routine( routine_name)

    ! The two endpoints in lat,lon
    pa_latlon = [76.74_dp, -74.79_dp]
    pb_latlon = [82.19_dp, -60.00_dp]

    ! The two endpoints in x,y
    call oblique_sg_projection( pa_latlon(2), pa_latlon(1), refgeo%grid%lambda_M, refgeo%grid%phi_M, refgeo%grid%alpha_stereo, xa, ya)
    call oblique_sg_projection( pb_latlon(2), pb_latlon(1), refgeo%grid%lambda_M, refgeo%grid%phi_M, refgeo%grid%alpha_stereo, xb, yb)

    pa = [xa,ya]
    pb = [xb,yb]

    do i = 1, refgeo%grid%nx

      yl_ab = pa(2) + (refgeo%grid%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))

      do j = 1, refgeo%grid%ny

        if (refgeo%grid%y(j) > pa(2) .and. refgeo%grid%y(j) > yl_ab .and. refgeo%grid%x(i) < pb(1)) then

          refgeo%Hi_grid( i,j) = 0._dp

        end if

      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_mask_noice_GRL_remove_Ellesmere_grid

end module reference_fields_module