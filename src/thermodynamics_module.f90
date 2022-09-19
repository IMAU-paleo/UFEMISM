module thermodynamics_module
  ! All the routines for calculating the englacial temperature profile

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module,       only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module,          only : CC, T0, ice_density, pi, sec_per_year, R_gas, grav
  use parallel_module,            only : par, sync, ierr, cerr, partition_list
  use utilities_module,           only : check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                         check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                         vertical_average, interpolate_ocean_depth, tridiagonal_solve
  use netcdf_module,              only : debug
  use data_types_module,          only : type_mesh, type_ice_model, type_remapping_mesh_mesh, &
                                         type_climate_snapshot_regional, type_ocean_snapshot_regional, &
                                         type_SMB_model
  use mesh_operators_module,      only : apply_Neumann_BC_direct_a_3D, ddx_a_to_a_2D, ddy_a_to_a_2D, &
                                         ddx_a_to_b_3D, ddy_a_to_b_3D
  use mesh_help_functions_module, only : CROSS2, find_containing_vertex
  use mesh_mapping_module,        only : remap_field_dp_3D
  use mpi_module,                 only : allgather_array
  use reallocate_mod,             only : reallocate_bounds
  use zeta_module,                only : calculate_zeta_derivatives, p_zeta

  implicit none

contains

! ===== Main routines =====
! =========================

  subroutine run_thermo_model( mesh, ice, climate, ocean, SMB, time, do_solve_heat_equation)
    ! Run the thermodynamics model. If so specified, solve the heat equation;
    ! if not, only prescribe a vertically uniform temperature to newly ice-covered grid cells.

    implicit none

    ! In/output variables
    type(type_mesh),                      intent(in)    :: mesh
    type(type_ice_model),                 intent(inout) :: ice
    type(type_climate_snapshot_regional), intent(in)    :: climate
    type(type_ocean_snapshot_regional),   intent(in)    :: ocean
    type(type_SMB_model),                 intent(in)    :: SMB
    real(dp),                             intent(in)    :: time
    logical,                              intent(in)    :: do_solve_heat_equation

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'run_thermo_model'
    integer                                             :: vi, vvi, vj
    logical                                             :: found_source_neighbour
    integer                                             :: n_source_neighbours
    real(dp), dimension(C%nz)                           :: Ti_source_neighbours
    real(dp)                                            :: T_surf_annual
    real(dp), dimension(:,:), allocatable               :: Ti_a

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate full temperature array
    allocate(Ti_a(mesh%nV,C%nz))

    ! Get data from each process
    Ti_a(mesh%vi1:mesh%vi2,:) = ice%Ti_a

    ! Gather the temperature across all processes
    call allgather_array(Ti_a)

    ! === Pick a method ===
    ! =====================

    select case (C%choice_thermo_model)

    case ('none')
      ! Do nothing: choice_ice_rheology_model should be set to "uniform"!

    case ('3D_heat_equation')
      ! Solve the 3-D heat equation

      ! NOTE: solved asynchronously from the ice dynamical equations.
      !       Since newly ice-covered pixels won't have a temperature assigned
      !       until the heat equation is solved again, treat these separately every time step.

      ! Prescribe a simple temperature profile to newly ice-covered grid cells.
      do vi = mesh%vi1, mesh%vi2

        if (ice%mask_ice_a( vi) == 1) then

          if (ice%mask_ice_a_prev( vi) == 0) then
            ! This grid cell is newly ice-covered
            ! If one of its neighbours was already ice-covered, assume the temperature
            ! profile here is equal to the profile from the upstream neighbour (due to advection).
            ! If no neighbours were ice-covered, the new ice must come from accumulation;
            ! just set a simple linear profile instead.

            found_source_neighbour = .false.
            Ti_source_neighbours   = 0._dp
            n_source_neighbours    = 0
            do vvi = 1, mesh%nC( vi)

              vj = mesh%C( vi,vvi)

              if (ice%mask_ice_a_prev( vj) == 1) then
                found_source_neighbour = .true.
                n_source_neighbours    = n_source_neighbours  + 1
                Ti_source_neighbours   = Ti_source_neighbours + Ti_a( vj,:)
              end if

            end do

            if (found_source_neighbour) then
              ! Ice probably was advected from neighbouring grid cells; copy temperature profile from there

              Ti_source_neighbours = Ti_source_neighbours / real( n_source_neighbours,dp)
              ice%Ti_a( vi,:) = Ti_source_neighbours

            else
              ! Ice probably came from surface accumulation; set temperature profile to annual mean surface temperature

              T_surf_annual = min( sum( climate%T2m( vi,:)) / 12._dp, T0)

              ice%Ti_a( vi,:) = T_surf_annual

            end if

          else
            ! This grid cell was already ice-covered in the previous time step, no need to do anything
          end if ! (ice%mask_ice_a_prev( j,i) == 0)

        else
          ! This pixel is ice-free; set temperature profile to zero

          ice%Ti_a( vi,:) = 0._dp

        end if ! (ice%mask_ice_a( vi) == 1)

      end do ! vi = mesh%vi1, mesh%vi2

      deallocate(Ti_a)

      ! Calculate various physical terms
      call calc_heat_capacity(          mesh, ice)
      call calc_thermal_conductivity(   mesh, ice)
      call calc_pressure_melting_point( mesh, ice)

      ! If so specified, solve the heat equation
      if (do_solve_heat_equation) then
        call solve_3D_heat_equation( mesh, ice, climate, ocean, SMB)
      end if

      ! Safety
      call check_for_NaN_dp_2D( ice%Ti_a, 'ice%Ti_a')

    case default
      ! Unknown case
      call crash('unknown choice_thermo_model "' // &
                  trim( C%choice_thermo_model) // '"!')

    end select

    ! Calculate the ice flow factor for the new temperature solution
    call calc_ice_rheology( mesh, ice, time)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_thermo_model

  subroutine initialise_thermo_model( mesh, ice, climate, ocean, SMB, region_name)
    ! Initialise the thermodynamics

    implicit none

    ! In/output variables
    type(type_mesh),                      intent(in)    :: mesh
    type(type_ice_model),                 intent(inout) :: ice
    type(type_climate_snapshot_regional), intent(in)    :: climate
    type(type_ocean_snapshot_regional),   intent(in)    :: ocean
    type(type_SMB_model),                 intent(in)    :: SMB
    character(len=3),                     intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'initialise_thermo_model'

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise ice temperatures
    call initialise_ice_temperature( mesh, ice, climate, ocean, SMB, region_name)

    if (par%master) then
      write (*,"(3A)") '  Initialising ice rheology "', &
                          trim(C%choice_ice_rheology), '"...'
    end if
    call sync

    ! Initialise ice rheology
    call calc_ice_rheology( mesh, ice, C%start_time_of_run)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_thermo_model

! ===== 3D heat-equation solver =====
! ===================================

  subroutine solve_3D_heat_equation( mesh, ice, climate, ocean, SMB)

    implicit none

    ! In/output variables
    type(type_mesh),                      intent(in)    :: mesh
    type(type_ice_model),                 intent(inout) :: ice
    type(type_climate_snapshot_regional), intent(in)    :: climate
    type(type_ocean_snapshot_regional),   intent(in)    :: ocean
    type(type_SMB_model),                 intent(in)    :: SMB

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'solve_3D_heat_equation'
    integer                                             :: vi, k, n_unstable
    real(dp), dimension(:,:), allocatable               :: u_times_dT_dx_upwind_a,  v_times_dT_dy_upwind_a
    real(dp), dimension(:,:), allocatable               :: Ti_new
    real(dp), dimension(:  ), allocatable               :: T_ocean_at_shelf_base
    integer,  dimension(:  ), allocatable               :: is_unstable
    real(dp)                                            :: f1, f2, f3, depth
    real(dp), dimension(2:C%nz)                         :: alpha
    real(dp), dimension(C%nz)                           :: beta
    real(dp), dimension(C%nz-1)                         :: gamma
    real(dp), dimension(C%nz)                           :: delta
    logical                                             :: hasnan

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( u_times_dT_dx_upwind_a (mesh%nTri, C%nz) )
    allocate( v_times_dT_dy_upwind_a (mesh%nTri, C%nz) )
    allocate( is_unstable            (mesh%nV        ) )
    allocate( Ti_new                 (mesh%nV,   C%nz) )
    allocate( T_ocean_at_shelf_base  (mesh%nV        ) )

    ! === Derivatives ===
    ! ===================

    ! Calculate upwind heat flux
    call calc_upwind_heat_flux_derivatives( mesh, ice, u_times_dT_dx_upwind_a, v_times_dT_dy_upwind_a)

    ! Calculate zeta derivatives required for solving the heat equation
    call calculate_zeta_derivatives( mesh, ice)

    ! === Heating terms ===
    ! =====================

    ! Calculate heating terms
    call calc_internal_heating( mesh, ice)
    call calc_frictional_heating( mesh, ice)

    ! === Surface temperatures ===
    ! ============================

    ! Set ice surface temperature equal to annual mean 2m air temperature
    do vi = mesh%vi1, mesh%vi2
      ice%Ti_a( vi,1) = min( T0, sum( climate%T2m( vi,:)) / 12._dp)
    end do

    ! === Ice shelf basal temperatures ===
    ! ====================================

    if (C%choice_ocean_model == 'none') then
      ! No ocean data available; use local pressure melting point of ice

      do vi = mesh%vi1, mesh%vi2
        if (ice%mask_shelf_a( vi) == 1) then
          T_ocean_at_shelf_base( vi) = ice%Ti_pmp_a( vi,C%nz)
        else
          T_ocean_at_shelf_base( vi) = T0
        end if
      end do

    else
      ! Calculate shelf base temperature from ocean data

      do vi = mesh%vi1, mesh%vi2

        if (ice%mask_shelf_a( vi) == 1) then
          depth = max( 0.1_dp, ice%Hi_a( vi) - ice%Hs_a( vi))   ! Depth is positive when below the sea surface!
          call interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( vi,:), depth, T_ocean_at_shelf_base( vi))
        else
          T_ocean_at_shelf_base( vi) = 0._dp
        end if

        ! NOTE: ocean data gives temperature in Celsius, thermodynamics wants Kelvin!
        T_ocean_at_shelf_base( vi) = T_ocean_at_shelf_base( vi) + T0

      end do

    end if

    ! Solve the heat equation for all vertices
    is_unstable( mesh%vi1:mesh%vi2) = 0
    n_unstable                      = 0

    do vi = mesh%vi1, mesh%vi2

      ! Skip ice-free vertices
      if (ice%mask_ice_a( vi) == 0) then
        cycle
      end if

      ! Ice surface boundary condition
      beta(  1) = 1._dp
      gamma( 1) = 0._dp
      delta( 1) = ice%Ti_a( vi,1)

      ! Loop over the whole vertical domain but not the surface (k=1) and the bottom (k=NZ):
      do k = 2, C%nz-1

        f1 = (ice%Ki_a( vi,k) * ice%dzeta_dz_a( vi)**2) / (ice_density * ice%Cpi_a( vi,k))

        f2 = ice%dzeta_dt_a( vi,k) + ice%dzeta_dx_a( vi,k) * ice%u_3D_a( vi,k) + ice%dzeta_dy_a( vi,k) * ice%v_3D_a( vi,k) + ice%dzeta_dz_a( vi) * ice%w_3D_a( vi,k)

        f3 = ice%internal_heating_a( vi,k) + (u_times_dT_dx_upwind_a( vi,k) + v_times_dT_dy_upwind_a( vi,k)) - ice%Ti_a( vi,k) / C%dt_thermo

        alpha(k) = f1 * p_zeta%a_zetazeta(k) - f2 * p_zeta%a_zeta(k)
        beta (k) = f1 * p_zeta%b_zetazeta(k) - f2 * p_zeta%b_zeta(k) - 1._dp / C%dt_thermo
        gamma(k) = f1 * p_zeta%c_zetazeta(k) - f2 * p_zeta%c_zeta(k)
        delta(k) = f3

      end do ! k = 2, C%nz-1

      ! Boundary conditions at the surface: set ice temperature equal to annual mean surface temperature
      beta(  1) = 1._dp
      gamma( 1) = 0._dp
      delta( 1) =  min( T0, sum( climate%T2m( vi,:)) / 12._dp)

      ! Boundary conditions at the base
      if (ice%mask_shelf_a( vi) == 1) then
        ! Set ice bottom temperature equal to seawater temperature (limited to the PMP)
        alpha( C%nz) = 0._dp
        beta ( C%nz) = 1._dp
        delta( C%nz) = min( T0, min( ice%Ti_pmp_a( vi,C%nz), T_ocean_at_shelf_base( vi) ))
      else
        if (ice%Ti_a( vi,C%nz) >= ice%Ti_pmp_a( vi,C%nz)) then
          ! Ice is already at/above pressure melting point; set temperature equal to PMP
          alpha( C%nz) = 0._dp
          beta ( C%nz) = 1._dp
          delta( C%nz) = ice%Ti_pmp_a( vi,C%nz)
        else
          ! Set a Neumann BC so the temperature gradient at the base is equal to basal heating rate (= geothermal + friction)
          alpha( C%nz) = 1._dp
          beta ( C%nz) = -1._dp
          delta( C%nz) = (C%zeta(C%nz) - C%zeta(C%nz-1)) * (ice%GHF_a( vi) + ice%frictional_heating_a( vi)) / (ice%dzeta_dz_a( vi) * ice%Ki_a( vi,C%nz))
        end if
      end if

      ! Solve the tridiagonal matrix equation representing the heat equation for this grid cell
      Ti_new( vi,:) = tridiagonal_solve( alpha, beta, gamma, delta)

      ! Make sure ice temperature doesn't exceed pressure melting point
      do k = 1, C%nz-1
        Ti_new( vi,k) = min( Ti_new( vi,k), ice%Ti_pmp_a( vi,k))
      end do

      if (Ti_new( vi,C%nz) >= ice%Ti_pmp_a( vi,C%nz)) then
        Ti_new( vi,C%nz) = min( ice%Ti_pmp_a( vi,C%nz), ice%Ti_a( vi,C%nz-1) - (C%zeta(C%nz) - C%zeta(C%nz-1)) * &
          (ice%GHF_a( vi) + ice%frictional_heating_a( vi)) / (ice%dzeta_dz_a( vi) * ice%Ki_a( vi,C%nz)))
      end if

      ! Mark temperatures below 150 K or NaN as unstable, to be replaced with the Robin solution.
      hasnan = .false.
      do k = 1, C%nz
        if (Ti_new( vi,k) /= Ti_new( vi,k)) then
          hasnan = .true.
        end if
      end do
      if (minval(Ti_new( vi,:)) < 150._dp .or. hasnan) then
        is_unstable( vi) = 1
        n_unstable  = n_unstable + 1
        write(0,*) '3D heat equation - instability detected; Hi = ', ice%Hi_a( vi), ', dHi_dt = ', ice%dHi_dt_a( vi)
      end if

    end do ! vi = mesh%vi1, mesh%vi2

    ! Gather the new temperature across all processes
    call allgather_array(Ti_new)
    ! Apply Neumann boundary conditions to the temperature field
    call apply_Neumann_BC_direct_a_3D( mesh, Ti_new)

    ! Get total number of instabilities
    call MPI_ALLREDUCE( MPI_IN_PLACE, n_unstable, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Cope with instability
    if (n_unstable < ceiling( real( mesh%nV, dp) / 100._dp)) then
      ! Instability is limited to an acceptably small number (< 1%) of grid cells;
      ! replace the temperature profile in those cells with the Robin solution

      do vi = mesh%vi1, mesh%vi2
        if (is_unstable( vi) == 1) then
          call replace_Ti_with_robin_solution( ice, climate, ocean, SMB, Ti_new(vi,:), vi)
        end if
      end do

    else
      ! An unacceptably large number of grid cells was unstable; throw an error.

      ! ! Save to debug file
      ! debug%dp_2D_a_01  = ice%Hi_a
      ! debug%int_2D_a_01 = ice%mask_ice_a
      ! debug%dp_3D_a_01  = ice%Ti_a
      ! debug%dp_3D_a_02  = Ti_new
      ! call write_to_debug_file

      ! Stop
      call crash('heat equation solver unstable for more than 1% of vertices!')

    end if

    ! Move the new temperature field to the ice data structure
    ice%Ti_a( mesh%vi1:mesh%vi2,:) = Ti_new( mesh%vi1:mesh%vi2,:)

    ! Clean up after yourself
    deallocate( u_times_dT_dx_upwind_a)
    deallocate( v_times_dT_dy_upwind_a)
    deallocate( Ti_new                )
    deallocate( is_unstable           )
    deallocate( T_ocean_at_shelf_base )

    ! Safety
    call check_for_NaN_dp_2D( ice%Ti_a, 'ice%Ti_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_3D_heat_equation

  subroutine calc_upwind_heat_flux_derivatives( mesh, ice, u_times_dT_dx_upwind_a, v_times_dT_dy_upwind_a)
    ! Calculate upwind heat flux derivatives at vertex vi, vertical layer k

    implicit none

    ! In/output variables:
    type(type_mesh),            intent(in)  :: mesh
    type(type_ice_model),       intent(in)  :: ice
    real(dp), dimension(:,:),   intent(out) :: u_times_dT_dx_upwind_a, v_times_dT_dy_upwind_a

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calc_upwind_heat_flux_derivatives'
    real(dp), dimension(:,:), allocatable   :: dTi_dx_3D_b, dTi_dy_3D_b, u_3D_b, v_3D_b
    integer                                 :: vi, k, vti, ti, n1, n2, n3, vib, vic, ti_upwind
    real(dp), dimension(2)                  :: u_upwind, ab, ac

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Temperature gradients ===
    ! =============================

    ! Allocate shared memory
    allocate( dTi_dx_3D_b (1:mesh%nTri, C%nz) )
    allocate( dTi_dy_3D_b (1:mesh%nTri, C%nz) )

    ! Calculate temperature gradients on the b-grid
    call ddx_a_to_b_3D( mesh, ice%Ti_a, dTi_dx_3D_b(mesh%ti1:mesh%ti2,:))
    call ddy_a_to_b_3D( mesh, ice%Ti_a, dTi_dy_3D_b(mesh%ti1:mesh%ti2,:))

    ! Gather the temperature gradients across all processes
    call allgather_array(dTi_dx_3D_b)
    call allgather_array(dTi_dy_3D_b)

    ! === Full velocity array ===
    ! ===========================

    ! Allocate full velocity arrays
    allocate(u_3D_b (1:mesh%nTri, C%nz) )
    allocate(v_3D_b (1:mesh%nTri, C%nz) )

    ! Get data from each process
    u_3D_b( mesh%ti1:mesh%ti2,:) = ice%u_3D_b
    v_3D_b( mesh%ti1:mesh%ti2,:) = ice%v_3D_b

    ! Gather velocities across all processes
    call allgather_array(u_3D_b)
    call allgather_array(v_3D_b)

    ! === Upwind derivatives ===
    ! ==========================

    ! Initialise
    u_times_dT_dx_upwind_a( mesh%vi1:mesh%vi2,:) = 0._dp
    v_times_dT_dy_upwind_a( mesh%vi1:mesh%vi2,:) = 0._dp

    do vi = mesh%vi1, mesh%vi2

      if (ice%mask_ice_a( vi) == 1) then

        ! The upwind velocity vector
        u_upwind = [-ice%u_vav_a( vi), -ice%v_vav_a( vi)]

        ! Find the upwind triangle
        ti_upwind = 0
        do vti = 1, mesh%niTri( vi)

          ! Triangle ti is spanned counter-clockwise by vertices [vi,vib,vic]
          ti  = mesh%iTri( vi,vti)
          vib = 0
          vic = 0
          do n1 = 1, 3
            n2 = n1 + 1
            if (n2 == 4) n2 = 1
            n3 = n2 + 1
            if (n3 == 4) n3 = 1

            if (mesh%Tri( ti,n1) == vi) then
              vib = mesh%Tri( ti,n2)
              vic = mesh%Tri( ti,n3)
              exit
            end if
          end do

          ! Check if the upwind velocity vector points into this triangle
          ab = mesh%V( vib,:) - mesh%V( vi,:)
          ac = mesh%V( vic,:) - mesh%V( vi,:)

          if (CROSS2( ab, u_upwind) >= 0._dp .and. CROSS2( u_upwind, ac) >= 0._dp) then
            ti_upwind = ti
            exit
          end if

        end do ! iti = 1, mesh%niTri( vi)

        ! Check if triangle was found
        if (ti_upwind == 0) then
          ! No triangle found

          ! Check whether it matters
          if (mesh%edge_index( vi) > 0) THEN
            ! Border vertex: just forget about it
            u_times_dT_dx_upwind_a( vi,:) = 0._dp
            v_times_dT_dy_upwind_a( vi,:) = 0._dp
          else
            ! Interior vertex: uh oh...
            call crash('couldnt find upwind triangle!')
          end if

        else
          ! Triangle found: calculate u * dT/dx, v * dT/dy

          do k = 1, C%nz
            u_times_dT_dx_upwind_a( vi,k) = u_3D_b( ti_upwind,k) * dTi_dx_3D_b( ti_upwind,k)
            v_times_dT_dy_upwind_a( vi,k) = v_3D_b( ti_upwind,k) * dTi_dy_3D_b( ti_upwind,k)
          end do

        end if

      end if ! (ice%mask_ice_a( vi) == 1)

    end do ! vi = mesh%vi1, mesh%vi2

    ! Clean up after yourself
    deallocate( dTi_dx_3D_b)
    deallocate( dTi_dy_3D_b)
    deallocate( u_3D_b)
    deallocate( v_3D_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_upwind_heat_flux_derivatives

! ===== Initial temperatures =====
! ================================

  subroutine initialise_ice_temperature( mesh, ice, climate, ocean, SMB, region_name)
    ! Initialise the englacial ice temperature at the start of a simulation

    implicit none

    ! In/output variables
    type(type_mesh),                      intent(in)    :: mesh
    type(type_ice_model),                 intent(inout) :: ice
    type(type_climate_snapshot_regional), intent(in)    :: climate
    type(type_ocean_snapshot_regional),   intent(in)    :: ocean
    type(type_SMB_model),                 intent(in)    :: SMB
    character(len=3),                     intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'initialise_ice_temperature'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write (*,"(3A)") '  Initialising ice temperature profile "', &
                          trim(C%choice_initial_ice_temperature), '"...'
    end if
    call sync

    ! === Temperature profile ===
    ! ===========================

    select case (C%choice_initial_ice_temperature)

    case ('uniform')
      ! Simple uniform temperature
      call initialise_ice_temperature_uniform( mesh, ice)
    case ('linear')
      ! Simple linear temperature profile
      call initialise_ice_temperature_linear( mesh, ice, climate)
    case ('Robin')
      ! Initialise with the Robin solution
      call initialise_ice_temperature_Robin( mesh, ice, climate, ocean, SMB)
    case ('restart')
      ! Initialise with the temperature field from a restart file
      call crash('choice_initial_ice_temperature "' // &
           trim( C%choice_initial_ice_temperature) // &
           '" not implemented yet!')
    case default
      ! Unknown case
      call crash('unknown choice_initial_ice_temperature "' // &
           trim( C%choice_initial_ice_temperature) // '"!')
    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ice_temperature

  subroutine initialise_ice_temperature_uniform( mesh, ice)
    ! Simple uniform temperature

    implicit none

    ! In/output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_ice_temperature_uniform'
    integer                             :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      if (ice%Hi_a( vi) > 0._dp) then
        ice%Ti_a( vi,:) = C%uniform_ice_temperature
      else
        ice%Ti_a( vi,:) = 0._dp
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ice_temperature_uniform

  subroutine initialise_ice_temperature_linear( mesh, ice, climate)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Simple linear temperature profile

    implicit none

    ! In/output variables
    type(type_mesh),                      intent(in)    :: mesh
    type(type_ice_model),                 intent(inout) :: ice
    type(type_climate_snapshot_regional), intent(in)    :: climate

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'initialise_ice_temperature_linear'
    integer                                             :: vi
    real(dp)                                            :: T_surf_annual, T_PMP_base

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      if (ice%Hi_a( vi) > 0._dp) then
        T_surf_annual = min( sum( climate%T2m( vi,:)) / 12._dp, T0)
        T_PMP_base    = T0 - (ice%Hi_a( vi) * 8.7E-04_dp)
        ice%Ti_a( vi,:) = T_surf_annual - C%zeta * (T_surf_annual - T_PMP_base)
      else
        ice%Ti_a( vi,:) = 0._dp
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ice_temperature_linear

  subroutine initialise_ice_temperature_Robin( mesh, ice, climate, ocean, SMB)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Initialise with the Robin solution

    implicit none

    ! In/output variables
    type(type_mesh),                      intent(in)    :: mesh
    type(type_ice_model),                 intent(inout) :: ice
    type(type_climate_snapshot_regional), intent(in)    :: climate
    type(type_ocean_snapshot_regional),   intent(in)    :: ocean
    type(type_SMB_model),                 intent(in)    :: SMB

    ! Local variables:
    character(len=256), parameter                       :: routine_name = 'initialise_ice_temperature_Robin'
    integer                                             :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate Ti_pmp
    call calc_pressure_melting_point( mesh, ice)

    ! Initialise with the Robin solution
    do vi = mesh%vi1, mesh%vi2
      call replace_Ti_with_robin_solution( ice, climate, ocean, SMB, ice%Ti_a(vi,:), vi)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ice_temperature_Robin

! ===== Rheology =====
! ====================

  subroutine calc_ice_rheology( mesh, ice, time)
    ! Calculate the flow factor A in Glen's flow law

    implicit none

    ! In/output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice
    real(dp),             intent(in)    :: time

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_ice_rheology'
    integer                             :: vi,k
    real(dp), dimension(C%nZ)           :: prof
    real(dp), parameter                 :: A_low_temp  = 1.14E-05_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    real(dp), parameter                 :: A_high_temp = 5.47E+10_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    real(dp), parameter                 :: Q_low_temp  = 6.0E+04_dp    ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    real(dp), parameter                 :: Q_high_temp = 13.9E+04_dp   ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    real(dp)                            :: A_flow_MISMIP

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Flow factor ===
    ! ===================

    select case (C%choice_ice_rheology)

      case ('uniform')
        ! Apply a uniform value for the ice flow factor
        ice%A_flow_3D_a( mesh%vi1:mesh%vi2,:) = C%uniform_flow_factor

      case ('Huybrechts1992')
        ! Calculate the ice flow factor as a function of the ice temperature according to the Arrhenius relationship (Huybrechts, 1992)

        do vi = mesh%vi1, mesh%vi2

          do k = 1, C%nz

            if (ice%mask_ice_a( vi) == 1) then

              if (ice%Ti_a( vi,k) < 263.15_dp) then
                ice%A_flow_3D_a( vi,k) = A_low_temp  * exp(-Q_low_temp  / (R_gas * ice%Ti_a( vi,k)))
              else
                ice%A_flow_3D_a( vi,k) = A_high_temp * exp(-Q_high_temp / (R_gas * ice%Ti_a( vi,k)))
              end if

            else

              if (C%choice_ice_margin == 'BC') then
                ice%A_flow_3D_a( vi,k) = 0._dp
              elseif (C%choice_ice_margin == 'infinite_slab') then
                ! In the "infinite slab" case, calculate effective viscosity everywhere
                ! (even when there's technically no ice present)
                ice%A_flow_3D_a( vi,k) = A_low_temp  * exp(-Q_low_temp  / (R_gas * 263.15_dp))
              else
                call crash('unknown choice_ice_margin "' // trim( C%choice_ice_margin) // '"!')
              end if

            end if

          end do ! k = 1, C%nz

        end do ! vi = mesh%vi1, mesh%vi2

      case default
        ! Unknown case
        call crash('unknown choice_ice_margin "' // trim( C%choice_ice_margin) // '"!')

    end select

    ! === Enhancement factors ===
    ! ===========================

    do vi = mesh%vi1, mesh%vi2
      if     (ice%mask_sheet_a( vi) == 1) then
        ice%A_flow_3D_a( vi,:) = ice%A_flow_3D_a( vi,:) * C%m_enh_sheet
      elseif (ice%mask_shelf_a( vi) == 1) then
        ice%A_flow_3D_a( vi,:) = ice%A_flow_3D_a( vi,:) * C%m_enh_shelf
      end if
    end do

    ! === Vertical average ===
    ! ========================

    do vi = mesh%vi1, mesh%vi2
      prof = ice%A_flow_3D_a( vi,:)
      call vertical_average( prof, ice%A_flow_vav_a( vi))
    end do

    ! === Finalisation ===
    ! ====================

    ! Safety
    call check_for_NaN_dp_2D( ice%A_flow_3D_a , 'ice%A_flow_3D_a' )
    call check_for_NaN_dp_1D( ice%A_flow_vav_a, 'ice%A_flow_vav_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_rheology

! ===== Tools =====
! =================

  subroutine replace_Ti_with_robin_solution( ice, climate, ocean, SMB, Ti, vi)
    ! This function calculates for one horizontal grid point the temperature profiles
    ! using the surface temperature and the geothermal heat flux as boundary conditions.
    ! See Robin solution in: Cuffey & Paterson 2010, 4th ed, chapter 9, eq. (9.13) - (9.22).

    implicit none

    ! In/output variables:
    type(type_ice_model),                 intent(in)    :: ice
    type(type_climate_snapshot_regional), intent(in)    :: climate
    type(type_ocean_snapshot_regional),   intent(in)    :: ocean
    type(type_SMB_model),                 intent(in)    :: SMB
    real(dp), dimension(:),               intent(inout) :: Ti
    integer,                              intent(in)    :: vi

    ! Local variables:
    integer                                             :: k
    real(dp)                                            :: Ts
    real(dp)                                            :: thermal_length_scale
    real(dp)                                            :: distance_above_bed
    real(dp)                                            :: erf1
    real(dp)                                            :: erf2

    real(dp)                                            :: thermal_conductivity_robin
    real(dp)                                            :: thermal_diffusivity_robin
    real(dp)                                            :: bottom_temperature_gradient_robin

    real(dp), parameter                                 :: kappa_0_ice_conductivity     = 9.828_dp                   ! The linear constant in the thermal conductivity of ice [J m^-1 K^-1 s^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    real(dp), parameter                                 :: kappa_e_ice_conductivity     = 0.0057_dp                  ! The exponent constant in the thermal conductivity of ice [K^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    real(dp), parameter                                 :: c_0_specific_heat            = 2127.5_dp                  ! The constant in the specific heat capacity of ice [J kg^-1 K^-1], see equation (12.5), Zwinger (2007), Cuffey & Paterson (2010, p. 400)

    real(dp)                                            :: depth
    real(dp)                                            :: T_ocean_at_shelf_base

    thermal_conductivity_robin        = kappa_0_ice_conductivity * sec_per_year * exp(-kappa_e_ice_conductivity * T0) ! Thermal conductivity            [J m^-1 K^-1 y^-1]
    thermal_diffusivity_robin         = thermal_conductivity_robin / (ice_density * c_0_specific_heat)                ! Thermal diffusivity             [m^2 y^-1]
    bottom_temperature_gradient_robin = - ice%GHF_a( vi) / thermal_conductivity_robin                                 ! Temperature gradient at bedrock

    Ts = min( T0, sum(climate%T2m( vi,:)) / 12._dp)

    if (ice%mask_sheet_a( vi) == 1 ) then
      ! Ice sheet vertex

      if (SMB%SMB_year( vi) > 0._dp) then
        ! The Robin solution can be used to estimate the subsurface temperature profile in an accumulation area

        thermal_length_scale = sqrt(2._dp * thermal_diffusivity_robin * ice%Hi_a( vi) / SMB%SMB_year( vi))
        do k = 1, C%nz
          distance_above_bed = (1._dp - C%zeta(k)) * ice%Hi_a( vi)
          erf1 = erf( distance_above_bed / thermal_length_scale)
          erf2 = erf( ice%Hi_a( vi) / thermal_length_scale)
          Ti( k) = Ts + sqrt(pi) / 2._dp * thermal_length_scale * bottom_temperature_gradient_robin * (erf1 - erf2)
        end do

      else
        ! Ablation area: use linear temperature profile from Ts to (offset below) T_pmp
        Ti( :) = Ts + ((T0 - CC * ice%Hi_a( vi)) - Ts) * C%zeta(:)

      end if

    elseif( ice%mask_shelf_a(vi) == 1) then
      ! Ice shelf vertex

      if (C%choice_ocean_model == 'none') then
        ! No ocean data available; use local pressure melting point of ice

        T_ocean_at_shelf_base = ice%Ti_pmp_a( vi,C%nz)

      else
        ! Calculate shelf base temperature from ocean data

        depth = max( 0.1_dp, ice%Hi_a( vi) - ice%Hs_a( vi))   ! Depth is positive when below the sea surface!
        call interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( vi,:), depth, T_ocean_at_shelf_base)

      end if

      Ti( :) = Ts + C%zeta(:) * (T0 + T_ocean_at_shelf_base - Ts)

    else
      ! Ice-free vertex

      ! Simply use Ts
      Ti( :) = Ts

    end if

    ! Correct all temperatures above T_pmp:
    do k = 1, C%nz
      Ti( k) = min( Ti( k), ice%Ti_pmp_a( vi,k))
    end do

  end subroutine replace_Ti_with_robin_solution

  subroutine calc_pressure_melting_point( mesh, ice)
    ! Calculate the pressure melting point of the ice according to Huybrechts (1992)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),      INTENT(IN)    :: mesh
    TYPE(type_ice_model), INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER       :: routine_name = 'calc_pressure_melting_point'
    INTEGER                             :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      ice%Ti_pmp_a( vi,:) = T0 - CC * ice%Hi_a( vi) * C%zeta
    END DO

    ! Safety
    call check_for_NaN_dp_2D( ice%Ti_pmp_a, 'ice%Ti_pmp_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_pressure_melting_point

  subroutine calc_heat_capacity( mesh, ice)
    ! Calculate the heat capacity of the ice

    implicit none

    ! In/output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_heat_capacity'
    integer                             :: vi

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_ice_heat_capacity)

      case ('uniform')
        ! Apply a uniform value for the heat capacity
        ice%Cpi_a( mesh%vi1:mesh%vi2,:) = C%uniform_ice_heat_capacity

      case ('Pounder1965')
        ! Calculate the heat capacity of ice according to Pounder: The Physics of Ice (1965)
        do vi = mesh%vi1, mesh%vi2
          ice%Cpi_a( vi,:) = 2115.3_dp + 7.79293_dp * (ice%Ti_a( vi,:) - T0)
        end do

      case default
        ! Unknown case
        call crash('unknown choice_ice_heat_capacity "' // &
                    trim( C%choice_ice_heat_capacity) // '"!')

    end select

    ! Safety
    call check_for_NaN_dp_2D( ice%Cpi_a, 'ice%Cpi_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_heat_capacity

  subroutine calc_thermal_conductivity( mesh, ice)
    ! Calculate the thermal conductivity of the ice

    implicit none

    ! In/output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_thermal_conductivity'
    integer                             :: vi

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_ice_thermal_conductivity)

      case ('uniform')
        ! Apply a uniform value for the thermal conductivity
        ice%Ki_a( mesh%vi1:mesh%vi2,:) = C%uniform_ice_thermal_conductivity

      case ('Ritz1987')
        ! Calculate the thermal conductivity of ice according to Ritz (1987)
        do vi = mesh%vi1, mesh%vi2
          ice%Ki_a( vi,:) = 3.101E+08_dp * exp(-0.0057_dp * ice%Ti_a( vi,:))
        end do

      case default
        ! Unknown case
        call crash('unknown choice_ice_thermal_conductivity "' // &
                    trim( C%choice_ice_thermal_conductivity) // '"!')

    end select

    ! Safety
    call check_for_NaN_dp_2D( ice%Ki_a, 'ice%Ki_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_thermal_conductivity

  subroutine calc_internal_heating( mesh, ice)
    ! Calculate internal heating due to deformation

    implicit none

    ! In- and output variables
    type(type_ice_model), intent(inout) :: ice
    type(type_mesh),      intent(in)    :: mesh

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_internal_heating'
    integer                             :: vi, k
    real(dp), dimension(:), allocatable :: dHs_dx_a, dHs_dy_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( dHs_dx_a (mesh%vi1:mesh%vi2) )
    allocate( dHs_dy_a (mesh%vi1:mesh%vi2) )

    ! Calculate surface slopes
    call ddx_a_to_a_2D( mesh, ice%Hs_a, dHs_dx_a)
    call ddy_a_to_a_2D( mesh, ice%Hs_a, dHs_dy_a)

    ! Calculate internal heating
    do vi = mesh%vi1, mesh%vi2

      ice%internal_heating_a( vi,:) = 0._dp

      if (mesh%edge_index( vi) > 0) cycle ! Skip the domain boundary
      if (ice%mask_ice_a( vi) == 0) cycle ! Skip ice-less elements

      ! Loop over the whole vertical domain but not the surface (k=1) nor the bottom (k=nz):
      do k = 2, C%nz-1
        ice%internal_heating_a( vi,k) = ((- grav * C%zeta(k)) / ice%Cpi_a( vi,k)) * &
                                        ((p_zeta%a_zeta(k) * ice%u_3D_a( vi,k-1) + &
                                          p_zeta%b_zeta(k) * ice%u_3D_a( vi,k) + &
                                          p_zeta%c_zeta(k) * ice%u_3D_a( vi,k+1)) * dHs_dx_a( vi) + &
                                         (p_zeta%a_zeta(k) * ice%v_3D_a( vi,k-1) + &
                                          p_zeta%b_zeta(k) * ice%v_3D_a( vi,k) + &
                                          p_zeta%c_zeta(k) * ice%v_3D_a( vi,k+1)) * dHs_dy_a( vi))
      end do

    end do

    ! Clean up after yourself
    deallocate( dHs_dx_a)
    deallocate( dHs_dy_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_internal_heating

  subroutine calc_frictional_heating( mesh, ice)
    ! Calculate frictional heating at the base due to sliding

    implicit none

    ! In- and output variables
    type(type_ice_model), intent(inout) :: ice
    type(type_mesh),      intent(in)    :: mesh

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_frictional_heating'
    integer                             :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Exception for when no sliding can occur
    if (C%choice_ice_dynamics == 'SIA' .or. &
        C%choice_sliding_law == 'no_sliding') then

      ice%frictional_heating_a( mesh%vi1:mesh%vi2) = 0._dp

    else

      ! Calculate frictional heating
      do vi = mesh%vi1, mesh%vi2
        if (ice%mask_sheet_a( vi) == 1) then
          ice%frictional_heating_a( vi) = ice%beta_a( vi) * (ice%u_base_a( vi)**2 + ice%u_base_a( vi)**2)
        else
          ice%frictional_heating_a( vi) = 0._dp
        end if
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_frictional_heating

! ===== Remapping =====
! =====================

  subroutine remap_ice_temperature( mesh_old, mesh_new, map, ice)
    ! Remap englacial temperature

    ! In/output variables:
    type(type_mesh),                intent(in)    :: mesh_old
    type(type_mesh),                intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh), intent(in)    :: map
    type(type_ice_model),           intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'remap_ice_temperature'
    integer,  dimension(:  ), allocatable         :: mask_ice_a_new
    integer                                       :: vi, vvi, vj
    real(dp), dimension(:,:), allocatable         :: Ti_ext, Ti_ext_local
    integer,  dimension(:  ), allocatable         :: Vmap, Vstack1, Vstack2, Vtemp
    integer                                       :: VstackN1, VstackN2
    integer                                       :: sti, n, it
    real(dp), dimension(C%nz)                     :: Ti_av

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    CALL init_routine( routine_name)

    ! === Ice mask ===
    ! ================

    ! Allocate memory
    allocate( mask_ice_a_new(mesh_new%vi1:mesh_new%vi2))

    do vi = mesh_new%vi1, mesh_new%vi2
      if (ice%Hi_a( vi) > 0._dp) then
        mask_ice_a_new( vi) = 1
      else
        mask_ice_a_new( vi) = 0
      end if
    end do

    ! === Temperature extrapolation ===
    ! =================================

    ! Allocation
    allocate(Ti_ext   ( 1:mesh_old%nV, C%nz) )
    allocate( Vmap    ( 1:mesh_old%nV      ) )
    allocate( Vstack1 ( 1:mesh_old%nV      ) )
    allocate( Vstack2 ( 1:mesh_old%nV      ) )

    ! Initial values
    Ti_ext(mesh_old%vi1:mesh_old%vi2,:) = ice%Ti_a

    ! Gather info from all processes
    call allgather_array(Ti_ext)

    ! Initialise the stack with all ice-free-next-to-ice-covered
    ! vertices (and also initialise the map)
    Vmap     = 0
    Vstack2  = 0
    VstackN2 = 0

    do vi = 1, mesh_old%nV

      if (ice%mask_ice_a( vi) == 1) then

        ! Ice vertex
        Vmap( vi) = 2

      else

        ! Ice-free vertex
        do vvi = 1, mesh_old%nC( vi)

          vj = mesh_old%C( vi,vvi)

          if (ice%mask_ice_a( vj) == 1) then
            ! Vertex vi is ice-free, but adjacent to ice-covered vertex vj
            VMap( vi) = 1
            VstackN2 = VstackN2 + 1
            Vstack2(   VstackN2) = vi
            exit
          end if

        end do

      end if

    end do

    ! Perform a flood-fill-style extrapolation
    it = 0
    do while (VstackN2 > 0)

      it = it + 1

      ! Cycle stacks, efficiently
      call move_alloc(Vstack1, Vtemp)
      call move_alloc(Vstack2, Vstack1)
      VstackN1 = VstackN2
      call move_alloc(Vtemp, Vstack2)
      VstackN2 = 0

      ! Extrapolate temperature values into data-less-next-to-data-filled pixels
      do sti = 1, VstackN1

        vi = Vstack1( sti)

        n     = 0
        Ti_av = 0._dp

        do vvi = 1, mesh_old%nC( vi)

          vj = mesh_old%C( vi,vvi)

          if (VMap( vj) == 2) then
            n     = n     + 1
            Ti_av = Ti_av + Ti_ext( vj,:)
          end if

        end do ! vvi = 1, mesh_old%nC( vi)

        ! Extrapolate temperature by averaging over data-filled neighbours
        Ti_av = Ti_av / real( n,dp)
        Ti_ext( vi,:) = Ti_av

      end do ! sti = 1: VstackN1

      ! Create new stack of data-less-next-to-data-filled pixels
      do sti = 1, VstackN1

        vi = Vstack1( sti)

        ! Mark this pixel as data-filled on the map
        Vmap( vi) = 2

        ! Add its data-less neighbours to the stack
        do vvi = 1, mesh_old%nC( vi)

          vj = mesh_old%C( vi,vvi)

          if (Vmap( vj) == 0) then
            Vmap( vj) = 1
            VstackN2 = VstackN2 + 1
            Vstack2(   VstackN2) = vj
          end if

        end do ! vvi = 1, mesh_old%nC( vi)

      end do ! sti = 1: VstackN1

    end do ! while (VstackN2 > 0)

    ! Clean up after yourself
    deallocate( Vmap   )
    deallocate( Vstack1)
    deallocate( Vstack2)

    ! === Remapping ===
    ! =================

    ! Remap the extrapolated temperature field
    allocate(Ti_ext_local(mesh_old%vi1:mesh_old%vi2,C%nz))
    Ti_ext_local = Ti_ext(mesh_old%vi1:mesh_old%vi2,:)
    deallocate(Ti_ext)

    call remap_field_dp_3D( mesh_old, mesh_new, map, Ti_ext_local, 'cons_2nd_order')

    ! Reallocate ice temperature field
    call reallocate_bounds( ice%Ti_a, mesh_new%vi1, mesh_new%vi2, C%nz )

    ! Copy remapped data only for ice-covered pixels
    do vi = mesh_new%vi1, mesh_new%vi2
      if (mask_ice_a_new( vi) == 1) then
        ice%Ti_a( vi,:) = Ti_ext_local( vi,:)
      end if
    end do

    ! Reallocate mask_ice_a_prev
    call reallocate_bounds( ice%mask_ice_a_prev,1,mesh_new%nV)

    ! Fill it in (needed for the generic temperature update)
    ice%mask_ice_a_prev(mesh_new%vi1:mesh_new%vi2) = mask_ice_a_new
    call allgather_array(ice%mask_ice_a_prev)

    ! === Finalisation ===
    ! ====================

    ! Clean up after yourself
    deallocate( mask_ice_a_new)
    deallocate( Ti_ext_local  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_ice_temperature

end module thermodynamics_module