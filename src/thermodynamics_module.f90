MODULE thermodynamics_module

  ! All the routines for calculating the englacial temperature profile

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE petscksp
  USE configuration_module,                ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                        ONLY: perr, petscmat_checksum
  USE parallel_module,                     ONLY: par, sync, ierr, cerr, partition_list, &
                                                 allocate_shared_int_0D,   allocate_shared_dp_0D, &
                                                 allocate_shared_int_1D,   allocate_shared_dp_1D, &
                                                 allocate_shared_int_2D,   allocate_shared_dp_2D, &
                                                 allocate_shared_int_3D,   allocate_shared_dp_3D, &
                                                 allocate_shared_bool_0D,  allocate_shared_bool_1D, &
                                                 reallocate_shared_int_0D, reallocate_shared_dp_0D, &
                                                 reallocate_shared_int_1D, reallocate_shared_dp_1D, &
                                                 reallocate_shared_int_2D, reallocate_shared_dp_2D, &
                                                 reallocate_shared_int_3D, reallocate_shared_dp_3D, &
                                                 deallocate_shared
  USE utilities_module,                    ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                                 check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                                 checksum_dp_1D, checksum_dp_2D, checksum_dp_3D
  USE netcdf_debug_module,                 ONLY: debug, write_to_debug_file, &
                                                 save_variable_as_netcdf_int_1D, save_variable_as_netcdf_dp_1D, &
                                                 save_variable_as_netcdf_int_2D, save_variable_as_netcdf_dp_2D, &
                                                 save_variable_as_netcdf_int_3D, save_variable_as_netcdf_dp_3D, &
                                                 write_PETSc_matrix_to_NetCDF, write_CSR_matrix_to_NetCDF

  ! Import specific functionality
  USE data_types_module,                   ONLY: type_mesh, type_ice_model, type_BMB_model, type_climate_snapshot_regional, &
                                                 type_SMB_model, type_restart_data
  USE utilities_module,                    ONLY: tridiagonal_solve, vertical_average
  USE mesh_operators_module,               ONLY: ddx_a_to_b_3D, ddy_a_to_b_3D, calc_zeta_gradients
  USE mesh_help_functions_module,          ONLY: CROSS2
  USE mesh_mapping_module,                 ONLY: remap_field_dp_3D
  USE ice_velocity_main_module,            ONLY: calc_vertical_velocities

  IMPLICIT NONE

CONTAINS

! == Run the chosen thermodynamics model
  SUBROUTINE run_thermo_model( mesh, ice, climate, SMB, BMB, time, do_solve_heat_equation)
    ! Run the thermodynamics model. If so specified, solve the heat equation;
    ! if not, only prescribe a vertically uniform temperature to newly ice-covered grid cells.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_snapshot_regional), INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                 INTENT(IN)    :: BMB
    REAL(dp),                             INTENT(IN)    :: time
    LOGICAL,                              INTENT(IN)    :: do_solve_heat_equation

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'run_thermo_model'
    INTEGER                                             :: vi, vvi, vj
    LOGICAL                                             :: found_source_neighbour
    INTEGER                                             ::  n_source_neighbours
    REAL(dp), DIMENSION(C%nz)                           :: Ti_source_neighbours
    REAL(dp)                                            :: T_surf_annual

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_thermo_model == 'none') THEN
      ! No need to do anything
      ! NOTE: choice_ice_rheology should be set to "uniform"!

    ELSEIF (C%choice_thermo_model == '3D_heat_equation') THEN
      ! Solve the 3-D heat equation

      ! NOTE: solved asynchronously from the ice dynamical equations.
      !       Since newly ice-covered pixels won't have a temperature assigned
      !       until the heat equation is solved again, treat these separately every time step.

      ! Prescribe a simple temperature profile to newly ice-covered grid cells.
      DO vi = mesh%vi1, mesh%vi2

        T_surf_annual = SUM( climate%T2m( vi,:)) / REAL( SIZE( climate%T2m,2),dp)

        IF (ice%mask_ice_a( vi) == 1) THEN
          ! This grid cell is now ice-covered

          IF (ice%mask_ice_a_prev( vi) == 0) THEN
            ! This grid cell is newly ice-covered; set temperature profile to annual mean surface temperature

            ice%Ti_a( vi,:) = MIN( T0, T_surf_annual)

          ELSE ! IF (ice%mask_ice_a_prev( vi) == 0) THEN
            ! This grid cell was already ice-covered in the previous time step, no need to do anything
          END IF ! IF (ice%mask_ice_a_prev( vi) == 0) THEN

        ELSE ! IF (ice%mask_ice_a( vi) == 1) THEN
          ! This pixel is ice-free; set temperature profile to annual mean surface temperature

          ice%Ti_a( vi,:) = MIN( T0, T_surf_annual)

        END IF ! IF (ice%mask_ice_a( vi) == 1) THEN

      END DO
      CALL sync

      ! Calculate various physical terms
      CALL calc_heat_capacity(          mesh, ice)
      CALL calc_thermal_conductivity(   mesh, ice)
      CALL calc_pressure_melting_point( mesh, ice)

      ! If so specified, solve the heat equation
      IF (do_solve_heat_equation) CALL solve_3D_heat_equation( mesh, ice, climate, SMB, BMB, C%dt_thermo)

      ! Safety
      CALL check_for_NaN_dp_2D( ice%Ti_a, 'ice%Ti_a')

    ELSE
      CALL crash('unknown choice_thermo_model "' // TRIM( C%choice_thermo_model) // '"!')
    END IF

    ! Calculate the ice flow factor for the new temperature solution
    CALL calc_ice_rheology( mesh, ice, time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_thermo_model

! == Solve the 3-D heat equation
  SUBROUTINE solve_3D_heat_equation( mesh, ice, climate, SMB, BMB, dt)
    ! Solve the three-dimensional heat equation
    !
    ! (See solve_1D_heat_equation for the derivation)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_snapshot_regional), INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                 INTENT(IN)    :: BMB
    REAL(dp),                             INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_3D_heat_equation'
    INTEGER                                            :: vi, k
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  u_times_dTdxp_upwind,  v_times_dTdyp_upwind
    INTEGER                                            :: wu_times_dTdxp_upwind, wv_times_dTdyp_upwind
    REAL(dp), DIMENSION(:    ), POINTER                ::  T_surf_annual,  Q_base_grnd,  T_base_float
    INTEGER                                            :: wT_surf_annual, wQ_base_grnd, wT_base_float
    REAL(dp)                                           :: dt_applied
    INTEGER                                            :: it_dt, it_it_dt
    LOGICAL                                            :: found_stable_solution
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  Ti_tplusdt
    INTEGER                                            :: wTi_tplusdt
    INTEGER,  DIMENSION(:    ), POINTER                ::  is_unstable
    INTEGER                                            :: wis_unstable
    INTEGER                                            :: n_unstable

    REAL(dp), DIMENSION( C%nz)                         :: icecol_Ti                   ! Vertical profile of ice temperature at time t                     [K]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_u                    !   "         "    "  horizontal ice velocity in the x-direction    [m yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_v                    !   "         "    "  horizontal ice velocity in the y-direction    [m yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_w                    !   "         "    "  vertical   ice velocity                       [m yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_u_times_dTdxp_upwind !   "         "    "  u * dT/dxp in the upwind direction            [K yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_v_times_dTdyp_upwind !   "         "    "  u * dT/dxp in the upwind direction            [K yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Ti_pmp               !   "         "    "  pressure melting point temperature            [K]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Ki                   !   "         "    "  thermal conductivity of ice                   [J yr^-1 m^-1 K^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Cpi                  !   "         "    "  specific heat of ice                          [J kg^-1 K^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_dzeta_dx             !   "         "    "  x-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_dzeta_dy             !   "         "    "  y-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_dzeta_dz             !   "         "    "  z-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_dzeta_dt             !   "         "    "  time-derivative of the scaled coordinate zeta [yr^-1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Phi                  !   "         "    "  internal heat production                      [J kg^1 yr^1]
    REAL(dp), DIMENSION( C%nz)                         :: icecol_Ti_tplusdt           ! Vertical profile of ice temperature at time t + dt                [K]

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(  mesh%nV,   C%nz, u_times_dTdxp_upwind, wu_times_dTdxp_upwind)
    CALL allocate_shared_dp_2D(  mesh%nV,   C%nz, v_times_dTdyp_upwind, wv_times_dTdyp_upwind)
    CALL allocate_shared_dp_1D(  mesh%nV,         T_surf_annual       , wT_surf_annual       )
    CALL allocate_shared_dp_1D(  mesh%nV,         Q_base_grnd         , wQ_base_grnd         )
    CALL allocate_shared_dp_1D(  mesh%nV,         T_base_float        , wT_base_float        )
    CALL allocate_shared_dp_2D(  mesh%nV  , C%nz, Ti_tplusdt          , wTi_tplusdt          )
    CALL allocate_shared_int_1D( mesh%nV  ,       is_unstable         , wis_unstable         )

    ! Calculate zeta gradients
    CALL calc_zeta_gradients( mesh, ice)

    ! Calculate vertical ice velocities
    CALL calc_vertical_velocities( mesh, ice, BMB)

    ! Calculate upwind velocity times temperature gradients
    CALL calc_upwind_heat_flux_derivatives( mesh, ice, u_times_dTdxp_upwind, v_times_dTdyp_upwind)

    ! Calculate heating terms
    CALL calc_strain_heating(     mesh, ice)
    CALL calc_frictional_heating( mesh, ice)

    ! Calculate annual mean surface temperature
    DO vi = mesh%vi1, mesh%vi2
      T_surf_annual( vi) = SUM( climate%T2m( vi,:)) / REAL( SIZE( climate%T2m,2),dp)
    END DO

    ! For floating ice, basal temperatures are assumed to be always at
    ! the pressure melting point (since ocean water cannot be colder than
    ! this, and the ice itself cannot be warmer than this).
    DO vi = mesh%vi1, mesh%vi2
     T_base_float( vi) = ice%Ti_pmp_a( vi,C%nz)
    END DO

    ! Calculate heat flux at the base of the grounded ice
    DO vi = mesh%vi1, mesh%vi2
      Q_base_grnd( vi) = ice%frictional_heating_a( vi) + ice%GHF_a( vi)
    END DO

    ! Solve the heat equation for all vertices
    is_unstable( mesh%vi1:mesh%vi2) = 0
    n_unstable                      = 0
    DO vi = mesh%vi1, mesh%vi2

      ! Set temperature for ice-free vertices to annual mean surface temperature
      IF (ice%mask_ice_a( vi) == 0) THEN
        is_unstable( vi) = 0
        Ti_tplusdt( vi,:) = T_surf_annual( vi)
        CYCLE
      END IF

      ! For very thin ice, just let the profile equal the surface temperature
      IF (ice%Hi_a( vi) < 10._dp) THEN
        is_unstable( vi) = 0
        Ti_tplusdt( vi,:) = T_surf_annual( vi)
        CYCLE
      END IF

      ! Gather all the data needed to solve the 1-D heat equation
      icecol_Ti                   = ice%Ti_a(               vi,:)
      icecol_u                    = ice%u_3D_a(             vi,:)
      icecol_v                    = ice%v_3D_a(             vi,:)
      icecol_w                    = ice%w_3D_a(             vi,:)
      icecol_u_times_dTdxp_upwind = u_times_dTdxp_upwind(   vi,:)
      icecol_v_times_dTdyp_upwind = v_times_dTdyp_upwind(   vi,:)
      icecol_Ti_pmp               = ice%Ti_pmp_a(           vi,:)
      icecol_Ki                   = ice%Ki_a(               vi,:)
      icecol_Cpi                  = ice%Cpi_a(              vi,:)
      icecol_dzeta_dx             = ice%dzeta_dx_ak(        vi,:)
      icecol_dzeta_dy             = ice%dzeta_dy_ak(        vi,:)
      icecol_dzeta_dz             = ice%dzeta_dz_ak(        vi,:)
      icecol_dzeta_dt             = ice%dzeta_dt_ak(        vi,:)
      icecol_Phi                  = ice%internal_heating_a( vi,:)

      ! Solve the 1-D heat equation in the vertical column;
      ! if the solution is unstable, try again with a smaller timestep

      found_stable_solution = .FALSE.
      it_dt = 0
      DO WHILE ((.NOT. found_stable_solution) .AND. it_dt < 10)

        it_dt = it_dt + 1
        dt_applied = dt * (0.5_dp**(REAL( it_dt-1,dp)))   ! When it_dt = 0, dt_applied = dt; when it_dt = 1, dt_applied = dt/2, etc.

        ! If dt_applied = dt, solve it once; if dt_applied = dt/2, solve it twice; etc.
        DO it_it_dt = 1, 2**(it_dt-1)

          ! Solve the heat equation in the vertical column
          IF     (ice%mask_sheet_a( vi) == 1) THEN
            ! Grounded ice: use Q_base_grnd as boundary condition

            CALL solve_1D_heat_equation( vi, mesh, ice%Hi_a( vi), icecol_Ti, icecol_u, icecol_v, icecol_w, &
              icecol_u_times_dTdxp_upwind, icecol_v_times_dTdyp_upwind, T_surf_annual( vi), &
              icecol_Ti_pmp, icecol_Ki, icecol_Cpi, icecol_dzeta_dx, icecol_dzeta_dy, icecol_dzeta_dz, icecol_dzeta_dt, &
              icecol_Phi, dt_applied, icecol_Ti_tplusdt, Q_base_grnd = Q_base_grnd( vi))

          ELSEIF (ice%mask_shelf_a( vi) == 1) THEN
            ! Floating ice: use T_base_float as boundary condition

            CALL solve_1D_heat_equation( vi, mesh, ice%Hi_a( vi), icecol_Ti, icecol_u, icecol_v, icecol_w, &
              icecol_u_times_dTdxp_upwind, icecol_v_times_dTdyp_upwind, T_surf_annual( vi), &
              icecol_Ti_pmp, icecol_Ki, icecol_Cpi, icecol_dzeta_dx, icecol_dzeta_dy, icecol_dzeta_dz, icecol_dzeta_dt, &
              icecol_Phi, dt_applied, icecol_Ti_tplusdt, T_base_float = T_base_float( vi))

          ELSE
            CALL crash('mask_ice = 1 but mask_sheet and mask_shelf are 0!')
          END IF

          ! Update temperature solution for next semi-time-step
          icecol_Ti = icecol_Ti_tplusdt

        END DO ! DO it_it_dt = 1, 2**(it_dt-1)

        ! Check if we found a stable solution
        found_stable_solution = .TRUE.
        DO k = 1, C%nz
          IF (icecol_Ti( k) /= icecol_Ti( k)) found_stable_solution = .FALSE.   ! If we found NaN, the solution is unstable
          IF (icecol_Ti( k) < 180._dp)        found_stable_solution = .FALSE.   ! If we found temperatures below 180 K, the solution is unstable
          IF (icecol_Ti( k) > T0)             found_stable_solution = .FALSE.   ! If we found temperatures above freezing point, the solution is unstable
        END DO

      END DO ! DO WHILE ((.NOT. found_stable_solution) .AND. it_dt < 10)

      ! If the solution is still unstable, set the vertical temperature profile
      ! here to the Robin solution. Not the most accurate, but at least it will
      ! keep the model running. If too many grid cells need this, crash the model.
      IF (.NOT. found_stable_solution) THEN
        is_unstable( vi) = 1
        n_unstable = n_unstable + 1
      END IF

      ! Copy temperature solution
      Ti_tplusdt( vi,:) = icecol_Ti

    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Cope with instability
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_unstable, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    IF (n_unstable < CEILING( REAL( mesh%nV) / 100._dp)) THEN
      ! Instability is limited to an acceptably small number (< 1%) of grid cells;
      ! replace the temperature profile in those cells with the Robin solution

      DO vi = mesh%vi1, mesh%vi2
        IF (is_unstable( vi) == 1) CALL replace_Ti_with_robin_solution( ice, climate, SMB, Ti_tplusdt, vi)
      END DO
      CALL sync

    ELSE
      ! An unacceptably large number of grid cells was unstable; throw an error.

      CALL save_variable_as_netcdf_int_1D( ice%mask_ice_a          , 'mask_ice_a'          )
      CALL save_variable_as_netcdf_dp_1D(  ice%Hi_a                , 'Hi_a'                )
      CALL save_variable_as_netcdf_dp_1D(  ice%Hs_a                , 'Hs_a'                )
      CALL save_variable_as_netcdf_dp_2D(  ice%dzeta_dx_ak         , 'dzeta_dx_ak'         )
      CALL save_variable_as_netcdf_dp_2D(  ice%dzeta_dy_ak         , 'dzeta_dy_ak'         )
      CALL save_variable_as_netcdf_dp_2D(  ice%dzeta_dz_ak         , 'dzeta_dz_ak'         )
      CALL save_variable_as_netcdf_dp_2D(  ice%dzeta_dt_ak         , 'dzeta_dt_ak'         )
      CALL save_variable_as_netcdf_dp_2D(  ice%u_3D_a              , 'u_3D_a'              )
      CALL save_variable_as_netcdf_dp_2D(  ice%v_3D_a              , 'v_3D_a'              )
      CALL save_variable_as_netcdf_dp_2D(  ice%w_3D_a              , 'w_3D_a'              )
      CALL save_variable_as_netcdf_dp_1D(  ice%frictional_heating_a, 'frictional_heating_a')
      CALL save_variable_as_netcdf_dp_2D(  ice%internal_heating_a  , 'internal_heating_a'  )
      CALL save_variable_as_netcdf_dp_1D(  T_surf_annual           , 'T_surf_annual'       )
      CALL save_variable_as_netcdf_dp_1D(  Q_base_grnd             , 'Q_base_grnd'         )
      CALL save_variable_as_netcdf_dp_1D(  T_base_float            , 'T_base_float'        )
      CALL save_variable_as_netcdf_dp_2D(  ice%Ti_a                , 'Ti_a'                )
      CALL save_variable_as_netcdf_dp_2D(  Ti_tplusdt              , 'Ti_tplusdt'          )
      CALL save_variable_as_netcdf_int_1D( is_unstable             , 'is_unstable'         )
      CALL write_to_debug_file
      CALL crash('heat equation solver unstable for more than 1% of vertices!')

    END IF

    ! Move the new temperature field to the ice data structure
    ice%Ti_a( mesh%vi1:mesh%vi2,:) = Ti_tplusdt( mesh%vi1:mesh%vi2,:)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wu_times_dTdxp_upwind)
    CALL deallocate_shared( wv_times_dTdyp_upwind)
    CALL deallocate_shared( wT_surf_annual       )
    CALL deallocate_shared( wQ_base_grnd         )
    CALL deallocate_shared( wT_base_float        )
    CALL deallocate_shared( wTi_tplusdt          )
    CALL deallocate_shared( wis_unstable         )

    ! Safety
    CALL check_for_NaN_dp_2D( ice%Ti_a, 'ice%Ti_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_3D_heat_equation

  SUBROUTINE solve_1D_heat_equation( vi, mesh, Hi, Ti, u, v, w, u_times_dTdxp_upwind, v_times_dTdyp_upwind, T_surf, &
    Ti_pmp, Ki, Cpi, dzeta_dx, dzeta_dy, dzeta_dz, dzeta_dt, Phi, dt, Ti_tplusdt, Q_base_grnd, T_base_float)
    ! Solve the heat equation in the vertical column, i.e. using implicit discretisation of dT/dz, but explicit for dT/dx, dT/dy
    !
    ! The general heat equation (i.e. conservation of energy) inside the ice reads:
    !
    !   dT/dt = k / (rho cp) grad^2 T - u dT/dx - v dT/dy - w dT/dz + Phi / (rho cp)
    !
    ! With the following quantities:
    !
    !   T     - ice temperature
    !   k     - thermal conductivity of the ice
    !   rho   - density of the ice
    !   cp    - specific heat of the ice
    !   u,v,w - velocity of the ice
    !   Phi   - internal heating of the ice (due to strain rates)
    !
    ! By neglecting horizontal diffusion of heat (which can be done because of the flat
    ! geometry of the ice sheet), the grad^2 operator reduces to d2/dz2:
    !
    !   dT/dt = k / (rho cp) d2T/dz2 - u dT/dx - v dT/dy - w dT/dz + Phi / (rho cp)
    !
    ! Transforming into [xp,yp,zeta]-coordinates (xp = x, yp = y, zeta = (Hs(x,y) - z) / Hi(x,y)) yields:
    !
    !   dT/dt + dT/dzeta dzeta/dt = k / (rho cp) d2T/dzeta2 (dzeta/dz)^2 ...
    !                                - u (dT/dxp + dT/dzeta dzeta/dx) ...
    !                                - v (dT/dyp + dT/dzeta dzeta/dy) ...
    !                                - w (         dT/dzeta dzeta/dz) ...
    !                                + Phi / (rho cp)
    !
    ! The horizontal temperature gradients dT/dxp, dT/dyp are discretised explicitly in time,
    ! whereas the vertical gradients dT/dzeta, d2T/dzeta2 are discretised implicitly, yielding:
    !
    !   (T( t+dt) - T( t)) / dt + dzeta/dt d/dzeta( T( t+dt)) = ...
    !       k / (rho cp) (dzeta/dz)^2 d2/dzeta2( T( t+dt)) ...
    !     - u (dT( t)/dxp + dzeta/dx  d /dzeta ( T( t+dt))) ...
    !     - v (dT( t)/dyp + dzeta/dy  d /dzeta ( T( t+dt))) ...
    !     - w (             dzeta/dz  d /dzeta ( T( t+dt))) ...
    !     + Phi / (rho cp)
    !
    ! Moving all terms involving the unknown T( t+dt) to the left-hand side,
    ! and all other terms to the right-hand side, yields:
    !
    !   [ 1 / dt + dzeta/dt d/dzeta ...
    !     - k / (rho cp) (dzeta/dz)^2 d2/dzeta2 ...
    !     + u dzeta/dx d/dzeta ...
    !     + v dzeta/dy d/dzeta ...
    !     + w dzeta/dz d/dzeta ] T( t+dt) = ...
    !      T( t) / dt - u dT( t)/dxp - v dT( t)/dyp + Phi / (rho cp)
    !
    ! This can be further rearranged to read:
    !
    !   [ 1/dt + (dzeta/dt + u dzeta/dx + v dzeta/dy + w dzeta/dz) d/dzeta - k / (rho cp) (dzeta/dz)^2 d2/dzeta2 ] T( t+dt) = ...
    !      T( t) / dt - u dT( t)/dxp - v dT( t)/dyp + Phi / (rho cp)
    !
    ! Discretising this expression in space, the d/dzeta and d2/dzeta2 operators on the
    ! left-hand side become matrices, while all other terms become vectors. The equation
    ! thus describes a system of linear equations, with T( t+dt) the solution.
    !
    ! In order to solve this, boundary conditions must be applied at the surface and
    ! base of the ice. At the surface, the ice temperature is assumed to equal the
    ! annual mean surface temperature. At the base of grounded ice, the temperature
    ! gradient dT/dz must follow the heat flux at the base: dT/dz = -Q_base / k, and
    ! the temperature must also not exceed the pressure melting point temperature.
    ! At the base of floating ice, the temperature is assumed to equal the temperature
    ! of the ocean, limited by the pressure melting point temperature.
    !
    ! NOTE: with the current vertical discretisation scheme, d/dzeta and d2/dzeta2 are
    !       tridiagonal matrices, so that the entire matrix is tridiagonal. The LAPACK
    !       solver depends on this fact; changing the vertical discretisation to a
    !       higher-order scheme will require a different matrix solver.

    IMPLICIT NONE

    ! In/output variables
    INTEGER, INTENT(IN) :: vi
    TYPE(type_mesh),                     INTENT(IN)    :: mesh                 ! Contains the d/dzeta, d2/dzeta2 operators
    REAL(dp),                            INTENT(IN)    :: Hi                   ! Thickness of the ice column                                       [m]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: Ti                   ! Vertical profile of ice temperature at time t                     [K]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: u                    !   "         "    "  horizontal ice velocity in the x-direction    [m yr^-1]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: v                    !   "         "    "  horizontal ice velocity in the y-direction    [m yr^-1]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: w                    !   "         "    "  vertical   ice velocity                       [m yr^-1]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: u_times_dTdxp_upwind !   "         "    "  u * dT/dxp in the upwind direction            [K yr^-1]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: v_times_dTdyp_upwind !   "         "    "  u * dT/dxp in the upwind direction            [K yr^-1]
    REAL(dp),                            INTENT(IN)    :: T_surf               ! Annual mean surface temperature                                   [K]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: Ti_pmp               !   "         "    "  pressure melting point temperature            [K]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: Ki                   !   "         "    "  thermal conductivity of ice                   [J yr^-1 m^-1 K^-1]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: Cpi                  !   "         "    "  specific heat of ice                          [J kg^-1 K^-1]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: dzeta_dx             !   "         "    "  x-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: dzeta_dy             !   "         "    "  y-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: dzeta_dz             !   "         "    "  z-derivative of the scaled coordinate zeta    [m^-1]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: dzeta_dt             !   "         "    "  time-derivative of the scaled coordinate zeta [yr^-1]
    REAL(dp), DIMENSION( C%nz),          INTENT(IN)    :: Phi                  !   "         "    "  internal heat production                      [J kg^1 yr^1]
    REAL(dp),                            INTENT(IN)    :: dt                   ! Time step                                                         [yr]
    REAL(dp), DIMENSION( C%nz),          INTENT(OUT)   :: Ti_tplusdt           ! Vertical profile of ice temperature at time t + dt                [K]
    REAL(dp),                            INTENT(IN), OPTIONAL :: Q_base_grnd   ! Heat flux at the base of grounded ice                             [J m^-2 yr^-1]
    REAL(dp),                            INTENT(IN), OPTIONAL :: T_base_float  ! Heat flux at the ice base                                         [J m^-2 yr^-1]

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_1D_heat_equation'
    REAL(dp), DIMENSION( C%nz-1)                       :: AA_ldiag           ! Lower diagonal of A
    REAL(dp), DIMENSION( C%nz  )                       :: AA_diag            !       Diagonal of A
    REAL(dp), DIMENSION( C%nz-1)                       :: AA_udiag           ! Upper diagonal of A
    REAL(dp), DIMENSION( C%nz)                         :: bb                 ! Right-hand side b
    INTEGER                                            :: k
    REAL(dp)                                           :: dzeta
    REAL(dp)                                           :: c_ddzeta, c_d2dzeta2

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Initialise
    AA_ldiag = 0._dp
    AA_diag  = 0._dp
    AA_udiag = 0._dp
    bb       = 0._dp

    ! Ice surface boundary conditions: T = MIN( T_surf, T0)
    AA_diag(  1) = 1._dp
    bb(       1) = MIN( T_surf, T0)

    ! Ice base: either dT( t+dt)/dz = -Q_base / k, T( t+dt) <= T_pmp  for grounded ice, or:
    !                   T           = MIN( T_base_float, T_pmp)       for floating ice

    IF (PRESENT( Q_base_grnd)) THEN
      ! For grounded ice, let dT/dz = -Q_base / k

      ! Safety
      IF (PRESENT( T_base_float)) CALL crash('must provide either Q_base_grnd or T_base_float, but not both!')

      AA_diag(  C%nz) = 1._dp
      bb(       C%nz) = MIN( Ti_pmp( C%nz), Ti( C%nz-1) - (C%zeta( C%nz) - C%zeta( C%nz-1)) * Q_base_grnd / (dzeta_dz( C%nz) * Ki( C%nz)))

    ELSEIF (PRESENT( T_base_float)) THEN
      ! For floating ice, let T = MIN( T_base_float, T_pmp)

      ! Safety
      IF (PRESENT( Q_base_grnd)) CALL crash('must provide either Q_base_grnd or T_base_float, but not both!')

      AA_diag(  C%nz) = 1._dp
      bb(       C%nz) = MIN( T_base_float, Ti_pmp( C%nz))

    ELSE
      CALL crash('must provide either Q_base_grnd or T_base_float, but not both!')
    END IF ! IF (PRESENT( Q_base_grnd)) THEN

    ! Ice column
    DO k = 2, C%nz-1

      ! Calculate matrix coefficients
      c_ddzeta   = dzeta_dt( k) + (u( k) * dzeta_dx( k)) + (v( k) * dzeta_dy( k)) + (w( k) * dzeta_dz( k))
      c_d2dzeta2 = -Ki( k) / (ice_density * Cpi( k)) * dzeta_dz( k)**2

      AA_ldiag( k-1) =              (c_ddzeta * mesh%M_ddzeta_k_k_ldiag( k-1)) + (c_d2dzeta2 * mesh%M_d2dzeta2_k_k_ldiag( k-1))
      AA_diag(  k  ) = 1._dp / dt + (c_ddzeta * mesh%M_ddzeta_k_k_diag(  k  )) + (c_d2dzeta2 * mesh%M_d2dzeta2_k_k_diag(  k  ))
      AA_udiag( k  ) =              (c_ddzeta * mesh%M_ddzeta_k_k_udiag( k  )) + (c_d2dzeta2 * mesh%M_d2dzeta2_k_k_udiag( k  ))

      ! Calculate right-hand side
      bb(       k) = Ti( k) / dt - u_times_dTdxp_upwind( k) - v_times_dTdyp_upwind( k) + Phi( k) / (ice_density * Cpi( k))

    END DO ! DO k = 1, C%nz

    ! Solve the tridiagonal matrix equation representing the heat equation for this grid cell
    Ti_tplusdt = tridiagonal_solve( AA_ldiag, AA_diag, AA_udiag, bb)

    ! Make sure ice temperature doesn't exceed pressure melting point
    DO k = 1, C%nz
      Ti_tplusdt( k) = MIN( Ti_tplusdt( k), Ti_pmp( k))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_1D_heat_equation

! == The Robin temperature solution
  SUBROUTINE replace_Ti_with_robin_solution( ice, climate, SMB, Ti, vi)
    ! This function calculates for one horizontal grid point the temperature profiles
    ! using the surface temperature and the geothermal heat flux as boundary conditions.
    ! See Robin solution in: Cuffey & Paterson 2010, 4th ed, chapter 9, eq. (9.13) - (9.22).

    USE parameters_module, ONLY: pi, sec_per_year

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_snapshot_regional), INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(:,:  ),           INTENT(INOUT) :: Ti
    INTEGER,                              INTENT(IN)    :: vi

    ! Local variables:
    INTEGER                                            :: k
    REAL(dp)                                           :: Ts
    REAL(dp)                                           :: thermal_length_scale
    REAL(dp)                                           :: distance_above_bed
    REAL(dp)                                           :: erf1
    REAL(dp)                                           :: erf2

    REAL(dp)                                           :: thermal_conductivity_robin
    REAL(dp)                                           :: thermal_diffusivity_robin
    REAL(dp)                                           :: bottom_temperature_gradient_robin

    REAL(dp), PARAMETER                                :: kappa_0_ice_conductivity     = 9.828_dp                   ! The linear constant in the thermal conductivity of ice [J m^-1 K^-1 s^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    REAL(dp), PARAMETER                                :: kappa_e_ice_conductivity     = 0.0057_dp                  ! The exponent constant in the thermal conductivity of ice [K^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    REAL(dp), PARAMETER                                :: c_0_specific_heat            = 2127.5_dp                  ! The constant in the specific heat capacity of ice [J kg^-1 K^-1], see equation (12.5), Zwinger (2007), Cuffey & Paterson (2010, p. 400)

    thermal_conductivity_robin        = kappa_0_ice_conductivity * sec_per_year * EXP(-kappa_e_ice_conductivity * T0) ! Thermal conductivity            [J m^-1 K^-1 y^-1]
    thermal_diffusivity_robin         = thermal_conductivity_robin / (ice_density * c_0_specific_heat)                ! Thermal diffusivity             [m^2 y^-1]
    bottom_temperature_gradient_robin = - ice%GHF_a( vi) / thermal_conductivity_robin                                 ! Temperature gradient at bedrock

    Ts = MIN( T0, SUM( climate%T2m( vi,:)) / REAL( SIZE( climate%T2m,2),dp))

    IF (ice%mask_sheet_a( vi) == 1 ) THEN

      IF (SMB%SMB_year( vi) > 0._dp) THEN
        ! The Robin solution can be used to estimate the subsurface temperature profile in an accumulation area

        thermal_length_scale = SQRT(2._dp * thermal_diffusivity_robin * ice%Hi_a( vi) / SMB%SMB_year( vi))
        DO k = 1, C%nz
          distance_above_bed = (1._dp - C%zeta(k)) * ice%Hi_a( vi)
          erf1 = erf( distance_above_bed / thermal_length_scale)
          erf2 = erf( ice%Hi_a( vi) / thermal_length_scale)
          Ti( vi,k) = Ts + SQRT(pi) / 2._dp * thermal_length_scale * bottom_temperature_gradient_robin * (erf1 - erf2)
        END DO

      ELSE

        ! Ablation area: use linear temperature profile from Ts to (offset below) T_pmp
        Ti( vi,:) = Ts + ((T0 - CC * ice%Hi_a( vi)) - Ts) * C%zeta(:)

      END IF

    ELSEIF( ice%mask_shelf_a(vi) == 1) THEN
      ! Use a linear profile between T_surf and Ti_pmp_base

      Ti( vi,:) = Ts + C%zeta * (ice%Ti_pmp_a( vi,C%nz) - Ts)

    ELSE

      ! No ice present: use Ts everywhere
      Ti( vi,:) = Ts

    END IF

    ! Correct all temperatures above T_pmp:
    DO k = 1, C%nz
      Ti( vi,k) = MIN( Ti( vi,k), ice%Ti_pmp_a( vi,k))
    END DO

  END SUBROUTINE replace_Ti_with_robin_solution

! == Calculate various physical terms
  SUBROUTINE calc_strain_heating( mesh, ice)
    ! Calculate internal heating due to strain rates
    !
    ! Bueler and Brown (2009), Eq. 8 (though they use Sigma instead of Phi):
    !
    !   Phi = 2 B(T*) D^(1/n + 1) = 2 A(T*)^(-1/n) D^(1/n + 1)
    !
    ! From the text just after their Eq. 6:
    !
    !   2D^2 = Dij Dij, so: D = SQRT( Dij Dij / 2)
    !
    ! Here, Dij is the strain rate tensor: Dij = 1/2 (dui/dxj + duj/dxi)
    !
    !         |          du/dx          1/2 (du/dy + dv/dx)     1/2 (du/dz + dw/dx) |
    !         |                                                                     |
    !   Dij = | 1/2 (du/dy + dv/dx)              dv/dy          1/2 (dv/dz + dw/dy) |
    !         |                                                                     |
    !         | 1/2 (du/dz + dw/dx)     1/2 (dv/dz + dw/dy)              dw/dz      |
    !
    ! So:
    !
    !   D = SQRT( 1/2 [ ...
    !                   (du/dx)^2     + 1/4 (du/dy + dv/dx)^2 + 1/4 (du/dz + dw/dx)^2 + ...
    !           1/4 (du/dy + dv/dx)^2 +         (dv/dy)^2     + 1/4 (dv/dz + dw/dy)^2 + ...
    !           1/4 (du/dz + dw/dx)^2 + 1/4 (dv/dx + dw/dy)^2 +         (dw/dz)^2 ])
    !
    !     = SQRT( 1/2 [ (du/dx)^2 + (dv/dy)^2 + (dv/dz)^2 + ...
    !                    1/2 (du/dy + dv/dx)^2 + 1/2 (du/dz + dw/dx)^2 + 1/2 (dv/dz + dw/dy)^2])

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_strain_heating'
    INTEGER                                            :: vi, k
    REAL(dp)                                           :: D

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
    DO k = 1, C%nz

      ! No ice means no heating
      IF (ice%mask_ice_a( vi) == 0) THEN
        ice%internal_heating_a( vi,k) = 0._dp
        CYCLE
      END IF

      ! Calculate the total strain rate D
      D = SQRT( 0.5_dp * (ice%du_dx_3D_a( vi,k)**2 + ice%dv_dy_3D_a( vi,k)**2 + ice%dw_dz_3D_a( vi,k)**2 + &
                          0.5_dp * (ice%du_dy_3D_a( vi,k) + ice%dv_dx_3D_a( vi,k))**2 + &
                          0.5_dp * (ice%du_dz_3D_a( vi,k) + ice%dw_dx_3D_a( vi,k))**2 + &
                          0.5_dp * (ice%dv_dz_3D_a( vi,k) + ice%dw_dy_3D_a( vi,k))**2 ))

      ! Calculate the strain heating rate Phi
      ice%internal_heating_a( vi,k) = 2._dp * ice%A_flow_3D_a( vi,k)**(-1._dp / C%n_flow) * D**(1._dp / C%n_flow + 1._dp)

    END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_strain_heating

  SUBROUTINE calc_frictional_heating( mesh, ice)
    ! Calculate frictional heating at the base due to sliding

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_frictional_heating'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! No sliding means no friction
    IF (C%choice_sliding_law == 'no_sliding') THEN
      ice%frictional_heating_a( mesh%vi1:mesh%vi2) = 0._dp
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Calculate frictional heating
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_sheet_a( vi) == 1) THEN
        ice%frictional_heating_a( vi) = ice%beta_b_a( vi) * ice%uabs_base_a( vi)
      ELSE
        ice%frictional_heating_a( vi) = 0._dp
      END IF
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_frictional_heating

  SUBROUTINE calc_heat_capacity( mesh, ice)
    ! Calculate the heat capacity of the ice

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_heat_capacity'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_heat_capacity == 'uniform') THEN
      ! Apply a uniform value for the heat capacity

      ice%Cpi_a( mesh%vi1:mesh%vi2,:) = C%uniform_ice_heat_capacity
      CALL sync

    ELSEIF (C%choice_ice_heat_capacity == 'Pounder1965') THEN
      ! Calculate the heat capacity of ice according to Pounder: The Physics of Ice (1965)

      DO vi = mesh%vi1, mesh%vi2
        ice%Cpi_a( vi,:) = 2115.3_dp + 7.79293_dp * (ice%Ti_a( vi,:) - T0)
      END DO
      CALL sync

    ELSE
      CALL crash('unknown choice_ice_heat_capacity "' // TRIM( C%choice_ice_heat_capacity) // '"!')
    END IF

    ! Safety
    CALL check_for_NaN_dp_2D( ice%Cpi_a, 'ice%Cpi_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_heat_capacity

  SUBROUTINE calc_thermal_conductivity( mesh, ice)
    ! Calculate the thermal conductivity of the ice

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_thermal_conductivity'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_thermal_conductivity == 'uniform') THEN
      ! Apply a uniform value for the thermal conductivity

      ice%Ki_a( mesh%vi1:mesh%vi2,:) = C%uniform_ice_thermal_conductivity
      CALL sync

    ELSEIF (C%choice_ice_thermal_conductivity == 'Ritz1987') THEN
      ! Calculate the thermal conductivity of ice according to Ritz (1987)

      DO vi = mesh%vi1, mesh%vi2
        ice%Ki_a( vi,:) = 3.101E+08_dp * EXP(-0.0057_dp * ice%Ti_a( vi,:))
      END DO
      CALL sync

    ELSE
      CALL crash('unknown choice_ice_thermal_conductivity "' // TRIM( C%choice_ice_thermal_conductivity) // '"!')
    END IF

    ! Safety
    CALL check_for_NaN_dp_2D( ice%Ki_a, 'ice%Ki_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_thermal_conductivity

  SUBROUTINE calc_pressure_melting_point( mesh, ice)
    ! Calculate the pressure melting point of the ice according to Huybrechts (1992)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pressure_melting_point'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      ice%Ti_pmp_a( vi,:) = T0 - CC * ice%Hi_a( vi) * C%zeta
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%Ti_pmp_a, 'ice%Ti_pmp_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pressure_melting_point

  SUBROUTINE calc_upwind_heat_flux_derivatives( mesh, ice, u_times_dTdxp_upwind, v_times_dTdyp_upwind)
    ! Calculate upwind heat flux derivatives at vertex vi, vertical layer k

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: u_times_dTdxp_upwind, v_times_dTdyp_upwind

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_upwind_heat_flux_derivatives'
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dTi_dxp_3D_b,  dTi_dyp_3D_b
    INTEGER                                            :: wdTi_dxp_3D_b, wdTi_dyp_3D_b
    INTEGER                                            :: vi, k, vti, ti, n1, n2, n3, vib, vic, ti_upwind
    REAL(dp), DIMENSION(2)                             :: u_upwind, ab, ac

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(  mesh%nTri, C%nz, dTi_dxp_3D_b, wdTi_dxp_3D_b)
    CALL allocate_shared_dp_2D(  mesh%nTri, C%nz, dTi_dyp_3D_b, wdTi_dyp_3D_b)

    ! Calculate dT/dxp, dT/dyp on the b-grid
    CALL ddx_a_to_b_3D( mesh, ice%Ti_a, dTi_dxp_3D_b)
    CALL ddy_a_to_b_3D( mesh, ice%Ti_a, dTi_dyp_3D_b)

    DO vi = mesh%vi1, mesh%vi2

      ! Exception for the trivial case of no ice
      IF (ice%mask_ice_a( vi) == 0) THEN
        u_times_dTdxp_upwind( vi,:) = 0._dp
        v_times_dTdyp_upwind( vi,:) = 0._dp
        CYCLE
      END IF

      ! The upwind velocity vector
      u_upwind = [-ice%u_vav_a( vi), -ice%v_vav_a( vi)]

      ! Find the upwind triangle
      ti_upwind = 0
      DO vti = 1, mesh%niTri( vi)

        ! Triangle ti is spanned counter-clockwise by vertices [vi,vib,vic]
        ti  = mesh%iTri( vi,vti)
        vib = 0
        vic = 0
        DO n1 = 1, 3
          n2 = n1 + 1
          IF (n2 == 4) n2 = 1
          n3 = n2 + 1
          IF (n3 == 4) n3 = 1

          IF (mesh%Tri( ti,n1) == vi) THEN
            vib = mesh%Tri( ti,n2)
            vic = mesh%Tri( ti,n3)
            EXIT
          END IF
        END DO

        ! Check if the upwind velocity vector points into this triangle
        ab = mesh%V( vib,:) - mesh%V( vi,:)
        ac = mesh%V( vic,:) - mesh%V( vi,:)

        IF (CROSS2( ab, u_upwind) >= 0._dp .AND. CROSS2( u_upwind, ac) >= 0._dp) THEN
          ti_upwind = ti
          EXIT
        END IF

      END DO ! DO iti = 1, mesh%niTri( vi)

      ! Safety
      IF (ti_upwind == 0) THEN
        CALL crash('could not find upwind triangle!')
      END IF

      ! Calculate u * dT/dx, v * dT/dy
      DO k = 1, C%nz
        u_times_dTdxp_upwind( vi,k) = ice%u_3D_b( ti_upwind,k) * dTi_dxp_3D_b( ti_upwind,k)
        v_times_dTdyp_upwind( vi,k) = ice%v_3D_b( ti_upwind,k) * dTi_dyp_3D_b( ti_upwind,k)
      END DO

      ! DENK DROM
      debug%int_2D_a_01( vi) = ti_upwind

    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wdTi_dxp_3D_b)
    CALL deallocate_shared( wdTi_dyp_3D_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_upwind_heat_flux_derivatives

  SUBROUTINE calc_ice_rheology( mesh, ice, time)
    ! Calculate the flow factor A in Glen's flow law

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_ice_rheology'
    INTEGER                                            :: vi,k
    REAL(dp), DIMENSION(C%nZ)                          :: prof
    REAL(dp), PARAMETER                                :: A_low_temp  = 1.14E-05_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: A_high_temp = 5.47E+10_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_low_temp  = 6.0E+04_dp    ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_high_temp = 13.9E+04_dp   ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp)                                           :: A_flow_MISMIP

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_rheology == 'uniform') THEN
      ! Apply a uniform value for the ice flow factor

      ice%A_flow_3D_a( mesh%vi1:mesh%vi2,:) = C%uniform_flow_factor
      CALL sync

    ELSEIF (C%choice_ice_rheology == 'Huybrechts1992') THEN

      ! Calculate the ice flow factor as a function of the ice temperature according to the Arrhenius relationship (Huybrechts, 1992)
      DO vi = mesh%vi1, mesh%vi2

        DO k = 1, C%nz
          IF (ice%mask_ice_a( vi) == 1) THEN
            IF (ice%Ti_a( vi,k) < 263.15_dp) THEN
              ice%A_flow_3D_a( vi,k) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti_a( vi,k)))
            ELSE
              ice%A_flow_3D_a( vi,k) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti_a( vi,k)))
            END IF
          ELSE
            IF (C%choice_ice_margin == 'BC') THEN
              ice%A_flow_3D_a( vi,k) = 0._dp
            ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
              ! In the "infinite slab" case, calculate effective viscosity everywhere
              ! (even when there's technically no ice present)
              ice%A_flow_3D_a( vi,k) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * 263.15_dp))
            ELSE
              CALL crash('unknown choice_ice_margin "' // TRIM( C%choice_ice_margin) // '"!')
            END IF
          END IF
        END DO ! DO k = 1, C%nz

      END DO
      CALL sync

    ELSEIF (C%choice_ice_rheology == 'MISMIP_mod') THEN
      ! The time-dependent, step-wise changing uniform flow factor in the MISMIP_mod experiment

      A_flow_MISMIP = 1.0E-16_dp
      IF     (time < 15000._dp) THEN
        A_flow_MISMIP = 1.0E-16_dp
      ELSEIF (time < 30000._dp) THEN
        A_flow_MISMIP = 1.0E-17_dp
      ELSEIF (time < 45000._dp) THEN
        A_flow_MISMIP = 1.0E-16_dp
      END IF

      ice%A_flow_3D_a(  mesh%vi1:mesh%vi2,:) = A_flow_MISMIP
      CALL sync

    ELSE
      CALL crash('unknown choice_ice_rheology "' // TRIM( C%choice_ice_rheology) // '"!')
    END IF

    ! Apply the flow enhancement factors
    DO vi = mesh%vi1, mesh%vi2
      IF     (ice%mask_sheet_a( vi) == 1) THEN
        ice%A_flow_3D_a( vi,:) = ice%A_flow_3D_a( vi,:) * C%m_enh_sheet
      ELSEIF (ice%mask_shelf_a( vi) == 1) THEN
        ice%A_flow_3D_a( vi,:) = ice%A_flow_3D_a( vi,:) * C%m_enh_shelf
      END IF
    END DO
    CALL sync

    ! Calculate vertical average
    DO vi = mesh%vi1, mesh%vi2
      prof = ice%A_flow_3D_a( vi,:)
      CALL vertical_average( C%zeta, prof, ice%A_flow_vav_a( vi))
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%A_flow_3D_a , 'ice%A_flow_3D_a' )
    CALL check_for_NaN_dp_1D( ice%A_flow_vav_a, 'ice%A_flow_vav_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_ice_rheology

! == Initialise the englacial ice temperature at the start of a simulation
  SUBROUTINE initialise_ice_temperature( mesh, ice, climate, SMB, region_name, restart)
    ! Initialise the englacial ice temperature at the start of a simulation

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_snapshot_regional), INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name
    TYPE(type_restart_data),              INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_temperature'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '  Initialising ice temperature profile "', TRIM(C%choice_initial_ice_temperature), '"...'

    IF     (C%choice_initial_ice_temperature == 'uniform') THEN
      ! Simple uniform temperature
      CALL initialise_ice_temperature_uniform( mesh, ice)
    ELSEIF (C%choice_initial_ice_temperature == 'linear') THEN
      ! Simple linear temperature profile
      CALL initialise_ice_temperature_linear( mesh, ice, climate)
    ELSEIF (C%choice_initial_ice_temperature == 'Robin') THEN
      ! Initialise with the Robin solution
      CALL initialise_ice_temperature_Robin( mesh, ice, climate, SMB)
    ELSEIF (C%choice_initial_ice_temperature == 'restart') THEN
      ! Initialise with the temperature field from the provided restart file
      CALL initialise_ice_temperature_restart( mesh, ice, restart)
    ELSE
      CALL crash('unknown choice_initial_ice_temperature "' // TRIM( C%choice_initial_ice_temperature) // '"!')
    END IF

    ! Initialise mask_ice_a_prev
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%Hi_a( vi) > 0._dp) THEN
        ice%mask_ice_a_prev( vi) = 1
      ELSE
        ice%mask_ice_a_prev( vi) = 0
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature

  SUBROUTINE initialise_ice_temperature_uniform( mesh, ice)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Simple uniform temperature

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_temperature_uniform'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      IF (ice%Hi_a( vi) > 0._dp) THEN
        ice%Ti_a( vi,:) = C%uniform_ice_temperature
      ELSE
        ice%Ti_a( vi,:) = 0._dp
      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_uniform

  SUBROUTINE initialise_ice_temperature_linear( mesh, ice, climate)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Simple linear temperature profile

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_snapshot_regional), INTENT(IN)    :: climate

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_ice_temperature_linear'
    INTEGER                                             :: vi,k
    REAL(dp)                                            :: T_surf_annual, T_PMP_base

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      T_surf_annual = MIN( T0, SUM( climate%T2m( vi,:)) / REAL( SIZE( climate%T2m( vi,:),1),dp))
      T_PMP_base    = T0 - CC * ice%Hi_a( vi)

      IF (ice%Hi_a( vi) > 0._dp) THEN
        DO k = 1, C%nz
          ice%Ti_a( vi,k) = ((1._dp - C%zeta( k)) * T_surf_annual) + (C%zeta( k) * T_PMP_base)
        END DO
      ELSE
        ice%Ti_a( vi,:) = T_surf_annual
      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_linear

  SUBROUTINE initialise_ice_temperature_Robin( mesh, ice, climate, SMB)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Initialise with the Robin solution

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_snapshot_regional), INTENT(IN)    :: climate
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_ice_temperature_Robin'
    INTEGER                                             :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate Ti_pmp
    CALL calc_pressure_melting_point( mesh, ice)

    ! Initialise with the Robin solution
    DO vi = mesh%vi1, mesh%vi2
      CALL replace_Ti_with_robin_solution( ice, climate, SMB, ice%Ti_a, vi)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_Robin

  SUBROUTINE initialise_ice_temperature_restart( mesh, ice, restart)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Initialise with the temperature field from the provided restart file

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_restart_data),              INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_ice_temperature_restart'
    INTEGER                                             :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '   Initialising ice temperatures using data read from restart file...'

    DO vi = mesh%vi1, mesh%vi2
      ice%Ti_a( vi,:) = restart%Ti( vi,:)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_restart

! == Remap englacial temperature
  SUBROUTINE remap_ice_temperature( mesh_old, mesh_new, ice)
    ! Remap englacial temperature

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_old
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_new
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_ice_temperature'
    INTEGER,  DIMENSION(:    ), POINTER                ::  mask_ice_a_old,  mask_ice_a_new
    INTEGER                                            :: wmask_ice_a_old, wmask_ice_a_new
    INTEGER                                            :: vi, vvi, vj
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  Ti_ext
    INTEGER                                            :: wTi_ext
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: Vmap, Vstack1, Vstack2
    INTEGER                                            :: VstackN1, VstackN2
    INTEGER                                            :: sti, n, it
    REAL(dp), DIMENSION(C%nz)                          :: Ti_av

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_1D( mesh_old%nV,       mask_ice_a_old, wmask_ice_a_old)
    CALL allocate_shared_int_1D( mesh_new%nV,       mask_ice_a_new, wmask_ice_a_new)
    CALL allocate_shared_dp_2D(  mesh_old%nV, C%nz, Ti_ext        , wTi_ext       )

    ! Fill in the old and new ice masks
    DO vi = mesh_old%vi1, mesh_old%vi2
      mask_ice_a_old( vi) = ice%mask_ice_a( vi)
    END DO
    DO vi = mesh_new%vi1, mesh_new%vi2
      IF (ice%Hi_a( vi) > 0._dp) THEN
        mask_ice_a_new( vi) = 1
      ELSE
        mask_ice_a_new( vi) = 0
      END IF
    END DO

  ! Extrapolate old ice temperature outside the ice to fill the entire domain
  ! =========================================================================

    ! Initialise
    Ti_ext( mesh_old%vi1:mesh_old%vi2,:) = ice%Ti_a( mesh_old%vi1:mesh_old%vi2,:)
    CALL sync

    IF (par%master) THEN

      ! Allocate map and stacks for extrapolation
      ALLOCATE( Vmap(    mesh_old%nV))
      ALLOCATE( Vstack1( mesh_old%nV))
      ALLOCATE( Vstack2( mesh_old%nV))

      ! Initialise the stack with all ice-free-next-to-ice-covered vertices
      ! (and also initialise the map)
      Vmap     = 0
      Vstack2  = 0
      VstackN2 = 0

      DO vi = 1, mesh_old%nV
        IF (mask_ice_a_old( vi) == 1) THEN
          Vmap( vi) = 2
        ELSE
          DO vvi = 1, mesh_old%nC( vi)
            vj = mesh_old%C( vi,vvi)
            IF (mask_ice_a_old( vj) == 1) THEN
              ! Vertex vi is ice-free, but adjacent to ice-covered vertex vj
              VMap( vi) = 1
              VstackN2 = VstackN2 + 1
              Vstack2(   VstackN2) = vi
              EXIT
            END IF
          END DO
        END IF
      END DO

      ! Perform a flood-fill-style extrapolation
      it = 0
      DO WHILE (VstackN2 > 0)

        it = it + 1

        ! Cycle stacks
        Vstack1( 1:VstackN2) = Vstack2( 1:VstackN2)
        VstackN1 = VstackN2
        Vstack2( 1:VstackN2) = 0
        VstackN2 = 0

        ! Extrapolate temperature values into data-less-next-to-data-filled pixels
        DO sti = 1, VstackN1

          vi = Vstack1( sti)

          n     = 0
          Ti_av = 0._dp

          DO vvi = 1, mesh_old%nC( vi)

            vj = mesh_old%C( vi,vvi)

            IF (VMap( vj) == 2) THEN
              n     = n     + 1
              Ti_av = Ti_av + Ti_ext( vj,:)
            END IF

          END DO ! DO vvi = 1, mesh_old%nC( vi)

          ! Extrapolate temperature by averaging over data-filled neighbours
          Ti_av = Ti_av / REAL( n,dp)
          Ti_ext( vi,:) = Ti_av

        END DO ! DO sti = 1: VstackN1

        ! Create new stack of data-less-next-to-data-filled pixels
        DO sti = 1, VstackN1

          vi = Vstack1( sti)

          ! Mark this pixel as data-filled on the Map
          Vmap( vi) = 2

          ! Add its data-less neighbours to the Stack
          DO vvi = 1, mesh_old%nC( vi)

            vj = mesh_old%C( vi,vvi)

            IF (Vmap( vj) == 0) THEN
              Vmap( vj) = 1
              VstackN2 = VstackN2 + 1
              Vstack2(   VstackN2) = vj
            END IF

          END DO ! DO vvi = 1, mesh_old%nC( vi)
        END DO ! DO sti = 1: VstackN1

      END DO ! DO WHILE (VstackN2 > 0)

      ! Clean up after yourself
      DEALLOCATE( Vmap   )
      DEALLOCATE( Vstack1)
      DEALLOCATE( Vstack2)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Remap the extrapolated temperature field
    CALL remap_field_dp_3D( mesh_old, mesh_new, Ti_ext, wTi_ext, 'cons_2nd_order')

    ! Reallocate ice temperature field, copy remapped data only for ice-covered pixels
    CALL reallocate_shared_dp_2D( mesh_new%nV, C%nz, ice%Ti_a, ice%wTi_a)

    DO vi = mesh_new%vi1, mesh_new%vi2
      IF (mask_ice_a_new( vi) == 1) ice%Ti_a( vi,:) = Ti_ext( vi,:)
    END DO

    ! Reallocate mask_ice_a_prev, fill it in (needed for the generic temperature update)
    CALL reallocate_shared_int_1D( mesh_new%nV, ice%mask_ice_a_prev, ice%wmask_ice_a_prev)
    ice%mask_ice_a_prev( mesh_new%vi1:mesh_new%vi2) = mask_ice_a_new( mesh_new%vi1:mesh_new%vi2)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wmask_ice_a_old)
    CALL deallocate_shared( wmask_ice_a_new)
    CALL deallocate_shared( wTi_ext        )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ice_temperature

END MODULE thermodynamics_module
