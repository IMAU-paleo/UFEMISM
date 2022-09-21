MODULE ice_dynamics_module

  ! Contains all the routines needed to calculate ice-sheet geometry at
  ! the next time step, including routines to determine said time step.

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module,                only : dp, C, routine_path, init_routine, &
                                                  finalise_routine, crash, warning
  use parallel_module,                     only : par, sync, ierr, cerr, partition_list
  use data_types_module,                   only : type_mesh, type_ice_model, type_model_region, &
                                                  type_reference_geometry, type_remapping_mesh_mesh
  use utilities_module,                    only : surface_elevation, thickness_above_floatation
  use reallocate_mod,                      only : reallocate_bounds
  use mesh_mapping_module,                 only : remap_field_dp_2D
  use general_ice_model_data_module,       only : update_general_ice_model_data, determine_masks
  use ice_velocity_module,                 only : solve_SIA, solve_SSA, solve_DIVA, remap_velocities, &
                                                  map_velocities_b_to_c_2D, initialise_velocity_solver
  use ice_thickness_module,                only : calc_dHi_dt
  use basal_conditions_and_sliding_module, only : initialise_basal_conditions, remap_basal_conditions
  use thermodynamics_module,               only : calc_ice_rheology, remap_ice_temperature
  use mpi_module,                          only : allgather_array
  use forcing_module,                      only : forcing

  implicit none

contains

! ===== Main =====
! ================

  SUBROUTINE run_ice_model( region, t_end)
    ! Calculate ice velocities and the resulting change in ice geometry

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ice_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_timestepping == 'direct') THEN
      ! CALL run_ice_dynamics_direct( region, t_end)
      call crash('choice_timestepping "direct" not implemented yet!')
    ELSEIF (C%choice_timestepping == 'pc') THEN
      CALL run_ice_dynamics_pc( region, t_end)
    ELSE
      CALL crash('unknown choice_timestepping "' // TRIM( C%choice_timestepping) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_model

  subroutine initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD)
    ! Allocate shared memory for all the data fields of the ice dynamical module, and
    ! initialise some of them

    implicit none

    ! In- and output variables
    type(type_mesh),               intent(in)    :: mesh
    type(type_ice_model),          intent(inout) :: ice
    type(type_reference_geometry), intent(in)    :: refgeo_init
    type(type_reference_geometry), intent(in)    :: refgeo_PD

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'initialise_ice_model'
    integer                                      :: vi

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write(*,"(A)") '  Initialising ice dynamics model...'
    end if
    call sync

    ! === Memory allocation ===
    ! =========================

    ! Allocate shared memory
    call allocate_ice_model( mesh, ice)

    ! === Predictor-corrector method ===
    ! ==================================

    ! Initialise some numbers for the predictor/corrector ice thickness update method
    ice%pc_zeta        = 1._dp
    ice%pc_eta         = C%pc_epsilon
    ice%pc_eta_prev    = C%pc_epsilon

    ! === Sea level ===
    ! =================

    if (C%is_restart) then

      ! Get sea level from restart data
      call crash('Sea level initialisation: restart not implement yet!')
      ! ice%SL_a( mesh%vi1:mesh%vi2) = restart%SL( mesh%vi1:mesh%vi2)

    else

      ! Get sea level at initial time
      select case (C%choice_sealevel_model)

        case ('fixed')
          ! Fixed sea level
          ice%SL_a( mesh%vi1:mesh%vi2) = C%fixed_sealevel

        case ('prescribed')
          ! Sea-level prescribed from external record file
          ice%SL_a( mesh%vi1:mesh%vi2) = forcing%sealevel_obs

        case ('eustatic')
          ! Eustatic sea level
          call crash('Sea level initialisation: eustatic method not implement yet!')
          ! ice%SL_a( mesh%vi1:mesh%vi2) = C%initial_guess_sealevel

        case ('SELEN')
          ! Sea level from SELEN
          call crash('Sea level initialisation: SELEN method not implement yet!')
          ! ice%SL_a( mesh%vi1:mesh%vi2) = C%initial_guess_sealevel

        case default
          ! Unknown case
          call crash('unknown choice_sealevel_model "' // &
                      TRIM( C%choice_sealevel_model) // '"!')

      end select

    end if

    ! Initialise with data from initial file
    do vi = mesh%vi1, mesh%vi2
      ! Main quantities
      ice%Hi_a( vi)   = refgeo_init%Hi( vi)
      ice%Hb_a( vi)   = refgeo_init%Hb( vi)

      if (C%is_restart) then
        ice%Hs_a( vi) = refgeo_init%Hs( vi)
      else
        ice%Hs_a( vi) = surface_elevation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))
      end if

      ice%TAF_a( vi)  = thickness_above_floatation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))

      ! Differences w.r.t. present-day
      ice%dHi_a( vi)  = ice%Hi_a( vi) - refgeo_PD%Hi( vi)
      ice%dHb_a( vi)  = ice%Hb_a( vi) - refgeo_PD%Hb( vi)
      ice%dHs_a( vi)  = ice%Hs_a( vi) - refgeo_PD%Hs( vi)
    end do

    ! Determine masks
    call determine_masks( mesh, ice)

    ice%dHb_dt_a = 0._dp
    ice%dHi_dt_a = 0._dp
    ice%dHs_dt_a = 0._dp

    ! Initialise the "previous ice mask", so that the first call to thermodynamics works correctly
    ice%mask_ice_a_prev( mesh%vi1:mesh%vi2) = ice%mask_ice_a( mesh%vi1:mesh%vi2)
    call allgather_array(ice%mask_ice_a_prev)

    ! Allocate and initialise basal conditions
    call initialise_basal_conditions( mesh, ice)

    ! Geothermal heat flux
    select case (C%choice_geothermal_heat_flux)

      case ('constant')
        ! Uniform value over whole domain
        ice%GHF_a( mesh%vi1:mesh%vi2) = C%constant_geothermal_heat_flux

      case ('spatial')
        ! Spatially variable field
        call crash ('spatially variable GHF not yet implemented!')
        ! call map_geothermal_heat_flux_to_mesh( mesh, ice)

      case default
        ! Unknown case
        call crash('unknown choice_geothermal_heat_flux "' // &
                    trim( C%choice_geothermal_heat_flux) // '"!')

    end select

    ! Initialise data and matrices for the velocity solver(s)
    call initialise_velocity_solver( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = huge( 1))

  end subroutine initialise_ice_model

! ===== Direct method =====
! =========================

  SUBROUTINE run_ice_dynamics_direct( region, t_end)
    ! Ice dynamics and time-stepping with the "direct" method
    ! NOTE: this does not work for DIVA ice dynamics!

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ice_dynamics_direct'
    REAL(dp)                                           :: dt_crit_SIA, dt_crit_SSA

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Calculate ice velocities with the selected ice-dynamical approximation
    ! ! ======================================================================

    ! IF     (C%choice_ice_dynamics == 'none') THEN
    !   ! Fixed ice geometry

    ! ELSEIF (C%choice_ice_dynamics == 'SIA') THEN
    !   ! Shallow ice approximation

    !   IF (region%time == region%t_next_SIA) THEN

    !     ! Calculate new ice velocities
    !     CALL solve_SIA( region%mesh, region%ice)

    !     ! Calculate critical time step
    !     CALL calc_critical_timestep_SIA( region%mesh, region%ice, dt_crit_SIA)
    !     IF (par%master) region%dt_crit_SIA = dt_crit_SIA

    !     ! Update timer
    !     IF (par%master) region%t_last_SIA = region%time
    !     IF (par%master) region%t_next_SIA = region%time + region%dt_crit_SIA
    !     CALL sync

    !   END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN

    ! ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
    !   ! Shallow shelf approximation

    !   IF (region%time == region%t_next_SSA) THEN

    !     ! Calculate new ice velocities
    !     CALL solve_SSA( region%mesh, region%ice)

    !     ! Calculate critical time step
    !     CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_SSA)
    !     IF (par%master) region%dt_crit_SSA = dt_crit_SSA

    !     ! Update timer
    !     IF (par%master) region%t_last_SSA = region%time
    !     IF (par%master) region%t_next_SSA = region%time + region%dt_crit_SSA
    !     CALL sync

    !   END IF ! IF (ABS(region%time - region%t_next_SSA) < dt_tol) THEN

    ! ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
    !   ! Hybrid SIA/SSA (Bueler and Brown, 2009)

    !   IF (region%time == region%t_next_SIA) THEN

    !     ! Calculate new ice velocities
    !     CALL solve_SIA( region%mesh, region%ice)

    !     ! Calculate critical time step
    !     CALL calc_critical_timestep_SIA( region%mesh, region%ice, dt_crit_SIA)
    !     IF (par%master) region%dt_crit_SIA = dt_crit_SIA

    !     ! Update timer
    !     IF (par%master) region%t_last_SIA = region%time
    !     IF (par%master) region%t_next_SIA = region%time + region%dt_crit_SIA
    !     CALL sync

    !   END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN

    !   IF (region%time == region%t_next_SSA) THEN

    !     ! Calculate new ice velocities
    !     CALL solve_SSA( region%mesh, region%ice)

    !     ! Calculate critical time step
    !     CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_SSA)
    !     IF (par%master) region%dt_crit_SSA = dt_crit_SSA

    !     ! Update timer
    !     IF (par%master) region%t_last_SSA = region%time
    !     IF (par%master) region%t_next_SSA = region%time + region%dt_crit_SSA
    !     CALL sync

    !   END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN

    ! ELSE ! IF     (C%choice_ice_dynamics == 'SIA') THEN
    !   CALL crash('"direct" time stepping works only with SIA, SSA, or SIA/SSA ice dynamics, not with DIVA!')
    ! END IF ! IF     (C%choice_ice_dynamics == 'SIA') THEN

    ! ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    ! CALL determine_timesteps_and_actions( region, t_end)

    ! !IF (par%master) WRITE(0,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_SIA = ', dt_crit_SIA, ', dt_crit_SSA = ', dt_crit_SSA, ', dt = ', region%dt

    ! ! Calculate new ice geometry
    ! ! ==========================

    ! IF (C%choice_ice_dynamics == 'none') THEN
    !   ! Fixed ice geometry
    !   region%ice%dHi_dt_a( region%mesh%vi1 : region%mesh%vi2) = 0._dp
    !   CALL sync
    ! ELSE
    !   CALL calc_dHi_dt( region%mesh, region%ice, region%SMB, region%BMB, region%dt, region%mask_noice, region%refgeo_PD)
    ! END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_direct

! ===== Predictor-corrector method =====
! ======================================

  SUBROUTINE run_ice_dynamics_pc( region, t_end)
    ! Ice dynamics and time-stepping with the predictor/correct method
    ! (adopted from Yelmo, originally based on Cheng et al., 2017)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ice_dynamics_pc'
    INTEGER                                            :: vi1,vi2
    LOGICAL                                            :: do_update_ice_velocity
    REAL(dp)                                           :: dt_from_pc, dt_crit_adv

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Abbreviations for cleaner code
    vi1 = region%mesh%vi1
    vi2 = region%mesh%vi2

    ! Determine whether or not we need to update ice velocities
    do_update_ice_velocity = .FALSE.
    IF     (C%choice_ice_dynamics == 'none') THEN
      region%ice%dHi_dt_a( vi1:vi2) = 0._dp
    ELSEIF (C%choice_ice_dynamics == 'SIA') THEN
      IF (region%time == region%t_next_SIA ) do_update_ice_velocity = .TRUE.
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      IF (region%time == region%t_next_SSA ) do_update_ice_velocity = .TRUE.
    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      IF (region%time == region%t_next_SIA ) do_update_ice_velocity = .TRUE.
    ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      IF (region%time == region%t_next_DIVA) do_update_ice_velocity = .TRUE.
    ELSE
      CALL crash('unknown choice_ice_dynamics "' // TRIM( C%choice_ice_dynamics) // '"!')
    END IF

    IF (do_update_ice_velocity) THEN

      ! Calculate time step based on the truncation error in ice thickness (Robinson et al., 2020, Eq. 33)
      CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_adv)

      region%dt_crit_ice_prev = region%dt_crit_ice
      region%ice%pc_eta_prev  = region%ice%pc_eta
      dt_from_pc              = (C%pc_epsilon / region%ice%pc_eta)**(C%pc_k_I + C%pc_k_p) * (C%pc_epsilon / region%ice%pc_eta_prev)**(-C%pc_k_p) * region%dt

      region%dt_crit_ice      = MAX(C%dt_min, MINVAL((/ C%dt_max, 2._dp * region%dt_crit_ice_prev, dt_crit_adv, dt_from_pc/)))

      region%ice%pc_zeta      = region%dt_crit_ice / region%dt_crit_ice_prev
      region%ice%pc_beta1     = 1._dp + region%ice%pc_zeta / 2._dp
      region%ice%pc_beta2     =       - region%ice%pc_zeta / 2._dp

      ! Predictor step
      ! ==============

      ! Calculate new ice geometry
      region%ice%pc_f2 = region%ice%dHi_dt_a
      CALL calc_dHi_dt( region%mesh, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%mask_noice, region%refgeo_PD)
      region%ice%pc_f1 = region%ice%dHi_dt_a
      region%ice%Hi_pred = MAX(0._dp, region%ice%Hi_a + region%dt_crit_ice * region%ice%dHi_dt_a)

      ! Update step
      ! ===========

      ! Calculate velocities for predicted geometry
      region%ice%Hi_old = region%ice%Hi_a
      region%ice%Hi_a   = region%ice%Hi_pred
      CALL update_general_ice_model_data( region%mesh, region%ice)

      IF     (C%choice_ice_dynamics == 'SIA') THEN

        ! Calculate velocities
       call solve_SIA( region%mesh, region%ice)

       ! Update timer
       region%t_last_SIA = region%time
       region%t_next_SIA = region%time + region%dt_crit_ice

      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN

        ! Calculate velocities
       CALL solve_SSA( region%mesh, region%ice)

       ! Update timer
       region%t_last_SSA = region%time
       region%t_next_SSA = region%time + region%dt_crit_ice


      ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN

        ! Calculate velocities
       CALL solve_SIA(  region%mesh, region%ice)
       CALL solve_SSA(  region%mesh, region%ice)

       ! Update timer
       region%t_last_SIA = region%time
       region%t_last_SSA = region%time
       region%t_next_SIA = region%time + region%dt_crit_ice
       region%t_next_SSA = region%time + region%dt_crit_ice

      ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN

        ! Calculate velocities
        CALL solve_DIVA( region%mesh, region%ice)

        ! Update timer
        region%t_last_DIVA = region%time
        region%t_next_DIVA = region%time + region%dt_crit_ice

      ELSE
        CALL crash('unknown choice_ice_dynamics "' // TRIM( C%choice_ice_dynamics) // '"!')
      END IF

      ! Corrector step
      ! ==============

      ! Go back to old ice thickness. Run all the other modules (climate, SMB, BMB, thermodynamics, etc.)
      ! and only go to new (corrected) ice thickness at the end of this time loop.
      region%ice%Hi_a = region%ice%Hi_old
      CALL update_general_ice_model_data( region%mesh, region%ice)

      ! Calculate "corrected" ice thickness based on new velocities
      CALL calc_dHi_dt( region%mesh, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%mask_noice, region%refgeo_PD)
      region%ice%pc_f3 = region%ice%dHi_dt_a
      region%ice%pc_f4 = region%ice%pc_f1
      region%ice%Hi_corr = MAX(0._dp, region%ice%Hi_a + 0.5_dp * region%dt_crit_ice * (region%ice%pc_f3 + region%ice%pc_f4))

      ! Determine truncation error
      CALL calc_pc_truncation_error( region%mesh, region%ice, region%dt_crit_ice, region%dt_prev)

    END IF ! IF (do_update_ice_velocity) THEN

    ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    CALL determine_timesteps_and_actions( region, t_end)

    ! ! Calculate ice thickness at the end of this model loop
    region%ice%Hi_tplusdt_a = MAX( 0._dp, region%ice%Hi_a + region%dt * region%ice%dHi_dt_a)

    ! IF (par%master) WRITE(*,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_adv = ', dt_crit_adv, ', dt_from_pc = ', dt_from_pc, ', dt = ', region%dt
    ! CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_pc

  SUBROUTINE calc_pc_truncation_error( mesh, ice, dt, dt_prev)
    ! Calculate the truncation error in the ice thickness rate of change (Robinson et al., 2020, Eq. 32)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: dt, dt_prev

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pc_truncation_error'
    INTEGER                                            :: vi, ci, vc
    LOGICAL                                            :: has_GL_neighbour
    REAL(dp)                                           :: zeta, eta_proc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Ratio of time steps
    ! zeta = dt / dt_prev

    ! ! Find maximum truncation error
    ! eta_proc = C%pc_eta_min

    ! DO vi = mesh%vi1, mesh%vi2

    !   ! Calculate truncation error (Robinson et al., 2020, Eq. 32)
    !   ice%pc_tau( vi) = ABS( zeta * (ice%Hi_corr( vi) - ice%Hi_pred( vi)) / ((3._dp * zeta + 3._dp) * dt))

    !   IF (ice%mask_sheet_a( vi) == 1) THEN

    !     has_GL_neighbour = .FALSE.
    !     DO ci = 1, mesh%nC( vi)
    !       vc = mesh%C( vi,ci)
    !       IF (ice%mask_gl_a( vc) == 1) THEN
    !         has_GL_neighbour = .TRUE.
    !         EXIT
    !       END IF
    !     END DO

    !     IF (.NOT.has_GL_neighbour) eta_proc = MAX( eta_proc, ice%pc_tau( vi))

    !   END IF

    ! END DO
    ! CALL MPI_REDUCE( eta_proc, ice%pc_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pc_truncation_error

! ===== Ice thickness update =====
! ================================

  subroutine update_ice_thickness( mesh, ice, refgeo_PD)
    ! Update the ice thickness at the end of a model time loop

    implicit none

    ! In- and output variables:
    type(type_mesh),               intent(in)    :: mesh
    type(type_ice_model),          intent(inout) :: ice
    type(type_reference_geometry), intent(in)    :: refgeo_PD

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'update_ice_thickness'

    ! Add routine to path
    call init_routine( routine_name)

    ! Save the previous ice mask, for use in thermodynamics
    ice%mask_ice_a_prev( mesh%vi1:mesh%vi2) = ice%mask_ice_a( mesh%vi1:mesh%vi2)
    call allgather_array(ice%mask_ice_a_prev)

    ! Set ice thickness to new value
    ice%Hi_a( mesh%vi1:mesh%vi2) = max( 0._dp, ice%Hi_tplusdt_a( mesh%vi1:mesh%vi2))

    call update_general_ice_model_data(      mesh, ice)

    ! Compute ice thickness/elevation difference w.r.t PD
    ice%dHi_a( mesh%vi1:mesh%vi2) = ice%Hi_a( mesh%vi1:mesh%vi2) - refgeo_PD%Hi( mesh%vi1:mesh%vi2)
    ice%dHs_a( mesh%vi1:mesh%vi2) = ice%Hs_a( mesh%vi1:mesh%vi2) - refgeo_PD%Hs( mesh%vi1:mesh%vi2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_ice_thickness

! ===== Time stepping =====
! =========================

  SUBROUTINE calc_critical_timestep_SIA( mesh, ice, dt_crit_SIA)
    ! Calculate the critical time step for advective ice flow (CFL criterion)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp),                            INTENT(OUT)   :: dt_crit_SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_critical_timestep_SIA'
    INTEGER                                            :: aci, vi ,vj, k
    REAL(dp), DIMENSION(:    ), POINTER                ::  u_c,  v_c,  Hi_c,  dHs_dx_c,  dHs_dy_c
    INTEGER                                            :: wu_c, wv_c, wHi_c, wdHs_dx_c, wdHs_dy_c
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  u_3D_c,  v_3D_c
    INTEGER                                            :: wu_3D_c, wv_3D_c
    REAL(dp)                                           :: D_SIA, dist, dt
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Allocate shared memory
    ! CALL allocate_shared_dp_1D( mesh%nAc,       u_c     , wu_c     )
    ! CALL allocate_shared_dp_1D( mesh%nAc,       v_c     , wv_c     )
    ! CALL allocate_shared_dp_1D( mesh%nAc,       Hi_c    , wHi_c    )
    ! CALL allocate_shared_dp_1D( mesh%nAc,       dHs_dx_c, wdHs_dx_c)
    ! CALL allocate_shared_dp_1D( mesh%nAc,       dHs_dy_c, wdHs_dy_c)
    ! CALL allocate_shared_dp_2D( mesh%nAc, C%nz, u_3D_c  , wu_3D_c  )
    ! CALL allocate_shared_dp_2D( mesh%nAc, C%nz, v_3D_c  , wv_3D_c  )

    ! ! Calculate ice velocity and thickness, and surface slopes on the staggered grid
    ! CALL map_velocities_b_to_c_2D( mesh, ice%u_vav_b, ice%v_vav_b, u_c, v_c)
    ! CALL map_a_to_c_2D( mesh, ice%Hi_a, Hi_c    )
    ! CALL ddx_a_to_c_2D( mesh, ice%Hs_a, dHs_dx_c)
    ! CALL ddy_a_to_c_2D( mesh, ice%Hs_a, dHs_dy_c)
    ! CALL map_velocities_b_to_c_3D( mesh, ice%u_3D_b, ice%v_3D_b, u_3D_c, v_3D_c)

    ! ! Initialise time step with maximum allowed value
    ! dt_crit_SIA = C%dt_max

    ! DO aci = mesh%ci1, mesh%ci2

    !   ! Calculate the SIA ice diffusivity
    !   D_SIA = 1E-9_dp
    !   D_SIA = MAX( D_SIA, ABS( u_c( aci) * Hi_c( aci) / dHs_dx_c( aci) ))
    !   D_SIA = MAX( D_SIA, ABS( v_c( aci) * Hi_c( aci) / dHs_dy_c( aci) ))

    !   vi = mesh%Aci( aci,1)
    !   vj = mesh%Aci( aci,2)
    !   dist = NORM2( mesh%V( vi,:) - mesh%V( vj,:))
    !   dt = dist**2 / (6._dp * D_SIA)
    !   dt_crit_SIA = MIN(dt_crit_SIA, dt)

    !   ! Also check the 3D advective time step
    !   DO k = 1, C%nz
    !     dt = dist / (ABS( u_3D_c( aci,k)) + ABS( v_3D_c( aci,k)))
    !     dt_crit_SIA = MIN( dt_crit_SIA, dt)
    !   END DO

    ! END DO

    ! CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_SIA, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    ! dt_crit_SIA = dt_crit_SIA * dt_correction_factor

    ! ! Clean up after yourself
    ! CALL deallocate_shared( wu_c     )
    ! CALL deallocate_shared( wv_c     )
    ! CALL deallocate_shared( wHi_c    )
    ! CALL deallocate_shared( wdHs_dx_c)
    ! CALL deallocate_shared( wdHs_dy_c)
    ! CALL deallocate_shared( wu_3D_c  )
    ! CALL deallocate_shared( wv_3D_c  )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_critical_timestep_SIA

  SUBROUTINE calc_critical_timestep_adv( mesh, ice, dt_crit_adv)
    ! Calculate the critical time step for advective ice flow (CFL criterion)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp),                            INTENT(OUT)   :: dt_crit_adv

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_critical_timestep_adv'
    INTEGER                                            :: aci, vi, vj
    REAL(dp), DIMENSION(:    ), allocatable            :: u_c,  v_c
    REAL(dp), DIMENSION(:    ), allocatable            :: u_vav_b, v_vav_b
    REAL(dp)                                           :: dist, dt
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    allocate(u_c(mesh%ci1:mesh%ci2))
    allocate(v_c(mesh%ci1:mesh%ci2))

    allocate(u_vav_b(mesh%nTri))
    allocate(v_vav_b(mesh%nTri))
    u_vav_b(mesh%ti1:mesh%ti2) = ice%u_vav_b
    v_vav_b(mesh%ti1:mesh%ti2) = ice%v_vav_b
    call allgather_array(u_vav_b)
    call allgather_array(v_vav_b)

    ! Calculate ice velocity on the staggered grid
    CALL map_velocities_b_to_c_2D( mesh, u_vav_b, v_vav_b, u_c, v_c)

    dt_crit_adv = 2._dp * C%dt_max

    DO aci = mesh%ci1, mesh%ci2

      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      dist = NORM2( mesh%V( vi,:) - mesh%V( vj,:))
      dt = dist / (ABS( u_c( aci)) + ABS( v_c( aci)))
      dt_crit_adv = MIN( dt_crit_adv, dt)

    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_adv = MIN(C%dt_max, dt_crit_adv * dt_correction_factor)

    ! Clean up after yourself
    deallocate( u_c)
    deallocate( v_c)
    deallocate( u_vav_b)
    deallocate( v_vav_b)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_critical_timestep_adv

  SUBROUTINE determine_timesteps_and_actions( region, t_end)
    ! Determine how long we can run just ice dynamics before another "action" (thermodynamics,
    ! GIA, output writing, inverse routine, etc.) has to be performed, and adjust the time step accordingly.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_timesteps_and_actions'
    REAL(dp)                                           :: t_next

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine when each model components should be updated

    t_next = MIN(t_end, region%time + C%dt_max)

    ! First the ice dynamics
    ! ======================

    ! IF     (C%choice_ice_dynamics == 'none') THEN
    !   ! Just stick to the maximum time step
    ! ELSEIF (C%choice_ice_dynamics == 'SIA') THEN
    !   t_next = MIN( t_next, region%t_next_SIA)
    ! ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
    !   t_next = MIN( t_next, region%t_next_SSA)
    ! ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
    !   t_next = MIN( t_next, region%t_next_SIA)
    !   t_next = MIN( t_next, region%t_next_SSA)
    ! ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
    !   t_next = MIN( t_next, region%t_next_DIVA)
    ! ELSE
    !   CALL crash('unknown choice_ice_dynamics "' // TRIM( C%choice_ice_dynamics) // '"!')
    ! END IF ! IF (C%choice_ice_dynamics == 'SIA') THEN

    ! Then the other model components
    ! ===============================

    region%do_thermo  = .FALSE.
    IF (region%time == region%t_next_thermo) THEN
      region%do_thermo      = .TRUE.
      region%t_last_thermo  = region%time
      region%t_next_thermo  = region%t_last_thermo + C%dt_thermo
    END IF
    t_next = MIN( t_next, region%t_next_thermo)

    region%do_climate = .FALSE.
    IF (region%time == region%t_next_climate) THEN
      region%do_climate     = .TRUE.
      region%t_last_climate = region%time
      region%t_next_climate = region%t_last_climate + C%dt_climate
    END IF
    t_next = MIN( t_next, region%t_next_climate)

    region%do_ocean   = .FALSE.
    IF (region%time == region%t_next_ocean) THEN
      region%do_ocean       = .TRUE.
      region%t_last_ocean   = region%time
      region%t_next_ocean   = region%t_last_ocean + C%dt_ocean
    END IF
    t_next = MIN( t_next, region%t_next_ocean)

    region%do_SMB     = .FALSE.
    IF (region%time == region%t_next_SMB) THEN
      region%do_SMB         = .TRUE.
      region%t_last_SMB     = region%time
      region%t_next_SMB     = region%t_last_SMB + C%dt_SMB
    END IF
    t_next = MIN( t_next, region%t_next_SMB)

    region%do_BMB     = .FALSE.
    IF (region%time == region%t_next_BMB) THEN
      region%do_BMB         = .TRUE.
      region%t_last_BMB     = region%time
      region%t_next_BMB     = region%t_last_BMB + C%dt_BMB
    END IF
    t_next = MIN( t_next, region%t_next_BMB)

    region%do_ELRA    = .FALSE.
    IF (region%time == region%t_next_ELRA) THEN
      region%do_ELRA        = .TRUE.
      region%t_last_ELRA    = region%time
      region%t_next_ELRA    = region%t_last_ELRA + C%dt_bedrock_ELRA
    END IF
    t_next = MIN( t_next, region%t_next_ELRA)

    region%do_output  = .FALSE.
    IF (region%time == region%t_next_output) THEN
      region%do_output      = .TRUE.
      region%t_last_output  = region%time
      region%t_next_output  = region%t_last_output + C%dt_output
    END IF
    t_next = MIN( t_next, region%t_next_output)

    ! Set time step so that we move forward to the next action
    region%dt = t_next - region%time

    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_timesteps_and_actions

! ===== Allocation and remapping =====
! ====================================

  subroutine allocate_ice_model( mesh, ice)
    ! Allocate ice model variables

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'allocate_ice_model'

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory
    ! ===============

    ! Basic data - ice thickness, bedrock height, surface height, mask, 3D ice velocities and temperature
    allocate(  ice%Hi_a                (mesh%vi1:mesh%vi2               ))
    allocate(  ice%Hb_a                (mesh%vi1:mesh%vi2               ))
    allocate(  ice%Hs_a                (mesh%vi1:mesh%vi2               ))
    allocate(  ice%SL_a                (mesh%vi1:mesh%vi2               ))
    ice%SL_a = 0d0
    allocate(  ice%TAF_a               (mesh%vi1:mesh%vi2               ))
    allocate(  ice%Ti_a                (mesh%vi1:mesh%vi2 , C%nz        ))

    ! Ice velocities
    allocate(  ice%u_3D_a              (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%u_3D_a = 0d0
    allocate(  ice%v_3D_a              (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%v_3D_a = 0d0
    allocate(  ice%u_3D_b              (mesh%ti1:mesh%ti2 , C%nz        ))
    ice%u_3D_b = 0d0
    allocate(  ice%v_3D_b              (mesh%ti1:mesh%ti2 , C%nz        ))
    ice%v_3D_b = 0d0
    allocate(  ice%w_3D_a              (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%w_3D_a = 0d0

    allocate(  ice%u_vav_a             (mesh%vi1:mesh%vi2               ))
    ice%u_vav_a = 0d0
    allocate(  ice%v_vav_a             (mesh%vi1:mesh%vi2               ))
    ice%v_vav_a = 0d0
    allocate(  ice%u_vav_b             (mesh%ti1:mesh%ti2               ))
    ice%u_vav_b = 0d0
    allocate(  ice%v_vav_b             (mesh%ti1:mesh%ti2               ))
    ice%v_vav_b = 0d0
    allocate(  ice%uabs_vav_a          (mesh%vi1:mesh%vi2               ))
    ice%uabs_vav_a = 0d0
    allocate(  ice%uabs_vav_b          (mesh%ti1:mesh%ti2               ))
    ice%uabs_vav_b = 0d0

    allocate(  ice%u_surf_a            (mesh%vi1:mesh%vi2               ))
    ice%u_surf_a = 0d0
    allocate(  ice%v_surf_a            (mesh%vi1:mesh%vi2               ))
    ice%v_surf_a = 0d0
    allocate(  ice%u_surf_b            (mesh%ti1:mesh%ti2               ))
    ice%u_surf_b = 0d0
    allocate(  ice%v_surf_b            (mesh%ti1:mesh%ti2               ))
    ice%v_surf_b = 0d0
    allocate(  ice%uabs_surf_a         (mesh%vi1:mesh%vi2               ))
    ice%uabs_surf_a = 0d0
    allocate(  ice%uabs_surf_b         (mesh%ti1:mesh%ti2               ))
    ice%uabs_surf_b = 0d0

    allocate(  ice%u_base_a            (mesh%vi1:mesh%vi2               ))
    ice%u_base_a = 0d0
    allocate(  ice%v_base_a            (mesh%vi1:mesh%vi2               ))
    ice%v_base_a = 0d0
    allocate(  ice%u_base_b            (mesh%ti1:mesh%ti2               ))
    ice%u_base_b = 0d0
    allocate(  ice%v_base_b            (mesh%ti1:mesh%ti2               ))
    ice%v_base_b = 0d0
    allocate(  ice%uabs_base_a         (mesh%vi1:mesh%vi2               ))
    ice%uabs_base_a = 0d0
    allocate(  ice%uabs_base_b         (mesh%ti1:mesh%ti2               ))
    ice%uabs_base_b = 0d0

    ! Different masks
    allocate(  ice%mask_land_a         (       1:mesh%nV                ))
    allocate(  ice%mask_ocean_a        (       1:mesh%nV                ))
    allocate(  ice%mask_lake_a         (       1:mesh%nV                ))
    allocate(  ice%mask_ice_a          (       1:mesh%nV                ))
    allocate(  ice%mask_sheet_a        (       1:mesh%nV                ))
    allocate(  ice%mask_shelf_a        (       1:mesh%nV                ))
    allocate(  ice%mask_coast_a        (       1:mesh%nV                ))
    allocate(  ice%mask_margin_a       (       1:mesh%nV                ))
    allocate(  ice%mask_gl_a           (       1:mesh%nV                ))
    allocate(  ice%mask_cf_a           (       1:mesh%nV                ))
    allocate(  ice%mask_a              (       1:mesh%nV                ))
    allocate(  ice%f_grnd_a            (mesh%vi1:mesh%vi2               ))
    allocate(  ice%f_grnd_b            (mesh%ti1:mesh%ti2               ))

    ! Ice physical properties
    allocate(  ice%A_flow_3D_a         (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%A_flow_3D_a = 0d0
    allocate(  ice%A_flow_vav_a        (mesh%vi1:mesh%vi2               ))
    ice%A_flow_vav_a = 0d0
    allocate(  ice%Ti_pmp_a            (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%Ti_pmp_a = 0d0
    allocate(  ice%Cpi_a               (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%Cpi_a = 0d0
    allocate(  ice%Ki_a                (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%Ki_a = 0d0

    ! Zeta derivatives
    allocate(  ice%dzeta_dt_a          (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%dzeta_dt_a = 0d0
    allocate(  ice%dzeta_dx_a          (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%dzeta_dx_a = 0d0
    allocate(  ice%dzeta_dy_a          (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%dzeta_dy_a = 0d0
    allocate(  ice%dzeta_dz_a          (mesh%vi1:mesh%vi2               ))
    ice%dzeta_dz_a = 0d0

    ! Ice dynamics - ice thickness calculation
    allocate(  ice%dVi_in              (       1:mesh%nV  , mesh%nC_mem ))
    ice%dVi_in = 0d0
    allocate(  ice%dHi_dt_a            (mesh%vi1:mesh%vi2               ))
    ice%dHi_dt_a = 0d0
    allocate(  ice%Hi_tplusdt_a        (mesh%vi1:mesh%vi2               ))
    ice%Hi_tplusdt_a = 0d0
    allocate(  ice%dHs_dt_a            (mesh%vi1:mesh%vi2               ))
    ice%dHs_dt_a = 0d0

    ! Ice dynamics - predictor/corrector ice thickness update
    allocate(  ice%pc_tau              (mesh%vi1:mesh%vi2               ))
    ice%pc_tau = 0d0
    allocate(  ice%pc_fcb              (mesh%vi1:mesh%vi2               ))
    ice%pc_fcb = 0d0
    allocate(  ice%pc_f1               (mesh%vi1:mesh%vi2               ))
    ice%pc_f1 = 0d0
    allocate(  ice%pc_f2               (mesh%vi1:mesh%vi2               ))
    ice%pc_f2 = 0d0
    allocate(  ice%pc_f3               (mesh%vi1:mesh%vi2               ))
    ice%pc_f3 = 0d0
    allocate(  ice%pc_f4               (mesh%vi1:mesh%vi2               ))
    ice%pc_f4 = 0d0
    allocate(  ice%Hi_old              (mesh%vi1:mesh%vi2               ))
    ice%Hi_old = 0d0
    allocate(  ice%Hi_pred             (mesh%vi1:mesh%vi2               ))
    ice%Hi_pred = 0d0
    allocate(  ice%Hi_corr             (mesh%vi1:mesh%vi2               ))
    ice%Hi_corr = 0d0

    ! Thermodynamics
    allocate(  ice%mask_ice_a_prev     (       1:mesh%nV                ))
    ice%mask_ice_a_prev = 0d0
    allocate(  ice%internal_heating_a  (mesh%vi1:mesh%vi2 , C%nz        ))
    ice%internal_heating_a = 0d0
    allocate(  ice%frictional_heating_a(mesh%vi1:mesh%vi2               ))
    ice%frictional_heating_a = 0d0
    allocate(  ice%GHF_a               (mesh%vi1:mesh%vi2               ))
    ice%GHF_a = 0d0

    ! Mesh adaptation data
    allocate(  ice%surf_curv           (       1:mesh%nV                ))
    ice%surf_curv = 0d0
    allocate(  ice%log_velocity        (mesh%vi1:mesh%vi2               ))
    ice%log_velocity = 0d0

    ! GIA
    allocate(  ice%dHb_a               (mesh%vi1:mesh%vi2               ))
    ice%dHb_a = 0d0
    allocate(  ice%dHb_dt_a            (mesh%vi1:mesh%vi2               ))
    ice%dHb_dt_a = 0d0
    allocate(  ice%dSL_dt_a            (mesh%vi1:mesh%vi2               ))
    ice%dSL_dt_a = 0d0

    ! Useful extra stuff
    allocate(  ice%dHi_a               (mesh%vi1:mesh%vi2               ))
    ice%dHi_a = 0d0
    allocate(  ice%dHs_a               (mesh%vi1:mesh%vi2               ))
    ice%dHs_a = 0d0

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = huge( 1))

  end subroutine allocate_ice_model

  SUBROUTINE remap_ice_model( mesh_old, mesh_new, map, ice, refgeo_PD, time)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_ice_model'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! === Remapping ===
    ! =================

    ! The only fields that actually need to be mapped. The rest only needs memory reallocation.
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%Hi_a        , 'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%dHi_dt_a    , 'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%Hi_tplusdt_a, 'cons_2nd_order')

    ! Remove very thin ice resulting from remapping errors
    DO vi = mesh_new%vi1, mesh_new%vi2
      IF (ice%Hi_a( vi) < 0.1_dp) THEN
        ice%Hi_a( vi) = 0._dp
      END IF
    END DO

    ! Remove artificial ice at domain border (messes up with thermodynamics)
    DO vi = mesh_new%vi1, mesh_new%vi2
      IF ( .NOT. mesh_new%edge_index( vi) == 0) THEN
        ice%Hi_a( vi) = 0._dp
      END IF
    END DO

    ! Remap sea level
    call remap_field_dp_2D( mesh_old, mesh_new, map, ice%SL_a, 'cons_2nd_order')

    ! Remap englacial temperature (needs some special attention because of the discontinuity at the ice margin)
    CALL remap_ice_temperature( mesh_old, mesh_new, map, ice)

    ! Remap bedrock change and add up to PD to prevent accumulation of numerical diffusion
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%dHb_a, 'cons_2nd_order')
    call reallocate_bounds( ice%Hb_a, mesh_new%vi1, mesh_new%vi2 )
    ice%Hb_a = refgeo_PD%Hb + ice%dHb_a

    ! Geothermal heat flux
    CALL reallocate_bounds( ice%GHF_a, mesh_new%vi1, mesh_new%vi2 )
    IF     (C%choice_geothermal_heat_flux == 'constant') THEN
      ice%GHF_a = C%constant_geothermal_heat_flux
    ! ELSEIF (C%choice_geothermal_heat_flux == 'spatial') THEN
    !   CALL map_geothermal_heat_flux_to_mesh( mesh_new, ice)
    ELSE
      CALL crash('unknown choice_geothermal_heat_flux "' // TRIM( C%choice_geothermal_heat_flux) // '"!')
    END IF

    ! Basal conditions
    CALL remap_basal_conditions( mesh_old, mesh_new, map, ice)

    ! Update the three elevation differences
    CALL reallocate_bounds( ice%Hs_a,  mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( ice%dHi_a, mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( ice%dHs_a, mesh_new%vi1, mesh_new%vi2)

    DO vi = mesh_new%vi1, mesh_new%vi2
      ice%Hs_a( vi) = surface_elevation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))
    END DO

    ice%dHi_a( mesh_new%vi1:mesh_new%vi2) = ice%Hi_a( mesh_new%vi1:mesh_new%vi2) - refgeo_PD%Hi( mesh_new%vi1:mesh_new%vi2)
    ice%dHs_a( mesh_new%vi1:mesh_new%vi2) = ice%Hs_a( mesh_new%vi1:mesh_new%vi2) - refgeo_PD%Hs( mesh_new%vi1:mesh_new%vi2)
    ice%dHb_a( mesh_new%vi1:mesh_new%vi2) = ice%Hb_a( mesh_new%vi1:mesh_new%vi2) - refgeo_PD%Hb( mesh_new%vi1:mesh_new%vi2)

    ! Simple memory reallocation for all the rest
    ! ===========================================

    ! Basic data - surface height, regional sea level, and thickness above floatation
    CALL reallocate_bounds ( ice%Hs_a         , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%SL_a         , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%TAF_a        , mesh_new%vi1, mesh_new%vi2 )

    ! Different masks
    CALL reallocate_bounds ( ice%mask_land_a  ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%mask_ocean_a ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%mask_lake_a  ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%mask_ice_a   ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%mask_sheet_a ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%mask_shelf_a ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%mask_coast_a ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%mask_margin_a,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%mask_gl_a    ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%mask_cf_a    ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%mask_a       ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%f_grnd_a     , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%f_grnd_b     , mesh_new%ti1, mesh_new%ti2 )

    ! Ice physical properties
    CALL reallocate_bounds ( ice%A_flow_3D_a  , mesh_new%vi1, mesh_new%vi2, C%nz )
    CALL reallocate_bounds ( ice%A_flow_vav_a , mesh_new%vi1, mesh_new%vi2       )
    CALL reallocate_bounds ( ice%Ti_pmp_a     , mesh_new%vi1, mesh_new%vi2, C%nz )
    CALL reallocate_bounds ( ice%Cpi_a        , mesh_new%vi1, mesh_new%vi2, C%nz )
    CALL reallocate_bounds ( ice%Ki_a         , mesh_new%vi1, mesh_new%vi2, C%nz )

    ! Zeta derivatives
    CALL reallocate_bounds ( ice%dzeta_dt_a   , mesh_new%vi1, mesh_new%vi2, C%nz )
    CALL reallocate_bounds ( ice%dzeta_dx_a   , mesh_new%vi1, mesh_new%vi2, C%nz )
    CALL reallocate_bounds ( ice%dzeta_dy_a   , mesh_new%vi1, mesh_new%vi2, C%nz )
    CALL reallocate_bounds ( ice%dzeta_dz_a   , mesh_new%vi1, mesh_new%vi2       )

    ! Ice dynamics - ice thickness calculation
    CALL reallocate_bounds ( ice%dVi_in       ,            1, mesh_new%nV , mesh_new%nC_mem )
    CALL reallocate_bounds ( ice%dHs_dt_a     , mesh_new%vi1, mesh_new%vi2       )

    ! Ice dynamics - predictor/corrector ice thickness update
    CALL reallocate_bounds ( ice%pc_tau       , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%pc_fcb       , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%pc_f1        , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%pc_f2        , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%pc_f3        , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%pc_f4        , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%Hi_old       , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%Hi_pred      , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%Hi_corr      , mesh_new%vi1, mesh_new%vi2 )

    ! Thermodynamics
    CALL reallocate_bounds ( ice%internal_heating_a  , mesh_new%vi1, mesh_new%vi2, C%nz )
    CALL reallocate_bounds ( ice%frictional_heating_a, mesh_new%vi1, mesh_new%vi2       )

    ! Mesh adaptation data
    CALL reallocate_bounds ( ice%surf_curv    ,            1, mesh_new%nV  )
    CALL reallocate_bounds ( ice%log_velocity , mesh_new%vi1, mesh_new%vi2 )

    ! GIA
    CALL reallocate_bounds ( ice%dHb_a        , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%dHb_dt_a     , mesh_new%vi1, mesh_new%vi2 )
    CALL reallocate_bounds ( ice%dSL_dt_a     , mesh_new%vi1, mesh_new%vi2 )

    ! Remap velocities
    ! ================

    ! If we're running with choice_ice_dynamics == "none",
    ! remap velocities pretending we are using the DIVA, just like we did
    ! during initialisation to get reasonable velocities for the thermodynamics.

    IF (C%choice_ice_dynamics == 'none') THEN
      C%choice_ice_dynamics = 'DIVA'
      CALL remap_velocities( mesh_old, mesh_new, map, ice)
      C%choice_ice_dynamics = 'none'
    ELSE
      CALL remap_velocities( mesh_old, mesh_new, map, ice)
    END IF

    ! Ice flow factor
    ! ===============

    ! ! Recalculate ice flow factor
    CALL calc_ice_rheology( mesh_new, ice, time)

    ! Update masks
    CALL update_general_ice_model_data( mesh_new, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ice_model

  SUBROUTINE map_geothermal_heat_flux_to_mesh( mesh, ice)

    USE data_types_module,          ONLY: type_remapping_latlon2mesh
    USE forcing_module,             ONLY: forcing
    USE mesh_mapping_module,        ONLY: create_remapping_arrays_glob_mesh, map_latlon2mesh_2D, deallocate_remapping_arrays_glob_mesh

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_geothermal_heat_flux_to_mesh'
    TYPE(type_remapping_latlon2mesh)                   :: map

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Calculate mapping arrays
    ! CALL create_remapping_arrays_glob_mesh( mesh, forcing%grid_ghf, map)

    ! ! Map global climate data to the mesh
    ! CALL map_latlon2mesh_2D( mesh, map, forcing%ghf_ghf, ice%GHF_a)

    ! ! Deallocate mapping arrays
    ! CALL deallocate_remapping_arrays_glob_mesh( map)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_geothermal_heat_flux_to_mesh

END MODULE ice_dynamics_module
